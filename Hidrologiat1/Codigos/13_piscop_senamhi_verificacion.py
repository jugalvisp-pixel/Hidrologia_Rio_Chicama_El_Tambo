"""
PISCOp oficial de SENAMHI (figshare, con DOI) para la cuenca El Tambo: descarga con verificación
MD5 y comprobación de que el promedio de cuenca reproduce la `prec` de CAMELS-PE, antes de usar la
malla para el mapa de precipitación media anual.

Productos (archivos mensuales; bastan para totales mensuales y anuales):
  v2.1 update : Gutierrez, Aybar, Endara, Lavado-Casimiro (2022), 1981–2019
                doi:10.6084/m9.figshare.21127423.v2 — PISCOp_monthly.nc
  v3.0        : Gutierrez, Lavado-Casimiro (2026), 1981–2025
                doi:10.6084/m9.figshare.32411886.v1 — PISCOp_m.nc y PISCOp_d.nc (diario, 1.5 GB)
El PISCOp del IRI se descartó (decisión del grupo, 2026-10-07): no reproduce CAMELS-PE.

Salidas: Tarea 1/El_Tambo/datos/piscop_senamhi/<version>/  (archivos descargados y PROCEDENCIA.json)
"""
from pathlib import Path
import hashlib
import io
import json
import sys
import time
import zipfile

import numpy as np
import pandas as pd
import geopandas as gpd
import xarray as xr
import requests
from shapely.geometry import box

ROOT = Path(__file__).resolve().parents[2] / "Tarea 1" / "El_Tambo"
DL = ROOT / "datos" / "piscop_senamhi"
OUT = ROOT / "resultados" / "precipitacion"
ZIP = ROOT.parent / "CAMELS_PE" / "CAMELS-PE_v1.0.1.zip"
BASIN = ROOT / "resultados" / "delimitacion" / "cuenca_el_tambo_glo30.gpkg"
RES = 0.1

PRODUCTOS = {
    "v2p1_update": {
        "doi": "10.6084/m9.figshare.21127423.v2",
        "archivos": {"PISCOp_monthly.nc": ("https://ndownloader.figshare.com/files/37488715",
                                           "038d625efbec91048ea1fe1210503d16")},
    },
    "v3p0": {
        "doi": "10.6084/m9.figshare.32411886.v1",
        "archivos": {"PISCOp_m.nc": ("https://ndownloader.figshare.com/files/64968111",
                                     "fbf4f19b5d537183a75d08615c404010"),
                     "readme.txt": ("https://ndownloader.figshare.com/files/64968015",
                                    "f8940f5331d9710a95e135cd0d2c4b56"),
                     # diario nacional (1.5 GB): solo para extraer el recorte de la cuenca (script 14);
                     # no va en el ZIP de entrega
                     "PISCOp_d.nc": ("https://ndownloader.figshare.com/files/64968183",
                                     "01a7259f4a4fa38158496e3670bb1794")},
    },
}


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(2**20), b""):
            h.update(b)
    return h.hexdigest()


def descargar(url, out, md5_ok, reintentos=3):
    if out.exists() and md5(out) == md5_ok:
        print(f"ya existe y MD5 correcto: {out.name}")
        return
    for k in range(1, reintentos + 1):
        try:
            t0 = time.time()
            with requests.get(url, stream=True, timeout=120) as r:
                r.raise_for_status()
                with open(out, "wb") as f:
                    for b in r.iter_content(2**20):
                        f.write(b)
            got = md5(out)
            if got != md5_ok:
                raise IOError(f"MD5 {got} ≠ {md5_ok}")
            print(f"{out.name}: {out.stat().st_size / 2**20:.1f} MB en {time.time() - t0:.0f} s, MD5 verificado")
            return
        except (requests.RequestException, IOError) as e:     # la falla queda visible
            print(f"  intento {k}/{reintentos} falló para {out.name}: {e}")
            out.unlink(missing_ok=True)
            time.sleep(5 * k)
    sys.exit(f"No se pudo descargar {url}")


for ver, info in PRODUCTOS.items():
    d = DL / ver
    d.mkdir(parents=True, exist_ok=True)
    for nombre, (url, h) in info["archivos"].items():
        descargar(url, d / nombre, h)
    (d / "PROCEDENCIA.json").write_text(json.dumps(
        {"producto": f"PISCOp {ver}", "doi": info["doi"], "proveedor": "SENAMHI (figshare)",
         "archivos": {n: {"url": u, "md5": h} for n, (u, h) in info["archivos"].items()},
         "descargado": time.strftime("%Y-%m-%d")}, indent=2, ensure_ascii=False), encoding="utf-8")


# ------------------------------------------------------------------ lectura de las mallas
def abrir(ver):
    """Malla mensual (mm/mes) con índice de fechas; el tiempo viene como n.º de mes desde 1981-01."""
    if ver == "v2p1_update":
        ds = xr.open_dataset(DL / ver / "PISCOp_monthly.nc", decode_times=False)
        da = ds["pc"].rename(z="time")
    else:
        ds = xr.open_dataset(DL / ver / "PISCOp_m.nc", decode_times=False)
        da = ds["precipitation"].rename(Z1="time", latitude="lat", longitude="lon")
    fechas = pd.date_range("1981-01-01", periods=da.sizes["time"], freq="MS")
    return da.assign_coords(time=fechas).sortby("lat")


with zipfile.ZipFile(ZIP) as z:
    cat = gpd.read_file(io.BytesIO(z.read("CAMELS-PE/04_geospatial/camels_pe_catchments.gpkg")))
    cat = cat.set_index("gauge_id").to_crs(4326)
    cam = {i: pd.read_csv(io.BytesIO(z.read(f"CAMELS-PE/03_timeseries/by_catchment/{i}.csv")),
                          parse_dates=["date"], index_col="date")["prec"].resample("MS").sum(min_count=28)
           for i in ("PE_47E0A2A8", "PE_472E74CC")}
POLIGONOS = {"GLO-30 (este trabajo)": ("PE_47E0A2A8", gpd.read_file(BASIN).to_crs(4326).geometry.iloc[0]),
             "CAMELS-PE El Tambo": ("PE_47E0A2A8", cat.loc["PE_47E0A2A8", "geometry"]),
             "CAMELS-PE Puente Coina": ("PE_472E74CC", cat.loc["PE_472E74CC", "geometry"])}


def pesos(da, geom, metodo):
    """Pesos por celda: área de intersección (UTM 17S) o 1 si el centro cae dentro del polígono."""
    x0, y0, x1, y1 = geom.bounds
    sub = da.sel(lon=slice(x0 - RES, x1 + RES), lat=slice(y0 - RES, y1 + RES))
    lons, lats = sub.lon.values, sub.lat.values
    g = gpd.GeoDataFrame({"lon": np.tile(lons, len(lats)), "lat": np.repeat(lats, len(lons))},
                         geometry=[box(lo - RES / 2, la - RES / 2, lo + RES / 2, la + RES / 2)
                                   for la in lats for lo in lons], crs=4326)
    if metodo == "área":
        g = g.to_crs(32717)
        w = g.intersection(gpd.GeoSeries([geom], crs=4326).to_crs(32717).iloc[0]).area.values
    else:
        w = gpd.points_from_xy(g["lon"], g["lat"]).within(geom).astype(float)
    return sub, xr.DataArray(w.reshape(len(lats), len(lons)), coords={"lat": lats, "lon": lons},
                             dims=("lat", "lon"))


def media_cuenca(sub, W):
    return ((sub * W).sum(("lat", "lon")) / W.where(sub.notnull()).sum(("lat", "lon"))).to_series()


filas = []
for ver in PRODUCTOS:
    da = abrir(ver)
    for pol, (gid, geom) in POLIGONOS.items():
        for metodo in ("área", "centro"):
            p = media_cuenca(*pesos(da, geom, metodo))
            for ini, fin in (("1998", "2019"), ("1998", "2025"), ("1981", "2019")):
                x = pd.concat([p.rename("pisco"), cam[gid].rename("camels")], axis=1).loc[ini:fin].dropna()
                if x.index.max().year < int(fin):
                    continue
                d = x["pisco"] - x["camels"]
                filas.append({"producto": ver, "poligono": pol, "pesos": metodo, "periodo": f"{ini}-{fin}",
                              "n_meses": len(x), "P_camels_mm_anio": x["camels"].mean() * 12,
                              "P_pisco_mm_anio": x["pisco"].mean() * 12,
                              "dif_pct": 100 * d.mean() / x["camels"].mean(),
                              "MAE_mm_mes": d.abs().mean(), "max_abs_dif_mm_mes": d.abs().max(),
                              "r": x.corr().iloc[0, 1]})
res = pd.DataFrame(filas)
OUT.mkdir(parents=True, exist_ok=True)
res.round(3).to_csv(OUT / "verificacion_piscop_senamhi_vs_camels.csv", index=False)
pd.set_option("display.width", 250)
print("\n=== Promedio de cuenca de PISCOp (SENAMHI) vs `prec` de CAMELS-PE ===")
print(res.round(2).to_string(index=False))
