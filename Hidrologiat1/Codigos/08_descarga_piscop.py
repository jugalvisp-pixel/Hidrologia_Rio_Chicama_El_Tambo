"""
Descarga optimizada de PISCOp v2.1 (SENAMHI) para la cuenca El Tambo, 1998–2025, y promedio
de cuenca ponderado por área.

Fuente: IRI/LDEO Climate Data Library, espejo oficial de PISCOp v2.1
  https://iridl.ldeo.columbia.edu/SOURCES/.SENAMHI/.HSR/.PISCO/.Prec/.v2p1/
Optimización: el recorte espacial (rectángulo de la cuenca + 1 celda) y temporal se hace en el
servidor con la sintaxis de la Data Library (X/Y/T ... RANGEEDGES), de modo que cada producto
se descarga en una sola petición de pocos cientos de kB en lugar del NetCDF nacional (~GB).

Productos (verificado 2026-10-07):
  - stable/monthly   1981-01 .. 2016-12   mm/month  (versión publicada, Aybar et al. 2020)
  - unstable/monthly 1981-01 .. 2026-09   mm/month  (versión operacional; única que cubre 2017+)
  - stable/daily     1981-01-01 .. 2016-12-31  mm/day (no existe "unstable" diario en IRI)

Salidas:
  Tarea 1/El_Tambo/datos/piscop/*.nc                       recortes descargados
  Tarea 1/El_Tambo/datos/piscop/pesos_piscop_cuenca.csv    fracción de cada celda dentro de la cuenca
  Tarea 1/El_Tambo/resultados/precipitacion/piscop_cuenca_mensual.csv
"""
from pathlib import Path
import time
import urllib.request

import numpy as np
import pandas as pd
import geopandas as gpd
import xarray as xr
from shapely.geometry import box

ROOT = Path(__file__).resolve().parents[2] / "Tarea 1"
BASIN = ROOT / "El_Tambo" / "resultados" / "delimitacion" / "cuenca_el_tambo_glo30.gpkg"
CAMELS_TS = ROOT / "CAMELS_PE" / "CAMELS-PE" / "03_timeseries" / "by_catchment" / "PE_47E0A2A8.csv"
DL = ROOT / "El_Tambo" / "datos" / "piscop"
OUT = ROOT / "El_Tambo" / "resultados" / "precipitacion"
DL.mkdir(parents=True, exist_ok=True)
OUT.mkdir(parents=True, exist_ok=True)

IRI = "https://iridl.ldeo.columbia.edu/SOURCES/.SENAMHI/.HSR/.PISCO/.Prec/.v2p1"
RES = 0.1
PRODUCTS = {  # nombre: (ruta IRI, rango temporal en sintaxis IRI)
    "piscop_v2p1_stable_monthly": (".stable/.monthly", "(Jan 1998)/(Dec 2016)"),
    "piscop_v2p1_unstable_monthly": (".unstable/.monthly", "(Jan 1998)/(Dec 2025)"),
    "piscop_v2p1_stable_daily": (".stable/.daily", "(1 Jan 1998)/(31 Dec 2016)"),
}

basin = gpd.read_file(BASIN).to_crs(4326)
geom = basin.geometry.iloc[0]
x0, y0, x1, y1 = geom.bounds
# Rectángulo de descarga: límites de la cuenca + una celda de margen, alineado a la malla
X0, X1 = np.floor(x0 / RES) * RES - RES, np.ceil(x1 / RES) * RES + RES
Y0, Y1 = np.floor(y0 / RES) * RES - RES, np.ceil(y1 / RES) * RES + RES


def iri_url(path, trange):
    q = (f"{IRI}/{path}/.Prec/X/({X0:.2f})/({X1:.2f})/RANGEEDGES/"
         f"Y/({Y0:.2f})/({Y1:.2f})/RANGEEDGES/T/{trange}/RANGEEDGES/data.nc")
    return q.replace(" ", "%20").replace("(", "%28").replace(")", "%29")


def download(name, path, trange, retries=3):
    out = DL / f"{name}_1998.nc"
    if out.exists():
        print(f"ya existe: {out.name}")
        return out
    url = iri_url(path, trange)
    for k in range(1, retries + 1):
        try:
            t0 = time.time()
            urllib.request.urlretrieve(url, out)
            print(f"{out.name}: {out.stat().st_size/1024:.0f} kB en {time.time()-t0:.1f} s")
            return out
        except Exception as e:           # la falla queda visible; no se sustituye por datos sintéticos
            print(f"intento {k}/{retries} falló para {name}: {e}")
            out.unlink(missing_ok=True)
            time.sleep(5 * k)
    raise RuntimeError(f"No se pudo descargar {name}:\n{url}")


def open_piscop(f):
    """Abre el recorte decodificando el tiempo (los mensuales usan calendario 360 días)."""
    ds = xr.open_dataset(f, decode_times=False)
    t = ds["T"]
    units = t.attrs["units"]
    if units.startswith("months since 1960-01-01"):
        m = np.floor(t.values).astype(int)
        time_index = pd.to_datetime({"year": 1960 + m // 12, "month": m % 12 + 1, "day": 1})
    elif units.startswith("days since 1960-01-01"):
        time_index = pd.Timestamp("1960-01-01") + pd.to_timedelta(np.floor(t.values), unit="D")
    else:
        raise ValueError(f"Unidades de tiempo no previstas: {units}")
    da = ds["Prec"].assign_coords(T=pd.DatetimeIndex(time_index)).rename(T="time", X="lon", Y="lat")
    return da


files = {n: download(n, p, t) for n, (p, t) in PRODUCTS.items()}

# ------------------------------------------------------------------ pesos por fracción de área
ref = open_piscop(files["piscop_v2p1_unstable_monthly"])
lons, lats = ref["lon"].values, ref["lat"].values
cells = [(lo, la, box(lo - RES / 2, la - RES / 2, lo + RES / 2, la + RES / 2))
         for la in lats for lo in lons]
g = gpd.GeoDataFrame({"lon": [c[0] for c in cells], "lat": [c[1] for c in cells]},
                     geometry=[c[2] for c in cells], crs=4326).to_crs(32717)
basin_utm = basin.to_crs(32717).geometry.iloc[0]
g["area_celda_km2"] = g.area / 1e6
g["area_en_cuenca_km2"] = g.intersection(basin_utm).area / 1e6
g["fraccion"] = g["area_en_cuenca_km2"] / g["area_celda_km2"]
w = g[g["area_en_cuenca_km2"] > 0].drop(columns="geometry")
w.round(5).to_csv(DL / "pesos_piscop_cuenca.csv", index=False)
print(f"\nCeldas PISCOp que tocan la cuenca: {len(w)} | completas (≥99 %): {(w['fraccion'] >= .99).sum()}"
      f" | área ponderada {w['area_en_cuenca_km2'].sum():.1f} km2")
W = xr.DataArray(
    g.pivot(index="lat", columns="lon", values="area_en_cuenca_km2").reindex(index=lats, columns=lons).values,
    coords={"lat": lats, "lon": lons}, dims=("lat", "lon"))


def basin_mean(da):
    """Promedio ponderado por área de intersección; ignora celdas sin dato renormalizando pesos."""
    valid = da.notnull()
    return (da * W).sum(("lat", "lon")) / (W.where(valid)).sum(("lat", "lon"))


# ------------------------------------------------------------------ series de cuenca
st_m = basin_mean(open_piscop(files["piscop_v2p1_stable_monthly"])).to_series()
un_m = basin_mean(open_piscop(files["piscop_v2p1_unstable_monthly"])).to_series()
st_d = basin_mean(open_piscop(files["piscop_v2p1_stable_daily"])).to_series()
st_d2m = st_d.resample("MS").sum(min_count=1)

cam = pd.read_csv(CAMELS_TS, parse_dates=["date"]).set_index("date")["prec"].loc["1998":"2025"]
cam_m = cam.resample("MS").sum(min_count=1)

df = pd.DataFrame({"P_piscop_unstable_mensual": un_m, "P_piscop_stable_mensual": st_m,
                   "P_piscop_stable_diario_agregado": st_d2m, "P_camels_pe": cam_m})
df.index.name = "fecha"
df = df.loc["1998-01":"2025-12"]
df.round(2).to_csv(OUT / "piscop_cuenca_mensual.csv")

# ------------------------------------------------------------------ verificaciones
print("\n=== Cobertura (meses con dato, 1998-01..2025-12 = 336) ===")
print(df.notna().sum().to_string())


def comp(a, b, name):
    x = df[[a, b]].dropna()
    d = x[a] - x[b]
    print(f"{name:52s} n={len(x):3d} | sesgo={d.mean():6.2f} mm/mes | MAE={d.abs().mean():6.2f}"
          f" | máx|dif|={d.abs().max():6.2f} | r={x[a].corr(x[b]):.4f}")


print("\n=== Comparaciones de consistencia ===")
comp("P_piscop_stable_diario_agregado", "P_piscop_stable_mensual", "diario agregado vs mensual (stable)")
comp("P_piscop_unstable_mensual", "P_piscop_stable_mensual", "unstable vs stable (1998-2016)")
comp("P_piscop_unstable_mensual", "P_camels_pe", "unstable (este trabajo) vs CAMELS-PE")
comp("P_piscop_stable_mensual", "P_camels_pe", "stable (este trabajo) vs CAMELS-PE")
ann = df.resample("YS").sum(min_count=12).round(0)
ann.index = ann.index.year
print("\n=== Totales anuales (mm) ===")
print(ann.to_string())
print(f"\nSalidas: {DL} y {OUT / 'piscop_cuenca_mensual.csv'}")
