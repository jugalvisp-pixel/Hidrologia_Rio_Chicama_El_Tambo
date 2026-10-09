"""
IMERG V07 Final mensual (GPM_3IMERGM.07, NASA GES DISC) para la cuenca El Tambo:
descarga optimizada por OPeNDAP, promedio ponderado por área, conversión a mm/mes y gráfica
de la serie mensual frente a P_L de CAMELS-PE en el periodo común.

Producto : GPM IMERG Final Precipitation L3 1 month 0.1° x 0.1° V07 (GPM_3IMERGM.07)
           doi:10.5067/GPM/IMERG/3B-MONTH/07 — variable `precipitation` (mm/hr, tasa media
           del mes; la versión Final incorpora pluviómetros GPCC).
Acceso manual (preferido): descargar desde https://disc.gsfc.nasa.gov/datasets/GPM_3IMERGM_07/summary
           → "Subset / Get Data" → OPeNDAP, 1998-01-01..2025-12-31, caja W -78.9 S -8.1
           E -78.1 N -7.2, variable `precipitation`, formato netCDF; copiar los archivos a
           Tarea 1/El_Tambo/datos/imerg/manual/ (también acepta los HDF5 globales originales).
Acceso automático (si la carpeta manual está vacía): requiere cuenta NASA Earthdata con la aplicación "NASA GESDISC DATA ARCHIVE"
           autorizada y un archivo ~/.netrc (Windows: %USERPROFILE%\\_netrc y \\.netrc) con
               machine urs.earthdata.nasa.gov login <usuario> password <clave>
           Las credenciales NO se guardan en el proyecto ni en el ZIP de entrega.
Optimiza : por OPeNDAP se pide solo el recorte de la cuenca (índices de lon/lat) de cada mes,
           unos pocos kB por archivo, en lugar del HDF5 global (~30 MB por mes).

Conversión: P[mm/mes] = tasa[mm/h] · 24 h/día · n_días del mes (tasa media mensual).
Promedio  : media ponderada por el área de intersección celda–cuenca (polígono GLO-30, UTM 17S).

Salidas:
  Tarea 1/El_Tambo/datos/imerg/GPM_3IMERGM.07_El_Tambo_recorte.nc
  Tarea 1/El_Tambo/datos/imerg/pesos_imerg_cuenca.csv
  Tarea 1/El_Tambo/resultados/series_mensuales/imerg_cuenca_mensual.csv
  Tarea 1/El_Tambo/resultados/series_mensuales/fig_P_L_vs_IMERG_mensual.png
"""
from pathlib import Path
import netrc
import os
import re
import sys
import time

import numpy as np
import pandas as pd
import geopandas as gpd
import xarray as xr
import requests
from shapely.geometry import box
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

ROOT = Path(__file__).resolve().parents[2] / "Tarea 1" / "El_Tambo"
BASIN = ROOT / "resultados" / "delimitacion" / "cuenca_el_tambo_glo30.gpkg"
SERIES = ROOT / "resultados" / "series_mensuales" / "series_mensuales_camels_pe_1998_2025.csv"
DL = ROOT / "datos" / "imerg"
MES_DIR = DL / "meses"
OUT = ROOT / "resultados" / "series_mensuales"
for p in (DL, MES_DIR, OUT):
    p.mkdir(parents=True, exist_ok=True)

COLL = "https://gpm1.gesdisc.eosdis.nasa.gov/data/GPM_L3/GPM_3IMERGM.07"
DAP = "https://gpm1.gesdisc.eosdis.nasa.gov/opendap/GPM_L3/GPM_3IMERGM.07"
URS = "urs.earthdata.nasa.gov"
RES = 0.1
INICIO, FIN = 1998, 2025


# ------------------------------------------------------------------ sesión Earthdata
class EarthdataSession(requests.Session):
    """Conserva la autenticación al ser redirigido a URS (receta de NASA Earthdata)."""

    def rebuild_auth(self, prepared_request, response):
        headers = prepared_request.headers
        url = prepared_request.url
        if "Authorization" in headers:
            orig = requests.utils.urlparse(response.request.url).hostname
            redir = requests.utils.urlparse(url).hostname
            if orig != redir and redir != URS and orig != URS:
                del headers["Authorization"]


def sesion():
    rutas = [Path.home() / ".netrc", Path.home() / "_netrc"]
    for r in rutas:
        if r.exists():
            try:
                auth = netrc.netrc(r).authenticators(URS)
            except netrc.NetrcParseError as e:
                sys.exit(f"Error leyendo {r}: {e}")
            if auth:
                s = EarthdataSession()
                s.auth = (auth[0], auth[2])
                return s
    sys.exit("No se encontraron credenciales de NASA Earthdata para urs.earthdata.nasa.gov en "
             f"{' ni '.join(map(str, rutas))}.\nCrear el archivo según la cabecera de este script "
             "y autorizar 'NASA GESDISC DATA ARCHIVE' en el perfil de Earthdata.")


# ------------------------------------------------------------------ recorte espacial
basin = gpd.read_file(BASIN).to_crs(4326)
x0, y0, x1, y1 = basin.total_bounds
# Malla IMERG: centros lon = -179.95 + 0.1 i (i = 0..3599), lat = -89.95 + 0.1 j (j = 0..1799)
i0, i1 = int(np.floor((x0 + 180) / RES)) - 1, int(np.floor((x1 + 180) / RES)) + 1
j0, j1 = int(np.floor((y0 + 90) / RES)) - 1, int(np.floor((y1 + 90) / RES)) + 1


def archivos_disponibles(s):
    """Lista los archivos mensuales V07 por año desde el directorio público de GES DISC."""
    files = {}
    for y in range(INICIO, FIN + 1):
        r = s.get(f"{COLL}/{y}/", timeout=60)
        r.raise_for_status()
        for name in sorted(set(re.findall(r'href="(3B-MO\.MS\.MRG\.3IMERG\.(\d{8})[^"]*\.HDF5)"', r.text))):
            files[pd.Timestamp(name[1][:6] + "01")] = (y, name[0])
    return files


def descargar_mes(s, fecha, y, name, reintentos=3):
    out = MES_DIR / f"imerg_{fecha:%Y%m}.nc4"
    if out.exists() and out.stat().st_size > 0:
        return out
    ce = (f"precipitation[0:0][{i0}:{i1}][{j0}:{j1}],lat[{j0}:{j1}],lon[{i0}:{i1}],time[0:0]")
    url = f"{DAP}/{y}/{name}.nc4?" + requests.utils.quote(ce, safe=",")   # restricción DAP2
    for k in range(1, reintentos + 1):
        try:
            r = s.get(url, timeout=120)
            if r.status_code == 401:
                sys.exit("Earthdata rechazó las credenciales (401). Revisar usuario/clave y que la "
                         "aplicación 'NASA GESDISC DATA ARCHIVE' esté autorizada.")
            r.raise_for_status()
            out.write_bytes(r.content)
            return out
        except requests.RequestException as e:      # la falla queda visible
            print(f"  {fecha:%Y-%m}: intento {k}/{reintentos} falló ({e})")
            time.sleep(5 * k)
    raise RuntimeError(f"No se pudo descargar {name}")


# ------------------------------------------------------------------ origen de los archivos
# Modo manual (preferido): archivos descargados desde el portal de GES DISC y copiados a
# Tarea 1/El_Tambo/datos/imerg/manual/. Se aceptan recortes netCDF (*.nc4, *.nc) generados con
# "Subset / Get Data" o los HDF5 globales originales (*.HDF5). Si la carpeta está vacía, se
# intenta la descarga automática por OPeNDAP con ~/.netrc.
MANUAL = DL / "manual"
MANUAL.mkdir(parents=True, exist_ok=True)
manual = sorted(p for p in MANUAL.iterdir() if p.suffix.lower() in (".nc4", ".nc", ".hdf5", ".h5"))
if manual:
    paths = manual
    print(f"Modo manual: {len(paths)} archivos en {MANUAL}")
else:
    s = sesion()
    disp = archivos_disponibles(s)
    print(f"Archivos IMERG V07 mensuales disponibles {INICIO}–{FIN}: {len(disp)} "
          f"({min(disp):%Y-%m} a {max(disp):%Y-%m})")
    paths = [descargar_mes(s, fecha, y, name) for fecha, (y, name) in sorted(disp.items())]
print(f"Tamaño total de los archivos de entrada: {sum(p.stat().st_size for p in paths)/2**20:.1f} MB")


def fecha_de(p):
    """Fecha del mes a partir del nombre GES DISC (3B-MO...3IMERG.YYYYMMDD-...) o imerg_YYYYMM."""
    m = re.search(r"3IMERG\.(\d{6})\d{2}-S", p.name) or re.search(r"imerg_(\d{6})", p.name)
    if not m:
        raise ValueError(f"No se reconoce la fecha en el nombre de archivo: {p.name}")
    if p.name.startswith("3B-") and not p.name.startswith("3B-MO"):
        raise ValueError(f"{p.name} no es un archivo mensual (3B-MO); revisar la colección descargada.")
    return pd.Timestamp(m.group(1) + "01")


# ------------------------------------------------------------------ armado y verificación
lon_min, lon_max = -179.95 + RES * i0 - 1e-3, -179.95 + RES * i1 + 1e-3
lat_min, lat_max = -89.95 + RES * j0 - 1e-3, -89.95 + RES * j1 + 1e-3
das, fechas = [], []
for p in paths:
    # Los HDF5 originales y los recortes "_subsetted.nc4" del portal guardan las variables en
    # el grupo "Grid"; los recortes OPeNDAP del modo automático las tienen en la raíz.
    ds = xr.open_dataset(p, decode_times=False)
    if "precipitation" not in ds:
        ds = xr.open_dataset(p, group="Grid", decode_times=False)
    da = ds["precipitation"]
    units = da.attrs.get("units", da.attrs.get("Units", ""))
    if isinstance(units, bytes):
        units = units.decode()
    if units.replace(" ", "").lower() not in ("mm/hr", "mm/h", "mmhr-1", "mm/hour"):
        raise ValueError(f"{p.name}: unidades inesperadas '{units}'; revisar conversión.")
    da = da.sel(lon=slice(lon_min, lon_max), lat=slice(lat_min, lat_max))
    if da.sizes["lon"] < (i1 - i0 + 1) or da.sizes["lat"] < (j1 - j0 + 1):
        raise ValueError(f"{p.name}: el recorte no cubre la cuenca (+1 celda). Descargar con la "
                         f"caja W {lon_min:.2f}, S {lat_min:.2f}, E {lon_max:.2f}, N {lat_max:.2f} o mayor.")
    f = fecha_de(p)
    fechas.append(f)
    das.append(da.isel(time=0, drop=True).expand_dims(time=[f]).transpose("time", "lat", "lon").load())
dup = pd.Series(fechas)[pd.Series(fechas).duplicated()]
if len(dup):
    raise ValueError(f"Meses duplicados en la entrada: {sorted(set(dup.dt.strftime('%Y-%m')))}")
P_rate = xr.concat(das, "time").sortby("time")
esperados = pd.date_range(P_rate.time.values.min(), P_rate.time.values.max(), freq="MS")
faltan = esperados.difference(pd.DatetimeIndex(P_rate.time.values))
print(f"IMERG: {P_rate.sizes['time']} meses, {pd.Timestamp(P_rate.time.values.min()):%Y-%m} a "
      f"{pd.Timestamp(P_rate.time.values.max()):%Y-%m} | meses faltantes dentro del rango: "
      f"{', '.join(f'{t:%Y-%m}' for t in faltan) if len(faltan) else 'ninguno'}")
# Coordenadas en float32: el residuo puede quedar cerca de 0 o de 1, se mide la distancia al entero.
for eje, c0 in (("lon", 179.95), ("lat", 89.95)):
    k = (P_rate[eje].values.astype(float) + c0) / RES
    assert np.allclose(k, np.round(k), atol=1e-2), f"malla {eje} no alineada a 0.1°"
P_rate.to_dataset(name="precipitation").to_netcdf(DL / "GPM_3IMERGM.07_El_Tambo_recorte.nc")

# Pesos: área de intersección de cada celda con la cuenca
lons, lats = P_rate.lon.values, P_rate.lat.values
g = gpd.GeoDataFrame(
    [{"lon": lo, "lat": la} for la in lats for lo in lons],
    geometry=[box(lo - RES / 2, la - RES / 2, lo + RES / 2, la + RES / 2) for la in lats for lo in lons],
    crs=4326).to_crs(32717)
cuenca = basin.to_crs(32717).geometry.iloc[0]
g["area_en_cuenca_km2"] = g.intersection(cuenca).area / 1e6
g["fraccion_celda"] = g["area_en_cuenca_km2"] / (g.area / 1e6)
g[g["area_en_cuenca_km2"] > 0].drop(columns="geometry").round(5).to_csv(DL / "pesos_imerg_cuenca.csv", index=False)
W = xr.DataArray(g.pivot(index="lat", columns="lon", values="area_en_cuenca_km2")
                 .reindex(index=lats, columns=lons).values,
                 coords={"lat": lats, "lon": lons}, dims=("lat", "lon"))
rate_cuenca = (P_rate * W).sum(("lat", "lon")) / W.where(P_rate.notnull()).sum(("lat", "lon"))
rate_cuenca = rate_cuenca.to_series()
dias = rate_cuenca.index.days_in_month
P_I = (rate_cuenca * 24 * dias).rename("P_I_mm_mes")              # mm/h → mm/mes
imerg = pd.DataFrame({"P_I_tasa_mm_h": rate_cuenca, "n_dias": dias, "P_I_mm_mes": P_I})
imerg.index.name = "fecha"
imerg.round(4).to_csv(OUT / "imerg_cuenca_mensual.csv")

# ------------------------------------------------------------------ periodo común y gráfica
cam = pd.read_csv(SERIES, index_col="fecha", parse_dates=True)["P_L_mm_mes"]
df = pd.concat([cam.rename("P_L"), P_I.rename("P_I")], axis=1)
comun = df.dropna()
ini, fin = comun.index.min(), comun.index.max()
df = df.loc[ini:fin]
n = len(comun)
d = comun["P_I"] - comun["P_L"]
print(f"\nPeriodo común: {ini:%Y-%m} a {fin:%Y-%m} | {n} pares mensuales")
print(f"Medias: P_L {comun['P_L'].mean():.1f} mm/mes | P_I {comun['P_I'].mean():.1f} mm/mes | "
      f"sesgo P_I−P_L {d.mean():.1f} | MAE {d.abs().mean():.1f} | RMSE {np.sqrt((d**2).mean()):.1f} | "
      f"r Pearson {comun.corr().iloc[0, 1]:.3f}")

fig, a1 = plt.subplots(figsize=(15, 5.5))
a1.plot(df.index, df["P_L"], color="#1b6ca8", lw=1.3, label="P$_L$ – CAMELS-PE (PISCOp v2.1)")
a1.plot(df.index, df["P_I"], color="#e4572e", lw=1.3, alpha=0.9,
        label="P$_I$ – IMERG V07 Final (promedio ponderado por área)")
a1.set_ylabel("Precipitación (mm/mes)")
a1.set_ylim(-10, df[["P_L", "P_I"]].max().max() * 1.18)        # margen para leyenda y métricas
a1.set_title("Cuenca del río Chicama hasta El Tambo – precipitación mensual en el periodo común "
             f"({ini:%Y-%m} a {fin:%Y-%m}, n = {n} meses)", loc="left")
a1.legend(loc="upper right", frameon=False)
a1.grid(alpha=0.3)
a1.text(0.01, 0.95, f"sesgo medio P$_I$−P$_L$ = {d.mean():+.1f} mm/mes · MAE = {d.abs().mean():.1f} · "
        f"RMSE = {np.sqrt((d**2).mean()):.1f} · r = {comun.corr().iloc[0, 1]:.2f}",
        transform=a1.transAxes, va="top", fontsize=9, bbox=dict(fc="white", alpha=0.8, lw=0))
a1.xaxis.set_major_locator(mdates.YearLocator(2))
a1.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
a1.set_xlim(ini - pd.Timedelta(days=20), fin + pd.Timedelta(days=20))
fig.text(0.01, 0.005, "Meses faltantes se muestran como vacíos (sin interpolar). P$_I$ convertido de "
         "mm/h a mm/mes con el número de días de cada mes.", fontsize=8)
fig.savefig(OUT / "fig_P_L_vs_IMERG_mensual.png", dpi=200, bbox_inches="tight")
print(f"Figura: {OUT / 'fig_P_L_vs_IMERG_mensual.png'}")
