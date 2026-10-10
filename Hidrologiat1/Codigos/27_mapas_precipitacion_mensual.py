"""
Precipitación media de cada mes calendario sobre la cuenca El Tambo (12 mapas, escala discretizada).

Fuente: CHIRPS v2.0 mensual 0.05° (script 26), periodo común P_L–P_I (periodo_comun.py). Para cada mes
calendario se promedian todos los años del periodo en cada celda (mm/mes). El campo se interpola de forma
bilineal a una malla de 0.0025° (~280 m) solo para representarlo de forma continua; el detalle real
sigue siendo de 0.05°. La suma de los 12 mapas es la precipitación media anual (panel (d) del script 07b).
Escala: intervalos discretos, más finos en valores bajos para distinguir los meses secos; colores de
naranja (seco) a azul oscuro (húmedo), en el estilo del ejemplo del curso.
En cada mapa se indica el promedio de cuenca de CHIRPS (ponderado por área, script 26) y el de P_L
(CAMELS-PE, sin malla espacial), que es mayor: el mapa muestra el patrón espacial, no la magnitud de P_L.

Salida: Tarea 1/El_Tambo/resultados/climatologia/mapas_precipitacion_mensual_chirps.png
Requiere 06, 09, 10 y 26.
"""
from pathlib import Path

import numpy as np
import pandas as pd
import geopandas as gpd
import xarray as xr
import rasterio
from rasterio.enums import Resampling
from rasterio.features import rasterize
from rasterio.transform import from_origin
from scipy.interpolate import RegularGridInterpolator
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, LinearSegmentedColormap, LightSource

from periodo_comun import periodo_comun

ROOT = Path(__file__).resolve().parents[2] / "Tarea 1" / "El_Tambo"
DATOS = ROOT / "datos"
DEL = ROOT / "resultados" / "delimitacion"
SMES = ROOT / "resultados" / "series_mensuales"
OUT = ROOT / "resultados" / "climatologia"
MESES = ["Enero", "Febrero", "Marzo", "Abril", "Mayo", "Junio", "Julio", "Agosto", "Septiembre", "Octubre",
         "Noviembre", "Diciembre"]
CLASES = [0, 5, 10, 15, 20, 30, 40, 50, 60, 80, 100, 125, 150, 175, 200, 250]   # mm/mes
COLORES = ["#e8761c", "#f0a21f", "#f5cf1f", "#eef01c", "#c9e02a", "#8fbf3a", "#5e9f4d", "#4b7f5f",
           "#4d6466", "#4a4f6e", "#3b3a63"]
RES = 0.0025                                                     # malla de representación (°)
PAD = 0.03

INI, FIN = periodo_comun()
basin = gpd.read_file(DEL / "cuenca_el_tambo_glo30.gpkg").to_crs(4326)
main = gpd.read_file(DEL / "rio_principal_el_tambo.gpkg").to_crs(4326)
outlet = gpd.read_file(DEL / "punto_cierre_el_tambo.gpkg").to_crs(4326).iloc[[1]]
x0, y0, x1, y1 = basin.total_bounds
x0, y0, x1, y1 = x0 - PAD, y0 - PAD, x1 + PAD, y1 + PAD

# ------------------------------------------------------------------ climatología mensual por celda
ch = xr.open_dataset(DATOS / "chirps" / "CHIRPS_v2p0_mensual_El_Tambo.nc")["precipitation"].sel(time=slice(INI, FIN))
assert ch.notnull().all(), "CHIRPS: faltantes"
clim = ch.groupby("time.month").mean()                          # (month, lat, lon), mm/mes
n_anios = ch.groupby("time.month").count().isel(lat=0, lon=0).values

# malla fina y máscara de cuenca
lons = np.arange(x0 + RES / 2, x1, RES)
lats = np.arange(y1 - RES / 2, y0, -RES)
T = from_origin(lons[0] - RES / 2, lats[0] + RES / 2, RES, RES)
mask = rasterize([(basin.geometry.iloc[0], 1)], out_shape=(lats.size, lons.size), transform=T, fill=0).astype(bool)
LA, LO = np.meshgrid(lats, lons, indexing="ij")
pts = np.column_stack([LA.ravel(), LO.ravel()])
campos = []
for m in range(1, 13):
    f = RegularGridInterpolator((clim.lat.values, clim.lon.values), clim.sel(month=m).values, method="linear")
    campos.append(np.where(mask, f(pts).reshape(LA.shape), np.nan))
campos = np.array(campos)

# relieve sombreado (DEM GLO-30 remuestreado a la misma malla)
with rasterio.open(DATOS / "dem" / "GLO30_El_Tambo_EPSG4326.tif") as src:
    win = rasterio.windows.from_bounds(x0, y0, x1, y1, src.transform)
    dem = src.read(1, window=win, out_shape=(lats.size, lons.size), resampling=Resampling.average).astype(float)
dx = RES * 111320 * np.cos(np.deg2rad((y0 + y1) / 2))
hs = LightSource(315, 45).hillshade(dem, vert_exag=2, dx=dx, dy=RES * 110570)

# promedios de cuenca por mes
cam = pd.read_csv(SMES / "series_mensuales_camels_pe.csv", index_col="fecha", parse_dates=True)["P_L_mm_mes"].loc[INI:FIN]
chs = pd.read_csv(SMES / "chirps_cuenca_mensual.csv", index_col="fecha", parse_dates=True)["P_chirps_mm_mes"].loc[INI:FIN]
m_PL = cam.groupby(cam.index.month).mean()
m_CH = chs.groupby(chs.index.month).mean()

# ------------------------------------------------------------------ figura
cmap = LinearSegmentedColormap.from_list("seco_humedo", COLORES, N=len(CLASES) - 1)
assert np.nanmax(campos) < CLASES[-1], f"máximo {np.nanmax(campos):.0f} mm fuera de la escala"
norm = BoundaryNorm(CLASES, cmap.N)
ext = (lons[0] - RES / 2, lons[-1] + RES / 2, lats[-1] - RES / 2, lats[0] + RES / 2)
aspect = 1 / np.cos(np.deg2rad((y0 + y1) / 2))
# Filas = trimestres (DEF, MAM, JJA, SON); columnas = sus tres meses en orden
TRIM = [("DEF", [12, 1, 2]), ("MAM", [3, 4, 5]), ("JJA", [6, 7, 8]), ("SON", [9, 10, 11])]
fig, axs = plt.subplots(4, 3, figsize=(12, 17.5))
fig.subplots_adjust(left=0.09, right=0.80, top=0.93, bottom=0.055, wspace=0.06, hspace=0.12)
for fila, (trim, meses_t) in enumerate(TRIM):
    pl_t = m_PL[meses_t].sum()
    pos = axs[fila, 0].get_position()
    fig.text(pos.x0 - 0.015, (pos.y0 + pos.y1) / 2, f"{trim}\nP$_L$ {pl_t:.0f} mm ({100 * pl_t / m_PL.sum():.0f} %)",
             rotation=90, ha="center", va="center", fontsize=11, fontweight="bold")
ubic = {m: (f, c) for f, (_, ms) in enumerate(TRIM) for c, m in enumerate(ms)}
for m in range(12):
    ax = axs[ubic[m + 1]]
    ax.imshow(hs, cmap="gray", extent=ext, alpha=0.45)
    im = ax.imshow(campos[m], cmap=cmap, norm=norm, extent=ext, alpha=0.88, interpolation="nearest")
    basin.boundary.plot(ax=ax, color="black", linewidth=1.0)
    main.plot(ax=ax, color="#08306b", linewidth=0.9)
    outlet.plot(ax=ax, color="red", edgecolor="k", markersize=18, zorder=5)
    ax.set_xlim(x0, x1); ax.set_ylim(y0, y1); ax.set_aspect(aspect)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel(""); ax.set_ylabel("")
    ax.set_title(MESES[m], fontsize=11)
    ax.text(0.98, 0.02, f"CHIRPS {m_CH[m + 1]:.0f} mm\nP$_L$ {m_PL[m + 1]:.0f} mm", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=8, bbox=dict(fc="white", alpha=0.85, lw=0))
cax = fig.add_axes([0.84, 0.25, 0.03, 0.5])
cb = fig.colorbar(im, cax=cax, ticks=CLASES, spacing="uniform")
cb.ax.set_yticklabels([f"{c}" for c in CLASES], fontsize=9)
cb.ax.set_title("[mm/mes]", fontsize=10, pad=8)
fig.suptitle(f"Precipitación media mensual – cuenca del río Chicama hasta El Tambo\n"
             f"CHIRPS v2.0 (0.05°), periodo común {INI:%Y-%m} a {FIN:%Y-%m}", fontsize=13)
fig.text(0.03, 0.012,
         f"Media de cada mes calendario sobre {int(n_anios.min())}–{int(n_anios.max())} años; interpolación bilineal de "
         f"las celdas de 0.05° solo para la representación.\nRecuadro: promedio de cuenca de CHIRPS (ponderado por área) "
         f"y de P$_L$ (CAMELS-PE). En el total anual CHIRPS es menor ({m_CH.sum():.0f} vs {m_PL.sum():.0f} mm): el mapa muestra "
         f"el patrón espacial, no la magnitud de P$_L$.\nFilas: trimestres DEF, MAM, JJA y SON, con el acumulado de P$_L$ "
         f"de cuenca y su fracción del total anual. Suma de los 12 meses = precipitación media anual. "
         f"Script 27_mapas_precipitacion_mensual.py.",
         fontsize=7.5, color="#444", va="bottom")
fig.savefig(OUT / "mapas_precipitacion_mensual_chirps.png", dpi=200)
plt.close(fig)

for m in range(12):
    v = campos[m][mask]
    print(f"{MESES[m]:10s} píxeles {np.nanmin(v):6.1f}–{np.nanmax(v):6.1f} mm | cuenca CHIRPS {m_CH[m + 1]:6.1f} | "
          f"P_L {m_PL[m + 1]:6.1f}")
print(f"Suma de los promedios de cuenca: CHIRPS {m_CH.sum():.0f} | P_L {m_PL.sum():.0f} mm/año")

