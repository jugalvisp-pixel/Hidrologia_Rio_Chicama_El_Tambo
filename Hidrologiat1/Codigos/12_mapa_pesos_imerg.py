"""
Mapa de las celdas IMERG (0.1°) sobre la cuenca El Tambo con el peso asignado a cada una en el
promedio de cuenca, y relieve del terreno (GLO-30) dentro de la porción de cada celda en la cuenca.

  (a) Peso de cada celda = área de intersección celda–cuenca / área de la cuenca (%).
  (b) Desnivel interno (P95 − P5 de la elevación) y elevación media de cada celda dentro de la
      cuenca: muestra cuánta variación topográfica queda dentro de un solo valor de IMERG.

El peso es solo por área (no considera la topografía); ver 10_imerg_mensual_y_grafica.py.
Requiere 06 y 10.

Salidas (Tarea 1/El_Tambo/resultados/series_mensuales/):
  mapa_celdas_imerg_pesos.png
  relieve_celdas_imerg.csv
"""
from pathlib import Path

import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
from rasterio.features import rasterize
from rasterio.windows import from_bounds
from shapely.geometry import box
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource
from matplotlib.lines import Line2D

from utils_mapas import formato_grados, escala_y_norte, m_por_grado_lon, M_PER_DEG_Y

ROOT = Path(__file__).resolve().parents[2] / "Tarea 1" / "El_Tambo"
DATOS = ROOT / "datos"
DEL = ROOT / "resultados" / "delimitacion"
OUT = ROOT / "resultados" / "series_mensuales"
RES = 0.1

basin = gpd.read_file(DEL / "cuenca_el_tambo_glo30.gpkg").to_crs(4326)
main = gpd.read_file(DEL / "rio_principal_el_tambo.gpkg").to_crs(4326)
outlet = gpd.read_file(DEL / "punto_cierre_el_tambo.gpkg").to_crs(4326).iloc[[1]]
w = pd.read_csv(DATOS / "imerg" / "pesos_imerg_cuenca.csv")
w["peso_pct"] = 100 * w["area_en_cuenca_km2"] / w["area_en_cuenca_km2"].sum()
w["completa"] = w["fraccion_celda"] > 0.99          # completa: ≥ 99 % de la celda en la cuenca
cells = gpd.GeoDataFrame(w, crs=4326, geometry=[box(lo - RES / 2, la - RES / 2, lo + RES / 2, la + RES / 2)
                                                for lo, la in zip(w["lon"], w["lat"])])

# Mapa con margen de media celda alrededor de la malla usada
xmin, ymin, xmax, ymax = cells.total_bounds + np.array([-0.02, -0.02, 0.02, 0.02])
LAT0 = (ymin + ymax) / 2

with rasterio.open(DATOS / "dem" / "GLO30_El_Tambo_EPSG4326.tif") as src:
    win = from_bounds(xmin, ymin, xmax, ymax, transform=src.transform).round_offsets().round_lengths()
    dem = src.read(1, window=win).astype(float)
    T = src.window_transform(win)
    if src.nodata is not None:
        dem[dem == src.nodata] = np.nan
EXT = (T.c, T.c + T.a * dem.shape[1], T.f + T.e * dem.shape[0], T.f)
DX, DY = abs(T.a) * m_por_grado_lon(LAT0), abs(T.e) * M_PER_DEG_Y
HS = LightSource(315, 45).hillshade(np.nan_to_num(dem, nan=np.nanmean(dem)), vert_exag=2, dx=DX, dy=DY)
in_basin = rasterize([(basin.geometry.iloc[0], 1)], out_shape=dem.shape, transform=T, fill=0).astype(bool)

# Relieve de la porción de cada celda dentro de la cuenca
rel = []
for k, c in cells.iterrows():
    m = rasterize([(c.geometry, 1)], out_shape=dem.shape, transform=T, fill=0).astype(bool) & in_basin
    z = dem[m & ~np.isnan(dem)]
    rel.append({"lon": c["lon"], "lat": c["lat"], "z_media_m": z.mean(), "z_min_m": z.min(),
                "z_max_m": z.max(), "desnivel_p5_p95_m": np.percentile(z, 95) - np.percentile(z, 5)})
rel = pd.DataFrame(rel)
cells = cells.merge(rel, on=["lon", "lat"])
cells.drop(columns="geometry").round(3).to_csv(OUT / "relieve_celdas_imerg.csv", index=False)


def marco(ax, title):
    ax.imshow(HS, cmap="gray", extent=EXT, alpha=0.55)
    ax.set_xlim(xmin, xmax); ax.set_ylim(ymin, ymax)
    ax.set_title(title, fontsize=11, loc="left", fontweight="bold")
    formato_grados(ax, LAT0)
    ax.tick_params(labelsize=7.5)


def contorno(ax):
    cells.boundary.plot(ax=ax, color="#444444", linewidth=0.6)
    main.plot(ax=ax, color="#08306b", linewidth=1.4)
    basin.boundary.plot(ax=ax, color="black", linewidth=1.6)
    outlet.plot(ax=ax, color="red", edgecolor="k", markersize=45, zorder=6)
    escala_y_norte(ax, LAT0)


def rotulos(ax, col, fmt, oscuro):
    """Valor de cada celda centrado en ella; texto blanco sobre fondos oscuros."""
    for _, c in cells.iterrows():
        ax.text(c["lon"], c["lat"], fmt(c), ha="center", va="center", fontsize=7, zorder=7,
                color="white" if c[col] > oscuro else "#1a1a1a",
                fontweight="bold" if c["completa"] else "normal")


fig, (a1, a2) = plt.subplots(1, 2, figsize=(14, 10))
fig.subplots_adjust(left=0.05, right=0.98, top=0.87, bottom=0.12, wspace=0.08)

marco(a1, "(a) Peso de cada celda en el promedio de cuenca")
cells.plot(ax=a1, column="peso_pct", cmap="Blues", vmin=0, vmax=cells["peso_pct"].max(),
           alpha=0.8, edgecolor="none")
contorno(a1)
rotulos(a1, "peso_pct", lambda c: f"{c['peso_pct']:.1f}", 0.6 * cells["peso_pct"].max())
cb = fig.colorbar(a1.collections[0], ax=a1, orientation="horizontal", fraction=0.04, pad=0.08)
cb.set_label("Peso (%) = área de la celda dentro de la cuenca / área de la cuenca", fontsize=8.5)
nb = (~cells["completa"]).sum()
pb = cells.loc[~cells["completa"], "peso_pct"].sum()
a1.text(0.01, 0.01, f"{len(cells)} celdas · {cells['completa'].sum()} completas (rótulo en negrita) · "
        f"{nb} de borde = {pb:.1f} % del peso\nCeldas equivalentes (Σ fracciones) = "
        f"{cells['fraccion_celda'].sum():.1f}", transform=a1.transAxes, fontsize=8, va="bottom",
        bbox=dict(fc="white", alpha=0.85, lw=0), zorder=9)

marco(a2, "(b) Desnivel del terreno dentro de cada celda (GLO-30)")
cells.plot(ax=a2, column="desnivel_p5_p95_m", cmap="Oranges", vmin=0,
           vmax=cells["desnivel_p5_p95_m"].max(), alpha=0.8, edgecolor="none")
contorno(a2)
rotulos(a2, "desnivel_p5_p95_m", lambda c: f"{c['desnivel_p5_p95_m']:.0f}\n({c['z_media_m']:.0f})",
        0.6 * cells["desnivel_p5_p95_m"].max())
cb = fig.colorbar(a2.collections[0], ax=a2, orientation="horizontal", fraction=0.04, pad=0.08)
cb.set_label("Desnivel P95 − P5 de la elevación dentro de la cuenca (m)\n"
             "entre paréntesis: elevación media de la celda en la cuenca (m s.n.m.)", fontsize=8.5)

fig.legend(handles=[Line2D([], [], color="black", lw=1.6, label="Límite de cuenca (GLO-30)"),
                    Line2D([], [], color="#444444", lw=0.6, label="Celda IMERG V07 (0.1°)"),
                    Line2D([], [], color="#08306b", lw=1.4, label="Río Chicama"),
                    Line2D([], [], color="red", marker="o", lw=0, markeredgecolor="k",
                           label="Punto de cierre El Tambo")],
           loc="lower center", ncol=4, fontsize=9, frameon=False, bbox_to_anchor=(0.5, 0.0))
fig.suptitle("Celdas de IMERG V07 sobre la cuenca del río Chicama hasta El Tambo y pesos del "
             "promedio de cuenca", fontsize=13, y=0.985)
fig.text(0.5, 0.945, "Coordenadas geográficas WGS 84 (EPSG:4326) · áreas de intersección calculadas "
         "en UTM 17S (EPSG:32717)\nEl peso solo considera el área, no la topografía",
         ha="center", va="top", fontsize=9)
fig.savefig(OUT / "mapa_celdas_imerg_pesos.png", dpi=200)

print(cells[["lon", "lat", "fraccion_celda", "peso_pct", "z_media_m", "desnivel_p5_p95_m"]]
      .round(2).sort_values("peso_pct", ascending=False).to_string(index=False))
print(f"\nCeldas: {len(cells)} | completas: {cells['completa'].sum()} | peso en celdas de borde: {pb:.1f} %")
print(f"Desnivel interno P95−P5: mediana {cells['desnivel_p5_p95_m'].median():.0f} m, "
      f"máx {cells['desnivel_p5_p95_m'].max():.0f} m")
print(f"Figura: {OUT / 'mapa_celdas_imerg_pesos.png'}")
