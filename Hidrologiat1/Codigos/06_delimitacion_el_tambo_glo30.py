"""
Delimitación de la cuenca del río Chicama hasta la estación El Tambo (CAMELS-PE PE_47E0A2A8)
a partir del modelo digital de elevación Copernicus GLO-30 (DSM, 1 arcsec ~ 30 m).

Pasos:
  1. Mosaico de tiles GLO-30 y recorte a la zona de la cuenca (+ margen).
  2. Reproyección a WGS 84 / UTM 17S (EPSG:32717) a 30 m (bilineal), para trabajar en metros.
  3. Acondicionamiento hidrológico: relleno de sumideros y depresiones, resolución de zonas planas.
  4. Direcciones de flujo D8 y acumulación de flujo (pysheds).
  5. Ajuste del punto de cierre (coordenadas SENAMHI) a la celda del cauce cuya área
     drenada se aproxima al área reportada por CAMELS-PE.
  6. Cuenca, red de drenaje (umbral de área), orden de Strahler y río principal (camino de flujo más largo).
  7. Parámetros morfométricos, comparación con el polígono CAMELS-PE y mapa.

Entradas : Tarea 1/El_Tambo/datos/dem/GLO30_tiles/*.tif
Salidas  : Tarea 1/El_Tambo/resultados/delimitacion/
"""
from pathlib import Path
import json
import re
import urllib.request

import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
from rasterio.merge import merge
from rasterio.warp import calculate_default_transform, reproject, Resampling
from rasterio.features import shapes, rasterize
from rasterio.windows import from_bounds
from shapely.geometry import shape, Point, LineString
from shapely.ops import unary_union
from pysheds.grid import Grid
from pysheds.view import Raster
# pysheds 0.5 usa np.in1d, eliminado en NumPy 2.4+; np.isin es el reemplazo equivalente.
if not hasattr(np, "in1d"):
    np.in1d = lambda a, b, **kw: np.isin(np.ravel(a), b, **kw)
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource
from matplotlib_scalebar.scalebar import ScaleBar
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.ticker import FuncFormatter, MultipleLocator

# ------------------------------------------------------------------ configuración
GID = "PE_47E0A2A8"
GAUGE_NAME = "El Tambo"
GAUGE_LON, GAUGE_LAT = -78.707, -7.574          # SENAMHI (stations.csv de CAMELS-PE)
AREA_CAMELS = 2171.585                           # km2, topographic_attributes.csv
BBOX = (-78.95, -8.15, -78.05, -7.20)            # lon_min, lat_min, lon_max, lat_max (cuenca + margen)
CRS_UTM = "EPSG:32717"
RES = 30.0                                       # m
STREAM_THR_KM2 = 5.0                             # área mínima drenada para definir cauce
SNAP_RADIUS_M = 1500.0

ROOT = Path(__file__).resolve().parents[2] / "Tarea 1"
TILES = ROOT / "El_Tambo" / "datos" / "dem" / "GLO30_tiles"
DEM_DIR = ROOT / "El_Tambo" / "datos" / "dem"
OUT = ROOT / "El_Tambo" / "resultados" / "delimitacion"
OUT.mkdir(parents=True, exist_ok=True)
CAMELS_GPKG = ROOT / "CAMELS_PE" / "CAMELS-PE" / "04_geospatial" / "camels_pe_catchments.gpkg"
NE_URL = ("https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/"
          "geojson/ne_50m_admin_0_countries.geojson")
NE_FILE = ROOT / "El_Tambo" / "datos" / "ne_50m_admin_0_countries.geojson"
SENAMHI_URL = "https://www.senamhi.gob.pe/mapas/mapa-estaciones-2/"

# ------------------------------------------------------------------ 1. mosaico y recorte
tiles = sorted(TILES.glob("Copernicus_DSM_COG_10_*_DEM.tif"))
srcs = [rasterio.open(t) for t in tiles]
for s in srcs:
    print(f"{s.name.split('/')[-1][-40:]}: bounds={tuple(round(b, 3) for b in s.bounds)} crs={s.crs}")
mosaic, mtrans = merge(srcs, bounds=BBOX)
meta = srcs[0].meta.copy()
meta.update(height=mosaic.shape[1], width=mosaic.shape[2], transform=mtrans, driver="GTiff",
            compress="deflate")
for s in srcs:
    s.close()
dem_geo = DEM_DIR / "GLO30_El_Tambo_EPSG4326.tif"
with rasterio.open(dem_geo, "w", **meta) as dst:
    dst.write(mosaic)

# ------------------------------------------------------------------ 2. reproyección a UTM 17S
dem_utm = DEM_DIR / "GLO30_El_Tambo_UTM17S_30m.tif"
with rasterio.open(dem_geo) as src:
    tr, w, h = calculate_default_transform(src.crs, CRS_UTM, src.width, src.height,
                                           *src.bounds, resolution=RES)
    m = src.meta.copy()
    m.update(crs=CRS_UTM, transform=tr, width=w, height=h, nodata=-9999.0, dtype="float32",
             compress="deflate")
    with rasterio.open(dem_utm, "w", **m) as dst:
        reproject(rasterio.band(src, 1), rasterio.band(dst, 1), src_nodata=src.nodata,
                  dst_nodata=-9999.0, resampling=Resampling.bilinear)
print(f"DEM UTM: {w} x {h} celdas de {RES:.0f} m")

# ------------------------------------------------------------------ 3-4. acondicionamiento, D8
grid = Grid.from_raster(str(dem_utm))
dem = grid.read_raster(str(dem_utm))
dirmap = (64, 128, 1, 2, 4, 8, 16, 32)          # N, NE, E, SE, S, SW, W, NW
filled = grid.resolve_flats(grid.fill_depressions(grid.fill_pits(dem)))
fdir = grid.flowdir(filled, dirmap=dirmap)
acc = grid.accumulation(fdir, dirmap=dirmap)
cell_km2 = RES * RES / 1e6
acc_km2 = np.asarray(acc) * cell_km2

# ------------------------------------------------------------------ 5. punto de cierre
gauge_utm = gpd.GeoSeries([Point(GAUGE_LON, GAUGE_LAT)], crs=4326).to_crs(CRS_UTM).iloc[0]
aff = grid.affine
col0, row0 = ~aff * (gauge_utm.x, gauge_utm.y)
col0, row0 = int(col0), int(row0)
rad = int(SNAP_RADIUS_M / RES)
rr, cc = np.mgrid[row0 - rad:row0 + rad + 1, col0 - rad:col0 + rad + 1]
dist = np.hypot(rr - row0, cc - col0) * RES
win = acc_km2[rr, cc]
# Celdas del cauce dentro del radio con área drenada dentro de ±15 % del área CAMELS
ok = (dist <= SNAP_RADIUS_M) & (np.abs(win / AREA_CAMELS - 1) <= 0.15)
if not ok.any():
    raise RuntimeError("No hay celda de cauce con área compatible cerca de la estación; revisar.")
i = np.argmin(np.where(ok, dist, np.inf))
r_out, c_out = rr.flat[i], cc.flat[i]
x_out, y_out = aff * (c_out + 0.5, r_out + 0.5)
print(f"Punto de cierre ajustado: desplazamiento {dist.flat[i]:.0f} m, "
      f"área drenada {acc_km2[r_out, c_out]:.1f} km2 (CAMELS {AREA_CAMELS:.1f} km2)")

# ------------------------------------------------------------------ 6. cuenca
catch = grid.catchment(x=c_out, y=r_out, fdir=fdir, dirmap=dirmap, xytype="index")
cmask = np.asarray(catch).astype(bool)
polys = [shape(g) for g, v in shapes(cmask.astype("uint8"), mask=cmask, transform=aff) if v == 1]
basin = max(polys, key=lambda p: p.area)
basin = basin.buffer(0)
basin_gdf = gpd.GeoDataFrame({"gauge_id": [GID], "nombre": [GAUGE_NAME]}, geometry=[basin],
                             crs=CRS_UTM)
basin_gdf.to_file(OUT / "cuenca_el_tambo_glo30.gpkg", layer="cuenca")
basin_gdf.to_crs(4326).to_file(OUT / "cuenca_el_tambo_glo30_wgs84.geojson", driver="GeoJSON")

# Red de drenaje y orden de Strahler (solo dentro de la cuenca)
thr_cells = STREAM_THR_KM2 / cell_km2
streams_mask = (acc > thr_cells) & catch
order = np.asarray(grid.stream_order(fdir, streams_mask, dirmap=dirmap))
branches = grid.extract_river_network(fdir, streams_mask, dirmap=dirmap)
lines, orders = [], []
for feat in branches["features"]:
    ln = LineString(feat["geometry"]["coordinates"])
    if ln.length == 0:
        continue
    pts = np.array(ln.coords)
    mid = pts[len(pts) // 2]
    c, r = ~aff * (mid[0], mid[1])
    lines.append(ln)
    orders.append(int(order[int(r), int(c)]))
rivers = gpd.GeoDataFrame({"strahler": orders}, geometry=lines, crs=CRS_UTM)
rivers["long_km"] = rivers.length / 1000
rivers.to_file(OUT / "red_drenaje_el_tambo.gpkg", layer="red_drenaje")

# Río principal: camino de flujo más largo hasta el cierre (distancia D8 ponderada)
OFF = {64: (-1, 0), 128: (-1, 1), 1: (0, 1), 2: (1, 1), 4: (1, 0), 8: (1, -1), 16: (0, -1),
       32: (-1, -1)}
fd = np.asarray(fdir)
w = Raster(np.where(np.isin(fd, [2, 8, 32, 128]), RES * np.sqrt(2), RES),
           viewfinder=fdir.viewfinder)
dto = np.asarray(grid.distance_to_outlet(x=c_out, y=r_out, fdir=fdir, dirmap=dirmap,
                                         xytype="index", weights=w))
dto = np.where(cmask, dto, np.nan)
r, c = np.unravel_index(np.nanargmax(dto), dto.shape)
path = []
while True:
    path.append(aff * (c + 0.5, r + 0.5))
    if (r, c) == (r_out, c_out):
        break
    dr, dc = OFF[int(fd[r, c])]
    r, c = r + dr, c + dc
    if not cmask[r, c]:
        break
main = LineString(path)
main_gdf = gpd.GeoDataFrame({"nombre": ["Río principal (Chicama)"]}, geometry=[main], crs=CRS_UTM)
main_gdf.to_file(OUT / "rio_principal_el_tambo.gpkg", layer="rio_principal")

outlet = gpd.GeoDataFrame(
    {"tipo": ["Estación SENAMHI", "Cierre ajustado (GLO-30)"],
     "gauge_id": [GID, GID]},
    geometry=[gauge_utm, Point(x_out, y_out)], crs=CRS_UTM)
outlet.to_file(OUT / "punto_cierre_el_tambo.gpkg", layer="punto_cierre")

# ------------------------------------------------------------------ 7. morfometría
z = np.asarray(dem)
zc = np.where(cmask, z, np.nan)
gy, gx = np.gradient(np.where(np.isfinite(zc), zc, np.nanmean(zc)), RES)
slope = np.degrees(np.arctan(np.hypot(gx, gy)))
slope = np.where(cmask, slope, np.nan)
A = basin.area / 1e6
P_raster = basin.length / 1000
# El perímetro del polígono raster tiene "escalones" de 30 m que lo inflan; se suaviza
# con cierre/apertura morfológica (±60 m) y simplificación (30 m) antes de calcular Kc.
P = basin.buffer(60).buffer(-60).simplify(RES).length / 1000
L = main.length / 1000
z_src = z[np.unravel_index(np.nanargmax(dto), dto.shape)]
z_out = z[r_out, c_out]
Ltot = rivers.length.sum() / 1000
camels = gpd.read_file(CAMELS_GPKG).query("gauge_id == @GID").to_crs(CRS_UTM).geometry.iloc[0]
iou = basin.intersection(camels).area / unary_union([basin, camels]).area
hyps = np.nanpercentile(zc, np.arange(0, 101, 5))

morfo = pd.Series({
    "Área A (km2)": A,
    "Área CAMELS-PE (km2)": AREA_CAMELS,
    "Diferencia de área vs CAMELS (%)": 100 * (A / AREA_CAMELS - 1),
    "IoU con polígono CAMELS-PE": iou,
    "Perímetro raster (km)": P_raster,
    "Perímetro suavizado P (km)": P,
    "Coef. de compacidad Kc = 0.282 P/sqrt(A)": 0.282 * P / np.sqrt(A),
    "Longitud río principal L (km)": L,
    "Factor de forma Kf = A/L^2": A / L**2,
    "Cota mínima (m s.n.m.)": np.nanmin(zc),
    "Cota máxima (m s.n.m.)": np.nanmax(zc),
    "Cota media (m s.n.m.)": np.nanmean(zc),
    "Cota mediana (m s.n.m.)": np.nanmedian(zc),
    "Pendiente media de la cuenca (°)": np.nanmean(slope),
    "Cota naciente río principal (m)": z_src,
    "Cota cierre (m)": z_out,
    "Pendiente media río principal (m/m)": (z_src - z_out) / (L * 1000),
    f"Longitud red de drenaje, umbral {STREAM_THR_KM2} km2 (km)": Ltot,
    "Densidad de drenaje Dd (km/km2)": Ltot / A,
    "Orden de Strahler máximo": rivers["strahler"].max(),
    "Coord. cierre ajustado UTM17S E (m)": x_out,
    "Coord. cierre ajustado UTM17S N (m)": y_out,
}).round(3)
morfo.to_csv(OUT / "parametros_morfometricos.csv", header=["valor"])
pd.DataFrame({"percentil_area_por_debajo": np.arange(0, 101, 5), "cota_m": hyps}).to_csv(
    OUT / "curva_hipsometrica.csv", index=False)
print("\n=== Parámetros morfométricos ===")
print(morfo.to_string())

# ------------------------------------------------------------------ 8. mapa
if not NE_FILE.exists():
    urllib.request.urlretrieve(NE_URL, NE_FILE)
countries = gpd.read_file(NE_FILE)

# El mapa se dibuja en coordenadas geográficas (EPSG:4326) para rotular la cuadrícula en
# grados; los cálculos hidrológicos y morfométricos se hicieron en UTM 17S.
basin_ll = basin_gdf.to_crs(4326)
rivers_ll = rivers.to_crs(4326)
main_ll = main_gdf.to_crs(4326)
out_ll = gpd.GeoSeries([Point(x_out, y_out)], crs=CRS_UTM).to_crs(4326).iloc[0]
bminx, bminy, bmaxx, bmaxy = basin_ll.total_bounds
PAD = 0.03                                       # margen del mapa (°)
xmin, xmax, ymin, ymax = bminx - PAD, bmaxx + PAD, bminy - PAD, bmaxy + PAD
lat0 = np.radians((ymin + ymax) / 2)
M_PER_DEG_X = 111320 * np.cos(lat0)              # m por grado de longitud a la latitud media
M_PER_DEG_Y = 110574                             # m por grado de latitud

# DEM geográfico original (sin re-muestreo) recortado al polígono de la cuenca
with rasterio.open(dem_geo) as src:
    win = from_bounds(xmin, ymin, xmax, ymax, transform=src.transform)
    win = win.round_offsets().round_lengths()
    zg = src.read(1, window=win).astype(float)
    tg = src.window_transform(win)
    if src.nodata is not None:
        zg[zg == src.nodata] = np.nan
inside_g = rasterize([(basin_ll.geometry.iloc[0], 1)], out_shape=zg.shape, transform=tg,
                     fill=0).astype(bool)
hs = LightSource(azdeg=315, altdeg=45).hillshade(
    np.nan_to_num(zg, nan=np.nanmean(zg)), vert_exag=2,
    dx=abs(tg.a) * M_PER_DEG_X, dy=abs(tg.e) * M_PER_DEG_Y)
ext_g = (tg.c, tg.c + tg.a * zg.shape[1], tg.f + tg.e * zg.shape[0], tg.f)

# Estaciones CAMELS-PE dentro del rectángulo del mapa (dentro de la cuenca y vecinas)
st = pd.read_csv(CAMELS_GPKG.parents[1] / "01_metadata" / "stations.csv")
st = gpd.GeoDataFrame(st, geometry=gpd.points_from_xy(st["gauge_lon"], st["gauge_lat"]), crs=4326)
st = st.cx[xmin:xmax, ymin:ymax].copy()
# Las estaciones de aforo están sobre el cauce, en el borde de su propia cuenca: se usa un
# margen de 500 m para decidir si están dentro de la cuenca de El Tambo.
basin_buf = basin_gdf.buffer(500).to_crs(4326).iloc[0]
st["en_cuenca"] = st.within(basin_buf)
st[["gauge_id", "gauge_name", "gauge_lat", "gauge_lon", "gauge_elev", "gauge_record_start",
    "gauge_record_end", "gauge_perc_obs", "en_cuenca"]].to_csv(OUT / "estaciones_camels_mapa.csv",
                                                              index=False)
print("\n=== Estaciones CAMELS-PE en el rectángulo del mapa ===")
print(st[["gauge_id", "gauge_name", "gauge_record_start", "gauge_record_end", "gauge_perc_obs",
          "en_cuenca"]].to_string(index=False))

# Red de estaciones SENAMHI: listado publicado en el mapa web de estaciones de SENAMHI
# (arreglo JSON embebido en la página). Se guarda una copia fechada para reproducibilidad.
SEN_DIR = ROOT / "El_Tambo" / "datos" / "estaciones"
SEN_DIR.mkdir(parents=True, exist_ok=True)
sen_files = sorted(SEN_DIR.glob("senamhi_estaciones_*.csv"))
if sen_files:
    sen_all = pd.read_csv(sen_files[-1], dtype={"cod": str, "cod_old": str})
else:
    html = urllib.request.urlopen(SENAMHI_URL, timeout=60).read().decode("utf-8", "ignore")
    arr = re.search(r"var\s+PruebaTest\s*=\s*(\[.*?\]);", html, re.S).group(1)
    arr = re.sub(r":\s*(-?)\.(\d)", r": \g<1>0.\2", arr)      # corrige números como "-.1172"
    sen_all = pd.DataFrame(json.loads(arr))
    sen_all.to_csv(SEN_DIR / f"senamhi_estaciones_{pd.Timestamp.today():%Y%m%d}.csv", index=False)
sen = gpd.GeoDataFrame(sen_all, geometry=gpd.points_from_xy(sen_all["lon"], sen_all["lat"]), crs=4326)
sen = sen.cx[xmin:xmax, ymin:ymax].copy()
sen["en_cuenca"] = sen.within(basin_buf)
sen.drop(columns="geometry").to_csv(OUT / "estaciones_senamhi_mapa.csv", index=False)
# Las estaciones hidrológicas que coinciden con aforos CAMELS-PE (< 300 m) ya se dibujan como CAMELS.
st_utm = st.to_crs(CRS_UTM)
sen_utm = sen.to_crs(CRS_UTM)
sen["es_camels"] = [st_utm.distance(p).min() < 300 for p in sen_utm.geometry]
# Une estaciones convencionales y automáticas del mismo nombre ubicadas a menos de 1 km.
sen["nom"] = sen["nom"].str.split().str.join(" ").str.title()
sen_plot = []
for nom, g in sen[~sen["es_camels"]].groupby("nom"):
    gu = g.to_crs(CRS_UTM)
    if gu.geometry.distance(gu.geometry.iloc[0]).max() < 1000:
        sen_plot.append(dict(nom=nom, cate=" + ".join(sorted(g["cate"])), tipo=g["ico"].iloc[0],
                             en_cuenca=bool(g["en_cuenca"].any()), geometry=g.geometry.iloc[0]))
    else:
        sen_plot += [dict(nom=nom, cate=r["cate"], tipo=r["ico"], en_cuenca=bool(r["en_cuenca"]),
                          geometry=r.geometry) for _, r in g.iterrows()]
sen_plot = gpd.GeoDataFrame(sen_plot, crs=4326)
print("\n=== Estaciones SENAMHI en el rectángulo del mapa (sin duplicar aforos CAMELS) ===")
print(sen_plot.drop(columns="geometry").to_string(index=False))

fig = plt.figure(figsize=(13.5, 11))
gs = GridSpec(3, 2, figure=fig, width_ratios=[3.6, 1], height_ratios=[1, 1, 1.15],
              wspace=0.06, hspace=0.12)
ax = fig.add_subplot(gs[:, 0])
im = ax.imshow(zg, cmap="terrain", extent=ext_g, alpha=0.9, interpolation="nearest",
               vmin=np.nanmin(zg), vmax=np.nanmax(zg))
ax.imshow(hs, cmap="gray", extent=ext_g, alpha=0.3, interpolation="nearest")
for o in sorted(rivers_ll["strahler"].unique()):
    rivers_ll[rivers_ll["strahler"] == o].plot(ax=ax, color="#1f5fbf",
                                               linewidth=0.3 + 0.45 * (o - 1))
main_ll.plot(ax=ax, color="#08306b", linewidth=2.6)
basin_ll.boundary.plot(ax=ax, color="black", linewidth=1.8)
ax.plot(out_ll.x, out_ll.y, "o", color="red", markeredgecolor="k", ms=11, zorder=7)
ax.annotate(f"{GAUGE_NAME} ({GID})\npunto de cierre", (out_ll.x, out_ll.y),
            xytext=(45, -45), textcoords="offset points", fontsize=8.5,
            arrowprops=dict(arrowstyle="->"), bbox=dict(boxstyle="round", fc="white", alpha=0.85))
for _, s in st[st["gauge_id"] != GID].iterrows():
    color = "#ffd92f" if s["en_cuenca"] else "#bdbdbd"
    ax.plot(s.geometry.x, s.geometry.y, "^", color=color, markeredgecolor="k", ms=10, zorder=7)
    ax.annotate(f"{s['gauge_name']} ({s['gauge_id']})\n{s['gauge_record_start'][:4]}–"
                f"{s['gauge_record_end'][:4]}, {s['gauge_perc_obs']:.0f} % obs.",
                (s.geometry.x, s.geometry.y), xytext=(8, 6), textcoords="offset points",
                fontsize=7.5, bbox=dict(boxstyle="round", fc="white", alpha=0.8, lw=0.3))
SEN_STYLE = {"M": dict(marker="s", color="#4daf4a"), "H": dict(marker="D", color="#00bfff")}
for _, s in sen_plot.iterrows():
    ax.plot(s.geometry.x, s.geometry.y, ls="", markeredgecolor="k", ms=7.5, zorder=7,
            **SEN_STYLE[s["tipo"]])
    below_ok = s.geometry.y > ymin + 0.03           # cerca del borde inferior: rótulo arriba
    ax.annotate(f"{s['nom']} ({s['cate']})", (s.geometry.x, s.geometry.y),
                xytext=(6, -10 if below_ok else 6),
                textcoords="offset points", fontsize=6.8,
                bbox=dict(boxstyle="round,pad=0.15", fc="white", alpha=0.75, lw=0))

# Fuentes de humedad (esquema; ver explicación y bibliografía en el informe)
ax.annotate("", xy=(0.80, 0.60), xytext=(0.995, 0.60), xycoords="axes fraction",
            arrowprops=dict(arrowstyle="simple,head_width=2.2,head_length=1.6,tail_width=1.0",
                            fc="#2166ac", ec="k", alpha=0.9), zorder=8)
ax.text(0.985, 0.635, "Humedad amazónica / atlántica\nflujo del E en niveles medios-altos\n"
        "(principal, dic–abr)", transform=ax.transAxes, ha="right", va="bottom", fontsize=8,
        color="#08306b", fontweight="bold", bbox=dict(fc="white", alpha=0.75, lw=0), zorder=8)
ax.annotate("", xy=(0.17, 0.78), xytext=(0.005, 0.78), xycoords="axes fraction",
            arrowprops=dict(arrowstyle="simple,head_width=1.6,head_length=1.2,tail_width=0.6",
                            fc="white", ec="#b2182b", ls="--", lw=1.5, alpha=0.9), zorder=8)
ax.text(0.01, 0.80, "Humedad del Pacífico\n(secundaria; relevante en El Niño)",
        transform=ax.transAxes, ha="left", va="bottom", fontsize=8, color="#b2182b",
        fontweight="bold", bbox=dict(fc="white", alpha=0.75, lw=0), zorder=8)
ax.set_xlim(xmin, xmax); ax.set_ylim(ymin, ymax)
ax.set_aspect(1 / np.cos(lat0))                  # evita distorsión E-O en grados
ax.xaxis.set_major_locator(MultipleLocator(0.1))
ax.yaxis.set_major_locator(MultipleLocator(0.1))
ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{abs(v):.2f}° {'O' if v < 0 else 'E'}"))
ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{abs(v):.2f}° {'S' if v < 0 else 'N'}"))
ax.set_xlabel("Longitud"); ax.set_ylabel("Latitud")
ax.tick_params(axis="y", labelrotation=90)
ax.grid(color="grey", alpha=0.35, linewidth=0.5, linestyle="--")
ax.add_artist(ScaleBar(M_PER_DEG_X, units="m", location="lower left", box_alpha=0.8,
                     rotation="horizontal-only"))
# Flecha de norte
ax.annotate("N", xy=(0.94, 0.95), xytext=(0.94, 0.86), xycoords="axes fraction",
            ha="center", va="center", fontsize=16, fontweight="bold",
            arrowprops=dict(facecolor="black", width=6, headwidth=16))
# Leyenda
handles = [Line2D([], [], color="black", lw=1.8, label="Límite de cuenca"),
           Line2D([], [], color="#08306b", lw=2.6, label="Río Chicama"),
           Line2D([], [], color="#1f5fbf", lw=1, label=f"Red de drenaje (A > {STREAM_THR_KM2:g} km², Strahler)"),
           Line2D([], [], color="red", marker="o", lw=0, markeredgecolor="k", ms=9,
                  label="Punto de cierre (El Tambo)"),
           Line2D([], [], color="#ffd92f", marker="^", lw=0, markeredgecolor="k", ms=9,
                  label="Estación CAMELS-PE dentro de la cuenca")]
if (~st["en_cuenca"]).any():
    handles.append(Line2D([], [], color="#bdbdbd", marker="^", lw=0, markeredgecolor="k", ms=9,
                          label="Estación CAMELS-PE vecina"))
handles += [Line2D([], [], ls="", markeredgecolor="k", ms=7.5, label="SENAMHI meteorológica",
                   **SEN_STYLE["M"]),
            Line2D([], [], ls="", markeredgecolor="k", ms=7.5, label="SENAMHI hidrológica",
                   **SEN_STYLE["H"])]
CATE_TXT = {"CO": "climática ordinaria", "CP": "climática principal", "PLU": "pluviométrica",
            "PE": "propósito específico", "EMA": "meteorológica automática",
            "HLM": "hidrológica limnimétrica", "HLG": "hidrológica limnigráfica",
            "EHMA": "hidrometeorológica automática", "EHA": "hidrológica automática",
            "EAMA": "agrometeorológica automática", "MAP": "meteorológica agrícola principal"}
cats = sorted({c for v in sen_plot["cate"] for c in v.split(" + ")})
cats_txt = "Categorías SENAMHI: " + "; ".join(f"{c} = {CATE_TXT.get(c, c)}" for c in cats)
ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.075), ncol=3,
          fontsize=8.5, frameon=False)
ax.set_title(f"Cuenca del río Chicama hasta la estación {GAUGE_NAME} ({GID})\n"
             f"Área = {A:,.0f} km²  ·  DEM Copernicus GLO-30 (30 m)", fontsize=12)
ax.text(0.5, -0.165, cats_txt + "\nCoordenadas geográficas WGS 84 (EPSG:4326)  ·  Delimitación "
        "y morfometría calculadas en WGS 84 / UTM 17S (EPSG:32717)\nFuente DEM: Copernicus "
        "GLO-30 DSM (ESA/Airbus)  ·  Red de drenaje derivada con D8 (pysheds)  ·  Estaciones: "
        "CAMELS-PE v1.0.1 y red SENAMHI (mapa web, consulta 2026-10)  ·  Límites: Natural Earth",
        transform=ax.transAxes, ha="center", va="top", fontsize=7.5)

# Mapas de ubicación: Sudamérica y Perú
peru = countries[countries["ADMIN"] == "Peru"]
PERU_BOX = (-82, -18.6, -68, 0.5)
ax_sa = fig.add_subplot(gs[0, 1])
countries[countries["CONTINENT"] == "South America"].plot(ax=ax_sa, color="#eeeeee",
                                                          edgecolor="#999999", linewidth=0.3)
peru.plot(ax=ax_sa, color="#d9c9a3", edgecolor="#555555", linewidth=0.5)
ax_sa.add_patch(Rectangle(PERU_BOX[:2], PERU_BOX[2] - PERU_BOX[0], PERU_BOX[3] - PERU_BOX[1],
                          fill=False, edgecolor="red", linewidth=1))
ax_sa.plot(out_ll.x, out_ll.y, "o", color="red", ms=4)
ax_sa.annotate("", xy=(-77, -7.5), xytext=(-52, -3), arrowprops=dict(
    arrowstyle="simple,head_width=0.9,head_length=0.8,tail_width=0.35", fc="#2166ac", ec="k"))
ax_sa.annotate("", xy=(-79.5, -7.6), xytext=(-83, -6.5), arrowprops=dict(
    arrowstyle="simple,head_width=0.6,head_length=0.6,tail_width=0.2", fc="white",
    ec="#b2182b", ls="--"))
ax_sa.set_xlim(-83, -33); ax_sa.set_ylim(-56, 14)
ax_sa.set_title("Ubicación en Sudamérica", fontsize=9)

ax_pe = fig.add_subplot(gs[1, 1])
countries.cx[-83:-66, -20:2].plot(ax=ax_pe, color="#eeeeee", edgecolor="#999999", linewidth=0.4)
peru.plot(ax=ax_pe, color="#d9c9a3", edgecolor="#555555", linewidth=0.6)
basin_ll.plot(ax=ax_pe, color="red")
ax_pe.add_patch(Rectangle((xmin, ymin), xmax - xmin, ymax - ymin, fill=False,
                          edgecolor="red", linewidth=0.8))
ax_pe.set_xlim(PERU_BOX[0], PERU_BOX[2]); ax_pe.set_ylim(PERU_BOX[1], PERU_BOX[3])
ax_pe.set_title("Ubicación en Perú", fontsize=9)
for a in (ax_sa, ax_pe):
    a.set_xticks([]); a.set_yticks([]); a.set_xlabel(""); a.set_ylabel("")
    a.set_facecolor("#cfe3f3"); a.set_aspect("equal")

# Barra de colores en la celda inferior derecha
cell = fig.add_subplot(gs[2, 1]); cell.axis("off")
cax = cell.inset_axes([0.30, 0.04, 0.12, 0.92])
cb = fig.colorbar(im, cax=cax)
cb.set_label("Elevación (m s.n.m.)\nCopernicus GLO-30", fontsize=9)
fig.savefig(OUT / "mapa_cuenca_el_tambo.png", dpi=200, bbox_inches="tight")
plt.close(fig)

# Comparación con CAMELS-PE
fig, ax = plt.subplots(figsize=(6, 6))
gpd.GeoSeries([camels], crs=CRS_UTM).boundary.plot(ax=ax, color="tab:orange", lw=1.5,
                                                   label="CAMELS-PE")
basin_gdf.boundary.plot(ax=ax, color="k", lw=1, label="GLO-30 (este trabajo)")
ax.legend(); ax.set_title(f"Comparación de delimitaciones · IoU = {iou:.3f}")
ax.ticklabel_format(style="plain"); ax.tick_params(axis="y", labelrotation=90)
fig.savefig(OUT / "comparacion_camels.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"\nSalidas en {OUT}")
