"""
Mapas temáticos de contexto físico de la cuenca El Tambo (río Chicama) y ficha de atributos.

Todos los mapas se presentan en coordenadas geográficas WGS 84 (EPSG:4326), con cuadrícula en
grados. Las capas se llevan a la malla del DEM GLO-30 geográfico (1", ~30 m); pendiente y
orientación se calculan con el espaciamiento de la malla convertido a metros a la latitud media.

  (a) Pendiente (GLO-30), escala continua
  (b) Orientación de laderas (GLO-30)
  (c) Cobertura del suelo (ESA WorldCover 2021, 10 m)
  (d) Unidades geológicas (INGEMMET, Geología 1:100 000 integrada)
  (e) Precipitación media anual del periodo común P_L–P_I, PISCOp v3.0 (SENAMHI, fuente de apoyo). CAMELS-PE
      (fuente principal de P_L) solo da series promediadas por cuenca, sin malla; el campo
      espacial se toma de PISCOp v3.0 diario, en sus celdas de 0.1° sin interpolar.
Además: curva hipsométrica y ficha de fracciones por clase (ficha_contexto_fisico.csv).

Requiere haber ejecutado 06_delimitacion_el_tambo_glo30.py, 07a_descarga_contexto_fisico.py y
14_piscop_v3_cuenca.py (recorte de PISCOp v3.0 y serie de cuenca).
"""
from pathlib import Path

import numpy as np
import pandas as pd
import geopandas as gpd
import xarray as xr
import rasterio
from rasterio.warp import reproject, Resampling
from rasterio.features import rasterize
from rasterio.windows import from_bounds
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm, LightSource
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

from utils_mapas import formato_grados, escala_y_norte, m_por_grado_lon, M_PER_DEG_Y
from periodo_comun import periodo_comun, fin_diario

ROOT = Path(__file__).resolve().parents[2] / "Tarea 1" / "El_Tambo"
DATOS = ROOT / "datos"
DEL = ROOT / "resultados" / "delimitacion"
OUT = ROOT / "resultados" / "contexto_fisico"
OUT.mkdir(parents=True, exist_ok=True)
CRS_TXT = "Coordenadas geográficas WGS 84 (EPSG:4326)"
PAD = 0.03  # margen del mapa alrededor de la cuenca (°)

# ------------------------------------------------------------------ vectores
basin = gpd.read_file(DEL / "cuenca_el_tambo_glo30.gpkg").to_crs(4326)
rivers = gpd.read_file(DEL / "red_drenaje_el_tambo.gpkg").to_crs(4326)
main = gpd.read_file(DEL / "rio_principal_el_tambo.gpkg").to_crs(4326)
outlet = gpd.read_file(DEL / "punto_cierre_el_tambo.gpkg").to_crs(4326).iloc[[1]]
A_km2 = basin.to_crs(32717).area.iloc[0] / 1e6
bx0, by0, bx1, by1 = basin.total_bounds
xmin, xmax, ymin, ymax = bx0 - PAD, bx1 + PAD, by0 - PAD, by1 + PAD
LAT0 = (ymin + ymax) / 2

# ------------------------------------------------------------------ malla geográfica (DEM)
with rasterio.open(DATOS / "dem" / "GLO30_El_Tambo_EPSG4326.tif") as src:
    win = from_bounds(xmin, ymin, xmax, ymax, transform=src.transform).round_offsets().round_lengths()
    dem = src.read(1, window=win).astype(float)
    T = src.window_transform(win)
    if src.nodata is not None:
        dem[dem == src.nodata] = np.nan
SHAPE = dem.shape
EXT = (T.c, T.c + T.a * SHAPE[1], T.f + T.e * SHAPE[0], T.f)
mask = rasterize([(basin.geometry.iloc[0], 1)], out_shape=SHAPE, transform=T, fill=0).astype(bool)
DX, DY = abs(T.a) * m_por_grado_lon(LAT0), abs(T.e) * M_PER_DEG_Y   # tamaño de celda (m)

# ------------------------------------------------------------------ derivados del DEM
gy, gx = np.gradient(np.where(np.isnan(dem), np.nanmean(dem), dem), DY, DX)
slope = np.degrees(np.arctan(np.hypot(gx, gy)))
# Orientación: dirección hacia la que mira la ladera (0 = N, 90 = E); las filas crecen al sur
aspect = (np.degrees(np.arctan2(-gx, gy)) + 360) % 360
asp_cls = (((aspect + 22.5) % 360) // 45).astype(int)          # 0..7 = N, NE, ..., NO
asp_cls[slope < 2] = 8                                         # plano
HS = LightSource(315, 45).hillshade(np.nan_to_num(dem, nan=np.nanmean(dem)), vert_exag=2,
                                    dx=DX, dy=DY)

# ------------------------------------------------------------------ cobertura (WorldCover)
WC = {10: ("Bosque", "#006400"), 20: ("Arbustal", "#ffbb22"), 30: ("Pastizal / pajonal", "#ffff4c"),
      40: ("Agricultura", "#f096ff"), 50: ("Área urbana", "#fa0000"),
      60: ("Suelo desnudo / vegetación escasa", "#b4b4b4"), 70: ("Nieve y hielo", "#f0f0f0"),
      80: ("Cuerpos de agua", "#0064c8"), 90: ("Humedal herbáceo (bofedal)", "#0096a0"),
      95: ("Manglar", "#00cf75"), 100: ("Musgo y liquen", "#fae6a0")}
lc = np.full(SHAPE, 0, dtype="uint8")
with rasterio.open(DATOS / "cobertura" / "ESA_WorldCover_10m_2021_v200_El_Tambo.tif") as src:
    reproject(rasterio.band(src, 1), lc, dst_transform=T, dst_crs="EPSG:4326",
              src_nodata=src.nodata, dst_nodata=0, resampling=Resampling.mode)

# ------------------------------------------------------------------ geología (INGEMMET 1:100k)
geo = gpd.read_file(DATOS / "geologia" / "ingemmet_geologia_100k.geojson").to_crs(4326)
geo["NAME"] = geo["NAME"].fillna("").str.strip()
geo.loc[geo["UNIDAD"].str.contains("Laguna", na=False), "NAME"] = "Laguna"
geo_in = gpd.overlay(geo, basin[["geometry"]], how="intersection").to_crs(32717)
units = (geo_in.assign(km2=geo_in.area / 1e6)
         .groupby("NAME").agg(UNIDAD=("UNIDAD", "first"), DESCRIP=("DESCRIP", "first"),
                              km2=("km2", "sum"))
         .sort_values("km2", ascending=False))
units["pct"] = 100 * units["km2"] / A_km2

# Agrupación hidrogeológica (propuesta a contrastar con INGEMMET, Boletines Serie H, y ANA):
# se agrupan unidades por tipo de porosidad/permeabilidad esperada, no por edad.
GEO_CLASSES = {
    1: ("Depósitos cuaternarios no consolidados (aluviales, fluviales, glaciares, lacustres)"
        " – acuífero poroso", "#ffd400"),
    2: ("Areniscas cuarzosas y cuarcitas (Fm. Chimú, Fm. Farrat) – acuífero fisurado "
        "sedimentario", "#ff6d00"),
    3: ("Secuencias de lutitas, limolitas y areniscas (Fm. Chicama, Santa, Carhuaz, Inca) y "
        "conglomerados (Fm. Huaylas) – acuitardo sedimentario", "#2ca02c"),
    4: ("Calizas y margas (Fm. Chúlec, Pariatambo, Pulluicana, Quilquiñán, Cajamarca) – "
        "acuífero kárstico-fisurado / acuitardo", "#1f78ff"),
    5: ("Rocas volcánicas (Gr. Calipuy, Fm. Llama, San Pablo, Chilete, Gr. Casma) – "
        "acuitardo / acuífero fisurado volcánico", "#e7298a"),
    6: ("Rocas intrusivas (granodiorita, granito, diorita) – acuífugo", "#d50000"),
    7: ("Lagunas", "#0064c8"),
}


def hidro_grupo(name, unidad):
    n, u = str(name), str(unidad).lower()
    if n == "Laguna":
        return 7
    if n.startswith("Q"):
        return 1
    if any(k in u for k in ["granodiorita", "granito", "diorita", "batolito"]):
        return 6
    if any(k in u for k in ["calipuy", "llama", "san pablo", "chilete", "andesita", "dacita",
                            "casma"]):
        return 5
    if any(k in u for k in ["chúlec", "pariatambo", "pulluicana", "quilquiñán", "cajamarca"]):
        return 4
    if any(k in u for k in ["chimú", "farrat"]):
        return 2
    return 3   # Chicama, Santa, Carhuaz, Inca, Santa-Carhuaz, Huaylas y otras clásticas


units["grupo_hidrogeologico"] = [hidro_grupo(n, u) for n, u in zip(units.index, units["UNIDAD"])]
units.round(3).to_csv(OUT / "unidades_geologicas_cuenca.csv")
print("\nUnidades INGEMMET y grupo hidrogeológico asignado:")
print(units[["UNIDAD", "pct", "grupo_hidrogeologico"]].round(2).to_string())
geo["cod"] = [hidro_grupo(n, u) for n, u in zip(geo["NAME"], geo["UNIDAD"])]
geo_r = rasterize([(g, v) for g, v in zip(geo.geometry, geo["cod"])], out_shape=SHAPE,
                  transform=T, fill=0).astype("uint8")

# ------------------------------------------------------------------ precipitación (PISCOp v3.0)
# Periodo común P_L–P_I. Como el último año está incompleto, la media anual se calcula como la suma de
# las 12 medias mensuales del periodo (usa todos los meses sin sesgar hacia los meses del año parcial).
INI, FIN = periodo_comun()
P_PER = f"{INI:%Y-%m} a {FIN:%Y-%m}"


def media_anual_clim(x):
    """Σ de las medias de cada mes calendario (x mensual: Series o DataArray con eje 'time')."""
    if isinstance(x, pd.Series):
        return x.groupby(x.index.month).mean().sum()
    return x.groupby("time.month").mean().sum("month")


pv3 = xr.open_dataset(DATOS / "piscop_senamhi" / "v3p0" / "PISCOp_v3_d_El_Tambo_recorte.nc")["precipitation"]
pv3 = pv3.sel(time=slice(INI, fin_diario(FIN)))
assert pv3.notnull().all(), "PISCOp v3.0: faltantes en el recorte"
P_CELDA = media_anual_clim(pv3.resample(time="MS").sum())          # mm/año por celda
# Cada píxel del DEM toma el valor de la celda de 0.1° que lo contiene (sin interpolar)
px_lon = T.c + T.a * (np.arange(SHAPE[1]) + 0.5)
px_lat = T.f + T.e * (np.arange(SHAPE[0]) + 0.5)
celda = P_CELDA.sel(lon=xr.DataArray(px_lon, dims="x"), lat=xr.DataArray(px_lat, dims="y"),
                    method="nearest", tolerance=0.05 + 1e-6)
precip = np.where(mask, celda.transpose("y", "x").values, np.nan)
SMES = ROOT / "resultados" / "series_mensuales"
P_MEDIA = media_anual_clim(pd.read_csv(SMES / "piscop_v3_cuenca_mensual.csv", parse_dates=["fecha"],
                                       index_col="fecha")["P_piscop_v3_mm_mes"].loc[INI:FIN])   # ponderado por área
P_CAMELS = media_anual_clim(pd.read_csv(SMES / "series_mensuales_camels_pe.csv", parse_dates=["fecha"],
                                        index_col="fecha")["P_L_mm_mes"].loc[INI:FIN])
print(f"\nP media anual {P_PER} (Σ medias mensuales): PISCOp v3.0 cuenca {P_MEDIA:.0f} mm/año (celdas "
      f"{np.nanmin(precip):.0f}–{np.nanmax(precip):.0f}) | CAMELS-PE prec {P_CAMELS:.0f} mm/año")

# ------------------------------------------------------------------ utilidades de mapa


def base(ax, title):
    ax.imshow(HS, cmap="gray", extent=EXT, alpha=0.5)
    ax.set_xlim(xmin, xmax); ax.set_ylim(ymin, ymax)
    ax.set_title(title, fontsize=11, loc="left", fontweight="bold")
    formato_grados(ax, LAT0)
    ax.tick_params(labelsize=7.5)
    ax.xaxis.label.set_size(8.5); ax.yaxis.label.set_size(8.5)
    ax.text(0.01, 0.99, "WGS 84 (EPSG:4326)", transform=ax.transAxes, ha="left", va="top",
            fontsize=7, bbox=dict(fc="white", alpha=0.8, lw=0), zorder=9)


def overlay(ax):
    rivers[rivers["strahler"] >= 3].plot(ax=ax, color="#1f5fbf", linewidth=0.5, alpha=0.8)
    main.plot(ax=ax, color="#08306b", linewidth=1.6)
    basin.boundary.plot(ax=ax, color="black", linewidth=1.3)
    outlet.plot(ax=ax, color="red", edgecolor="k", markersize=45, zorder=6)
    escala_y_norte(ax, LAT0)


def cat_layer(ax, arr, classes, title, legend_title):
    """Capa categórica: classes = {código: (etiqueta, color)}; leyenda con % de la cuenca."""
    base(ax, title)
    codes = sorted(classes)
    idx = np.full(arr.shape, np.nan)
    for i, c in enumerate(codes):
        idx[arr == c] = i
    cmap = ListedColormap([classes[c][1] for c in codes])
    ax.imshow(np.where(mask, idx, np.nan), cmap=cmap, extent=EXT, alpha=0.9,
              norm=BoundaryNorm(np.arange(len(codes) + 1) - 0.5, len(codes)),
              interpolation="nearest")
    overlay(ax)
    pct = {c: 100 * np.sum(arr[mask] == c) / mask.sum() for c in codes}
    handles = [Patch(facecolor=classes[c][1], edgecolor="grey", label=f"{classes[c][0]} ({pct[c]:.1f} %)")
               for c in codes if pct[c] >= 0.05]
    ax.legend(handles=handles, title=legend_title, fontsize=7, title_fontsize=8,
              loc="upper left", bbox_to_anchor=(0, -0.08), frameon=False)
    return pct


def slope_layer(ax):
    base(ax, "(a) Pendiente – Copernicus GLO-30")
    im = ax.imshow(np.where(mask, slope, np.nan), cmap="turbo", vmin=0, vmax=60,
                   extent=EXT, alpha=0.9, interpolation="nearest")
    overlay(ax)
    cb = plt.colorbar(im, ax=ax, orientation="horizontal", fraction=0.045, pad=0.08,
                      extend="max")
    s = slope[mask]
    cb.set_label(f"Pendiente (°)  ·  media {s.mean():.1f}°, mediana {np.median(s):.1f}°, "
                 f"P90 {np.percentile(s, 90):.1f}°", fontsize=8)
    cb.ax.tick_params(labelsize=7.5)
    return {}


def precip_layer(ax):
    base(ax, f"(e) Precipitación media anual {P_PER} – PISCOp v3.0 (SENAMHI)")
    im = ax.imshow(precip, cmap="Blues", vmin=0, vmax=np.ceil(np.nanmax(precip) / 100) * 100,
                   extent=EXT, alpha=0.85, interpolation="nearest")
    overlay(ax)
    cb = plt.colorbar(im, ax=ax, orientation="horizontal", fraction=0.045, pad=0.08)
    cb.set_label(f"P (mm/año) en celdas de 0.1° · promedio de cuenca ponderado por área "
                 f"{P_MEDIA:.0f} mm/año", fontsize=8)
    cb.ax.tick_params(labelsize=7.5)
    ax.text(0.01, -0.2, "Fuente de apoyo: PISCOp v3.0 diario (doi:10.6084/m9.figshare.32411886.v1),\n"
            "producto interpolado de pluviómetros y satélite; celdas sin interpolar.\n"
            f"La P$_L$ del análisis es CAMELS-PE ({P_CAMELS:.0f} mm/año), sin malla espacial.\n"
            "Periodo común P$_L$–P$_I$; media anual = Σ de las 12 medias mensuales.",
            transform=ax.transAxes, fontsize=7.5, va="top")
    return {}


ASP = {i: (n, c) for i, (n, c) in enumerate(zip(
    ["N", "NE", "E", "SE", "S", "SO", "O", "NO", "Plano (< 2°)"],
    ["#3288bd", "#66c2a5", "#abdda4", "#e6f598", "#fee08b", "#fdae61", "#f46d43", "#d53e4f",
     "#dddddd"]))}
themes = {
    "a_pendiente": slope_layer,
    "b_orientacion": lambda ax: cat_layer(ax, asp_cls, ASP, "(b) Orientación de laderas – GLO-30",
                                          "Ladera orientada hacia (% del área)"),
    "c_cobertura": lambda ax: cat_layer(ax, lc, {k: v for k, v in WC.items() if (lc[mask] == k).any()},
                                        "(c) Cobertura – ESA WorldCover 2021 (10 m)",
                                        "Cobertura (% del área)"),
    "d_geologia": lambda ax: cat_layer(ax, geo_r, GEO_CLASSES,
                                       "(d) Unidades hidrogeológicas – agrupación de "
                                       "INGEMMET 1:100 000",
                                       "Grupo hidrogeológico (% del área)"),
    "e_precipitacion": precip_layer,
}
COMMON = [Line2D([], [], color="black", lw=1.3, label="Límite de cuenca"),
          Line2D([], [], color="#08306b", lw=1.6, label="Río Chicama"),
          Line2D([], [], color="#1f5fbf", lw=0.6, label="Red de drenaje (Strahler ≥ 3)"),
          Line2D([], [], color="red", marker="o", lw=0, markeredgecolor="k",
                 label="Punto de cierre El Tambo")]

# Elimina mapas de versiones anteriores que ya no forman parte del producto
for old in ["mapa_a_elevacion.png", "mapa_b_pendiente.png", "mapa_c_orientacion.png",
            "mapa_d_cobertura.png", "mapa_e_geologia.png", "mapa_f_suelos.png"]:
    (OUT / old).unlink(missing_ok=True)

pcts = {}
for name, fn in themes.items():
    fig, ax = plt.subplots(figsize=(8.5, 10.5))
    pcts[name] = fn(ax)
    fig.legend(handles=COMMON, loc="lower right", fontsize=7, frameon=False,
               bbox_to_anchor=(0.98, 0.0))
    fig.savefig(OUT / f"mapa_{name}.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

# Panel compuesto 3 x 2 (cinco mapas; la sexta casilla queda vacía)
fig, axs = plt.subplots(3, 2, figsize=(17, 37))
fig.subplots_adjust(left=0.05, right=0.98, top=0.955, bottom=0.05, hspace=0.36, wspace=0.12)
for ax, fn in zip(axs.flat, themes.values()):
    fn(ax)
axs.flat[-1].axis("off")
fig.legend(handles=COMMON, loc="upper center", ncol=4, fontsize=10, frameon=False,
           bbox_to_anchor=(0.5, 0.983))
fig.suptitle("Contexto físico de la cuenca del río Chicama hasta la estación El Tambo "
             f"(CAMELS-PE PE_47E0A2A8, {A_km2:,.0f} km²)", fontsize=15, y=0.995)
fig.text(0.5, 0.975, CRS_TXT + "  ·  Límite, red y río Chicama derivados de Copernicus GLO-30 (D8)",
         ha="center", va="top", fontsize=9)
fig.savefig(OUT / "panel_contexto_fisico.png", dpi=150)
plt.close(fig)

# ------------------------------------------------------------------ curva hipsométrica
z = np.sort(dem[mask])[::-1]
frac = np.arange(1, z.size + 1) / z.size
bands = np.arange(500, 4500, 500)
band_pct = np.histogram(dem[mask], bins=bands)[0] / mask.sum() * 100
fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.5))
a1.plot(frac * 100, z, color="k")
a1.set_xlabel("Área por encima de la cota (%)"); a1.set_ylabel("Elevación (m s.n.m.)")
a1.set_title("Curva hipsométrica"); a1.grid(alpha=0.3)
hi = (np.nanmean(dem[mask]) - z.min()) / (z.max() - z.min())
a1.text(0.55, 0.85, f"Integral hipsométrica = {hi:.2f}", transform=a1.transAxes)
a2.barh(bands[:-1] + 250, band_pct, height=450, color="#7a9e7e")
a2.set_xlabel("% del área de la cuenca"); a2.set_ylabel("Banda altitudinal (m s.n.m.)")
a2.set_title("Distribución del área por bandas de 500 m"); a2.grid(alpha=0.3)
fig.tight_layout(); fig.savefig(OUT / "curva_hipsometrica.png", dpi=150); plt.close(fig)

# ------------------------------------------------------------------ ficha
s = slope[mask]
rows = [("elevacion", f"Banda {b}-{b+500} m", p) for b, p in zip(bands[:-1], band_pct)]
rows += [("pendiente", f"{a}-{b}°", 100 * np.mean((s >= a) & (s < b)))
         for a, b in [(0, 5), (5, 15), (15, 30), (30, 45), (45, 90)]]
for key, cls in [("b_orientacion", ASP), ("c_cobertura", WC), ("d_geologia", GEO_CLASSES)]:
    rows += [(key[2:], cls[c][0], p) for c, p in pcts[key].items()]
p_in = precip[mask]
rows += [("precipitacion_piscop_v3", f"{a}-{b} mm/año", 100 * np.mean((p_in >= a) & (p_in < b)))
         for a, b in [(0, 500), (500, 750), (750, 1000), (1000, 1500)]]
ficha = pd.DataFrame(rows, columns=["tema", "clase", "pct_area"]).round(2)
ficha.to_csv(OUT / "ficha_contexto_fisico.csv", index=False)
pd.set_option("display.width", 220, "display.max_colwidth", 110)
print(f"Área {A_km2:.1f} km2 | integral hipsométrica {hi:.2f} | pendiente media {s.mean():.1f}°")
print(ficha.to_string(index=False))
print(f"Salidas en {OUT}")
