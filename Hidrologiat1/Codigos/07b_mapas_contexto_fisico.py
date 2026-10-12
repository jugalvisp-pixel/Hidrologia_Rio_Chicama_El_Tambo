"""
Mapas temáticos de contexto físico de la cuenca El Tambo (río Chicama) y ficha de atributos.

Todos los mapas se presentan en coordenadas geográficas WGS 84 (EPSG:4326), con cuadrícula en
grados. Las capas se llevan a la malla del DEM GLO-30 geográfico (1", ~30 m); pendiente y
orientación se calculan con el espaciamiento de la malla convertido a metros a la latitud media.

  (a) Orientación de laderas (GLO-30)
  (b) Cobertura del suelo (ESA WorldCover 2021, 10 m)
  (c) Unidades geológicas (INGEMMET, Geología 1:100 000 integrada)
  (d) Precipitación media anual del periodo común P_L–P_I, CHIRPS v2.0 (0.05°, fuente de apoyo), interpolada
      de forma bilineal a la malla del DEM para una representación continua, con isoyetas cada 100 mm.
      La interpolación suaviza la presentación pero no agrega detalle: la información sigue siendo de 0.05°.
      CAMELS-PE (fuente principal de P_L) solo da series promediadas por cuenca, sin malla.
La pendiente se sigue calculando (orientación de zonas planas y ficha), pero ya no se dibuja (decisión del
grupo, 2026-10-10). Además: curva hipsométrica y ficha de fracciones por clase (ficha_contexto_fisico.csv).

Requiere haber ejecutado 06_delimitacion_el_tambo_glo30.py, 07a_descarga_contexto_fisico.py,
09_datos_camels_pe_el_tambo.py y 26_chirps_cuenca.py (recorte de CHIRPS y serie de cuenca).
"""
from scipy.interpolate import RegularGridInterpolator
from pathlib import Path
import textwrap

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
import matplotlib.patheffects as pe

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
        " – acuífero poroso", "#2b7bba"),
    2: ("Areniscas cuarzosas y cuarcitas (Fm. Chimú, Fm. Farrat) – acuífero fisurado "
        "sedimentario", "#2e9e44"),
    3: ("Secuencias de lutitas, limolitas y areniscas (Fm. Chicama, Santa, Carhuaz, Inca) y "
        "conglomerados (Fm. Huaylas) – acuitardo sedimentario", "#e6cfa3"),
    4: ("Calizas y margas (Fm. Chúlec, Pariatambo, Pulluicana, Quilquiñán, Cajamarca) – "
        "acuífero kárstico-fisurado / acuitardo", "#8fd18a"),
    5: ("Rocas volcánicas (Gr. Calipuy, Fm. Llama, San Pablo, Chilete, Gr. Casma) – "
        "acuitardo / acuífero fisurado volcánico", "#c49a6c"),
    6: ("Rocas intrusivas (granodiorita, granito, diorita) – acuífugo", "#7a5c48"),
    7: ("Lagunas", "#9fd3f0"),
}
# Colores según la convención de los mapas hidrogeológicos (IAH/UNESCO, usada por INGEMMET): azul =
# acuífero poroso; verdes = acuíferos fisurados y kársticos; marrones = acuitardos y acuífugos.
# Orden de la leyenda: de mayor a menor capacidad de almacenar y transmitir agua.
ORDEN_HIDRO = [1, 2, 4, 5, 3, 6, 7]


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
geo_in["grupo"] = [hidro_grupo(n, u) for n, u in zip(geo_in["NAME"], geo_in["UNIDAD"])]
units.round(3).to_csv(OUT / "unidades_geologicas_cuenca.csv")
print("\nUnidades INGEMMET y grupo hidrogeológico asignado:")
print(units[["UNIDAD", "pct", "grupo_hidrogeologico"]].round(2).to_string())
geo["cod"] = [hidro_grupo(n, u) for n, u in zip(geo["NAME"], geo["UNIDAD"])]
geo_r = rasterize([(g, v) for g, v in zip(geo.geometry, geo["cod"])], out_shape=SHAPE,
                  transform=T, fill=0).astype("uint8")

# Colores oficiales de INGEMMET: simbología de la capa Geología 1:100 000 integrada de GEOCATMIN
# (renderer de valores únicos por el código NAME). Se guarda una copia para reproducibilidad.
SIMB = DATOS / "geologia" / "ingemmet_simbologia_geologia100k.json"
if not SIMB.exists():
    import json, urllib.request
    url = ("https://geocatmin.ingemmet.gob.pe/arcgis/rest/services/SERV_GEOLOGIA_100K_INTEGRADA/MapServer/6"
           "?f=json")
    capa = json.loads(urllib.request.urlopen(url, timeout=60).read().decode("utf-8"))
    SIMB.write_text(json.dumps({"fuente": url, "consulta": pd.Timestamp.today().strftime("%Y-%m-%d"),
                                "renderer": capa["drawingInfo"]["renderer"]}, ensure_ascii=False), encoding="utf-8")
import json
_rend = json.loads(SIMB.read_text(encoding="utf-8"))["renderer"]
COLOR_INGEMMET = {v["value"].strip(): tuple(c / 255 for c in v["symbol"]["color"][:3])
                  for v in _rend["uniqueValueInfos"]}
COLOR_INGEMMET["Laguna"] = COLOR_INGEMMET.get("", (0.81, 1.0, 1.0))   # en GEOCATMIN las lagunas no tienen código
UID = {n: i + 1 for i, n in enumerate(sorted(geo["NAME"].unique()))}
geo_u = rasterize([(g, UID[n]) for g, n in zip(geo.geometry, geo["NAME"])], out_shape=SHAPE,
                  transform=T, fill=0).astype("int32")
ROTULO_KM2 = 20                                     # área mínima de un polígono para rotular su código
sin_color = [n for n in units.index if n not in COLOR_INGEMMET]
assert not sin_color, f"Unidades sin color en la simbología de INGEMMET: {sin_color}"

# ------------------------------------------------------------------ precipitación (CHIRPS v2.0)
# Periodo común P_L–P_I en años completos. La media anual se calcula como la suma de las 12 medias mensuales
# del periodo (con años completos equivale a la media de los totales anuales).
INI, FIN = periodo_comun()
P_PER = f"{INI:%Y-%m} a {FIN:%Y-%m}"


def media_anual_clim(x):
    """Σ de las medias de cada mes calendario (x mensual: Series o DataArray con eje 'time')."""
    if isinstance(x, pd.Series):
        return x.groupby(x.index.month).mean().sum()
    return x.groupby("time.month").mean().sum("month")


chirps = xr.open_dataset(DATOS / "chirps" / "CHIRPS_v2p0_mensual_El_Tambo.nc")["precipitation"]
chirps = chirps.sel(time=slice(INI, FIN))
assert chirps.notnull().all(), "CHIRPS: faltantes en el recorte"
P_CELDA = media_anual_clim(chirps)                                  # mm/año por celda de 0.05°
# Interpolación bilineal entre centros de celda a cada píxel del DEM (representación continua)
px_lon = T.c + T.a * (np.arange(SHAPE[1]) + 0.5)
px_lat = T.f + T.e * (np.arange(SHAPE[0]) + 0.5)
interp = RegularGridInterpolator((P_CELDA.lat.values, P_CELDA.lon.values), P_CELDA.values, method="linear")
LAT2, LON2 = np.meshgrid(px_lat, px_lon, indexing="ij")
P_CONT = interp(np.column_stack([LAT2.ravel(), LON2.ravel()])).reshape(SHAPE)
precip = np.where(mask, P_CONT, np.nan)
SMES = ROOT / "resultados" / "series_mensuales"
P_MEDIA = media_anual_clim(pd.read_csv(SMES / "chirps_cuenca_mensual.csv", parse_dates=["fecha"],
                                       index_col="fecha")["P_chirps_mm_mes"].loc[INI:FIN])   # ponderado por área
P_CAMELS = media_anual_clim(pd.read_csv(SMES / "series_mensuales_camels_pe.csv", parse_dates=["fecha"],
                                        index_col="fecha")["P_L_mm_mes"].loc[INI:FIN])
print(f"\nP media anual {P_PER} (Σ medias mensuales): CHIRPS cuenca {P_MEDIA:.0f} mm/año (píxeles "
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


def contornos(ax, campo):
    """Límites de los polígonos disueltos por `campo` (±25 m cierra las rendijas entre hojas 1:100 000)."""
    pol = geo_in[[campo, "geometry"]].assign(geometry=geo_in.buffer(25)).dissolve(campo).reset_index()
    pol = pol.assign(geometry=pol.buffer(-25)).explode(index_parts=False)
    pol = pol[pol.area > 0]
    pol.to_crs(4326).boundary.plot(ax=ax, color="#3a3a3a", linewidth=0.35, alpha=0.8, zorder=3)
    return pol


def hidro_layer(ax):
    """Unidades de INGEMMET unificadas por comportamiento hídrico (grupos hidrogeológicos)."""
    base(ax, "(c) Unidades hidrogeológicas – geología INGEMMET 1:100 000 agrupada")
    idx = np.full(SHAPE, np.nan)
    for i, c in enumerate(ORDEN_HIDRO):
        idx[geo_r == c] = i
    ax.imshow(np.where(mask, idx, np.nan), cmap=ListedColormap([GEO_CLASSES[c][1] for c in ORDEN_HIDRO]),
              norm=BoundaryNorm(np.arange(len(ORDEN_HIDRO) + 1) - 0.5, len(ORDEN_HIDRO)), extent=EXT,
              alpha=0.95, interpolation="nearest")
    contornos(ax, "grupo")
    overlay(ax)
    pct_g = {c: 100 * np.sum(geo_r[mask] == c) / mask.sum() for c in sorted(GEO_CLASSES)}
    miembros = {c: ", ".join(units.index[units["grupo_hidrogeologico"] == c].str.replace("Laguna", "—"))
                for c in GEO_CLASSES}
    handles = [Patch(facecolor=GEO_CLASSES[c][1], edgecolor="grey",
                     label="\n".join(textwrap.wrap(f"{GEO_CLASSES[c][0]} ({pct_g[c]:.1f} %)", 95))
                     + f"\n      unidades: {miembros[c]}")
               for c in ORDEN_HIDRO if pct_g[c] >= 0.005]
    ax.legend(handles=handles, title="Comportamiento hidrogeológico (% del área; de más a menos permeable)",
              fontsize=7, title_fontsize=8, loc="upper left", bbox_to_anchor=(0, -0.08), frameon=False,
              labelspacing=0.9)
    ax.text(0.0, -0.47, "Agrupación propia por porosidad y permeabilidad esperadas de cada unidad; a contrastar con "
            "el mapa hidrogeológico de INGEMMET y la ANA.\nColores según la convención IAH/UNESCO de los mapas "
            "hidrogeológicos. Unidades individuales: mapa_geologia_unidades_ingemmet.png.",
            transform=ax.transAxes, fontsize=7, va="top", color="#444")
    return pct_g


def geo_layer(ax):
    """Unidades geológicas con los colores oficiales de INGEMMET; el grupo hidrogeológico va en la leyenda."""
    base(ax, "Geología – INGEMMET 1:100 000 (colores oficiales) y grupo hidrogeológico")
    orden = list(units.index)                                          # de mayor a menor área
    idx = np.full(SHAPE, np.nan)
    for i, n in enumerate(orden):
        idx[geo_u == UID[n]] = i
    ax.imshow(np.where(mask, idx, np.nan), cmap=ListedColormap([COLOR_INGEMMET[n] for n in orden]),
              norm=BoundaryNorm(np.arange(len(orden) + 1) - 0.5, len(orden)), extent=EXT, alpha=0.95,
              interpolation="nearest")
    # Como en los mapas de INGEMMET: límites de unidad y código rotulado (unidades de la misma edad
    # comparten casi el mismo color, p. ej. las formaciones del Cretácico inferior)
    # Se disuelve por unidad para que no aparezcan los bordes entre hojas 1:100 000 (no son contactos)
    pol = contornos(ax, "NAME")
    grandes = pol[pol.area / 1e6 >= ROTULO_KM2]
    for _, r in grandes.to_crs(4326).iterrows():
        p = r.geometry.representative_point()
        ax.text(p.x, p.y, r["NAME"] if r["NAME"] != "Laguna" else "", fontsize=6.5, ha="center", va="center",
                zorder=4, path_effects=[pe.withStroke(linewidth=2, foreground="white")])
    overlay(ax)
    h_u = [Patch(facecolor=COLOR_INGEMMET[n], edgecolor="grey",
                 label=f"{n if n != 'Laguna' else '—'}  {r.UNIDAD} ({r.pct:.1f} %) · G{r.grupo_hidrogeologico}")
           for n, r in units.iterrows()]
    leg = ax.legend(handles=h_u, title="Unidad INGEMMET (% del área) · grupo hidrogeológico", ncol=2,
                    fontsize=6.5, title_fontsize=8, loc="upper left", bbox_to_anchor=(0, -0.08), frameon=False,
                    columnspacing=1.0, handlelength=1.4)
    ax.add_artist(leg)
    pct_g = {c: 100 * np.sum(geo_r[mask] == c) / mask.sum() for c in sorted(GEO_CLASSES)}
    h_g = [Line2D([], [], lw=0, label=f"G{c}: {GEO_CLASSES[c][0]} ({pct_g[c]:.1f} %)")
           for c in sorted(GEO_CLASSES) if pct_g[c] >= 0.005]
    ax.legend(handles=h_g, title="Grupos hidrogeológicos (agrupación propia, a contrastar con INGEMMET y ANA)",
              fontsize=6.5, title_fontsize=8, loc="upper left", bbox_to_anchor=(0, -0.36), frameon=False,
              handlelength=0, handletextpad=0)
    return pct_g


def precip_layer(ax):
    base(ax, f"(d) Precipitación media anual {P_PER} – CHIRPS v2.0")
    vmin = np.floor(np.nanmin(precip) / 100) * 100
    vmax = np.ceil(np.nanmax(precip) / 100) * 100
    im = ax.imshow(precip, cmap="Blues", vmin=vmin, vmax=vmax, extent=EXT, alpha=0.85, interpolation="bilinear")
    cs = ax.contour(px_lon, px_lat, precip, levels=np.arange(vmin, vmax + 1, 100), colors="#08306b",
                    linewidths=0.6, alpha=0.8)
    ax.clabel(cs, fmt="%d", fontsize=7, inline=True)
    overlay(ax)
    # Barra en un eje insertado: el mapa conserva el mismo tamaño que los demás del panel
    cb = plt.colorbar(im, cax=ax.inset_axes([0.0, -0.13, 1.0, 0.025]), orientation="horizontal")
    cb.set_label(f"P (mm/año) · isoyetas cada 100 mm · promedio de cuenca ponderado por área "
                 f"{P_MEDIA:.0f} mm/año", fontsize=8)
    cb.ax.tick_params(labelsize=7.5)
    ax.text(0.0, -0.22, "Fuente de apoyo: CHIRPS v2.0 mensual, 0.05° (Funk et al. 2015; doi:10.15780/G2RP4Q),\n"
            "interpolado de forma bilineal a la malla del DEM: la representación es continua, pero el\n"
            "detalle real es de 0.05° (~5.5 km). CHIRPS subestima frente a la P$_L$ del análisis,\n"
            f"CAMELS-PE ({P_CAMELS:.0f} mm/año, sin malla espacial): el mapa muestra el patrón espacial,\n"
            "no la magnitud. Periodo común P$_L$–P$_I$; media anual = Σ de las 12 medias mensuales.",
            transform=ax.transAxes, fontsize=7.5, va="top")
    return {}


ASP = {i: (n, c) for i, (n, c) in enumerate(zip(
    ["N", "NE", "E", "SE", "S", "SO", "O", "NO", "Plano (< 2°)"],
    ["#3288bd", "#66c2a5", "#abdda4", "#e6f598", "#fee08b", "#fdae61", "#f46d43", "#d53e4f",
     "#dddddd"]))}
themes = {
    "a_orientacion": lambda ax: cat_layer(ax, asp_cls, ASP, "(a) Orientación de laderas – GLO-30",
                                          "Ladera orientada hacia (% del área)"),
    "b_cobertura": lambda ax: cat_layer(ax, lc, {k: v for k, v in WC.items() if (lc[mask] == k).any()},
                                        "(b) Cobertura – ESA WorldCover 2021 (10 m)",
                                        "Cobertura (% del área)"),
    "c_hidrogeologia": hidro_layer,
    "d_precipitacion": precip_layer,
}
COMMON = [Line2D([], [], color="black", lw=1.3, label="Límite de cuenca"),
          Line2D([], [], color="#08306b", lw=1.6, label="Río Chicama"),
          Line2D([], [], color="#1f5fbf", lw=0.6, label="Red de drenaje (Strahler ≥ 3)"),
          Line2D([], [], color="red", marker="o", lw=0, markeredgecolor="k",
                 label="Punto de cierre El Tambo")]

# Elimina mapas de versiones anteriores que ya no forman parte del producto
for old in ["mapa_a_elevacion.png", "mapa_b_pendiente.png", "mapa_c_orientacion.png",
            "mapa_d_cobertura.png", "mapa_e_geologia.png", "mapa_f_suelos.png",
            "mapa_a_pendiente.png", "mapa_b_orientacion.png", "mapa_c_cobertura.png",
            "mapa_d_geologia.png", "mapa_e_precipitacion.png", "mapa_c_geologia.png"]:
    (OUT / old).unlink(missing_ok=True)

pcts = {}
for name, fn in themes.items():
    fig, ax = plt.subplots(figsize=(8.5, 10.5))
    pcts[name] = fn(ax)
    fig.legend(handles=COMMON, loc="lower right", fontsize=7, frameon=False,
               bbox_to_anchor=(0.98, 0.0))
    fig.savefig(OUT / f"mapa_{name}.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

# Mapa de apoyo (fuera del panel): unidades individuales con los colores oficiales de INGEMMET
fig, ax = plt.subplots(figsize=(8.5, 10.5))
geo_layer(ax)
fig.legend(handles=COMMON, loc="lower right", fontsize=7, frameon=False, bbox_to_anchor=(0.98, 0.0))
fig.savefig(OUT / "mapa_geologia_unidades_ingemmet.png", dpi=200, bbox_inches="tight")
plt.close(fig)

# Panel compuesto 1 x 4 (leyendas de cada mapa debajo de él)
fig, axs = plt.subplots(1, 4, figsize=(31, 14.5))
fig.subplots_adjust(left=0.03, right=0.99, top=0.87, bottom=0.30, wspace=0.10)
for ax, fn in zip(axs.flat, themes.values()):
    fn(ax)
fig.legend(handles=COMMON, loc="upper center", ncol=4, fontsize=11, frameon=False,
           bbox_to_anchor=(0.5, 0.935))
fig.suptitle("Contexto físico de la cuenca del río Chicama hasta la estación El Tambo "
             f"(CAMELS-PE PE_47E0A2A8, {A_km2:,.0f} km²)", fontsize=17, y=0.99)
fig.text(0.5, 0.955, CRS_TXT + "  ·  Límite, red y río Chicama derivados de Copernicus GLO-30 (D8)",
         ha="center", va="top", fontsize=10)
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
for key, cls in [("a_orientacion", ASP), ("b_cobertura", WC), ("c_hidrogeologia", GEO_CLASSES)]:
    rows += [(key[2:], cls[c][0], p) for c, p in pcts[key].items()]
p_in = precip[mask]
rows += [("precipitacion_chirps", f"{a}-{b} mm/año", 100 * np.mean((p_in >= a) & (p_in < b)))
         for a, b in [(0, 500), (500, 750), (750, 1000), (1000, 1500)]]
ficha = pd.DataFrame(rows, columns=["tema", "clase", "pct_area"]).round(2)
ficha.to_csv(OUT / "ficha_contexto_fisico.csv", index=False)
pd.set_option("display.width", 220, "display.max_colwidth", 110)
print(f"Área {A_km2:.1f} km2 | integral hipsométrica {hi:.2f} | pendiente media {s.mean():.1f}°")
print(ficha.to_string(index=False))
print(f"Salidas en {OUT}")
