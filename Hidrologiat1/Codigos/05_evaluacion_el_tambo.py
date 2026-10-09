"""
Evaluación detallada de la cuenca El Tambo (río Chicama, CAMELS-PE PE_47E0A2A8)
frente a los criterios de selección de la Tarea 1:
  1. Registro común >= 25 años de caudal observado y precipitación diarios,
     <= 10 % de días faltantes por variable (antes de rellenar) y su distribución temporal.
  2. Área 100-10000 km2 y representatividad frente a la malla IMERG (0.1°).
  3. Disponibilidad de temperatura (PISCOt) y necesidad de ERA5-Land.
  4. Controles de calidad preliminares (duplicados, negativos, ceros, constantes, extremos,
     balance anual P-R).
Salidas en Tarea 1/seleccion_PE/el_tambo/.
"""
from pathlib import Path
import warnings

import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from shapely.geometry import box

warnings.filterwarnings("ignore", message="Geometry is in a geographic CRS")

GID = "PE_47E0A2A8"
ROOT = Path(__file__).resolve().parents[2] / "Tarea 1"
BASE = ROOT / "CAMELS_PE" / "CAMELS-PE"
OUT = ROOT / "seleccion_PE" / "el_tambo"
OUT.mkdir(parents=True, exist_ok=True)
UTM = 32717  # WGS84 / UTM 17S

d = pd.read_csv(BASE / "03_timeseries" / "by_catchment" / f"{GID}.csv", parse_dates=["date"])

# ---------------------------------------------------------------- 1. Controles básicos
print("=== Controles básicos ===")
print("Rango de fechas:", d["date"].min().date(), "->", d["date"].max().date())
print("Fechas duplicadas:", int(d["date"].duplicated().sum()))
full = pd.date_range(d["date"].min(), d["date"].max(), freq="D")
print("Días ausentes del índice:", len(full.difference(d["date"])))
d = d.set_index("date")
for v in ["prec", "flow_obs", "tmean"]:
    s = d[v]
    print(f"{v:9s} NA={int(s.isna().sum()):6d}  neg={int((s < 0).sum())}  ceros={int((s == 0).sum())}"
          f"  min={s.min():.3f}  max={s.max():.3f}")

# ---------------------------------------------------------------- 2. Huecos de caudal
q = d["flow_obs"]
na = q.isna()
grp = (na != na.shift()).cumsum()
gaps = (pd.DataFrame({"na": na, "g": grp}).reset_index()
        .query("na").groupby("g")["date"].agg(["min", "max", "count"]))
gaps.columns = ["inicio", "fin", "dias"]
gaps = gaps[gaps["inicio"] >= "1998-01-01"].sort_values("dias", ascending=False)
gaps.to_csv(OUT / "huecos_caudal.csv", index=False)
print("\n=== Huecos de caudal desde 1998 (10 más largos) ===")
print(gaps.head(10).to_string(index=False))
print(f"Total huecos: {len(gaps)}; huecos de 1 día: {(gaps['dias'] == 1).sum()}")

# ---------------------------------------------------------------- 3. Ventanas candidatas
print("\n=== Ventanas candidatas ===")
wins = [("1998-01-01", "2025-12-31"), ("1998-01-01", "2024-12-31"), ("1998-01-01", "2025-09-30"), ("1998-01-01", "2020-12-31")]
for a, b in wins:
    w = d.loc[a:b]
    m = w["flow_obs"].isna().resample("MS").mean()
    print(f"{a}..{b}: {len(w)/365.25:.1f} años | faltantes Q={100*w['flow_obs'].isna().mean():.2f}% "
          f"P={100*w['prec'].isna().mean():.2f}% T={100*w['tmean'].isna().mean():.1f}% | "
          f"meses Q con >10% faltante: {int((m > 0.10).sum())}/{len(m)}; "
          f"con 0% faltante: {int((m == 0).sum())}")

# ---------------------------------------------------------------- 4. Disponibilidad año x mes
W = d.loc["1998":"2025"]
valid = W["flow_obs"].notna().groupby([W.index.year, W.index.month]).sum().unstack()
ndays = W["flow_obs"].groupby([W.index.year, W.index.month]).size().unstack()
frac = valid / ndays
valid.to_csv(OUT / "dias_validos_caudal_anio_mes.csv")
yr = (100 * W["flow_obs"].isna().groupby(W.index.year).mean()).round(1)
print("\n=== % faltante de caudal por año ===")
print(" ".join(f"{y}:{v}" for y, v in yr.items()))

fig, ax = plt.subplots(figsize=(7, 8))
im = ax.imshow(frac.values, aspect="auto", cmap="viridis", vmin=0, vmax=1)
ax.set_xticks(range(12), ["E", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"])
ax.set_yticks(range(len(frac)), frac.index)
ax.set_title(f"El Tambo ({GID}) – fracción de días con caudal observado")
fig.colorbar(im, ax=ax, label="fracción de días válidos")
fig.tight_layout(); fig.savefig(OUT / "disponibilidad_caudal.png", dpi=150); plt.close(fig)

# ---------------------------------------------------------------- 5. Calidad del caudal
runs = (q.diff() != 0).cumsum()
rl = q.dropna().groupby(runs.dropna()).size()
print("\n=== Calidad de caudal ===")
print("Secuencias constantes >= 7 días:", int((rl >= 7).sum()), "| máx:", int(rl.max()))
top = q.nlargest(5)
print("5 máximos diarios (mm/d):", ", ".join(f"{t.date()}={v:.1f}" for t, v in top.items()))

# Balance anual en años con >= 95 % de días de caudal
ann = W.groupby(W.index.year).agg(P=("prec", "sum"), R=("flow_obs", "sum"),
                                   nQ=("flow_obs", "count"), n=("prec", "size"))
ann = ann[ann["nQ"] / ann["n"] >= 0.95]
ann["R/P"] = (ann["R"] / ann["P"]).round(2)
ann.to_csv(OUT / "balance_anual.csv")
print("\n=== Balance anual (años con >=95 % de caudal) ===")
print(ann[["P", "R"]].round(0).join(ann["R/P"]).to_string())

# Climatología rápida (solo orientativa; meses con <=10 % faltantes)
M = pd.DataFrame({"P": W["prec"].resample("MS").sum(),
                  "R": W["flow_obs"].resample("MS").sum(),
                  "fQ": W["flow_obs"].isna().resample("MS").mean(),
                  "T": W["tmean"].resample("MS").mean()})
M.loc[M["fQ"] > 0.10, "R"] = np.nan
clim = M.groupby(M.index.month)[["P", "R", "T"]].mean().round(1)
print("\n=== Climatología orientativa 1998-2025 (mm/mes, °C) ===")
print(clim.T.to_string())

# ---------------------------------------------------------------- 6. Escalas espaciales
cat = gpd.read_file(BASE / "04_geospatial" / "camels_pe_catchments.gpkg")
geom = cat.loc[cat["gauge_id"] == GID].geometry.iloc[0]
A = gpd.GeoSeries([geom], crs=4326).to_crs(UTM).area.iloc[0] / 1e6
minx, miny, maxx, maxy = geom.bounds
print(f"\n=== Escalas espaciales ===\nÁrea del polígono (UTM 17S): {A:.1f} km2")
print(f"Extensión: lon {minx:.2f}..{maxx:.2f}, lat {miny:.2f}..{maxy:.2f}")


def grid_stats(res, origin=0.0, name=""):
    xs = np.arange(np.floor((minx - origin) / res) * res + origin, maxx, res)
    ys = np.arange(np.floor((miny - origin) / res) * res + origin, maxy, res)
    cells = gpd.GeoSeries([box(x, y, x + res, y + res) for x in xs for y in ys], crs=4326)
    f = (cells.intersection(geom).to_crs(UTM).area / cells.to_crs(UTM).area).values
    cell_km2 = cells.to_crs(UTM).area.mean() / 1e6
    inside = f[f > 0]
    print(f"{name:28s} celda≈{cell_km2:6.0f} km2 | área/celda={A/cell_km2:6.1f} | "
          f"celdas que tocan={len(inside):3d} | completas(>=99%)={int((inside >= .99).sum()):3d} | "
          f"peso en celdas de borde={100*(1 - inside[inside >= .99].sum()/inside.sum()):.0f}%")
    return cells, f


imerg_cells, imerg_f = grid_stats(0.1, 0.0, "IMERG V07 (0.1°)")
grid_stats(0.1, 0.0, "ERA5-Land (0.1°)")
grid_stats(0.25, 0.125, "ERA5 (0.25°)")
grid_stats(0.1, 0.0, "PISCOp v2.1 (0.1°)")

# Mapa de la cuenca con la malla IMERG
gdf = gpd.GeoDataFrame({"frac": imerg_f}, geometry=imerg_cells.values, crs=4326)
gdf = gdf[gdf["frac"] > 0]
gauges = gpd.read_file(BASE / "04_geospatial" / "camels_pe_gauges.gpkg")
fig, ax = plt.subplots(figsize=(7, 6))
gdf.plot(ax=ax, column="frac", cmap="Blues", edgecolor="grey", linewidth=0.4, legend=True,
         legend_kwds={"label": "fracción de la celda IMERG dentro de la cuenca"})
gpd.GeoSeries([geom], crs=4326).boundary.plot(ax=ax, color="k", linewidth=1.2)
g = gauges[gauges["gauge_id"] == GID]
g.plot(ax=ax, color="red", markersize=40, zorder=5)
ax.set_title(f"El Tambo ({GID}), {A:.0f} km² – malla IMERG 0.1°")
ax.set_xlabel("Longitud (°)"); ax.set_ylabel("Latitud (°)")
fig.tight_layout(); fig.savefig(OUT / "mapa_malla_imerg.png", dpi=150); plt.close(fig)
gdf.to_file(OUT / "pesos_imerg.gpkg")
print(f"\nFiguras y tablas en {OUT}")
