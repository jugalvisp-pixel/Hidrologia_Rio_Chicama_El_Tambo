"""
Tamizado de cuencas CAMELS-PE v1.0.1 según los requisitos de la Tarea 1:
  - Área de drenaje entre 100 y 10000 km2 (incluidos).
  - Caudal observado (flow_obs) y precipitación (prec) diarios con un periodo
    común de al menos 25 años y como máximo 10 % de días faltantes por variable
    (calculado sobre los días esperados, antes de cualquier relleno).
  - Superposición con IMERG V07 Final (desde 2000-06-01) y número de celdas
    IMERG de 0.1° que cubren la cuenca.
"""
from pathlib import Path

import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import box

BASE = Path(__file__).resolve().parents[2] / "Tarea 1" / "CAMELS_PE" / "CAMELS-PE"
OUT = Path(__file__).resolve().parents[2] / "Tarea 1" / "seleccion_PE"
OUT.mkdir(exist_ok=True)

MIN_YEARS = 25
MAX_MISS = 0.10
AREA_MIN, AREA_MAX = 100, 10000
IMERG_START = pd.Timestamp("2000-06-01")
IMERG_RES = 0.1

st = pd.read_csv(BASE / "01_metadata" / "stations.csv")
topo = pd.read_csv(BASE / "02_attributes" / "topographic_attributes.csv")
clim = pd.read_csv(BASE / "02_attributes" / "climatic_indices.csv")
hum = pd.read_csv(BASE / "02_attributes" / "human_intervention_attributes.csv")
hyd = pd.read_csv(BASE / "02_attributes" / "hydrological_signatures.csv")
meta = (st.merge(topo, on="gauge_id").merge(clim, on="gauge_id")
          .merge(hum, on="gauge_id").merge(hyd, on="gauge_id"))

ts = pd.read_csv(BASE / "03_timeseries" / "timeseries.csv",
                 usecols=["date", "gauge_id", "prec", "flow_obs"],
                 parse_dates=["date"])


def best_window(g):
    """Ventana de años calendario completos más larga (y luego más reciente)
    en la que prec y flow_obs tienen <= 10 % de días faltantes cada una."""
    g = g.set_index("date").sort_index()
    yrs = g.index.year
    miss_q = g["flow_obs"].isna().groupby(yrs).sum()
    miss_p = g["prec"].isna().groupby(yrs).sum()
    ndays = g.groupby(yrs).size()
    years = ndays.index.values
    best = None
    for i in range(len(years)):
        for j in range(i + MIN_YEARS - 1, len(years)):
            sl = slice(years[i], years[j])
            n = ndays.loc[sl].sum()
            fq, fp = miss_q.loc[sl].sum() / n, miss_p.loc[sl].sum() / n
            if fq <= MAX_MISS and fp <= MAX_MISS:
                L = years[j] - years[i] + 1
                if best is None or (L, years[j]) > (best[1] - best[0] + 1, best[1]):
                    best = (years[i], years[j], fq, fp)
    return best


def longest_gap(s):
    na = s.isna().values
    run = mx = 0
    for v in na:
        run = run + 1 if v else 0
        mx = max(mx, run)
    return mx


meta_ok = meta[(meta["area"] >= AREA_MIN) & (meta["area"] <= AREA_MAX)]
rows = []
for gid in meta_ok["gauge_id"]:
    g = ts[ts["gauge_id"] == gid]
    w = best_window(g)
    if w is None:
        continue
    y0, y1, fq, fp = w
    gw = g[(g["date"].dt.year >= y0) & (g["date"].dt.year <= y1)].set_index("date")
    # Meses con <= 10 % de días faltantes de caudal (criterio provisional de completitud)
    mfrac = gw["flow_obs"].isna().resample("MS").mean()
    imerg = gw[gw.index >= IMERG_START]
    rows.append(dict(gauge_id=gid, y0=y0, y1=y1, n_years=y1 - y0 + 1,
                     miss_q_pct=100 * fq, miss_p_pct=100 * fp,
                     max_gap_q_days=longest_gap(gw["flow_obs"]),
                     months_q_ok=int((mfrac <= 0.10).sum()), months_total=len(mfrac),
                     imerg_overlap_years=round(len(imerg) / 365.25, 1),
                     miss_q_imerg_pct=100 * imerg["flow_obs"].isna().mean()))

res = pd.DataFrame(rows).merge(meta_ok, on="gauge_id")

# Celdas IMERG (malla 0.1°, bordes en múltiplos de 0.1°) que intersectan cada cuenca
cat = gpd.read_file(BASE / "04_geospatial" / "camels_pe_catchments.gpkg")
cat = cat[cat["gauge_id"].isin(res["gauge_id"])]
ncell, nfull = {}, {}
for _, r in cat.iterrows():
    minx, miny, maxx, maxy = r.geometry.bounds
    xs = np.arange(np.floor(minx / IMERG_RES) * IMERG_RES, maxx, IMERG_RES)
    ys = np.arange(np.floor(miny / IMERG_RES) * IMERG_RES, maxy, IMERG_RES)
    cells = gpd.GeoSeries([box(x, y, x + IMERG_RES, y + IMERG_RES) for x in xs for y in ys], crs=4326)
    frac = cells.intersection(r.geometry).area / cells.area
    ncell[r["gauge_id"]] = int((frac > 0).sum())
    nfull[r["gauge_id"]] = int((frac >= 0.99).sum())
res["imerg_cells_any"] = res["gauge_id"].map(ncell)
res["imerg_cells_full"] = res["gauge_id"].map(nfull)

cols = ["gauge_id", "gauge_name", "name_cat", "gauge_region", "gauge_lat", "gauge_lon",
        "area", "elev_min", "elev_max", "elev_mean", "y0", "y1", "n_years",
        "miss_q_pct", "miss_p_pct", "max_gap_q_days", "months_q_ok", "months_total",
        "imerg_overlap_years", "miss_q_imerg_pct", "imerg_cells_any", "imerg_cells_full",
        "p_mean", "aridity", "runoff_ratio", "baseflow_index", "reservoir_n", "reservoir_vol",
        "surface_adh_vol", "is_nested", "nested_group_id"]
res = res[cols].sort_values(["n_years", "miss_q_pct"], ascending=[False, True])
res.to_csv(OUT / "candidatas_PE.csv", index=False)
pd.set_option("display.width", 250, "display.max_columns", 40)
print(f"Cuencas en rango de área: {len(meta_ok)} de {len(meta)}")
print(f"Cuencas que cumplen >= {MIN_YEARS} años y <= {MAX_MISS:.0%} faltantes: {len(res)}")
print(res.round(2).to_string(index=False))
