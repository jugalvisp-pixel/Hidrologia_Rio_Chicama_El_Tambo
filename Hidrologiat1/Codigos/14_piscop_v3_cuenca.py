"""
PISCOp v3.0 diario (SENAMHI, doi:10.6084/m9.figshare.32411886.v1) para la cuenca El Tambo:
recorte de la malla nacional, promedio de cuenca ponderado por área y control de requisitos de
la guía (registro diario, ≥ 25 años comunes con el caudal, ≤ 10 % de días faltantes).

Producto : PISCOp v3.0, precipitación diaria (mm/día), malla 0.1°, 1981-01-01 .. 2025-12-31.
           Producto interpolado: pluviómetros (Perú y países vecinos) con control de calidad y
           relleno, interpolados con predictores satelitales y climatologías radar-pluviómetro;
           el diario se corrige para sumar el total mensual. No es una observación directa.
Promedio : igual que IMERG (script 10): media ponderada por el área de intersección celda–cuenca
           (polígono GLO-30, UTM 17S); pesos renormalizados si una celda no tiene dato.
Requiere : 06 (polígono) y 13 (descarga de PISCOp_d.nc y PISCOp_m.nc con MD5).

Salidas:
  Tarea 1/El_Tambo/datos/piscop_senamhi/v3p0/PISCOp_v3_d_El_Tambo_recorte.nc   (va en el ZIP)
  Tarea 1/El_Tambo/datos/piscop_senamhi/v3p0/pesos_piscop_v3_cuenca.csv
  Tarea 1/El_Tambo/resultados/series_mensuales/piscop_v3_cuenca_diaria.csv
  Tarea 1/El_Tambo/resultados/series_mensuales/piscop_v3_cuenca_mensual.csv
"""
from pathlib import Path

import numpy as np
import pandas as pd
import geopandas as gpd
import xarray as xr
from shapely.geometry import box

ROOT = Path(__file__).resolve().parents[2] / "Tarea 1" / "El_Tambo"
V3 = ROOT / "datos" / "piscop_senamhi" / "v3p0"
BASIN = ROOT / "resultados" / "delimitacion" / "cuenca_el_tambo_glo30.gpkg"
CAMELS = ROOT / "datos" / "camels_pe" / "PE_47E0A2A8.csv"
OUT = ROOT / "resultados" / "series_mensuales"
RES = 0.1
INI, FIN = "1981-01-01", "2025-12-31"
COMUN_Q = ("1998-01-01", "2025-12-31")        # periodo con caudal observado (CAMELS-PE flow_obs)


def estandarizar(ds):
    """Nombres comunes (time, lat, lon) y variable de precipitación, sea cual sea el archivo."""
    ren = {k: v for k, v in {"Z1": "time", "z": "time", "T": "time", "latitude": "lat",
                             "longitude": "lon", "X": "lon", "Y": "lat"}.items() if k in ds.dims}
    ds = ds.rename(ren)
    var = [v for v in ds.data_vars if ds[v].dims == ("time", "lat", "lon")]
    if len(var) != 1:
        raise ValueError(f"Variable de precipitación ambigua: {list(ds.data_vars)}")
    return ds[var[0]]


basin = gpd.read_file(BASIN).to_crs(4326)
geom = basin.geometry.iloc[0]
x0, y0, x1, y1 = geom.bounds

# ------------------------------------------------------------------ recorte de la malla diaria
d = estandarizar(xr.open_dataset(V3 / "PISCOp_d.nc", decode_times=False))
print(f"PISCOp_d.nc: variable '{d.name}', unidades '{d.attrs.get('units', 's/d')}', "
      f"{d.sizes['time']} pasos, malla {d.sizes['lat']}×{d.sizes['lon']}")
fechas = pd.date_range(INI, FIN, freq="D")
if d.sizes["time"] != len(fechas):
    raise ValueError(f"Se esperaban {len(fechas)} días {INI}..{FIN} y hay {d.sizes['time']}")
d = d.assign_coords(time=fechas).sortby("lat")
d = d.sel(lon=slice(x0 - RES, x1 + RES), lat=slice(y0 - RES, y1 + RES)).load()
rec = V3 / "PISCOp_v3_d_El_Tambo_recorte.nc"
d.to_dataset(name="precipitation").assign_attrs(
    fuente="PISCOp v3.0 (SENAMHI), doi:10.6084/m9.figshare.32411886.v1, archivo PISCOp_d.nc",
    recorte="rectángulo de la cuenca El Tambo (GLO-30) + 1 celda", unidades="mm/día").to_netcdf(rec)
print(f"Recorte guardado: {rec.name} ({rec.stat().st_size / 2**20:.1f} MB), "
      f"lon {d.lon.values.min():.2f}..{d.lon.values.max():.2f}, lat {d.lat.values.min():.2f}..{d.lat.values.max():.2f}")

# ------------------------------------------------------------------ pesos por área de intersección
lons, lats = d.lon.values, d.lat.values
g = gpd.GeoDataFrame({"lon": np.tile(lons, len(lats)), "lat": np.repeat(lats, len(lons))},
                     geometry=[box(lo - RES / 2, la - RES / 2, lo + RES / 2, la + RES / 2)
                               for la in lats for lo in lons], crs=4326).to_crs(32717)
g["area_en_cuenca_km2"] = g.intersection(basin.to_crs(32717).geometry.iloc[0]).area / 1e6
g["fraccion_celda"] = g["area_en_cuenca_km2"] / (g.area / 1e6)
g[g["area_en_cuenca_km2"] > 0].drop(columns="geometry").round(5).to_csv(V3 / "pesos_piscop_v3_cuenca.csv", index=False)
W = xr.DataArray(g["area_en_cuenca_km2"].values.reshape(len(lats), len(lons)),
                 coords={"lat": lats, "lon": lons}, dims=("lat", "lon"))
print(f"Celdas en la cuenca: {(W > 0).sum().item()} | completas (≥ 99 %): "
      f"{(g['fraccion_celda'] >= 0.99).sum()} | área {W.sum().item():.1f} km²")

# ------------------------------------------------------------------ controles de calidad
dc = d.where(W > 0)
nan_celda_dia = int(dc.isnull().sum().item() - (W == 0).sum().item() * d.sizes["time"])
print(f"\nControles (celdas dentro de la cuenca):")
print(f"  valores faltantes: {nan_celda_dia} | negativos: {int((dc < 0).sum())} | "
      f"máximo diario en una celda: {float(dc.max()):.1f} mm")
P = ((d * W).sum(("lat", "lon")) / W.where(d.notnull()).sum(("lat", "lon"))).to_series()
P[W.where(d.notnull()).sum(("lat", "lon")).to_series() == 0] = np.nan     # día sin ninguna celda
P.name = "P_piscop_v3_mm_d"
faltan = P.isna()
c0, c1 = COMUN_Q
print(f"  días sin dato en el promedio de cuenca: {int(faltan.sum())} de {len(P)} en {INI[:4]}–{FIN[:4]}; "
      f"{int(faltan.loc[c0:c1].sum())} de {len(P.loc[c0:c1])} "
      f"({100 * faltan.loc[c0:c1].mean():.2f} %) en el periodo común con el caudal {c0[:4]}–{c1[:4]}")
const = (P.diff() == 0) & (P > 0)
print(f"  días con valor positivo idéntico al anterior: {int(const.sum())}")

# ------------------------------------------------------------------ agregación mensual y contraste
n_validos = P.notna().resample("MS").sum()
Pm = P.resample("MS").sum(min_count=1).where(n_validos == P.resample("MS").size())   # solo meses completos
m = estandarizar(xr.open_dataset(V3 / "PISCOp_m.nc", decode_times=False))
m = m.assign_coords(time=pd.date_range(INI, periods=m.sizes["time"], freq="MS")).sortby("lat")
m = m.sel(lon=slice(lons.min() - 1e-3, lons.max() + 1e-3), lat=slice(lats.min() - 1e-3, lats.max() + 1e-3))
Pm_arch = ((m * W).sum(("lat", "lon")) / W.where(m.notnull()).sum(("lat", "lon"))).to_series()
dif = (Pm - Pm_arch).dropna()
print(f"\nΣ diario vs archivo mensual de v3.0: máx |dif| = {dif.abs().max():.2f} mm/mes, "
      f"media = {dif.mean():+.3f} mm/mes ({len(dif)} meses)")

cam = pd.read_csv(CAMELS, parse_dates=["date"], index_col="date")["prec"]
cmp_ = pd.concat([Pm.rename("v3"), cam.resample("MS").sum().rename("camels")], axis=1).loc[c0:c1].dropna()
print(f"Media anual {c0[:4]}–{c1[:4]}: PISCOp v3.0 {cmp_['v3'].mean() * 12:.0f} mm/año | "
      f"CAMELS-PE prec {cmp_['camels'].mean() * 12:.0f} mm/año | r mensual {cmp_.corr().iloc[0, 1]:.3f}")

# ------------------------------------------------------------------ salidas
P.round(3).to_frame().rename_axis("fecha").to_csv(OUT / "piscop_v3_cuenca_diaria.csv")
pd.DataFrame({"dias_validos": n_validos, "n_dias": P.resample("MS").size(),
              "P_piscop_v3_mm_mes": Pm.round(2)}).rename_axis("fecha").to_csv(OUT / "piscop_v3_cuenca_mensual.csv")
anual = P.resample("YS").sum(min_count=365)
print(f"\nP anual (cuenca) {INI[:4]}–{FIN[:4]}: media {anual.mean():.0f} mm, mín {anual.min():.0f} "
      f"({anual.idxmin().year}), máx {anual.max():.0f} ({anual.idxmax().year})")
print(f"Salidas: {OUT / 'piscop_v3_cuenca_diaria.csv'}, {OUT / 'piscop_v3_cuenca_mensual.csv'}")
