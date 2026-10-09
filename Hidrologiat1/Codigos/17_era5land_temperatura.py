"""
Temperatura del aire a 2 m de ERA5-Land para la cuenca El Tambo y comparación con la temperatura
de CAMELS-PE (`tmean`, PISCOt v1.2), sin unir las dos fuentes en una sola serie (guía 3.1).

Producto : ERA5-Land monthly averaged data from 1950 to present (Copernicus Climate Change Service,
           doi:10.24381/cds.68d2bb30), product_type monthly_averaged_reanalysis, variable
           2m_temperature (K), malla 0.1° (~9 km). Es un reanálisis, no una observación.
Acceso   : Copernicus Climate Data Store (CDS) con la librería cdsapi. Requiere cuenta en
           https://cds.climate.copernicus.eu, aceptar la licencia de ERA5-Land y un archivo
           ~/.cdsapirc (Windows: %USERPROFILE%\\.cdsapirc) con
               url: https://cds.climate.copernicus.eu/api
               key: <token personal>
           Las credenciales NO se guardan en el proyecto ni en el ZIP de entrega.
Recorte  : se pide solo el rectángulo de la cuenca + 1 celda (pocos MB).
Promedio : media ponderada por el área de intersección celda–cuenca (polígono GLO-30, UTM 17S),
           igual que IMERG y PISCOp v3.0. T[°C] = T[K] − 273.15 (media mensual, no suma).

Salidas:
  Tarea 1/El_Tambo/datos/era5land/ERA5-Land_t2m_mensual_El_Tambo.nc  (+ PROCEDENCIA.json)
  Tarea 1/El_Tambo/datos/era5land/pesos_era5land_cuenca.csv
  Tarea 1/El_Tambo/resultados/series_mensuales/era5land_t2m_cuenca_mensual.csv
  Tarea 1/El_Tambo/resultados/series_mensuales/temperatura_piscot_extendida_era5land.csv
      (PISCOt hasta 2020-12 y, en columna aparte, 2021–2025 estimada con ERA5-Land + sesgo mensual)
  Tarea 1/El_Tambo/resultados/series_mensuales/fig_T_camels_vs_era5land.png
"""
from pathlib import Path
import json
import sys
import time
import zipfile

import numpy as np
import pandas as pd
import geopandas as gpd
import xarray as xr
from shapely.geometry import box
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

from periodo_comun import periodo_comun

ROOT = Path(__file__).resolve().parents[2] / "Tarea 1" / "El_Tambo"
BASIN = ROOT / "resultados" / "delimitacion" / "cuenca_el_tambo_glo30.gpkg"
SERIES = ROOT / "resultados" / "series_mensuales" / "series_mensuales_camels_pe.csv"
DL = ROOT / "datos" / "era5land"
OUT = ROOT / "resultados" / "series_mensuales"
DL.mkdir(parents=True, exist_ok=True)
NC = DL / "ERA5-Land_t2m_mensual_El_Tambo.nc"
RES = 0.1
INICIO, FIN = 1981, 2025
DATASET = "reanalysis-era5-land-monthly-means"

basin = gpd.read_file(BASIN).to_crs(4326)
geom = basin.geometry.iloc[0]
x0, y0, x1, y1 = geom.bounds
# Rectángulo de descarga alineado a la malla de 0.1° con una celda de margen
N, W = np.ceil(y1 / RES) * RES + RES, np.floor(x0 / RES) * RES - RES
S, E = np.floor(y0 / RES) * RES - RES, np.ceil(x1 / RES) * RES + RES

# ------------------------------------------------------------------ descarga (CDS)
if not NC.exists():
    try:
        import cdsapi
    except ImportError:
        sys.exit("Falta la librería cdsapi (pip install cdsapi).")
    if not (Path.home() / ".cdsapirc").exists():
        sys.exit(f"No existe {Path.home() / '.cdsapirc'}. Crear el archivo con url y key del CDS "
                 "(ver cabecera del script) y aceptar la licencia de ERA5-Land en el portal.")
    req = {"product_type": ["monthly_averaged_reanalysis"], "variable": ["2m_temperature"],
           "year": [str(y) for y in range(INICIO, FIN + 1)], "month": [f"{m:02d}" for m in range(1, 13)],
           "time": ["00:00"], "area": [round(N, 2), round(W, 2), round(S, 2), round(E, 2)],
           "data_format": "netcdf", "download_format": "unarchived"}
    t0 = time.time()
    tmp = DL / "descarga_cds.tmp"
    cdsapi.Client().retrieve(DATASET, req).download(str(tmp))     # una falla queda visible
    if zipfile.is_zipfile(tmp):                                     # por si el CDS entrega ZIP
        with zipfile.ZipFile(tmp) as z:
            nc = [n for n in z.namelist() if n.endswith(".nc")]
            if len(nc) != 1:
                raise RuntimeError(f"ZIP del CDS con contenido inesperado: {z.namelist()}")
            NC.write_bytes(z.read(nc[0]))
        tmp.unlink()
    else:
        tmp.rename(NC)
    (DL / "PROCEDENCIA.json").write_text(json.dumps(
        {"producto": "ERA5-Land monthly averaged data from 1950 to present", "dataset_cds": DATASET,
         "doi": "10.24381/cds.68d2bb30", "proveedor": "Copernicus Climate Change Service (C3S), ECMWF",
         "solicitud": req, "descargado": time.strftime("%Y-%m-%d")}, indent=2, ensure_ascii=False),
        encoding="utf-8")
    print(f"Descargado {NC.name} ({NC.stat().st_size / 2**20:.1f} MB) en {time.time() - t0:.0f} s")

# ------------------------------------------------------------------ lectura y promedio de cuenca
ds = xr.open_dataset(NC)
tdim = "valid_time" if "valid_time" in ds.dims else "time"
da = ds["t2m"].rename({tdim: "time", "latitude": "lat", "longitude": "lon"}).sortby("lat")
if "expver" in da.dims:
    da = da.max("expver")
units = da.attrs.get("units", "")
if units != "K":
    raise ValueError(f"Unidades inesperadas de t2m: '{units}'")
fechas = pd.DatetimeIndex(da.time.values).to_period("M").to_timestamp()
da = da.assign_coords(time=fechas)
esperados = pd.date_range(f"{INICIO}-01-01", f"{FIN}-12-01", freq="MS")
faltan = esperados.difference(fechas)
print(f"ERA5-Land: {da.sizes['time']} meses {fechas.min():%Y-%m} a {fechas.max():%Y-%m} | faltantes: "
      f"{', '.join(f'{t:%Y-%m}' for t in faltan) if len(faltan) else 'ninguno'} | malla {da.sizes['lat']}×{da.sizes['lon']}")

lons, lats = da.lon.values, da.lat.values
g = gpd.GeoDataFrame({"lon": np.tile(lons, len(lats)), "lat": np.repeat(lats, len(lons))},
                     geometry=[box(lo - RES / 2, la - RES / 2, lo + RES / 2, la + RES / 2)
                               for la in lats for lo in lons], crs=4326).to_crs(32717)
g["area_en_cuenca_km2"] = g.intersection(basin.to_crs(32717).geometry.iloc[0]).area / 1e6
g["fraccion_celda"] = g["area_en_cuenca_km2"] / (g.area / 1e6)
g[g["area_en_cuenca_km2"] > 0].drop(columns="geometry").round(5).to_csv(DL / "pesos_era5land_cuenca.csv", index=False)
W_ = xr.DataArray(g["area_en_cuenca_km2"].values.reshape(len(lats), len(lons)),
                  coords={"lat": lats, "lon": lons}, dims=("lat", "lon"))
print(f"Celdas en la cuenca: {int((W_ > 0).sum())} | completas (≥ 99 %): {(g['fraccion_celda'] >= 0.99).sum()} | "
      f"área {float(W_.sum()):.1f} km²")
dentro = da.where(W_ > 0)
n_nan = int(dentro.isnull().sum()) - int((W_ == 0).sum()) * da.sizes["time"]
print(f"Valores faltantes en celdas de la cuenca: {n_nan}")
T5 = ((da * W_).sum(("lat", "lon")) / W_.where(da.notnull()).sum(("lat", "lon"))).to_series() - 273.15
T5.name = "T_era5land_C"
T5.rename_axis("fecha").round(3).to_frame().to_csv(OUT / "era5land_t2m_cuenca_mensual.csv")

# ------------------------------------------------------------------ comparación con CAMELS-PE (PISCOt)
TC = pd.read_csv(SERIES, index_col="fecha", parse_dates=True)["T_C"].rename("T_camels_C")
PC_INI, PC_FIN = periodo_comun()                   # periodo común P_L–P_I (todas las figuras y ajustes)
df = pd.concat([TC, T5], axis=1).loc[PC_INI:PC_FIN]
c = df.dropna()
d = c["T_era5land_C"] - c["T_camels_C"]
clim = c.groupby(c.index.month).mean()
an = c - c.groupby(c.index.month).transform("mean")
r_an = an.corr().iloc[0, 1]
print(f"\nPeriodo común P_L–P_I {PC_INI:%Y-%m} a {PC_FIN:%Y-%m}; meses con ambas temperaturas: "
      f"{c.index.min():%Y-%m} a {c.index.max():%Y-%m} ({len(c)} meses)")
print(f"Medias: CAMELS-PE {c['T_camels_C'].mean():.2f} °C | ERA5-Land {c['T_era5land_C'].mean():.2f} °C | "
      f"sesgo ERA5L−CAMELS {d.mean():+.2f} °C | MAE {d.abs().mean():.2f} | RMSE {np.sqrt((d**2).mean()):.2f}")
print(f"r (series) {c.corr().iloc[0, 1]:.3f} | r (anomalías respecto al ciclo anual) {r_an:.3f} | "
      f"amplitud del ciclo anual: CAMELS {np.ptp(clim['T_camels_C']):.2f} °C, ERA5-Land {np.ptp(clim['T_era5land_C']):.2f} °C")
print("Sesgo por mes calendario (°C):", (clim["T_era5land_C"] - clim["T_camels_C"]).round(2).to_dict())

# ------------------------------------------------------------------ extensión de PISCOt a 2021–2025
# T_est(t) = T_ERA5L(t) + sesgo_j, con sesgo_j = media(T_PISCOt − T_ERA5L) del mes calendario j en
# los meses del periodo común con PISCOt (1998–2020). Validación fuera de muestra: sesgo ajustado en
# 1998–2010 y evaluado en 2011–2020.
def sesgo_mensual(x):
    return (x["T_camels_C"] - x["T_era5land_C"]).groupby(x.index.month).mean()


aj, ev = c.loc[:"2010"], c.loc["2011":]
est_ev = ev["T_era5land_C"] + ev.index.month.map(sesgo_mensual(aj)).values
e = est_ev - ev["T_camels_C"]
clim_aj = aj["T_camels_C"].groupby(aj.index.month).mean()
e_clim = ev.index.month.map(clim_aj).values - ev["T_camels_C"]
print(f"\nValidación 2011–2020 (sesgo ajustado en {aj.index.min():%Y}–2010, {len(ev)} meses): sesgo {e.mean():+.2f} °C | "
      f"MAE {e.abs().mean():.2f} | RMSE {np.sqrt((e**2).mean()):.2f} | r {np.corrcoef(est_ev, ev['T_camels_C'])[0, 1]:.2f}"
      f" · referencia climatología {aj.index.min():%Y}–2010: MAE {np.abs(e_clim).mean():.2f} °C")
pd.Series({"ajuste": f"{aj.index.min():%Y-%m} a {aj.index.max():%Y-%m}",
           "evaluacion": f"{ev.index.min():%Y-%m} a {ev.index.max():%Y-%m}", "n_eval": len(ev),
           "sesgo_C": round(e.mean(), 3), "MAE_C": round(e.abs().mean(), 3), "RMSE_C": round(np.sqrt((e**2).mean()), 3),
           "r": round(np.corrcoef(est_ev, ev["T_camels_C"])[0, 1], 3), "MAE_climatologia_C": round(np.abs(e_clim).mean(), 3)},
          name="valor").rename_axis("indicador").to_csv(OUT / "validacion_extension_T.csv")
b = sesgo_mensual(c)
ext = pd.DataFrame({"T_piscot_C": df["T_camels_C"]})
ext["T_estimada_era5land_C"] = (df["T_era5land_C"] + df.index.month.map(b).values).where(df["T_camels_C"].isna())
ext["fuente"] = np.where(ext["T_piscot_C"].notna(), "PISCOt v1.2 (CAMELS-PE)",
                         np.where(ext["T_estimada_era5land_C"].notna(),
                                  f"estimada: ERA5-Land + sesgo mensual {c.index.min():%Y}-{c.index.max():%Y}", ""))
ext.rename_axis("fecha").round(3).to_csv(OUT / "temperatura_piscot_extendida_era5land.csv")
print(f"Meses estimados: {ext['T_estimada_era5land_C'].notna().sum()} "
      f"({ext['T_estimada_era5land_C'].first_valid_index():%Y-%m} a {ext['T_estimada_era5land_C'].last_valid_index():%Y-%m}); "
      f"sesgo aplicado {b.min():+.2f} a {b.max():+.2f} °C")

# ------------------------------------------------------------------ ¿explica la altitud el sesgo?
# Relieve del modelo (geopotencial invariante de ERA5-Land / g) frente a la elevación media de la
# cuenca (atributo elev_mean de CAMELS-PE, FABDEM; ≈ GLO-30). El archivo trae la longitud en 0–360°.
GP = DL / "ERA5-Land_geopotencial_El_Tambo.nc"
if not GP.exists():
    import cdsapi
    cdsapi.Client().retrieve(DATASET, {"product_type": ["monthly_averaged_reanalysis"], "variable": ["geopotential"],
                                       "year": ["2000"], "month": ["01"], "time": ["00:00"],
                                       "area": [round(N, 2), round(W, 2), round(S, 2), round(E, 2)],
                                       "data_format": "netcdf", "download_format": "unarchived"}).download(str(GP))
zm = xr.open_dataset(GP)["z"].squeeze(drop=True).rename(latitude="lat", longitude="lon") / 9.80665
zm = zm.assign_coords(lon=((zm.lon + 180) % 360) - 180).sortby(["lat", "lon"])
zm = zm.reindex(lat=lats, lon=lons, method="nearest", tolerance=1e-3)
z_era = float((zm * W_).sum() / W_.sum())
attrs = pd.read_csv(ROOT / "datos" / "camels_pe" / "PE_47E0A2A8_atributos.csv", index_col=0)["valor"]
z_real = float(attrs["elev_mean"])
GRAD = -6.5e-3                                                    # gradiente térmico estándar (°C/m)
print(f"\nElevación media: ERA5-Land {z_era:.0f} m | cuenca (CAMELS-PE elev_mean) {z_real:.0f} m | diferencia "
      f"{z_era - z_real:+.0f} m → con {GRAD * 1e3:.1f} °C/km explica {GRAD * (z_era - z_real):+.2f} °C "
      f"de los {d.mean():+.2f} °C de sesgo")

# ------------------------------------------------------------------ figura
ROJO_T, MORADO = "#c8553d", "#6a4c93"
fig, (a1, a2) = plt.subplots(1, 2, figsize=(16, 5.2), gridspec_kw={"width_ratios": [3.2, 1], "wspace": 0.15})
a1.plot(df.index, df["T_camels_C"], color=ROJO_T, lw=1.1, label="CAMELS-PE `tmean` (PISCOt v1.2)")
a1.plot(df.index, df["T_era5land_C"], color=MORADO, lw=1.1, alpha=0.9, label="ERA5-Land T 2 m (promedio ponderado por área)")
a1.set_ylabel("Temperatura media mensual (°C)")
a1.grid(alpha=0.3)
a1.xaxis.set_major_locator(mdates.YearLocator(2))
a1.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
a1.set_xlim(df.index.min() - pd.Timedelta(days=60), df.index.max() + pd.Timedelta(days=60))
a1.legend(loc="lower left", bbox_to_anchor=(0, 1.0), frameon=False, ncol=2)
a1.set_title("(a) Serie mensual", loc="left", fontweight="bold", pad=24)
RESUMEN = (f"Meses con ambas fuentes {c.index.min():%Y-%m} a {c.index.max():%Y-%m} (n = {len(c)}): sesgo ERA5-Land − "
           f"CAMELS = {d.mean():+.2f} °C · MAE {d.abs().mean():.2f} · RMSE {np.sqrt((d**2).mean()):.2f} °C · "
           f"r = {c.corr().iloc[0, 1]:.2f} · r anomalías = {r_an:.2f}")      # va en la nota (no tapa datos)
meses = np.arange(1, 13)
a2.plot(meses, clim["T_camels_C"], color=ROJO_T, marker="o", ms=4, lw=1.4)
a2.plot(meses, clim["T_era5land_C"], color=MORADO, marker="o", ms=4, lw=1.4)
a2.set_xticks(meses, ["E", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"])
a2.set_ylabel("°C")
a2.grid(alpha=0.3)
a2.set_title(f"(b) Ciclo anual medio {c.index.min():%Y}–{c.index.max():%Y}", loc="left", fontweight="bold", pad=24)
fig.suptitle("Cuenca del río Chicama hasta El Tambo – temperatura media mensual: CAMELS-PE vs ERA5-Land",
             fontsize=13, y=1.03)
fig.text(0.01, -0.06, "ERA5-Land (C3S, doi:10.24381/cds.68d2bb30): reanálisis a 0.1°, T[°C] = T[K] − 273.15. CAMELS-PE "
         f"`tmean`: PISCOt v1.2 (SENAMHI), disponible hasta 2020-12. Periodo común P$_L$–P$_I$ {PC_INI:%Y-%m} a "
         f"{PC_FIN:%Y-%m}. Las fuentes se comparan, no se unen.\n" + RESUMEN,
         fontsize=8)
fig.savefig(OUT / "fig_T_camels_vs_era5land.png", dpi=200, bbox_inches="tight")
print(f"Figura: {OUT / 'fig_T_camels_vs_era5land.png'}")
