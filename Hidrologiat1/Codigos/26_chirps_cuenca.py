"""
CHIRPS v2.0 mensual (0.05°) para la cuenca El Tambo: campo de apoyo para el mapa continuo de
precipitación media anual del panel de contexto físico (script 07b) y serie de cuenca.

Producto : Climate Hazards group InfraRed Precipitation with Station data, v2.0, mensual, 0.05°
           (Funk et al. 2015, doi:10.1038/sdata.2015.66; dataset doi:10.15780/G2RP4Q). Combina una
           climatología de alta resolución que tiene en cuenta el relieve (CHPclim), estimaciones
           infrarrojas y estaciones. Es también la componente satelital de PISCOp.
Acceso   : recorte por OPeNDAP/NetCDF del IRI Data Library (SOURCES/.UCSB/.CHIRPS/.v2p0/.monthly/.global),
           mismo contenido que los GeoTIFF de UCSB. En el IRI la coordenada X marca el BORDE OESTE de la
           celda: verificado contra chirps-v2.0.2017.03.tif de data.chc.ucsb.edu (diferencia 0.000 mm
           con X + 0.025°; ~94 mm si no se corrige). Aquí se lleva X al centro de la celda.
Periodo  : periodo común P_L–P_I (periodo_comun.py).
Promedio de cuenca: ponderado por el área de intersección celda–polígono (UTM 17S), como IMERG y PISCOp.

Salidas:
  Tarea 1/El_Tambo/datos/chirps/CHIRPS_v2p0_mensual_El_Tambo.nc (+ PROCEDENCIA.json, pesos_chirps_cuenca.csv)
  Tarea 1/El_Tambo/resultados/series_mensuales/chirps_cuenca_mensual.csv
"""
from pathlib import Path
import json
import time
import urllib.request

import numpy as np
import pandas as pd
import geopandas as gpd
import xarray as xr
from shapely.geometry import box

from periodo_comun import periodo_comun

ROOT = Path(__file__).resolve().parents[2] / "Tarea 1" / "El_Tambo"
DL = ROOT / "datos" / "chirps"
DL.mkdir(parents=True, exist_ok=True)
SMES = ROOT / "resultados" / "series_mensuales"
NC = DL / "CHIRPS_v2p0_mensual_El_Tambo.nc"
INI, FIN = periodo_comun()
X0, X1, Y0, Y1 = -79.0, -78.0, -8.2, -7.2                     # recorte (cuenca + margen)
MESES_EN = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
URL = ("https://iridl.ldeo.columbia.edu/SOURCES/.UCSB/.CHIRPS/.v2p0/.monthly/.global/.precipitation/"
       f"X/{X0}/{X1}/RANGEEDGES/Y/{Y0}/{Y1}/RANGEEDGES/"
       f"T/({MESES_EN[INI.month - 1]}%20{INI.year})/({MESES_EN[FIN.month - 1]}%20{FIN.year})/RANGEEDGES/data.nc")

if not NC.exists():
    tmp = DL / "descarga.tmp"
    urllib.request.urlretrieve(URL, tmp)
    raw = xr.open_dataset(tmp, decode_times=False)
    t = raw["T"].values                                        # meses desde 1960-01 (centro del mes)
    fechas = pd.to_datetime([f"{1960 + int(m) // 12}-{int(m) % 12 + 1:02d}-01" for m in np.floor(t)])
    ds = xr.Dataset({"precipitation": (("time", "lat", "lon"), raw["precipitation"].values.astype("float32"))},
                    coords={"time": fechas, "lat": raw["Y"].values.astype("float64").round(3),
                            "lon": (raw["X"].values.astype("float64") + 0.025).round(3)}).sortby("lat")
    ds["precipitation"].attrs = {"units": "mm/month", "long_name": "CHIRPS v2.0 monthly precipitation"}
    ds.attrs = {"fuente": "CHIRPS v2.0 (UCSB/CHC) vía IRI Data Library; lon corregida al centro de celda",
                "doi": "10.15780/G2RP4Q"}
    raw.close()
    ds.to_netcdf(NC, encoding={"precipitation": {"zlib": True, "complevel": 5}})
    tmp.unlink()
    (DL / "PROCEDENCIA.json").write_text(json.dumps(
        {"producto": "CHIRPS v2.0 monthly 0.05°", "doi_dataset": "10.15780/G2RP4Q",
         "referencia": "Funk et al. 2015, Scientific Data, doi:10.1038/sdata.2015.66",
         "url": URL, "correccion": "X del IRI = borde oeste de la celda; se suma 0.025° (verificado contra "
                                   "chirps-v2.0.2017.03.tif de UCSB)",
         "descargado": time.strftime("%Y-%m-%d")}, indent=2, ensure_ascii=False), encoding="utf-8")

da = xr.open_dataset(NC)["precipitation"].sel(time=slice(INI, FIN))
n_meses = len(pd.period_range(INI, FIN, freq="M"))
assert da.sizes["time"] == n_meses, f"CHIRPS: {da.sizes['time']} meses de {n_meses}"
assert da.notnull().all(), "CHIRPS: faltantes en el recorte"

# ------------------------------------------------------------------ pesos por área de intersección
basin = gpd.read_file(ROOT / "resultados" / "delimitacion" / "cuenca_el_tambo_glo30.gpkg").to_crs(32717)
h = 0.025
celdas = gpd.GeoDataFrame(
    [{"lat": la, "lon": lo, "geometry": box(lo - h, la - h, lo + h, la + h)}
     for la in da.lat.values for lo in da.lon.values], crs=4326).to_crs(32717)
celdas["area_en_cuenca_km2"] = celdas.intersection(basin.geometry.iloc[0]).area / 1e6
celdas["fraccion_celda"] = celdas["area_en_cuenca_km2"] / (celdas.area / 1e6)
pesos = celdas[celdas["area_en_cuenca_km2"] > 0].drop(columns="geometry")
pesos.round(5).to_csv(DL / "pesos_chirps_cuenca.csv", index=False)
W = xr.DataArray(celdas["area_en_cuenca_km2"].values.reshape(da.sizes["lat"], da.sizes["lon"]),
                 coords={"lat": da.lat, "lon": da.lon}, dims=("lat", "lon"))
serie = (da * W).sum(["lat", "lon"]) / W.sum()
pd.DataFrame({"P_chirps_mm_mes": serie.to_series().round(2)}).rename_axis("fecha").to_csv(
    SMES / "chirps_cuenca_mensual.csv")

# ------------------------------------------------------------------ resumen
def media_anual(s):
    return s.groupby(s.index.month).mean().sum()


cam = pd.read_csv(SMES / "series_mensuales_camels_pe.csv", index_col="fecha", parse_dates=True)["P_L_mm_mes"].loc[INI:FIN]
pv3 = pd.read_csv(SMES / "piscop_v3_cuenca_mensual.csv", index_col="fecha", parse_dates=True)["P_piscop_v3_mm_mes"].loc[INI:FIN]
s = serie.to_series()
print(f"Celdas en la cuenca: {len(pesos)} | completas (≥ 99 %): {(pesos.fraccion_celda >= 0.99).sum()} | "
      f"celdas equivalentes {pesos.fraccion_celda.sum():.1f}")
print(f"P media anual {INI:%Y-%m}..{FIN:%Y-%m} (Σ medias mensuales): CHIRPS {media_anual(s):.0f} | "
      f"PISCOp v3 {media_anual(pv3):.0f} | CAMELS-PE {media_anual(cam):.0f} mm/año")
print(f"r mensual CHIRPS vs CAMELS-PE: {np.corrcoef(s, cam)[0, 1]:.3f}")
