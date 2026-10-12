"""
Descarga de perfiles verticales para analizar la inversión térmica sobre la cuenca El Tambo.

1. Radiosondeos observados de Trujillo (IGRA2, NOAA NCEI)
   Estación PEM00084516 «RS TRUJILLO» (WMO 84516), −8.0439°, −79.0593°, 127 m s.n.m., a ~40 km al SO de la
   salida de la cuenca. Lanzamientos a las 12 UTC (7 a.m. hora local), desde 2024-10.
   Durre, I., Yin, X., Vose, R. S., Applequist, S., Arnfield, J. (2018). Enhancing the Data Coverage in the
   Integrated Global Radiosonde Archive. J. Atmos. Oceanic Technol., 35, 1753–1770. doi:10.1175/JTECH-D-17-0223.1
   Archivo: https://www.ncei.noaa.gov/data/integrated-global-radiosonde-archive/access/data-por/

2. ERA5 mensual por niveles de presión, promedio de las 12 UTC (misma hora que los sondeos)
   reanalysis-era5-pressure-levels-monthly-means, product_type monthly_averaged_reanalysis_by_hour_of_day,
   doi:10.24381/cds.6860a573. Temperatura, humedad relativa y geopotencial en 1000–300 hPa, malla 0.25°,
   recuadro 7.0°S–8.3°S × 79.5°O–78.0°O (costa de Trujillo + cuenca), periodo común.
   Requiere cuenta en el CDS (~/.cdsapirc, fuera del repositorio).

Salidas: Tarea 1/El_Tambo/datos/radiosondeo/  y  Tarea 1/El_Tambo/datos/era5_perfiles/
"""
from pathlib import Path
import json
import sys
import time
import urllib.request
import zipfile

from periodo_comun import periodo_comun

ROOT = Path(__file__).resolve().parents[2] / "Tarea 1" / "El_Tambo" / "datos"
RS = ROOT / "radiosondeo"
E5 = ROOT / "era5_perfiles"
RS.mkdir(parents=True, exist_ok=True)
E5.mkdir(parents=True, exist_ok=True)
INI, FIN = periodo_comun()

ESTACION = "PEM00084516"
IGRA_URL = f"https://www.ncei.noaa.gov/data/integrated-global-radiosonde-archive/access/data-por/{ESTACION}-data.txt.zip"
NIVELES = ["300", "350", "400", "450", "500", "550", "600", "650", "700", "750", "775", "800", "825", "850", "875",
           "900", "925", "950", "975", "1000"]
AREA = [-7.0, -79.5, -8.3, -78.0]                            # N, O, S, E

# ------------------------------------------------------------------ 1. IGRA2
zp = RS / f"{ESTACION}-data.txt.zip"
urllib.request.urlretrieve(IGRA_URL, zp)                    # siempre la versión vigente (se actualiza a diario)
with zipfile.ZipFile(zp) as z:
    z.extractall(RS)
txt = RS / f"{ESTACION}-data.txt"
n = sum(1 for linea in txt.open(encoding="ascii") if linea.startswith("#"))
(RS / "PROCEDENCIA.json").write_text(json.dumps({
    "producto": "Integrated Global Radiosonde Archive v2 (IGRA2), NOAA NCEI", "estacion": ESTACION,
    "nombre": "RS TRUJILLO (WMO 84516)", "lat": -8.0439, "lon": -79.0593, "elevacion_m": 127,
    "url": IGRA_URL, "referencia": "doi:10.1175/JTECH-D-17-0223.1", "sondeos": n,
    "descargado": time.strftime("%Y-%m-%d")}, indent=2, ensure_ascii=False), encoding="utf-8")
print(f"IGRA2 {ESTACION}: {n} sondeos")

# ------------------------------------------------------------------ 2. ERA5 12 UTC
nc = E5 / "ERA5_perfiles_12UTC_mensual.nc"
if not nc.exists():
    if not (Path.home() / ".cdsapirc").exists():
        sys.exit(f"No existe {Path.home() / '.cdsapirc'} (ver cabecera del script 17).")
    import cdsapi
    req = {"product_type": ["monthly_averaged_reanalysis_by_hour_of_day"],
           "variable": ["temperature", "relative_humidity", "geopotential"], "pressure_level": NIVELES,
           "year": [str(y) for y in range(INI.year, FIN.year + 1)], "month": [f"{m:02d}" for m in range(1, 13)],
           "time": ["12:00"], "area": AREA, "data_format": "netcdf", "download_format": "unarchived"}
    t0 = time.time()
    tmp = E5 / "descarga.tmp"
    cdsapi.Client().retrieve("reanalysis-era5-pressure-levels-monthly-means", req).download(str(tmp))
    if zipfile.is_zipfile(tmp):
        with zipfile.ZipFile(tmp) as z:
            miembros = [m for m in z.namelist() if m.endswith(".nc")]
            nc.write_bytes(z.read(miembros[0]))
        tmp.unlink()
    else:
        tmp.rename(nc)
    (E5 / "PROCEDENCIA.json").write_text(json.dumps({
        "producto": "ERA5 monthly averaged reanalysis by hour of day on pressure levels",
        "dataset_cds": "reanalysis-era5-pressure-levels-monthly-means", "doi": "10.24381/cds.6860a573",
        "proveedor": "Copernicus Climate Change Service (C3S), ECMWF", "solicitud": req,
        "descargado": time.strftime("%Y-%m-%d")}, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"{nc.name}: {nc.stat().st_size / 1e6:.1f} MB en {time.time() - t0:.0f} s")
else:
    print(f"{nc.name} ya existe")
