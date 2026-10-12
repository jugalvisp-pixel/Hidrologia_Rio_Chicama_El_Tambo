"""
Punto 5.1: campos climáticos globales mensuales para los mapas de correlación.

Tres campos, malla común de 2° (moderada, adecuada para patrones de gran escala), 1981-01 a 2025-12:
  1. SST  – NOAA ERSST v5 (Huang et al. 2017, doi:10.7289/V5T72FNM), 2°, mensual, °C. Descarga directa del
            catálogo de NOAA PSL (archivo sst.mnmean.nc, sin credenciales). Tierra = NaN en el archivo original.
  2. SLP  – presión media al nivel del mar de ERA5 (msl, Pa → hPa), reanalysis-era5-single-levels-monthly-means
            (doi:10.24381/cds.f17050d7). Representa el anticiclón del Pacífico sur (subsidencia e inversión que
            frenan la lluvia en la estación seca, script 29) y el dipolo de presión del ENSO.
  3. U200 – componente zonal del viento en 200 hPa de ERA5 (u, m/s; positiva hacia el este),
            reanalysis-era5-pressure-levels-monthly-means (doi:10.24381/cds.6860a573). Los vientos del este
            en altura sobre los Andes acompañan el transporte de humedad desde la Amazonía (scripts 24–25). A
            200 hPa ninguna celda queda bajo el terreno del modelo, así que no hace falta máscara.
ERA5 se pide al CDS ya interpolado a 2° (parámetro "grid"); ERSST v5 ya está en 2°.

Acceso ERA5: cdsapi con ~/.cdsapirc (ver script 17). Las credenciales nunca se guardan en el repositorio.

Salidas (Tarea 1/El_Tambo/datos/campos_globales/):
  ERSST_v5_sst_mensual_1981_2025.nc, ERA5_msl_mensual_2grados.nc, ERA5_u200_mensual_2grados.nc, PROCEDENCIA.json
"""
from pathlib import Path
import json
import sys
import time
import urllib.request
import zipfile

import xarray as xr

ROOT = Path(__file__).resolve().parents[2] / "Tarea 1" / "El_Tambo"
DL = ROOT / "datos" / "campos_globales"
DL.mkdir(parents=True, exist_ok=True)
ANIOS = [str(y) for y in range(1981, 2026)]
GRID = [2.0, 2.0]
URL_ERSST = "https://downloads.psl.noaa.gov/Datasets/noaa.ersst.v5/sst.mnmean.nc"
PEDIDOS = {
    "msl": ("reanalysis-era5-single-levels-monthly-means", "10.24381/cds.f17050d7",
            {"variable": ["mean_sea_level_pressure"]}),
    "u200": ("reanalysis-era5-pressure-levels-monthly-means", "10.24381/cds.6860a573",
             {"variable": ["u_component_of_wind"], "pressure_level": ["200"]}),
}
PROC = DL / "PROCEDENCIA.json"


def registrar(clave, info):
    meta = json.loads(PROC.read_text(encoding="utf-8")) if PROC.exists() else {}
    meta[clave] = {**info, "descargado": time.strftime("%Y-%m-%d")}
    PROC.write_text(json.dumps(meta, indent=2, ensure_ascii=False), encoding="utf-8")


def ersst():
    nc = DL / "ERSST_v5_sst_mensual_1981_2025.nc"
    if nc.exists():
        return nc
    tmp = DL / "sst.mnmean.tmp"
    print(f"Descargando {URL_ERSST} …")
    urllib.request.urlretrieve(URL_ERSST, tmp)                  # una falla queda visible (excepción)
    ds = xr.open_dataset(tmp).sel(time=slice("1981-01-01", "2025-12-01")).load()
    ds.to_netcdf(nc)
    ds.close()
    tmp.unlink()
    registrar("sst", {"producto": "NOAA Extended Reconstructed SST v5 (ERSST v5)", "doi": "10.7289/V5T72FNM",
                      "url": URL_ERSST, "proveedor": "NOAA PSL", "variable": "sst (°C)", "malla": "2°",
                      "recorte": "1981-01 a 2025-12", "tamano_MB": round(nc.stat().st_size / 1e6, 1)})
    print(f"{nc.name}: {nc.stat().st_size / 1e6:.0f} MB")
    return nc


def era5(nombre):
    nc = DL / f"ERA5_{nombre}_mensual_2grados.nc"
    if nc.exists():
        return nc
    try:
        import cdsapi
    except ImportError:
        sys.exit("Falta la librería cdsapi (pip install cdsapi).")
    if not (Path.home() / ".cdsapirc").exists():
        sys.exit(f"No existe {Path.home() / '.cdsapirc'} (ver cabecera del script 17).")
    dataset, doi, extra = PEDIDOS[nombre]
    req = {"product_type": ["monthly_averaged_reanalysis"], **extra, "year": ANIOS,
           "month": [f"{m:02d}" for m in range(1, 13)], "time": ["00:00"], "grid": GRID,
           "data_format": "netcdf", "download_format": "unarchived"}
    t0 = time.time()
    tmp = DL / f"descarga_{nombre}.tmp"
    cdsapi.Client().retrieve(dataset, req).download(str(tmp))
    if zipfile.is_zipfile(tmp):                                  # por si el CDS entrega ZIP
        with zipfile.ZipFile(tmp) as z:
            miembros = [m for m in z.namelist() if m.endswith(".nc")]
            if len(miembros) != 1:
                sys.exit(f"ZIP inesperado: {z.namelist()}")
            nc.write_bytes(z.read(miembros[0]))
        tmp.unlink()
    else:
        tmp.rename(nc)
    registrar(nombre, {"producto": "ERA5 monthly averaged reanalysis", "dataset_cds": dataset, "doi": doi,
                       "proveedor": "Copernicus Climate Change Service (C3S), ECMWF", "solicitud": req,
                       "tamano_MB": round(nc.stat().st_size / 1e6, 1)})
    print(f"{nc.name}: {nc.stat().st_size / 1e6:.0f} MB en {time.time() - t0:.0f} s")
    return nc


if __name__ == "__main__":
    for f in (ersst(), era5("msl"), era5("u200")):
        print(xr.open_dataset(f))
