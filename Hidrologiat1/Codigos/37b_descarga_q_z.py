"""
Punto 5.1 (alternativa): campos de humedad específica y altura geopotencial para contrastar la elección de campos
atmosféricos del script 37a (presión al nivel del mar y viento zonal en 200 hPa).

Campos (ERA5 mensual, pedidos al CDS ya interpolados a 2°, 1981-01 a 2025-12, mismo formato que el 37a):
  · q   – humedad específica en 700 y 500 hPa (kg/kg → g/kg en el script 37). 700 hPa (~3 km) queda cerca de la cota
          mediana de la cuenca (2 578 m): es la humedad disponible a la altura de la cuenca. Sobre los Andes y el
          Altiplano ese nivel queda bajo el terreno del modelo; esas celdas se enmascaran con la presión en superficie.
  · z   – geopotencial en 500 y 700 hPa (m²/s²). En el script 37 se convierte a ALTURA geopotencial Z = z / g0
          (g0 = 9.80665 m/s²), en metros geopotenciales; no se confunden geopotencial y altura geopotencial.
  · sp  – presión en superficie (Pa), solo para la máscara de niveles bajo el terreno.
reanalysis-era5-pressure-levels-monthly-means (doi:10.24381/cds.6860a573) y
reanalysis-era5-single-levels-monthly-means (doi:10.24381/cds.f17050d7). Credenciales en ~/.cdsapirc (script 17).

Salidas (Tarea 1/El_Tambo/datos/campos_globales/):
  ERA5_q_z_500_700_mensual_2grados.nc, ERA5_sp_mensual_2grados.nc, PROCEDENCIA.json (se agrega la entrada)
"""
from pathlib import Path
import json
import sys
import time
import zipfile

import xarray as xr

ROOT = Path(__file__).resolve().parents[2] / "Tarea 1" / "El_Tambo"
DL = ROOT / "datos" / "campos_globales"
DL.mkdir(parents=True, exist_ok=True)
ANIOS = [str(y) for y in range(1981, 2026)]
PEDIDOS = {
    "q_z_500_700": ("reanalysis-era5-pressure-levels-monthly-means", "10.24381/cds.6860a573",
                    {"variable": ["specific_humidity", "geopotential"], "pressure_level": ["500", "700"]}),
    "sp": ("reanalysis-era5-single-levels-monthly-means", "10.24381/cds.f17050d7",
           {"variable": ["surface_pressure"]}),
}
PROC = DL / "PROCEDENCIA.json"


def descargar(nombre):
    nc = DL / f"ERA5_{nombre}_mensual_2grados.nc"
    if nc.exists():
        return nc
    if not (Path.home() / ".cdsapirc").exists():
        sys.exit(f"No existe {Path.home() / '.cdsapirc'} (ver cabecera del script 17).")
    import cdsapi
    dataset, doi, extra = PEDIDOS[nombre]
    req = {"product_type": ["monthly_averaged_reanalysis"], **extra, "year": ANIOS,
           "month": [f"{m:02d}" for m in range(1, 13)], "time": ["00:00"], "grid": [2.0, 2.0],
           "data_format": "netcdf", "download_format": "unarchived"}
    t0 = time.time()
    tmp = DL / f"descarga_{nombre}.tmp"
    cdsapi.Client().retrieve(dataset, req).download(str(tmp))
    if zipfile.is_zipfile(tmp):
        with zipfile.ZipFile(tmp) as z:
            miembros = [m for m in z.namelist() if m.endswith(".nc")]
            if len(miembros) != 1:
                sys.exit(f"ZIP inesperado: {z.namelist()}")
            nc.write_bytes(z.read(miembros[0]))
        tmp.unlink()
    else:
        tmp.replace(nc)
    meta = json.loads(PROC.read_text(encoding="utf-8")) if PROC.exists() else {}
    meta[nombre] = {"producto": "ERA5 monthly averaged reanalysis", "dataset_cds": dataset, "doi": doi,
                    "proveedor": "Copernicus Climate Change Service (C3S), ECMWF", "solicitud": req,
                    "script": Path(__file__).name, "tamano_MB": round(nc.stat().st_size / 1e6, 1),
                    "descargado": time.strftime("%Y-%m-%d")}
    PROC.write_text(json.dumps(meta, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"{nc.name}: {nc.stat().st_size / 1e6:.0f} MB en {time.time() - t0:.0f} s")
    return nc


if __name__ == "__main__":
    for n in PEDIDOS:
        print(xr.open_dataset(descargar(n)))
