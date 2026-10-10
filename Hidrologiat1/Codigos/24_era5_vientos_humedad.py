"""
Viento y transporte de humedad regional (ERA5 mensual) para el mapa de ubicación en Sudamérica.

Producto : ERA5 monthly averaged reanalysis (Copernicus C3S / ECMWF), dos conjuntos del CDS:
  - reanalysis-era5-pressure-levels-monthly-means (doi:10.24381/cds.6860a573):
      u, v (componentes del viento, m/s) y q (humedad específica, kg/kg) en 850, 700, 500 y 200 hPa.
  - reanalysis-era5-single-levels-monthly-means (doi:10.24381/cds.f17050d7):
      componentes este y norte del flujo de vapor de agua integrado en la vertical (IVT, kg m⁻¹ s⁻¹) y
      agua precipitable (tcwv, kg/m²).
Malla 0.5° (se pide al CDS interpolada desde 0.25°), años completos que cubren el periodo común
(periodo_comun.py); el análisis se limita a ese periodo.
Dominio: 20°N–40°S, 110°O–20°O. Cubre las regiones que pueden aportar humedad a la cuenca, que el
grupo debe sustentar con bibliografía (1.5.c):
  (1) Atlántico tropical: alisios que llevan humedad al continente;
  (2) cuenca amazónica: flujo del este, reciclaje y chorro de bajo nivel al este de los Andes;
  (3) Pacífico oriental tropical y costa peruana: Niño 1+2 (0–10°S, 90–80°O) y parte de Niño 3.
Niveles: 850 hPa (capa límite y chorro de bajo nivel), 700 hPa (cerca de la elevación media de la cuenca,
2557 m), 500 hPa (flujo medio sobre los Andes) y 200 hPa (alta troposfera, Alta de Bolivia).

Acceso: cdsapi con ~/.cdsapirc (ver script 17); hay que aceptar la licencia de ERA5 en el portal del CDS.
Los NetCDF mensuales completos (~cientos de MB) quedan fuera del repositorio; el script guarda además
una climatología estacional pequeña (1°) que sí se versiona y es la que usa el mapa del script 06.

Salidas:
  Tarea 1/El_Tambo/datos/era5/ERA5_niveles_mensual.nc, ERA5_unico_mensual.nc, PROCEDENCIA.json
  Tarea 1/El_Tambo/datos/era5/era5_viento_ivt_climatologia_estacional.nc
    (medias del periodo común para dic–abr, may–sep, oct–nov y anual)
"""
from pathlib import Path
import json
import sys
import time
import zipfile

import numpy as np
import xarray as xr

from periodo_comun import periodo_comun

ROOT = Path(__file__).resolve().parents[2] / "Tarea 1" / "El_Tambo"
DL = ROOT / "datos" / "era5"
DL.mkdir(parents=True, exist_ok=True)
INI, FIN = periodo_comun()
AREA = [20, -110, -40, -20]                                 # N, O, S, E
GRID = [0.5, 0.5]
NIVELES = ["200", "500", "700", "850"]
ESTACIONES = {"dic-abr": [12, 1, 2, 3, 4], "may-sep": [5, 6, 7, 8, 9], "oct-nov": [10, 11],
              "anual": list(range(1, 13))}
PEDIDOS = {
    "niveles": ("reanalysis-era5-pressure-levels-monthly-means", "10.24381/cds.6860a573",
                {"variable": ["u_component_of_wind", "v_component_of_wind", "specific_humidity"],
                 "pressure_level": NIVELES}),
    "unico": ("reanalysis-era5-single-levels-monthly-means", "10.24381/cds.f17050d7",
              {"variable": ["vertical_integral_of_eastward_water_vapour_flux",
                            "vertical_integral_of_northward_water_vapour_flux",
                            "total_column_water_vapour"]}),
    "superficie": ("reanalysis-era5-single-levels-monthly-means", "10.24381/cds.f17050d7",
                   {"variable": ["surface_pressure"]}),
}


def descargar(nombre):
    nc = DL / f"ERA5_{nombre}_mensual.nc"
    if nc.exists():
        return nc
    try:
        import cdsapi
    except ImportError:
        sys.exit("Falta la librería cdsapi (pip install cdsapi).")
    if not (Path.home() / ".cdsapirc").exists():
        sys.exit(f"No existe {Path.home() / '.cdsapirc'} (ver cabecera del script 17).")
    dataset, doi, extra = PEDIDOS[nombre]
    req = {"product_type": ["monthly_averaged_reanalysis"], **extra,
           "year": [str(y) for y in range(INI.year, FIN.year + 1)], "month": [f"{m:02d}" for m in range(1, 13)],
           "time": ["00:00"], "area": AREA, "grid": GRID, "data_format": "netcdf", "download_format": "unarchived"}
    t0 = time.time()
    tmp = DL / f"descarga_{nombre}.tmp"
    cdsapi.Client().retrieve(dataset, req).download(str(tmp))
    if zipfile.is_zipfile(tmp):                              # por si el CDS entrega ZIP
        with zipfile.ZipFile(tmp) as z:
            miembros = [m for m in z.namelist() if m.endswith(".nc")]
            if len(miembros) != 1:
                sys.exit(f"ZIP inesperado: {z.namelist()}")
            nc.write_bytes(z.read(miembros[0]))
        tmp.unlink()
    else:
        tmp.rename(nc)
    proc = DL / "PROCEDENCIA.json"
    meta = json.loads(proc.read_text(encoding="utf-8")) if proc.exists() else {}
    meta[nombre] = {"producto": "ERA5 monthly averaged reanalysis", "dataset_cds": dataset, "doi": doi,
                    "proveedor": "Copernicus Climate Change Service (C3S), ECMWF", "solicitud": req,
                    "descargado": time.strftime("%Y-%m-%d"), "tamano_MB": round(nc.stat().st_size / 1e6, 1)}
    proc.write_text(json.dumps(meta, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"{nc.name}: {nc.stat().st_size / 1e6:.0f} MB en {time.time() - t0:.0f} s")
    return nc


def abrir(nc):
    ds = xr.open_dataset(nc)
    ren = {k: v for k, v in {"valid_time": "time", "latitude": "lat", "longitude": "lon",
                             "pressure_level": "nivel"}.items() if k in ds.dims or k in ds.coords}
    ds = ds.rename(ren).sortby("lat")
    ds = ds.drop_vars([v for v in ("number", "expver") if v in ds.variables])
    ds["time"] = ds["time"].values.astype("datetime64[M]").astype("datetime64[ns]")
    return ds.sel(time=slice(INI, FIN))


niv = abrir(descargar("niveles"))
uni = abrir(descargar("unico"))
sup = abrir(descargar("superficie"))
print(niv)
print(uni)
n_meses = niv.sizes["time"]
assert n_meses == uni.sizes["time"], "los dos conjuntos no cubren los mismos meses"
print(f"Meses en el periodo común: {n_meses} ({INI:%Y-%m} a {FIN:%Y-%m})")

# ------------------------------------------------------------------ climatología estacional (1°)
iv_e = [v for v in uni.data_vars if "eastward" in uni[v].attrs.get("long_name", "").lower() or v.startswith("p71")][0]
iv_n = [v for v in uni.data_vars if "northward" in uni[v].attrs.get("long_name", "").lower() or v.startswith("p72")][0]
out = xr.Dataset({
    "u": niv["u"], "v": niv["v"], "q": niv["q"] * 1000,
    "ivt_e": uni[iv_e], "ivt_n": uni[iv_n], "tcwv": uni["tcwv"], "sp": sup["sp"] / 100})
out = out.coarsen(lat=2, lon=2, boundary="trim").mean()      # 0.5° → 1° (suficiente para el mapa)
clim = xr.concat([out.sel(time=out.time.dt.month.isin(m)).mean("time").expand_dims(estacion=[e])
                  for e, m in ESTACIONES.items()], dim="estacion")
# Nivel bajo el terreno del modelo: presión en superficie menor que la del nivel → valores extrapolados.
clim["bajo_terreno"] = (clim["sp"] < clim["nivel"]).astype("int8")
for k in ("u", "v", "q"):
    clim[k] = clim[k].where(clim["bajo_terreno"] == 0)
clim["vel"] = np.hypot(clim["u"], clim["v"])
clim["ivt"] = np.hypot(clim["ivt_e"], clim["ivt_n"])
unidades = {"u": "m/s", "v": "m/s", "vel": "m/s", "q": "g/kg", "ivt_e": "kg m-1 s-1", "ivt_n": "kg m-1 s-1",
            "ivt": "kg m-1 s-1", "tcwv": "kg m-2", "sp": "hPa"}
for k, u in unidades.items():
    clim[k].attrs["units"] = u
clim.attrs = {"fuente": "ERA5 monthly averaged reanalysis (C3S/ECMWF); script 24_era5_vientos_humedad.py",
              "periodo": f"{INI:%Y-%m} a {FIN:%Y-%m} ({n_meses} meses)",
              "estaciones": json.dumps(ESTACIONES), "malla": "1° (promedio de 2×2 celdas de 0.5°)"}
clim["bajo_terreno"].attrs["description"] = "1 = el nivel queda bajo la superficie del modelo (u, v, q enmascarados)"
for k in clim.data_vars:
    if k != "bajo_terreno":
        clim[k] = clim[k].astype("float32")
f = DL / "era5_viento_ivt_climatologia_estacional.nc"
clim.to_netcdf(f, encoding={v: {"zlib": True, "complevel": 5} for v in clim.data_vars})
print(f"{f.name}: {f.stat().st_size / 1e6:.1f} MB")

# Resumen en el entorno de la cuenca (−7.6°, −78.5°), para la descripción del mapa
pt = clim.sel(lat=-7.5, lon=-78.5, method="nearest")
for e in ESTACIONES:
    p = pt.sel(estacion=e)
    txt = "; ".join(f"{int(n)} hPa: {float(p.vel.sel(nivel=int(n))):.1f} m/s desde "
                    f"{(np.degrees(np.arctan2(-float(p.u.sel(nivel=int(n))), -float(p.v.sel(nivel=int(n))))) % 360):.0f}°"
                    for n in NIVELES)
    print(f"{e:8s} IVT {float(p.ivt):.0f} kg/m/s | {txt}")
