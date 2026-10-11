"""
Inversión térmica, subsidencia y convergencia de humedad sobre la cuenca El Tambo (ERA5 mensual).

Hipótesis a contrastar (1.5.c): el anticiclón del Pacífico Sur mantiene subsidencia de gran escala frente a la
costa; el aire que desciende se calienta por compresión y, sobre la capa marina enfriada por el afloramiento,
forma una inversión térmica que impide la convección. La subsidencia implica divergencia en niveles bajos
(continuidad de masa), así que la humedad que llega no se acumula. Si la hipótesis es correcta, la inversión
y la subsidencia deberían ser más intensas en may–sep (divergencia del IVT, estación seca) que en dic–abr
(convergencia, estación húmeda), y los años con inversión más débil deberían ser más lluviosos.

Datos (CDS, credenciales en ~/.cdsapirc, ver script 17):
  - reanalysis-era5-pressure-levels-monthly-means (doi:10.24381/cds.6860a573): temperatura y velocidad
    vertical ω (Pa/s; ω < 0 = ascenso) en 23 niveles de 1000 a 200 hPa.
  - reanalysis-era5-single-levels-monthly-means (doi:10.24381/cds.f17050d7): presión en superficie y
    temperatura a 2 m; geopotencial de superficie (relieve del modelo, un solo mes).
  Dominio 2°S–12°S × 88°O–72°O (Pacífico frente a la costa, cuenca y vertiente amazónica), malla 0.25°,
  años completos del periodo común (se recorta a periodo_comun.py).
Diagnósticos:
  (a) Corte vertical oeste–este promediado en 7°S–8°S (latitud de la cuenca): temperatura potencial θ,
      ω y capas de inversión (dT/dz ≥ 0, temperatura que aumenta con la altura), para dic–abr y may–sep.
  (b) Ciclo anual: estabilidad de la baja troposfera LTS = θ(700 hPa) − θ(superficie) (Klein & Hartmann 1993)
      sobre el océano costero (83°O–80°O × 6.5°S–9°S); ω en 500 hPa sobre la cuenca; convergencia del IVT
      sobre la cuenca (ERA5 del script 24, mm/d); P_L de CAMELS-PE.
  (c) Anomalías respecto al ciclo anual en dic–abr: LTS costera y ω500 sobre la cuenca frente a P_L.

Salidas:
  Tarea 1/El_Tambo/datos/era5_inversion/ (NetCDF, fuera de Git) y climatologia_inversion.nc (versionada)
  Tarea 1/El_Tambo/resultados/contexto_fisico/fig_inversion_termica_era5.png, ciclo_inversion_convergencia.csv,
  correlaciones_inversion_lluvia.csv
"""
from pathlib import Path
import json
import sys
import time
import zipfile

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from periodo_comun import periodo_comun

ROOT = Path(__file__).resolve().parents[2] / "Tarea 1" / "El_Tambo"
DL = ROOT / "datos" / "era5_inversion"
DL.mkdir(parents=True, exist_ok=True)
OUT = ROOT / "resultados" / "contexto_fisico"
INI, FIN = periodo_comun()
AREA = [-2, -88, -12, -72]                                   # N, O, S, E
NIVELES = [1000, 975, 950, 925, 900, 875, 850, 825, 800, 775, 750, 700, 650, 600, 550, 500, 450, 400, 350, 300,
           250, 225, 200]
BASE = {"product_type": ["monthly_averaged_reanalysis"], "time": ["00:00"], "area": AREA, "grid": [0.25, 0.25],
        "data_format": "netcdf", "download_format": "unarchived"}
ANIOS = [str(y) for y in range(INI.year, FIN.year + 1)]
MESES = [f"{m:02d}" for m in range(1, 13)]
PEDIDOS = {
    "niveles": ("reanalysis-era5-pressure-levels-monthly-means", "10.24381/cds.6860a573",
                {"variable": ["temperature", "vertical_velocity"], "pressure_level": [str(p) for p in NIVELES],
                 "year": ANIOS, "month": MESES}),
    "superficie": ("reanalysis-era5-single-levels-monthly-means", "10.24381/cds.f17050d7",
                   {"variable": ["surface_pressure", "2m_temperature"], "year": ANIOS, "month": MESES}),
    "relieve": ("reanalysis-era5-single-levels-monthly-means", "10.24381/cds.f17050d7",
                {"variable": ["geopotential"], "year": ["2000"], "month": ["01"]}),
}
CAJA_COSTA = dict(lon=slice(-83, -80), lat=slice(-9, -6.5))       # océano frente a la costa
CAJA_CUENCA = dict(lon=slice(-79, -78), lat=slice(-8, -7.25))     # celdas que cubren la cuenca
LAT_CORTE = slice(-8, -7)
HUMEDA, SECA = [12, 1, 2, 3, 4], [5, 6, 7, 8, 9]


ANIOS_POR_PEDIDO = 7          # el CDS limita el tamaño de cada solicitud: niveles × variables × meses


def _bajar(dataset, req, destino):
    import cdsapi
    tmp = destino.with_suffix(".tmp")
    cdsapi.Client().retrieve(dataset, req).download(str(tmp))
    if zipfile.is_zipfile(tmp):
        with zipfile.ZipFile(tmp) as z:
            miembros = [m for m in z.namelist() if m.endswith(".nc")]
            assert len(miembros) == 1, z.namelist()
            destino.write_bytes(z.read(miembros[0]))
        tmp.unlink()
    else:
        tmp.rename(destino)


def descargar(nombre):
    nc = DL / f"ERA5_{nombre}.nc"
    if nc.exists():
        return nc
    try:
        import cdsapi  # noqa: F401
    except ImportError:
        sys.exit("Falta la librería cdsapi (pip install cdsapi).")
    dataset, doi, extra = PEDIDOS[nombre]
    req = {**BASE, **extra}
    t0 = time.time()
    if nombre == "niveles":                                   # por bloques de años y luego se unen
        partes = []
        for i in range(0, len(ANIOS), ANIOS_POR_PEDIDO):
            bloque = ANIOS[i:i + ANIOS_POR_PEDIDO]
            parte = DL / f"ERA5_{nombre}_{bloque[0]}_{bloque[-1]}.part.nc"
            if not parte.exists():
                _bajar(dataset, {**req, "year": bloque}, parte)
                print(f"  {parte.name}: {parte.stat().st_size / 1e6:.0f} MB")
            partes.append(parte)
        ds_partes = [xr.load_dataset(p) for p in partes]               # en memoria: libera los archivos
        xr.concat(ds_partes, dim="valid_time").to_netcdf(nc)
        for p in partes:
            p.unlink()
    else:
        _bajar(dataset, req, nc)
    proc = DL / "PROCEDENCIA.json"
    meta = json.loads(proc.read_text(encoding="utf-8")) if proc.exists() else {}
    meta[nombre] = {"dataset_cds": dataset, "doi": doi, "solicitud": req, "descargado": time.strftime("%Y-%m-%d"),
                    "tamano_MB": round(nc.stat().st_size / 1e6, 1)}
    proc.write_text(json.dumps(meta, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"{nc.name}: {nc.stat().st_size / 1e6:.0f} MB en {time.time() - t0:.0f} s")
    return nc


def abrir(nc, recortar=True):
    ds = xr.open_dataset(nc)
    ren = {k: v for k, v in {"valid_time": "time", "latitude": "lat", "longitude": "lon",
                             "pressure_level": "nivel"}.items() if k in ds.dims or k in ds.coords}
    ds = ds.rename(ren).sortby("lat")
    ds = ds.drop_vars([v for v in ("number", "expver") if v in ds.variables])
    if "time" in ds.dims:
        ds["time"] = ds["time"].values.astype("datetime64[M]").astype("datetime64[ns]")
        if recortar:
            ds = ds.sel(time=slice(INI, FIN))
    return ds


niv = abrir(descargar("niveles"))
sup = abrir(descargar("superficie"))
zsf = abrir(descargar("relieve"), recortar=False)["z"].squeeze(drop=True) / 9.80665     # m
print(niv)
assert niv.sizes["time"] == sup.sizes["time"] == len(pd.period_range(INI, FIN, freq="M"))

# ------------------------------------------------------------------ variables derivadas
P = niv["nivel"]                                                     # hPa
theta = niv["t"] * (1000 / P) ** 0.2857                              # temperatura potencial (K)
bajo = sup["sp"] / 100 < P                                           # nivel bajo la superficie del modelo
theta_s = sup["t2m"] * (1000 / (sup["sp"] / 100)) ** 0.2857
LTS = theta.sel(nivel=700) - theta_s                                 # K


def media_caja(da, caja):
    d = da.sel(**caja)
    return d.weighted(np.cos(np.deg2rad(d.lat))).mean(["lat", "lon"])


lts_costa = media_caja(LTS, CAJA_COSTA).to_series()
w500 = media_caja(niv["w"].sel(nivel=500), CAJA_CUENCA).to_series()          # Pa/s

# Convergencia del IVT sobre la cuenca (ERA5 del script 24, 0.5°): −∇·IVT en mm/d
uni = abrir(ROOT / "datos" / "era5" / "ERA5_unico_mensual.nc")
R_T = 6.371e6
fl = np.deg2rad(uni.lat)
div = (uni["viwve"].differentiate("lon") + (uni["viwvn"] * np.cos(fl)).differentiate("lat")) \
      * np.rad2deg(1) / (R_T * np.cos(fl))
conv = media_caja(-div * 86400, dict(lon=slice(-79, -78), lat=slice(-8.25, -7.25))).to_series()

cam = pd.read_csv(ROOT / "resultados" / "series_mensuales" / "series_mensuales_camels_pe.csv", index_col="fecha",
                  parse_dates=True)["P_L_mm_mes"].loc[INI:FIN]
PLd = cam / cam.index.days_in_month
tabla = pd.DataFrame({"LTS_costa_K": lts_costa, "omega500_cuenca_Pa_s": w500, "conv_IVT_cuenca_mm_d": conv,
                      "P_L_mm_d": PLd})
ciclo = tabla.groupby(tabla.index.month).mean()
ciclo.index.name = "mes"
ciclo.round(4).to_csv(OUT / "ciclo_inversion_convergencia.csv")
print("\nCiclo anual:\n", ciclo.round(3).to_string())

# ------------------------------------------------------------------ anomalías (dic–abr)
an = tabla - tabla.groupby(tabla.index.month).transform("mean")
hum = an[an.index.month.isin(HUMEDA)].dropna()
filas = []
for x in ("LTS_costa_K", "omega500_cuenca_Pa_s", "conv_IVT_cuenca_mm_d"):
    r = hum[x].corr(hum["P_L_mm_d"])
    n = len(hum)
    # tamaño efectivo aproximado por autocorrelación de orden 1 de ambas series
    r1 = hum[x].autocorr() * hum["P_L_mm_d"].autocorr()
    n_ef = n * (1 - r1) / (1 + r1)
    t = r * np.sqrt((n_ef - 2) / (1 - r ** 2))
    filas.append({"predictor": x, "r_anomalias_dic_abr": round(r, 3), "n": n, "n_efectivo": round(n_ef, 1),
                  "t": round(t, 2)})
corr = pd.DataFrame(filas)
corr.to_csv(OUT / "correlaciones_inversion_lluvia.csv", index=False)
print("\nCorrelación de anomalías en dic–abr con P_L:\n", corr.to_string(index=False))

# ------------------------------------------------------------------ climatologías estacionales del corte
def estacion(da, meses):
    return da.sel(time=da.time.dt.month.isin(meses)).mean("time")


corte = {}
for nombre, meses in (("dic-abr", HUMEDA), ("may-sep", SECA)):
    th = estacion(theta, meses).sel(lat=LAT_CORTE).mean("lat")
    tt = estacion(niv["t"], meses).sel(lat=LAT_CORTE).mean("lat")
    om = estacion(niv["w"], meses).sel(lat=LAT_CORTE).mean("lat")
    spc = estacion(sup["sp"] / 100, meses).sel(lat=LAT_CORTE).mean("lat")
    corte[nombre] = xr.Dataset({"theta": th, "t": tt, "w": om, "sp": spc})
clim = xr.concat([c.expand_dims(estacion=[k]) for k, c in corte.items()], dim="estacion")
clim.attrs = {"fuente": "ERA5 monthly (C3S/ECMWF); script 29", "periodo": f"{INI:%Y-%m} a {FIN:%Y-%m}",
              "corte": "media 7°S–8°S"}
clim.to_netcdf(DL / "climatologia_inversion.nc")

# altura aproximada de cada nivel (atmósfera estándar) para el eje secundario
def p_a_z(p):
    return 44330.8 * (1 - (np.asarray(p) / 1013.25) ** 0.190263)


def z_a_p(z):
    return 1013.25 * (1 - np.asarray(z) / 44330.8) ** (1 / 0.190263)


# ------------------------------------------------------------------ figura
fig = plt.figure(figsize=(15, 15))
gs_top = fig.add_gridspec(1, 2, top=0.94, bottom=0.68, wspace=0.28)
gs = fig.add_gridspec(2, 2, top=0.52, bottom=0.06, hspace=0.38, wspace=0.28)
lon = clim.lon.values
pniv = clim.nivel.values
CUENCA_LON = (-78.75, -78.25)
ejes_corte = []
for j, (nombre, titulo) in enumerate((("dic-abr", "(a) Estación húmeda, dic–abr"),
                                      ("may-sep", "(b) Estación seca, may–sep"))):
    ax = fig.add_subplot(gs_top[0, j])
    ejes_corte.append(ax)
    c = clim.sel(estacion=nombre)
    mascara = c["sp"].values[None, :] < pniv[:, None]
    th = np.ma.masked_where(mascara, c["theta"].values)
    w = np.ma.masked_where(mascara, c["w"].values)
    cf = ax.contourf(lon, pniv, th, levels=np.arange(290, 350, 2), cmap="RdYlBu_r", extend="both")
    ax.contour(lon, pniv, th, levels=np.arange(290, 350, 2), colors="k", linewidths=0.25, alpha=0.5)
    # ω: contornos azules (ascenso, ω < 0) y rojos (descenso, ω > 0)
    ax.contour(lon, pniv, w, levels=[-0.08, -0.06, -0.04, -0.02], colors="#08306b", linewidths=1.3,
               linestyles="solid")
    ax.contour(lon, pniv, w, levels=[0.02, 0.04, 0.06, 0.08], colors="#a50f15", linewidths=1.3, linestyles="--")
    # capas de inversión: T aumenta con la altura (dT/dp > 0 hacia arriba ⇒ T(p_arriba) ≥ T(p_abajo))
    T = np.ma.masked_where(mascara, c["t"].values)
    inv = np.zeros_like(T, dtype=bool)
    inv[1:, :] = (T[1:, :] - T[:-1, :]) >= 0                         # nivel k (más alto) más cálido que k−1
    ax.contourf(lon, pniv, np.ma.masked_where(mascara, inv.astype(float)), levels=[0.5, 1.5], colors="none",
                hatches=["////"])
    ax.fill_between(lon, c["sp"].values, 1010, color="#5b4636", zorder=5)                  # relieve del modelo
    ax.fill_between(CUENCA_LON, z_a_p(684), z_a_p(4317), color="#2ca02c", alpha=0.45, zorder=6, lw=0)
    ax.annotate("cuenca", xy=(np.mean(CUENCA_LON), z_a_p(4317)), xytext=(-76.6, 380), fontsize=9,
                arrowprops=dict(arrowstyle="->", lw=0.8), zorder=7,
                bbox=dict(fc="white", alpha=0.8, lw=0, pad=1))
    ax.set_ylim(1010, 200); ax.set_xlim(lon.min(), lon.max())
    ax.set_yscale("log"); ax.set_yticks([1000, 850, 700, 500, 300, 200], ["1000", "850", "700", "500", "300", "200"])
    ax.minorticks_off()
    ax.set_ylabel("Presión (hPa)")
    ax.set_xticks(range(-88, -71, 4), [f"{abs(x)}° O" for x in range(-88, -71, 4)])
    sec = ax.secondary_yaxis("right", functions=(p_a_z, z_a_p))
    sec.set_yticks([0, 1000, 2000, 3000, 5000, 7000, 10000])
    sec.set_yticklabels(["0", "1", "2", "3", "5", "7", "10"]); sec.set_ylabel("Altura aprox. (km)")
    ax.set_title(f"{titulo} – corte a 7–8°S", loc="left", fontsize=11)
    ax.text(-87.6, 960, "Pacífico", fontsize=9, color="#08306b", zorder=7)
    ax.text(-74.5, 960, "Amazonía", fontsize=9, color="white", zorder=7)
pa, pb = ejes_corte[0].get_position(), ejes_corte[1].get_position()
cax = fig.add_axes([pa.x0, pa.y0 - 0.06, pb.x1 - pa.x0, 0.011])
cb = fig.colorbar(cf, cax=cax, orientation="horizontal")
cb.set_label("Temperatura potencial θ (K): las isolíneas se juntan donde el aire es más estable", fontsize=9)
leg = [Line2D([], [], color="#08306b", lw=1.3, label="ω < 0: ascenso (−0.02 a −0.08 Pa/s)"),
       Line2D([], [], color="#a50f15", lw=1.3, ls="--", label="ω > 0: descenso o subsidencia (0.02 a 0.08 Pa/s)"),
       Patch(fc="white", ec="k", hatch="////", label="Inversión: la temperatura aumenta con la altura"),
       Patch(fc="#5b4636", label="Relieve del modelo ERA5 (~30 km, suavizado)"),
       Patch(fc="#2ca02c", alpha=0.45, label="Cuenca: de 684 a 4 317 m (cotas reales)")]
fig.legend(handles=leg, loc="upper center", ncol=3, fontsize=8.5, frameon=False,
           bbox_to_anchor=(0.5, pa.y0 - 0.085))

# (c) ciclo anual
MES = list("EFMAMJJASOND")
ax1 = fig.add_subplot(gs[0, 0])
ax1.bar(range(12), ciclo["P_L_mm_d"], color="#1b6ca8", alpha=0.35, label="P$_L$ (mm/d)")
ax1.plot(range(12), ciclo["conv_IVT_cuenca_mm_d"], "-o", color="#0b7a75", label="Convergencia del IVT sobre la cuenca (mm/d)")
ax1.axhline(0, color="#555", lw=0.8)
ax1.set_xticks(range(12), MES); ax1.set_ylabel("mm/d")
ax1.set_title("(c) Lluvia y convergencia de humedad sobre la cuenca", loc="left", fontsize=11)
ax1.legend(fontsize=8, frameon=False, loc="upper right"); ax1.grid(alpha=0.3)
ax2 = fig.add_subplot(gs[0, 1])
ax2.plot(range(12), ciclo["LTS_costa_K"], "-o", color="#a50f15", label="LTS frente a la costa (K): intensidad de la inversión")
ax2.set_ylabel("LTS (K)", color="#a50f15"); ax2.tick_params(axis="y", colors="#a50f15")
ax2.set_xticks(range(12), MES); ax2.grid(alpha=0.3)
ax2b = ax2.twinx()
ax2b.plot(range(12), ciclo["omega500_cuenca_Pa_s"], "-s", color="#08306b", label="ω500 sobre la cuenca (Pa/s)")
ax2b.axhline(0, color="#08306b", lw=0.6, ls=":")
ax2b.set_ylabel("ω500 (Pa/s; > 0 = descenso)", color="#08306b"); ax2b.tick_params(axis="y", colors="#08306b")
ax2.set_title("(d) Inversión costera y movimiento vertical sobre la cuenca", loc="left", fontsize=11)
h1, l1 = ax2.get_legend_handles_labels(); h2, l2 = ax2b.get_legend_handles_labels()
ax2.legend(h1 + h2, l1 + l2, fontsize=8, frameon=True, framealpha=0.9, loc="lower center")

# (e) anomalías en dic–abr
for k, (x, etq, col) in enumerate((("LTS_costa_K", "Anomalía de LTS costera (K)", "#a50f15"),
                                   ("omega500_cuenca_Pa_s", "Anomalía de ω500 sobre la cuenca (Pa/s)", "#08306b"))):
    ax = fig.add_subplot(gs[1, k])
    ax.scatter(hum[x], hum["P_L_mm_d"], s=14, color=col, alpha=0.55)
    for a in (1998, 2017):
        s = hum[(hum.index.year == a) & hum.index.month.isin([2, 3, 4])]
        ax.scatter(s[x], s["P_L_mm_d"], s=40, facecolor="none", edgecolor="k", lw=1.2)
        for f_, r_ in s.iterrows():
            ax.annotate(f"{'FMA'[f_.month - 2]}-{str(a)[2:]}", (r_[x], r_["P_L_mm_d"]), fontsize=7.5,
                        xytext=(4, 3), textcoords="offset points")
    rr = corr.set_index("predictor").loc[x]
    ax.text(0.02, 0.96, f"r = {rr.r_anomalias_dic_abr:.2f} (n = {int(rr.n)}, n efectivo ≈ {rr.n_efectivo:.0f})",
            transform=ax.transAxes, va="top", fontsize=9)
    ax.axhline(0, color="#999", lw=0.6); ax.axvline(0, color="#999", lw=0.6)
    ax.set_xlabel(etq); ax.set_ylabel("Anomalía de P$_L$ (mm/d)")
    ax.set_title(f"({'ef'[k]}) Anomalías mensuales en dic–abr (círculos: feb–abr de 1998 y 2017)", loc="left",
                 fontsize=10.5)
    ax.grid(alpha=0.3)
fig.suptitle(f"Inversión térmica, subsidencia y convergencia de humedad – ERA5, periodo común {INI:%Y-%m} a {FIN:%Y-%m}",
             fontsize=13, y=0.995)
fig.text(0.5, 0.005, "LTS = θ(700 hPa) − θ(superficie), media 83–80°O × 9–6.5°S (Klein & Hartmann 1993). ω500 y convergencia "
         "del IVT: celdas sobre la cuenca. Corte: media 7–8°S; niveles bajo la superficie del modelo, en marrón.\n"
         "Fuente: ERA5 monthly averaged reanalysis (C3S/ECMWF; doi:10.24381/cds.6860a573, 10.24381/cds.f17050d7). "
         "Script 29_era5_inversion_termica.py.", ha="center", fontsize=8, color="#444")
fig.savefig(OUT / "fig_inversion_termica_era5.png", dpi=200, bbox_inches="tight")
plt.close(fig)
print(f"\nFigura: {OUT / 'fig_inversion_termica_era5.png'}")
