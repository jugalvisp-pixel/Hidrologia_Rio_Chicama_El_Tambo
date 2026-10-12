"""
Contraste de los radiosondeos de Trujillo con ERA5, en dos comparaciones separadas.

Los sondeos de Trujillo (IGRA2 84516) solo existen desde 2024-10, casi todos fuera del periodo común
(1998-01 a 2024-12). Para saber si la inversión que muestran representa el periodo común se separan dos preguntas:

  Comparación 1 · Trujillo (sondeos) frente a ERA5 en los MISMOS meses (2024-10 a 2026-09):
      ¿reproduce ERA5 el perfil observado? Mide el error del reanálisis, sin diferencia de periodo.
  Comparación 2 · ERA5 en los meses de los sondeos (2024-10 a 2026-09) frente a ERA5 en el periodo común:
      ¿fueron esos dos años normales? Misma fuente, solo cambia el periodo.

Datos:
  - Perfiles de cada sondeo en una malla de 25 m (script 34b, perfiles_por_sondeo.nc), 12 UTC.
  - ERA5 mensual por niveles de presión, promedio de las 12 UTC (reanalysis-era5-pressure-levels-monthly-means,
    monthly_averaged_reanalysis_by_hour_of_day; doi:10.24381/cds.6860a573), celda de 0.25° más cercana a
    Trujillo (−8.0°, −79.0°). Periodo común: archivo del script 34a. 2025-01 a 2026-09: se descarga aquí
    (mismo recuadro y niveles que 34a); 2024-10 a 2024-12 se toman del archivo del 34a.

Método: temporadas dic–abr (húmeda) y may–sep (seca), las de los scripts 34b y 35. ERA5 es un promedio
mensual, así que de los sondeos se usa la MEDIA por nivel (no la mediana) para comparar lo mismo; la banda es el
rango intercuartílico de los sondeos. El perfil de ERA5 se dibuja en la altura geopotencial media de cada nivel
y se interpola a la malla de 25 m para las diferencias. Estabilidad de la capa baja: T(1 500 m) − T(500 m);
positiva = inversión neta en esa capa.

Salidas (Tarea 1/El_Tambo/resultados/inversion/):
  fig_contraste_1_trujillo_vs_era5.png, fig_contraste_2_era5_periodos.png, contraste_sondeos_era5.csv
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
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

from periodo_comun import periodo_comun
from estilo_figuras import EST

ROOT = Path(__file__).resolve().parents[2] / "Tarea 1" / "El_Tambo"
E5 = ROOT / "datos" / "era5_perfiles"
NC_COMUN = E5 / "ERA5_perfiles_12UTC_mensual.nc"                       # script 34a
NC_SOND = E5 / "ERA5_perfiles_12UTC_mensual_2025_2026.nc"               # este script
INV = ROOT / "resultados" / "inversion"
MORF = ROOT / "resultados" / "delimitacion" / "parametros_morfometricos.csv"
INI, FIN = periodo_comun()
S_INI, S_FIN = pd.Timestamp("2024-10-01"), pd.Timestamp("2026-09-30")  # meses completos con sondeos
G = 9.80665
NIVELES = ["300", "350", "400", "450", "500", "550", "600", "650", "700", "750", "775", "800", "825", "850", "875",
           "900", "925", "950", "975", "1000"]                         # los de 34a
AREA = [-7.0, -79.5, -8.3, -78.0]
TEMP = {"may–sep (seca)": ([5, 6, 7, 8, 9], EST["seca"]), "dic–abr (húmeda)": ([12, 1, 2, 3, 4], EST["húmeda"])}
Z_TOP = 5000

# ------------------------------------------------------------------ 1. ERA5 2025-01 a 2026-09
if not NC_SOND.exists():
    if not (Path.home() / ".cdsapirc").exists():
        sys.exit(f"No existe {Path.home() / '.cdsapirc'} (ver cabecera del script 17).")
    import cdsapi
    req = {"product_type": ["monthly_averaged_reanalysis_by_hour_of_day"],
           "variable": ["temperature", "relative_humidity", "geopotential"], "pressure_level": NIVELES,
           "year": ["2025", "2026"], "month": [f"{m:02d}" for m in range(1, 13)], "time": ["12:00"], "area": AREA,
           "data_format": "netcdf", "download_format": "unarchived"}
    tmp = E5 / "descarga_2025_2026.tmp"
    cdsapi.Client().retrieve("reanalysis-era5-pressure-levels-monthly-means", req).download(str(tmp))
    if zipfile.is_zipfile(tmp):
        with zipfile.ZipFile(tmp) as z:
            NC_SOND.write_bytes(z.read([m for m in z.namelist() if m.endswith(".nc")][0]))
        tmp.unlink()
    else:
        tmp.replace(NC_SOND)
    proc = E5 / "PROCEDENCIA.json"
    p = json.loads(proc.read_text(encoding="utf-8")) if proc.exists() else {}
    p["extension_2025_2026"] = {"archivo": NC_SOND.name, "solicitud": req, "script": Path(__file__).name,
                                "descargado": time.strftime("%Y-%m-%d")}
    proc.write_text(json.dumps(p, indent=2, ensure_ascii=False), encoding="utf-8")


def abrir(nc):
    ds = xr.open_dataset(nc)
    ds = ds.rename({k: v for k, v in {"valid_time": "time", "latitude": "lat", "longitude": "lon",
                                       "pressure_level": "level"}.items() if k in ds.dims or k in ds.coords})
    return ds.sel(lat=-8.0, lon=-79.0, method="nearest")[["t", "r", "z"]].load()


pc, ps = abrir(NC_COMUN), abrir(NC_SOND)
CELDA = (float(pc.lat), float(pc.lon))
e5_comun = pc.sel(time=slice(INI, FIN + pd.offsets.MonthEnd(0) + pd.Timedelta(hours=23)))
e5_sond = xr.concat([pc.sel(time=slice(S_INI, "2024-12-31T23")), ps.sel(time=slice("2025-01-01", S_FIN + pd.Timedelta(hours=23)))],
                    "time")
n_esp = len(pd.date_range(S_INI, S_FIN, freq="MS"))
if e5_sond.time.size != n_esp:
    sys.exit(f"ERA5 2024-10 a 2026-09: {e5_sond.time.size} meses de {n_esp}; ¿falta publicar algún mes en el CDS?")
print(f"ERA5 celda {CELDA}: periodo común {e5_comun.time.size} meses; meses de los sondeos {e5_sond.time.size}")

# ------------------------------------------------------------------ 2. perfiles por temporada
snd = xr.open_dataset(INV / "perfiles_por_sondeo.nc").load()
snd = snd.sel(fecha=slice(S_INI, S_FIN + pd.Timedelta(hours=23)))
MALLA = snd.z_m.values


def era5_temporada(ds, meses):
    s = ds.sel(time=ds.time.dt.month.isin(meses))
    m = s.mean("time")
    z = (m["z"] / G).values
    o = np.argsort(z)
    return dict(z=z[o], T=(m["t"] - 273.15).values[o], HR=m["r"].values[o], n=int(s.time.size))


def sondeos_temporada(meses):
    s = snd.sel(fecha=snd.fecha.dt.month.isin(meses))
    ok = np.isfinite(s.T_C).sum("fecha").values >= max(10, 0.5 * s.fecha.size)
    out = {"n": int(s.fecha.size)}
    for v in ("T_C", "HR"):
        a = s[v].values
        out[v[:2].rstrip("_")] = np.where(ok, np.nanmean(a, 0), np.nan)
        out[v[:2].rstrip("_") + "_q1"] = np.where(ok, np.nanpercentile(a, 25, 0), np.nan)
        out[v[:2].rstrip("_") + "_q3"] = np.where(ok, np.nanpercentile(a, 75, 0), np.nan)
    return out


def en_malla(p, v):
    return np.interp(MALLA, p["z"], p[v], left=np.nan, right=np.nan)


def estab(z, t):
    return np.interp(1500, z, t) - np.interp(500, z, t)


S, EC, ES = {}, {}, {}
for t, (meses, _) in TEMP.items():
    S[t] = sondeos_temporada(meses)
    EC[t] = era5_temporada(e5_comun, meses)
    ES[t] = era5_temporada(e5_sond, meses)

# tabla resumen
filas = []
for t in TEMP:
    for fuente, z, T, HR, n in (("Trujillo, sondeos IGRA2 (media) 2024-10 a 2026-09", MALLA, S[t]["T"], S[t]["HR"], S[t]["n"]),
                                ("ERA5 12 UTC 2024-10 a 2026-09", ES[t]["z"], ES[t]["T"], ES[t]["HR"], ES[t]["n"]),
                                (f"ERA5 12 UTC {INI:%Y-%m} a {FIN:%Y-%m} (periodo común)", EC[t]["z"], EC[t]["T"],
                                 EC[t]["HR"], EC[t]["n"])):
        ok = np.isfinite(T)
        filas.append({"temporada": t, "fuente": fuente, "n": n,
                      "T_500m_C": np.interp(500, z[ok], T[ok]), "T_1500m_C": np.interp(1500, z[ok], T[ok]),
                      "estabilidad_T1500_menos_T500_C": estab(z[ok], T[ok]),
                      "HR_500m": np.interp(500, z[ok], HR[ok]), "HR_2000m": np.interp(2000, z[ok], HR[ok])})
tab = pd.DataFrame(filas)
tab.round(2).to_csv(INV / "contraste_sondeos_era5.csv", index=False)
print(tab.round(1).to_string())

# ------------------------------------------------------------------ 3. figuras
morf = pd.read_csv(MORF, index_col=0)["valor"]
Z_REF = [(float(morf["Cota cierre (m)"]), "El Tambo (cierre)"),
         (float(morf["Cota mediana (m s.n.m.)"]), "cota mediana de la cuenca"),
         (float(morf["Cota máxima (m s.n.m.)"]), "cota máxima")]
plt.rcParams.update({"font.size": 9.5, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.edgecolor": "#898781", "axes.grid": True, "grid.color": "#e1e0d9", "grid.linewidth": 0.6})


def referencias(ax, x, ha):
    for z, txt in Z_REF:
        ax.axhline(z, color="#52514e", lw=0.8, ls=(0, (1, 2)), zorder=1)
        ax.text(x, z + 40, f"{txt} · {z:,.0f} m".replace(",", " "), fontsize=7.3, color="#52514e", ha=ha)


def figura(titulo, sub, curvas, etq_dif, nombre, nota):
    """curvas[t] = (A, B): A continua, B discontinua; cada una dict(z, T, HR, [T_q1, T_q3, HR_q1, HR_q3]).
    Panel (c): B − A en la malla de 25 m."""
    fig, axs = plt.subplots(1, 3, figsize=(13.5, 7.2), sharey=True, gridspec_kw={"width_ratios": [1.2, 1, 0.9]})
    for t, (meses, col) in TEMP.items():
        A, B = curvas[t]
        for ax, v in ((axs[0], "T"), (axs[1], "HR")):
            if f"{v}_q1" in A:
                ax.fill_betweenx(A["z"], A[f"{v}_q1"], A[f"{v}_q3"], color=col, alpha=0.15, lw=0)
            ax.plot(A[v], A["z"], color=col, lw=2.4)
            ax.plot(B[v], B["z"], color=col, lw=1.4, ls="--", marker="o", ms=3.5)
        dif = np.interp(MALLA, B["z"], B["T"], left=np.nan, right=np.nan) - np.interp(
            MALLA, A["z"], A["T"], left=np.nan, right=np.nan)
        axs[2].plot(dif, MALLA, color=col, lw=2)
    axs[2].axvline(0, color="#0b0b0b", lw=0.8)
    axs[0].set(xlim=(-2, 26), ylim=(0, Z_TOP), xlabel="Temperatura (°C)", ylabel="Altitud (m s.n.m.)")
    axs[1].set(xlim=(0, 102), xlabel="Humedad relativa (%)")
    axs[2].set(xlim=(-4, 4), xlabel=f"Diferencia de temperatura (°C)\n{etq_dif}")
    axs[2].xaxis.label.set_size(8.5)
    referencias(axs[0], -1.6, "left")
    referencias(axs[1], 2, "left")
    for z, _ in Z_REF:
        axs[2].axhline(z, color="#52514e", lw=0.8, ls=(0, (1, 2)))
    axs[0].set_title("(a) Temperatura", loc="left", fontsize=10, fontweight="bold")
    axs[1].set_title("(b) Humedad relativa", loc="left", fontsize=10, fontweight="bold")
    axs[2].set_title("(c) Diferencia de temperatura", loc="left", fontsize=10, fontweight="bold")
    h = [Line2D([], [], color="#0b0b0b", lw=2.4), Line2D([], [], color="#0b0b0b", lw=1.4, ls="--", marker="o", ms=3.5)]
    l = list(sub)
    if "T_q1" in curvas["may–sep (seca)"][0]:
        h.append(Patch(color="#898781", alpha=0.3)); l.append("Rango intercuartílico de los sondeos")
    h += [Patch(color=c) for _, c in TEMP.values()]
    l += [f"{t} (n = {curvas[t][0].get('n', '')} / {curvas[t][1].get('n', '')})" for t in TEMP]
    axs[0].legend(h, l, loc="lower left", fontsize=7.6, framealpha=0.94)
    fig.suptitle(titulo, fontsize=11.5, fontweight="bold", x=0.02, ha="left")
    fig.text(0.02, 0.008, nota, fontsize=7.5, color="#52514e")
    fig.subplots_adjust(left=0.06, right=0.99, top=0.89, bottom=0.17, wspace=0.08)
    fig.savefig(INV / nombre, dpi=200)
    plt.close(fig)
    print("Figura:", INV / nombre)


def sond_dict(t):
    s = S[t]
    return dict(z=MALLA, T=s["T"], HR=s["HR"], T_q1=s["T_q1"], T_q3=s["T_q3"], HR_q1=s["HR_q1"], HR_q3=s["HR_q3"],
                n=s["n"])


PS = f"{S_INI:%Y-%m} a {S_FIN:%Y-%m}"
PC = f"{INI:%Y-%m} a {FIN:%Y-%m}"
figura(f"Comparación 1 · Trujillo (sondeos IGRA2) frente a ERA5 en los mismos meses, {PS}",
       [f"Trujillo: media de los sondeos 12 UTC ({PS})", f"ERA5: promedio de las 12 UTC ({PS})"],
       {t: (sond_dict(t), ES[t]) for t in TEMP}, "ERA5 − Trujillo\n> 0: ERA5 más cálido",
       "fig_contraste_1_trujillo_vs_era5.png",
       f"Misma hora (12 UTC = 7 a.m.) y mismos meses: la diferencia es error de ERA5, no del periodo. n = sondeos / meses de "
       f"ERA5. ERA5: celda {CELDA[0]:.2f}°, {CELDA[1]:.2f}° (0.25°), niveles de\npresión dibujados en su altura geopotencial "
       f"media. Sondeos: IGRA2 84516 (Trujillo, 127 m s.n.m.) interpolados cada 25 m. Scripts 34a, 34b y 34c.")
figura(f"Comparación 2 · ERA5: meses de los sondeos ({PS}) frente al periodo común ({PC})",
       [f"ERA5 periodo común ({PC})", f"ERA5 meses de los sondeos ({PS})"],
       {t: (EC[t], ES[t]) for t in TEMP}, f"({PS}) − ({PC})\n> 0: más cálido que lo normal",
       "fig_contraste_2_era5_periodos.png",
       f"Una sola fuente (ERA5, promedio mensual de las 12 UTC, celda {CELDA[0]:.2f}°, {CELDA[1]:.2f}°): la diferencia "
       f"muestra si los años con sondeos fueron típicos del periodo común. n = meses de ERA5\n(periodo común / meses de los "
       f"sondeos). Los meses 2024-10 a 2024-12 están en ambos. Scripts 34a y 34c.")
