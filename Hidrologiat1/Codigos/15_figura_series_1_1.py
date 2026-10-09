"""
Punto 1.1: series mensuales de precipitación, caudal y temperatura de la cuenca El Tambo en
paneles alineados (mismo eje de tiempo, unidades explícitas y vacíos visibles).

Series construidas a partir de los registros diarios de CAMELS-PE en 09_datos_camels_pe_el_tambo.py:
  P_L,m  = Σ P_d                       [mm/mes]  solo meses completos
  Q̄_m    = (1/n_m) Σ Q_d               [m³/s]    media de días válidos si faltan ≤ 10 %
  R_m    = Σ q_d (q ya en mm/d)        [mm/mes]  mes completo: suma; ≤ 10 % faltantes: escalada
  T_m    = media de T_d                [°C]      si faltan ≤ 10 %
Periodo: periodo común de análisis, meses con P_L e IMERG (periodo_comun.py). Meses sin dato o
excluidos quedan como vacíos (sin interpolar) y se sombrean. T 2021 en adelante: estimación del
script 17 (ERA5-Land + sesgo mensual), dibujada aparte.

Salidas (Tarea 1/El_Tambo/resultados/series_mensuales/):
  fig_series_mensuales_1_1.png
  extremos_series_mensuales_1_1.csv
"""
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

from periodo_comun import periodo_comun

ROOT = Path(__file__).resolve().parents[2] / "Tarea 1" / "El_Tambo"
OUT = ROOT / "resultados" / "series_mensuales"
ini, fin = periodo_comun()
m = pd.read_csv(OUT / "series_mensuales_camels_pe.csv", parse_dates=["fecha"], index_col="fecha").loc[ini:fin]
EXT_T = OUT / "temperatura_piscot_extendida_era5land.csv"      # script 17
VAL_T = pd.read_csv(OUT / "validacion_extension_T.csv", index_col=0)["valor"]
COLOR_EST = "#6a4c93"                                            # tramo estimado de T

PANELES = [  # columna, etiqueta, unidades, color
    ("P_L_mm_mes", "P$_L$ – precipitación acumulada mensual (CAMELS-PE, PISCOp v2.1)", "mm/mes", "#1b6ca8"),
    ("Q_m3_s", "$\\bar{Q}$ – caudal medio mensual (CAMELS-PE, SENAMHI)", "m³/s", "#0b7a75"),
    ("R_mm_mes", "R – escorrentía mensual en lámina", "mm/mes", "#5b3f8c"),
    ("T_C", "T – temperatura media mensual (CAMELS-PE, PISCOt v1.2)", "°C", "#c8553d"),
]


def sombrear(ax, faltan, color, alpha):
    """Sombrea cada mes sin dato como un bloque del primer al último día del mes."""
    for t in faltan:
        ax.axvspan(t, t + pd.offsets.MonthBegin(1), color=color, alpha=alpha, lw=0, zorder=0)


fig, axs = plt.subplots(len(PANELES), 1, figsize=(15, 11.5), sharex=True,
                        gridspec_kw={"hspace": 0.12, "height_ratios": [1.15, 1, 1, 0.8]})
for ax, (col, titulo, unid, color) in zip(axs, PANELES):
    s = m[col]
    ax.plot(s.index, s.values, color=color, lw=1.2, marker="o", ms=1.8)   # NaN corta la línea
    ax.set_ylabel(unid)
    ax.set_title(titulo, loc="left", fontsize=10.5)
    ax.grid(alpha=0.3)
    faltan = s.index[s.isna()]
    if col != "T_C":
        sombrear(ax, faltan, "#9e9e9e", 0.35)
    v = s.dropna()
    # Resumen encima del panel (a la altura del título) para no tapar datos
    ax.text(1.0, 1.01, f"{len(v)}/{len(s)} meses válidos · media {v.mean():.1f} · máx {v.max():.1f} "
            f"({v.idxmax():%Y-%m})", transform=ax.transAxes, ha="right", va="bottom", fontsize=8.5)
    if col in ("Q_m3_s", "R_mm_mes"):
        esc = m.index[m["R_flag"] == "escalado"]
        ax.plot(esc, s.loc[esc], ls="none", marker="o", ms=4, mfc="none", mec="#e08e0b", mew=0.9)
    if col == "T_C" and len(faltan):
        # Tramo estimado (script 17): ERA5-Land + sesgo mensual medio frente a PISCOt en el periodo
        # común. Se dibuja aparte (discontinuo, fondo propio) y no reemplaza a T_C en los cálculos.
        te = pd.read_csv(EXT_T, parse_dates=["fecha"], index_col="fecha")["T_estimada_era5land_C"].loc[ini:fin].dropna()
        ax.axvspan(te.index.min(), te.index.max() + pd.offsets.MonthBegin(1), color=COLOR_EST, alpha=0.12,
                   lw=0, zorder=0)
        enlace = pd.concat([s.loc[[v.index.max()]], te])                 # une visualmente 2020-12 → 2021-01
        ax.plot(enlace.index, enlace.values, color=COLOR_EST, lw=1.2, ls="--", marker="o", ms=1.8)
        ax.set_title(titulo + f" · {te.index.min():%Y-%m} a {te.index.max():%Y-%m} estimada (fondo morado, "
                     "ver leyenda)", loc="left", fontsize=10.5)
        ax.set_ylim(top=max(s.max(), te.max()) + 0.25)
axs[0].set_ylim(bottom=0)
axs[-1].xaxis.set_major_locator(mdates.YearLocator(2))
axs[-1].xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
axs[-1].xaxis.set_minor_locator(mdates.YearLocator(1))
axs[-1].set_xlim(ini - pd.Timedelta(days=15), fin + pd.offsets.MonthBegin(1) + pd.Timedelta(days=15))

fig.legend(handles=[Patch(fc="#9e9e9e", alpha=0.35, label="Mes sin dato o excluido (> 10 % de días faltantes)"),
                    Line2D([], [], ls="none", marker="o", mfc="none", mec="#e08e0b",
                           label="Mes con ≤ 10 % de días faltantes (Q̄: media de días válidos; R escalada a n$_m$)"),
                    Line2D([], [], color=COLOR_EST, ls="--", lw=1.2,
                           label=f"T desde 2021 estimada: ERA5-Land + sesgo mensual vs PISCOt (no es PISCOt; "
                                 f"validación {VAL_T['evaluacion']}: MAE {float(VAL_T['MAE_C']):.2f} °C)")],
           loc="lower center", ncol=2, fontsize=8.5, frameon=False, bbox_to_anchor=(0.5, 0.0))
fig.suptitle("Cuenca del río Chicama hasta El Tambo (CAMELS-PE PE_47E0A2A8, 2171.6 km²) – series "
             f"mensuales {ini:%Y-%m} a {fin:%Y-%m}", fontsize=13, y=0.985)
fig.text(0.5, 0.955, "Agregadas desde registros diarios: P$_L$ = ΣP$_d$; $\\bar{Q}$ = media de Q$_d$; "
         "R = Σq$_d$ (q en mm/d, sin repetir 86.4/A); T = media de T$_d$. Vacíos sin interpolar.",
         ha="center", fontsize=9)
fig.subplots_adjust(left=0.06, right=0.985, top=0.93, bottom=0.085)
fig.savefig(OUT / "fig_series_mensuales_1_1.png", dpi=200)

# ------------------------------------------------------------------ extremos (apoyo para describir 1.1)
filas = []
for col, _, unid, _ in PANELES:
    v = m[col].dropna()
    for tipo, f in (("máximo", v.idxmax()), ("mínimo", v.idxmin())):
        filas.append({"variable": col, "unidades": unid, "tipo": tipo, "fecha": f"{f:%Y-%m}", "valor": v[f]})
ext = pd.DataFrame(filas)
ext.round(2).to_csv(OUT / "extremos_series_mensuales_1_1.csv", index=False)
print(ext.round(2).to_string(index=False))
rp = m[["P_L_mm_mes", "R_mm_mes"]].dropna()
print(f"\nMeses con R > P_L: {(rp['R_mm_mes'] > rp['P_L_mm_mes']).sum()} de {len(rp)} → "
      + ", ".join(f"{t:%Y-%m}" for t in rp.index[rp['R_mm_mes'] > rp['P_L_mm_mes']]))
print(f"Figura: {OUT / 'fig_series_mensuales_1_1.png'}")
