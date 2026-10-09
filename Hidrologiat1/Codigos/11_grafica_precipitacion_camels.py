"""
Serie temporal de precipitación de CAMELS-PE (P_L, PISCOp v2.1 promediado sobre la cuenca)
en el periodo común de análisis (meses con P_L e IMERG; ver periodo_comun.py).

(a) Acumulado mensual (Σ diario) y media móvil centrada de 12 meses.
(b) Totales anuales (solo años con los 12 meses dentro del periodo común) y media.
Mismo formato que 16_grafica_precipitacion_imerg.py.

Requiere 09 y 10.
"""
from pathlib import Path

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

from periodo_comun import periodo_comun

ROOT = Path(__file__).resolve().parents[2] / "Tarea 1" / "El_Tambo"
SERIES = ROOT / "resultados" / "series_mensuales" / "series_mensuales_camels_pe.csv"
OUT = ROOT / "resultados" / "series_mensuales"

AZUL = "#2a78d6"           # serie 1 (paleta de referencia)
TINTA = "#0b0b0b"
TINTA_2 = "#52514e"
GRILLA = "#d9d8d4"

INI, FIN = periodo_comun()
m = pd.read_csv(SERIES, index_col="fecha", parse_dates=True)
P = m["P_L_mm_mes"].loc[INI:FIN]
mm12 = P.rolling(12, center=True, min_periods=12).mean()
anual = P.resample("YS").sum(min_count=12)
anual = anual[P.resample("YS").count() == 12]                             # años incompletos fuera
anual.index = anual.index.year
incompletos = [y for y in P.index.year.unique() if y not in anual.index]
media_anual = anual.mean()
top = P.nlargest(4)

print(f"Periodo común: {INI:%Y-%m} a {FIN:%Y-%m} | meses válidos {P.notna().sum()}/{len(P)}")
print(f"Media mensual {P.mean():.1f} mm | mediana {P.median():.1f} | máx {P.max():.1f} ({P.idxmax():%Y-%m})"
      f" | meses con P = 0: {(P == 0).sum()}")
print(f"Total anual medio {media_anual:.0f} mm ({anual.index.min()}–{anual.index.max()}) | mín {anual.min():.0f} "
      f"({anual.idxmin()}) | máx {anual.max():.0f} ({anual.idxmax()}) | años incompletos excluidos: {incompletos}")
print("Meses más lluviosos:", ", ".join(f"{t:%Y-%m} {v:.0f} mm" for t, v in top.items()))

plt.rcParams.update({"font.size": 10, "axes.edgecolor": TINTA_2, "axes.labelcolor": TINTA,
                     "xtick.color": TINTA_2, "ytick.color": TINTA_2})
fig, (a1, a2) = plt.subplots(2, 1, figsize=(15, 8.5), sharex=True,
                             gridspec_kw={"height_ratios": [3, 1.6], "hspace": 0.1})
for ax in (a1, a2):
    ax.grid(axis="y", color=GRILLA, lw=0.6)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

# (a) mensual
a1.plot(P.index, P.values, color=AZUL, lw=1.1, label="Acumulado mensual")
a1.plot(mm12.index, mm12.values, color=TINTA, lw=2, label="Media móvil centrada de 12 meses")
# Rótulos selectivos: el máximo de cada uno de los dos años más extremos (sin solapamientos)
picos = P.groupby(P.index.year).max().nlargest(2)
for anio in picos.index:
    t = P[str(anio)].idxmax()
    a1.annotate(f"{t:%m/%Y}: {P[t]:.0f} mm", (t, P[t]), xytext=(8, 0),
                textcoords="offset points", fontsize=8.5, color=TINTA_2, va="center")
a1.set_ylabel("Precipitación (mm/mes)")
a1.set_ylim(0, P.max() * 1.12)
a1.legend(loc="lower left", bbox_to_anchor=(0, 1.0), frameon=False, ncol=2)
a1.set_title("Precipitación de la cuenca del río Chicama hasta El Tambo – CAMELS-PE "
             f"(PISCOp v2.1), periodo común {INI:%Y-%m} a {FIN:%Y-%m}", loc="left", color=TINTA, pad=26)
a1.text(0.005, 0.98, "(a)", transform=a1.transAxes, fontweight="bold", va="top")

# (b) anual
fechas_anual = pd.to_datetime([f"{y}-07-01" for y in anual.index])
a2.bar(fechas_anual, anual.values, width=300, color=AZUL, edgecolor="white", lw=0)
a2.axhline(media_anual, color=TINTA_2, lw=1.2, ls="--")
x_fin = pd.Timestamp(f"{FIN.year + 1}-01-01")
a2.text(x_fin, media_anual, f" media\n {media_anual:.0f} mm/año", fontsize=8.5, color=TINTA_2,
        va="center", ha="left")
for y in incompletos:
    a2.text(pd.Timestamp(f"{y}-07-01"), anual.max() * 0.04, f"{y}\nincompleto", ha="center",
            va="bottom", fontsize=7.5, color=TINTA_2)
a2.set_ylim(0, anual.max() * 1.05)
a2.set_ylabel("Total anual (mm/año)")
a2.text(0.005, 0.97, "(b)", transform=a2.transAxes, fontweight="bold", va="top")
a2.xaxis.set_major_locator(mdates.YearLocator(2))
a2.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
a2.set_xlim(pd.Timestamp(f"{INI.year}-01-01") - pd.Timedelta(days=120), x_fin + pd.Timedelta(days=90))

fig.text(0.01, -0.005, "Fuente: CAMELS-PE v1.0.1 (doi:10.5281/zenodo.21195425), variable `prec` "
         "(PISCOp v2.1 promediado sobre la cuenca). Acumulado mensual = Σ precipitación diaria; todos los meses "
         "están completos.\nPeriodo común = meses con P$_L$ e IMERG; "
         f"{', '.join(map(str, incompletos))} incompleto (hasta {FIN:%m/%Y}) y excluido de (b).", fontsize=8, color=TINTA_2)
out = OUT / f"fig_precipitacion_camels_{INI:%Y}_{FIN:%Y}.png"
fig.savefig(out, dpi=200, bbox_inches="tight", facecolor="white")
print(f"Figura: {out}")
