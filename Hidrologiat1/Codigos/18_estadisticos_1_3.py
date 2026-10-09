"""
Punto 1.3: distribución de los datos mensuales de P_L, P_I, Q̄, R y T en la cuenca El Tambo.

Tablas (estadísticos descriptivos):
  (A) registro completo de cada variable: P_L 1981–2025, T 1981–2020 (solo PISCOt; la estimación
      2021–2025 con ERA5-Land NO entra), Q̄ y R 1998–2025, P_I 1998-01..2025-08. Solo tabla: la guía
      pide distinguirla de las estadísticas del periodo común.
  (B) periodo común de análisis (meses con P_L e IMERG; periodo_comun.py): cada variable con sus meses
      válidos dentro de ese periodo. P_L y P_I tienen exactamente los mismos meses.
Figuras y meses extremos: siempre sobre el periodo común (B).
Convenciones:
  - Faltantes: meses sin dato o excluidos por completitud (script 09) se descartan; no se rellenan.
    Un cero es un valor medido (P = 0 mm), distinto de un faltante (NaN).
  - Percentiles y cuartiles: interpolación lineal entre estadísticos de orden (Hyndman & Fan
    tipo 7, método por defecto de numpy). IQR = Q3 − Q1.
  - Desviación estándar muestral (ddof = 1). Asimetría: coeficiente de Fisher–Pearson ajustado (G1).
Histogramas:
  - Ancho de clase por la regla de Freedman–Diaconis, h = 2·IQR·n^(−1/3), redondeado a un valor
    cómodo (1, 2, 2.5, 5 × 10^k); límites desde 0 (o desde un múltiplo de h para T).
  - P_L y P_I: mismos meses (periodo común) y mismos límites de clase, calculados con la muestra
    conjunta.
  - Eje vertical en frecuencia relativa (% de meses): los anchos de clase son iguales, así que la
    frecuencia relativa es proporcional a la densidad y se lee directamente como fracción de meses.
Diagramas de caja: caja = Q1–Q3, línea = mediana, rombo = media, bigotes hasta el dato más extremo
dentro de 1.5·IQR (Tukey); puntos fuera de los bigotes se muestran, no se eliminan.

Requiere 09 (series CAMELS-PE) y 10 (serie IMERG).
Salidas (Tarea 1/El_Tambo/resultados/series_mensuales/):
  estadisticos_1_3_registro_completo.csv, estadisticos_1_3_periodo_comun.csv, extremos_1_3.csv,
  fig_histogramas_1_3.png, fig_cajas_1_3.png
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import skew
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from periodo_comun import periodo_comun

ROOT = Path(__file__).resolve().parents[2] / "Tarea 1" / "El_Tambo"
OUT = ROOT / "resultados" / "series_mensuales"

cam = pd.read_csv(OUT / "series_mensuales_camels_pe.csv", index_col="fecha", parse_dates=True)
imerg = pd.read_csv(OUT / "imerg_cuenca_mensual.csv", index_col="fecha", parse_dates=True)["P_I_mm_mes"]
df = pd.DataFrame({"P_L": cam["P_L_mm_mes"], "P_I": imerg, "Q": cam["Q_m3_s"], "R": cam["R_mm_mes"],
                   "T": cam["T_C"]})
VARS = {  # clave: (nombre, unidades, color)
    "P_L": ("P$_L$ (CAMELS-PE)", "mm/mes", "#1b6ca8"),
    "P_I": ("P$_I$ (IMERG V07)", "mm/mes", "#e4572e"),
    "Q": ("$\\bar{Q}$", "m³/s", "#0b7a75"),
    "R": ("R", "mm/mes", "#5b3f8c"),
    "T": ("T (PISCOt)", "°C", "#c8553d"),
}
PCT = [5, 10, 25, 50, 75, 90, 95]


def estadisticos(s):
    v = s.dropna().values
    p = np.percentile(v, PCT)                                 # tipo 7 (interpolación lineal)
    return {"n_meses_validos": len(v), "media": v.mean(), "mediana": p[3], "desv_est": v.std(ddof=1),
            "minimo": v.min(), "maximo": v.max(), "rango": v.max() - v.min(), "Q1": p[2], "Q3": p[4],
            "IQR": p[4] - p[2], "P5": p[0], "P10": p[1], "P90": p[5], "P95": p[6],
            "asimetria_G1": skew(v, bias=False), "media_menos_mediana": v.mean() - p[3],
            "meses_igual_0": int((v == 0).sum())}


def tabla(d, periodos):
    filas = []
    for k, s in d.items():
        st = estadisticos(s)
        filas.append({"variable": k, "unidades": VARS[k][1], "periodo": periodos[k],
                      "meses_faltantes_en_periodo": int(s.isna().sum()), **st})
    return pd.DataFrame(filas).set_index("variable")


def rango_valido(s):
    return s.loc[s.first_valid_index():s.last_valid_index()]


# (A) registro completo de cada variable (desde su primer hasta su último mes válido)
reg = {k: rango_valido(df[k]) for k in VARS}
per_A = {k: f"{s.index.min():%Y-%m} a {s.index.max():%Y-%m}" for k, s in reg.items()}
tA = tabla(reg, per_A)
# (B) periodo común P_L–P_I: cada variable con sus meses válidos dentro de él
INI, FIN = periodo_comun()
pc = {k: df[k].loc[INI:FIN] for k in VARS}
per_B = {k: f"{INI:%Y-%m} a {FIN:%Y-%m} (periodo común)" for k in VARS}
tB = tabla(pc, per_B)
pp = df.loc[INI:FIN, ["P_L", "P_I"]].dropna()              # mismos meses para las dos precipitaciones
assert len(pp) == len(pd.date_range(INI, FIN, freq="MS")), "P_L y P_I deben coincidir en todo el periodo común"

for t, nombre in ((tA, "registro_completo"), (tB, "periodo_comun")):
    t.round(3).to_csv(OUT / f"estadisticos_1_3_{nombre}.csv")
pd.set_option("display.width", 250, "display.max_columns", 30)
cols = ["periodo", "n_meses_validos", "media", "mediana", "desv_est", "minimo", "maximo", "Q1", "Q3", "IQR",
        "P5", "P10", "P90", "P95", "asimetria_G1", "meses_igual_0"]
print("=== (A) Registro completo de cada variable ===\n", tA[cols].round(2).to_string())
print(f"\n=== (B) Periodo común {INI:%Y-%m} a {FIN:%Y-%m} ===\n", tB[cols].round(2).to_string())


# ------------------------------------------------------------------ meses extremos y coherencia (periodo común)
filas = []
for k in VARS:
    s = pc[k].dropna()
    for tipo, sel in (("máximo", s.nlargest(3)), ("mínimo", s.nsmallest(3))):
        for rango_, (f, val) in enumerate(sel.items(), 1):
            fila = {"variable": k, "tipo": tipo, "orden": rango_, "fecha": f"{f:%Y-%m}", "valor": val}
            fila.update({f"{o}_mismo_mes": df.loc[f, o] for o in VARS if o != k})
            filas.append(fila)
ext = pd.DataFrame(filas)
ext.round(2).to_csv(OUT / "extremos_1_3.csv", index=False)
print("\n=== Meses extremos (valor de las otras variables en el mismo mes) ===\n", ext.round(1).to_string(index=False))
print(f"\nPeriodo común – meses con P_L < 1 mm: {(pc['P_L'] < 1).sum()} | P_I < 1 mm: {(pc['P_I'] < 1).sum()} | "
      f"Q̄ < 1 m³/s: {(pc['Q'] < 1).sum()} (ceros exactos: P_L {tB.loc['P_L', 'meses_igual_0']}, "
      f"P_I {tB.loc['P_I', 'meses_igual_0']}, Q̄ {tB.loc['Q', 'meses_igual_0']})")


# ------------------------------------------------------------------ histogramas
def ancho_bonito(x):
    v = np.asarray(x)
    h = 2 * (np.percentile(v, 75) - np.percentile(v, 25)) * len(v) ** (-1 / 3)
    e = 10 ** np.floor(np.log10(h))
    return min((m * e for m in (1, 2, 2.5, 5, 10)), key=lambda c: abs(c - h))


def bordes(x, h, desde_cero=True):
    lo = 0 if desde_cero else np.floor(np.min(x) / h) * h
    return np.arange(lo, np.max(x) + h, h)


fig, axs = plt.subplots(2, 2, figsize=(15, 9.5))
fig.subplots_adjust(hspace=0.38, wspace=0.18, top=0.9, bottom=0.08, left=0.06, right=0.98)
# (a) P_L vs P_I con los mismos meses y límites
hP = ancho_bonito(np.r_[pp["P_L"], pp["P_I"]])
bP = bordes(np.r_[pp["P_L"], pp["P_I"]], hP)
ax = axs[0, 0]
for k in ("P_L", "P_I"):
    w = np.full(len(pp), 100 / len(pp))
    ax.hist(pp[k], bins=bP, weights=w, histtype="stepfilled", alpha=0.3, color=VARS[k][2])
    ax.hist(pp[k], bins=bP, weights=w, histtype="step", lw=1.6, color=VARS[k][2],
            label=f"{VARS[k][0]}: media {pp[k].mean():.1f}, mediana {pp[k].median():.1f}")
ax.set_title(f"(a) Precipitación mensual (n = {len(pp)} meses, mismos meses)", loc="left", fontweight="bold",
             fontsize=10.5)
ax.set_xlabel(f"mm/mes (clases de {hP:g} mm)")
ax.legend(frameon=False, fontsize=9)
# (b)–(d) Q̄, R y T con sus meses válidos dentro del periodo común
for ax, k, letra in ((axs[0, 1], "Q", "b"), (axs[1, 0], "R", "c"), (axs[1, 1], "T", "d")):
    s = pc[k].dropna()
    h = ancho_bonito(s)
    b = bordes(s, h, desde_cero=(k != "T"))
    ax.hist(s, bins=b, weights=np.full(len(s), 100 / len(s)), color=VARS[k][2], alpha=0.85,
            edgecolor="white", lw=0.5)
    ax.axvline(s.mean(), color="k", lw=1.2, ls="--", label=f"media {s.mean():.1f}")
    ax.axvline(s.median(), color="k", lw=1.2, ls=":", label=f"mediana {s.median():.1f}")
    ax.set_title(f"({letra}) {VARS[k][0]} (n = {len(s)} meses válidos"
                 + (f", {s.index.min():%Y-%m} a {s.index.max():%Y-%m}" if k == "T" else "") + ")",
                 loc="left", fontweight="bold", fontsize=10.5)
    ax.set_xlabel(f"{VARS[k][1]} (clases de {h:g} {VARS[k][1]})")
    ax.legend(frameon=False, fontsize=9)
for ax in axs.flat:
    ax.set_ylabel("Frecuencia relativa (% de meses)")
    ax.grid(alpha=0.3)
fig.suptitle("Cuenca del río Chicama hasta El Tambo – distribución de los valores mensuales en el periodo común "
             f"{INI:%Y-%m} a {FIN:%Y-%m}", fontsize=13, y=0.975)
fig.text(0.5, 0.935, "Ancho de clase: regla de Freedman–Diaconis redondeada; P$_L$ y P$_I$ con los mismos meses y "
         "límites de clase. Frecuencia relativa con clases de igual ancho. Faltantes excluidos, sin relleno; "
         "T solo PISCOt (sin la estimación desde 2021).", ha="center", fontsize=9)
fig.savefig(OUT / "fig_histogramas_1_3.png", dpi=200)
plt.close(fig)

# ------------------------------------------------------------------ diagramas de caja
fig, axs = plt.subplots(1, 4, figsize=(15, 5.8), gridspec_kw={"width_ratios": [1.6, 1, 1, 1], "wspace": 0.35})
grupos = [(axs[0], [("P_L", pp["P_L"]), ("P_I", pp["P_I"])], f"(a) P (n = {len(pp)})", "mm/mes"),
          (axs[1], [("Q", pc["Q"].dropna())], f"(b) $\\bar{{Q}}$ (n = {pc['Q'].notna().sum()})", "m³/s"),
          (axs[2], [("R", pc["R"].dropna())], f"(c) R (n = {pc['R'].notna().sum()})", "mm/mes"),
          (axs[3], [("T", pc["T"].dropna())], f"(d) T PISCOt (n = {pc['T'].notna().sum()})", "°C")]
for ax, series, titulo, unid in grupos:
    bp = ax.boxplot([s.values for _, s in series], whis=1.5, widths=0.55, patch_artist=True, showmeans=True,
                    meanprops=dict(marker="D", markerfacecolor="white", markeredgecolor="k", ms=5),
                    medianprops=dict(color="k", lw=1.5), flierprops=dict(marker="o", ms=3, alpha=0.6))
    for patch, (k, _) in zip(bp["boxes"], series):
        patch.set_facecolor(VARS[k][2]); patch.set_alpha(0.55)
    ax.set_xticks(range(1, len(series) + 1), [VARS[k][0] for k, _ in series])
    ax.set_ylabel(unid)
    ax.set_title(titulo, loc="left", fontsize=10, fontweight="bold")
    ax.grid(axis="y", alpha=0.3)
fig.suptitle(f"Diagramas de caja de los valores mensuales, periodo común {INI:%Y-%m} a {FIN:%Y-%m} (caja Q1–Q3, "
             "línea = mediana, rombo = media, bigotes = 1.5·IQR de Tukey)", fontsize=11.5, y=1.0)
fig.savefig(OUT / "fig_cajas_1_3.png", dpi=200, bbox_inches="tight")
print(f"\nFiguras: {OUT / 'fig_histogramas_1_3.png'}, {OUT / 'fig_cajas_1_3.png'}")
