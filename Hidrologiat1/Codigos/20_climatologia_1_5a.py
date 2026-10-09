"""
Punto 1.5.a: climatología de 12 meses (ciclo anual) de P_L, P_I, Q̄, R y T en la cuenca El Tambo,
sobre el periodo común de análisis (meses con P_L e IMERG; periodo_comun.py).

Para cada mes calendario: n.º de años válidos, media, mediana, desviación estándar (ddof = 1),
cuartiles, P10 y P90 (percentiles tipo 7, interpolación lineal). P_L y P_I usan los mismos meses
(coinciden en todo el periodo común). Q̄ y R usan sus meses válidos (criterio de completitud del 09).
T: PISCOt (CAMELS-PE) hasta 2020-12 y, desde 2021-01, la estimación ERA5-Land + sesgo mensual del
script 17 (columna aparte); se informa cuántos años de cada mes son estimados.

Figuras:
  fig_climatologia_1_5a.png      media y mediana con bandas IQR (Q1–Q3) y P10–P90. Las bandas describen
                                 la dispersión entre años; no son intervalos de confianza de la media.
  fig_anio_mes_1_5a.png          mapas año × mes que conservan los años individuales.
Tabla: climatologia_1_5a.csv (formato largo: variable × mes).
Requiere 09, 10 y 17.
"""
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from periodo_comun import periodo_comun

ROOT = Path(__file__).resolve().parents[2] / "Tarea 1" / "El_Tambo"
SM = ROOT / "resultados" / "series_mensuales"
OUT = ROOT / "resultados" / "climatologia"
OUT.mkdir(parents=True, exist_ok=True)

INI, FIN = periodo_comun()
cam = pd.read_csv(SM / "series_mensuales_camels_pe.csv", index_col="fecha", parse_dates=True).loc[INI:FIN]
imerg = pd.read_csv(SM / "imerg_cuenca_mensual.csv", index_col="fecha", parse_dates=True)["P_I_mm_mes"].loc[INI:FIN]
text = pd.read_csv(SM / "temperatura_piscot_extendida_era5land.csv", index_col="fecha", parse_dates=True).loc[INI:FIN]
T = text["T_piscot_C"].fillna(text["T_estimada_era5land_C"])          # PISCOt + estimación (no se mezclan filas)
T_est = text["T_piscot_C"].isna() & text["T_estimada_era5land_C"].notna()

df = pd.DataFrame({"P_L": cam["P_L_mm_mes"], "P_I": imerg, "Q": cam["Q_m3_s"], "R": cam["R_mm_mes"], "T": T})
VARS = {  # clave: (etiqueta, unidades, color)
    "P_L": ("P$_L$ (CAMELS-PE)", "mm/mes", "#1b6ca8"),
    "P_I": ("P$_I$ (IMERG V07)", "mm/mes", "#e4572e"),
    "Q": ("$\\bar{Q}$", "m³/s", "#0b7a75"),
    "R": ("R", "mm/mes", "#5b3f8c"),
    "T": ("T (PISCOt + estimación ERA5-Land desde 2021)", "°C", "#c8553d"),
}
MESES = ["Ene", "Feb", "Mar", "Abr", "May", "Jun", "Jul", "Ago", "Sep", "Oct", "Nov", "Dic"]
assert df[["P_L", "P_I"]].notna().all().all(), "P_L y P_I deben tener los mismos meses en el periodo común"

# ------------------------------------------------------------------ tabla de climatología
filas = []
for k in VARS:
    for mes in range(1, 13):
        s = df.loc[df.index.month == mes, k]
        v = s.dropna().values
        p = np.percentile(v, [10, 25, 50, 75, 90])
        fila = {"variable": k, "unidades": VARS[k][1], "mes": mes, "mes_nombre": MESES[mes - 1],
                "n_anios_validos": len(v), "n_anios_periodo": len(s), "media": v.mean(), "mediana": p[2],
                "desv_est": v.std(ddof=1), "Q1": p[1], "Q3": p[3], "P10": p[0], "P90": p[4],
                "minimo": v.min(), "maximo": v.max()}
        if k == "T":
            fila["n_anios_estimados"] = int(T_est[T_est.index.month == mes].sum())
        filas.append(fila)
clim = pd.DataFrame(filas)
clim.round(3).to_csv(OUT / "climatologia_1_5a.csv", index=False)
pd.set_option("display.width", 250, "display.max_columns", 30)
print(f"Periodo común {INI:%Y-%m} a {FIN:%Y-%m}")
for k in VARS:
    c = clim[clim.variable == k].set_index("mes_nombre")
    cols = ["n_anios_validos", "media", "mediana", "desv_est", "Q1", "Q3", "P10", "P90"] + \
           (["n_anios_estimados"] if k == "T" else [])
    print(f"\n=== {k} ({VARS[k][1]}) ===\n{c[cols].round(2).T.to_string()}")

# Resumen del ciclo (apoyo para 1.5.b)
for k in ("P_L", "P_I", "Q", "R"):
    c = clim[clim.variable == k].set_index("mes")["media"]
    tot = c.sum()
    djfma = c.loc[[12, 1, 2, 3, 4]].sum() / tot * 100
    print(f"{k}: máx medio {MESES[c.idxmax() - 1]} ({c.max():.1f}), mín {MESES[c.idxmin() - 1]} ({c.min():.1f}); "
          f"{djfma:.0f} % del total anual medio en dic–abr")
cT = clim[clim.variable == "T"].set_index("mes")["media"]
print(f"T: máx {MESES[cT.idxmax() - 1]} ({cT.max():.2f} °C), mín {MESES[cT.idxmin() - 1]} ({cT.min():.2f} °C), "
      f"amplitud {cT.max() - cT.min():.2f} °C")

# ------------------------------------------------------------------ figura 1: ciclo anual con bandas
x = np.arange(1, 13)


def ciclo(ax, k, desplaz=0.0):
    c = clim[clim.variable == k].set_index("mes")
    col = VARS[k][2]
    ax.fill_between(x + desplaz, c["P10"], c["P90"], color=col, alpha=0.13, lw=0)
    ax.fill_between(x + desplaz, c["Q1"], c["Q3"], color=col, alpha=0.30, lw=0)
    ax.plot(x + desplaz, c["media"], color=col, lw=2, marker="o", ms=4, label=f"{VARS[k][0]} – media")
    ax.plot(x + desplaz, c["mediana"], color=col, lw=1.5, ls="--", label=f"{VARS[k][0]} – mediana")


fig, axs = plt.subplots(2, 2, figsize=(15, 10))
fig.subplots_adjust(hspace=0.32, wspace=0.16, top=0.88, bottom=0.1, left=0.06, right=0.98)
paneles = [(axs[0, 0], ["P_L", "P_I"], "(a) Precipitación", "mm/mes"),
           (axs[0, 1], ["Q"], "(b) Caudal medio mensual", "m³/s"),
           (axs[1, 0], ["R"], "(c) Escorrentía en lámina", "mm/mes"),
           (axs[1, 1], ["T"], "(d) Temperatura media", "°C")]
for ax, ks, titulo, unid in paneles:
    for k in ks:
        ciclo(ax, k)
    n = clim[clim.variable == ks[0]].set_index("mes")["n_anios_validos"]
    sub = f"años válidos por mes: {n.min()}–{n.max()}"
    if ks == ["T"]:
        ne = clim[clim.variable == "T"].set_index("mes")["n_anios_estimados"].astype(int)
        sub = f"{n.min()}–{n.max()} años por mes, {ne.min()}–{ne.max()} estimados"
    ax.set_title(f"{titulo} · {sub}", loc="left", fontsize=10.5, fontweight="bold")
    ax.set_xticks(x, MESES)
    ax.set_ylabel(unid)
    ax.grid(alpha=0.3)
    ax.set_xlim(0.6, 12.4)
    if unid != "°C":
        ax.set_ylim(bottom=0)
    ax.legend(frameon=False, fontsize=8.5, loc="upper right" if ks != ["T"] else "lower left")
fig.suptitle(f"Cuenca del río Chicama hasta El Tambo – ciclo anual en el periodo común {INI:%Y-%m} a {FIN:%Y-%m}",
             fontsize=13, y=0.97)
fig.text(0.5, 0.925, "Línea continua = media; discontinua = mediana; banda oscura = rango intercuartílico (Q1–Q3); "
         "banda clara = P10–P90. Las bandas describen la dispersión entre años, no son intervalos de confianza.",
         ha="center", fontsize=9)
fig.text(0.01, 0.02, "P$_L$ y P$_I$ con los mismos meses. $\\bar{Q}$ y R con sus meses válidos (≤ 10 % de días "
         "faltantes). T: PISCOt v1.2 hasta 2020-12 y ERA5-Land + sesgo mensual vs PISCOt desde 2021-01 (script 17). "
         "Percentiles tipo 7.", fontsize=8)
fig.savefig(OUT / "fig_climatologia_1_5a.png", dpi=200)
plt.close(fig)

# ------------------------------------------------------------------ figura 2: mapas año × mes
anios = np.arange(INI.year, FIN.year + 1)


def matriz(s):
    g = s.to_frame("v").assign(a=s.index.year, m=s.index.month)
    return g.pivot(index="a", columns="m", values="v").reindex(index=anios, columns=range(1, 13))


fig, axs = plt.subplots(1, 4, figsize=(16, 9), sharey=True, gridspec_kw={"wspace": 0.1})
pmax = df[["P_L", "P_I"]].max().max()
specs = [("P_L", "Blues", dict(vmin=0, vmax=pmax), "mm/mes"),
         ("P_I", "Blues", dict(vmin=0, vmax=pmax), "mm/mes"),
         ("Q", "GnBu", dict(norm=LogNorm(vmin=max(df["Q"].min(), 0.1), vmax=df["Q"].max())), "m³/s (escala log)"),
         ("T", "OrRd", dict(vmin=df["T"].min(), vmax=df["T"].max()), "°C")]
for ax, (k, cmap, kw, unid) in zip(axs, specs):
    M = matriz(df[k])
    cm = plt.get_cmap(cmap).copy()
    cm.set_bad("#d9d9d9")
    im = ax.imshow(M.values, cmap=cm, aspect="auto", interpolation="nearest",
                   extent=(0.5, 12.5, anios[-1] + 0.5, anios[0] - 0.5), **kw)
    ax.set_xticks(range(1, 13), [m[0] for m in MESES])
    ax.set_title(VARS[k][0].replace(" (PISCOt + estimación ERA5-Land desde 2021)", " (PISCOt; ▨ estimada)"), fontsize=10)
    cb = fig.colorbar(im, ax=ax, orientation="horizontal", fraction=0.035, pad=0.05)
    cb.set_label(unid, fontsize=8.5)
    if k == "T":                                           # años estimados con trama
        for f in T_est.index[T_est]:
            ax.add_patch(plt.Rectangle((f.month - 0.5, f.year - 0.5), 1, 1, fill=False, hatch="////",
                                       edgecolor="#555555", lw=0))
axs[0].set_yticks(anios[::2])
axs[0].set_ylabel("Año")
fig.suptitle(f"Valores mensuales por año y mes calendario – periodo común {INI:%Y-%m} a {FIN:%Y-%m} "
             "(gris: sin dato o fuera del periodo)", fontsize=12, y=0.97)
fig.text(0.01, 0.01, "P$_L$ y P$_I$ con la misma escala de color. $\\bar{Q}$ en escala logarítmica por su fuerte "
         "asimetría. T con trama: estimación ERA5-Land + sesgo mensual (no es PISCOt).", fontsize=8.5)
fig.subplots_adjust(top=0.92, bottom=0.08, left=0.05, right=0.99)
fig.savefig(OUT / "fig_anio_mes_1_5a.png", dpi=200)
print(f"\nSalidas en {OUT}")
