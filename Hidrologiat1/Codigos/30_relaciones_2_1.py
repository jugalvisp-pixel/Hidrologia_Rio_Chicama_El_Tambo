"""
Punto 2.1: diagramas de dispersión y relaciones entre series mensuales (cuenca El Tambo).

Pares (meses coincidentes y válidos del periodo común, periodo_comun.py):
  (a) P_I (IMERG) frente a P_L (CAMELS-PE), con línea 1:1 y escalas iguales.
  (b) P_L frente a Q̄ (m³/s); apoyo en lámina: P_L frente a R (mm/mes) con la línea R = P.
  (c) P_I frente a Q̄.
Puntos coloreados por estación climática (DEF, MAM, JJA, SON). Q̄ y R usan la serie rellenada (script 09);
como sensibilidad se repiten las métricas solo con los meses observados completos (R_flag = completo).
Métricas: n, Pearson r (y p), Spearman ρ; para P_I − P_L: sesgo medio, MAE, RMSE, también por estación.
Rezagos: correlación de Q̄_t con P_(t−k), k = 0..6 meses (la lluvia antecede al caudal: nunca se usa lluvia
futura), y con la lluvia acumulada de los últimos k+1 meses (P_t + … + P_(t−k)).
Anomalías: a_t = X_t − μ_j(t), con μ_j la media del mes calendario j en el periodo común (misma climatología
del apartado 1.5; definición del punto 3.2). Se repiten las correlaciones y los rezagos con anomalías.

Salidas (Tarea 1/El_Tambo/resultados/relaciones/):
  metricas_2_1.csv, metricas_IMERG_por_estacion_2_1.csv, rezagos_2_1.csv,
  fig_dispersion_2_1.png, fig_dispersion_anomalias_2_1.png, fig_rezagos_2_1.png
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from periodo_comun import periodo_comun

ROOT = Path(__file__).resolve().parents[2] / "Tarea 1" / "El_Tambo"
SM = ROOT / "resultados" / "series_mensuales"
OUT = ROOT / "resultados" / "relaciones"
OUT.mkdir(parents=True, exist_ok=True)
INI, FIN = periodo_comun()
EST = {12: "DEF", 1: "DEF", 2: "DEF", 3: "MAM", 4: "MAM", 5: "MAM", 6: "JJA", 7: "JJA", 8: "JJA",
       9: "SON", 10: "SON", 11: "SON"}
COL_EST = {"DEF": "#1b6ca8", "MAM": "#0b7a75", "JJA": "#e4572e", "SON": "#c9a227"}
REZ = range(0, 7)

cam = pd.read_csv(SM / "series_mensuales_camels_pe.csv", index_col="fecha", parse_dates=True)
imerg = pd.read_csv(SM / "imerg_cuenca_mensual.csv", index_col="fecha", parse_dates=True)["P_I_mm_mes"]
# Rezagos: se toma la lluvia desde 12 meses antes del inicio para que P_(t−k) exista desde el primer mes
df_l = pd.DataFrame({"P_L": cam["P_L_mm_mes"], "P_I": imerg}).loc[INI - pd.DateOffset(months=12):FIN]
df = pd.DataFrame({"P_L": cam["P_L_mm_mes"], "P_I": imerg, "Q": cam["Q_m3_s"], "R": cam["R_mm_mes"],
                   "flag": cam["R_flag"]}).loc[INI:FIN]
df["est"] = df.index.month.map(EST)
assert df[["P_L", "P_I", "Q", "R"]].notna().all().all(), "faltan datos en el periodo común"
num = ["P_L", "P_I", "Q", "R"]
clim = df[num].groupby(df.index.month).transform("mean")
an = df[num] - clim
an["est"], an["flag"] = df["est"], df["flag"]
n_anios = df.groupby(df.index.month).size()


def correl(x, y):
    r, p = stats.pearsonr(x, y)
    rho, prho = stats.spearmanr(x, y)
    return r, p, rho, prho


# ------------------------------------------------------------------ métricas
filas = []
for nombre, (xv, yv) in {"(a) P_I vs P_L": ("P_L", "P_I"), "(b) P_L vs Q̄": ("P_L", "Q"),
                         "(b') P_L vs R": ("P_L", "R"), "(c) P_I vs Q̄": ("P_I", "Q")}.items():
    for serie, d in (("original", df), ("anomalías", an)):
        for muestra, dd in (("todos los meses", d), ("solo Q observado completo", d[d["flag"] == "completo"])):
            if muestra != "todos los meses" and "Q" not in (xv, yv) and "R" not in (xv, yv):
                continue
            r, p, rho, prho = correl(dd[xv], dd[yv])
            filas.append({"par": nombre, "serie": serie, "muestra": muestra, "n": len(dd), "pearson_r": r,
                          "pearson_p": p, "spearman_rho": rho, "spearman_p": prho})
met = pd.DataFrame(filas)
met.round(4).to_csv(OUT / "metricas_2_1.csv", index=False)


def errores(d):
    e = d["P_I"] - d["P_L"]
    return pd.Series({"n": len(d), "media_P_L": d["P_L"].mean(), "media_P_I": d["P_I"].mean(), "sesgo": e.mean(),
                      "sesgo_pct": 100 * e.sum() / d["P_L"].sum(), "MAE": e.abs().mean(),
                      "RMSE": np.sqrt((e ** 2).mean()), "pearson_r": d[["P_I", "P_L"]].corr().iloc[0, 1]})


err = pd.concat({"todo el año": errores(df), **{e: errores(df[df.est == e]) for e in COL_EST}}, axis=1).T
err.round(3).to_csv(OUT / "metricas_IMERG_por_estacion_2_1.csv")
print(f"Periodo {INI:%Y-%m}..{FIN:%Y-%m}, n = {len(df)} meses ({n_anios.min()}–{n_anios.max()} años por mes)\n")
print(met.round(3).to_string(index=False))
print("\nP_I − P_L por estación:\n", err.round(2).to_string())

# ------------------------------------------------------------------ rezagos
clim_l = df_l.groupby(df_l.index.month).transform("mean")         # climatología del mismo periodo extendido
an_l = df_l - df_l.loc[INI:FIN].groupby(df_l.loc[INI:FIN].index.month).mean().reindex(df_l.index.month).values
filas = []
for serie, base, lluvia in (("original", df, df_l), ("anomalías", an, an_l)):
    for p in ("P_L", "P_I"):
        for k in REZ:
            x_rez = lluvia[p].shift(k).reindex(base.index)
            x_acu = lluvia[p].rolling(k + 1).sum().reindex(base.index)
            for tipo, x in (("rezago P(t−k)", x_rez), ("acumulada P(t−k..t)", x_acu)):
                ok = x.notna()
                r, _, rho, _ = correl(x[ok], base.loc[ok, "Q"])
                filas.append({"serie": serie, "lluvia": p, "tipo": tipo, "k": k, "n": int(ok.sum()),
                              "pearson_r": r, "spearman_rho": rho})
rez = pd.DataFrame(filas)
rez.round(4).to_csv(OUT / "rezagos_2_1.csv", index=False)
print("\nMejor rezago/acumulación (Spearman) por serie y lluvia:")
for (s, p, t), g in rez.groupby(["serie", "lluvia", "tipo"]):
    b = g.loc[g.spearman_rho.idxmax()]
    print(f"  {s:9s} {p} {t:20s}: k = {int(b.k)} ρ = {b.spearman_rho:.3f} (k=0: {g[g.k == 0].spearman_rho.iloc[0]:.3f})")


# ------------------------------------------------------------------ figuras de dispersión
def dispersion(d, titulo_fig, archivo, anom):
    fig, axs = plt.subplots(2, 2, figsize=(13, 11.5))
    pares = [("P_L", "P_I", "(a) IMERG frente a lluvia local", "P$_L$ (mm/mes)", "P$_I$ (mm/mes)"),
             ("P_L", "Q", "(b) Lluvia local frente a caudal", "P$_L$ (mm/mes)", "$\\bar{Q}$ (m³/s)"),
             ("P_I", "Q", "(c) IMERG frente a caudal", "P$_I$ (mm/mes)", "$\\bar{Q}$ (m³/s)"),
             ("P_L", "R", "(b') Lluvia local frente a escorrentía en lámina", "P$_L$ (mm/mes)", "R (mm/mes)")]
    for ax, (xv, yv, tit, xl, yl) in zip(axs.flat, pares):
        if anom:
            xl, yl = "Anomalía de " + xl, "Anomalía de " + yl
        for e, c in COL_EST.items():
            s = d[d.est == e]
            ax.scatter(s[xv], s[yv], s=18, color=c, alpha=0.7, edgecolor="white", lw=0.4, label=e)
        if yv in ("P_I", "R"):                                    # línea 1:1 y escalas iguales
            lo = min(d[xv].min(), d[yv].min()); hi = max(d[xv].max(), d[yv].max())
            pad = 0.04 * (hi - lo)
            ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], color="#3d3d3d", lw=1,
                    label="1:1" if yv == "P_I" else "R = P$_L$")
            ax.set_xlim(lo - pad, hi + pad); ax.set_ylim(lo - pad, hi + pad); ax.set_aspect("equal")
        if anom:
            ax.axhline(0, color="#999", lw=0.6); ax.axvline(0, color="#999", lw=0.6)
        elif yv == "Q":                              # caudal: varía 3 órdenes de magnitud → eje logarítmico
            ax.set_yscale("log")
            yl += " – escala logarítmica"
        m = met[(met.par.str.startswith(tit[:3].replace("(b'", "(b'"))) & (met.serie == ("anomalías" if anom else "original"))
                & (met.muestra == "todos los meses")]
        if len(m):
            r = m.iloc[0]
            txt = f"n = {int(r.n)}\nPearson r = {r.pearson_r:.2f}\nSpearman ρ = {r.spearman_rho:.2f}"
            if yv == "P_I" and not anom:
                e = err.loc["todo el año"]
                txt += f"\nsesgo = {e.sesgo:+.1f} mm/mes ({e.sesgo_pct:+.0f} %)\nMAE = {e.MAE:.1f} · RMSE = {e.RMSE:.1f}"
            ax.text(0.03, 0.97, txt, transform=ax.transAxes, va="top", fontsize=9,
                    bbox=dict(fc="white", alpha=0.85, lw=0))
        ax.set_xlabel(xl); ax.set_ylabel(yl); ax.set_title(tit, loc="left", fontsize=11)
        ax.grid(alpha=0.3)
    axs[0, 0].legend(loc="lower right", fontsize=8.5, frameon=False, title="Estación", title_fontsize=8.5)
    fig.suptitle(titulo_fig, fontsize=12.5)
    fig.text(0.5, 0.01, f"Periodo común {INI:%Y-%m} a {FIN:%Y-%m}, {len(d)} meses coincidentes. DEF = dic–feb, MAM = mar–may, "
             "JJA = jun–ago, SON = sep–nov. Q̄ y R con faltantes rellenados (script 09); ver la sensibilidad sin "
             "relleno en metricas_2_1.csv." + (" Anomalías: valor menos la media de su mes calendario en el periodo "
                                                "común." if anom else ""), ha="center", fontsize=8, color="#444", wrap=True)
    fig.tight_layout(rect=(0, 0.03, 1, 0.97))
    fig.savefig(OUT / archivo, dpi=200)
    plt.close(fig)


dispersion(df, "Punto 2.1 – Relaciones mensuales entre lluvia y caudal (series originales)", "fig_dispersion_2_1.png", False)
dispersion(an, "Punto 2.1 – Relaciones mensuales en anomalías (sin el ciclo anual compartido)",
           "fig_dispersion_anomalias_2_1.png", True)

# ------------------------------------------------------------------ figura de rezagos
fig, axs = plt.subplots(1, 2, figsize=(13, 4.8), sharey=True)
for ax, serie in zip(axs, ("original", "anomalías")):
    for p, col in (("P_L", "#1b6ca8"), ("P_I", "#e4572e")):
        for tipo, ls, mk in (("rezago P(t−k)", "-", "o"), ("acumulada P(t−k..t)", "--", "s")):
            g = rez[(rez.serie == serie) & (rez.lluvia == p) & (rez.tipo == tipo)]
            ax.plot(g.k, g.spearman_rho, ls=ls, marker=mk, color=col, ms=5,
                    label=f"{'P$_L$' if p == 'P_L' else 'P$_I$'} – {'lluvia del mes t−k' if tipo.startswith('rez') else 'lluvia acumulada de t−k a t'}")
    ax.axhline(0, color="#999", lw=0.7)
    ax.set_xlabel("k (meses de anticipación de la lluvia respecto al caudal)")
    ax.set_title(f"({'ab'[serie == 'anomalías']}) {'Series originales' if serie == 'original' else 'Anomalías'}",
                 loc="left", fontsize=11)
    ax.grid(alpha=0.3); ax.set_xticks(list(REZ))
axs[0].set_ylabel("Spearman ρ con $\\bar{Q}_t$")
axs[0].legend(fontsize=8, frameon=False, loc="lower left")
fig.suptitle("Punto 2.1 – ¿La lluvia de meses anteriores ayuda a explicar el caudal?", fontsize=12.5)
fig.text(0.5, -0.02, "La lluvia siempre antecede al caudal (k ≥ 0). Acumulada: suma de los últimos k+1 meses, incluido el "
         f"mes t. Periodo común {INI:%Y-%m} a {FIN:%Y-%m}.", ha="center", fontsize=8, color="#444")
fig.tight_layout()
fig.savefig(OUT / "fig_rezagos_2_1.png", dpi=200, bbox_inches="tight")
plt.close(fig)
print(f"\nSalidas en {OUT}")
