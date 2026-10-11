"""
Puntos 2.2 y 2.3: modelos estadísticos mensuales y su evaluación fuera del periodo de ajuste.

Objetivos: (i) estimar P_L (CAMELS-PE) a partir de P_I (IMERG); (ii) estimar Q̄ (m³/s) a partir de P_L o de P_I.

Partición por bloques de años completos (nunca meses al azar, que son dependientes):
  A (principal): ajuste 1998–2015 (18 años), evaluación 2016–2024 (9 años; incluye El Niño costero 2017 y los
    años con más relleno de caudal, 2019–2021).
  B (estabilidad): ajuste 2007–2024, evaluación 1998–2006 (incluye El Niño 1997–98).
Selección dentro del ajuste: validación cruzada dejando un año fuera (bloques anuales) entre los años de ajuste;
se elige el candidato y el rezago K de menor MAE de validación cruzada (el MAE representa el error típico y pesa
menos que el RMSE los pocos meses de crecida extrema). Las climatologías (referencia de caudal y anomalías)
se calculan SOLO con los años de ajuste.

Candidatos (i), P_L a partir de P_I [mm/mes]:
  ref   P̂_L = P_I                               IMERG sin corrección (referencia)
  lin   P̂_L = a + b·P_I                          regresión lineal (referencia de modelos)
  prop  P̂_L = b·P_I                              proporcional (sin intercepto)
  mes   P̂_L = b_j·P_I                            factor por mes calendario j (Σ P_L / Σ P_I del mes en el ajuste)
  pot   P̂_L = exp(a)·(P_I+1)^b·D − 1             potencial ajustada en ln(·+1), D = corrección de Duan
Candidatos (ii), Q̄ a partir de P (P_L o P_I) [m³/s; P en mm/mes; rezagos k = 0..K, la lluvia antecede]:
  ref   Q̂ = Q̄_clim,j                            climatología mensual del caudal (referencia)
  lin0  Q̂ = a + b·P_t                            regresión lineal simple
  linK  Q̂ = a + Σ b_k·P_(t−k)                    lineal con rezagos
  logK  ln Q̂ = a + Σ b_k·ln(P_(t−k)+1)           log-lineal con rezagos (Duan)
  logmK ln Q̂ = a_j + Σ b_k·ln(P_(t−k)+1)         log-lineal con intercepto por mes calendario (Duan)
  anoK  Q̂ = Q̄_clim,j + Σ b_k·(P_(t−k) − P̄_clim,j(t−k))   climatología + regresión de anomalías
Ceros: P_L e IMERG tienen meses de casi 0 mm → se usa ln(P+1); Q̄ > 0 siempre (mín ≈ 0.4 m³/s). Predicciones
negativas de los modelos lineales se cuentan como físicamente implausibles (no se recortan en las métricas).
Métricas fuera del ajuste: sesgo, MAE, RMSE, NSE y habilidad = 1 − MAE/MAE_ref. Sensibilidad: evaluación solo
con meses de caudal observado completo (R_flag = completo).

Salidas (Tarea 1/El_Tambo/resultados/modelos/):
  seleccion_cv_2_2.csv, evaluacion_2_3.csv, coeficientes_2_2.csv, predicciones_2_3.csv,
  fig_modelo_PL_2_3.png, fig_modelo_Q_2_3.png, fig_evaluacion_resumen_2_3.png
"""
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from periodo_comun import periodo_comun

ROOT = Path(__file__).resolve().parents[2] / "Tarea 1" / "El_Tambo"
SM = ROOT / "resultados" / "series_mensuales"
OUT = ROOT / "resultados" / "modelos"
OUT.mkdir(parents=True, exist_ok=True)
INI, FIN = periodo_comun()
PARTICIONES = {"A": (range(1998, 2016), range(2016, 2025)), "B": (range(2007, 2025), range(1998, 2007))}
KMAX = 4
MESES = ["Ene", "Feb", "Mar", "Abr", "May", "Jun", "Jul", "Ago", "Sep", "Oct", "Nov", "Dic"]

cam = pd.read_csv(SM / "series_mensuales_camels_pe.csv", index_col="fecha", parse_dates=True)
imerg = pd.read_csv(SM / "imerg_cuenca_mensual.csv", index_col="fecha", parse_dates=True)["P_I_mm_mes"]
base = pd.DataFrame({"P_L": cam["P_L_mm_mes"], "P_I": imerg, "Q": cam["Q_m3_s"], "flag": cam["R_flag"]})
base = base.loc[INI - pd.DateOffset(months=KMAX):FIN]            # meses previos para los rezagos
for p in ("P_L", "P_I"):
    for k in range(KMAX + 1):
        base[f"{p}_{k}"] = base[p].shift(k)
df = base.loc[INI:FIN].copy()
df["mes"], df["anio"] = df.index.month, df.index.year
assert df[["P_L", "P_I", "Q"]].notna().all().all()
# IMERG empieza en 1998-01: sus rezagos no existen en los primeros KMAX meses. Para Q̄ ← P_I se excluyen esos
# meses (se reporta n); P_L (CAMELS-PE) sí tiene datos previos, así que Q̄ ← P_L usa todos los meses.
DATOS = {"P_L": df, "P_I": df.dropna(subset=[f"P_I_{k}" for k in range(KMAX + 1)])}
print(f"Meses: {len(df)} (Q̄ ← P_I: {len(DATOS['P_I'])}, desde {DATOS['P_I'].index.min():%Y-%m})")


# ------------------------------------------------------------------ utilidades
def ols(X, y):
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    return b


def con1(*cols):
    return np.column_stack([np.ones(len(cols[0])), *cols])


def metricas(o, p):
    e = p - o
    return {"n": len(o), "sesgo": e.mean(), "MAE": e.abs().mean(), "RMSE": np.sqrt((e ** 2).mean()),
            "NSE": 1 - (e ** 2).sum() / ((o - o.mean()) ** 2).sum(), "r": np.corrcoef(o, p)[0, 1],
            "n_negativos": int((p < 0).sum())}


# ---- (i) P_L a partir de P_I
def fit_P(tr, nombre):
    x, y = tr["P_I"].values, tr["P_L"].values
    if nombre == "ref":
        return {}
    if nombre == "lin":
        return {"b": ols(con1(x), y)}
    if nombre == "prop":
        return {"b": np.array([x @ y / (x @ x)])}
    if nombre == "mes":
        return {"f": tr.groupby("mes")["P_L"].sum() / tr.groupby("mes")["P_I"].sum()}
    if nombre == "pot":
        ly = np.log(y + 1)
        b = ols(con1(np.log(x + 1)), ly)
        res = ly - con1(np.log(x + 1)) @ b
        return {"b": b, "D": np.mean(np.exp(res))}


def pred_P(m, te, nombre):
    x = te["P_I"].values
    if nombre == "ref":
        return x
    if nombre == "lin":
        return con1(x) @ m["b"]
    if nombre == "prop":
        return m["b"][0] * x
    if nombre == "mes":
        return x * te["mes"].map(m["f"]).values
    if nombre == "pot":
        return np.exp(con1(np.log(x + 1)) @ m["b"]) * m["D"] - 1


# ---- (ii) Q̄ a partir de P
def cols(p, K):
    return [f"{p}_{k}" for k in range(K + 1)]


def fit_Q(tr, nombre, p, K):
    y = tr["Q"].values
    climQ = tr.groupby("mes")["Q"].mean()
    if nombre == "ref":
        return {"climQ": climQ}
    if nombre in ("lin0", "linK"):
        X = tr[cols(p, 0 if nombre == "lin0" else K)].values
        return {"b": ols(con1(*X.T), y)}
    if nombre in ("logK", "logmK"):
        X = np.log(tr[cols(p, K)].values + 1)
        if nombre == "logmK":
            D = pd.get_dummies(tr["mes"]).reindex(columns=range(1, 13), fill_value=0).values.astype(float)
            A = np.column_stack([D, X])
        else:
            A = con1(*X.T)
        b = ols(A, np.log(y))
        return {"b": b, "D": np.mean(np.exp(np.log(y) - A @ b))}
    if nombre == "anoK":
        climP = tr.groupby("mes")[p].mean()
        A = np.column_stack([tr[f"{p}_{k}"].values - climP.reindex(((tr["mes"] - k - 1) % 12) + 1).values
                             for k in range(K + 1)])
        b = ols(A, y - climQ.reindex(tr["mes"]).values)
        return {"b": b, "climQ": climQ, "climP": climP}


def pred_Q(m, te, nombre, p, K):
    if nombre == "ref":
        return m["climQ"].reindex(te["mes"]).values
    if nombre in ("lin0", "linK"):
        return con1(*te[cols(p, 0 if nombre == "lin0" else K)].values.T) @ m["b"]
    if nombre in ("logK", "logmK"):
        X = np.log(te[cols(p, K)].values + 1)
        if nombre == "logmK":
            D = pd.get_dummies(te["mes"]).reindex(columns=range(1, 13), fill_value=0).values.astype(float)
            A = np.column_stack([D, X])
        else:
            A = con1(*X.T)
        return np.exp(A @ m["b"]) * m["D"]
    if nombre == "anoK":
        A = np.column_stack([te[f"{p}_{k}"].values - m["climP"].reindex(((te["mes"] - k - 1) % 12) + 1).values
                             for k in range(K + 1)])
        return m["climQ"].reindex(te["mes"]).values + A @ m["b"]


def cv_anual(tr, fit, pred, obs):
    """Validación cruzada dejando un año fuera dentro del periodo de ajuste."""
    p = pd.Series(np.nan, index=tr.index)
    for a in tr["anio"].unique():
        m = fit(tr[tr.anio != a])
        p[tr.anio == a] = pred(m, tr[tr.anio == a])
    return metricas(tr[obs], p)


# ------------------------------------------------------------------ selección y evaluación
sel, ev, coef, preds = [], [], [], {}
CAND_P = ["ref", "lin", "prop", "mes", "pot"]
CAND_Q = [("ref", 0), ("lin0", 0)] + [(n, K) for n in ("linK", "logK", "logmK", "anoK") for K in range(1, KMAX + 1)]
for part, (a_aj, a_ev) in PARTICIONES.items():
    tr, te = df[df.anio.isin(a_aj)], df[df.anio.isin(a_ev)]
    # (i) P_L a partir de P_I
    for n in CAND_P:
        cv = cv_anual(tr, lambda t, n=n: fit_P(t, n), lambda m, t, n=n: pred_P(m, t, n), "P_L")
        sel.append({"particion": part, "objetivo": "P_L ← P_I", "modelo": n, "K": 0, **{f"cv_{k}": v for k, v in cv.items()}})
    # (ii) Q̄ a partir de P_L y de P_I
    for p in ("P_L", "P_I"):
        tr = DATOS[p][DATOS[p].anio.isin(a_aj)]
        for n, K in CAND_Q:
            cv = cv_anual(tr, lambda t, n=n, K=K: fit_Q(t, n, p, K), lambda m, t, n=n, K=K: pred_Q(m, t, n, p, K), "Q")
            sel.append({"particion": part, "objetivo": f"Q̄ ← {p}", "modelo": n, "K": K,
                        **{f"cv_{k}": v for k, v in cv.items()}})
sel = pd.DataFrame(sel)
sel.round(4).to_csv(OUT / "seleccion_cv_2_2.csv", index=False)

elegidos = {}
for (part, obj), g in sel.groupby(["particion", "objetivo"]):
    g2 = g[g.modelo != "ref"]
    elegidos[(part, obj)] = g2.loc[g2.cv_MAE.idxmin(), ["modelo", "K"]].tolist()
print("Elegidos por validación cruzada en el ajuste:", elegidos)

for part, (a_aj, a_ev) in PARTICIONES.items():
    tr, te = df[df.anio.isin(a_aj)], df[df.anio.isin(a_ev)]
    obs_ok = te["flag"] == "completo"
    tareas = [("P_L ← P_I", "P_L", n, 0, None) for n in CAND_P]
    tareas += [(f"Q̄ ← {p}", "Q", n, K, p) for p in ("P_L", "P_I") for n, K in CAND_Q]
    for obj, y, n, K, p in tareas:
        d_ = DATOS[p] if p else df
        tr, te = d_[d_.anio.isin(a_aj)], d_[d_.anio.isin(a_ev)]
        obs_ok = te["flag"] == "completo"
        if obj.startswith("P_L"):
            m = fit_P(tr, n); pr, pa = pred_P(m, te, n), pred_P(m, tr, n)
        else:
            m = fit_Q(tr, n, p, K); pr, pa = pred_Q(m, te, n, p, K), pred_Q(m, tr, n, p, K)
        es_sel = elegidos[(part, obj)] == [n, K]
        for muestra, mask in (("todos", np.ones(len(te), bool)), ("solo Q observado", obs_ok.values)):
            if muestra != "todos" and y != "Q":
                continue
            met = metricas(te[y][mask], pd.Series(pr, index=te.index)[mask])
            ev.append({"particion": part, "evaluacion": f"{min(a_ev)}–{max(a_ev)}", "objetivo": obj, "modelo": n,
                       "K": K, "elegido_cv": es_sel, "muestra": muestra, **met})
        if es_sel or n == "ref" or n in ("lin", "lin0"):
            preds[(part, obj, n, K)] = (pd.Series(pa, index=tr.index), pd.Series(pr, index=te.index))
        if es_sel or n in ("lin", "lin0"):
            coef.append({"particion": part, "objetivo": obj, "modelo": n, "K": K,
                         "parametros": "; ".join(f"{k}={np.round(v, 4).tolist()}" if not isinstance(v, pd.Series)
                                                 else f"{k}=" + ", ".join(f"{MESES[i - 1]}:{x:.3f}" for i, x in v.items())
                                                 for k, v in m.items())})
ev = pd.DataFrame(ev)
# habilidad respecto a la referencia de la misma partición, objetivo y muestra
ref = ev[ev.modelo == "ref"].set_index(["particion", "objetivo", "muestra"])["MAE"]
ev["habilidad_MAE"] = 1 - ev["MAE"] / ev.set_index(["particion", "objetivo", "muestra"]).index.map(ref)
ev.round(4).to_csv(OUT / "evaluacion_2_3.csv", index=False)
pd.DataFrame(coef).to_csv(OUT / "coeficientes_2_2.csv", index=False)

pd.set_option("display.width", 220)
cols_v = ["particion", "evaluacion", "objetivo", "modelo", "K", "muestra", "n", "sesgo", "MAE", "RMSE", "NSE",
          "habilidad_MAE", "n_negativos"]
print("\nEvaluación fuera del ajuste (referencia, lineal y elegido):")
mostrar = ev[(ev.modelo == "ref") | ev.elegido_cv | ev.modelo.isin(["lin", "lin0"])]
print(mostrar[cols_v].round(3).to_string(index=False))

# ------------------------------------------------------------------ figuras
COL = {"obs": "#3d3d3d", "ref": "#9e9e9e", "lin": "#c9a227", "sel": "#1b6ca8"}


def figura(obj, y, unidades, archivo, titulo, logy=False):
    part = "A"
    a_aj, a_ev = PARTICIONES[part]
    te = df[df.anio.isin(a_ev)]
    n_sel, K_sel = elegidos[(part, obj)]
    lin = "lin" if obj.startswith("P_L") else "lin0"
    _, p_ref = preds[(part, obj, "ref", 0)]
    _, p_lin = preds[(part, obj, lin, 0)]
    p_aj, p_sel = preds[(part, obj, n_sel, K_sel)]
    tr = df[df.anio.isin(a_aj)]
    fig = plt.figure(figsize=(15, 10.5))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.1, 1], hspace=0.35, wspace=0.28)
    ax = fig.add_subplot(gs[0, :])
    ax.plot(te.index, te[y], color=COL["obs"], lw=1.6, label="Observado")
    ax.plot(te.index, p_ref, color=COL["ref"], lw=1.1, ls=":", label="Referencia: " + ("IMERG sin corregir" if y == "P_L"
                                                                                      else "climatología mensual (años de ajuste)"))
    ax.plot(te.index, p_lin, color=COL["lin"], lw=1.1, ls="--", label="Regresión lineal simple")
    ax.plot(te.index, p_sel, color=COL["sel"], lw=1.4, label=f"Modelo elegido: {n_sel}" + (f", K = {K_sel}" if K_sel else ""))
    if logy:
        ax.set_yscale("log")
    for _, r in ev[(ev.particion == part) & (ev.objetivo == obj) & (ev.muestra == "todos")
                   & ((ev.modelo == "ref") | ev.elegido_cv | (ev.modelo == lin))].iterrows():
        pass
    e_s = ev[(ev.particion == part) & (ev.objetivo == obj) & (ev.muestra == "todos") & ev.elegido_cv].iloc[0]
    e_r = ev[(ev.particion == part) & (ev.objetivo == obj) & (ev.muestra == "todos") & (ev.modelo == "ref")].iloc[0]
    e_l = ev[(ev.particion == part) & (ev.objetivo == obj) & (ev.muestra == "todos") & (ev.modelo == lin)].iloc[0]
    ax.text(0.01, 0.97, f"Evaluación {min(a_ev)}–{max(a_ev)} (n = {int(e_s.n)}):\n"
            f"referencia  MAE {e_r.MAE:.1f} · RMSE {e_r.RMSE:.1f} · sesgo {e_r.sesgo:+.1f}\n"
            f"lineal      MAE {e_l.MAE:.1f} · RMSE {e_l.RMSE:.1f} · sesgo {e_l.sesgo:+.1f}\n"
            f"elegido     MAE {e_s.MAE:.1f} · RMSE {e_s.RMSE:.1f} · sesgo {e_s.sesgo:+.1f} · "
            f"habilidad {100 * e_s.habilidad_MAE:+.0f} %", transform=ax.transAxes, va="top", fontsize=8.5,
            family="monospace", bbox=dict(fc="white", alpha=0.9, lw=0))
    ax.set_ylabel(unidades); ax.grid(alpha=0.3); ax.legend(fontsize=8.5, frameon=False, loc="upper right", ncol=2)
    ax.set_title(f"(a) Serie fuera del ajuste ({min(a_ev)}–{max(a_ev)}); ajuste {min(a_aj)}–{max(a_aj)}", loc="left", fontsize=11)
    res_te = p_sel - te[y]
    res_tr = p_aj - tr[y]
    b1 = fig.add_subplot(gs[1, 0])
    b1.scatter(tr.index, res_tr, s=8, color="#bbbbbb", label="ajuste")
    b1.scatter(te.index, res_te, s=12, color=COL["sel"], label="evaluación")
    b1.axhline(0, color="k", lw=0.8); b1.set_title("(b) Residuo (estimado − observado) en el tiempo", loc="left", fontsize=10)
    b1.set_ylabel(unidades); b1.legend(fontsize=8, frameon=False); b1.grid(alpha=0.3)
    b2 = fig.add_subplot(gs[1, 1])
    b2.scatter(p_sel, res_te, s=12, color=COL["sel"])
    b2.axhline(0, color="k", lw=0.8); b2.set_xlabel(f"Valor estimado ({unidades})")
    b2.set_title("(c) Residuo frente al valor estimado", loc="left", fontsize=10); b2.grid(alpha=0.3)
    if logy:
        b2.set_xscale("log")
    b3 = fig.add_subplot(gs[1, 2])
    b3.boxplot([res_te[te.mes == m].values for m in range(1, 13)], tick_labels=[m[0] for m in MESES],
               medianprops=dict(color=COL["sel"]), flierprops=dict(ms=3))
    b3.axhline(0, color="k", lw=0.8); b3.set_title("(d) Residuo por mes calendario (evaluación)", loc="left", fontsize=10)
    b3.grid(alpha=0.3)
    fig.suptitle(titulo, fontsize=12.5)
    fig.savefig(OUT / archivo, dpi=200, bbox_inches="tight")
    plt.close(fig)


figura("P_L ← P_I", "P_L", "mm/mes", "fig_modelo_PL_2_3.png", "Punto 2.2–2.3 – Estimar la lluvia local (P$_L$) a partir de IMERG")
figura("Q̄ ← P_L", "Q", "m³/s", "fig_modelo_Q_2_3.png",
       "Punto 2.2–2.3 – Estimar el caudal medio mensual a partir de la lluvia local (P$_L$)", logy=True)

# resumen: habilidad de todos los candidatos en ambas particiones
fig, axs = plt.subplots(1, 3, figsize=(15, 5), sharey=False)
for ax, obj in zip(axs, ["P_L ← P_I", "Q̄ ← P_L", "Q̄ ← P_I"]):
    g = ev[(ev.objetivo == obj) & (ev.muestra == "todos") & (ev.modelo != "ref")].copy()
    g["etq"] = g.modelo + np.where(g.K > 0, " K=" + g.K.astype(str), "")
    piv = g.pivot_table(index="etq", columns="particion", values="habilidad_MAE")
    piv = piv.loc[g.drop_duplicates("etq").etq]
    y = np.arange(len(piv))
    ax.barh(y - 0.2, 100 * piv["A"], height=0.38, color="#1b6ca8", label="A: evaluación 2016–2024")
    ax.barh(y + 0.2, 100 * piv["B"], height=0.38, color="#c9a227", label="B: evaluación 1998–2006")
    sels = {elegidos[(p, obj)][0] + (f" K={elegidos[(p, obj)][1]}" if elegidos[(p, obj)][1] else "") for p in "AB"}
    ax.set_yticks(y, [("★ " if e in sels else "") + e for e in piv.index], fontsize=8)
    ax.axvline(0, color="k", lw=0.8); ax.invert_yaxis(); ax.grid(axis="x", alpha=0.3)
    ax.set_xlabel("Habilidad = 1 − MAE/MAE$_{ref}$ (%)")
    ax.set_title(obj.replace("P_L", "P$_L$").replace("P_I", "P$_I$"), loc="left", fontsize=11)
axs[0].legend(fontsize=8, frameon=False, loc="lower right")
fig.suptitle("Punto 2.3 – Mejora fuera del ajuste respecto a la referencia (★ = elegido por validación cruzada en el ajuste)",
             fontsize=12)
fig.text(0.5, -0.03, "Referencias: IMERG sin corregir (para P$_L$) y climatología mensual del caudal calculada con los años de "
         "ajuste (para Q̄). Habilidad > 0: el modelo mejora a la referencia.", ha="center", fontsize=8.5, color="#444")
fig.tight_layout()
fig.savefig(OUT / "fig_evaluacion_resumen_2_3.png", dpi=200, bbox_inches="tight")
plt.close(fig)

out = []
for (part, obj, n, K), (pa, pr) in preds.items():
    out.append(pd.DataFrame({"particion": part, "objetivo": obj, "modelo": n, "K": K, "conjunto": "evaluacion",
                             "estimado": pr}))
pd.concat(out).rename_axis("fecha").round(3).to_csv(OUT / "predicciones_2_3.csv")
print(f"\nSalidas en {OUT}")
