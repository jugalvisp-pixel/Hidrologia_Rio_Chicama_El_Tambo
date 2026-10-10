"""
Validación del método de relleno del caudal (relleno_caudal.py, aplicado en el script 09).

1. Huecos diarios sintéticos: en meses observados completos (1998-01 en adelante) se borran L días seguidos
   (L = 1–15; 600 repeticiones por L, semilla fija) y se compara el total mensual reconstruido con el
   observado, con dos métodos: interpolación lineal de ln(q) (el elegido) y el escalado anterior
   (media de los días válidos · n_días).
2. Modelo mensual: validación cruzada dejando un año fuera sobre los meses observados completos, para
   cinco variantes de regresión log-lineal (con dummies de mes) y la razón obs/sim por mes calendario.
   Se elige la de menor error con predictores disponibles en todos los meses a rellenar.
3. Resumen de qué se rellenó en el periodo común y con qué método.

Salidas (Tarea 1/El_Tambo/resultados/series_mensuales/):
  validacion_relleno_huecos_diarios.csv, validacion_relleno_modelo_mensual.csv,
  resumen_relleno_caudal.csv, fig_validacion_relleno_caudal.png
"""
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import relleno_caudal as rc
from periodo_comun import periodo_comun

ROOT = Path(__file__).resolve().parents[2] / "Tarea 1" / "El_Tambo"
OUT = ROOT / "resultados" / "series_mensuales"
INI, FIN = periodo_comun()
rng = np.random.default_rng(42)
d = pd.read_csv(ROOT / "datos" / "camels_pe" / "PE_47E0A2A8.csv", parse_dates=["date"], index_col="date")
q = d["flow_obs"].loc[d["flow_obs"].first_valid_index():]

# ------------------------------------------------------------------ 1. huecos diarios sintéticos
meses_ok = [p for p, g in q.groupby(q.index.to_period("M")) if g.notna().all()]
filas = []
for L in (1, 2, 3, 4, 5, 7, 10, 12, 15):
    e_i, e_s = [], []
    while len(e_i) < 600:
        g = q[str(meses_ok[rng.integers(len(meses_ok))])]
        i0 = rng.integers(1, len(g) - L - 1)
        s = g.copy()
        s.iloc[i0:i0 + L] = np.nan
        e_i.append(100 * (rc.interpolar_log(s, L)[0].sum() / g.sum() - 1))
        e_s.append(100 * (s.mean() * len(g) / g.sum() - 1))
    for metodo, e in (("interpolacion_ln_q", np.array(e_i)), ("escalado", np.array(e_s))):
        filas.append({"dias_borrados": L, "metodo": metodo, "error_abs_mediano_pct": np.median(np.abs(e)),
                      "error_abs_P90_pct": np.percentile(np.abs(e), 90), "sesgo_medio_pct": e.mean()})
hue = pd.DataFrame(filas).round(2)
hue.to_csv(OUT / "validacion_relleno_huecos_diarios.csv", index=False)
print(hue.to_string(index=False))

# ------------------------------------------------------------------ 2. modelo mensual (dejando un año fuera)
m = pd.read_csv(OUT / "series_mensuales_camels_pe.csv", index_col="fecha", parse_dates=True)
g = d.resample("MS")
t = pd.DataFrame({"R": m["R_obs_mm_mes"].where(m["R_flag_obs"] == "completo"), "S": g["flow_sim"].sum(min_count=1),
                  "P": m["P_L_mm_mes"]}).loc["1998":]
t["R1"], t["P1"] = t["R"].shift(1), t["P"].shift(1)
VARS = {"P, P(m−1)": ["P", "P1"], "S": ["S"], "S, R(m−1)": ["S", "R1"], "P, P(m−1), R(m−1)": ["P", "P1", "R1"],
        "S, P, R(m−1)": ["S", "P", "R1"]}


def X(tt, vars_):
    meses = pd.get_dummies(tt.index.month).reindex(columns=range(1, 13), fill_value=0).values.astype(float)
    return np.column_stack([meses, *[np.log(tt[c].clip(lower=0.01)).values for c in vars_]])


filas, preds = [], {}
obs = t.dropna(subset=["R"])
for nombre, vars_ in VARS.items():
    pr = []
    for a in obs.index.year.unique():
        tr, te = obs[obs.index.year != a].dropna(subset=vars_), obs[obs.index.year == a].dropna(subset=vars_)
        if te.empty:
            continue
        y = np.log(tr["R"].clip(lower=0.01))
        b, *_ = np.linalg.lstsq(X(tr, vars_), y, rcond=None)
        s2 = np.var(y - X(tr, vars_) @ b)
        pr.append(pd.Series(np.exp(X(te, vars_) @ b + s2 / 2), index=te.index))
    pr = pd.concat(pr)
    o = obs["R"].reindex(pr.index)
    rel = (100 * (pr / o - 1)).abs()
    filas.append({"modelo": f"ln R = a_mes + {nombre}", "n": len(pr),
                  "NSE": 1 - ((pr - o) ** 2).sum() / ((o - o.mean()) ** 2).sum(),
                  "NSE_log": 1 - ((np.log(pr) - np.log(o)) ** 2).sum() / ((np.log(o) - np.log(o).mean()) ** 2).sum(),
                  "error_rel_mediano_pct": rel.median(), "sesgo_pct": 100 * (pr.sum() / o.sum() - 1)})
    preds[nombre] = (o, pr)
mod = pd.DataFrame(filas).round(3)
mod.to_csv(OUT / "validacion_relleno_modelo_mensual.csv", index=False)
print("\n", mod.to_string(index=False))

# ------------------------------------------------------------------ 3. resumen de lo rellenado
pc = m.loc[INI:FIN]
res = pd.crosstab(pc["R_flag_obs"], pc["R_flag"])
res.to_csv(OUT / "resumen_relleno_caudal.csv")
print("\nPeriodo común: condición observada (filas) vs método final (columnas)\n", res)
cambio = pc[pc["R_flag_obs"] != "completo"]
print(f"Meses con algún relleno: {len(cambio)} | días interpolados: {int(pc['dias_interpolados_Q'].sum())} | "
      f"volumen de R en meses antes excluidos: {pc.loc[pc['R_flag_obs'].isin(['excluido', 'sin_dato']), 'R_mm_mes'].sum():.0f} mm "
      f"de {pc['R_mm_mes'].sum():.0f} mm ({100 * pc.loc[pc['R_flag_obs'].isin(['excluido', 'sin_dato']), 'R_mm_mes'].sum() / pc['R_mm_mes'].sum():.1f} %)")

# ------------------------------------------------------------------ figura
fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 5.2))
for metodo, col, etq in (("interpolacion_ln_q", "#1b6ca8", "Interpolación de ln(q) (método elegido)"),
                         ("escalado", "#e4572e", "Escalado: media de días válidos × n$_m$ (método anterior)")):
    h = hue[hue.metodo == metodo]
    a1.plot(h.dias_borrados, h.error_abs_mediano_pct, "-o", color=col, label=f"{etq} – mediana")
    a1.plot(h.dias_borrados, h.error_abs_P90_pct, "--", color=col, alpha=0.7, label="  percentil 90")
a1.set_xlabel("Días seguidos faltantes en el mes")
a1.set_ylabel("Error absoluto del total mensual (%)")
a1.set_title("(a) Huecos diarios sintéticos en meses completos", loc="left", fontsize=10.5)
a1.grid(alpha=0.3); a1.legend(fontsize=8, frameon=False)
o, pr = preds["S, R(m−1)"]
a2.loglog(o, pr, "o", ms=3.5, color="#5b3f8c", alpha=0.6)
lim = [min(o.min(), pr.min()) * 0.8, max(o.max(), pr.max()) * 1.2]
a2.plot(lim, lim, color="#555", lw=1)
fila = mod[mod.modelo.str.endswith("S, R(m−1)")].iloc[0]
a2.set_title("(b) Modelo mensual ln R = a$_{mes}$ + b·ln S + c·ln R$_{m-1}$\n(validación dejando un año fuera)",
             loc="left", fontsize=10.5)
a2.text(0.03, 0.97, f"NSE {fila.NSE:.2f} · NSE(log) {fila.NSE_log:.2f}\nerror rel. mediano {fila.error_rel_mediano_pct:.0f} %"
        f" · n = {int(fila.n)}", transform=a2.transAxes, va="top", fontsize=9)
a2.set_xlabel("R observado (mm/mes)"); a2.set_ylabel("R predicho (mm/mes)")
a2.grid(alpha=0.3, which="both")
fig.suptitle("Validación del relleno de caudal – cuenca El Tambo (CAMELS-PE flow_obs; S = flow_sim, PISCO-ARNOVIC)",
             fontsize=12)
fig.tight_layout()
fig.savefig(OUT / "fig_validacion_relleno_caudal.png", dpi=200)
plt.close(fig)
