"""
Punto 1.5.b: descripción y contraste del ciclo anual de la cuenca El Tambo en el periodo común
(meses con P_L e IMERG; periodo_comun.py). Variables: P_L, P_I, Q̄ (como R en mm/mes para comparar
formas con la lluvia) y T (PISCOt + estimación ERA5-Land desde 2021, script 17).

Indicadores (criterios explícitos; ver resumen al final):
  - Máximo, mínimo, amplitud (máx − mín de la media mensual) y amplitud relativa (amplitud / media).
  - Índice de estacionalidad SI (Walsh & Lawler 1981): SI = (1/X) Σ |x_m − X/12|, X = total anual medio.
    Clases: < 0.20 muy uniforme; 0.20–0.39 uniforme con estación más húmeda; 0.40–0.59 algo estacional
    con estación seca corta; 0.60–0.79 estacional; 0.80–0.99 marcadamente estacional con estación seca
    larga; 1.00–1.19 casi toda la lluvia en ≤ 3 meses; ≥ 1.20 extremo.
  - Índice de concentración PCI (Oliver 1980): PCI = 100 Σ x_m² / (Σ x_m)². < 10 uniforme;
    11–15 estacionalidad moderada; 16–20 estacional; > 20 fuertemente estacional.
  - Centroide circular y concentración (Markham 1970): cada mes es un vector de ángulo (m − 0.5)·30°
    y módulo x_m; el ángulo de la resultante da el mes "centro de masa" y |resultante|/Σx (0–1) la
    concentración. El desfase Q–P es la diferencia de centroides (meses).
  - Picos: máximos locales del ciclo (vecinos circulares) con media ≥ X/12 y prominencia ≥ 10 % de la
    amplitud. Armónicos: fracción de la varianza del ciclo de 12 medias explicada por el 1.º (anual)
    y el 2.º (semianual). Criterio: UNIMODAL si hay un solo pico y H1 ≥ 2·H2; BIMODAL si hay dos picos
    y H2 ≥ 0.5·H1; ESTACIONALIDAD DÉBIL si SI < 0.40.
  - Meses húmedos: media ≥ X/12. Meses secos: criterio de Bagnouls–Gaussen P [mm] < 2·T [°C].
Variabilidad por mes: desviación estándar, CV (solo P y Q; inestable si la media es pequeña: se marca
cuando la media < 10 % del máximo mensual), asimetría G1 y peso de cada año en la media (jackknife:
máximo cambio de la media al quitar un año, en %).
Picos por año: mes del máximo de cada año hidrológico (sep–ago) completo.
Subperiodos: 1998-01..2011-12 y 2012-01..fin del periodo común; mismas métricas.
Anomalías: a_t = X_t − μ_j (μ_j: media del mes j en el periodo común); totales/medias por año
hidrológico (sep–ago, estación húmeda dic–abr) y por año calendario. Los años contrastantes se eligen
por año calendario porque el año hidrológico 1997-98 (El Niño) queda incompleto en el periodo común.

Salidas (Tarea 1/El_Tambo/resultados/climatologia/):
  indicadores_estacionalidad_1_5b.csv, variabilidad_por_mes_1_5b.csv, mes_pico_por_anio_1_5b.csv,
  subperiodos_1_5b.csv, anomalias_anio_hidrologico_1_5b.csv, fig_ciclo_1_5b.png, fig_anomalias_1_5b.png
Requiere 09, 10 y 17.
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import skew
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

from periodo_comun import periodo_comun

ROOT = Path(__file__).resolve().parents[2] / "Tarea 1" / "El_Tambo"
SM = ROOT / "resultados" / "series_mensuales"
OUT = ROOT / "resultados" / "climatologia"
OUT.mkdir(parents=True, exist_ok=True)

INI, FIN = periodo_comun()
cam = pd.read_csv(SM / "series_mensuales_camels_pe.csv", index_col="fecha", parse_dates=True).loc[INI:FIN]
imerg = pd.read_csv(SM / "imerg_cuenca_mensual.csv", index_col="fecha", parse_dates=True)["P_I_mm_mes"].loc[INI:FIN]
text = pd.read_csv(SM / "temperatura_piscot_extendida_era5land.csv", index_col="fecha", parse_dates=True).loc[INI:FIN]
df = pd.DataFrame({"P_L": cam["P_L_mm_mes"], "P_I": imerg, "R": cam["R_mm_mes"], "Q": cam["Q_m3_s"],
                   "T": text["T_piscot_C"].fillna(text["T_estimada_era5land_C"])})
COL = {"P_L": "#1b6ca8", "P_I": "#e4572e", "R": "#5b3f8c", "Q": "#0b7a75", "T": "#c8553d"}
ETQ = {"P_L": "P$_L$", "P_I": "P$_I$", "R": "R ($\\bar{Q}$ en lámina)", "Q": "$\\bar{Q}$", "T": "T"}
MESES = ["Ene", "Feb", "Mar", "Abr", "May", "Jun", "Jul", "Ago", "Sep", "Oct", "Nov", "Dic"]
INICIO_AH = 9                                              # año hidrológico sep–ago
SUB = {"1998–2011": (INI, pd.Timestamp("2011-12-01")), f"2012–{FIN:%Y-%m}": (pd.Timestamp("2012-01-01"), FIN)}


# ------------------------------------------------------------------ funciones de indicadores
def clim_media(s):
    return s.groupby(s.index.month).mean().reindex(range(1, 13))


def centroide(c):
    """Mes centroide (1–12, decimal; 1.0 = inicio de ene) y concentración de Markham."""
    ang = np.deg2rad((np.arange(1, 13) - 0.5) * 30)
    z = np.sum(c.values * np.exp(1j * ang))
    mes = (np.rad2deg(np.angle(z)) % 360) / 30 + 0.5
    return mes, np.abs(z) / c.sum()


def armonicos(c):
    """Fracción de la varianza del ciclo de 12 medias explicada por el 1.º y 2.º armónico."""
    f = np.fft.rfft(c.values - c.mean())
    pot = np.abs(f[1:]) ** 2
    pot[-1] /= 2                                           # el armónico de Nyquist (k = 6) no se duplica
    return pot[0] / pot.sum(), pot[1] / pot.sum()


def picos(c, umbral):
    amp = c.max() - c.min()
    v = c.values
    out = []
    for i in range(12):
        izq, der = v[i - 1], v[(i + 1) % 12]
        if v[i] > izq and v[i] >= der and v[i] >= umbral:
            # prominencia: altura sobre el mínimo más alto de los dos lados (recorrido circular)
            lados = []
            for paso in (-1, 1):
                j, mn = i, v[i]
                while True:
                    j = (j + paso) % 12
                    mn = min(mn, v[j])
                    if v[j] > v[i] or j == i:
                        break
                lados.append(mn)
            if v[i] - max(lados) >= 0.1 * amp:
                out.append(i + 1)
    return out


def indicadores(sub, nombre):
    filas = []
    cT = clim_media(sub["T"])
    for k in ("P_L", "P_I", "R", "T"):
        c = clim_media(sub[k])
        fila = {"muestra": nombre, "variable": k, "mes_max": MESES[c.idxmax() - 1], "max": c.max(),
                "mes_min": MESES[c.idxmin() - 1], "min": c.min(), "amplitud": c.max() - c.min()}
        h1, h2 = armonicos(c)
        fila.update({"var_H1_pct": 100 * h1, "var_H2_pct": 100 * h2})
        if k == "T":
            filas.append(fila)
            continue
        X = c.sum()
        cm, conc = centroide(c)
        pk = picos(c, X / 12)
        si = np.abs(c - X / 12).sum() / X
        humedos = [MESES[m - 1] for m in range(1, 13) if c[m] >= X / 12]
        fila.update({"total_anual_medio": X, "amplitud_relativa": (c.max() - c.min()) / (X / 12),
                     "SI": si, "PCI": 100 * (c ** 2).sum() / X ** 2, "centroide_mes": cm,
                     "centroide": f"{MESES[int(cm) % 12 - 1 if int(cm) % 12 else 11]} ({cm:.2f})",
                     "concentracion": conc, "n_picos": len(pk), "picos": ", ".join(MESES[m - 1] for m in pk),
                     "meses_humedos": f"{len(humedos)}: {', '.join(humedos)}",
                     "pct_en_dic_abr": 100 * c.loc[[12, 1, 2, 3, 4]].sum() / X})
        if k in ("P_L", "P_I"):
            secos = [MESES[m - 1] for m in range(1, 13) if c[m] < 2 * cT[m]]
            fila["meses_secos_gaussen"] = f"{len(secos)}: {', '.join(secos)}"
        fila["regimen"] = ("estacionalidad débil" if si < 0.40 else
                           "unimodal" if len(pk) == 1 and h1 >= 2 * h2 else
                           "bimodal" if len(pk) == 2 and h2 >= 0.5 * h1 else "indeterminado")
        filas.append(fila)
    return filas


# ------------------------------------------------------------------ 1. indicadores del periodo común y subperiodos
ind = pd.DataFrame(indicadores(df, f"periodo común {INI:%Y-%m} a {FIN:%Y-%m}")
                   + sum((indicadores(df.loc[a:b], f"subperiodo {n}") for n, (a, b) in SUB.items()), []))
ind.round(3).to_csv(OUT / "indicadores_estacionalidad_1_5b.csv", index=False)
pd.set_option("display.width", 260, "display.max_columns", 40, "display.max_colwidth", 60)
cols = ["muestra", "variable", "mes_max", "max", "mes_min", "min", "amplitud", "SI", "PCI", "centroide",
        "concentracion", "n_picos", "picos", "var_H1_pct", "var_H2_pct", "meses_humedos", "meses_secos_gaussen",
        "pct_en_dic_abr", "regimen"]
print(f"=== Indicadores de estacionalidad (periodo común {INI:%Y-%m} a {FIN:%Y-%m} y subperiodos) ===")
print(ind[cols].round(2).to_string(index=False))

# ------------------------------------------------------------------ 2. desfases y fase
base = ind[ind.muestra.str.startswith("periodo")].set_index("variable")
lag = {}
for p in ("P_L", "P_I"):
    d_c = (base.loc["R", "centroide_mes"] - base.loc[p, "centroide_mes"]) % 12
    d_pk = (MESES.index(base.loc["R", "mes_max"]) - MESES.index(base.loc[p, "mes_max"])) % 12
    d_min = (MESES.index(base.loc["R", "mes_min"]) - MESES.index(base.loc[p, "mes_min"])) % 12
    lag[p] = (d_c, d_pk, d_min)
    print(f"\nDesfase Q–{p}: centroide {d_c:.2f} meses | mes del máximo {d_pk} | mes del mínimo {d_min}")
# Correlación entre ciclos normalizados con desplazamiento (forma y fase)
for p in ("P_L", "P_I"):
    cp, cr = clim_media(df[p]), clim_media(df["R"])
    r = {s: np.corrcoef(np.roll(cp.values, s), cr.values)[0, 1] for s in range(0, 4)}
    print(f"r(ciclo R, ciclo {p} desplazado s meses): " + ", ".join(f"s={s}: {v:.2f}" for s, v in r.items()))

# ------------------------------------------------------------------ 3. variabilidad y asimetría por mes, peso de cada año
filas = []
for k in ("P_L", "P_I", "Q", "R", "T"):
    cmax = clim_media(df[k]).max()
    for mes in range(1, 13):
        s = df.loc[df.index.month == mes, k].dropna()
        jk = pd.Series({y: s.drop(s.index[s.index.year == y]).mean() for y in s.index.year})
        cambio = (jk - s.mean()) / s.mean() * 100
        fila = {"variable": k, "mes": MESES[mes - 1], "n": len(s), "media": s.mean(), "desv_est": s.std(ddof=1),
                "asimetria_G1": skew(s, bias=False),
                "jackknife_max_cambio_pct": cambio.abs().max(), "anio_mas_influyente": int(cambio.abs().idxmax())}
        if k != "T":                                         # CV no se usa en °C (cero convencional)
            fila["CV"] = s.std(ddof=1) / s.mean()
            fila["CV_inestable"] = s.mean() < 0.1 * cmax
        filas.append(fila)
var = pd.DataFrame(filas)
var.round(3).to_csv(OUT / "variabilidad_por_mes_1_5b.csv", index=False)
print("\n=== Variabilidad por mes (CV, asimetría, peso del año más influyente) ===")
for k in ("P_L", "P_I", "Q"):
    v = var[var.variable == k].set_index("mes")
    print(f"{k}: CV " + " ".join(f"{m}:{c:.2f}{'*' if u else ''}" for m, c, u in zip(v.index, v.CV, v.CV_inestable)))
    print(f"    G1 " + " ".join(f"{m}:{g:.1f}" for m, g in zip(v.index, v.asimetria_G1)))
    print(f"    jackknife máx " + " ".join(f"{m}:{j:.0f}%({a})" for m, j, a in
                                          zip(v.index, v.jackknife_max_cambio_pct, v.anio_mas_influyente)))
print("(* = CV inestable: media del mes < 10 % del máximo mensual)")

# ------------------------------------------------------------------ 4. ¿aparece el pico medio en la mayoría de los años?
ah = df.copy()
ah["AH"] = np.where(ah.index.month >= INICIO_AH, ah.index.year + 1, ah.index.year)      # año hidrológico (año final)
completos = ah.groupby("AH").size()
ah = ah[ah["AH"].isin(completos[completos == 12].index)]
pico = []
for y, g in ah.groupby("AH"):
    fila = {"anio_hidrologico": f"{y - 1}-{str(y)[2:]}"}
    for k in ("P_L", "P_I", "Q"):
        s = g[k]
        fila[f"mes_max_{k}"] = MESES[s.idxmax().month - 1] if s.notna().all() else None
    pico.append(fila)
pico = pd.DataFrame(pico)
pico.to_csv(OUT / "mes_pico_por_anio_1_5b.csv", index=False)
print(f"\n=== Mes del máximo por año hidrológico (sep–ago), {len(pico)} años completos ===")
frec = {}
for k in ("P_L", "P_I", "Q"):
    s = pico[f"mes_max_{k}"].dropna()
    frec[k] = s.value_counts().reindex(MESES, fill_value=0)
    pk_medio = base.loc["R" if k == "Q" else k, "mes_max"]
    feb_abr = s.isin(["Feb", "Mar", "Abr"]).mean() * 100
    print(f"{k}: n = {len(s)} | máximo en {pk_medio} (mes del pico medio) en {(s == pk_medio).mean() * 100:.0f} % "
          f"de los años | en feb–abr en {feb_abr:.0f} % | {dict(frec[k][frec[k] > 0])}")

# ------------------------------------------------------------------ 5. subperiodos (tabla de medias mensuales)
sub_rows = []
for n, (a, b) in SUB.items():
    for k in ("P_L", "P_I", "Q", "T"):
        c = clim_media(df.loc[a:b, k])
        sub_rows.append({"subperiodo": n, "variable": k, **{MESES[m - 1]: c[m] for m in range(1, 13)},
                         "anual": c.sum() if k != "T" and k != "Q" else c.mean()})
subt = pd.DataFrame(sub_rows)
subt.round(2).to_csv(OUT / "subperiodos_1_5b.csv", index=False)
print("\n=== Medias mensuales por subperiodo (anual: suma para P; media para Q y T) ===")
print(subt.round(1).to_string(index=False))

# ------------------------------------------------------------------ 6. anomalías respecto al ciclo anual
mu = {k: clim_media(df[k]) for k in df}
anom = pd.DataFrame({k: df[k] - df.index.month.map(mu[k]).values for k in df}, index=df.index)
anom["AH"] = np.where(anom.index.month >= INICIO_AH, anom.index.year + 1, anom.index.year)
anom["humeda"] = anom.index.month.isin([12, 1, 2, 3, 4])
completo_ah = anom.groupby("AH").size() == 12
filas = []
for y in completo_ah[completo_ah].index:
    g = anom[anom["AH"] == y]
    gh = g[g["humeda"]]
    fila = {"anio_hidrologico": f"{y - 1}-{str(y)[2:]}"}
    for k in ("P_L", "P_I", "R"):
        fila[f"{k}_anom_total_mm"] = g[k].sum() if g[k].notna().all() else np.nan
        fila[f"{k}_anom_dic_abr_mm"] = gh[k].sum() if gh[k].notna().all() else np.nan
    fila["T_anom_media_C"] = g["T"].mean()
    filas.append(fila)
aah = pd.DataFrame(filas)
aah.round(2).to_csv(OUT / "anomalias_anio_hidrologico_1_5b.csv", index=False)
print(f"\nAños hidrológicos completos en el periodo común: {len(aah)} ({aah.anio_hidrologico.iloc[0]} a "
      f"{aah.anio_hidrologico.iloc[-1]}); el año 1997-98 (El Niño) queda incompleto porque empieza en 1997-09.")
# Años calendario completos (incluyen 1998): base para elegir los años contrastantes
ac = anom.assign(anio=anom.index.year)
comp = ac.groupby("anio").size() == 12
filas = []
for y in comp[comp].index:
    g = ac[ac.anio == y]
    fila = {"anio": y}
    for k in ("P_L", "P_I", "R"):
        fila[f"{k}_anom_total_mm"] = g[k].sum() if g[k].notna().all() else np.nan
        fila[f"{k}_anom_ene_abr_mm"] = g.loc[g.index.month <= 4, k].sum() if g.loc[g.index.month <= 4, k].notna().all() else np.nan
    fila["T_anom_media_C"] = g["T"].mean()
    filas.append(fila)
acal = pd.DataFrame(filas)
acal.round(2).to_csv(OUT / "anomalias_anio_calendario_1_5b.csv", index=False)
o = acal.sort_values("P_L_anom_total_mm")
print(f"\n=== Años calendario contrastantes ({len(acal)} años completos; anomalía total de P_L, mm) ===")
print("Más secos:\n", o.head(4).round(1).to_string(index=False))
print("Más húmedos:\n", o.tail(4).round(1).to_string(index=False))
print(f"r entre anomalías anuales (año calendario): P_L–P_I {acal['P_L_anom_total_mm'].corr(acal['P_I_anom_total_mm']):.2f}"
      f" | P_L–R {acal['P_L_anom_total_mm'].corr(acal['R_anom_total_mm']):.2f} (n = {acal['R_anom_total_mm'].notna().sum()}) | "
      f"ene–abr: P_L–R {acal['P_L_anom_ene_abr_mm'].corr(acal['R_anom_ene_abr_mm']):.2f} "
      f"(n = {acal['R_anom_ene_abr_mm'].notna().sum()}), P_I–R {acal['P_I_anom_ene_abr_mm'].corr(acal['R_anom_ene_abr_mm']):.2f}")
signo = (np.sign(acal["P_L_anom_total_mm"]) != np.sign(acal["P_I_anom_total_mm"]))
print("Años con signo opuesto de anomalía P_L vs P_I:", acal.loc[signo, "anio"].tolist())

# ------------------------------------------------------------------ figura 1: forma, fase, picos por año y subperiodos
x = np.arange(1, 13)
fig, axs = plt.subplots(2, 2, figsize=(15, 10))
fig.subplots_adjust(hspace=0.36, wspace=0.18, top=0.87, bottom=0.08, left=0.06, right=0.98)
ax = axs[0, 0]
for k in ("P_L", "P_I", "R"):
    c = clim_media(df[k])
    ax.plot(x, 100 * c / c.sum(), color=COL[k], lw=2, marker="o", ms=4,
            label=f"{ETQ[k]}: centroide {base.loc[k, 'centroide_mes']:.1f} · SI {base.loc[k, 'SI']:.2f} · "
                  f"PCI {base.loc[k, 'PCI']:.0f}")
    ax.axvline(base.loc[k, "centroide_mes"], color=COL[k], lw=1, ls=":")
ax.axhline(100 / 12, color="k", lw=0.8, ls="--")
ax.text(12.3, 100 / 12, "1/12", fontsize=8, va="center")
ax.set_title("(a) Forma y fase: % del total anual medio en cada mes\n"
             f"desfase del centroide R–P$_L$ = {lag['P_L'][0]:.1f} meses, R–P$_I$ = {lag['P_I'][0]:.1f} meses "
             "(líneas punteadas: centroides)", loc="left", fontweight="bold", fontsize=10)
ax.set_ylabel("% del total anual")
ax.legend(frameon=False, fontsize=8.5)
ax.set_xticks(x, MESES)
ax = axs[0, 1]
w = 0.27
for i, k in enumerate(("P_L", "P_I", "Q")):
    f = frec[k]
    ax.bar(x + (i - 1) * w, 100 * f.values / f.sum(), width=w, color=COL[k], label=f"{ETQ[k]} (n = {f.sum()})")
ax.set_xticks(x, MESES)
ax.set_ylabel("% de años hidrológicos")
ax.set_title("(b) Mes del máximo de cada año hidrológico (sep–ago)", loc="left", fontweight="bold", fontsize=10.5)
ax.legend(frameon=False, fontsize=8.5)
ax = axs[1, 0]
for (n, (a, b)), ls in zip(SUB.items(), ("-", "--")):
    for k in ("P_L", "P_I", "R"):
        c = clim_media(df.loc[a:b, k])
        ax.plot(x, c, color=COL[k], lw=1.8, ls=ls, label=f"{ETQ[k]} {n}")
ax.set_xticks(x, MESES)
ax.set_ylabel("mm/mes")
ax.set_title("(c) Ciclo medio por subperiodo", loc="left", fontweight="bold", fontsize=10.5)
ax.legend(frameon=False, fontsize=8, ncol=2)
ax = axs[1, 1]
for k in ("P_L", "P_I", "Q"):
    v = var[var.variable == k]
    ax.plot(x, v["CV"], color=COL[k], lw=1.8, marker="o", ms=4, label=ETQ[k])
    inst = v["CV_inestable"].astype(bool).values
    ax.plot(x[inst], v["CV"].values[inst], ls="none", marker="o", ms=9, mfc="none", mec=COL[k])
ax.set_xticks(x, MESES)
ax.set_ylabel("CV (desv. est. / media)")
ax.set_title("(d) Variabilidad entre años por mes\n(círculo: CV inestable, media del mes < 10 % del máximo)",
             loc="left", fontweight="bold", fontsize=10)
ax.legend(frameon=False, fontsize=8.5)
for a_ in axs.flat:
    a_.grid(alpha=0.3)
fig.suptitle(f"Cuenca del río Chicama hasta El Tambo – rasgos del ciclo anual, periodo común {INI:%Y-%m} a {FIN:%Y-%m}",
             fontsize=13, y=0.975)
fig.text(0.5, 0.94, "SI: índice de estacionalidad (Walsh & Lawler 1981); PCI: índice de concentración (Oliver 1980); "
         "centroide y concentración: Markham (1970). R = caudal en lámina (mm/mes).", ha="center", fontsize=9)
fig.savefig(OUT / "fig_ciclo_1_5b.png", dpi=200)
plt.close(fig)

# ------------------------------------------------------------------ figura 2: anomalías
fig, axs = plt.subplots(4, 1, figsize=(15, 11), sharex=True, gridspec_kw={"hspace": 0.15})
for ax, k, unid in zip(axs, ("P_L", "P_I", "R", "T"), ("mm/mes", "mm/mes", "mm/mes", "°C")):
    a = anom[k]
    ax.bar(a.index, a.values, width=25, color=np.where(a >= 0, COL[k], "#9e9e9e"), lw=0)
    ax.plot(a.index, a.rolling(12, center=True, min_periods=10).mean(), color="k", lw=1.4)
    ax.axhline(0, color="k", lw=0.6)
    ax.set_ylabel(unid)
    ax.set_title(f"{ETQ[k]}: anomalía respecto a la media del mes calendario (barras) y media móvil de 12 meses",
                 loc="left", fontsize=10)
    ax.grid(alpha=0.3)
    if k == "T":
        ax.axvspan(pd.Timestamp("2021-01-01"), FIN + pd.offsets.MonthBegin(1), color="#6a4c93", alpha=0.1, lw=0)
        ax.text(pd.Timestamp("2021-02-01"), ax.get_ylim()[1] * 0.85, "T estimada (ERA5-Land + sesgo)", fontsize=8,
                color="#6a4c93")
for y in o.tail(3)["anio"].tolist() + o.head(3)["anio"].tolist():
    y0 = pd.Timestamp(f"{y}-01-01")
    for ax in axs:
        ax.axvspan(y0, y0 + pd.DateOffset(years=1), color="#fde68a", alpha=0.25, lw=0, zorder=0)
    axs[0].text(y0 + pd.Timedelta(days=182), axs[0].get_ylim()[1] * 0.92, str(y), ha="center", fontsize=7.5)
axs[-1].xaxis.set_major_locator(mdates.YearLocator(2))
axs[-1].xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
fig.suptitle(f"Anomalías mensuales respecto al ciclo anual (referencia: medias por mes del periodo común "
             f"{INI:%Y-%m} a {FIN:%Y-%m}); sombreado amarillo: 3 años calendario más húmedos y 3 más secos según P$_L$",
             fontsize=11, y=0.93)
fig.savefig(OUT / "fig_anomalias_1_5b.png", dpi=200, bbox_inches="tight")
print(f"\nSalidas en {OUT}")
