"""
Punto 4: análisis de frecuencias mediante Fourier de las series mensuales de la cuenca El Tambo.

Series (secuencia mensual cronológica, Δt = 1 mes; nunca las 12 medias climatológicas):
  P_L      CAMELS-PE `prec` (PISCOp v2.1)   1981-01 a 2025-12  N = 540, sin vacíos
  P_I      IMERG V07 Final                  1998-01 a 2025-08  N = 332, sin vacíos
  Q        caudal El Tambo, m³/s            1998-01 a 2025-12  N = 336, continuo gracias al relleno documentado
                                            (80 meses: 76 interpolados ≤ 15 d, 3 escalados, 1 modelo). Efecto de los
                                            vacíos evaluado con Lomb–Scargle sobre los 256 meses observados completos
                                            (estimador para muestreo irregular; no se concatena ni se rellena con ceros).
  T_PISCOt CAMELS-PE `tmean`                1981-01 a 2020-12  N = 480
  T_ERA5L  ERA5-Land T2m (reanálisis)       1981-01 a 2025-12  N = 540
Periodos: registro completo de cada variable (distinta longitud → distinta resolución) y periodo común
1998-01 a 2024-12 (N = 324 = 27 años; T_PISCOt y su contraste con ERA5-Land: 1998–2020, N = 276).

Transformaciones (se reporta cuál se aplica a cada espectro):
  centrada            x − media global (conserva el ciclo anual)
  anomalia            x − μ_j con μ_j del periodo de referencia fijo 1998–2020 (como en el punto 3); luego se
                      resta su media global (solo afecta f = 0)
  *_sin_tendencia     además se resta la recta OLS (tendencia del punto 3)
Estimadores:
  Periodograma (FFT, ventana rectangular, sin relleno con ceros): f_k = k/(NΔt), k = 1..⌊N/2⌋, Δf = 1/(NΔt),
    T_k = 1/f_k. Densidad espectral unilateral S_k = 2Δt|X_k|²/N (sin duplicar en Nyquist), de modo que
    Σ S_k Δf = varianza (Parseval; se verifica). Unidades: (unidad)²/(ciclo/mes).
  Periodograma con ventana de Hann (reduce la fuga espectral de los picos fuertes, como el anual, a costa de
    ensanchar el lóbulo principal ~2 bins); normalizado por la potencia de la ventana.
  Welch: segmentos de 120 meses (10 años), Hann, 50 % de traslape → Δf = 1/120 ciclos/mes; pierde periodos
    > 10 años y gana estabilidad (grados de libertad ν calculados con la correlación entre segmentos).
Fondo de ruido: AR(1) ajustado a la serie (ρ = autocorrelación de rezago 1, varianza igual a la de la serie):
  S_AR(f) = σ²(1 − ρ²)/(1 + ρ² − 2ρ cos 2πf)·(2Δt). Ordenadas del periodograma ~ S_AR·χ²₂/2. Umbral puntual
  95 % y umbral global (búsqueda entre M = ⌊N/2⌋ frecuencias, Bonferroni: 1 − 0.05/M). Solo se aplica a las
  anomalías: en las series centradas el ciclo anual es determinístico y no se prueba contra ruido rojo.
Estabilidad (4.2.4): ventana (rectangular, Hann, Welch), tendencia, extremos (anomalías recortadas en los
  percentiles 1 y 99), subperiodo (mitades del registro), registro completo frente a periodo común, fuente.

Salidas: Tarea 1/El_Tambo/resultados/espectros/
"""
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import signal, stats

from periodo_comun import periodo_comun
from estilo_figuras import VAR, TINTA, TINTA2, GRIS_CLARO, colores_ordinales, aplicar_estilo

ROOT = Path(__file__).resolve().parents[2] / "Tarea 1" / "El_Tambo"
SM = ROOT / "resultados" / "series_mensuales"
OUT = ROOT / "resultados" / "espectros"
OUT.mkdir(parents=True, exist_ok=True)
PC_INI, PC_FIN = periodo_comun()
REF = ("1998-01-01", "2020-12-01")
ALFA = 0.05
SEG = 120
BANDAS = [  # (nombre, periodo mínimo, periodo máximo) en meses
    ("intraestacional 2–5 m", 2, 5.5), ("semianual 5.5–7 m", 5.5, 7), ("7–11 m", 7, 11),
    ("anual 11–13 m", 11, 13.5), ("13 m–2 años", 13.5, 24), ("interanual 2–7 años", 24, 84),
    ("> 7 años", 84, np.inf)]
PERIODOS_A = [0.25, 0.5, 1, 2, 5, 10, 20, 50]                        # marcas del eje de periodo (años)
FRECUENCIAS = [1 / 240, 1 / 60, 1 / 24, 1 / 12, 1 / 6, 1 / 3, 0.5]   # marcas del eje de frecuencia (ciclos/mes)

cam = pd.read_csv(SM / "series_mensuales_camels_pe.csv", index_col="fecha", parse_dates=True)
imerg = pd.read_csv(SM / "imerg_cuenca_mensual.csv", index_col="fecha", parse_dates=True)["P_I_mm_mes"]
t5l = pd.read_csv(SM / "era5land_t2m_cuenca_mensual.csv", index_col="fecha", parse_dates=True)["T_era5land_C"]
SERIES = {
    "P_L": (cam["P_L_mm_mes"], "P_L (CAMELS-PE)", "mm/mes", VAR["P_L"]),     # colores comunes (estilo_figuras)
    "P_I": (imerg, "P_I (IMERG)", "mm/mes", VAR["P_I"]),
    "Q": (cam["Q_m3_s"], "Q̄ El Tambo", "m³/s", VAR["Q"]),
    "T_PISCOt": (cam["T_C"], "T (PISCOt)", "°C", VAR["T"]),
    "T_ERA5L": (t5l, "T2m (ERA5-Land)", "°C", VAR["T_ERA5L"]),
}
VARS = list(SERIES)
SER = {k: v[0].dropna().sort_index() for k, v in SERIES.items()}
for k, s in SER.items():                       # muestreo regular: sin meses faltantes dentro del tramo
    full = pd.date_range(s.index.min(), s.index.max(), freq="MS")
    assert len(full) == len(s), f"{k} tiene vacíos internos: no se puede usar la FFT directa"
Q_OBS = cam["Q_m3_s"].where(cam["R_flag"] == "completo").dropna()


# ------------------------------------------------------------------ transformaciones
def anomalia(s, k):
    ref = SER[k].loc[REF[0]:REF[1]]
    mu = ref.groupby(ref.index.month).mean()
    a = s - s.index.month.map(mu).values
    return a - a.mean()


def sin_tendencia(s):
    t = np.arange(len(s), dtype=float)
    b = np.polyfit(t, s.values, 1)
    return pd.Series(s.values - np.polyval(b, t), s.index)


def transformar(k, s, tr):
    if tr.startswith("centrada"):
        x = s - s.mean()
    else:
        x = anomalia(s, k)
    if tr.endswith("sin_tendencia"):
        x = sin_tendencia(x)
    return x


# ------------------------------------------------------------------ estimadores
def periodograma(x, ventana="rect"):
    x = np.asarray(x, float)
    N = len(x)
    w = np.ones(N) if ventana == "rect" else signal.windows.hann(N, sym=False)
    U = np.mean(w ** 2)
    X = np.fft.rfft(x * w)
    f = np.fft.rfftfreq(N, d=1.0)
    S = 2.0 * np.abs(X) ** 2 / (N * U)
    if N % 2 == 0:
        S[-1] /= 2
    return f[1:], S[1:]                         # se descarta f = 0 (media)


def welch(x):
    f, S = signal.welch(np.asarray(x, float), fs=1.0, window="hann", nperseg=SEG, noverlap=SEG // 2,
                        detrend=False, scaling="density", return_onesided=True)
    paso = SEG // 2
    K = (len(x) - SEG) // paso + 1
    w = signal.windows.hann(SEG, sym=False)
    rho = [(np.sum(w[:SEG - m * paso] * w[m * paso:]) / np.sum(w ** 2)) ** 2 if m * paso < SEG else 0.0
           for m in range(1, K)]                   # solo se solapan segmentos vecinos
    nu = 2 * K / (1 + 2 * sum((1 - m / K) * r for m, r in enumerate(rho, start=1)))
    return f[1:], S[1:], K, nu


def ar1_fondo(x, f):
    x = np.asarray(x, float)
    rho = np.corrcoef(x[:-1], x[1:])[0, 1]
    var = np.var(x)
    S = 2 * var * (1 - rho ** 2) / (1 + rho ** 2 - 2 * rho * np.cos(2 * np.pi * f))
    return rho, S


def fraccion_bandas(f, S):
    T = 1 / f
    tot = np.sum(S)
    return {b: float(np.sum(S[(T >= lo) & (T < hi)]) / tot) for b, lo, hi in BANDAS}


def picos(f, S, N, n_max, fondo=None, nu=2, M=None):
    loc = signal.argrelmax(np.concatenate([[0], S, [0]]))[0] - 1
    loc = loc[np.argsort(S[loc])[::-1]][:n_max]
    filas = []
    for i in loc:
        k = int(round(f[i] * N))
        fila = {"f_ciclos_mes": f[i], "k_ciclos_en_registro": k, "periodo_meses": 1 / f[i],
                "periodo_min_meses": N / (k + 0.5), "periodo_max_meses": N / (k - 0.5) if k > 0.5 else np.inf,
                "fraccion_varianza": S[i] / S.sum()}
        if fondo is not None:
            fila["cociente_fondo_AR1"] = S[i] / fondo[i]
            fila["supera_95_puntual"] = bool(S[i] > fondo[i] * stats.chi2.ppf(1 - ALFA, nu) / nu)
            fila["supera_95_global"] = bool(S[i] > fondo[i] * stats.chi2.ppf(1 - ALFA / M, nu) / nu)
        filas.append(fila)
    return filas


# ------------------------------------------------------------------ cálculo principal
PERIODOS = {"completo": lambda k: SER[k],
            "comun": lambda k: SER[k].loc[PC_INI:min(PC_FIN, SER[k].index.max())],
            "1998-2020": lambda k: SER[k].loc["1998":"2020"]}
TRANS = ["centrada", "anomalia", "centrada_sin_tendencia", "anomalia_sin_tendencia"]
meta, pk, bd, esp, ESPEC = [], [], [], [], {}
for k in VARS:
    u = SERIES[k][2]
    for per, fn in PERIODOS.items():
        if per == "1998-2020" and not k.startswith("T"):
            continue
        s = fn(k)
        N = len(s)
        for tr in TRANS:
            x = transformar(k, s, tr)
            f, S = periodograma(x)
            fh, Sh = periodograma(x, "hann")
            fw, Sw, K, nu = welch(x)
            parseval = np.sum(S) * (1 / N) / np.var(x.values)
            rho, fondo = ar1_fondo(x.values, f)
            ESPEC[(k, per, tr)] = dict(f=f, S=S, fh=fh, Sh=Sh, fw=fw, Sw=Sw, nu=nu, K=K, rho=rho, fondo=fondo, N=N,
                                       x=x)
            meta.append({"variable": k, "periodo": per, "inicio": f"{s.index.min():%Y-%m}", "fin": f"{s.index.max():%Y-%m}",
                         "N_meses": N, "dt_meses": 1, "df_ciclos_mes": 1 / N, "periodo_max_resoluble_meses": N,
                         "nyquist_ciclos_mes": 0.5, "transformacion": tr,
                         "estimador": "periodograma FFT (rect.), Hann y Welch", "normalizacion":
                         "densidad unilateral, Σ S·Δf = varianza", "unidades_densidad": f"({u})²/(ciclo/mes)",
                         "varianza": float(np.var(x.values)), "parseval_cociente": parseval,
                         "anual_en_bin_exacto": N % 12 == 0, "welch_segmentos": K, "welch_nu": nu,
                         "welch_df": 1 / SEG, "rho1_AR1": rho,
                         "vacios": "ninguno" if k != "Q" else "80 meses rellenados (R_flag); ver Lomb–Scargle"})
            M = len(f)
            prueba = tr.startswith("anomalia")
            for p in picos(f, S, N, 6 if prueba else 4, fondo if prueba else None, 2, M):
                pk.append({"variable": k, "periodo": per, "transformacion": tr, "estimador": "periodograma",
                           "N_meses": N, "df": 1 / N} | p)
            for b, v in fraccion_bandas(f, S).items():
                bd.append({"variable": k, "periodo": per, "transformacion": tr, "banda": b, "fraccion": v})
            for i in range(len(f)):
                esp.append({"variable": k, "periodo": per, "transformacion": tr, "estimador": "periodograma",
                            "f": f[i], "S": S[i], "fondo_AR1": fondo[i]})
            for i in range(len(fw)):
                esp.append({"variable": k, "periodo": per, "transformacion": tr, "estimador": "welch",
                            "f": fw[i], "S": Sw[i]})
META = pd.DataFrame(meta)
PK = pd.DataFrame(pk)
BD = pd.DataFrame(bd)
META.to_csv(OUT / "metadatos_espectros_4_1.csv", index=False, float_format="%.5g")
PK.to_csv(OUT / "picos_4_2.csv", index=False, float_format="%.4g")
BD.to_csv(OUT / "bandas_varianza_4_2.csv", index=False, float_format="%.4f")
pd.DataFrame(esp).to_csv(OUT / "espectros_4.csv", index=False, float_format="%.6g")
print("=== 4.1 Parseval (Σ S·Δf / varianza): ", f"{META['parseval_cociente'].min():.6f}–{META['parseval_cociente'].max():.6f}")
print(META[META.transformacion == "anomalia"][["variable", "periodo", "inicio", "fin", "N_meses", "df_ciclos_mes",
                                                "anual_en_bin_exacto", "welch_segmentos", "welch_nu", "rho1_AR1"]]
      .round(4).to_string(index=False))

# Verificación del estimador propio contra scipy.signal.periodogram
x = ESPEC[("P_L", "completo", "anomalia")]["x"].values
f_s, S_s = signal.periodogram(x, fs=1.0, window="boxcar", detrend=False, scaling="density")
assert np.allclose(S_s[1:], ESPEC[("P_L", "completo", "anomalia")]["S"]), "estimador propio ≠ scipy"
print("Estimador propio = scipy.signal.periodogram ✓")

# ------------------------------------------------------------------ fracciones de varianza por banda
piv = BD[BD.transformacion.isin(["centrada", "anomalia"])].pivot_table(
    index=["variable", "periodo", "transformacion"], columns="banda", values="fraccion")[[b for b, *_ in BANDAS]]
print("\n=== 4.2 Fracción de varianza por banda (periodograma) ===")
print((piv * 100).round(1).to_string())

# ------------------------------------------------------------------ picos principales
print("\n=== Picos de las anomalías (registro completo), con fondo AR(1) ===")
v = PK[(PK.periodo == "completo") & (PK.transformacion == "anomalia")]
print(v[["variable", "periodo_meses", "periodo_min_meses", "periodo_max_meses", "k_ciclos_en_registro",
         "fraccion_varianza", "cociente_fondo_AR1", "supera_95_puntual", "supera_95_global"]].round(3).to_string(index=False))

# ------------------------------------------------------------------ estabilidad (4.2.4)
est = []


def registrar(k, variante, x, ventana="rect", detalle=""):
    f, S = periodograma(x, ventana)
    fb = fraccion_bandas(f, S)
    T = 1 / f
    inter = (T >= 24) & (T <= 96)
    i = np.argmax(np.where(inter, S, -np.inf))
    est.append({"variable": k, "variante": variante, "N": len(x),
                "detalle": detalle or f"{x.index[0]:%Y-%m}–{x.index[-1]:%Y-%m}", "pico_interanual_meses": T[i],
                "rango_pico_meses": f"{len(x) / (round(f[i] * len(x)) + 0.5):.0f}–{len(x) / (round(f[i] * len(x)) - 0.5):.0f}",
                "frac_anual": fb["anual 11–13 m"], "frac_semianual": fb["semianual 5.5–7 m"],
                "frac_interanual_2_7a": fb["interanual 2–7 años"], "frac_mas_7a": fb["> 7 años"],
                "frac_alta_2_5m": fb["intraestacional 2–5 m"]})


for k in VARS:
    s = SER[k]
    a = transformar(k, s, "anomalia")
    registrar(k, "anomalía, rectangular (base)", a)
    registrar(k, "anomalía, Hann", a, "hann")
    registrar(k, "anomalía sin tendencia", sin_tendencia(a))
    lo, hi = np.percentile(a, [1, 99])
    registrar(k, "anomalía recortada P1–P99 (extremos)", a.clip(lo, hi) - a.clip(lo, hi).mean())
    n2 = len(a) // 2
    n2 -= n2 % 12
    registrar(k, "primera mitad del registro", a.iloc[:n2] - a.iloc[:n2].mean())
    registrar(k, "segunda mitad del registro", a.iloc[-n2:] - a.iloc[-n2:].mean())
    c = PERIODOS["comun"](k)
    registrar(k, "periodo común", transformar(k, c, "anomalia"))
    if k == "Q":                               # caudal muy asimétrico: crecidas de 1998 y 2017 dominan la varianza
        lq = np.log(s)
        mu = lq.loc[REF[0]:REF[1]].groupby(lq.loc[REF[0]:REF[1]].index.month).mean()
        al = lq - lq.index.month.map(mu).values
        registrar(k, "anomalía de ln Q (extremos)", al - al.mean())
EST = pd.DataFrame(est)
EST.to_csv(OUT / "estabilidad_4_2.csv", index=False, float_format="%.4f")
print("\n=== Estabilidad (anomalías) ===")
print(EST.round(3).to_string(index=False))

# ------------------------------------------------------------------ vacíos de caudal: Lomb–Scargle
qa_f = transformar("Q", SER["Q"].loc[PC_INI:PC_FIN], "anomalia")
qo = Q_OBS.loc[PC_INI:PC_FIN]
ref = SER["Q"].loc[REF[0]:REF[1]]
qo_a = qo - qo.index.month.map(ref.groupby(ref.index.month).mean()).values
qo_a = qo_a - qo_a.mean()
tq = (qo_a.index.year - PC_INI.year) * 12 + qo_a.index.month - 1
fq, Sq = periodograma(qa_f)
ls = signal.lombscargle(tq.values.astype(float), qo_a.values, 2 * np.pi * fq, precenter=True)
tr_ = np.arange(len(qa_f), dtype=float)        # control: Lomb–Scargle de la serie rellenada en malla regular
ls_reg = signal.lombscargle(tr_, qa_f.values, 2 * np.pi * fq, precenter=True)
ls_n, fq_n = ls / ls.sum(), Sq / Sq.sum()
fb_ls, fb_fill = fraccion_bandas(fq, ls), fraccion_bandas(fq, Sq)
LS = pd.DataFrame({"banda": list(fb_fill), "relleno_FFT": list(fb_fill.values()),
                   "relleno_LombScargle_malla_regular": list(fraccion_bandas(fq, ls_reg).values()),
                   "observado_LombScargle": list(fb_ls.values())})
print(f"Control: FFT y Lomb–Scargle en malla regular coinciden: "
      f"{np.allclose(ls_reg / ls_reg.sum(), fq_n, atol=1e-6)}")
LS.to_csv(OUT / "vacios_caudal_lombscargle_4_1.csv", index=False, float_format="%.4f")
print(f"\n=== Vacíos del caudal ({len(qo_a)} meses observados de {len(qa_f)}) ===")
print((LS.set_index("banda") * 100).round(1).to_string())
print(f"r entre espectros normalizados (log): {np.corrcoef(np.log(fq_n), np.log(ls_n + 1e-12))[0, 1]:.2f}")

# ------------------------------------------------------------------ figuras
aplicar_estilo()
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
GRIS = TINTA2


def eje_doble(ax, fmin=1 / 600, fmax=0.5, rotulos=True):
    """Datos en frecuencia f (ciclos/mes), escala log. Eje inferior: periodo en años (= 1/(12 f)), creciente hacia la
    derecha; eje superior: la frecuencia en ciclos/mes que pide la guía (Nyquist 0.5)."""
    ax.set_xscale("log")
    ax.set_xlim(fmax, fmin)                                   # periodo crece hacia la derecha
    ax.set_xticks([1 / (12 * a) for a in PERIODOS_A], [f"{a:g}" for a in PERIODOS_A])
    ax.set_xticks([], minor=True)
    sec = ax.secondary_xaxis("top")
    sec.set_xticks([f for f in FRECUENCIAS if fmin <= f <= fmax],
                   [f"1/{round(1 / f)}" if f < 0.5 else "0.5" for f in FRECUENCIAS if fmin <= f <= fmax])
    sec.set_xticks([], minor=True)
    sec.tick_params(labelsize=7, colors=GRIS)
    if rotulos:
        sec.set_xlabel("frecuencia (ciclos/mes)", fontsize=7.5, color=GRIS, labelpad=2)
    for m in (12, 6):
        ax.axvline(1 / m, color=GRIS, lw=0.6, ls=":")
    ax.axvspan(1 / 84, 1 / 24, color=GRIS_CLARO, alpha=0.3, lw=0)
    ax.grid(alpha=0.2, lw=0.5, which="both")


# Figura 1: originales centradas vs anomalías (registro completo), periodograma + Welch + fondo AR(1).
# Misma escala vertical en las dos columnas de cada variable: se ve cuánto cae la potencia al quitar la climatología.
fig, axs = plt.subplots(len(VARS), 2, figsize=(14, 3.3 * len(VARS)))
for i, k in enumerate(VARS):
    col, u = SERIES[k][3], SERIES[k][2]
    ee = [ESPEC[(k, "completo", tr)] for tr in ("centrada", "anomalia")]
    ylim = (min(np.percentile(e["S"], 1) for e in ee) / 3, max(e["S"].max() for e in ee) * 3)
    for j, (tr, e) in enumerate(zip(("centrada", "anomalia"), ee)):
        ax = axs[i, j]
        ax.plot(e["f"], e["S"], color=col, lw=0.8, alpha=0.75, label="periodograma (rect.)")
        ax.plot(e["fw"], e["Sw"], color=TINTA, lw=1.4, label=f"Welch (Hann, {SEG} m, K={e['K']}, ν≈{e['nu']:.0f})")
        if tr == "anomalia":
            M = len(e["f"])
            ax.plot(e["f"], e["fondo"], color=GRIS, lw=1.0, ls="-", label=f"AR(1) ρ={e['rho']:.2f}")
            ax.plot(e["f"], e["fondo"] * stats.chi2.ppf(1 - ALFA, 2) / 2, color=TINTA2, lw=0.9, ls="--",
                    label="95 % puntual")
            ax.plot(e["f"], e["fondo"] * stats.chi2.ppf(1 - ALFA / M, 2) / 2, color=TINTA2, lw=0.9, ls=":",
                    label=f"95 % global (M={M})")
        ax.set_yscale("log")
        eje_doble(ax, rotulos=(i == 0))
        ax.set_ylim(*ylim)
        ax.set_ylabel(f"S(f) [({u})²·mes]")
        ax.set_title(f"{SERIES[k][1]} · {'centrada (con ciclo anual)' if tr == 'centrada' else 'anomalía (sin climatología)'}"
                     f" · {e['x'].index[0]:%Y-%m}–{e['x'].index[-1]:%Y-%m}, N={e['N']}, Δf=1/{e['N']}", fontsize=8.5,
                     loc="left", pad=28 if i == 0 else 16)
        if i == 0:
            ax.legend(fontsize=6.8, frameon=False, loc="lower right", ncol=1 if j == 0 else 2)
        if i == len(VARS) - 1:
            ax.set_xlabel("periodo (años)")
fig.suptitle("Punto 4 · Espectros de densidad unilaterales (Σ S·Δf = varianza) del registro completo de cada variable\n"
             "Eje inferior: periodo (años); superior: frecuencia (ciclos/mes) · punteadas verticales: 1 año y 6 meses · "
             "sombreado: banda 2–7 años · misma escala vertical en las dos columnas de cada variable", fontsize=9.5)
fig.tight_layout(rect=(0, 0, 1, 0.975))
fig.savefig(OUT / "fig_espectros_4_2.png", dpi=150)
plt.close(fig)

# Figura 1b: los mismos periodogramas en % de la varianza, eje de periodo lineal (años). Lectura directa de "qué
# periodos dominan" (idea del visor de la compañera). % de cada frecuencia de Fourier = S_k·Δf / varianza. Los picos se
# rotulan con el rango que permite la resolución (±½ bin) y solo si caben al menos 3 ciclos en el registro.
T_MAX = 20                                                       # años mostrados; periodos más largos ≈ tendencia


def picos_pct(T, pct, k, ymax, n_max=4, t_min=0.4):
    """Máximos locales con t_min ≤ T ≤ T_MAX y k ≥ 3 ciclos, del mayor al menor. Se descarta un pico si su rótulo
    caería encima de otro ya elegido (cerca en periodo y en altura)."""
    i = np.arange(1, len(pct) - 1)
    loc = i[(pct[i] > pct[i - 1]) & (pct[i] >= pct[i + 1]) & (T[i] >= t_min) & (T[i] <= T_MAX) & (k[i] >= 3)]
    def caja(j, w=2.5, h=0.14 * ymax):                      # extensión aproximada del rótulo (años, %)
        if T[j] < 2:                                         # a la derecha del pico, centrado en vertical
            return T[j], T[j] + w, pct[j] - h / 2, pct[j] + h / 2
        return T[j] - w / 2, T[j] + w / 2, pct[j], pct[j] + h  # encima del pico

    def chocan(a, b):
        return a[0] < b[1] and b[0] < a[1] and a[2] < b[3] and b[2] < a[3]

    elegidos = []
    for j in loc[np.argsort(pct[loc])[::-1]]:
        if not any(chocan(caja(j), caja(m)) or abs(T[j] - T[m]) < 0.15 * T[m] for m in elegidos):
            elegidos.append(j)
        if len(elegidos) == n_max:
            break
    return elegidos


def fmt_a(t):
    return f"{t:.1f}" if t < 10 else f"{t:.0f}"


fig, axs = plt.subplots(len(VARS), 2, figsize=(14, 3.0 * len(VARS)), sharex=True)
for i, k_ in enumerate(VARS):
    col = SERIES[k_][3]
    for j, tr in enumerate(("centrada", "anomalia")):
        ax = axs[i, j]
        e = ESPEC[(k_, "completo", tr)]
        var, N = np.var(e["x"].values), e["N"]
        kk = np.rint(e["f"] * N)                                    # índice de Fourier = ciclos en el registro
        T = 1 / (12 * e["f"])
        pct = 100 * e["S"] * (1 / N) / var
        M = len(e["f"])
        umb = 100 * e["fondo"] * stats.chi2.ppf(1 - ALFA, 2) / 2 * (1 / N) / var
        umb_g = 100 * e["fondo"] * stats.chi2.ppf(1 - ALFA / M, 2) / 2 * (1 / N) / var
        ok = T <= T_MAX
        ax.vlines(T[ok], 0, pct[ok], color=col, lw=1.3)
        ax.plot(T[ok], umb[ok], color=TINTA2, lw=0.9, ls="--")
        ymax = pct[ok].max() * 1.45                                  # margen para los rótulos
        for p in picos_pct(T, pct, kk, ymax):
            lo_, hi_ = N / (12 * (kk[p] + 0.5)), N / (12 * (kk[p] - 0.5))
            marca = " **" if pct[p] > umb_g[p] else (" *" if pct[p] > umb[p] else "")
            ax.plot(T[p], pct[p], "o", color=col, ms=4)
            corto = T[p] < 2                                         # cerca del eje: rótulo a la derecha del pico
            ax.annotate(f"≈{fmt_a(T[p])} a ({fmt_a(lo_)}–{fmt_a(hi_)}){marca}\n{pct[p]:.1f} % · {int(kk[p])} ciclos",
                        (T[p], pct[p]), xytext=(5, 0) if corto else (0, 4), textcoords="offset points",
                        ha="left" if corto else "center", va="center" if corto else "bottom", fontsize=7,
                        bbox=dict(fc="white", ec="none", alpha=0.8, pad=0.5))
        fuera = 100 * (e["S"][~ok] * (1 / N)).sum() / var
        ax.text(0.99, 0.95, f"periodos > {T_MAX} años: {fuera:.1f} % (fuera del eje)", transform=ax.transAxes,
                ha="right", va="top", fontsize=7, color=TINTA2)
        ax.set_xlim(0, T_MAX)
        ax.set_ylim(0, ymax)
        ax.set_title(f"{SERIES[k_][1]} · {'centrada (con ciclo anual)' if tr == 'centrada' else 'anomalía (sin climatología)'}"
                     f" · {e['x'].index[0]:%Y-%m}–{e['x'].index[-1]:%Y-%m}, N = {N}", fontsize=8.5, loc="left")
        if j == 0:
            ax.set_ylabel("% de la varianza total")
        ax.grid(alpha=0.25, lw=0.5)
for ax in axs[-1]:
    ax.set_xlabel("periodo (años)")
    ax.set_xticks(range(0, T_MAX + 1, 2))
fig.suptitle("Punto 4 · Periodograma en % de la varianza total por periodo (registro completo; escala vertical propia de "
             "cada panel)\nCada línea es una frecuencia de Fourier k/N. Rótulo: periodo ≈ T con el rango de ±½ bin de "
             "frecuencia, % de la varianza y ciclos en el registro (solo picos con ≥ 3 ciclos).\nDiscontinua: 95 % puntual "
             "de un fondo AR(1). * supera ese umbral; ** supera también el umbral global (búsqueda entre todas las "
             "frecuencias).",
             fontsize=9.5)
fig.tight_layout(rect=(0, 0, 1, 0.955))
fig.savefig(OUT / "fig_espectros_porcentaje_4_2.png", dpi=150)
plt.close(fig)

# Figura 2: fracción de varianza por banda
fig, axs = plt.subplots(1, 2, figsize=(14, 5.2), sharey=True)
colb = colores_ordinales(len(BANDAS))                 # de periodos cortos (oscuro) a largos (claro)
for ax, tr in zip(axs, ("centrada", "anomalia")):
    filas = [(k, per) for k in VARS for per in ("completo", "comun")]
    etiquetas = []
    base = np.zeros(len(filas))
    for i_b, ((b, *_), c) in enumerate(zip(BANDAS, colb)):
        vals = np.array([BD[(BD.variable == k) & (BD.periodo == per) & (BD.transformacion == tr) & (BD.banda == b)]
                         ["fraccion"].iloc[0] for k, per in filas])
        ax.barh(range(len(filas)), vals * 100, left=base.copy() * 100, color=c, edgecolor="white", lw=1, label=b)
        for y, (v_, b0) in enumerate(zip(vals, base)):
            if v_ > 0.08:
                ax.text((b0 + v_ / 2) * 100, y, f"{v_ * 100:.0f}", ha="center", va="center", fontsize=7,
                        color="white" if i_b < 3 else "k")
        base += vals
    for k, per in filas:
        m = META[(META.variable == k) & (META.periodo == per) & (META.transformacion == tr)].iloc[0]
        etiquetas.append(f"{k} · {m['inicio'][:4]}–{m['fin'][:4]} (N={m['N_meses']})")
    ax.set_yticks(range(len(filas)), etiquetas, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.set_xlabel("% de la varianza")
    ax.set_title("Serie centrada (con ciclo anual)" if tr == "centrada" else "Anomalías (sin climatología)",
                 fontsize=9.5, loc="left")
axs[1].legend(fontsize=7.5, frameon=False, loc="lower right", bbox_to_anchor=(1.0, -0.32), ncol=4)
fig.suptitle("Punto 4 · Reparto de la varianza por bandas de periodo (periodograma)", fontsize=10.5)
fig.tight_layout(rect=(0, 0.04, 1, 0.96))
fig.savefig(OUT / "fig_bandas_varianza_4_2.png", dpi=150)
plt.close(fig)

# Figura 3: comparación entre fuentes y entre lluvia y caudal (periodo común, Welch normalizado por la varianza)
fig, axs = plt.subplots(2, 2, figsize=(13.5, 8.5))
pares = [(("P_L", "comun"), ("P_I", "comun"), "Precipitación: CAMELS-PE vs IMERG"),
         (("T_PISCOt", "1998-2020"), ("T_ERA5L", "1998-2020"), "Temperatura: PISCOt vs ERA5-Land (1998–2020)"),
         (("P_L", "comun"), ("Q", "comun"), "Lluvia local vs caudal")]
for ax, (a_, b_, tit) in zip(axs.flat[:3], pares):
    for kk, per in (a_, b_):
        for tr, ls_ in (("centrada", "-"), ("anomalia", "--")):
            e = ESPEC[(kk, per, tr)]
            var = np.var(e["x"].values)
            ax.plot(e["fw"], e["Sw"] / var, ls=ls_, color=SERIES[kk][3], lw=1.5,
                    label=f"{kk} {'centrada' if tr == 'centrada' else 'anomalía'}")
    ax.set_yscale("log")
    eje_doble(ax, fmin=1 / SEG)
    ax.set_ylabel("S(f)/varianza [mes]")
    ax.set_title(tit, fontsize=9.5, loc="left", pad=28)
    ax.legend(fontsize=7, frameon=False, loc="lower left")
ax = axs.flat[3]
for tr, ls_ in (("centrada", "-"), ("anomalia", "--")):
    ep, eq = ESPEC[("P_L", "comun", tr)], ESPEC[("Q", "comun", tr)]
    r = (eq["Sw"] / np.var(eq["x"].values)) / (ep["Sw"] / np.var(ep["x"].values))
    ax.plot(ep["fw"], r, ls=ls_, color=TINTA, lw=1.5, label=f"{'centrada' if tr == 'centrada' else 'anomalía'}")
ax.axhline(1, color=GRIS, lw=0.7)
ax.set_yscale("log")
eje_doble(ax, fmin=1 / SEG)
ax.set_ylabel("cociente de espectros normalizados Q / P_L")
ax.set_title("¿Filtra la cuenca las variaciones rápidas? (> 1: el caudal concentra más\nvarianza relativa en esa escala que la lluvia; no contiene fases)",
             fontsize=9, loc="left", pad=28)
ax.legend(fontsize=7.5, frameon=False)
for a in axs[1]:
    a.set_xlabel("periodo (años)")
fig.suptitle(f"Punto 4 · Formas espectrales normalizadas por la varianza (Welch, Hann, {SEG} meses, 50 % traslape); "
             "las alturas no son variabilidad absoluta", fontsize=10)
fig.tight_layout(rect=(0, 0, 1, 0.97))
fig.savefig(OUT / "fig_comparacion_espectros_4_2.png", dpi=150)
plt.close(fig)

# Figura 4: estabilidad
ORDEN_V = list(dict.fromkeys(EST["variante"]))
fig, axs = plt.subplots(1, len(VARS), figsize=(19, 5.2), sharey=True)
for ax, k in zip(axs, VARS):
    d = EST[EST.variable == k].set_index("variante").reindex(ORDEN_V).reset_index()
    y = np.arange(len(d))
    ax.barh(y - 0.2, d["frac_anual"] * 100, height=0.38, color=colb[3], label="anual 11–13 m")
    ax.barh(y + 0.2, d["frac_interanual_2_7a"] * 100, height=0.38, color=colb[5], label="interanual 2–7 a")
    for yi, r in zip(y, d.itertuples()):
        if np.isnan(r.frac_anual):
            ax.text(0.3, yi, "no aplica", va="center", fontsize=7, color=GRIS)
            continue
        ax.text(max(r.frac_anual, r.frac_interanual_2_7a) * 100 + 1, yi, f"pico {r.pico_interanual_meses / 12:.1f} a",
                va="center", fontsize=7, color=GRIS)
    ax.set_yticks(y, d["variante"], fontsize=7.5)
    ax.invert_yaxis()
    ax.set_title(f"{SERIES[k][1]} (anomalías, N base={int(d.N.iloc[0])})", fontsize=9, loc="left")
    ax.set_xlabel("% de la varianza")
    ax.grid(alpha=0.2, axis="x")
axs[0].legend(fontsize=7.5, frameon=False, loc="lower right")
fig.suptitle("Punto 4.2 · Estabilidad del reparto de varianza y del pico interanual (2–8 años) frente a ventana, tendencia, "
             "extremos y subperiodo", fontsize=10)
fig.tight_layout(rect=(0, 0, 1, 0.95))
fig.savefig(OUT / "fig_estabilidad_4_2.png", dpi=150)
plt.close(fig)

# Figura 5: vacíos de caudal
fig, ax = plt.subplots(figsize=(9, 4.2))
ax.plot(fq, fq_n, color=SERIES["Q"][3], lw=1.1, label=f"FFT de la serie rellenada (N={len(qa_f)})")
ax.plot(fq, ls_n, color=TINTA, lw=1.0, ls="--", label=f"Lomb–Scargle solo meses observados completos (n={len(qo_a)})")
ax.set_yscale("log")
eje_doble(ax, fmin=1 / len(qa_f))
ax.set_ylabel("potencia / potencia total")
ax.set_xlabel("periodo (años)")
ax.legend(fontsize=7.5, frameon=False, loc="lower right")
ax.set_title(f"Punto 4.1 · Efecto de los vacíos del caudal sobre el espectro (anomalías, {PC_INI:%Y}–{PC_FIN:%Y})",
             fontsize=9.5, loc="left", pad=28)
fig.tight_layout()
fig.savefig(OUT / "fig_vacios_caudal_4_1.png", dpi=150)
plt.close(fig)

print("\nSalidas en", OUT)
for f_ in sorted(OUT.iterdir()):
    print("  ", f_.name)
