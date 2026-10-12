"""
Punto 4 (complemento): análisis de onditas (wavelets) para ver CUÁNDO aparece la variabilidad fuera del ciclo anual
y si la banda de 2–7 años de la cuenca coincide en el tiempo con la del ENSO.

Método (Torrence & Compo 1998, BAMS 79:61–78; implementación propia con FFT):
  · Series: anomalías mensuales (climatología 1998–2020, como en los puntos 3 y 4) de P_L (1981–2025), Q̄ (1998–2025)
    y T ERA5-Land (1981–2025), estandarizadas (media 0, varianza 1). ENSO: anomalías Niño 1+2 y Niño 3.4 de ERSST v5
    (CPC/NOAA, base 1991–2020), mismo tratamiento.
  · Ondita de Morlet, ω0 = 6; Δt = 1 mes; escalas s_j = 2 meses · 2^(j/12) hasta ~20 años; periodo de Fourier
    = 1.033 s. Se rellena con ceros hasta la siguiente potencia de 2 (procedimiento estándar de la transformada; los
    datos no se rellenan): las zonas afectadas por los bordes quedan fuera del cono de influencia (COI), sombreado.
  · Corrección de sesgo (Liu, Liang y Weisberg 2007, J. Atmos. Oceanic Technol. 24:2093–2102): la potencia |W|² sin
    corregir infla los periodos largos; se divide por la escala. Expresada como 100·|W|²/(C_δ·s)/σ² es el % de la
    varianza por octava de periodo (su integral en escala recupera la varianza, TC98 ec. 14), la misma unidad que el
    periodograma de Fourier por octava, 100·S_k·f_k·ln2/σ², que se superpone al espectro global para compararlos.
  · Significancia: fondo de ruido rojo AR(1) con la autocorrelación de rezago 1 de cada serie; 95 % puntual con
    χ²₂/2 (contorno negro). Espectro global (promedio en el tiempo) con sus grados de libertad (TC98 ec. 23).
  · Varianza promediada en la banda 2–7 años (TC98 ec. 24) con su umbral de 95 % (ec. 25–28), y reconstrucción de la
    componente de 2–7 años (ec. 29) para correlacionar cuenca y ENSO en esa banda.
  · Con series filtradas a 2–7 años quedan muy pocos grados de libertad (≈ 1 por cada 2–4 años de registro): las
    correlaciones entre componentes son descriptivas y se reportan con ese aviso, sin valor p.

Salidas (Tarea 1/El_Tambo/resultados/espectros/):
  fig_onditas_4_3.png, fig_onditas_banda_4_3.png, onditas_resumen_4_3.csv, onditas_banda_4_3.csv
"""
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from scipy import stats

from estilo_figuras import VAR, TINTA, TINTA2, GRIS, GRIS_CLARO, CMAP_P, aplicar_estilo

aplicar_estilo()
ROOT = Path(__file__).resolve().parents[2] / "Tarea 1" / "El_Tambo"
SM = ROOT / "resultados" / "series_mensuales"
OUT = ROOT / "resultados" / "espectros"
REF = ("1998-01-01", "2020-12-01")
W0, DJ, S0 = 6.0, 1 / 12, 2.0                                    # Morlet, resolución en escala, escala mínima (meses)
BANDA = (24, 84)                                                 # 2–7 años, en meses
CDELTA, PSI0, DJ0, GAMMA = 0.776, np.pi ** -0.25, 0.60, 2.32     # constantes de Morlet ω0 = 6 (TC98, tabla 2)
NINOS = [1983, 1998, 2016, 2017, 2023]                           # años de El Niño fuerte o costero (marcas de lectura)

cam = pd.read_csv(SM / "series_mensuales_camels_pe.csv", index_col="fecha", parse_dates=True)
t5l = pd.read_csv(SM / "era5land_t2m_cuenca_mensual.csv", index_col="fecha", parse_dates=True)["T_era5land_C"]
nino = pd.read_csv(ROOT / "datos" / "indices_clima" / "ersst5.nino.mth.91-20.ascii", sep=r"\s+")
nino.index = pd.to_datetime(dict(year=nino.YR, month=nino.MON, day=1))


def anomalia(s):
    s = s.dropna().sort_index()
    ref = s.loc[REF[0]:REF[1]]
    return s - s.index.month.map(ref.groupby(ref.index.month).mean()).values


SERIES = {  # clave: (anomalías, nombre, color)
    "P_L": (anomalia(cam["P_L_mm_mes"]).loc["1981":"2025"], "P$_L$ (CAMELS-PE)", VAR["P_L"]),
    "Q": (anomalia(cam["Q_m3_s"]).loc["1998":"2025"], "$\\bar{Q}$ (El Tambo)", VAR["Q"]),
    "T_ERA5L": (anomalia(t5l).loc["1981":"2025"], "T (ERA5-Land)", VAR["T_ERA5L"]),
    "N12": (nino.iloc[:, 3].loc["1981":"2025"], "Niño 1+2 (SST)", TINTA),
    "N34": (nino.iloc[:, 9].loc["1981":"2025"], "Niño 3.4 (SST)", GRIS),
}
for k, (s, *_r) in SERIES.items():
    assert len(pd.date_range(s.index[0], s.index[-1], freq="MS")) == len(s), f"{k} tiene vacíos internos"


def onditas(x):
    """Transformada continua de Morlet (TC98). x estandarizada. Devuelve W (escalas × tiempo), escalas, periodos, COI."""
    n = len(x)
    npad = int(2 ** np.ceil(np.log2(n)))
    xp = np.concatenate([x, np.zeros(npad - n)])
    k = 2 * np.pi * np.fft.fftfreq(npad)                       # Δt = 1 mes
    J = int(np.log2(240 / S0) / DJ)                              # hasta 240 meses (20 años) de escala
    esc = S0 * 2 ** (DJ * np.arange(J + 1))
    xh = np.fft.fft(xp)
    W = np.empty((len(esc), n), dtype=complex)
    for j, s in enumerate(esc):
        psi = np.sqrt(2 * np.pi * s) * PSI0 * np.exp(-(s * k - W0) ** 2 / 2) * (k > 0)
        W[j] = np.fft.ifft(xh * psi)[:n]
    fac = 4 * np.pi / (W0 + np.sqrt(2 + W0 ** 2))
    per = fac * esc
    coi = fac / np.sqrt(2) * np.concatenate([[1e-5], np.arange(1, (n + 1) // 2), np.arange((n // 2) - 1, 0, -1), [1e-5]])[:n]
    return W, esc, per, coi


def fondo_rojo(alfa, per):
    f = 1 / per
    return (1 - alfa ** 2) / (1 + alfa ** 2 - 2 * alfa * np.cos(2 * np.pi * f))


RES, BAND, COMP = {}, {}, {}
filas = []
for k, (s, nom, col) in SERIES.items():
    x = ((s - s.mean()) / s.std()).values
    n = len(x)
    alfa = np.corrcoef(x[:-1], x[1:])[0, 1]
    W, esc, per, coi = onditas(x)
    pot = np.abs(W) ** 2
    teor = fondo_rojo(alfa, per)
    sig = pot / (teor[:, None] * stats.chi2.ppf(0.95, 2) / 2)    # > 1: supera el 95 % puntual
    glob = pot.mean(1)
    dof = 2 * np.sqrt(1 + (n / (GAMMA * esc)) ** 2)
    glob_sig = teor * stats.chi2.ppf(0.95, dof) / dof
    b = (per >= BANDA[0]) & (per <= BANDA[1])
    savg = DJ / CDELTA * (pot[b] / esc[b, None]).sum(0)          # varianza en la banda (unidades de varianza)
    Sb = (1 / esc[b]).sum(); Smid = np.exp(np.log(esc[b]).mean()); na = b.sum()
    dof_b = (na * Sb / Smid) * np.sqrt(1 + (na * DJ / DJ0) ** 2)
    sig_b = DJ / CDELTA * (teor[b] / esc[b]).sum() * stats.chi2.ppf(0.95, dof_b) / dof_b
    comp = DJ / (CDELTA * PSI0) * (W[b].real / np.sqrt(esc[b, None])).sum(0)   # componente 2–7 años (σ = 1)
    dentro = per[:, None] <= coi[None, :]
    var = np.var(x)
    rect = 100 * pot / (CDELTA * esc[:, None]) / var             # Liu et al. 2007: % de la varianza por octava
    glob_oct = rect.mean(1)
    glob_sig_oct = 100 * glob_sig / (CDELTA * esc) / var
    X = np.fft.rfft(x - x.mean())
    fk = np.arange(len(X)) / n
    Sk = 2 * np.abs(X) ** 2 / n                                  # densidad unilateral, Δt = 1 mes
    four_oct = 100 * Sk[1:] * fk[1:] * np.log(2) / var
    RES[k] = dict(t=s.index, pot=pot, sig=sig, per=per, coi=coi, glob=glob, glob_sig=glob_sig, alfa=alfa, nom=nom, col=col,
                  rect=rect, glob_oct=glob_oct, glob_sig_oct=glob_sig_oct, four_per=1 / fk[1:] / 12, four_oct=four_oct,
                  suma_oct=float(glob_oct.sum() * DJ))
    BAND[k] = pd.Series(savg, index=s.index)
    COMP[k] = pd.Series(comp, index=s.index)
    frac_sig_banda = float(((sig[b] > 1) & dentro[b]).any(0).mean())
    filas.append({"serie": k, "inicio": f"{s.index[0]:%Y-%m}", "fin": f"{s.index[-1]:%Y-%m}", "N": n, "alfa_AR1": alfa,
                  "var_componente_2_7a": float(np.var(comp)),
                  "frac_tiempo_banda_2_7a_significativa": frac_sig_banda,
                  "umbral_95_varianza_banda": float(sig_b),
                  "frac_tiempo_varianza_banda_sobre_umbral": float((savg > sig_b).mean()),
                  "periodo_max_espectro_global_2_7a_anios": float(per[b][np.argmax(glob[b])] / 12),
                  "global_supera_95_en_banda": bool((glob[b] > glob_sig[b]).any()),
                  "pico_global_corregido_anios": float(per[np.argmax(glob_oct)] / 12),
                  "pico_global_corregido_2_7a_anios": float(per[b][np.argmax(glob_oct[b])] / 12),
                  "pct_varianza_reconstruida": float(glob_oct.sum() * DJ)})
    RES[k]["sig_b"] = sig_b

# correlación entre componentes de 2–7 años (descriptiva)
for f in filas:
    k = f["serie"]
    for ref in ("N12", "N34"):
        a, b_ = COMP[k].align(COMP[ref], join="inner")
        f[f"r_comp_2_7a_con_{ref}"] = float(np.corrcoef(a, b_)[0, 1]) if k != ref else 1.0
        cc = [(lag, np.corrcoef(a.values[lag:], b_.values[:len(b_) - lag])[0, 1]) for lag in range(0, 13)]
        lag, r = max(cc, key=lambda z: z[1])
        f[f"r_max_{ref}_rezago_0_12m"] = float(r)
        f[f"rezago_meses_{ref}_antecede"] = int(lag)
        f["n_meses_comunes_" + ref] = len(a)
RS = pd.DataFrame(filas)
RS.round(3).to_csv(OUT / "onditas_resumen_4_3.csv", index=False)
pd.DataFrame({f"varianza_banda_{k}": v for k, v in BAND.items()} | {f"componente_2_7a_{k}": v for k, v in COMP.items()}
             ).round(4).to_csv(OUT / "onditas_banda_4_3.csv", index_label="fecha")
print(RS.round(2).T.to_string())

# ------------------------------------------------------------------ figura 1: potencia de onditas
ORDEN = ["P_L", "Q", "T_ERA5L", "N12"]
NIV = [0, 5, 10, 20, 40, 80, 160, 320]
NORM = BoundaryNorm(NIV, CMAP_P.N, extend="max")
fig = plt.figure(figsize=(15, 3.0 * len(ORDEN)))
gs = fig.add_gridspec(len(ORDEN), 2, width_ratios=[5, 1], hspace=0.45, wspace=0.04, bottom=0.12, top=0.89)
xlim = (pd.Timestamp("1981-01-01"), pd.Timestamp("2026-01-01"))
for i, k in enumerate(ORDEN):
    d = RES[k]
    ax = fig.add_subplot(gs[i, 0])
    T = d["t"]
    P = d["per"] / 12
    m = ax.contourf(T, P, d["rect"], levels=NIV, cmap=CMAP_P, norm=NORM, extend="max")
    ax.contour(T, P, d["sig"], levels=[1], colors=TINTA, linewidths=0.9)
    ax.fill_between(T, d["coi"] / 12, P.max(), color="white", alpha=0.55, hatch="xx", edgecolor=GRIS_CLARO, lw=0)
    ax.axhspan(BANDA[0] / 12, BANDA[1] / 12, color="none", ec=TINTA2, ls="--", lw=0.8)
    for y in NINOS:
        ax.axvline(pd.Timestamp(y, 3, 1), color=TINTA2, lw=0.6, ls=":")
    ax.set_yscale("log", base=2)
    ax.set_yticks([0.25, 0.5, 1, 2, 4, 8, 16], ["0.25", "0.5", "1", "2", "4", "8", "16"])
    ax.set_ylim(20, P.min())                                     # periodos cortos arriba (como Torrence y Compo 1998)
    ax.set_xlim(*xlim)
    ax.set_ylabel("periodo (años)")
    ax.set_title(f"{d['nom']} · anomalías estandarizadas {T[0]:%Y-%m}–{T[-1]:%Y-%m} · AR(1) α = {d['alfa']:.2f}",
                 fontsize=9.5, loc="left")
    g = fig.add_subplot(gs[i, 1], sharey=ax)
    fo = d["four_per"] <= 20
    g.plot(d["four_oct"][fo], d["four_per"][fo], color=GRIS_CLARO, lw=0.8, label="Fourier (periodograma)")
    g.plot(d["glob_oct"], P, color=d["col"], lw=1.8, label="onditas (global)")
    g.plot(d["glob_sig_oct"], P, color=TINTA2, lw=0.9, ls="--", label="95 % AR(1)")
    g.axhspan(BANDA[0] / 12, BANDA[1] / 12, color=GRIS_CLARO, alpha=0.3, lw=0)
    g.set_xlim(0, 1.15 * max(d["glob_oct"].max(), d["glob_sig_oct"][P <= 20].max()))
    g.tick_params(labelleft=False)
    ax.set_ylim(20, P.min())                                     # periodos cortos arriba (como Torrence y Compo 1998)                                     # el panel global comparte el eje: fijarlo al final
    g.set_title("global", fontsize=9, loc="left")
    if i == 0:
        g.legend(fontsize=7, frameon=False, loc="lower right")
    if i == len(ORDEN) - 1:
        g.set_xlabel("% de la varianza por octava")
cax = fig.add_axes([0.2, 0.06, 0.45, 0.012])
cb = fig.colorbar(m, cax=cax, orientation="horizontal", ticks=NIV)
cb.set_label("potencia corregida por sesgo (Liu et al. 2007): % de la varianza por octava de periodo")
fig.suptitle("Punto 4 · Onditas de Morlet: en qué años aparece cada periodo (potencia corregida por sesgo; a la derecha, "
             "espectro global de onditas sobre el periodograma de Fourier, en la misma unidad)\n"
             "Contorno negro: supera el 95 % frente a un ruido rojo AR(1) · rayado: cono de influencia (bordes, no interpretar) · "
             "líneas discontinuas: banda de 2–7 años · punteadas verticales: El Niño 1983, 1998, 2016, 2017, 2023",
             fontsize=10, y=0.94)
fig.savefig(OUT / "fig_onditas_4_3.png", dpi=150, bbox_inches="tight")
plt.close(fig)

# ------------------------------------------------------------------ figura 2: banda 2–7 años frente a ENSO
fig, axs = plt.subplots(2, 1, figsize=(14, 7.5), sharex=True)
ax = axs[0]
for k in ("P_L", "Q", "T_ERA5L", "N12"):
    d = RES[k]
    ax.plot(BAND[k].index, BAND[k] / d["sig_b"], color=d["col"], lw=1.6 if k != "N12" else 1.2,
            ls=":" if k == "N12" else "-", label=d["nom"])
ax.axhline(1, color=TINTA2, lw=0.9, ls="--", label="umbral 95 % (ruido rojo)")
ax.set_ylabel("varianza en 2–7 años / umbral 95 %")
ax.set_title("(a) Varianza promediada en la banda de 2–7 años: > 1 = la banda está activa por encima del ruido",
             fontsize=10, loc="left")
ax.legend(fontsize=8, frameon=False, ncol=5, loc="upper right")
ax = axs[1]
for k in ("P_L", "Q", "T_ERA5L"):
    r12 = RS.set_index("serie").loc[k, "r_comp_2_7a_con_N12"]
    ax.plot(COMP[k].index, COMP[k], color=RES[k]["col"], lw=1.5, label=f"{RES[k]['nom']} (r con Niño 1+2 = {r12:.2f})")
ax.plot(COMP["N12"].index, COMP["N12"], color=TINTA, lw=1.2, ls=":", label="Niño 1+2")
ax.axhline(0, color=GRIS, lw=0.6)
ax.set_ylabel("componente de 2–7 años (σ)")
ax.set_title("(b) Componente de 2–7 años reconstruida desde las onditas (descriptivo: pocos grados de libertad, sin p)",
             fontsize=10, loc="left")
ax.legend(fontsize=8, frameon=False, ncol=2, loc="upper right")
for a in axs:
    for y in NINOS:
        a.axvline(pd.Timestamp(y, 3, 1), color=GRIS_CLARO, lw=6, alpha=0.4, zorder=0)
    a.grid(alpha=0.25, lw=0.5)
    a.set_xlim(*xlim)
fig.suptitle("Punto 4 · ¿La variabilidad de 2–7 años de la cuenca va con el ENSO? (sombreado: años El Niño)", fontsize=11)
fig.tight_layout()
fig.savefig(OUT / "fig_onditas_banda_4_3.png", dpi=150)
plt.close(fig)
print(f"\nSalidas en {OUT}")
