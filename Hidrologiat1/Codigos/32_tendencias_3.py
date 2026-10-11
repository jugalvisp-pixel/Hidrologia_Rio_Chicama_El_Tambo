"""
Punto 3: tendencias hidroclimáticas mensuales y sus posibles causas (cuenca El Tambo, río Chicama).

Variables (siempre mensuales; ninguna fuente se une con otra):
  P_L      precipitación de referencia, CAMELS-PE `prec` (PISCOp v2.1)       mm/mes   1981-01 a 2025-12
  P_I      IMERG V07 Final                                                 mm/mes   1998-01 a 2025-08
  Q        caudal medio mensual El Tambo (CAMELS-PE, relleno documentado)  m³/s     1998-01 a 2025-12
  T_PISCOt temperatura media de CAMELS-PE `tmean` (PISCOt v1.2)            °C       1981-01 a 2020-12
  T_ERA5L  temperatura a 2 m de ERA5-Land (reanálisis, serie aparte)       °C       1981-01 a 2025-12
  Sensibilidad de fuente: PISCOp v3.0 (1981–2025), CHIRPS v2.0 (1998–2024) y caudal solo con meses observados
  completos (R_flag = completo).

Dos preguntas, dos periodos (guía 3.1):
  · registro completo de cada variable (cambio de largo plazo; IMERG no recorta a las demás);
  · periodo común 1998-01 a 2024-12 (periodo_comun.py) para comparar variables y fuentes; T_PISCOt termina en
    2020-12, así que su comparación con ERA5-Land se hace en 1998–2020 y 1981–2020.

Anomalías (3.2): a_t = X_t − μ_j, z_t = a_t / s_j con μ_j y s_j de un periodo de referencia FIJO, 1998-01 a 2020-12
(23 años): es el periodo más largo en que las cinco series tienen datos, así que todas las fuentes comparten la
misma referencia. Se reporta el número de años por mes. Una climatología fija (no móvil) no remueve tendencia.

Métodos (3.3–3.4), en las tres representaciones y en dos escalas (todos los meses; doce subseries mensuales):
  OLS      X = β0 + β1·t + ε con t en años decimales (fechas reales). Errores estándar HAC de Newey–West
           (12 rezagos en la serie mensual; regla 4(n/100)^(2/9) en las subseries anuales) además de los usuales.
           Para la serie original se compara con OLS con efectos de mes calendario (X = β1·t + Σ γ_j·D_j).
           Diagnóstico de residuos: autocorrelación de rezago 1, Ljung–Box (12), Breusch–Pagan, Jarque–Bera.
  MK/Sen   Mann–Kendall + pendiente de Theil–Sen con intervalo de Gilbert (1987). Serie original: MK estacional
           (Hirsch et al. 1982) y su versión con dependencia entre estaciones (Hirsch & Slack 1984); pendiente
           estacional de Sen (mediana de pendientes dentro de cada mes). Anomalías: MK con varianza corregida por
           autocorrelación (Hamed & Rao 1998). Subseries por mes: MK con Hamed & Rao (también se reporta el MK
           sin corregir). Implementación: pymannkendall (Hussain & Mahmud 2019) + cálculo propio con fechas reales.
  No lineal LOWESS (Cleveland 1979, iteraciones robustas) con ventana de 10 años y sensibilidad 5 y 15 años; banda
           por remuestreo de residuos en bloques móviles de 24 meses (curva EXPLORATORIA, sin inferencia formal).
           Serie original: modelo aditivo X = s_j + f(t) por retroajuste (backfitting) de efectos mensuales y LOWESS.
           Subseries por mes: LOWESS con 75 % de los años cuando hay ≥ 15 años.
           Regresión segmentada (una bisagra) sobre anomalías anuales, comparada con la lineal por BIC; prueba de
           Pettitt para salto de nivel.
Pendientes por década en unidades de la variable (mm/mes/década = cambio del acumulado MENSUAL, no anual;
m³/s/década; °C/década; unidades estandarizadas/década).
Significancia α = 0.05. Comparaciones múltiples: control FDR de Benjamini–Hochberg sobre la familia de las
12 subseries × 5 variables principales (registro completo), por método.
Sensibilidad (3.5): ventanas de inicio/fin (≥ 15 años), retiro de años extremos (1998, 2017, 2023: El Niño costero),
recorte de 5 años en cada extremo, fuente, relleno del caudal. 3.6: índice Niño 1+2 (ERSST v5, CPC/NOAA) para
evaluar si la variabilidad ENSO produce una tendencia aparente; coeficiente de escorrentía y residuo Q–P anual.

Salidas: Tarea 1/El_Tambo/resultados/tendencias/  (CSV y PNG listados al final)
"""
from pathlib import Path
import warnings

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.ticker
from scipy import stats
import statsmodels.api as sm
from statsmodels.nonparametric.smoothers_lowess import lowess
from statsmodels.stats.diagnostic import acorr_ljungbox, het_breuschpagan
from statsmodels.stats.stattools import jarque_bera
from statsmodels.stats.multitest import multipletests
import pymannkendall as pmk

from periodo_comun import periodo_comun

warnings.filterwarnings("ignore", category=RuntimeWarning)
ROOT = Path(__file__).resolve().parents[2] / "Tarea 1" / "El_Tambo"
SM = ROOT / "resultados" / "series_mensuales"
OUT = ROOT / "resultados" / "tendencias"
OUT.mkdir(parents=True, exist_ok=True)
NINO = ROOT / "datos" / "indices_clima" / "ersst5.nino.mth.91-20.ascii"
URL_NINO = "https://www.cpc.ncep.noaa.gov/data/indices/ersst5.nino.mth.91-20.ascii"
PC_INI, PC_FIN = periodo_comun()
REF_INI, REF_FIN = pd.Timestamp("1998-01-01"), pd.Timestamp("2020-12-01")
ALFA = 0.05
HAC_MENSUAL = 12
VENTANA_LOESS = 10                      # años
N_BOOT = 300
BLOQUE = 24                             # meses
ANIOS_EXTREMOS = [1998, 2017, 2023]
MESES = ["Ene", "Feb", "Mar", "Abr", "May", "Jun", "Jul", "Ago", "Sep", "Oct", "Nov", "Dic"]
RNG = np.random.default_rng(32)

# ------------------------------------------------------------------ datos
cam = pd.read_csv(SM / "series_mensuales_camels_pe.csv", index_col="fecha", parse_dates=True)
imerg = pd.read_csv(SM / "imerg_cuenca_mensual.csv", index_col="fecha", parse_dates=True)["P_I_mm_mes"]
t5l = pd.read_csv(SM / "era5land_t2m_cuenca_mensual.csv", index_col="fecha", parse_dates=True)["T_era5land_C"]
v3 = pd.read_csv(SM / "piscop_v3_cuenca_mensual.csv", index_col=0, parse_dates=True)["P_piscop_v3_mm_mes"]
chirps = pd.read_csv(SM / "chirps_cuenca_mensual.csv", index_col=0, parse_dates=True)["P_chirps_mm_mes"]

SERIES = {   # clave: (serie, nombre, unidad, color)
    "P_L": (cam["P_L_mm_mes"], "P_L (CAMELS-PE, PISCOp v2.1)", "mm/mes", "#1f6fb4"),
    "P_I": (imerg, "P_I (IMERG V07 Final)", "mm/mes", "#e4572e"),
    "Q": (cam["Q_m3_s"], "Q̄ El Tambo (relleno documentado)", "m³/s", "#2a9d8f"),
    "T_PISCOt": (cam["T_C"], "T (CAMELS-PE, PISCOt v1.2)", "°C", "#8c510a"),
    "T_ERA5L": (t5l, "T2m (ERA5-Land, reanálisis)", "°C", "#6a4c93"),
}
EXTRA = {
    "P_v3": (v3, "P (PISCOp v3.0)", "mm/mes"),
    "P_CHIRPS": (chirps, "P (CHIRPS v2.0)", "mm/mes"),
    "Q_obs": (cam["Q_m3_s"].where(cam["R_flag"] == "completo"), "Q̄ solo meses observados completos", "m³/s"),
}
PRINC = list(SERIES)
UNID = {k: v[2] for k, v in SERIES.items()} | {k: v[2] for k, v in EXTRA.items()}
SER = {k: v[0].dropna().sort_index() for k, v in SERIES.items()} | {k: v[0].dropna().sort_index() for k, v in EXTRA.items()}


def t_dec(idx):
    return np.asarray(idx.year + (idx.month - 0.5) / 12, float)


def recorte(s, ini=None, fin=None):
    return s.loc[(ini or s.index.min()):(fin or s.index.max())]


# ------------------------------------------------------------------ 3.1 registros y fuentes
def resumen_registro(k, s, src_txt):
    full = pd.date_range(s.index.min(), s.index.max(), freq="MS")
    falt = full.difference(s.index)
    return {"variable": k, "fuente": src_txt, "inicio": f"{s.index.min():%Y-%m}", "fin": f"{s.index.max():%Y-%m}",
            "meses_validos": len(s), "meses_faltantes_internos": len(falt),
            "faltantes": ", ".join(f"{t:%Y-%m}" for t in falt[:12]) + (" …" if len(falt) > 12 else "")}


META = {
    "P_L": "CAMELS-PE v1.0.1 `prec` = PISCOp v2.1 (SENAMHI; pluviómetros + CHIRPS), 0.1°, diaria → suma mensual; "
           "promedio de cuenca de CAMELS-PE",
    "P_I": "GPM IMERG V07 Final mensual (NASA GES DISC), 0.1°, mm/h × horas del mes; media ponderada por área (31 celdas)",
    "Q": "CAMELS-PE `flow_obs` (SENAMHI 211003 El Tambo), diario → media mensual; relleno: interpolación log ≤ 15 d, "
         "escalado ≥ 50 % de días, 1 mes con modelo (1998-08)",
    "T_PISCOt": "CAMELS-PE `tmean` = PISCOt v1.2 (SENAMHI), 0.01°, diaria → media mensual; termina 2020-12",
    "T_ERA5L": "ERA5-Land monthly averaged reanalysis, 2m_temperature (K → °C), 0.1° (~9 km), CDS "
               "doi:10.24381/cds.68d2bb30; media ponderada por área de intersección (UTM 17S)",
    "P_v3": "PISCOp v3.0 diaria (figshare doi:10.6084/m9.figshare.32411886.v1) → suma mensual, pesos por área",
    "P_CHIRPS": "CHIRPS v2.0 mensual 0.05° (IRI), pesos por área; recorte 1998–2024",
    "Q_obs": "Q̄ solo meses con todos los días observados (R_flag = completo)",
}
reg = pd.DataFrame([resumen_registro(k, SER[k], META[k]) for k in SER])
reg["cambios_de_fuente"] = ["ninguno en la serie (PISCOp v2.1 usa una red de estaciones variable en el tiempo)",
                            "TRMM-era (1998–2014, entradas TMPA/GPROF de TRMM) y GPM (2014→) dentro de la misma versión V07",
                            "relleno en 80 meses (76 interpolados, 3 escalados, 1 modelo)",
                            "ninguno; termina en 2020-12 (no se extiende con ERA5-Land en este punto)",
                            "reanálisis; cambia el sistema de observación asimilado por ERA5 a lo largo del tiempo",
                            "producto distinto de PISCOp v2.1", "producto distinto (satélite + estaciones GTS)",
                            "submuestra del caudal"]
# relieve de ERA5-Land frente a la cuenca (resolución y altura, guía 3.1)
try:
    pes = pd.read_csv(ROOT / "datos" / "era5land" / "pesos_era5land_cuenca.csv")
    zm = xr.open_dataset(ROOT / "datos" / "era5land" / "ERA5-Land_geopotencial_El_Tambo.nc")["z"].squeeze(drop=True)
    zm = (zm / 9.80665).rename(latitude="lat", longitude="lon")
    zm = zm.assign_coords(lon=((zm.lon + 180) % 360) - 180)
    zc = np.array([float(zm.sel(lat=la, lon=lo, method="nearest")) for la, lo in zip(pes["lat"], pes["lon"])])
    Z_ERA = float((zc * pes["area_en_cuenca_km2"]).sum() / pes["area_en_cuenca_km2"].sum())
except Exception as e:                                           # el NetCDF de relieve lo genera el script 17
    print("Aviso: no se pudo leer el relieve de ERA5-Land:", e)
    Z_ERA = np.nan
attrs = pd.read_csv(ROOT / "datos" / "camels_pe" / "PE_47E0A2A8_atributos.csv", index_col=0)["valor"]
Z_CUENCA, Z_MAX = float(attrs["elev_mean"]), float(attrs["elev_max"])
reg["nota_altura"] = ""
reg.loc[reg["variable"] == "T_ERA5L", "nota_altura"] = (
    f"relieve medio del modelo {Z_ERA:.0f} m vs cuenca {Z_CUENCA:.0f} m (Δ {Z_ERA - Z_CUENCA:+.0f} m; "
    f"con −6.5 °C/km ≈ {-6.5e-3 * (Z_ERA - Z_CUENCA):+.2f} °C); afecta el nivel, no necesariamente la tendencia")
reg.to_csv(OUT / "registros_fuentes_3_1.csv", index=False)
print("=== 3.1 Registros ===")
print(reg[["variable", "inicio", "fin", "meses_validos", "meses_faltantes_internos"]].to_string(index=False))

# ------------------------------------------------------------------ 3.2 anomalías con referencia fija
CLIM, AN, ZZ = {}, {}, {}
filas = []
for k, s in SER.items():
    ref = s.loc[REF_INI:REF_FIN]
    g = ref.groupby(ref.index.month)
    clim = pd.DataFrame({"n_anios": g.count(), "mu": g.mean(), "s": g.std(ddof=1)})
    CLIM[k] = clim
    mu = s.index.month.map(clim["mu"]).values
    sd = s.index.month.map(clim["s"]).values
    AN[k] = pd.Series(s.values - mu, s.index).dropna()        # meses sin climatología de referencia quedan fuera
    ZZ[k] = pd.Series((s.values - mu) / np.where(sd > 0, sd, np.nan), s.index).dropna()
    sin_ref = sorted(set(range(1, 13)) - set(clim.index[clim["n_anios"] >= 2]))
    if sin_ref:
        print(f"Aviso 3.2: {k} sin climatología de referencia en los meses {sin_ref}: sus anomalías no se definen")
    for m, r in clim.iterrows():
        filas.append({"variable": k, "mes": m, "n_anios": int(r["n_anios"]), "media": r["mu"], "desv_est": r["s"],
                      "cv": r["s"] / r["mu"] if r["mu"] else np.nan})
clim_tab = pd.DataFrame(filas)
clim_tab.round(4).to_csv(OUT / "climatologia_referencia_3_2.csv", index=False)
inest = clim_tab[(clim_tab["variable"].isin(PRINC)) & ((clim_tab["desv_est"] < 1) & clim_tab["variable"].str.startswith("P"))]
print(f"\n=== 3.2 Referencia {REF_INI:%Y-%m}–{REF_FIN:%Y-%m}; meses de P con s_j < 1 mm/mes (z inestable): "
      f"{', '.join(f'{r.variable}-{MESES[r.mes - 1]} ({r.desv_est:.2f})' for r in inest.itertuples()) or 'ninguno'}")
REPS = {"original": SER, "anomalia": AN, "estandarizada": ZZ}
UREP = {"original": None, "anomalia": None, "estandarizada": "z"}


def unidad(k, rep):
    return "z/década" if rep == "estandarizada" else f"{UNID[k]}/década"


# ------------------------------------------------------------------ utilidades estadísticas
def mk_s_var(x):
    n = len(x)
    i, j = np.triu_indices(n, 1)
    s = np.sign(x[j] - x[i]).sum()
    _, c = np.unique(x, return_counts=True)
    var = (n * (n - 1) * (2 * n + 5) - (c * (c - 1) * (2 * c + 5)).sum()) / 18
    return float(s), float(var)


def z_p(s, var):
    z = 0.0 if s == 0 else (s - np.sign(s)) / np.sqrt(var)
    return z, 2 * (1 - stats.norm.cdf(abs(z)))


def pendientes_pares(x, t):
    i, j = np.triu_indices(len(x), 1)
    dt = t[j] - t[i]
    ok = dt > 0
    return (x[j] - x[i])[ok] / dt[ok]


def sen_ci(sl, var):
    sl = np.sort(sl)
    N = len(sl)
    C = stats.norm.ppf(1 - ALFA / 2) * np.sqrt(var)
    m1, m2 = int(np.floor((N - C) / 2)), int(np.ceil((N + C) / 2))
    return np.median(sl), sl[max(m1 - 1, 0)], sl[min(m2, N - 1)]


def mk_sen(s, modificado=True):
    """MK (S y varianza propios, orden cronológico), varianza de Hamed–Rao, Sen con fechas reales [unid/año]."""
    x, t = s.values.astype(float), t_dec(s.index)
    S, var0 = mk_s_var(x)
    out = {"n": len(x), "tau": S / (len(x) * (len(x) - 1) / 2), "S": S}
    out["z_mk"], out["p_mk"] = z_p(S, var0)
    var = var0
    if modificado:
        hr = pmk.hamed_rao_modification_test(x)
        assert abs(hr.s - S) < 1e-6, f"S propio {S} ≠ pymannkendall {hr.s} (n={len(x)}, {s.index.min()})"
        var = hr.var_s if np.isfinite(hr.var_s) and hr.var_s > 0 else var0   # H–R indefinida → MK usual
        out["p_hr"] = z_p(S, var)[1]
        out["factor_var_hr"] = var / var0
    b, lo, hi = sen_ci(pendientes_pares(x, t), var)
    out.update(sen=b * 10, sen_lo=lo * 10, sen_hi=hi * 10)
    return out


def mk_estacional(s):
    """Serie original: MK estacional (Hirsch 1982) y con dependencia entre meses (Hirsch–Slack 1984);
    pendiente estacional de Sen (mediana de pendientes intra-mes) con IC de Gilbert usando la varianza H–S."""
    full = pd.date_range(f"{s.index.min().year}-01-01", f"{s.index.max().year}-12-01", freq="MS")
    x = s.reindex(full).values.astype(float)
    a = pmk.seasonal_test(x, period=12)
    try:                                   # Hirsch–Slack necesita años con todos los meses (falla con muchos vacíos)
        b = pmk.correlated_seasonal_test(x, period=12)
        p_hs, var_hs = b.p, b.var_s
    except ZeroDivisionError:
        p_hs, var_hs = np.nan, a.var_s
    sl = np.concatenate([pendientes_pares(g.values.astype(float), t_dec(g.index))
                         for _, g in s.groupby(s.index.month) if len(g) > 2])
    med, lo, hi = sen_ci(sl, var_hs)
    return {"n": len(s), "tau": a.Tau, "S": a.s, "p_mk": a.p, "p_hr": p_hs, "sen": med * 10, "sen_lo": lo * 10,
            "sen_hi": hi * 10, "factor_var_hr": var_hs / a.var_s}


def ols(s, hac, mes=False):
    t = t_dec(s.index)
    X = pd.DataFrame({"t": t - t.mean()}, index=s.index)
    if mes:
        X = X.join(pd.get_dummies(s.index.month, prefix="m", dtype=float).set_index(s.index))
    else:
        X = sm.add_constant(X)
    y = s.values.astype(float)
    r0 = sm.OLS(y, X).fit()
    r = sm.OLS(y, X).fit(cov_type="HAC", cov_kwds={"maxlags": hac})
    res = r0.resid
    ci0, ci = r0.conf_int(ALFA).loc["t"], r.conf_int(ALFA).loc["t"]
    out = {"n": len(y), "pend": r.params["t"] * 10, "lo_usual": ci0.iloc[0] * 10, "hi_usual": ci0.iloc[1] * 10,
           "p_usual": r0.pvalues["t"], "lo": ci.iloc[0] * 10, "hi": ci.iloc[1] * 10, "p": r.pvalues["t"],
           "r2": r0.rsquared, "r1_resid": pd.Series(np.asarray(res)).autocorr(1)}
    if len(y) > 30:
        out["p_ljungbox12"] = float(acorr_ljungbox(res, lags=[12])["lb_pvalue"].iloc[0])
    out["p_breuschpagan"] = het_breuschpagan(res, sm.add_constant(X["t"].values))[1]   # ¿cambia la dispersión con t?
    out["p_jarquebera"] = jarque_bera(res)[1]
    return out


def hac_anual(n):
    return int(np.floor(4 * (n / 100) ** (2 / 9)))


def loess(s, anios, it=3):
    t = t_dec(s.index)
    frac = min(1.0, anios * 12 / len(s)) if anios else 0.75
    return lowess(s.values, t, frac=frac, it=it, return_sorted=False)


def loess_banda(s, anios=VENTANA_LOESS):
    f = loess(s, anios)
    res = s.values - f
    n = len(res)
    reps = np.empty((N_BOOT, n))
    for b in range(N_BOOT):
        idx = np.concatenate([np.arange(st, st + BLOQUE) for st in RNG.integers(0, n - BLOQUE, n // BLOQUE + 1)])[:n]
        reps[b] = loess(pd.Series(f + res[idx], s.index), anios)
    return f, np.percentile(reps, 5, axis=0), np.percentile(reps, 95, axis=0)


def aditivo(s, anios=VENTANA_LOESS, n_it=10):
    """X = s_j + f(t): retroajuste de efecto mensual (media de X − f por mes) y LOWESS de X − s_j."""
    f = np.zeros(len(s))
    m = s.index.month
    for _ in range(n_it):
        sj = pd.Series(s.values - f).groupby(m).transform("mean").values
        f = loess(pd.Series(s.values - sj, s.index), anios)
        f = f - f.mean()
    return sj, f


def tramos(t, f, unidad_):
    """Intervalos de aumento/disminución de una curva suave (signo de la derivada)."""
    d = np.sign(np.diff(f))
    d[d == 0] = 1
    out, ini = [], 0
    for i in range(1, len(d) + 1):
        if i == len(d) or d[i] != d[ini]:
            out.append({"desde": t[ini], "hasta": t[i], "sentido": "aumento" if d[ini] > 0 else "disminución",
                        "cambio": f[i] - f[ini], "duracion_anios": t[i] - t[ini], "unidad": unidad_})
            ini = i
    return out


def anual(s, min_meses=10):
    g = s.groupby(s.index.year)
    a = g.mean()
    return a[g.count() >= min_meses]


def pettitt(x):
    x = np.asarray(x, float)
    n = len(x)
    U = np.array([np.sign(x[k + 1:][None, :] - x[:k + 1][:, None]).sum() for k in range(n - 1)])
    k = int(np.argmax(np.abs(U)))
    K = abs(U[k])
    p = min(1.0, 2 * np.exp(-6 * K ** 2 / (n ** 3 + n ** 2)))
    return k, p


def segmentada(y):
    """Lineal vs bisagra continua (una ruptura de pendiente) sobre anomalías anuales; BIC."""
    t = y.index.values.astype(float)
    v = y.values
    n = len(v)
    X1 = np.column_stack([np.ones(n), t])
    rss1 = np.sum((v - X1 @ np.linalg.lstsq(X1, v, rcond=None)[0]) ** 2)
    bic1 = n * np.log(rss1 / n) + 2 * np.log(n)
    mejor = (np.inf, None, None)
    for tb in t[5:-5]:
        X2 = np.column_stack([X1, np.clip(t - tb, 0, None)])
        b = np.linalg.lstsq(X2, v, rcond=None)[0]
        rss = np.sum((v - X2 @ b) ** 2)
        if rss < mejor[0]:
            mejor = (rss, tb, b)
    rss2, tb, b = mejor
    bic2 = n * np.log(rss2 / n) + 4 * np.log(n)               # 3 coeficientes + ruptura
    return {"anio_ruptura": int(tb), "pend_antes": b[1] * 10, "pend_despues": (b[1] + b[2]) * 10,
            "BIC_lineal": bic1, "BIC_segmentada": bic2, "prefiere": "segmentada" if bic2 < bic1 - 2 else "lineal"}


# ------------------------------------------------------------------ 3.3–3.4 todos los meses
PERIODOS = {"completo": (None, None), "comun": (PC_INI, PC_FIN)}
tod = []
for k in SER:
    for per, (i0, i1) in PERIODOS.items():
        for rep, D in REPS.items():
            s = recorte(D[k], i0, i1)
            fila0 = {"variable": k, "periodo": per, "inicio": f"{s.index.min():%Y-%m}", "fin": f"{s.index.max():%Y-%m}",
                     "representacion": rep, "unidad": unidad(k, rep)}
            o = ols(s, HAC_MENSUAL)
            tod.append(fila0 | {"metodo": "OLS (HAC 12)"} | o)
            if rep == "original":
                tod.append(fila0 | {"metodo": "OLS + efectos de mes (HAC 12)"} | ols(s, HAC_MENSUAL, mes=True))
                me = mk_estacional(s)
                tod.append(fila0 | {"metodo": "MK estacional + Sen estacional", "pend": me["sen"], "lo": me["sen_lo"],
                                    "hi": me["sen_hi"], "p": me["p_hr"], "p_usual": me["p_mk"], "n": me["n"],
                                    "tau": me["tau"], "factor_var_hr": me["factor_var_hr"]})
            else:
                m = mk_sen(s)
                tod.append(fila0 | {"metodo": "MK Hamed-Rao + Sen", "pend": m["sen"], "lo": m["sen_lo"], "hi": m["sen_hi"],
                                    "p": m["p_hr"], "p_usual": m["p_mk"], "n": m["n"], "tau": m["tau"],
                                    "factor_var_hr": m["factor_var_hr"]})
TOD = pd.DataFrame(tod)
TOD.to_csv(OUT / "tendencias_todos_meses_3_4.csv", index=False, float_format="%.5g")
print("\n=== 3.3–3.4 Todos los meses (anomalías; pendiente por década) ===")
v = TOD[(TOD["representacion"] == "anomalia")][["variable", "periodo", "metodo", "pend", "lo", "hi", "p", "n"]]
print(v.round(3).to_string(index=False))

# ------------------------------------------------------------------ 3.3–3.4 subseries mensuales
mes_f, verif = [], []
for k in SER:
    for per, (i0, i1) in PERIODOS.items():
        for m in range(1, 13):
            sub = {rep: recorte(D[k], i0, i1) for rep, D in REPS.items()}
            sub = {rep: s[s.index.month == m] for rep, s in sub.items()}
            n = len(sub["original"])
            sj = CLIM[k]["s"].get(m, np.nan)
            if n < 8 or not np.isfinite(sj):
                continue
            o = {rep: ols(s, hac_anual(n)) for rep, s in sub.items()}
            mk = {rep: mk_sen(s) for rep, s in sub.items()}
            verif.append({"variable": k, "periodo": per, "mes": m, "n": n,
                          "dif_pend_OLS_X_menos_a": o["original"]["pend"] - o["anomalia"]["pend"],
                          "dif_pend_OLS_z_menos_a_sobre_s": o["estandarizada"]["pend"] - o["anomalia"]["pend"] / sj,
                          "dif_p_OLS_X_z": o["original"]["p"] - o["estandarizada"]["p"],
                          "dif_tau_MK_X_z": mk["original"]["tau"] - mk["estandarizada"]["tau"],
                          "dif_sen_X_menos_a": mk["original"]["sen"] - mk["anomalia"]["sen"]})
            for rep in REPS:
                base = {"variable": k, "periodo": per, "mes": m, "representacion": rep, "unidad": unidad(k, rep),
                        "n": n, "inicio": sub[rep].index.min().year, "fin": sub[rep].index.max().year,
                        "s_j_ref": sj, "r1_subserie": sub[rep].autocorr(1)}
                mes_f.append(base | {"metodo": "OLS (HAC)", "pend": o[rep]["pend"], "lo": o[rep]["lo"], "hi": o[rep]["hi"],
                                     "p": o[rep]["p"], "p_usual": o[rep]["p_usual"], "r1_resid": o[rep]["r1_resid"],
                                     "p_breuschpagan": o[rep]["p_breuschpagan"]})
                mes_f.append(base | {"metodo": "MK Hamed-Rao + Sen", "pend": mk[rep]["sen"], "lo": mk[rep]["sen_lo"],
                                     "hi": mk[rep]["sen_hi"], "p": mk[rep]["p_hr"], "p_usual": mk[rep]["p_mk"],
                                     "tau": mk[rep]["tau"]})
MES = pd.DataFrame(mes_f)
VER = pd.DataFrame(verif)
VER.to_csv(OUT / "verificacion_transformaciones_3_3.csv", index=False, float_format="%.3g")
print(f"\nVerificación (máx |dif|): OLS X−a {VER['dif_pend_OLS_X_menos_a'].abs().max():.2e}, "
      f"OLS z − a/s {VER['dif_pend_OLS_z_menos_a_sobre_s'].abs().max():.2e}, p OLS X−z {VER['dif_p_OLS_X_z'].abs().max():.2e}, "
      f"τ MK X−z {VER['dif_tau_MK_X_z'].abs().max():.2e}, Sen X−a {VER['dif_sen_X_menos_a'].abs().max():.2e}")

# FDR: familia = 12 meses × 5 variables principales, registro completo, por método (las 3 representaciones dan
# la misma prueba, así que se usa una sola: la original).
MES["q_fdr"] = np.nan
fam = (MES["periodo"] == "completo") & MES["variable"].isin(PRINC) & (MES["representacion"] == "original")
for met in MES["metodo"].unique():
    sel = fam & (MES["metodo"] == met)
    MES.loc[sel, "q_fdr"] = multipletests(MES.loc[sel, "p"], alpha=ALFA, method="fdr_bh")[1]
fam_c = (MES["periodo"] == "comun") & MES["variable"].isin(PRINC) & (MES["representacion"] == "original")
for met in MES["metodo"].unique():
    sel = fam_c & (MES["metodo"] == met)
    MES.loc[sel, "q_fdr"] = multipletests(MES.loc[sel, "p"], alpha=ALFA, method="fdr_bh")[1]
for col in ("q_fdr",):                                         # copiar q a las otras representaciones
    key = ["variable", "periodo", "mes", "metodo"]
    q = MES[MES["representacion"] == "original"].set_index(key)["q_fdr"]
    MES["q_fdr"] = MES.set_index(key).index.map(q)
MES.to_csv(OUT / "tendencias_por_mes_3_4.csv", index=False, float_format="%.5g")
sig = MES[(MES["representacion"] == "original") & (MES["p"] < ALFA)]
print(f"\n=== Subseries con p < {ALFA} (antes de FDR) ===")
print(sig[["variable", "periodo", "mes", "metodo", "pend", "lo", "hi", "p", "q_fdr", "n"]].round(4).to_string(index=False))

# ------------------------------------------------------------------ 3.4 no lineal: LOWESS, aditivo, segmentada
curvas, tr, seg = {}, [], []
for k in PRINC:
    s_a, s_z, s_x = AN[k], ZZ[k], SER[k]
    fa, la, ha = loess_banda(s_a)
    fz = loess(s_z, VENTANA_LOESS)
    sj, fx = aditivo(s_x)
    sens = {w: loess(s_a, w) for w in (5, 15)}
    sens_lin = loess(s_a, VENTANA_LOESS, it=0)                   # sin iteraciones robustas (efecto de extremos)
    fo = ols(s_a, HAC_MENSUAL)
    curvas[k] = dict(fa=fa, la=la, ha=ha, fz=fz, fx=fx, sj=sj, sens=sens, nr=sens_lin)
    tr += [{"variable": k, "representacion": "anomalia", "ventana_anios": VENTANA_LOESS} | d
           for d in tramos(t_dec(s_a.index), fa, UNID[k])]
    ya = anual(s_a)
    sg = segmentada(ya)
    kp, pp = pettitt(ya.values)
    yrs = ya.index.values
    seg.append({"variable": k, "anios": f"{yrs[0]}–{yrs[-1]}", "n_anios": len(yrs)} | sg |
               {"pettitt_anio_cambio": int(yrs[kp + 1]), "pettitt_p": pp,
                "media_antes": ya.iloc[:kp + 1].mean(), "media_despues": ya.iloc[kp + 1:].mean(),
                "dif_lowess_vs_ols_max": float(np.max(np.abs(fa - (fo["pend"] / 10) * (t_dec(s_a.index) - t_dec(s_a.index).mean())
                                                               - fa.mean()))),
                "dif_aditivo_vs_lowess_anom_max": float(np.max(np.abs((fx - fx.mean()) - (fa - fa.mean())))),
                "amplitud_lowess": float(fa.max() - fa.min())})
TR = pd.DataFrame(tr)
TR = TR[TR["duracion_anios"] >= 1.0]
TR.to_csv(OUT / "tramos_lowess_3_4.csv", index=False, float_format="%.4g")
SEG = pd.DataFrame(seg)
SEG.to_csv(OUT / "segmentada_pettitt_3_5.csv", index=False, float_format="%.4g")
print("\n=== Segmentada (anomalías anuales) y Pettitt ===")
print(SEG[["variable", "anios", "anio_ruptura", "pend_antes", "pend_despues", "prefiere", "pettitt_anio_cambio",
           "pettitt_p"]].round(3).to_string(index=False))

# ------------------------------------------------------------------ 3.5 sensibilidad
sen_f = []


def fila_sens(k, caso, s, rep="anomalia"):
    if len(s) < 60:
        return
    m = mk_sen(s)
    o = ols(s, HAC_MENSUAL)
    sen_f.append({"variable": k, "caso": caso, "inicio": f"{s.index.min():%Y-%m}", "fin": f"{s.index.max():%Y-%m}",
                  "n": len(s), "unidad": unidad(k, rep), "sen": m["sen"], "sen_lo": m["sen_lo"], "sen_hi": m["sen_hi"],
                  "p_mk_hr": m["p_hr"], "ols": o["pend"], "ols_lo": o["lo"], "ols_hi": o["hi"], "p_ols_hac": o["p"]})


for k in list(SER):
    a = AN[k]
    for per, (i0, i1) in PERIODOS.items():
        fila_sens(k, f"base ({'común' if per == 'comun' else per})", recorte(a, i0, i1))
    fila_sens(k, "sin 1998, 2017, 2023 (extremos)", a[~a.index.year.isin(ANIOS_EXTREMOS)])
    y0, y1 = a.index.min().year, a.index.max().year
    fila_sens(k, "inicio +5 años", a.loc[f"{y0 + 5}":])
    fila_sens(k, "fin −5 años", a.loc[:f"{y1 - 5}"])
    fila_sens(k, "1998–2020 (periodo de T_PISCOt)", a.loc["1998":"2020"])
    fila_sens(k, "1981–2020 (registro de T_PISCOt)", a.loc["1981":"2020"])
    fila_sens(k, "z, registro completo", ZZ[k], rep="estandarizada")
SENS = pd.DataFrame(sen_f).drop_duplicates(subset=["variable", "inicio", "fin", "n", "unidad"])
SENS.to_csv(OUT / "sensibilidad_3_5.csv", index=False, float_format="%.4g")
print("\n=== Sensibilidad (anomalías, Sen y OLS-HAC por década) ===")
print(SENS[["variable", "caso", "inicio", "fin", "sen", "sen_lo", "sen_hi", "p_mk_hr", "ols", "p_ols_hac"]]
      .round(3).to_string(index=False))

# Ventanas de inicio y fin (anomalías anuales, Sen y MK; ≥ 15 años)
VENT = {}
for k in PRINC:
    ya = anual(AN[k])
    yrs = ya.index.values
    M = np.full((len(yrs), len(yrs)), np.nan)
    P = np.full_like(M, np.nan)
    for i, a0 in enumerate(yrs):
        for j, a1 in enumerate(yrs):
            if a1 - a0 + 1 >= 15:
                w = ya.loc[a0:a1]
                M[i, j] = np.median(pendientes_pares(w.values, w.index.values.astype(float))) * 10
                P[i, j] = z_p(*mk_s_var(w.values))[1]
    VENT[k] = (yrs, M, P)

# ------------------------------------------------------------------ 3.6 ENSO, escorrentía
if not NINO.exists():
    import requests
    NINO.parent.mkdir(parents=True, exist_ok=True)
    NINO.write_bytes(requests.get(URL_NINO, timeout=60).content)
n12 = pd.read_csv(NINO, sep=r"\s+")
n12.index = pd.to_datetime(dict(year=n12["YR"], month=n12["MON"], day=1))
n12 = n12.iloc[:, 3].rename("Nino12_anom")                      # columna ANOM de NINO1+2 (base 1991–2020)
enso = []
for k in list(SER):
    a = AN[k]
    d = pd.concat([a.rename("a"), n12, n12.shift(1).rename("n12_1")], axis=1).dropna()
    X = sm.add_constant(d[["Nino12_anom", "n12_1"]])
    r = sm.OLS(d["a"], X).fit(cov_type="HAC", cov_kwds={"maxlags": HAC_MENSUAL})
    resid = pd.Series(d["a"].values - r.fittedvalues.values + r.params["const"], d.index)
    m0, m1 = mk_sen(d["a"]), mk_sen(resid)
    mn = mk_sen(n12.loc[d.index.min():d.index.max()])
    enso.append({"variable": k, "inicio": f"{d.index.min():%Y-%m}", "fin": f"{d.index.max():%Y-%m}", "n": len(d),
                 "r_Nino12": d["a"].corr(d["Nino12_anom"]), "r_Nino12_rezago1": d["a"].corr(d["n12_1"]),
                 "R2_enso": r.rsquared, "sen_anomalia": m0["sen"], "p_anomalia": m0["p_hr"],
                 "sen_residuo_sin_enso": m1["sen"], "lo_residuo": m1["sen_lo"], "hi_residuo": m1["sen_hi"],
                 "p_residuo": m1["p_hr"], "sen_Nino12_mismo_periodo_C_dec": mn["sen"], "p_Nino12": mn["p_hr"],
                 "unidad": unidad(k, "anomalia")})
ENSO = pd.DataFrame(enso)
ENSO.to_csv(OUT / "enso_3_6.csv", index=False, float_format="%.4g")
print("\n=== ENSO (Niño 1+2) ===")
print(ENSO[["variable", "inicio", "fin", "r_Nino12", "R2_enso", "sen_anomalia", "p_anomalia", "sen_residuo_sin_enso",
            "p_residuo", "sen_Nino12_mismo_periodo_C_dec"]].round(3).to_string(index=False))

# Respuesta del caudal a la lluvia: año hidrológico sep–ago completo (R en mm, P_L en mm)
R = cam["R_mm_mes"]
ah = pd.DataFrame({"P": cam["P_L_mm_mes"], "R": R, "T": t5l})
ah["ah"] = ah.index.year + (ah.index.month >= 9)
g = ah.groupby("ah")
anu = pd.DataFrame({"P": g["P"].sum(min_count=12), "R": g["R"].sum(min_count=12), "T": g["T"].mean()}).dropna()
anu["C"] = anu["R"] / anu["P"]
rr = sm.OLS(anu["R"], sm.add_constant(anu["P"])).fit()
anu["resid_R_dado_P"] = rr.resid
esc = []
for c, u in (("P", "mm/año"), ("R", "mm/año"), ("C", "–"), ("resid_R_dado_P", "mm/año"), ("T", "°C")):
    sl = pendientes_pares(anu[c].values, anu.index.values.astype(float))
    S, var = mk_s_var(anu[c].values)
    b, lo, hi = sen_ci(sl, var)
    esc.append({"serie_anual": c, "anios_hidrologicos": f"{anu.index.min() - 1}/{anu.index.min()}–"
                f"{anu.index.max() - 1}/{anu.index.max()}", "n": len(anu), "sen_dec": b * 10, "lo": lo * 10,
                "hi": hi * 10, "p_mk": z_p(S, var)[1], "unidad": u + "/década"})
ESC = pd.DataFrame(esc)
ESC.to_csv(OUT / "escorrentia_anual_3_6.csv", index=False, float_format="%.4g")
anu.round(4).to_csv(OUT / "balance_anual_hidrologico_3_6.csv")
print(f"\n=== Año hidrológico (R = {rr.params['const']:.0f} + {rr.params['P']:.3f}·P, r² {rr.rsquared:.2f}) ===")
print(ESC.round(3).to_string(index=False))

# ------------------------------------------------------------------ figuras
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
GRIS = "#555555"


def fmt_eje_tiempo(ax):
    ax.xaxis.set_major_locator(mdates.YearLocator(5))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    ax.grid(alpha=0.25, lw=0.5)


# Figura 1: tres representaciones con OLS, Sen y LOWESS (registro completo)
fig, axs = plt.subplots(len(PRINC), 3, figsize=(17, 3.0 * len(PRINC)), sharex=True)
for i, k in enumerate(PRINC):
    col = SERIES[k][3]
    c = curvas[k]
    for j, rep in enumerate(REPS):
        ax = axs[i, j]
        s = REPS[rep][k]
        t = t_dec(s.index)
        ax.plot(s.index, s.values, color=col, lw=0.6, alpha=0.55)
        sub = TOD[(TOD["variable"] == k) & (TOD["periodo"] == "completo") & (TOD["representacion"] == rep)]
        if rep == "original":
            ax.plot(s.index, s.values.mean() + c["fx"], color="k", lw=1.6, label="aditivo s_j + LOWESS 10 a (tendencia f(t) + media)")
            o = sub[sub["metodo"] == "OLS + efectos de mes (HAC 12)"].iloc[0]
            ax.plot(s.index, s.values.mean() + o["pend"] / 10 * (t - t.mean()), color="#d62728", lw=1.3, ls="--",
                    label=f"OLS + efectos de mes: {o['pend']:+.2f} [{o['lo']:+.2f}, {o['hi']:+.2f}]")
            me = sub[sub["metodo"] == "MK estacional + Sen estacional"].iloc[0]
            ax.set_title(f"{SERIES[k][1]} — original\nSen estacional {me['pend']:+.2f} {unidad(k, rep)} "
                         f"(p H–S {me['p']:.2f})", fontsize=8.5, loc="left")
        else:
            f = c["fa"] if rep == "anomalia" else c["fz"]
            if rep == "anomalia":
                ax.fill_between(s.index, c["la"], c["ha"], color="k", alpha=0.15, lw=0, label="banda 90 % (bloques 24 m)")
            ax.plot(s.index, f, color="k", lw=1.6, label="LOWESS 10 años")
            o = sub[sub["metodo"] == "OLS (HAC 12)"].iloc[0]
            m = sub[sub["metodo"] == "MK Hamed-Rao + Sen"].iloc[0]
            ax.plot(s.index, o["pend"] / 10 * (t - t.mean()) + np.mean(s.values), color="#d62728", lw=1.3, ls="--",
                    label=f"OLS-HAC {o['pend']:+.2f} [{o['lo']:+.2f}, {o['hi']:+.2f}] p={o['p']:.2f}")
            med_t = np.median(t)
            ax.plot(s.index, m["pend"] / 10 * (t - med_t) + np.median(s.values), color="#1b9e77", lw=1.3, ls=":",
                    label=f"Sen {m['pend']:+.2f} [{m['lo']:+.2f}, {m['hi']:+.2f}] p H–R={m['p']:.2f}")
            ax.axhline(0, color=GRIS, lw=0.6)
            ax.set_title(f"{SERIES[k][1]} — {'anomalía a' if rep == 'anomalia' else 'anomalía estandarizada z'}"
                         f"\npendientes en {unidad(k, rep)}", fontsize=8.5, loc="left")
        if k == "Q" and rep != "estandarizada":   # la crecida de 1998 (~700 m³/s) aplana el panel: eje recortado
            lo_, hi_ = np.nanpercentile(s.values, [0, 99])
            fuera = int((s.values > hi_ * 1.15).sum())
            ax.set_ylim(lo_ - 0.05 * (hi_ - lo_), hi_ * 1.15 if hi_ > 0 else hi_ + 1)
            ax.text(0.99, 0.97, f"eje recortado: {fuera} meses fuera (máx {s.max():.0f})", transform=ax.transAxes,
                    ha="right", va="top", fontsize=6.8, color=GRIS)
        ax.legend(fontsize=6.6, loc="upper left", frameon=False)
        ax.set_ylabel(UNID[k] if rep != "estandarizada" else "z")
        fmt_eje_tiempo(ax)
fig.suptitle(f"Punto 3 · Series originales, anomalías y anomalías estandarizadas (registro completo de cada variable; "
             f"referencia {REF_INI:%Y}–{REF_FIN:%Y})", fontsize=11, y=0.995)
fig.tight_layout()
fig.savefig(OUT / "fig_series_tendencias_3_4.png", dpi=150)
plt.close(fig)

# Figura 2: pendientes por mes calendario e intervalos (registro completo), unidades originales y z
fig, axs = plt.subplots(len(PRINC), 2, figsize=(14, 2.7 * len(PRINC)))
x = np.arange(1, 13)
for i, k in enumerate(PRINC):
    for j, rep in enumerate(("original", "estandarizada")):
        ax = axs[i, j]
        for dx, met, mk_, col in ((-0.17, "OLS (HAC)", "s", "#d62728"), (0.17, "MK Hamed-Rao + Sen", "o", "#1b9e77")):
            d = MES[(MES["variable"] == k) & (MES["periodo"] == "completo") & (MES["representacion"] == rep)
                    & (MES["metodo"] == met)].set_index("mes").reindex(x)
            ax.errorbar(x + dx, d["pend"], yerr=[d["pend"] - d["lo"], d["hi"] - d["pend"]], fmt="none", ecolor=col,
                        elinewidth=1.1, capsize=0)
            q = d["q_fdr"] < ALFA
            p = (d["p"] < ALFA) & ~q
            ax.scatter((x + dx)[~(q | p)], d["pend"][~(q | p)], marker=mk_, s=26, facecolor="white", edgecolor=col, zorder=3,
                       label=f"{met}" if i == 0 else None)
            ax.scatter((x + dx)[p], d["pend"][p], marker=mk_, s=26, facecolor=col, edgecolor=col, alpha=0.45, zorder=3)
            ax.scatter((x + dx)[q], d["pend"][q], marker=mk_, s=40, facecolor=col, edgecolor="k", zorder=4)
        ax.axhline(0, color=GRIS, lw=0.7)
        ax.set_xticks(x, MESES)
        n = MES[(MES["variable"] == k) & (MES["periodo"] == "completo")]["n"]
        nn = f"{n.min()}" if n.min() == n.max() else f"{n.min()}–{n.max()}"
        ax.set_title(f"{SERIES[k][1]} · {unidad(k, rep)} · {nn} años por mes", fontsize=8.5, loc="left")
        ax.grid(alpha=0.25, lw=0.5, axis="y")
axs[0, 0].legend(fontsize=7, frameon=False, loc="upper left")
fig.text(0.5, 0.003, f"Barras: IC {100 * (1 - ALFA):.0f} % (OLS con HAC; Sen con Gilbert y varianza Hamed–Rao). Marcador vacío: "
         f"p ≥ {ALFA}; relleno claro: p < {ALFA} pero q FDR ≥ {ALFA};\nrelleno con borde negro: q FDR < {ALFA} "
         f"(familia 12 meses × 5 variables). Las anomalías a tienen las mismas pendientes que la serie original.",
         ha="center", fontsize=8)
fig.suptitle("Punto 3 · Pendiente de cada mes calendario a través de los años (registro completo)", fontsize=11)
fig.tight_layout(rect=(0, 0.03, 1, 0.98))
fig.savefig(OUT / "fig_pendientes_mensuales_3_5.png", dpi=150)
plt.close(fig)

# Figura 3: subseries mensuales de cada variable (12 paneles) con OLS, Sen y LOWESS
for k in PRINC:
    fig, axs = plt.subplots(4, 3, figsize=(13, 11), sharex=True)
    col = SERIES[k][3]
    for m in range(1, 13):
        ax = axs.flat[m % 12]                                      # filas DEF / MAM / JJA / SON
        s = SER[k][SER[k].index.month == m]
        yr = s.index.year.values.astype(float)
        ax.plot(yr, s.values, "o-", color=col, ms=3, lw=0.7, alpha=0.8)
        d = MES[(MES["variable"] == k) & (MES["periodo"] == "completo") & (MES["mes"] == m)
                & (MES["representacion"] == "original")].set_index("metodo")
        if len(d) == 0:
            continue
        o, mk_ = d.loc["OLS (HAC)"], d.loc["MK Hamed-Rao + Sen"]
        ax.plot(yr, s.mean() + o["pend"] / 10 * (yr - yr.mean()), color="#d62728", ls="--", lw=1.2)
        ax.plot(yr, np.median(s.values) + mk_["pend"] / 10 * (yr - np.median(yr)), color="#1b9e77", ls=":", lw=1.4)
        if len(s) >= 15:
            ax.plot(yr, lowess(s.values, yr, frac=0.75, return_sorted=False), color="k", lw=1.3)
        star = "**" if mk_["q_fdr"] < ALFA else ("*" if mk_["p"] < ALFA else "")
        ax.set_title(f"{MESES[m - 1]}: Sen {mk_['pend']:+.2f}{star} · OLS {o['pend']:+.2f} {UNID[k]}/déc",
                     fontsize=8, loc="left")
        ax.grid(alpha=0.25, lw=0.5)
    for i, tri in enumerate(("DEF", "MAM", "JJA", "SON")):
        axs[i, 0].set_ylabel(f"{tri}\n{UNID[k]}", fontsize=9)
    fig.suptitle(f"Punto 3 · {SERIES[k][1]}: cada mes calendario a través de los años ({UNID[k]}). "
                 f"\nRojo discontinuo OLS, verde punteado Sen, negro LOWESS (75 % de los años, exploratorio). "
                 f"* p MK–HR < {ALFA}; ** q FDR < {ALFA}", fontsize=9.5)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(OUT / f"fig_subseries_{k}_3_4.png", dpi=140)
    plt.close(fig)

# Figura 4: robustez — comparación de pendientes de anomalías (métodos, periodos, fuentes, casos)
grupos = [("Precipitación (mm/mes/década)", ["P_L", "P_I", "P_v3", "P_CHIRPS"]),
          ("Caudal (m³/s/década)", ["Q", "Q_obs"]), ("Temperatura (°C/década)", ["T_PISCOt", "T_ERA5L"])]
COL_EXTRA = {"P_v3": "#5b7d9a", "P_CHIRPS": "#b07aa1", "Q_obs": "#7a9e3a"}
nfil = [SENS["variable"].isin(ks).sum() - (SENS["variable"].isin(ks) & (SENS["unidad"] == "z/década")).sum()
        for _, ks in grupos]
fig, axs = plt.subplots(3, 1, figsize=(12, 0.3 * sum(nfil) + 2.5), gridspec_kw={"height_ratios": nfil})
for ax, (tit, ks) in zip(axs, grupos):
    rows = SENS[SENS["variable"].isin(ks) & (SENS["unidad"] != "z/década")].copy()
    rows["orden"] = rows["variable"].map({v: i for i, v in enumerate(ks)})
    rows = rows.sort_values(["orden", "caso"]).reset_index(drop=True)
    y = np.arange(len(rows))[::-1]
    for yi, r in zip(y, rows.itertuples()):
        col = SERIES[r.variable][3] if r.variable in SERIES else COL_EXTRA[r.variable]
        ax.plot([r.sen_lo, r.sen_hi], [yi + 0.15] * 2, color=col, lw=1.6)
        ax.plot(r.sen, yi + 0.15, "o", color=col, mfc=col if r.p_mk_hr < ALFA else "white", ms=5)
        ax.plot([r.ols_lo, r.ols_hi], [yi - 0.15] * 2, color=col, lw=1.0, alpha=0.6)
        ax.plot(r.ols, yi - 0.15, "s", color=col, mfc=col if r.p_ols_hac < ALFA else "white", ms=4, alpha=0.8)
    ax.set_yticks(y, [f"{r.variable} · {r.caso} · {r.inicio[:4]}–{r.fin[:4]}" for r in rows.itertuples()], fontsize=7)
    ax.axvline(0, color=GRIS, lw=0.7)
    ax.set_title(tit, fontsize=9.5, loc="left")
    ax.grid(alpha=0.25, lw=0.5, axis="x")
fig.suptitle("Punto 3.5 · Robustez de la tendencia de las anomalías: ● Sen + MK Hamed–Rao (arriba), ■ OLS-HAC (abajo); "
             f"relleno = p < {ALFA}; barras IC {100 * (1 - ALFA):.0f} %", fontsize=10)
fig.tight_layout(rect=(0, 0, 1, 0.96))
fig.savefig(OUT / "fig_robustez_3_5.png", dpi=150)
plt.close(fig)

# Figura 5: ventanas de inicio y fin
fig, axs = plt.subplots(1, len(PRINC), figsize=(4.1 * len(PRINC), 4.4))
for ax, k in zip(axs, PRINC):
    yrs, M, P = VENT[k]
    lim = np.nanmax(np.abs(M))
    im = ax.pcolormesh(yrs, yrs, M.T, cmap="RdBu_r" if k.startswith("T") else "BrBG", vmin=-lim, vmax=lim, shading="nearest")
    ii, jj = np.where(P < ALFA)
    ax.scatter(yrs[ii], yrs[jj], s=3, c="k", marker=".")
    ax.set_xlabel("año inicial")
    ax.set_ylabel("año final")
    ax.set_title(f"{k}: Sen de anomalías anuales\n({UNID[k]}/década; • p MK < {ALFA}, sin corregir)", fontsize=8.5)
    plt.colorbar(im, ax=ax, shrink=0.8)
    ax.set_xlim(yrs[0] - 0.5, yrs[-1] - 13.5)
    ax.set_ylim(yrs[0] + 13.5, yrs[-1] + 0.5)
    ax.grid(alpha=0.25, lw=0.5)
    ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True, nbins=5))
    ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True, nbins=6))
    ax.set_aspect("equal")
fig.suptitle("Punto 3.5 · Sensibilidad de la pendiente a la fecha inicial y final (ventanas ≥ 15 años)", fontsize=10.5)
fig.tight_layout()
fig.savefig(OUT / "fig_ventanas_inicio_fin_3_5.png", dpi=150)
plt.close(fig)

# Figura 6: 3.6 — anomalías anuales vs Niño 1+2, balance anual
fig, axs = plt.subplots(4, 1, figsize=(13, 12), sharex=True, gridspec_kw={"height_ratios": [0.7, 1.1, 1, 0.7]})
n12a = anual(n12.loc["1981":])
ax = axs[0]
ax.bar(n12a.index, n12a.values, color=np.where(n12a.values > 0, "#d6604d", "#4393c3"), width=0.8)
ax.axhline(0, color=GRIS, lw=0.6)
ax.set_ylabel("°C")
ax.set_title("(a) Niño 1+2 anual (ERSST v5, CPC/NOAA; anomalías respecto a 1991–2020)", fontsize=9.5, loc="left")
ax = axs[1]
for k in ("P_L", "Q", "T_ERA5L"):
    z = anual(ZZ[k])
    ax.plot(z.index, z.values, "o-", ms=3, lw=1.2, color=SERIES[k][3], label=k)
ax.axhline(0, color=GRIS, lw=0.6)
ax.set_ylabel("z (media anual)")
ax.legend(fontsize=7.5, frameon=False, ncol=3, loc="upper left")
ax.set_title("(b) Anomalías estandarizadas, media por año calendario", fontsize=9.5, loc="left")
ax = axs[2]
ax.plot(anu.index, anu["P"], "o-", color=SERIES["P_L"][3], ms=3, label="P_L")
ax.plot(anu.index, anu["R"], "o-", color=SERIES["Q"][3], ms=3, label="R")
ax.set_ylabel("mm/año")
ax.legend(fontsize=7.5, frameon=False, ncol=2, loc="upper left")
ax = axs[3]
ax.plot(anu.index, anu["C"], "s-", color=GRIS, ms=3, lw=0.9)
ax.set_ylabel("R/P")
ax.set_title("(d) Coeficiente de escorrentía del año hidrológico (etiqueta = año en que termina)", fontsize=9.5, loc="left")
for a_ in axs:
    a_.grid(alpha=0.25, lw=0.5)
ax = axs[2]
e = ESC.set_index("serie_anual")
ax.set_title(f"(c) Años hidrológicos sep–ago completos ({e.loc['P', 'anios_hidrologicos']}). Sen: P {e.loc['P', 'sen_dec']:+.0f}, "
             f"R {e.loc['R', 'sen_dec']:+.0f} mm/año/déc, R/P {e.loc['C', 'sen_dec']:+.3f}/déc, residuo R|P "
             f"{e.loc['resid_R_dado_P', 'sen_dec']:+.0f} mm/año/déc (p {e.loc['resid_R_dado_P', 'p_mk']:.2f})",
             fontsize=9, loc="left")
fig.suptitle("Punto 3.6 · Variabilidad ENSO y respuesta del caudal a la lluvia", fontsize=11)
fig.tight_layout()
fig.savefig(OUT / "fig_causas_3_6.png", dpi=150)
plt.close(fig)

print("\nSalidas en", OUT)
for f in sorted(OUT.iterdir()):
    print("  ", f.name)
