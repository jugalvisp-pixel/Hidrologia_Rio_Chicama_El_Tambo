"""
Relleno de los faltantes de caudal (CAMELS-PE `flow_obs`, mm/d) de la cuenca El Tambo.

Decisión del grupo (2026-10-10): rellenar los faltantes de Q y R aceptando una pérdida de calidad, con un
método explícito y validado (script 28). Se rellena solo dentro del registro observado (1998-01 en adelante).

Paso 1 – Interpolación diaria. Cada hueco de hasta MAX_HUECO_DIAS días con datos a ambos lados se rellena
  interpolando linealmente ln(q) entre el último día anterior y el primero posterior. En escala logarítmica
  la interpolación sigue la forma exponencial de la recesión. Validación (script 28): al borrar L días
  seguidos de meses completos, el error mediano del total mensual es 0.1–1.4 % para L ≤ 7 y 2–5 % para
  L = 10–15 (el escalado anterior, media de los días válidos · n_días, daba 3–4 veces más error).
Paso 2 – Meses con faltantes restantes (huecos > MAX_HUECO_DIAS días):
  - si después del paso 1 hay datos en al menos la mitad del mes, R = media de los días disponibles · n_días
    (escalado). Validación: con 15 días seguidos borrados el error mediano del escalado es ~11 %, la mitad
    del error del modelo mensual, así que conviene aprovechar los días observados;
  - si hay menos de la mitad de los días, el mes se reconstruye con un modelo mensual de regresión en escala logarítmica:
        ln R_m = a_mes + b·ln S_m + c·ln R_(m−1)
    con S = caudal simulado de CAMELS-PE (`flow_sim`, modelo PISCO-ARNOVIC v1.1) y R_(m−1) el mes anterior
    (observado o ya rellenado). Si R_(m−1) no existe se usa ln R_m = a_mes + b·ln S_m. Se ajusta con los
    meses observados completos y se corrige el sesgo de la retransformación (+σ²/2). Validación cruzada
    dejando un año fuera: NSE 0.72, error relativo mediano 23 %. No usa P_L directamente, pero S proviene
    de un modelo forzado con PISCOp, así que los meses reconstruidos no deben usarse para evaluar la
    relación lluvia–caudal (punto 2): se identifican con la bandera R_flag = "modelo_mensual".
"""
import numpy as np
import pandas as pd

MAX_HUECO_DIAS = 15
MAX_FALT_ESCALA = 0.50


def interpolar_log(q, max_hueco=MAX_HUECO_DIAS):
    """Rellena huecos de hasta max_hueco días (con datos a ambos lados) interpolando ln(q)."""
    na = q.isna()
    grupo = (na != na.shift()).cumsum()
    largo = na.groupby(grupo).transform("size")
    lq = np.log(q.clip(lower=1e-4)).interpolate(limit_area="inside")
    relleno = np.exp(lq).where(na & (largo <= max_hueco))
    return q.combine_first(relleno), relleno.notna()


def _X(t, vars_):
    meses = pd.get_dummies(t.index.month).reindex(columns=range(1, 13), fill_value=0).values.astype(float)
    return np.column_stack([meses, *[np.log(t[c].clip(lower=0.01)).values for c in vars_]])


def ajustar_modelo(R, S):
    """Ajusta los dos modelos log-lineales con meses observados completos. Devuelve {nombre: (coef, var)}."""
    t = pd.DataFrame({"R": R, "S": S, "R1": R.shift(1)})
    out = {}
    for nombre, vars_ in (("S_R1", ["S", "R1"]), ("S", ["S"])):
        tt = t.dropna(subset=["R", *vars_])
        y = np.log(tt["R"].clip(lower=0.01)).values
        b, *_ = np.linalg.lstsq(_X(tt, vars_), y, rcond=None)
        out[nombre] = (b, vars_, np.var(y - _X(tt, vars_) @ b))
    return out


def predecir(modelos, S_m, R1, fecha):
    """Predice R de un mes (fecha) con S_m y el R del mes anterior (puede ser NaN)."""
    nombre = "S_R1" if np.isfinite(R1) else "S"
    b, vars_, var = modelos[nombre]
    t = pd.DataFrame({"S": [S_m], "R1": [R1]}, index=pd.DatetimeIndex([fecha]))
    return float(np.exp(_X(t, vars_) @ b + var / 2)[0]), nombre
