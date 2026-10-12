"""
Visor HTML de los productos de la Tarea 1 (cuenca El Tambo, río Chicama).

Genera un único archivo autocontenido, `visor_resultados.html`, en la raíz del proyecto: las figuras
van incrustadas (JPEG comprimido en base64) y no se usa ninguna librería ni fuente externa, de modo
que se abre con doble clic en cualquier navegador, sin conexión y sin instalar nada. Si además se
clona el repositorio completo, cada figura enlaza a su PNG original en resolución completa.

Las cifras de los textos se leen de las tablas de resultados (no están escritas a mano), así que al
volver a ejecutar los scripts 06–21 y luego este, la página queda actualizada.
Requiere haber ejecutado los scripts 06–21.
"""
from pathlib import Path
from datetime import date
import base64
import html
import io
import json

import numpy as np
import pandas as pd
import xarray as xr
from PIL import Image

from periodo_comun import periodo_comun

PROY = Path(__file__).resolve().parents[2]                 # raíz del proyecto (Hidrologia/)
R = PROY / "Tarea 1" / "El_Tambo" / "resultados"
SM, CL, CC, CF, DL = R / "series_mensuales", R / "climatologia", R / "control_calidad", R / "contexto_fisico", R / "delimitacion"
SALIDA = PROY / "visor_resultados.html"
ANCHO_MAX, CALIDAD = 1800, 82                               # compresión de las figuras incrustadas

INI, FIN = periodo_comun()
PER = f"{INI:%Y-%m} a {FIN:%Y-%m}"
N_MESES = len(pd.date_range(INI, FIN, freq="MS"))


# ------------------------------------------------------------------ cifras leídas de los resultados
def leer(p, **kw):
    return pd.read_csv(p, **kw)


cam = leer(SM / "series_mensuales_camels_pe.csv", index_col="fecha", parse_dates=True)
imerg = leer(SM / "imerg_cuenca_mensual.csv", index_col="fecha", parse_dates=True)["P_I_mm_mes"]
pp = pd.concat([cam["P_L_mm_mes"], imerg], axis=1).loc[INI:FIN].dropna()
d = pp["P_I_mm_mes"] - pp["P_L_mm_mes"]
IM = dict(PL=pp["P_L_mm_mes"].mean(), PI=pp["P_I_mm_mes"].mean(), sesgo=d.mean(), mae=d.abs().mean(),
          rmse=np.sqrt((d ** 2).mean()), r=pp.corr().iloc[0, 1], ratio=pp["P_I_mm_mes"].sum() / pp["P_L_mm_mes"].sum())
pesos = leer(PROY / "Tarea 1" / "El_Tambo" / "datos" / "imerg" / "pesos_imerg_cuenca.csv")
PW = dict(n=len(pesos), comp=int((pesos.fraccion_celda > 0.99).sum()),
          borde=100 * pesos.loc[pesos.fraccion_celda <= 0.99, "area_en_cuenca_km2"].sum() / pesos.area_en_cuenca_km2.sum(),
          eq=pesos.fraccion_celda.sum())
est = leer(SM / "estadisticos_1_3_periodo_comun.csv", index_col="variable")
cl = leer(CL / "climatologia_1_5a.csv")
ind = leer(CL / "indicadores_estacionalidad_1_5b.csv")
ip = ind[ind.muestra.str.startswith("periodo")].set_index("variable")
val = leer(SM / "validacion_extension_T.csv", index_col=0)["valor"]
e5 = leer(SM / "era5land_t2m_cuenca_mensual.csv", index_col="fecha", parse_dates=True)["T_era5land_C"]
tt = pd.concat([cam["T_C"], e5], axis=1).loc[INI:FIN].dropna()
TB = dict(sesgo=(tt["T_era5land_C"] - tt["T_C"]).mean(), r=tt.corr().iloc[0, 1], n=len(tt))
acal = leer(CL / "anomalias_anio_calendario_1_5b.csv").sort_values("P_L_anom_total_mm")
humedos = ", ".join(str(int(a)) for a in acal.tail(4)["anio"][::-1])
secos = ", ".join(str(int(a)) for a in acal.head(4)["anio"])
reg = leer(CC / "registro_anomalias_1_4.csv")
RE = PROY / "Tarea 1" / "El_Tambo" / "resultados" / "relaciones"
MO = PROY / "Tarea 1" / "El_Tambo" / "resultados" / "modelos"
_m21 = leer(RE / "metricas_2_1.csv")
_m21 = _m21[(_m21.serie == "original") & (_m21.muestra == "todos los meses")].set_index("par")
R21 = dict(a=_m21.loc["(a) P_I vs P_L", "pearson_r"], b=_m21.loc["(b) P_L vs Q̄", "spearman_rho"],
           c=_m21.loc["(c) P_I vs Q̄", "spearman_rho"])
E21 = leer(RE / "metricas_IMERG_por_estacion_2_1.csv", index_col=0).loc["todo el año"]
_rz = leer(RE / "rezagos_2_1.csv").query("lluvia == 'P_L'")
RZ = dict(k0=_rz.query("serie == 'original' and tipo == 'rezago P(t−k)' and k == 0").spearman_rho.iloc[0],
          k1=_rz.query("serie == 'original' and tipo == 'rezago P(t−k)' and k == 1").spearman_rho.iloc[0],
          acu=_rz.query("serie == 'original' and tipo == 'acumulada P(t−k..t)'").spearman_rho.max(),
          an1=_rz.query("serie == 'anomalías' and tipo == 'rezago P(t−k)' and k == 1").spearman_rho.iloc[0])
# Inversión térmica: radiosondeos de Trujillo (scripts 34a–34b) y transporte de humedad sobre la cuenca (script 24)
IV = R / "inversion"
_snd = leer(IV / "inversion_por_sondeo.csv", parse_dates=["fecha"])


def _temporada(t):
    s = _snd[_snd.temporada == t]
    si = s[s.inversion]
    return dict(n=len(s), frec=100 * s.inversion.mean(), base=si.base_m.median(), tope=si.tope_m.median(),
                dT=si.intensidad_C.median())


SND = dict(seca=_temporada("may–sep (seca)"), hum=_temporada("dic–abr (húmeda)"), n=len(_snd),
           frec=100 * _snd.inversion.mean(), ini=_snd.fecha.min(), fin=_snd.fecha.max())
_ct = leer(IV / "contraste_sondeos_era5.csv")
_ct["clave"] = (_ct.temporada.str[:3].map({"may": "seca", "dic": "hum"}) + "_"
                + np.select([_ct.fuente.str.startswith("Trujillo"), _ct.fuente.str.contains("periodo común")],
                            ["snd", "e5c"], "e5s"))
CT = _ct.set_index("clave")["estabilidad_T1500_menos_T500_C"].to_dict()     # T(1 500 m) − T(500 m), script 34c
_cl = xr.open_dataset(PROY / "Tarea 1" / "El_Tambo" / "datos" / "era5" / "era5_viento_ivt_climatologia_estacional.nc")
_cu = _cl.sel(lat=-7.574, lon=-78.707, method="nearest")      # celda de 1° que contiene el cierre
IVT = {e: (float(_cu["ivt"].sel(estacion=e)), float(_cu["tcwv"].sel(estacion=e))) for e in ("dic-abr", "may-sep")}
_fl = cam.loc[INI:FIN]
RL = dict(interp=int((_fl.R_flag == "interpolado_diario").sum()), esc=int((_fl.R_flag == "escalado").sum()),
          mod=int((_fl.R_flag == "modelo_mensual").sum()),
          mes_mod=", ".join(f"{f:%Y-%m}" for f in _fl.index[_fl.R_flag == "modelo_mensual"]),
          obs=int((_fl.R_flag == "completo").sum()), excl=int(_fl.R_flag_obs.isin(["excluido", "sin_dato"]).sum()))
TE = PROY / "Tarea 1" / "El_Tambo" / "resultados" / "tendencias"
_td = leer(TE / "tendencias_todos_meses_3_4.csv")
_ta = _td[(_td.representacion == "anomalia") & _td.variable.isin(["P_L", "P_I", "Q", "T_PISCOt", "T_ERA5L"])]
t3 = _ta[["variable", "periodo", "inicio", "fin", "n", "metodo", "pend", "lo", "hi", "p", "unidad"]].copy()
t3["periodo"] = t3["periodo"].replace({"comun": "común"})
t3 = t3.round({"pend": 3, "lo": 3, "hi": 3, "p": 3})
t3.columns = ["Variable", "Periodo", "Inicio", "Fin", "n", "Método", "Pendiente", "IC inf.", "IC sup.", "p", "Unidad"]
TD = _ta.set_index(["variable", "periodo", "metodo"])


def pend3(v, per="completo", met="MK Hamed-Rao + Sen"):
    r = TD.loc[(v, per, met)]
    return f"{r.pend:+.2f} [{r.lo:+.2f}, {r.hi:+.2f}] (p = {r.p:.3f})"


_enso = leer(TE / "enso_3_6.csv").set_index("variable")
_esc = leer(TE / "escorrentia_anual_3_6.csv").set_index("serie_anual")
_seg = leer(TE / "segmentada_pettitt_3_5.csv").set_index("variable")
_mes = leer(TE / "tendencias_por_mes_3_4.csv")
_fdr = _mes[(_mes.representacion == "original") & (_mes.periodo == "completo") & (_mes.q_fdr < 0.05)]
FDR = {v: ", ".join(sorted({["Ene", "Feb", "Mar", "Abr", "May", "Jun", "Jul", "Ago", "Sep", "Oct", "Nov", "Dic"][m - 1]
                            for m in g.mes}, key=lambda s: ["Ene", "Feb", "Mar", "Abr", "May", "Jun", "Jul", "Ago", "Sep",
                                                             "Oct", "Nov", "Dic"].index(s)))
       for v, g in _fdr.groupby("variable")}
ES = PROY / "Tarea 1" / "El_Tambo" / "resultados" / "espectros"
_bd = leer(ES / "bandas_varianza_4_2.csv").set_index(["variable", "periodo", "transformacion", "banda"])["fraccion"]


def banda(v, tr, b, per="completo"):
    return 100 * _bd.loc[(v, per, tr, b)]


_me = leer(ES / "metadatos_espectros_4_1.csv")
RHO = _me[(_me.periodo == "completo") & (_me.transformacion == "anomalia")].set_index("variable")["rho1_AR1"]
_pk = leer(ES / "picos_4_2.csv")
_pk = _pk[(_pk.periodo == "completo") & (_pk.transformacion == "anomalia")]
N_GLOBAL = int(_pk["supera_95_global"].astype(str).eq("True").sum())
# mayor pico interanual (≥ 2 años, ≥ 3 ciclos) de las anomalías del registro completo: (periodo en años, % varianza)
_pi = _pk[(_pk.estimador == "periodograma") & (_pk.periodo_meses >= 24) & (_pk.k_ciclos_en_registro >= 3)]
PK = {v: (g.loc[g.fraccion_varianza.idxmax(), "periodo_meses"] / 12, 100 * g.fraccion_varianza.max())
      for v, g in _pi.groupby("variable")}
OND = leer(ES / "onditas_resumen_4_3.csv").set_index("serie")                 # script 38
t4 = _pk[_pk.groupby("variable").cumcount() < 3][["variable", "N_meses", "periodo_meses", "periodo_min_meses",
                                                   "periodo_max_meses", "k_ciclos_en_registro", "fraccion_varianza",
                                                   "cociente_fondo_AR1", "supera_95_puntual", "supera_95_global"]].copy()
t4["fraccion_varianza"] = (t4["fraccion_varianza"] * 100).round(1)
t4 = t4.round({"periodo_meses": 1, "periodo_min_meses": 0, "periodo_max_meses": 0, "cociente_fondo_AR1": 1})
t4.columns = ["Variable", "N", "Periodo (meses)", "Rango mín.", "Rango máx.", "Ciclos en el registro", "% varianza",
              "Cociente con AR(1)", "> 95 % puntual", "> 95 % global"]
_ls = leer(ES / "vacios_caudal_lombscargle_4_1.csv").set_index("banda")
pico = leer(CL / "mes_pico_por_anio_1_5b.csv")
pct_mar = (pico["mes_max_P_L"] == "Mar").mean() * 100
pct_fa = pico["mes_max_P_L"].isin(["Feb", "Mar", "Abr"]).mean() * 100
pct_qabr = (pico["mes_max_Q"].dropna() == "Abr").mean() * 100


def cmes(var, mes):
    return cl[(cl.variable == var) & (cl.mes == mes)].iloc[0]


def lag_centroide(p):
    return (ip.loc["R", "centroide_mes"] - ip.loc[p, "centroide_mes"]) % 12


# ------------------------------------------------------------------ utilidades HTML
def img64(p):
    im = Image.open(p).convert("RGB")
    if im.width > ANCHO_MAX:
        im = im.resize((ANCHO_MAX, round(im.height * ANCHO_MAX / im.width)), Image.LANCZOS)
    buf = io.BytesIO()
    im.save(buf, "JPEG", quality=CALIDAD, optimize=True, progressive=True)
    return "data:image/jpeg;base64," + base64.b64encode(buf.getvalue()).decode()


def rel(p):
    return Path(p).relative_to(PROY).as_posix().replace(" ", "%20")


# ------------------------------------------------------------------ datos de las gráficas interactivas (Plotly)
# Construcción tomada del 36_pagina_web.py del visor de la compañera; las gráficas están en visor_graficas.js.
MO = PROY / "Tarea 1" / "El_Tambo" / "resultados" / "modelos"
MESES = ["Ene", "Feb", "Mar", "Abr", "May", "Jun", "Jul", "Ago", "Sep", "Oct", "Nov", "Dic"]


def lista(s, nd=2):
    return [None if pd.isna(v) else round(float(v), nd) for v in s]


_im_full = leer(SM / "imerg_cuenca_mensual.csv", index_col="fecha", parse_dates=True)["P_I_mm_mes"]
_chirps = leer(SM / "chirps_cuenca_mensual.csv", index_col="fecha", parse_dates=True)["P_chirps_mm_mes"]
_piscop = leer(SM / "piscop_v3_cuenca_mensual.csv", index_col="fecha", parse_dates=True)["P_piscop_v3_mm_mes"]
_text = leer(SM / "temperatura_piscot_extendida_era5land.csv", index_col="fecha", parse_dates=True)
idx = pd.date_range(INI, FIN, freq="MS")
_anom = {}
for k, s in {"P_L": cam["P_L_mm_mes"], "P_I": _im_full, "R": cam["R_mm_mes"],
             "T": _text["T_piscot_C"].fillna(_text["T_estimada_era5land_C"])}.items():
    s = s.reindex(idx)
    a = s - s.index.month.map(s.groupby(s.index.month).mean()).values         # igual que el script 21
    _anom[k] = {"m": lista(a), "mm12": lista(a.rolling(12, center=True, min_periods=10).mean())}
_anios_comp = [y for y in range(INI.year, FIN.year + 1) if pp.index.year.tolist().count(y) == 12]
_anual = pd.DataFrame({"P_L": cam["P_L_mm_mes"], "P_I": _im_full}).groupby(lambda f: f.year).sum().loc[_anios_comp]
_acal = acal.sort_values("anio")
DATOS = {
    "fechas": [f"{f:%Y-%m}" for f in idx],
    "P_L": lista(cam["P_L_mm_mes"].reindex(idx)), "P_I": lista(_im_full.reindex(idx)),
    "chirps": lista(_chirps.reindex(idx)), "piscop": lista(_piscop.reindex(idx)),
    "Q": lista(cam["Q_m3_s"].reindex(idx), 3), "R": lista(cam["R_mm_mes"].reindex(idx)),
    "R_flag_obs": cam["R_flag_obs"].reindex(idx).fillna("sin_dato").tolist(),
    "T_piscot": lista(_text["T_piscot_C"].reindex(idx)), "T_est": lista(_text["T_estimada_era5land_C"].reindex(idx)),
    "T_era5": lista(e5.reindex(idx)),
    "anual": {"anios": [int(a) for a in _anual.index], "P_L": lista(_anual["P_L"], 1), "P_I": lista(_anual["P_I"], 1)},
    "clim": {v: {"unid": g.unidades.iloc[0], **{c: lista(g.sort_values("mes")[c], 2)
                                                for c in ("media", "mediana", "Q1", "Q3", "P10", "P90")}}
             for v, g in cl.groupby("variable")},
    "anom": _anom,
    "anom_anual": {"anios": _acal["anio"].astype(int).tolist(),
                   "P_L": lista(_acal["P_L_anom_total_mm"], 1), "P_I": lista(_acal["P_I_anom_total_mm"], 1),
                   "R": lista(_acal["R_anom_total_mm"], 1), "T": lista(_acal["T_anom_media_C"], 2)},
}
# diarios (control de calidad): caudal sin rellenar, días interpolados, precipitación y temperatura
_dia = leer(PROY / "Tarea 1" / "El_Tambo" / "datos" / "camels_pe" / "PE_47E0A2A8.csv", index_col="date",
            parse_dates=True).loc[INI:FIN + pd.offsets.MonthEnd(0)]
_qint = leer(SM / "caudal_diario_interpolado.csv", index_col="fecha", parse_dates=True).reindex(_dia.index)
A_KM2 = 2171.585                                                  # área CAMELS-PE, igual que el script 09
_cond = leer(CC / "disponibilidad_condicion_mes.csv", index_col="fecha", parse_dates=True).loc[INI:FIN]
_anios = list(range(INI.year, FIN.year + 1))


def _matriz(valido):
    """Fracción de días válidos por año × mes."""
    g = valido.groupby([valido.index.year, valido.index.month]).mean()
    return [[None if (y, m) not in g.index else round(float(g[(y, m)]), 3) for m in range(1, 13)] for y in _anios]


def _matriz_cond(col):
    return [[_cond[col].get(pd.Timestamp(y, m, 1), None) for m in range(1, 13)] for y in _anios]


DATOS["diario"] = {"ini": f"{_dia.index[0]:%Y-%m-%d}", "n": len(_dia),
                   "P_L": lista(_dia["prec"], 1), "T": lista(_dia["tmean"], 1),
                   "Q": lista(_dia["flow_obs"] * A_KM2 / 86.4, 2),
                   "Q_int": lista((_qint["q_mm_d"] * A_KM2 / 86.4).where(_qint["interpolado"] == 1), 2)}
DATOS["disp"] = {"anios": _anios,
                 "P_L": {"f": _matriz(_dia["prec"].notna()), "c": _matriz_cond("P_L")},
                 "Q": {"f": _matriz(_dia["flow_obs"].notna()), "c": _matriz_cond("Q")},
                 "T": {"f": _matriz(_dia["tmean"].notna()), "c": _matriz_cond("T")},
                 "P_I": {"f": _matriz(_im_full.reindex(_dia.index, method="ffill").notna()), "c": _matriz_cond("P_I")}}
# evaluación fuera de muestra (2.3)
_pr = leer(MO / "predicciones_2_3.csv", parse_dates=["fecha"])
_pr = _pr[_pr.conjunto == "evaluacion"]
_obs2 = {"P_L ← P_I": cam["P_L_mm_mes"], "Q̄ ← P_L": cam["Q_m3_s"]}
_mod2 = {"P_L ← P_I": [("ref", "Referencia: climatología"), ("lin", "Regresión lineal"),
                       ("mes", "Elegido: factor por mes calendario")],
         "Q̄ ← P_L": [("ref", "Referencia: climatología"), ("lin0", "Regresión lineal"),
                     ("logmK", "Elegido: ln Q̄ con término por mes y lluvia acumulada")]}
DATOS["modelos"] = {"fechas": [f"{f:%Y-%m}" for f in idx]}
for obj, mods in _mod2.items():
    d2 = {"obs": lista(_obs2[obj].reindex(idx), 2), "mods": []}
    for cod, nom in mods:
        s = _pr[(_pr.objetivo == obj) & (_pr.modelo == cod)].groupby("fecha")["estimado"].mean().reindex(idx)
        d2["mods"].append({"nombre": nom, "y": lista(s, 2)})
    DATOS["modelos"]["PL" if obj.startswith("P_L") else "Q"] = d2
# anomalías con la tendencia de Sen (3.4), registro completo
_clr = leer(TE / "climatologia_referencia_3_2.csv")
DATOS["tend"] = {}
for v, s in {"P_L": cam["P_L_mm_mes"], "P_I": _im_full, "Q": cam["Q_m3_s"], "T_PISCOt": cam["T_C"],
             "T_ERA5L": e5}.items():
    r = TD.loc[(v, "completo", "MK Hamed-Rao + Sen")]
    s = s.loc[pd.Timestamp(r.inicio + "-01"):pd.Timestamp(r.fin + "-01")]
    mu = _clr[_clr.variable == v].set_index("mes")["media"]
    a = (s - s.index.month.map(mu).values).dropna()
    t = a.index.year + (a.index.month - 0.5) / 12                      # años decimales, como el script 32
    b = r.pend / 10                                                    # pendiente por año
    c = float(np.median(a.values - b * t))                             # ordenada de Sen
    DATOS["tend"][v] = {"fechas": [f"{f:%Y-%m}" for f in a.index], "a": lista(a, 3),
                        "mm12": lista(a.rolling(12, center=True, min_periods=10).mean(), 3),
                        "sen": [round(c + b * t[0], 3), round(c + b * t[-1], 3)],
                        "sen_x": [f"{a.index[0]:%Y-%m}", f"{a.index[-1]:%Y-%m}"],
                        "txt": f"Sen {r.pend:+.2f} {r.unidad} [{r.lo:+.2f}, {r.hi:+.2f}], p = {r.p:.3f}",
                        "unid": r.unidad.split("/")[0] if v.startswith("T") else r.unidad.rsplit("/", 1)[0]}

# panel lateral del control de calidad (1.4): resumen, tratamiento y registro de cada variable
_rd = leer(CC / "resumen_disponibilidad_1_4.csv")
_rd = _rd[_rd.periodo.str.startswith("periodo")].set_index("variable")
_VQC = {"P_L": ("P<sub>L</sub> · precipitación", r"P_L|Todas"), "Q": ("Q̄ · caudal", r"\bQ\b|Todas"),
        "T": ("T · temperatura", r"\bT\b|Todas"), "P_I": ("P<sub>I</sub> · IMERG", r"P_I")}

# «Cómo se construyó» de cada figura del visor (obligatorio: si falta una entrada, el script se detiene).
COMO = {
    "mapa_cuenca_el_tambo.png":
        "Delimitación con el DEM Copernicus GLO-30 (30 m) reproyectado a UTM 17S, acondicionado (relleno de "
        "de pozos y depresiones, resolución de planos) y direcciones de flujo D8 con pysheds; punto de cierre reubicado, "
        "dentro de 1.5 km de la estación, sobre una celda cuya área drenada difiere menos de 15 % de la de CAMELS-PE. Red de drenaje = celdas con área aportante > 5 km², ordenada por Strahler. "
        "Relieve sombreado del mismo DEM. Estaciones de CAMELS-PE v1.0.1 y del mapa web de SENAMHI. Recuadros con "
        "límites de Natural Earth 1:50m.",
    "panel_contexto_fisico.png":
        "Capas llevadas a la malla del DEM GLO-30 geográfico (~30 m). (a) Orientación calculada del DEM. (b) ESA "
        "WorldCover 2021 (10 m) remuestreado por moda. (c) Geología 1:100 000 de INGEMMET (GEOCATMIN), unidades "
        "agrupadas en grupos hidrogeológicos según su litología. (d) CHIRPS v2.0 (0.05°): media anual del periodo "
        "común, interpolada de forma bilineal solo para dibujarla continua, con isoyetas cada 100 mm.",
    "fig_sondeo_tipico_inversion.png":
        "Sondeo observado de Trujillo (IGRA2 84516, NOAA NCEI; 12 UTC). Las alturas de los niveles significativos "
        "se calcularon con la ecuación hipsométrica desde la superficie y se ajustaron a la altura geopotencial de "
        "los niveles estándar; HR desde la depresión del punto de rocío (Magnus–Bolton). La inversión es el primer "
        "aumento continuo de T con la altura (≥ 1 °C, por debajo de 3 500 m); la capa de nubes, el tramo con "
        "HR ≥ 95 %. La parcela sube por la adiabática seca hasta el NCA y luego por la pseudoadiabática. Se eligió "
        "el sondeo de may–sep más cercano a las medianas de base, tope e intensidad de esa estación.",
    "fig_temporadas_humedad_inversion.png":
        "Izquierda: IVT de ERA5 mensual (0.5°, promediado a 1°), media de cada temporada en el periodo común; "
        "color = magnitud, flechas = dirección (una cada 3°). Derecha: radiosondeos de Trujillo interpolados a una "
        "malla vertical de 25 m sin extrapolar; mediana y rango intercuartílico de la temporada, adiabática seca de "
        "referencia desde la superficie y capa de inversión entre las medianas de base y tope de los sondeos con "
        "inversión. Mismas escalas en las dos filas.",
    "fig_contraste_1_trujillo_vs_era5.png":
        "Sondeos de Trujillo (IGRA2 84516, 12 UTC) de 2024-10 a 2026-09 interpolados a una malla de 25 m: media "
        "por nivel y rango intercuartílico de cada temporada. ERA5 mensual por niveles de presión, promedio de las "
        "12 UTC de los mismos meses, en la celda de 0.25° más cercana a Trujillo (−8.0°, −79.0°), con cada nivel en "
        "su altura geopotencial media. Se usa la media de los sondeos (no la mediana) porque ERA5 es un promedio. "
        "(c) ERA5 interpolado a la malla de 25 m menos la media de los sondeos.",
    "fig_contraste_2_era5_periodos.png":
        "ERA5 mensual por niveles de presión, promedio de las 12 UTC, misma celda que la comparación 1: media de "
        "cada temporada en el periodo común frente a la media en los meses con sondeos (2024-10 a 2026-09). "
        "(c) Diferencia entre ambos periodos en la malla de 25 m.",
    "fig_series_mensuales_1_1.png":
        "Desde los registros diarios de CAMELS-PE: P_L = suma diaria (solo meses completos); Q̄ = media diaria en "
        "m³/s; R = suma del caudal diario en lámina (mm); T = media diaria de PISCOt. Caudal rellenado: interpolación "
        "logarítmica de huecos ≤ 15 días, escalado si hay ≥ 50 % de días y modelo mensual en el resto. T desde 2021 "
        "= ERA5-Land + sesgo de cada mes (script 17).",
    "fig_T_camels_vs_era5land.png":
        "ERA5-Land T2m mensual (0.1°, CDS) promediado sobre la cuenca con pesos por área de intersección y "
        "convertido de K a °C, frente a `tmean` de CAMELS-PE (PISCOt v1.2). Comparación en los meses comunes; "
        "relieve del modelo calculado del geopotencial de ERA5-Land.",
    "mapa_celdas_imerg_pesos.png":
        "Celdas de IMERG (0.1°) intersecadas con el polígono de la cuenca en UTM 17S; peso = área dentro de la "
        "cuenca / área de la cuenca. Desnivel = percentil 95 − percentil 5 de la elevación GLO-30 dentro de cada celda.",
    "fig_P_L_vs_IMERG_mensual.png":
        "IMERG V07 Final mensual (tasa en mm/h × 24 × días del mes) promediado con los pesos por área, frente a "
        "P_L de CAMELS-PE, mes a mes en el periodo común.",
    f"fig_precipitacion_camels_{INI:%Y}_{FIN:%Y}.png":
        "P_L de CAMELS-PE: (a) acumulado mensual y media móvil centrada de 12 meses; (b) totales de los años con "
        "los 12 meses en el periodo común.",
    f"fig_precipitacion_imerg_{INI:%Y}_{FIN:%Y}.png":
        "IMERG con el mismo procedimiento y los mismos límites de ejes que la figura de P_L, para compararlas.",
    "fig_histogramas_1_3.png":
        "Valores mensuales del periodo común. Ancho de clase por Freedman–Diaconis (2·IQR·n^−1/3, redondeado); "
        "P_L y P_I con los mismos meses y límites de clase; eje en % de meses.",
    "fig_cajas_1_3.png":
        "Cuartiles con interpolación lineal (tipo 7), mediana, media y bigotes de Tukey hasta el dato más extremo "
        "dentro de 1.5·IQR; los valores fuera de los bigotes se muestran.",
    "fig_disponibilidad_1_4.png":
        "Días válidos de cada mes en los registros diarios de CAMELS-PE, clasificados con el mismo criterio de "
        "completitud del análisis (antes del relleno); IMERG y ERA5-Land son mensuales.",
    "fig_validacion_relleno_caudal.png":
        "(a) Se borraron huecos sintéticos de distinta longitud (600 por longitud) en meses completos y se compararon "
        "la interpolación logarítmica y el escalado con el valor real. (b) Cinco modelos mensuales evaluados con "
        "validación cruzada dejando un año fuera.",
    "fig_climatologia_1_5a.png":
        "Para cada mes calendario del periodo común: media, mediana, cuartiles y percentiles 10 y 90 (tipo 7) de "
        "los años válidos.",
    "mapas_precipitacion_mensual_chirps.png":
        "CHIRPS v2.0 (0.05°): media de cada mes calendario en el periodo común, celda por celda, interpolada de "
        "forma bilineal a ~280 m solo para dibujarla; escala discreta con clases más finas en valores bajos.",
    "fig_anio_mes_1_5a.png":
        "Matriz año × mes de cada variable en el periodo común; caudal en escala logarítmica; trama en los meses "
        "de temperatura estimados con ERA5-Land.",
    "fig_ciclo_1_5b.png":
        "Indicadores del ciclo de 12 medias: SI (Walsh & Lawler), PCI (Oliver), centroide circular de Markham, "
        "picos y armónicos; mes del máximo de cada año hidrológico sep–ago; subperiodos 1998–2011 y 2012–fin; "
        "CV y asimetría por mes.",
    "fig_anomalias_1_5b.png":
        "Anomalía = valor del mes − media de su mes calendario en el periodo común, con media móvil de 12 meses.",
    "fig_sintesis_ciclo_1_5c.png":
        "Ciclo de P_L frente a la evapotranspiración potencial (Penman–Monteith de CAMELS-PE hasta 2016 y "
        "Hargreaves–Samani con PISCOt); clases FAO (P vs PET), meses secos de Gaussen y aguas altas. Figura "
        "tentativa según el grupo.",
    "fig_dispersion_2_1.png":
        "Pares de meses del periodo común coloreados por estación; Pearson y Spearman; caudal en escala logarítmica; "
        "línea 1:1 en IMERG frente a P_L y línea R = P en lámina.",
    "fig_rezagos_2_1.png":
        "Correlación de Spearman del caudal del mes t con la lluvia de t − k (k = 0 a 6) y con la lluvia acumulada "
        "de t − k a t, en series originales y en anomalías.",
    "fig_evaluacion_resumen_2_3.png":
        "Habilidad = 1 − MAE / MAE de la referencia, fuera del periodo de ajuste, en dos particiones por bloques de "
        "años. Modelo y rezago elegidos con validación cruzada dejando un año fuera dentro del ajuste.",
    "fig_modelo_PL_2_3.png":
        "Modelos de P_L a partir de IMERG ajustados en 1998–2015 y evaluados en 2016–2024; residuos frente al tiempo, "
        "al valor estimado y al mes.",
    "fig_modelo_Q_2_3.png":
        "Modelos de caudal a partir de P_L ajustados en 1998–2015 y evaluados en 2016–2024; los logarítmicos usan la "
        "corrección de Duan al volver a m³/s.",
    "fig_series_tendencias_3_4.png":
        "Registro completo de cada variable. Anomalías con una climatología fija 1998–2020. OLS con errores HAC "
        "(Newey–West, 12 rezagos), Sen con Mann–Kendall (Hamed–Rao) y LOWESS de 10 años con banda por remuestreo de "
        "residuos en bloques de 24 meses; en la serie original, modelo aditivo de efecto mensual + LOWESS.",
    "fig_pendientes_mensuales_3_5.png":
        "Doce subseries (un mes calendario a través de los años). OLS con HAC e intervalo de Sen–Gilbert con varianza "
        "de Hamed–Rao; control FDR de Benjamini–Hochberg sobre 12 meses × 5 variables.",
    "fig_robustez_3_5.png":
        "Pendientes de las anomalías recalculadas con otros métodos, periodos (recorte de 5 años en cada extremo, "
        "sin 1998, 2017 y 2023) y fuentes (PISCOp v3, CHIRPS, caudal solo observado).",
    "fig_ventanas_inicio_fin_3_5.png":
        "Pendiente de Sen y prueba de Mann–Kendall de las medias anuales de anomalías para cada par de años inicial y "
        "final con 15 años o más.",
    "fig_causas_3_6.png":
        "Índice Niño 1+2 de ERSST v5 (CPC/NOAA); medias anuales de anomalías estandarizadas; balance por año "
        "hidrológico sep–ago con P_L y R, y regresión de R sobre P para el residuo.",
    "fig_espectros_4_2.png":
        "Periodograma de la FFT sin relleno con ceros: f<sub>k</sub> = k/N, Δf = 1/N ciclos/mes, densidad unilateral "
        "normalizada para que Σ S·Δf sea la varianza; Welch con Hann de 120 meses y 50 % de traslape; fondo AR(1) con "
        "umbrales χ² al 95 % puntual y global (Bonferroni). Eje inferior en periodo (años) = 1/(12 f), creciente hacia "
        "la derecha; eje superior en frecuencia (ciclos/mes). Las dos columnas de cada variable tienen la misma escala "
        "vertical.",
    "fig_espectros_porcentaje_4_2.png":
        "Los mismos periodogramas expresados como % de la varianza en cada frecuencia de Fourier: % = S<sub>k</sub>·Δf / "
        "varianza, con el eje de periodo lineal en años (hasta 20; lo que queda por encima se indica en cada panel). "
        "Cada pico rotulado lleva el rango de periodos que cabe en ±½ bin de frecuencia y el número de ciclos k en el "
        "registro; solo se rotulan picos con k ≥ 3. * = supera el 95 % puntual del fondo AR(1); ** = supera también el "
        "umbral global. Idea tomada del visor paralelo del grupo.",
    "fig_onditas_4_3.png":
        "Transformada continua con la ondita de Morlet (ω<sub>0</sub> = 6; Torrence y Compo 1998) de las anomalías "
        "estandarizadas (climatología 1998–2020), con la potencia corregida por sesgo (Liu et al. 2007: se divide por la "
        "escala) y expresada en % de la varianza por octava de periodo. Contorno: supera el 95 % de un ruido rojo AR(1); "
        "rayado: cono de influencia (bordes); a la derecha, el espectro global sobre el periodograma de Fourier en la "
        "misma unidad. El relleno con ceros es interno a la transformada; los datos no se rellenan.",
    "fig_onditas_banda_4_3.png":
        "(a) Varianza promediada en la banda de 2–7 años dividida por su umbral del 95 % (Torrence y Compo 1998, ec. "
        "24–28). (b) Componente de 2–7 años reconstruida (ec. 29) de la cuenca y de Niño 1+2 (ERSST v5, CPC/NOAA) y su "
        "correlación; es descriptiva, porque una serie filtrada a 2–7 años tiene muy pocos grados de libertad.",
    "fig_bandas_varianza_4_2.png":
        "Suma de las ordenadas del periodograma dentro de cada banda de periodos, dividida por la suma total.",
    "fig_comparacion_espectros_4_2.png":
        "Espectros de Welch divididos por la varianza de cada serie, en el periodo común; cociente de los espectros "
        "normalizados de caudal y P_L.",
    "fig_estabilidad_4_2.png":
        "Fracciones por banda y periodo del pico de 2–8 años, repetidos con ventana de Hann, sin tendencia, con "
        "anomalías recortadas en los percentiles 1 y 99, con ln Q y por mitades del registro.",
    "fig_vacios_caudal_4_1.png":
        "Periodograma de Lomb–Scargle con solo los meses de caudal observados completos, frente a la FFT de la serie "
        "rellenada; ambos normalizados por su potencia total.",
    # gráficas interactivas (clave "int:<nombre>", ver visor_graficas.js)
    "int:serie11":
        "P<sub>L</sub> = suma de la precipitación diaria (PISCOp v2.1 promediado en la cuenca por CAMELS-PE), solo meses "
        "completos. Q̄ = media del caudal diario observado (SENAMHI), convertido de mm/d a m³/s con Q̄ = q̄·A/86.4 y "
        "A = 2 171.6 km². R = suma de q diario en mm. T = media de la temperatura diaria (PISCOt v1.2) si faltan ≤ 10 % de "
        "los días. El relleno de caudal sigue tres pasos: interpolación de ln q en huecos ≤ 15 días, escalado de la media "
        "en meses con ≥ 50 % de días y modelo mensual con el caudal simulado en el resto (validación en el punto 1.4).",
    "int:temp":
        "ERA5-Land T2m mensual (0.1°) descargado del CDS y promediado sobre la cuenca ponderando por área (30 celdas). El "
        "sesgo de cada mes calendario frente a PISCOt se calculó en 1998–2010 y se evaluó fuera de muestra en 2011–2020; "
        "con él se estima T desde 2021: T<sub>est</sub> = ERA5-Land + sesgo del mes. (b) Media de cada mes calendario en "
        "1998–2020, solo con los meses en que existen las dos fuentes.",
    "int:precip":
        "IMERG mensual: la tasa media (mm/h) se convierte a acumulado con mm/mes = tasa · 24 · días del mes, y se promedia "
        "sobre la cuenca ponderando cada celda por su área de intersección con el polígono (UTM 17S). CHIRPS (0.05°) y "
        "PISCOp v3.0 (0.1°) se promediaron igual. P<sub>L</sub> es el promedio de cuenca que entrega CAMELS-PE.",
    "int:anual":
        "Suma de los 12 meses de cada año calendario del periodo común (años completos).",
    "int:qc":
        "Días válidos contados en el registro diario de CAMELS-PE (P<sub>L</sub>, caudal y T) y meses de IMERG. Un mes es "
        "completo si faltan ≤ 10 % de sus días, criterio aplicado igual a todas las variables; no se presentan sumas "
        "parciales como acumulados mensuales completos. Caudal convertido de mm/d a m³/s con A = 2 171.6 km². La condición "
        "de cada mes (completo, escalado, excluido, sin dato) aparece al pasar el cursor.",
    "int:clim":
        "Los meses del periodo común se agrupan por mes calendario y en cada grupo se calculan media, mediana, cuartiles "
        "y percentiles 10 y 90 (interpolación lineal). Para T se usa PISCOt y, desde 2021, la estimación con ERA5-Land.",
    "int:anom":
        "Anomalía = valor del mes − media del mismo mes calendario en el periodo común. La media móvil centrada exige al "
        "menos 10 de los 12 meses. Para T se usa PISCOt y, desde 2021, la estimación con ERA5-Land.",
    "int:anomanual":
        "Suma de las 12 anomalías mensuales de cada año calendario completo del periodo común (media para T); usa la "
        "variable elegida en la gráfica anterior.",
    "int:modelos":
        "Predicciones de la evaluación por bloques de años (script 31): partición B (ajuste 2007–2024, evaluación "
        "1998–2006) y partición A (ajuste 1998–2015, evaluación 2016–2024). Referencia = climatología del mes calculada "
        "con los años de ajuste. Modelos elegidos por validación cruzada dentro del ajuste.",
    "int:tend":
        "Anomalía = valor del mes − media del mismo mes en 1998–2020. Barras = anomalía mensual; línea gris = media móvil "
        "de 12 meses; recta = tendencia de Sen con prueba de Mann–Kendall corregida por autocorrelación (Hamed–Rao), del "
        "script 32. La ordenada de la recta es la mediana de a − b·t.",
}


def tarjeta(titulo, tipo, cuerpo, que, como, leer_, meta, clase="fig"):
    """Tarjeta de figura (diseño del visor de la compañera): título y tipo arriba, figura a todo el ancho y, debajo,
    «Qué muestra» y «Cómo se construyó» en dos columnas y «Cómo leerla» en una fila completa."""
    extra = f'<div class="f-leer"><h5>Cómo leerla</h5><p class="leer">{leer_}</p></div>' if leer_ else ""
    chip = '<span class="tipo int">interactiva</span>' if tipo == "int" else '<span class="tipo">imagen</span>'
    return f"""
<figure class="{clase}">
  <header class="fig-h"><h4>{titulo}</h4>{chip}</header>
  {cuerpo}
  <figcaption>
    <div class="ficha">
      <div class="f-que"><h5>Qué muestra</h5><p>{que}</p></div>
      <div class="f-como"><h5>Cómo se construyó</h5><p>{como}</p></div>{extra}
    </div>
    <p class="meta">{meta}</p>
  </figcaption>
</figure>"""


def figura(p, titulo, que, leer_=None, script=None):
    p = Path(p)
    if p.name not in COMO:
        raise KeyError(f"Falta «Cómo se construyó» para {p.name} en COMO")
    pie = f'<span class="script">{html.escape(script)}</span>' if script else ""
    cuerpo = (f'<div class="fig-cuerpo"><button class="zoom" aria-label="Ampliar figura: {html.escape(titulo)}">'
              f'<img loading="lazy" src="{img64(p)}" alt="{html.escape(titulo)}"></button></div>')
    return tarjeta(titulo, "img", cuerpo, que, COMO[p.name], leer_,
                   f'{pie}<a href="{rel(p)}">PNG original</a> <span class="nota">(disponible si se clona el repositorio)</span>')


def selector(grafica, opciones, multi=False, etiqueta="Variable"):
    """Botones para elegir la variable de una gráfica interactiva (multi: varias a la vez)."""
    b = "".join(f'<button type="button" data-v="{v}" aria-pressed="{"true" if multi or i == 0 else "false"}">{t}</button>'
                for i, (v, t) in enumerate(opciones))
    tit = f'<span class="sel-tit">{etiqueta}:</span>' if multi else ""
    return (f'<div class="selector{" multi" if multi else ""}" data-para="{grafica}" role="group" '
            f'aria-label="{etiqueta}">{tit}{b}</div>')


def grafica(nombre, titulo, que, leer_=None, script=None, controles="", lado=""):
    """Gráfica interactiva (Plotly, visor_graficas.js). «Cómo se construyó» en COMO["int:<nombre>"].
    lado: panel lateral junto a la gráfica (control de calidad)."""
    clave = f"int:{nombre}"
    if clave not in COMO:
        raise KeyError(f"Falta «Cómo se construyó» para la gráfica interactiva {nombre} en COMO")
    plot = f'<div class="plot" data-plot="{nombre}" role="img" aria-label="{html.escape(titulo)}"></div>'
    cuerpo = (f'{controles}<div class="qc-grid"><div class="fig-cuerpo">{plot}</div>{lado}</div>' if lado
              else f'<div class="fig-cuerpo">{controles}{plot}</div>')
    pie = f'<span class="script">{html.escape(script)}</span>' if script else ""
    return tarjeta(titulo, "int", cuerpo, que, COMO[clave], leer_, pie)


TRAT = {
    "P_L": "Sin días faltantes en 1981–2025. P<sub>L</sub> mensual = suma de los días del mes, solo con meses completos. "
           "No es una observación directa: es PISCOp v2.1 (pluviómetros de SENAMHI + CHIRPS) promediado sobre la cuenca "
           "por CAMELS-PE, y no se pudo reproducir con PISCOp v3.0 (−14.5 %). Los meses con menos de 1 mm son meses secos "
           "reales, no faltantes.",
    "Q": f"{_rd.loc['Q', 'pct_dias_faltantes']:.1f} % de días faltantes en el periodo común, concentrados en ago–oct y en "
         f"2019–2021. Con el criterio de completitud (≤ 10 % de días faltantes) quedaban {RL['excl']} meses excluidos; el "
         f"grupo decidió rellenarlos: (1) huecos ≤ 15 días, interpolación de ln q; (2) meses con ≥ 50 % de días, media × "
         f"días del mes; (3) el resto, modelo mensual con el caudal simulado de CAMELS-PE. Resultado: {RL['obs']} meses "
         f"observados completos, {RL['interp']} con días interpolados, {RL['esc']} escalados y {RL['mod']} reconstruido "
         f"({RL['mes_mod']}). Lo observado se conserva aparte (R_obs, R_flag_obs) y la bandera R_flag indica el origen de "
         f"cada mes. Validación en la figura de abajo.",
    "T": f"PISCOt (CAMELS-PE) termina en 2020-12: los {int(_rd.loc['T', 'sin_dato'])} meses desde 2021 no tienen dato "
         f"observado. Se estimaron con ERA5-Land + sesgo de cada mes calendario (MAE {float(val['MAE_C']):.2f} °C en la "
         f"validación {val['evaluacion']}) y se guardan en una columna aparte, sin unir las fuentes. En las estadísticas "
         f"del punto 1.3 solo entra PISCOt.",
    "P_I": f"Producto mensual (IMERG V07 Final): los {N_MESES} meses del periodo común están completos. La tasa media "
           f"(mm/h) se convirtió a mm/mes con los días de cada mes y se promedió sobre la cuenca ponderando por área. "
           f"Nunca baja de ≈ 2 mm/mes y sus picos son mucho menores que los de P<sub>L</sub> (limitación del producto, "
           f"no error de datos).",
}


def panel_qc(v, primero):
    nombre, patron = _VQC[v]
    r = _rd.loc[v]
    falt = "—" if pd.isna(r.pct_dias_faltantes) else f"{r.pct_dias_faltantes:.1f} %"
    chips = "".join(f'<div class="dato">{t}<b>{x}</b></div>' for t, x in [
        ("Días faltantes", falt), ("Meses completos", int(r.completo)), ("Escalados", int(r.escalado)),
        ("Excluidos", int(r.excluido)), ("Sin dato", int(r.sin_dato))])
    filas = reg[reg.variable.str.contains(patron, regex=True)]
    items = "".join(f'<li><span class="estado e-{html.escape(x.estado.replace(" ", "-"))}">{html.escape(x.estado)}</span> '
                    f'<b>{html.escape(x.id)}</b> · {html.escape(str(x.fecha_o_periodo))}: {html.escape(x.anomalia)} '
                    f'<i>→ {html.escape(x.decision)}</i></li>' for _, x in filas.iterrows())
    return (f'<div class="qc-var" data-var="{v}"{"" if primero else " hidden"}>'
            f'<h5>{nombre}</h5><p class="nota">Periodo común {PER}, antes del relleno.</p>'
            f'<div class="datos-clave">{chips}</div>'
            f'<h5>Qué se hizo</h5><p>{TRAT[v]}</p>'
            f'<h5>Comprobaciones y decisiones</h5><ul class="qc-lista">{items}</ul></div>')


def lista_morfometria():
    """Lista horizontal de los parámetros morfométricos (script 06)."""
    pm = leer(DL / "parametros_morfometricos.csv", index_col=0)["valor"]
    filas = [  # (clave en el CSV, etiqueta, formato, unidad)
        ("Área A (km2)", "Área A", "{:,.1f}", "km²"),
        ("Área CAMELS-PE (km2)", "Área CAMELS-PE", "{:,.1f}", "km²"),
        ("IoU con polígono CAMELS-PE", "Coincidencia con CAMELS-PE (IoU)", "{:.3f}", ""),
        ("Perímetro suavizado P (km)", "Perímetro P (suavizado)", "{:.1f}", "km"),
        ("Coef. de compacidad Kc = 0.282 P/sqrt(A)", "Coef. de compacidad Kc = 0.282·P/√A", "{:.2f}", ""),
        ("Factor de forma Kf = A/L^2", "Factor de forma Kf = A/L²", "{:.3f}", ""),
        ("Longitud río principal L (km)", "Longitud del río principal L", "{:.1f}", "km"),
        ("Pendiente media río principal (m/m)", "Pendiente media del río principal", "{:.3f}", "m/m"),
        ("Cota naciente río principal (m)", "Cota de la naciente", "{:,.0f}", "m s.n.m."),
        ("Cota mínima (m s.n.m.)", "Cota mínima (cierre)", "{:,.0f}", "m s.n.m."),
        ("Cota máxima (m s.n.m.)", "Cota máxima", "{:,.0f}", "m s.n.m."),
        ("Cota media (m s.n.m.)", "Cota media", "{:,.0f}", "m s.n.m."),
        ("Cota mediana (m s.n.m.)", "Cota mediana", "{:,.0f}", "m s.n.m."),
        ("Pendiente media de la cuenca (°)", "Pendiente media de la cuenca", "{:.1f}", "°"),
        ("Longitud red de drenaje, umbral 5.0 km2 (km)", "Longitud de la red (A > 5 km²)", "{:,.0f}", "km"),
        ("Densidad de drenaje Dd (km/km2)", "Densidad de drenaje Dd", "{:.2f}", "km/km²"),
        ("Orden de Strahler máximo", "Orden de Strahler máximo", "{:.0f}", ""),
    ]
    items = "".join(f'<div class="param"><span>{e}</span><b>{f.format(pm[k]).replace(",", " ")} {u}</b></div>'
                    for k, e, f, u in filas)
    return (f'<div class="tabla-cont"><h4>Parámetros morfométricos</h4><div class="params">{items}</div>'
            f'<p class="nota">Calculados en UTM 17S con la delimitación GLO-30 (script 06; '
            f'parametros_morfometricos.csv). Kc ≈ 1 indica una cuenca circular; valores mayores, más alargada. '
            f'La red de drenaje depende del umbral de área elegido (5 km²), así que Dd es relativa a ese umbral.</p></div>')


def tabla(df, titulo, nota=""):
    t = df.to_html(index=False, border=0, classes="tabla", na_rep="—", escape=True)
    n = f'<p class="nota">{nota}</p>' if nota else ""
    return f'<div class="tabla-cont"><h4>{titulo}</h4><div class="scroll">{t}</div>{n}</div>'


def seccion(id_, titulo, intro, cuerpo, punto=""):
    chip = f'<span class="chip">{punto}</span>' if punto else ""
    return f'<section id="{id_}"><h2>{chip}{titulo}</h2><p class="intro">{intro}</p>{cuerpo}</section>'


# ------------------------------------------------------------------ contenido
fuentes = pd.DataFrame([
    ["P_L – precipitación de referencia", "CAMELS-PE v1.0.1 (`prec`, PISCOp v2.1 promediado sobre la cuenca)", "doi:10.5281/zenodo.21195425", "diaria → mensual"],
    ["Q / R – caudal observado", "CAMELS-PE v1.0.1 (`flow_obs`, SENAMHI)", "doi:10.5281/zenodo.21195425", "diaria → mensual"],
    ["T – temperatura", "CAMELS-PE (`tmean`, PISCOt v1.2) hasta 2020; ERA5-Land T2m desde 2021 (estimación aparte)", "doi:10.24381/cds.68d2bb30", "diaria / mensual"],
    ["P_I – precipitación satelital", "GPM IMERG V07 Final (GPM_3IMERGM.07, NASA GES DISC)", "doi:10.5067/GPM/IMERG/3B-MONTH/07", "mensual, 0.1°"],
    ["Precipitación en malla (apoyo)", "PISCOp v3.0 (SENAMHI)", "doi:10.6084/m9.figshare.32411886.v1", "diaria, 0.1°"],
    ["Mapa de precipitación media anual", "CHIRPS v2.0 (UCSB/CHC), recorte del IRI Data Library", "doi:10.15780/G2RP4Q", "mensual, 0.05°"],
    ["Transporte de humedad y viento", "ERA5 (C3S/ECMWF), niveles de presión y nivel único", "doi:10.24381/cds.6860a573; 10.24381/cds.f17050d7", "mensual, 0.5°"],
    ["Perfiles verticales (inversión térmica)", "Radiosondeos IGRA2 84516 Trujillo (NOAA NCEI); ERA5 por niveles de presión, 12 UTC",
     "doi:10.1175/JTECH-D-17-0223.1; 10.24381/cds.6860a573", "diaria 12 UTC / mensual, 0.25°"],
    ["Relieve", "Copernicus GLO-30 DSM", "ESA / Airbus", "30 m"],
    ["Cobertura y geología", "ESA WorldCover 2021; INGEMMET 1:100 000 (GEOCATMIN)", "—", "10 m / vectorial"],
], columns=["Variable", "Fuente", "DOI / proveedor", "Resolución"])

t13 = est.reset_index()[["variable", "unidades", "n_meses_validos", "media", "mediana", "desv_est", "minimo",
                         "maximo", "Q1", "Q3", "P10", "P90", "asimetria_G1"]].round(2)
t13.columns = ["Variable", "Unidades", "Meses válidos", "Media", "Mediana", "Desv. est.", "Mín.", "Máx.",
               "Q1", "Q3", "P10", "P90", "Asimetría G1"]
t15 = ip.reset_index()[["variable", "mes_max", "mes_min", "SI", "PCI", "centroide", "concentracion", "picos",
                        "meses_humedos", "meses_secos_gaussen", "pct_en_dic_abr", "regimen"]].round(2)
t15.columns = ["Variable", "Mes máx.", "Mes mín.", "SI", "PCI", "Centroide", "Concentración", "Picos",
               "Meses húmedos", "Meses secos (Gaussen)", "% en dic–abr", "Régimen"]
t15c = leer(CL / "criterios_clasificacion_1_5c.csv").drop(columns="fuente")
t15c.columns = ["Rasgo", "Criterio", "Evidencia (P_L)", "Evidencia alternativa", "Lectura propuesta",
                "Robustez / incertidumbre"]
ari = leer(CL / "aridez_budyko_1_5c.csv")
ari = ari[ari.fuente_PET.str.startswith("Penman")].set_index("fuente_P")
t14 = reg[["id", "variable", "fecha_o_periodo", "anomalia", "decision", "estado"]]
t14.columns = ["Id", "Variable", "Fecha / periodo", "Anomalía", "Decisión", "Estado"]

secciones = [
    seccion("cuenca", "La cuenca y su delimitación",
            "Río Chicama hasta la estación hidrométrica El Tambo (SENAMHI 211003; CAMELS-PE PE_47E0A2A8), en la "
            "vertiente del Pacífico del norte del Perú (−7.574°, −78.707°; 684 m s.n.m.). La cuenca se delimitó "
            "con el modelo digital de elevación Copernicus GLO-30 y se comparó con el polígono oficial de CAMELS-PE.",
            figura(DL / "mapa_cuenca_el_tambo.png", "Mapa de la cuenca",
                   "Relieve (GLO-30), límite de la cuenca, río Chicama, red de drenaje, punto de cierre y estaciones de SENAMHI. "
                   "Área delimitada 2 170 km² (CAMELS-PE: 2 171.6 km²; índice de intersección sobre unión 0.994), con "
                   "cotas entre unos 680 y 4 300 m s.n.m. Los recuadros ubican la cuenca en Sudamérica y en el Perú, "
                   "con Lima como referencia.",
                   script="06_delimitacion_el_tambo_glo30.py")
            + lista_morfometria()),
    seccion("contexto", "Contexto físico",
            "Rasgos del terreno, de la cobertura y de la atmósfera que ayudan a explicar cómo la cuenca recibe, "
            "guarda y entrega el agua. Haga clic en una figura para ampliarla.",
            figura(CF / "panel_contexto_fisico.png", "Panel de contexto físico",
                   "(a) orientación de laderas, (b) cobertura del suelo, (c) unidades geológicas de INGEMMET "
                   "agrupadas por su comportamiento hídrico y (d) precipitación media anual de CHIRPS. En (c), los "
                   "acuíferos fisurados (areniscas y cuarcitas, 40 %) y los acuitardos sedimentarios (lutitas, 37 %) "
                   "alternan en franjas NO–SE; los acuíferos porosos se limitan a los fondos de valle (2.5 %).",
                   "Colores de (c) según la convención de los mapas hidrogeológicos (IAH/UNESCO): azul = acuífero "
                   "poroso, verdes = acuíferos fisurados y kársticos, marrones = acuitardos y acuífugos. La agrupación "
                   "es una propuesta del grupo (por porosidad y permeabilidad esperadas) que debe contrastarse con "
                   "INGEMMET y ANA.", script="07b_mapas_contexto_fisico.py")
            + f'<h3 class="sub2">Fuentes de humedad e inversión térmica</h3><p class="intro">Radiosondeos observados '
              f'de Trujillo (IGRA2 84516, a ~40 km de la salida de la cuenca): {SND["n"]} sondeos de las 7 a.m. entre '
              f'{SND["ini"]:%Y-%m} y {SND["fin"]:%Y-%m}, con inversión en el {SND["frec"]:.0f} % de ellos; y transporte '
              f'de humedad de ERA5 en el periodo común. Los sondeos solo existen desde {SND["ini"]:%Y-%m}, fuera del '
              f'periodo común ({PER}); al final de la sección se comprueba con ERA5 si representan ese periodo.</p>'
            + figura(IV / "fig_sondeo_tipico_inversion.png", "Sondeo típico de la estación seca",
                     "La temperatura baja con la altura cerca de la adiabática seca hasta el nivel de condensación (NCA), "
                     "donde se forman los estratocúmulos (HR ≈ 100 %). En la base de la inversión el aire se calienta de "
                     "golpe, unos 5 °C, y por encima la humedad cae a menos del 50 %. Debajo queda el aire marino frío y "
                     "húmedo; encima, aire de la troposfera libre que desciende y se calienta. La capa de nubes queda "
                     "atrapada a la altura del cierre de la cuenca y no llega a la cuenca alta.",
                     script="34a_descarga_perfiles.py, 34b_perfiles_inversion.py")
            + figura(IV / "fig_temporadas_humedad_inversion.png",
                     "De dónde llega la humedad y qué la frena: estación seca frente a estación húmeda",
                     f"<b>La diferencia entre temporadas no está en cuánto vapor pasa sobre la cuenca, sino en si puede "
                     f"subir.</b> En ERA5 el transporte sobre la cuenca va hacia el oeste en las dos temporadas y no es "
                     f"menor en la seca (≈ {IVT['may-sep'][0]:.0f} frente a {IVT['dic-abr'][0]:.0f} kg m⁻¹ s⁻¹), aunque "
                     f"la columna lleva menos agua precipitable (≈ {IVT['may-sep'][1]:.0f} frente a "
                     f"{IVT['dic-abr'][1]:.0f} mm). <b>Estación seca (may–sep):</b> la inversión es fuerte "
                     f"(+{SND['seca']['dT']:.1f} °C) y está casi siempre ({SND['seca']['frec']:.0f} % de los días), con "
                     f"base a {SND['seca']['base']:.0f} m y tope a {SND['seca']['tope']:.0f} m, justo sobre el cierre: "
                     f"casi toda la cuenca queda en aire cálido y seco que desciende, y el vapor pasa sin condensar. "
                     f"Frente a la costa, el flujo del anticiclón del Pacífico corre paralelo a ella, de sureste a "
                     f"noroeste. <b>Estación húmeda (dic–abr):</b> el flujo del Atlántico cruza la Amazonía y al este de "
                     f"los Andes gira hacia el sur junto a la cordillera; la columna es más húmeda y la inversión es más "
                     f"baja y débil (+{SND['hum']['dT']:.1f} °C, {SND['hum']['base']:.0f}–{SND['hum']['tope']:.0f} m, "
                     f"por debajo del cierre) y falta algunos días ({SND['hum']['frec']:.0f} % con inversión). Sobre la "
                     f"cuenca el perfil se enfría con la altura sin tapa, lo que permite la convección.",
                     "El transporte de humedad no equivale a lluvia: llueve donde el flujo converge y asciende. Los "
                     "sondeos cubren dos años y son de las 7 a.m. La explicación física (anticiclón del Pacífico sur, "
                     "corriente de Humboldt, convección andina) debe sustentarse con bibliografía.",
                     script="24_era5_vientos_humedad.py, 34b_perfiles_inversion.py, 35_temporadas_humedad_inversion.py")
            + f'<h3 class="sub2">¿Representan los sondeos el periodo común?</h3><p class="intro">Los sondeos de '
              f'Trujillo cubren {SND["ini"]:%Y-%m} a {SND["fin"]:%Y-%m}, casi todo fuera del periodo común ({PER}). '
              f'Se separan dos preguntas: (1) si ERA5 reproduce los sondeos cuando se comparan los mismos meses, y '
              f'(2) si esos meses fueron típicos del periodo común según ERA5. Estabilidad de la capa baja = '
              f'T(1 500 m) − T(500 m); positiva indica inversión neta en esa capa.</p>'
            + figura(IV / "fig_contraste_1_trujillo_vs_era5.png",
                     "Comparación 1 · Trujillo (sondeos) frente a ERA5 en los mismos meses (2024-10 a 2026-09)",
                     f"Misma hora y mismos meses, así que la diferencia es error del reanálisis. Por encima de unos "
                     f"1 500 m ERA5 sigue de cerca a los sondeos en las dos temporadas (diferencias menores de 1 °C). "
                     f"Por debajo, ERA5 es más cálido, hasta unos 2 °C cerca de la superficie, y no forma la inversión: "
                     f"en may–sep la estabilidad de la capa baja es {CT['seca_snd']:+.1f} °C en los sondeos y "
                     f"{CT['seca_e5s']:+.1f} °C en ERA5; en dic–abr, {CT['hum_snd']:+.1f} y {CT['hum_e5s']:+.1f} °C. "
                     f"La caída de la humedad por encima de la capa en may–sep sí aparece en ERA5, aunque más gradual.",
                     "ERA5 es un promedio mensual en niveles de presión (cada ~25–50 hPa) y una celda de 0.25° que mezcla "
                     "mar y tierra: suaviza la inversión y no sirve para medir su intensidad. Sí distingue las dos "
                     "temporadas por la forma del perfil: capa casi isoterma en may–sep, enfriamiento continuo en dic–abr.",
                     script="34a_descarga_perfiles.py, 34b_perfiles_inversion.py, 34c_contraste_sondeos_era5.py")
            + figura(IV / "fig_contraste_2_era5_periodos.png",
                     f"Comparación 2 · ERA5: meses de los sondeos (2024-10 a 2026-09) frente al periodo común ({PER})",
                     f"Una sola fuente, así que la diferencia es del periodo. En dic–abr los meses con sondeos son casi "
                     f"iguales al periodo común (diferencias de unas décimas de grado). En may–sep fueron más cálidos, unos "
                     f"0.5–1 °C cerca de la superficie y por encima de 1 500 m (casi sin diferencia hacia los 1 000 m), "
                     f"pero con la misma forma del perfil: estabilidad de la capa baja "
                     f"{CT['seca_e5s']:+.1f} °C frente a {CT['seca_e5c']:+.1f} °C en el periodo común.",
                     "Si la forma del perfil en los meses con sondeos es la del periodo común, la inversión que miden los "
                     "sondeos puede tomarse como representativa de 1998–2024 en su estructura (posición respecto a la "
                     "cuenca, contraste entre temporadas), no en sus valores exactos. Dos años son pocos; el grupo debe "
                     "decidir cuánto peso darle a esta comprobación.",
                     script="34a_descarga_perfiles.py, 34c_contraste_sondeos_era5.py")),
    seccion("p11", "Series mensuales de precipitación, caudal y temperatura",
            f"Series construidas desde los registros diarios: P<sub>L</sub> = suma diaria; Q̄ = media diaria (m³/s); "
            f"R = caudal en lámina (mm/mes); T = media diaria. Todas las figuras usan el <b>periodo común {PER}</b> "
            f"({N_MESES} meses): los meses con precipitación de CAMELS-PE y de IMERG a la vez. En las gráficas "
            f"interactivas, arrastre para hacer zoom, use la barra inferior para elegir un intervalo y haga doble clic "
            f"para volver a la vista completa.",
            grafica("serie11", "Series mensuales en paneles alineados",
                    f"Los cuatro paneles comparten el eje de tiempo: al pasar el cursor se leen los valores del mismo mes. "
                    f"El caudal y la escorrentía tienen sus faltantes rellenados ({RL['interp']} meses con días "
                    f"interpolados, {RL['esc']} escalados y {RL['mod']} reconstruido con modelo, {RL['mes_mod']}; el "
                    f"método de cada mes está en la columna R_flag del CSV); el fondo rojo claro marca los {RL['excl']} "
                    f"meses que sin relleno quedaban excluidos. Desde 2021 la temperatura es una estimación con ERA5-Land "
                    f"(línea discontinua), no PISCOt. Con los botones se eligen los paneles visibles.",
                    "Los máximos de 1998 y 2017 corresponden a eventos de El Niño; el caudal queda alto también en los "
                    "meses siguientes.", script="09_datos_camels_pe_el_tambo.py, relleno_caudal.py",
                    controles=selector("serie11", [("P_L", "P<sub>L</sub>"), ("Q", "Q̄"), ("R", "R"), ("T", "T")],
                                       multi=True, etiqueta="Paneles visibles"))
            + grafica("temp", "Temperatura: PISCOt frente a ERA5-Land",
                      f"PISCOt (la temperatura de CAMELS-PE) termina en 2020. En los {TB['n']} meses con ambas fuentes "
                      f"varían igual (r = {TB['r']:.2f}), pero ERA5-Land es {abs(TB['sesgo']):.1f} °C más fría de forma "
                      f"casi constante. La estimación desde 2021 (ERA5-Land + sesgo de cada mes) tuvo un error medio de "
                      f"{float(val['MAE_C']):.2f} °C al validarla en {val['evaluacion']}. Las dos fuentes no se unen en "
                      f"una sola serie. El ciclo anual (b) tiene la misma forma en las dos: casi plano, con el máximo en "
                      f"agosto–septiembre.",
                      "Solo unos 0.5 °C del desfase se explican por la altitud del relieve del modelo. Otra causa probable "
                      "es la definición: PISCOt usa (Tmáx + Tmín)/2, y ERA5-Land promedia las 24 horas.",
                      script="17_era5land_temperatura.py"), "Punto 1.1"),
    seccion("p12", "Precipitación satelital IMERG",
            f"IMERG V07 Final (0.1°) promediado sobre la cuenca, ponderando cada celda por su área dentro de la cuenca. "
            f"Frente a P<sub>L</sub>: media {IM['PI']:.1f} vs {IM['PL']:.1f} mm/mes, sesgo {IM['sesgo']:+.1f} mm/mes, "
            f"MAE {IM['mae']:.1f}, RMSE {IM['rmse']:.1f}, r = {IM['r']:.2f}. IMERG acumula el {100 * IM['ratio']:.0f} % "
            f"de la lluvia de P<sub>L</sub>. Ninguna de las dos es una observación directa: P<sub>L</sub> interpola "
            f"pluviómetros con información satelital, e IMERG Final también se ajusta con pluviómetros.",
            figura(SM / "mapa_celdas_imerg_pesos.png", "Celdas de IMERG y pesos en el promedio de cuenca",
                   f"(a) Peso de cada celda = área dentro de la cuenca / área de la cuenca: {PW['n']} celdas, {PW['comp']} "
                   f"completas, {PW['borde']:.0f} % del peso en celdas de borde (≈ {PW['eq']:.1f} celdas equivalentes). "
                   "(b) Desnivel del terreno dentro de cada celda: una sola celda de IMERG abarca más de 1 000 m de "
                   "desnivel; los pesos no consideran la topografía.", script="12_mapa_pesos_imerg.py")
            + grafica("precip", "P<sub>L</sub> frente a IMERG, mes a mes",
                      "IMERG sigue bien el ritmo de la lluvia, pero sus picos son mucho más bajos, sobre todo en los "
                      "meses más lluviosos. CHIRPS y PISCOp v3.0 están ocultas: haga clic en su nombre en la leyenda "
                      "para compararlas.",
                      script="10_imerg_mensual_y_grafica.py, 26_chirps_cuenca.py, 14_piscop_v3_cuenca.py")
            + grafica("anual", "Totales anuales de P<sub>L</sub> e IMERG",
                      f"IMERG queda por debajo de P<sub>L</sub> todos los años; las mayores diferencias son las de "
                      f"{', '.join(str(a) for a in (_anual.P_L - _anual.P_I).nlargest(3).index)} "
                      f"({(_anual.P_L - _anual.P_I).nlargest(3).min():.0f} mm/año o más).", script="10_imerg_mensual_y_grafica.py, 11_grafica_precipitacion_camels.py"),
            "Punto 1.2"),
    seccion("p13", "Distribución de los valores mensuales",
            "Histogramas, diagramas de caja y estadísticos descriptivos en el periodo común. P<sub>L</sub> y P<sub>I</sub> "
            "usan los mismos meses y los mismos límites de clase.",
            figura(SM / "fig_histogramas_1_3.png", "Histogramas",
                   "Ancho de clase por la regla de Freedman–Diaconis y frecuencia relativa (% de meses). La lluvia y "
                   "sobre todo el caudal son asimétricos: muchos meses bajos y pocos muy altos. La temperatura es casi "
                   "simétrica.", script="18_estadisticos_1_3.py")
            + figura(SM / "fig_cajas_1_3.png", "Diagramas de caja",
                     "Caja = cuartiles 1 a 3, línea = mediana, rombo = media; bigotes de Tukey (1.5 × IQR). Los puntos "
                     "fuera de los bigotes se muestran, no se eliminan.", script="18_estadisticos_1_3.py")
            + tabla(t13, f"Estadísticos en el periodo común ({PER})",
                    "Percentiles con interpolación lineal (tipo 7). T solo PISCOt, sin la estimación desde 2021. "
                    "La tabla del registro completo de cada variable está en estadisticos_1_3_registro_completo.csv."),
            "Punto 1.3"),
    seccion("p14", "Control de calidad",
            "Comprobaciones sobre los datos diarios (fechas, códigos de faltante, negativos, secuencias constantes, "
            "tramos interpolados, picos aislados), revisiones manuales y registro de anomalías con la decisión tomada.",
            grafica("qc", "Disponibilidad y control de calidad por variable",
                    "Arriba, la condición de cada mes (año × mes) antes de cualquier relleno: azul completo, celeste "
                    "escalado (≤ 10 % de días faltantes), rojo excluido (> 10 %) y gris sin dato; el número es la cantidad "
                    "de días con dato en los meses incompletos. Abajo, la serie diaria con los vacíos sin rellenar; en el "
                    "caudal, los días interpolados se marcan con puntos negros. A la derecha, lo que se comprobó y se "
                    "decidió con esa variable. Los faltantes de caudal se concentran en septiembre (y agosto–octubre) casi "
                    "todos los años, además de un periodo con más vacíos en 2019–2021.",
                    script="19_control_calidad_1_4.py, 09_datos_camels_pe_el_tambo.py",
                    controles=selector("qc", [("P_L", "P<sub>L</sub>"), ("Q", "Q̄ caudal"), ("T", "T"),
                                              ("P_I", "P<sub>I</sub> IMERG")]),
                    lado=f'<aside class="qc-lado">{"".join(panel_qc(v, i == 0) for i, v in enumerate(_VQC))}</aside>')
            + figura(SM / "fig_validacion_relleno_caudal.png", "Relleno de los faltantes de caudal y su validación",
                     f"Decisión del grupo: rellenar los faltantes de caudal. Método, en orden: (1) huecos diarios de hasta "
                     f"15 días, interpolando el logaritmo del caudal entre el día anterior y el posterior, lo que sigue la "
                     f"forma de la recesión; (2) meses con datos en al menos la mitad de los días, media de los días "
                     f"disponibles × días del mes; (3) meses con menos datos, modelo mensual con el caudal simulado de "
                     f"CAMELS-PE y el caudal del mes anterior. En el periodo común: {RL['interp']} meses con días "
                     f"interpolados, {RL['esc']} escalados y {RL['mod']} reconstruido(s) con el modelo ({RL['mes_mod']}).",
                     "(a) Al borrar días de meses completos, la interpolación da errores 3–8 veces menores que el escalado "
                     "usado antes (≤ 1.5 % de mediana hasta 7 días). (b) El modelo mensual, evaluado dejando un año "
                     "fuera, tiene un error típico del 23 %; por eso solo se usa donde no hay suficientes días medidos. "
                     "La bandera R_flag de series_mensuales_camels_pe.csv indica el origen de cada mes.",
                     script="relleno_caudal.py, 09_datos_camels_pe_el_tambo.py, 28_validacion_relleno_caudal.py")
            + tabla(t14, "Registro de anomalías",
                    "Estado: «sin problema» = comprobado, nada que corregir; «incierto» = la duda se declara como "
                    "limitación. El detalle de cada comprobación está en registro_anomalias_1_4.csv."),
            "Punto 1.4"),
    seccion("p15a", "Climatología mensual",
            f"Climatología de 12 meses en el periodo común: para cada mes, número de años, media, mediana, desviación "
            f"estándar, cuartiles y percentiles 10 y 90. P<sub>L</sub> llega a su máximo en marzo "
            f"({cmes('P_L', 3)['media']:.0f} mm de media) y a su mínimo en julio–agosto; el caudal llega a su máximo en "
            f"abril ({cmes('Q', 4)['media']:.0f} m³/s).",
            grafica("clim", "Ciclo anual medio con bandas de dispersión",
                    "Línea continua = media; discontinua = mediana; banda oscura = rango intercuartílico (Q1–Q3); banda "
                    "clara = percentiles 10–90. Elija la variable con los botones.",
                    "Las bandas describen cuánto cambian los años entre sí; no son intervalos de confianza.",
                    script="20_climatologia_1_5a.py",
                    controles=selector("clim", [("P_L", "P<sub>L</sub>"), ("P_I", "P<sub>I</sub>"), ("Q", "Q̄"),
                                                ("R", "R"), ("T", "T")]))
            + figura(CL / "mapas_precipitacion_mensual_chirps.png", "Precipitación media de cada mes sobre la cuenca",
                     "Doce mapas con la precipitación media de cada mes calendario (CHIRPS, celdas de 0.05°), con una "
                     "escala discreta de naranja (seco) a azul oscuro (húmedo). Cada fila es un trimestre: DEF y MAM "
                     "aportan cada uno cerca del 42–43 % de la lluvia anual de P<sub>L</sub>, SON el 14 % y JJA el 2 %. "
                     "La lluvia crece de dic a mar, cuando cubre toda la cuenca, y en jun–ago casi desaparece. En los "
                     "meses húmedos llueve más en las cabeceras del sureste y menos cerca de la salida.",
                     "El recuadro de cada mapa compara el promedio de cuenca de CHIRPS con el de P<sub>L</sub>. CHIRPS "
                     "es menor en el total anual, así que los mapas muestran el patrón espacial, no la magnitud.",
                     script="27_mapas_precipitacion_mensual.py")
            + figura(CL / "fig_anio_mes_1_5a.png", "Mapa año × mes",
                     "Cada celda es un mes de un año, de modo que se ven los años individuales. El caudal va en escala "
                     "logarítmica. En temperatura, la trama marca los valores estimados con ERA5-Land.",
                     script="20_climatologia_1_5a.py"),
            "Punto 1.5.a"),
    seccion("p15b", "Rasgos del ciclo anual y su variabilidad",
            f"Forma, fase, estacionalidad y variabilidad entre años del ciclo anual. Con los criterios propuestos, el "
            f"régimen de P<sub>L</sub> es <b>{ip.loc['P_L', 'regimen']}</b>.",
            figura(CL / "fig_ciclo_1_5b.png", "Rasgos del ciclo anual",
                     f"(a) Forma y fase: el caudal va {lag_centroide('P_L'):.1f} meses detrás de P<sub>L</sub> según el "
                     f"centroide. (b) El máximo anual de P<sub>L</sub> cae en marzo en el {pct_mar:.0f} % de los años "
                     f"(feb–abr en el {pct_fa:.0f} %), y el del caudal en abril en el {pct_qabr:.0f} %. (c) Los dos "
                     f"subperiodos mantienen el mismo régimen. (d) Variabilidad entre años por mes; en los meses secos "
                     f"el CV es inestable porque la media es casi cero.",
                     "SI = índice de estacionalidad (Walsh & Lawler, 1981); PCI = índice de concentración (Oliver, 1980); "
                     "centroide = Markham (1970).", script="21_ciclo_anual_1_5b.py")
            + grafica("anom", "Anomalías mensuales respecto al ciclo anual",
                      "Barras = diferencia de cada mes respecto a la media de su mes calendario (azul por encima de lo "
                      "normal, rojo por debajo; en T, rojo = más cálido). Línea = media móvil centrada de 12 meses, que "
                      "deja ver los años húmedos y secos sin el ciclo estacional.",
                      script="21_ciclo_anual_1_5b.py",
                      controles=selector("anom", [("P_L", "P<sub>L</sub>"), ("P_I", "P<sub>I</sub>"), ("R", "R"),
                                                  ("T", "T")]))
            + grafica("anomanual", "Anomalía anual (año calendario)",
                      f"Suma de las anomalías de cada año (media para T), para la variable elegida arriba. Años más "
                      f"húmedos según P<sub>L</sub>: {humedos}; más secos: {secos}.",
                      script="21_ciclo_anual_1_5b.py")
            + tabla(t15, "Indicadores de estacionalidad (periodo común)",
                    "R = caudal en lámina. Mes seco según Bagnouls–Gaussen (P < 2T). Criterios de régimen: unimodal si "
                    "hay un solo pico y el armónico anual domina; bimodal si hay dos picos y el semianual es comparable; "
                    "estacionalidad débil si SI < 0.40."),
            "Punto 1.5.b"),
    seccion("p15c", "Síntesis del ciclo anual",
            f"Evidencia reunida para proponer la clasificación hidroclimática mensual. Para la aridez se usa un "
            f"criterio climático (P/PET), no solo el total de lluvia. La clasificación final y su explicación física "
            f"las decide y sustenta el grupo con bibliografía.",
            figura(CL / "fig_sintesis_ciclo_1_5c.png", "Ciclo anual frente a la demanda de agua",
                   f"Barras: lluvia media de cada mes (P<sub>L</sub> e IMERG). Líneas: evapotranspiración potencial (PET, "
                   f"Penman–Monteith de CAMELS-PE, 1998–2016), la mitad de la PET y el caudal en lámina (R). Debajo, cada "
                   f"mes clasificado por varios criterios. Con P<sub>L</sub>, la lluvia supera a la PET solo en feb–abr, y "
                   f"de mayo a agosto el caudal medio supera a la lluvia del mismo mes.",
                   f"Fila FAO: húmedo si P ≥ PET, intermedio si P ≥ 0.5·PET, seco si es menor; el porcentaje es la "
                   f"fracción de años en que P<sub>L</sub> ≥ 0.5·PET. Aridez UNEP (P/PET): P<sub>L</sub> "
                   f"{ari.loc['P_L (CAMELS-PE)', 'AI_P_PET']:.2f} ({ari.loc['P_L (CAMELS-PE)', 'clase_UNEP']}), "
                   f"IMERG {ari.loc['P_I (IMERG)', 'AI_P_PET']:.2f} ({ari.loc['P_I (IMERG)', 'clase_UNEP']}): la clase "
                   f"depende de la fuente de lluvia.", script="23_sintesis_ciclo_1_5c.py")
            + tabla(t15c, "Criterios, evidencia e incertidumbre",
                    "«Lectura propuesta» es un apoyo cuantitativo para la discusión del grupo, no la conclusión. "
                    "Detalle en criterios_clasificacion_1_5c.csv, aridez_budyko_1_5c.csv y koppen_1_5c.csv."),
            "Punto 1.5.c"),
    seccion("p2", "Relaciones entre series y modelos estadísticos",
            f"Comparaciones mensuales con los meses coincidentes del periodo común ({PER}). El caudal usa la serie "
            f"rellenada; las métricas sin los meses rellenados están en los CSV como prueba de sensibilidad.",
            figura(RE / "fig_dispersion_2_1.png", "Dispersión entre lluvia y caudal",
                   f"(a) IMERG frente a P<sub>L</sub>: muy correlacionadas (r = {R21['a']:.2f}) pero no concordantes: "
                   f"IMERG queda bajo la línea 1:1 en los meses húmedos (sesgo {E21['sesgo']:+.1f} mm/mes, "
                   f"{E21['sesgo_pct']:+.0f} %) y algo por encima en jun–ago. (b, c) La relación de la lluvia del mismo "
                   f"mes con el caudal es débil (Spearman ρ = {R21['b']:.2f} con P<sub>L</sub> y {R21['c']:.2f} con IMERG) "
                   f"y muy dispersa: lluvias parecidas dan caudales que difieren hasta 100 veces según la época.",
                   "Colores por estación. El caudal va en escala logarítmica porque varía tres órdenes de magnitud. "
                   "(b') En lámina, casi todos los meses quedan bajo la línea R = P<sub>L</sub>; los que la superan "
                   "son de may–ago (aporte de almacenamiento) y crecidas de El Niño.", script="30_relaciones_2_1.py")
            + figura(RE / "fig_rezagos_2_1.png", "La lluvia de meses anteriores explica mejor el caudal",
                     f"El caudal se relaciona mejor con la lluvia del mes anterior (ρ = {RZ['k1']:.2f}) que con la del "
                     f"mismo mes ({RZ['k0']:.2f}), y aún mejor con la lluvia acumulada de los últimos 4 meses "
                     f"(ρ = {RZ['acu']:.2f}): la cuenca almacena agua y la entrega con retraso. En anomalías la asociación "
                     f"se debilita pero persiste (ρ = {RZ['an1']:.2f} con un mes de rezago).",
                     "La lluvia siempre antecede al caudal; nunca se usa lluvia futura.", script="30_relaciones_2_1.py")
            + grafica("modelos", "Observado frente a modelado, fuera del periodo de ajuste",
                      "Cada modelo se ajustó con unos años y se evaluó con otros que no vio; el hueco 2007–2015 no se "
                      "evalúa. Para la lluvia local, el factor por mes calendario corrige el sesgo de IMERG en los meses "
                      "húmedos; para el caudal, el modelo logarítmico sigue la recesión y el estiaje, pero subestima la "
                      "crecida de El Niño costero de 2017. Elija el objetivo con los botones; las series se activan o "
                      "desactivan en la leyenda.", script="31_modelos_2_2_2_3.py",
                      controles=selector("modelos", [("PL", "P<sub>L</sub> ← IMERG"), ("Q", "Q̄ ← P<sub>L</sub>")]))
            + figura(MO / "fig_evaluacion_resumen_2_3.png", "¿Mejoran los modelos a una referencia sencilla?",
                     "Habilidad fuera del periodo de ajuste en dos particiones por bloques de años (A: ajuste 1998–2015, "
                     "evaluación 2016–2024; B: ajuste 2007–2024, evaluación 1998–2006). Para P<sub>L</sub> a partir de "
                     "IMERG, un factor de corrección por mes calendario reduce el error típico cerca de 40 %. Para el "
                     "caudal, solo los modelos logarítmicos con término por mes y lluvia acumulada (2–3 meses de rezago) "
                     "mejoran la climatología; la regresión lineal simple es peor que la climatología.",
                     "★ = modelo elegido por validación cruzada dentro del ajuste. Habilidad = 1 − MAE/MAE de la referencia.",
                     script="31_modelos_2_2_2_3.py")
            + figura(MO / "fig_modelo_PL_2_3.png", "Modelo de lluvia local a partir de IMERG",
                     "Evaluación 2016–2024 del modelo elegido (factor por mes calendario) frente a IMERG sin corregir y a la "
                     "regresión lineal, con los residuos frente al tiempo, al valor estimado y al mes.",
                     script="31_modelos_2_2_2_3.py")
            + figura(MO / "fig_modelo_Q_2_3.png", "Modelo de caudal a partir de la lluvia local",
                     "Evaluación 2016–2024 del modelo elegido (ln Q̄ con intercepto por mes y lluvia de los últimos 4 meses). "
                     "Sigue bien la recesión y el estiaje, pero subestima la crecida de El Niño costero de 2017 (residuo de "
                     "−310 m³/s en abril). La regresión lineal da 30 caudales negativos (las caídas al borde inferior).",
                     "Una relación estadística lluvia–caudal no conserva masa ni sustituye el balance hídrico.",
                     script="31_modelos_2_2_2_3.py"),
            "Punto 2"),
    seccion("p3", "Tendencias hidroclimáticas",
            "Cada variable se analiza con su registro mensual completo (cambio de largo plazo) y, aparte, en el periodo "
            f"común {PER} (comparación entre fuentes). Las anomalías usan una climatología fija 1998–2020, el periodo más "
            "largo en que las cinco series tienen datos. Métodos: OLS con errores robustos a autocorrelación (HAC), "
            "Mann–Kendall con pendiente de Sen (estacional en la serie original; con corrección por autocorrelación de "
            "Hamed–Rao en las anomalías y en cada mes) y LOWESS como curva exploratoria. Pendientes por década; en "
            "precipitación son cambios del acumulado mensual, no del anual.",
            grafica("tend", "Anomalías y tendencia de Sen, registro completo",
                    f"La temperatura sube en las dos fuentes (PISCOt {pend3('T_PISCOt')} °C/década; ERA5-Land "
                    f"{pend3('T_ERA5L')} °C/década) de forma casi lineal. P<sub>L</sub> no muestra tendencia "
                    f"({pend3('P_L')} mm/mes/década) y el caudal tampoco ({pend3('Q')} m³/s/década). IMERG sí baja "
                    f"({pend3('P_I')} mm/mes/década), aunque la lluvia de referencia en el mismo periodo no lo hace. "
                    f"Elija la variable con los botones.", script="32_tendencias_3.py",
                    controles=selector("tend", [("P_L", "P<sub>L</sub>"), ("P_I", "P<sub>I</sub>"), ("Q", "Q̄"),
                                                ("T_PISCOt", "T PISCOt"), ("T_ERA5L", "T ERA5-Land")]))
            + figura(TE / "fig_series_tendencias_3_4.png", "Series originales, anomalías y anomalías estandarizadas",
                   f"Las mismas tendencias vistas en la serie original, en anomalías y en anomalías estandarizadas.",
                   "Columna izquierda: serie original con un modelo aditivo (efecto de cada mes + curva suave) y OLS con "
                   "efectos de mes. Centro: anomalía a = X − μ<sub>mes</sub> con banda de remuestreo por bloques. Derecha: "
                   "z = a/s<sub>mes</sub>. El panel de caudal tiene el eje recortado por la crecida de 1998.",
                   script="32_tendencias_3.py")
            + figura(TE / "fig_pendientes_mensuales_3_5.png", "Tendencia de cada mes calendario",
                     "Doce subseries por variable (todos los eneros, todos los febreros…). Tras controlar la tasa de falsos "
                     "descubrimientos (FDR) sobre 60 pruebas, las tendencias que se sostienen son de temperatura "
                     f"(ERA5-Land: {FDR.get('T_ERA5L', 'ninguno')}; PISCOt: {FDR.get('T_PISCOt', 'ninguno')}), con el "
                     "mayor calentamiento estandarizado al final de la estación seca y en la transición a la húmeda (ago–nov). En lluvia y caudal solo "
                     f"quedan casos aislados (P<sub>L</sub>: {FDR.get('P_L', 'ninguno')}; P<sub>I</sub>: "
                     f"{FDR.get('P_I', 'ninguno')}; Q̄: {FDR.get('Q', 'ninguno')}) que no son consistentes entre métodos.",
                     "Verificado: restar μ<sub>mes</sub> no cambia la pendiente OLS ni la de Sen, y dividir por s<sub>mes</sub> "
                     "solo cambia la escala; X, a y z no son evidencias independientes.", script="32_tendencias_3.py")
            + figura(TE / "fig_robustez_3_5.png", "Robustez frente a método, periodo y fuente",
                     "Las tendencias de temperatura persisten al cambiar el método, recortar los extremos del registro o "
                     "retirar los años de El Niño costero. En precipitación el signo depende de la fuente y del periodo: "
                     "P<sub>L</sub> y PISCOp v3 son ligeramente positivas en 1981–2025 y negativas, sin significancia, en "
                     "1998–2024, mientras IMERG y CHIRPS son negativas. La tendencia de IMERG desaparece si se empieza en 2003.",
                     "● Sen con MK Hamed–Rao; ■ OLS-HAC. Relleno = p < 0.05.", script="32_tendencias_3.py")
            + figura(TE / "fig_ventanas_inicio_fin_3_5.png", "Sensibilidad a la fecha inicial y final",
                     "Pendiente de Sen de las anomalías anuales para todas las ventanas de 15 años o más. En lluvia y caudal "
                     "el signo cambia según el año inicial, sobre todo si la ventana empieza o termina en un año de El Niño "
                     "(1983, 1998, 2017); en temperatura es positivo en casi todas las ventanas.",
                     script="32_tendencias_3.py")
            + figura(TE / "fig_causas_3_6.png", "¿A qué podrían deberse los cambios?",
                     f"Los años húmedos y de caudal alto coinciden con El Niño 1+2 cálido (r mensual con P<sub>L</sub> "
                     f"{_enso.loc['P_L', 'r_Nino12']:.2f}, con T PISCOt {_enso.loc['T_PISCOt', 'r_Nino12']:.2f}); retirar la "
                     f"señal ENSO no elimina el calentamiento. En años hidrológicos, el caudal dada la lluvia baja "
                     f"{_esc.loc['resid_R_dado_P', 'sen_dec']:+.0f} mm/año por década (p = {_esc.loc['resid_R_dado_P', 'p_mk']:.2f}), "
                     "algo compatible con más evapotranspiración o con usos del agua, pero que también depende de la crecida "
                     "de 1998 al inicio del registro. Es una hipótesis por contrastar, no una atribución.",
                     "Para atribuir faltan datos de uso del agua, de cambios de cobertura y metadatos de la estación.",
                     script="32_tendencias_3.py")
            + tabla(t3, "Tendencia de las anomalías, todos los meses",
                    "Pendiente por década con IC 95 %. OLS-HAC: Newey–West con 12 rezagos. MK Hamed–Rao + Sen: IC de Gilbert "
                    "con la varianza corregida. La serie original y las doce subseries están en tendencias_todos_meses_3_4.csv "
                    "y tendencias_por_mes_3_4.csv."),
            "Punto 3"),
    seccion("p4", "Escalas dominantes de variabilidad (Fourier)",
            "Espectros de la secuencia mensual cronológica (Δt = 1 mes) de cada variable, en su registro completo y en el "
            f"periodo común {PER}. Se comparan la serie centrada (conserva el ciclo anual) y las anomalías respecto a la "
            "climatología 1998–2020, con y sin la tendencia del punto 3. El estimador principal es el periodograma de la "
            "FFT sin relleno con ceros: f<sub>k</sub> = k/N, Δf = 1/N ciclos/mes, densidad unilateral normalizada para que "
            "Σ S·Δf sea la varianza (comprobado). Se contrasta con ventana de Hann y con Welch (segmentos de 10 años). "
            "Ninguna serie tiene vacíos internos; el caudal es continuo gracias al relleno documentado.",
            figura(ES / "fig_espectros_porcentaje_4_2.png", "Qué periodos dominan (% de la varianza)",
                   f"Es la lectura más directa del espectro. En las series centradas manda el ciclo anual "
                   f"({banda('P_L', 'centrada', 'anual 11–13 m'):.0f} % de la varianza de P<sub>L</sub> y "
                   f"{banda('Q', 'centrada', 'anual 11–13 m'):.0f} % del caudal en la banda de 11–13 meses), más el "
                   f"semianual por la estación húmeda concentrada. Al quitar la climatología ningún periodo domina: en "
                   f"la lluvia y el caudal el mayor pico interanual reúne apenas {PK['P_L'][1]:.1f} % "
                   f"(P<sub>L</sub>, ≈ {PK['P_L'][0]:.1f} años) y {PK['Q'][1]:.1f} % (caudal, ≈ {PK['Q'][0]:.1f} años), "
                   f"mientras que en la temperatura destaca un pico de ≈ {PK['T_PISCOt'][0]:.1f} años con "
                   f"{PK['T_PISCOt'][1]:.0f} % (PISCOt).",
                   "Los periodos largos se apoyan en pocos ciclos dentro del registro (un pico de 9 años cabe 3 veces en "
                   "28 años): por eso cada rótulo lleva el rango de periodos compatible con la resolución y el número de "
                   "ciclos. Un pico que supera el umbral puntual no es una periodicidad demostrada.",
                   script="33_espectros_4.py")
            + figura(ES / "fig_espectros_4_2.png", "Espectros del registro completo",
                   f"En las series centradas domina el ciclo anual en la lluvia (P<sub>L</sub> "
                   f"{banda('P_L', 'centrada', 'anual 11–13 m'):.0f} % de la varianza; IMERG "
                   f"{banda('P_I', 'centrada', 'anual 11–13 m'):.0f} %), con un armónico semianual que refleja la "
                   f"estación húmeda concentrada en pocos meses. En el caudal el ciclo anual pesa menos "
                   f"({banda('Q', 'centrada', 'anual 11–13 m'):.0f} %) porque las crecidas de El Niño (1998, 2017) reparten "
                   f"mucha varianza en todas las frecuencias. En la temperatura el ciclo anual es pequeño (PISCOt "
                   f"{banda('T_PISCOt', 'centrada', 'anual 11–13 m'):.0f} %) y dominan las escalas interanuales. Al retirar "
                   f"la climatología, la lluvia se parece a ruido casi blanco (AR(1): ρ = {RHO['P_L']:.2f} P<sub>L</sub>, "
                   f"{RHO['P_I']:.2f} IMERG); el caudal es más rojo (ρ = {RHO['Q']:.2f}) y la temperatura mucho más "
                   f"(ρ = {RHO['T_PISCOt']:.2f}–{RHO['T_ERA5L']:.2f}).",
                   f"Eje inferior en periodo (años), creciente hacia la derecha; eje superior en frecuencia (ciclos/mes): el "
                   f"ciclo anual está en 1/12 ≈ 0.083 ciclos/mes, el semianual en 1/6 y la frecuencia de Nyquist en 0.5 "
                   f"(periodo mínimo resoluble de 2 meses). Línea gris continua: fondo AR(1); discontinua y punteada: "
                   f"umbrales al 95 % puntual y global (corregido por buscar entre todas las frecuencias). "
                   + ("Ningún pico de las anomalías supera el umbral global" if N_GLOBAL == 0 else
                      f"Solo {N_GLOBAL} picos de las anomalías superan el umbral global")
                   + ": los picos interanuales no son periodicidades demostradas.",
                   script="33_espectros_4.py")
            + figura(ES / "fig_bandas_varianza_4_2.png", "¿En qué escalas está la varianza?",
                     f"Al quitar la climatología, la banda anual cae casi a cero en la lluvia, pero el caudal conserva "
                     f"{banda('Q', 'anomalia', 'anual 11–13 m'):.0f} % en esa banda: la amplitud del ciclo anual del caudal "
                     f"cambia mucho de un año a otro (años de El Niño). La banda interanual de 2–7 años reúne "
                     f"{banda('P_L', 'anomalia', 'interanual 2–7 años'):.0f} % de las anomalías de P<sub>L</sub> (1981–2025), "
                     f"{banda('Q', 'anomalia', 'interanual 2–7 años'):.0f} % del caudal y "
                     f"{banda('T_PISCOt', 'anomalia', 'interanual 2–7 años'):.0f} % de la temperatura PISCOt.",
                     "Fracciones del periodograma por bandas de periodo; cada barra suma 100 %. Registros de distinta "
                     "longitud tienen distinta resolución: compare las filas del mismo periodo.", script="33_espectros_4.py")
            + figura(ES / "fig_comparacion_espectros_4_2.png", "Fuentes y respuesta de la cuenca",
                     "Formas espectrales normalizadas por la varianza (las alturas no son variabilidad absoluta). CAMELS-PE e "
                     "IMERG tienen la misma forma en las series centradas; en anomalías IMERG tiene más varianza relativa en "
                     "escalas de 2–5 meses. PISCOt y ERA5-Land coinciden en casi todas las escalas. En las anomalías del "
                     f"periodo común, la escala de 2–5 meses reúne {banda('Q', 'anomalia', 'intraestacional 2–5 m', 'comun'):.0f} % "
                     f"de la varianza del caudal frente a {banda('P_L', 'anomalia', 'intraestacional 2–5 m', 'comun'):.0f} % de "
                     f"la lluvia, y la de 13 meses a 7 años {banda('Q', 'anomalia', '13 m–2 años', 'comun') + banda('Q', 'anomalia', 'interanual 2–7 años', 'comun'):.0f} % "
                     f"frente a {banda('P_L', 'anomalia', '13 m–2 años', 'comun') + banda('P_L', 'anomalia', 'interanual 2–7 años', 'comun'):.0f} %. "
                     "Es coherente con un almacenamiento que suaviza las variaciones rápidas y da persistencia.",
                     "Un espectro de potencia no contiene la fase: los desfases lluvia–caudal se estiman en el punto 2 "
                     "(rezagos), no aquí.", script="33_espectros_4.py")
            + figura(ES / "fig_estabilidad_4_2.png", "Estabilidad de los resultados",
                     "La banda anual y la interanual cambian poco con la tendencia o al recortar extremos, pero el periodo "
                     "del pico interanual cambia con la ventana y con la mitad del registro analizada (por ejemplo, P<sub>L</sub> "
                     "pasa de unos 5 a unos 2.4 años). Con unos 8–13 ciclos en el registro, los periodos interanuales solo se "
                     "pueden dar como bandas amplias.", script="33_espectros_4.py")
            + figura(ES / "fig_vacios_caudal_4_1.png", "Efecto de los vacíos del caudal",
                     f"Con solo los meses observados completos (Lomb–Scargle, estimador para muestreo irregular), la banda de "
                     f"2–5 meses pasa de {100 * _ls.loc['intraestacional 2–5 m', 'relleno_FFT']:.0f} % a "
                     f"{100 * _ls.loc['intraestacional 2–5 m', 'observado_LombScargle']:.0f} % de la varianza y bajan las bandas "
                     "largas. En la malla regular, la FFT y Lomb–Scargle dan el mismo resultado, así que la diferencia la "
                     "producen los vacíos, no el método: las conclusiones sobre escalas cortas del caudal son las más inciertas.",
                     script="33_espectros_4.py")
            + tabla(t4, "Picos principales de las anomalías (registro completo)",
                    "Rango = periodos que caben en ±½ bin de frecuencia (resolución Δf = 1/N). Ciclos = k, número de "
                    "oscilaciones completas en el registro. Todos los picos y transformaciones están en picos_4_2.csv; "
                    "N, Δf, ventana y normalización de cada espectro en metadatos_espectros_4_1.csv.")
            + '<h3 class="sub2">Complemento: onditas (en qué años aparece cada periodo)</h3><p class="intro">La guía '
              'pide espectros de Fourier, que promedian todo el registro. Como complemento, la transformada de onditas '
              'muestra si la variabilidad interanual está presente siempre o solo en algunos años, y si coincide en el '
              'tiempo con el ENSO. Es un análisis descriptivo: no atribuye la variabilidad de la cuenca al ENSO.</p>'
            + figura(ES / "fig_onditas_4_3.png", "Onditas: en qué años aparece cada periodo",
                     "La variabilidad fuera del ciclo anual no es constante en el tiempo. En la lluvia, la banda de 4–5 "
                     "años supera el ruido rojo entre 1983 y 2000 (El Niño 1983 y 1997–98) y después casi desaparece. El "
                     "caudal solo muestra pulsos cortos en 1998 y 2017. La temperatura mantiene una banda de 3–5 años "
                     "durante todo el registro, igual que Niño 1+2.",
                     "Lo que queda dentro del rayado (cono de influencia) está afectado por los bordes y no se interpreta. "
                     "El espectro global de la derecha es un Fourier suavizado: por eso sus picos coinciden con los del "
                     "periodograma.", script="38_onditas_4.py")
            + figura(ES / "fig_onditas_banda_4_3.png", "¿La banda de 2–7 años va con el ENSO?",
                     f"La componente de 2–7 años de la lluvia sigue a Niño 1+2 (r = {OND.loc['P_L', 'r_comp_2_7a_con_N12']:.2f}) "
                     f"y poco a Niño 3.4 (r = {OND.loc['P_L', 'r_comp_2_7a_con_N34']:.2f}): se parece más al ENSO costero "
                     f"que al del Pacífico central. El caudal la sigue menos (r = {OND.loc['Q', 'r_comp_2_7a_con_N12']:.2f}) "
                     f"y la temperatura mucho más (r = {OND.loc['T_ERA5L', 'r_comp_2_7a_con_N34']:.2f} con Niño 3.4). Pero "
                     f"esa banda es una parte pequeña de las anomalías de la cuenca: "
                     f"{100 * OND.loc['P_L', 'var_componente_2_7a']:.0f} % en la lluvia, "
                     f"{100 * OND.loc['Q', 'var_componente_2_7a']:.0f} % en el caudal y "
                     f"{100 * OND.loc['T_ERA5L', 'var_componente_2_7a']:.0f} % en la temperatura. Solo en 1997–98 la banda "
                     f"de la lluvia se acerca al umbral del 95 %.",
                     "Las correlaciones son descriptivas, sin valor p: una serie filtrada a 2–7 años tiene muy pocos grados "
                     "de libertad. La explicación física (El Niño costero, Pacífico central) debe sustentarse con "
                     "bibliografía.", script="38_onditas_4.py"),
            "Punto 4"),
    seccion("p5", "Mapas de correlación con el clima global",
            "Pendiente. Esta sección reunirá los mapas de correlación, mes a mes, entre la precipitación y el caudal de "
            "la cuenca y los campos globales de temperatura superficial del mar y de dos variables atmosféricas.",
            '<div class="tabla-cont"><h4>En construcción</h4><p>Aún no hay productos para este punto.</p></div>',
            "Punto 5"),
]

nav = [("cuenca", "Cuenca"), ("contexto", "Contexto físico"), ("p11", "1.1 Series"), ("p12", "1.2 IMERG"),
       ("p13", "1.3 Distribución"), ("p14", "1.4 Calidad"), ("p15a", "1.5a Climatología"),
       ("p15b", "1.5b Ciclo anual"), ("p15c", "1.5c Síntesis"), ("p2", "2. Relaciones"), ("p3", "3. Tendencias"),
       ("p4", "4. Fourier"), ("p5", "5. Clima global"), ("datos", "Datos")]

CSS = """
:root{--bg:#f6f5f2;--card:#ffffff;--ink:#1d1d1b;--ink2:#55534e;--line:#e2e0da;--acc:#1b6ca8;--chip:#e7f0f8;--code:#f0efeb;
  --muted:#898781;--axis:#c3c2b7;--pos:#2a78d6;--neg:#e34948;--excl:rgba(208,59,59,.12);
  --vPL:#2a78d6;--vPI:#eb6834;--vQ:#1baf7a;--vT:#c0392b;--vTE:#7d1f16;--apoyo:#898781}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){--bg:#151514;--card:#1f1f1d;--ink:#ecebe7;--ink2:#a9a7a1;--line:#34332f;--acc:#6aaee6;--chip:#1c2d3d;--code:#2a2a27;
  --muted:#898781;--axis:#383835;--pos:#3987e5;--neg:#e66767;--excl:rgba(230,103,103,.18);
  --vPL:#3987e5;--vPI:#d95926;--vQ:#199e70;--vT:#e66767;--vTE:#f0a39b;--apoyo:#a3a19b}}
:root[data-theme="dark"]{--bg:#151514;--card:#1f1f1d;--ink:#ecebe7;--ink2:#a9a7a1;--line:#34332f;--acc:#6aaee6;--chip:#1c2d3d;--code:#2a2a27;
  --muted:#898781;--axis:#383835;--pos:#3987e5;--neg:#e66767;--excl:rgba(230,103,103,.18);
  --vPL:#3987e5;--vPI:#d95926;--vQ:#199e70;--vT:#e66767;--vTE:#f0a39b;--apoyo:#a3a19b}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);font:16px/1.6 system-ui,-apple-system,"Segoe UI",Roboto,Ubuntu,sans-serif}
header.top{background:var(--card);border-bottom:1px solid var(--line)}
.wrap{max-width:1200px;margin:0 auto;padding:0 16px}
header.top .wrap{padding-top:28px;padding-bottom:22px}
h1{font-size:1.75rem;line-height:1.25;margin:0 0 6px}
.sub{color:var(--ink2);margin:0}
.datos-clave{display:grid;grid-template-columns:repeat(6,minmax(0,1fr));gap:10px;margin-top:16px}
.dato{background:var(--bg);border:1px solid var(--line);border-radius:8px;padding:8px 12px;font-size:.85rem;min-width:0}
.dato b{display:block;font-size:.98rem;line-height:1.3}
@media (max-width:900px){.datos-clave{grid-template-columns:repeat(3,minmax(0,1fr))}}
@media (max-width:520px){.datos-clave{grid-template-columns:repeat(2,minmax(0,1fr))}}
nav{position:sticky;top:0;z-index:5;background:var(--card);border-bottom:1px solid var(--line)}
nav .wrap{display:flex;gap:4px;overflow-x:auto;padding-top:8px;padding-bottom:8px}
nav a{white-space:nowrap;color:var(--ink2);text-decoration:none;padding:6px 10px;border-radius:6px;font-size:.92rem}
nav a:hover,nav a:focus-visible{background:var(--chip);color:var(--ink)}
nav a[aria-current="page"]{background:var(--chip);color:var(--acc);font-weight:600}
main .wrap{padding-top:8px;padding-bottom:40px}
section{padding-top:28px;display:none}
section.activa{display:block}
.pasos{display:flex;justify-content:space-between;gap:12px;margin:8px 0 0;flex-wrap:wrap}
.pasos a{color:var(--acc);text-decoration:none;border:1px solid var(--line);background:var(--card);border-radius:8px;padding:8px 14px;font-size:.92rem}
.pasos a:hover,.pasos a:focus-visible{background:var(--chip)}
.pasos .sig{margin-left:auto}
h2{font-size:1.4rem;margin:0 0 8px;display:flex;align-items:center;gap:10px;flex-wrap:wrap}
.chip{font-size:.75rem;font-weight:600;background:var(--chip);color:var(--acc);padding:3px 9px;border-radius:999px}
.intro{color:var(--ink2);margin:0 0 20px;text-align:justify;hyphens:auto}
/* Tarjeta de figura (diseño del visor de la compañera): título y tipo arriba; figura a todo el ancho de la
   columna (1200 px); debajo, "Qué muestra" y "Cómo se construyó" en dos columnas y "Cómo leerla" en una fila */
.fig{background:var(--card);border:1px solid var(--line);border-radius:12px;margin:0 0 24px;overflow:hidden}
.fig-h{display:flex;align-items:baseline;justify-content:space-between;gap:12px;padding:14px 18px 6px}
.fig-h h4{margin:0;font-size:1.08rem}
.tipo{font-size:.72rem;font-weight:600;letter-spacing:.02em;color:var(--ink2);border:1px solid var(--line);border-radius:999px;padding:1px 8px;white-space:nowrap}
.tipo.int{color:var(--acc);background:var(--chip);border-color:transparent}
.fig-cuerpo{padding:4px 12px 8px}
.zoom{all:unset;display:block;box-sizing:border-box;width:100%;cursor:zoom-in;background:#fff;border-radius:8px;padding:8px}
.zoom:focus-visible{outline:3px solid var(--acc)}
.zoom img{display:block;max-width:100%;max-height:80vh;width:auto;height:auto;margin:0 auto}
.plot{width:100%;min-height:420px}
.selector{display:flex;flex-wrap:wrap;gap:4px;padding:4px 6px 8px}
.fig>.selector{padding:4px 18px 8px}
.selector button{border:1px solid var(--line);background:var(--bg);color:var(--ink2);border-radius:6px;padding:3px 12px;font:inherit;font-size:.88rem;cursor:pointer}
.selector button[aria-pressed="true"]{background:var(--chip);color:var(--acc);border-color:var(--acc);font-weight:600}
.selector button:focus-visible{outline:2px solid var(--acc);outline-offset:1px}
.selector .sel-tit{align-self:center;color:var(--ink2);font-size:.85rem;margin-right:4px}
figcaption{padding:12px 18px 4px;border-top:1px solid var(--line)}
.ficha{display:grid;grid-template-columns:1fr 1fr;gap:6px 28px}
.ficha>div{min-width:0}
.ficha h5{margin:0 0 4px;font-size:.8rem;text-transform:uppercase;letter-spacing:.05em;color:var(--muted)}
.ficha p{margin:0 0 8px;font-size:.95rem}
.f-leer{grid-column:1 / -1;border-top:1px solid var(--line);padding-top:8px}
.leer{color:var(--ink2)}
.sub2{margin:40px 0 6px;padding-top:24px;border-top:1px solid var(--line);font-size:1.2rem}
/* control de calidad: gráfica y panel de la variable elegida lado a lado */
.qc-grid{display:grid;grid-template-columns:minmax(0,1.6fr) minmax(0,1fr);gap:4px 12px;align-items:start}
.qc-lado{padding:4px 18px 8px 4px;font-size:.92rem}
.qc-lado h5{margin:10px 0 4px;font-size:.8rem;text-transform:uppercase;letter-spacing:.05em;color:var(--muted)}
.qc-lado h5:first-child{margin-top:4px;font-size:1rem;text-transform:none;letter-spacing:0;color:var(--ink)}
.qc-lado p{margin:0 0 6px}
.qc-lado .datos-clave{display:flex;flex-wrap:wrap;margin-top:6px;gap:6px}.qc-lado .dato{padding:5px 9px;font-size:.78rem}.qc-lado .dato b{font-size:.95rem}
.qc-lista{margin:0;padding-left:0;list-style:none;max-height:300px;overflow-y:auto}
.qc-lista li{padding:6px 0;border-bottom:1px solid var(--line);font-size:.86rem}
.qc-lista i{color:var(--ink2)}
.estado{display:inline-block;font-size:.7rem;font-weight:600;border-radius:999px;padding:0 7px;border:1px solid var(--line);color:var(--ink2)}
.estado.e-corregido{color:var(--acc);border-color:var(--acc)}
.estado.e-incierto{color:var(--neg);border-color:var(--neg)}
@media (max-width:900px){.qc-grid{grid-template-columns:1fr}.qc-lado{padding:4px 12px 8px}}
@media (max-width:760px){
  .ficha{grid-template-columns:1fr}
  .fig-h{padding:12px 12px 4px}.fig-cuerpo{padding:4px 4px 8px}figcaption{padding:10px 12px 4px}
  .zoom{padding:4px}.plot{min-height:360px}
}
.params{display:grid;grid-template-columns:repeat(auto-fill,minmax(250px,1fr));gap:8px}
.param{background:var(--bg);border:1px solid var(--line);border-radius:8px;padding:7px 12px;font-size:.85rem;color:var(--ink2)}
.param b{display:block;color:var(--ink);font-size:1rem;font-variant-numeric:tabular-nums}
.meta{font-size:.8rem;color:var(--ink2);display:flex;flex-wrap:wrap;gap:10px;align-items:center;margin:0;padding:2px 0 12px}
.meta a{color:var(--acc)}
.script{font-family:ui-monospace,Consolas,monospace;background:var(--code);padding:1px 6px;border-radius:4px}
.nota{color:var(--ink2);font-size:.82rem}
.tabla-cont{background:var(--card);border:1px solid var(--line);border-radius:12px;padding:14px 18px;margin:0 0 22px}
.tabla-cont h4{margin:0 0 10px}
.scroll{overflow-x:auto}
table.tabla{border-collapse:collapse;font-size:.85rem;min-width:100%}
.tabla th,.tabla td{padding:6px 10px;border-bottom:1px solid var(--line);text-align:left;vertical-align:top}
.tabla th{background:var(--bg);position:sticky;top:0;font-weight:600}
.tabla td{font-variant-numeric:tabular-nums}
footer{border-top:1px solid var(--line);background:var(--card);color:var(--ink2);font-size:.88rem}
footer .wrap{padding-top:20px;padding-bottom:28px}
footer code{font-family:ui-monospace,Consolas,monospace;background:var(--code);padding:1px 5px;border-radius:4px}
#visor{position:fixed;inset:0;background:rgba(0,0,0,.88);display:none;align-items:center;justify-content:center;z-index:20;padding:16px;cursor:zoom-out}
#visor.abierto{display:flex}
#visor img{max-width:100%;max-height:100%;background:#fff;border-radius:6px}
#visor button{position:absolute;top:12px;right:16px;font-size:1.6rem;background:none;border:0;color:#fff;cursor:pointer}
.tema{margin-left:auto;border:1px solid var(--line);background:var(--bg);color:var(--ink);border-radius:6px;padding:4px 10px;cursor:pointer;font-size:.85rem}
@media (max-width:600px){h1{font-size:1.35rem}}
"""

JS = """
const secs=[...document.querySelectorAll('main section')],enl=[...document.querySelectorAll('nav a')];
const titulos=Object.fromEntries(enl.map(a=>[a.getAttribute('href').slice(1),a.textContent]));
secs.forEach((s,i)=>{const d=document.createElement('div');d.className='pasos';
  if(i>0){const a=document.createElement('a');a.href='#'+secs[i-1].id;a.textContent='← '+titulos[secs[i-1].id];d.appendChild(a);}
  if(i<secs.length-1){const a=document.createElement('a');a.className='sig';a.href='#'+secs[i+1].id;a.textContent=titulos[secs[i+1].id]+' →';d.appendChild(a);}
  s.appendChild(d);});
function mostrar(){const id=location.hash.slice(1);const s=secs.find(x=>x.id===id)||secs[0];
  secs.forEach(x=>x.classList.toggle('activa',x===s));
  enl.forEach(a=>{if(a.getAttribute('href')==='#'+s.id){a.setAttribute('aria-current','page');a.scrollIntoView({block:'nearest',inline:'nearest'});}else a.removeAttribute('aria-current');});
  document.title=titulos[s.id]+' · Cuenca El Tambo';window.scrollTo(0,0);
  if(window.Graficas)Graficas.dibujarEn(s);}
window.addEventListener('hashchange',mostrar);mostrar();
const visor=document.getElementById('visor'),vimg=visor.querySelector('img');
document.querySelectorAll('.zoom').forEach(b=>b.addEventListener('click',()=>{const i=b.querySelector('img');vimg.src=i.src;vimg.alt=i.alt;visor.classList.add('abierto');visor.querySelector('button').focus();}));
function cerrar(){visor.classList.remove('abierto');vimg.removeAttribute('src');}
visor.addEventListener('click',cerrar);
document.addEventListener('keydown',e=>{if(e.key==='Escape')cerrar();});
const raiz=document.documentElement,bt=document.querySelector('.tema');
function aplicar(t){if(t)raiz.setAttribute('data-theme',t);else raiz.removeAttribute('data-theme');}
try{aplicar(localStorage.getItem('tema'));}catch(e){}
bt.addEventListener('click',()=>{const oscuro=raiz.getAttribute('data-theme')==='dark'||(!raiz.getAttribute('data-theme')&&matchMedia('(prefers-color-scheme: dark)').matches);const t=oscuro?'light':'dark';aplicar(t);try{localStorage.setItem('tema',t);}catch(e){}
  if(window.Graficas)Graficas.redibujar();});
matchMedia('(prefers-color-scheme: dark)').addEventListener('change',()=>{if(window.Graficas)Graficas.redibujar();});
"""

# Plotly (incluido en el repositorio, sin CDN) + datos + gráficas: el visor sigue siendo un solo archivo sin conexión
PLOTLY = (Path(__file__).parent / "recursos" / "plotly-2.35.2.min.js").read_text(encoding="utf-8")
GRAFICAS_JS = (Path(__file__).parent / "visor_graficas.js").read_text(encoding="utf-8")
DATOS_JSON = json.dumps(DATOS, ensure_ascii=False, separators=(",", ":")).replace("</", "<\\/")

pagina = f"""<!doctype html>
<html lang="es">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Cuenca El Tambo – resultados</title>
<style>{CSS}</style>
</head>
<body>
<header class="top"><div class="wrap">
  <h1>Clasificación hidroclimática de la cuenca del río Chicama hasta El Tambo</h1>
  <p class="sub">Tarea 1 · Hidrología 2026-2S · Universidad Nacional de Colombia - Sede Medellín</p>
  <div class="datos-clave">
    <div class="dato">Estación<b>El Tambo · PE_47E0A2A8</b></div>
    <div class="dato">Área<b>2 170 km²</b></div>
    <div class="dato">Periodo común<b>{PER} ({N_MESES} meses)</b></div>
    <div class="dato">P<sub>L</sub> media anual<b>{cl[cl.variable == 'P_L'].media.sum():.0f} mm/año</b></div>
    <div class="dato">P<sub>I</sub> media anual<b>{cl[cl.variable == 'P_I'].media.sum():.0f} mm/año</b></div>
    <div class="dato">Régimen (propuesto)<b>{ip.loc['P_L', 'regimen']}, máximo en marzo</b></div>
  </div>
</div></header>
<nav aria-label="Secciones"><div class="wrap">
  {''.join(f'<a href="#{i}">{t}</a>' for i, t in nav)}
  <button class="tema" type="button" aria-label="Cambiar tema claro u oscuro">◐ Tema</button>
</div></nav>
<main><div class="wrap">
{''.join(secciones)}
<section id="datos"><h2>Datos y reproducibilidad</h2>
<p class="intro">Fuentes usadas y cómo regenerar todos los productos. Haga clic en cualquier figura para ampliarla (Esc para cerrar).</p>
{tabla(fuentes, "Fuentes de datos")}
<div class="tabla-cont"><h4>Cómo reproducir</h4>
<p>Con Python 3 y las dependencias de <code>requirements.txt</code>, ejecute en <code>Hidrologiat1/Codigos/</code>, en este orden:
<code>06 → 07a → 13 → 14 → 09 → 10 → 28 → 26 → 07b → 11 → 12 → 17 → 15 → 16 → 18 → 19 → 20 → 21 → 23 → 24 → 25 → 29 → 27 → 30 → 31 → 32 → 33 → 38 → 34a → 34b → 34c → 35 → 22</code>.
El periodo común se calcula en <code>periodo_comun.py</code>. Las descargas de IMERG y ERA5-Land requieren cuentas personales
(NASA Earthdata y Copernicus CDS), cuyas credenciales no se incluyen en el repositorio.</p></div>
</section>
</div></main>
<footer><div class="wrap">
Página generada el {date.today():%Y-%m-%d} con <code>Hidrologiat1/Codigos/22_visor_html.py</code> a partir de los archivos de
resultados. Archivo autocontenido: funciona sin conexión y sin instalar nada. Las descripciones distinguen lo que muestran los
datos de las hipótesis, que el grupo debe sustentar con bibliografía.
</div></footer>
<div id="visor" role="dialog" aria-modal="true" aria-label="Figura ampliada"><button type="button" aria-label="Cerrar">✕</button><img alt=""></div>
<script>{PLOTLY}</script>
<script>window.DATOS = {DATOS_JSON};</script>
<script>{GRAFICAS_JS}</script>
<script>{JS}</script>
</body>
</html>"""

SALIDA.write_text(pagina, encoding="utf-8")
print(f"Visor: {SALIDA} ({SALIDA.stat().st_size / 2**20:.1f} MB, {pagina.count('<figure')} figuras "
      f"({pagina.count('class=\"plot\"')} interactivas), {pagina.count('<table')} tablas)")
