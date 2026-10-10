"""
Punto 1.5.c (síntesis del ciclo anual): evidencia cuantitativa para proponer la clasificación
hidroclimática mensual de la cuenca El Tambo. La explicación física con bibliografía queda aparte.

La guía pide criterios explícitos, no copiar atributos de CAMELS, y que una categoría de aridez use un
criterio climático (no solo el total de lluvia). Se calculan:
  1. Aridez (UNEP; Middleton & Thomas 1997): AI = P/PET. < 0.05 hiperárido; 0.05–0.20 árido;
     0.20–0.50 semiárido; 0.50–0.65 subhúmedo seco; ≥ 0.65 húmedo. Se combinan 3 fuentes de P
     (P_L CAMELS-PE, P_I IMERG, PISCOp v3.0) con 2 de PET:
       - PET_PM: `pet` de CAMELS-PE (PISCOeo_pm v1.0, Penman–Monteith), disponible solo hasta 2016-12
         → periodo 1998-01..2016-12.
       - PET_HS: Hargreaves–Samani (FAO-56 ec. 52) con Tmin/Tmax diarias de PISCOt (CAMELS-PE) y Ra
         de la latitud del centroide de la cuenca → periodo 1998-01..2020-12. Es una alternativa de
         contraste: usa la temperatura media de cuenca (684–4317 m), no celda a celda.
     Balance de largo plazo frente a Budyko (1974): E/P = sqrt(φ·tanh(1/φ)·(1 − e^(−φ))), φ = PET/P;
     escorrentía esperada R/P = 1 − E/P frente a la observada (Σ de las 12 medias mensuales de R).
  2. Köppen–Geiger (criterios de Peel et al. 2007) con T media de cuenca y cada fuente de P
     (verano austral = oct–mar). Es una clase de cuenca agregada: el valle bajo y la puna difieren.
  3. Clasificación de cada mes calendario:
       - Humedad (criterio FAO de periodo de crecimiento): húmedo si P ≥ PET; intermedio si
         0.5·PET ≤ P < PET; seco si P < 0.5·PET (con PET_PM; mismos años 1998–2016).
       - Robustez: % de años (1998–2016) en que P_L ≥ 0.5·PET ese mes.
       - Seco de Bagnouls–Gaussen (P < 2T), periodo común.
       - Aguas altas: Q̄ medio del mes ≥ Q̄ medio anual, periodo común.
       - Frecuencia del mes como máximo anual de P_L y de Q̄ (años hidrológicos, de 1.5.b).
  4. Láminas por estación: P y R de la estación húmeda (dic–abr) y seca (may–nov), periodo común.
Las cifras del ciclo (SI, PCI, centroides, picos, subperiodos) vienen de 1.5.b y se resumen en la tabla
de criterios para tenerlas en un solo lugar.

Salidas (Tarea 1/El_Tambo/resultados/climatologia/):
  aridez_budyko_1_5c.csv, koppen_1_5c.csv, clasificacion_mensual_1_5c.csv, laminas_estacionales_1_5c.csv,
  criterios_clasificacion_1_5c.csv, fig_sintesis_ciclo_1_5c.png
Requiere 09, 10, 14, 17, 20 y 21.
"""
from pathlib import Path

import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

from periodo_comun import periodo_comun

ROOT = Path(__file__).resolve().parents[2] / "Tarea 1" / "El_Tambo"
SM = ROOT / "resultados" / "series_mensuales"
OUT = ROOT / "resultados" / "climatologia"

INI, FIN = periodo_comun()
FIN_PM = pd.Timestamp("2016-12-01")                        # último mes con PET de CAMELS-PE
FIN_HS = pd.Timestamp("2020-12-01")                        # último mes con Tmin/Tmax de PISCOt
MESES = ["Ene", "Feb", "Mar", "Abr", "May", "Jun", "Jul", "Ago", "Sep", "Oct", "Nov", "Dic"]
HUMEDA, SECA = [12, 1, 2, 3, 4], [5, 6, 7, 8, 9, 10, 11]
COL = {"P_L": "#1b6ca8", "P_I": "#e4572e", "R": "#5b3f8c", "PET": "#3d3d3d"}
CLASE_COL = {"húmedo": "#1b6ca8", "intermedio": "#9cc3e0", "seco": "#e8d9b5"}

# ------------------------------------------------------------------ datos mensuales
cam = pd.read_csv(SM / "series_mensuales_camels_pe.csv", index_col="fecha", parse_dates=True)
imerg = pd.read_csv(SM / "imerg_cuenca_mensual.csv", index_col="fecha", parse_dates=True)["P_I_mm_mes"]
pv3 = pd.read_csv(SM / "piscop_v3_cuenca_mensual.csv", index_col="fecha", parse_dates=True)["P_piscop_v3_mm_mes"]
text = pd.read_csv(SM / "temperatura_piscot_extendida_era5land.csv", index_col="fecha", parse_dates=True)
df = pd.DataFrame({"P_L": cam["P_L_mm_mes"], "P_I": imerg, "P_v3": pv3, "R": cam["R_mm_mes"], "Q": cam["Q_m3_s"],
                   "T": text["T_piscot_C"].fillna(text["T_estimada_era5land_C"])})

# PET diaria: Penman–Monteith de CAMELS-PE y Hargreaves–Samani desde Tmin/Tmax de PISCOt
dia = pd.read_csv(ROOT / "datos" / "camels_pe" / "PE_47E0A2A8.csv", parse_dates=["date"], index_col="date")
cuenca = gpd.read_file(ROOT / "resultados" / "delimitacion" / "cuenca_el_tambo_glo30.gpkg")
LAT = cuenca.to_crs(32717).centroid.to_crs(4326).y.iloc[0]
phi = np.deg2rad(LAT)
J = dia.index.dayofyear.values
dr = 1 + 0.033 * np.cos(2 * np.pi * J / 365)
dec = 0.409 * np.sin(2 * np.pi * J / 365 - 1.39)
ws = np.arccos(-np.tan(phi) * np.tan(dec))
Ra = 24 * 60 / np.pi * 0.0820 * dr * (ws * np.sin(phi) * np.sin(dec) + np.cos(phi) * np.cos(dec) * np.sin(ws))
pet_hs = 0.0023 * 0.408 * Ra * (dia["tmean"] + 17.8) * np.sqrt((dia["tmax"] - dia["tmin"]).clip(lower=0))


def mensual_completo(s):
    """Suma mensual solo de meses sin días faltantes."""
    g = s.resample("MS")
    return g.sum(min_count=1).where(g.count() == g.size())


df["PET_PM"] = mensual_completo(dia["pet"])
df["PET_HS"] = mensual_completo(pet_hs)


def clim(s):
    return s.groupby(s.index.month).mean().reindex(range(1, 13))


# ------------------------------------------------------------------ 1. aridez y Budyko
def budyko_R_P(ip):
    return 1 - np.sqrt(ip * np.tanh(1 / ip) * (1 - np.exp(-ip)))


def clase_unep(ai):
    for lim, nombre in [(0.05, "hiperárido"), (0.20, "árido"), (0.50, "semiárido"), (0.65, "subhúmedo seco")]:
        if ai < lim:
            return nombre
    return "húmedo"


filas = []
for pet, fin in (("PET_PM", FIN_PM), ("PET_HS", FIN_HS)):
    sub = df.loc[INI:fin]
    PET = clim(sub[pet]).sum()
    R = clim(sub["R"]).sum()
    for p, etq in (("P_L", "P_L (CAMELS-PE)"), ("P_I", "P_I (IMERG)"), ("P_v3", "PISCOp v3.0")):
        P = clim(sub[p]).sum()
        filas.append({"fuente_P": etq, "fuente_PET": "Penman–Monteith (PISCOeo_pm)" if pet == "PET_PM"
                      else "Hargreaves–Samani (PISCOt)", "periodo": f"{INI:%Y-%m}..{fin:%Y-%m}",
                      "P_mm_anio": round(P, 1), "PET_mm_anio": round(PET, 1), "AI_P_PET": round(P / PET, 3),
                      "clase_UNEP": clase_unep(P / PET), "PET_P": round(PET / P, 3),
                      "R_mm_anio": round(R, 1), "R_P_obs": round(R / P, 3),
                      "R_P_Budyko": round(budyko_R_P(PET / P), 3), "P_menos_R_mm": round(P - R, 1)})
aridez = pd.DataFrame(filas)
aridez.to_csv(OUT / "aridez_budyko_1_5c.csv", index=False, encoding="utf-8")


# ------------------------------------------------------------------ 2. Köppen–Geiger (Peel et al. 2007)
def koppen(P, T):
    """P y T: climatologías de 12 meses (índice 1–12). Hemisferio sur: verano = oct–mar."""
    ver, inv = [10, 11, 12, 1, 2, 3], [4, 5, 6, 7, 8, 9]
    MAP, MAT = P.sum(), T.mean()
    frac_ver = P[ver].sum() / MAP
    pthr = 2 * MAT + (28 if frac_ver >= 0.7 else 0 if frac_ver <= 0.3 else 14)
    d = {"MAP_mm": round(MAP, 1), "MAT_C": round(MAT, 2), "frac_P_verano": round(frac_ver, 3),
         "Pumbral_x10_mm": round(10 * pthr, 1), "Thot_C": round(T.max(), 2), "Tcold_C": round(T.min(), 2),
         "P_seco_invierno_mm": round(P[inv].min(), 1), "P_humedo_verano_mm": round(P[ver].max(), 1)}
    if MAP < 10 * pthr:
        clase = ("BW" if MAP < 5 * pthr else "BS") + ("h" if MAT >= 18 else "k")
    elif T.min() >= 18:
        clase = "A"
    else:
        ps_dry, pw_dry, ps_wet, pw_wet = P[ver].min(), P[inv].min(), P[ver].max(), P[inv].max()
        if ps_dry < 40 and ps_dry < pw_wet / 3:
            seg = "s"
        elif pw_dry < ps_wet / 10:
            seg = "w"
        else:
            seg = "f"
        ter = "a" if T.max() >= 22 else ("b" if (T >= 10).sum() >= 4 else "c")
        clase = "C" + seg + ter
    d["clase"] = clase
    return d


comun = df.loc[INI:FIN]
Tc = clim(comun["T"])
kop = pd.DataFrame([{"fuente_P": etq, "periodo": f"{INI:%Y-%m}..{FIN:%Y-%m}", **koppen(clim(comun[p]), Tc)}
                    for p, etq in (("P_L", "P_L (CAMELS-PE)"), ("P_I", "P_I (IMERG)"), ("P_v3", "PISCOp v3.0"))])
kop.to_csv(OUT / "koppen_1_5c.csv", index=False, encoding="utf-8")


# ------------------------------------------------------------------ 3. clasificación de cada mes
def clase_fao(p, pet):
    return "húmedo" if p >= pet else ("intermedio" if p >= 0.5 * pet else "seco")


pm = df.loc[INI:FIN_PM]
cP_L, cP_I, cPET = clim(pm["P_L"]), clim(pm["P_I"]), clim(pm["PET_PM"])
ok = pm[["P_L", "PET_PM"]].dropna()
pct_moist = (ok["P_L"] >= 0.5 * ok["PET_PM"]).groupby(ok.index.month).mean().reindex(range(1, 13)) * 100
picos = pd.read_csv(OUT / "mes_pico_por_anio_1_5b.csv")
n_pico_P = picos["mes_max_P_L"].value_counts().reindex(MESES, fill_value=0).values
n_pico_Q = picos["mes_max_Q"].value_counts().reindex(MESES, fill_value=0).values
n_anios_Q = picos["mes_max_Q"].notna().sum()
cQ = clim(comun["Q"])
cPLc, cPIc, cR = clim(comun["P_L"]), clim(comun["P_I"]), clim(comun["R"])

mens = pd.DataFrame({
    "mes": MESES,
    "P_L_mm": cPLc.values.round(1), "P_I_mm": cPIc.values.round(1), "R_mm": cR.values.round(1),
    "Q_m3_s": cQ.values.round(2), "T_C": Tc.values.round(2),
    "PET_PM_mm_1998_2016": cPET.values.round(1), "PET_HS_mm_1998_2020": clim(df.loc[INI:FIN_HS, "PET_HS"]).values.round(1),
    "P_L_PET_1998_2016": (cP_L / cPET).values.round(2), "P_I_PET_1998_2016": (cP_I / cPET).values.round(2),
    "clase_FAO_P_L": [clase_fao(a, b) for a, b in zip(cP_L, cPET)],
    "clase_FAO_P_I": [clase_fao(a, b) for a, b in zip(cP_I, cPET)],
    "pct_anios_P_L_ge_media_PET": pct_moist.values.round(0),
    "seco_Gaussen_P_L": (cPLc < 2 * Tc).values, "seco_Gaussen_P_I": (cPIc < 2 * Tc).values,
    "aguas_altas_Q_ge_media": (cQ >= cQ.mean()).values,
    "R_mayor_que_P_L": (cR > cPLc).values,
    "n_anios_pico_P_L": n_pico_P, "n_anios_pico_Q": n_pico_Q,
})
mens.to_csv(OUT / "clasificacion_mensual_1_5c.csv", index=False, encoding="utf-8")

# ------------------------------------------------------------------ 4. láminas por estación
est = []
for nombre, meses in (("húmeda dic–abr", HUMEDA), ("seca may–nov", SECA)):
    idx = [m - 1 for m in meses]
    P, PI, R = cPLc.values[idx].sum(), cPIc.values[idx].sum(), cR.values[idx].sum()
    est.append({"estacion": nombre, "P_L_mm": round(P, 1), "P_I_mm": round(PI, 1), "R_mm": round(R, 1),
                "pct_P_L_anual": round(100 * P / cPLc.sum(), 1), "pct_R_anual": round(100 * R / cR.sum(), 1),
                "R_P_L": round(R / P, 2), "R_P_I": round(R / PI, 2)})
estac = pd.DataFrame(est)
estac.to_csv(OUT / "laminas_estacionales_1_5c.csv", index=False, encoding="utf-8")

# ------------------------------------------------------------------ 5. tabla de criterios
ind = pd.read_csv(OUT / "indicadores_estacionalidad_1_5b.csv")
pc = ind[ind["muestra"].str.startswith("periodo")].set_index("variable")
sp = ind[ind["muestra"].str.startswith("subperiodo")]


def walsh(si):
    for lim, nombre in [(0.2, "muy uniforme"), (0.4, "uniforme con estación más húmeda"),
                        (0.6, "algo estacional"), (0.8, "estacional"),
                        (1.0, "marcadamente estacional con estación seca larga"), (1.2, "casi toda en ≤ 3 meses")]:
        if si < lim:
            return nombre
    return "extremo"


def oliver(pci):
    return "uniforme" if pci <= 10 else "moderada" if pci <= 15 else "estacional" if pci <= 20 else "fuerte"


def rango_sub(v, col, fmt="{:.2f}"):
    s = sp[sp["variable"] == v][col]
    return f"{fmt.format(s.min())}–{fmt.format(s.max())}"


def meses_de(cond):
    return ", ".join(m for m, c in zip(MESES, cond) if c) or "ninguno"


pm_L = aridez[(aridez["fuente_PET"].str.startswith("Penman"))].set_index("fuente_P")
hs_L = aridez[(aridez["fuente_PET"].str.startswith("Harg"))].set_index("fuente_P")
pctP = (100 * np.array(n_pico_P) / len(picos))
sdT = comun["T"].groupby(comun.index.month).std()
crit = [
    ("Número de picos", "picos del ciclo (prominencia ≥ 10 % amplitud) y armónicos (H1 ≥ 2·H2 → unimodal)",
     f"P_L: {int(pc.loc['P_L','n_picos'])} pico ({pc.loc['P_L','picos']}); H1 {pc.loc['P_L','var_H1_pct']:.0f} %, H2 {pc.loc['P_L','var_H2_pct']:.0f} %",
     f"P_I: {int(pc.loc['P_I','n_picos'])} pico ({pc.loc['P_I','picos']}); R: {int(pc.loc['R','n_picos'])} pico ({pc.loc['R','picos']})",
     "unimodal", "igual en ambos subperiodos", "fig_ciclo_1_5b, indicadores_estacionalidad_1_5b"),
    ("Mes del máximo de P", "mes de mayor media; % de años hidrológicos con el máximo en feb–abr",
     f"P_L: {pc.loc['P_L','mes_max']} ({pc.loc['P_L','max']:.0f} mm); máximo en feb–abr en {pctP[1:4].sum():.0f} % de los años",
     f"P_I: {pc.loc['P_I','mes_max']} ({pc.loc['P_I','max']:.0f} mm)", "máximo en marzo",
     "máximo fuera de feb–abr: " + ", ".join(f"{a} ({m})" for a, m in picos.loc[~picos["mes_max_P_L"].isin(
         ["Feb", "Mar", "Abr"]), ["anio_hidrologico", "mes_max_P_L"]].values),
     "mes_pico_por_anio_1_5b, fig_anio_mes_1_5a"),
    ("Concentración estacional", "SI (Walsh & Lawler 1981), PCI (Oliver 1980), % de P en dic–abr",
     f"P_L: SI {pc.loc['P_L','SI']:.2f} ({walsh(pc.loc['P_L','SI'])}); PCI {pc.loc['P_L','PCI']:.1f} ({oliver(pc.loc['P_L','PCI'])}); {pc.loc['P_L','pct_en_dic_abr']:.0f} % en dic–abr",
     f"P_I: SI {pc.loc['P_I','SI']:.2f} ({walsh(pc.loc['P_I','SI'])}); PCI {pc.loc['P_I','PCI']:.1f} ({oliver(pc.loc['P_I','PCI'])}); {pc.loc['P_I','pct_en_dic_abr']:.0f} % en dic–abr",
     "estacional a marcadamente estacional",
     f"SI P_L en subperiodos {rango_sub('P_L', 'SI')} y PCI {rango_sub('P_L', 'PCI', '{:.1f}')}: en los límites "
     "entre clases (SI 0.80; PCI 15–16)",
     "indicadores_estacionalidad_1_5b"),
    ("Estación húmeda y seca (lluvia)", "meses ≥ media mensual; Gaussen P < 2T",
     f"P_L húmedos: {pc.loc['P_L','meses_humedos']}; secos (Gaussen): {pc.loc['P_L','meses_secos_gaussen']}",
     f"P_I: mismos meses húmedos y secos", "húmeda dic–abr; seca may–sep; oct–nov transición",
     "en 2012–2025 P_I suma nov como seco", "clasificacion_mensual_1_5c"),
    ("Balance de humedad mensual", "FAO: húmedo P ≥ PET; intermedio 0.5·PET ≤ P < PET; seco P < 0.5·PET (PET_PM 1998–2016)",
     f"P_L húmedos: {meses_de(mens['clase_FAO_P_L'] == 'húmedo')}; intermedios: {meses_de(mens['clase_FAO_P_L'] == 'intermedio')}",
     f"P_I húmedos: {meses_de(mens['clase_FAO_P_I'] == 'húmedo')}; intermedios: {meses_de(mens['clase_FAO_P_I'] == 'intermedio')}",
     "déficit hídrico climático la mayor parte del año",
     "el número de meses húmedos depende de la fuente de P", "clasificacion_mensual_1_5c, fig_sintesis_ciclo_1_5c"),
    ("Aridez anual", "UNEP AI = P/PET (2 fuentes de PET)",
     f"P_L: AI {pm_L.loc['P_L (CAMELS-PE)','AI_P_PET']:.2f} ({pm_L.loc['P_L (CAMELS-PE)','clase_UNEP']}) con PM; {hs_L.loc['P_L (CAMELS-PE)','AI_P_PET']:.2f} ({hs_L.loc['P_L (CAMELS-PE)','clase_UNEP']}) con HS",
     f"P_I: AI {pm_L.loc['P_I (IMERG)','AI_P_PET']:.2f} ({pm_L.loc['P_I (IMERG)','clase_UNEP']}); PISCOp v3: {pm_L.loc['PISCOp v3.0','AI_P_PET']:.2f} ({pm_L.loc['PISCOp v3.0','clase_UNEP']})",
     "subhúmedo seco a semiárido", "la clase cambia con la fuente de P y de PET", "aridez_budyko_1_5c"),
    ("Clima (Köppen–Geiger)", "Peel et al. (2007) con T y P de cuenca",
     f"P_L: {kop.loc[0,'clase']} (MAP {kop.loc[0,'MAP_mm']:.0f} vs umbral B {kop.loc[0,'Pumbral_x10_mm']:.0f} mm)",
     f"P_I: {kop.loc[1,'clase']}; PISCOp v3: {kop.loc[2,'clase']}", f"{kop.loc[0,'clase']} con P_L",
     "frontera C/B según la fuente de P; T media de cuenca oculta el gradiente 684–4317 m", "koppen_1_5c"),
    ("Respuesta del caudal", "centroide de Markham R − P; mes de máximo Q̄; aguas altas Q̄ ≥ media",
     f"desfase R−P_L {pc.loc['R','centroide_mes'] - pc.loc['P_L','centroide_mes']:.1f} meses; máximo Q̄ {MESES[int(cQ.idxmax()) - 1]}; aguas altas: {meses_de(mens['aguas_altas_Q_ge_media'])}",
     f"desfase R−P_I {pc.loc['R','centroide_mes'] - pc.loc['P_I','centroide_mes']:.1f} meses; máximo Q̄ en abr en {n_pico_Q[3]} de {n_anios_Q} años",
     "caudal unimodal desfasado ~1–1.5 meses", f"SI de R {pc.loc['R','SI']:.2f} > SI de P_L: caudal más concentrado que la lluvia",
     "fig_ciclo_1_5b, clasificacion_mensual_1_5c"),
    ("Partición estacional", "láminas P y R por estación (periodo común)",
     f"seca may–nov: P_L {estac.loc[1,'P_L_mm']:.0f} mm, R {estac.loc[1,'R_mm']:.0f} mm (R/P_L {estac.loc[1,'R_P_L']:.2f})",
     f"húmeda dic–abr: P_L {estac.loc[0,'P_L_mm']:.0f} mm, R {estac.loc[0,'R_mm']:.0f} mm (R/P_L {estac.loc[0,'R_P_L']:.2f})",
     "caudal de estación seca sostenido por almacenamiento", f"R > P_L en la climatología de: {meses_de(mens['R_mayor_que_P_L'])}",
     "laminas_estacionales_1_5c"),
    ("Coherencia del balance", "R/P observado vs Budyko (1974)",
     f"P_L: R/P {pm_L.loc['P_L (CAMELS-PE)','R_P_obs']:.2f} vs Budyko {pm_L.loc['P_L (CAMELS-PE)','R_P_Budyko']:.2f}",
     f"P_I: R/P {pm_L.loc['P_I (IMERG)','R_P_obs']:.2f} vs Budyko {pm_L.loc['P_I (IMERG)','R_P_Budyko']:.2f}",
     "R observado mayor que el esperado con cualquier P", "P subestimada, PET sobrestimada o aportes no representados",
     "aridez_budyko_1_5c"),
    ("Temperatura", "amplitud del ciclo de 12 medias",
     f"T: amplitud {pc.loc['T','amplitud']:.2f} °C (máx {pc.loc['T','mes_max']}, mín {pc.loc['T','mes_min']})",
     f"desv. est. entre años de cada mes {sdT.min():.2f}–{sdT.max():.2f} °C, del orden de la amplitud del ciclo",
     "isotermal (régimen térmico tropical)",
     "T de dos fuentes (PISCOt y estimación ERA5-Land desde 2021)", "climatologia_1_5a"),
]
crit = pd.DataFrame(crit, columns=["rasgo", "criterio", "evidencia_principal", "evidencia_alternativa",
                                   "lectura_propuesta", "robustez_o_incertidumbre", "fuente"])
crit.to_csv(OUT / "criterios_clasificacion_1_5c.csv", index=False, encoding="utf-8")

# ------------------------------------------------------------------ figura de síntesis
x = np.arange(12)
fig = plt.figure(figsize=(12, 8.2))
gs = fig.add_gridspec(2, 1, height_ratios=[3, 1.55], hspace=0.08)
ax = fig.add_subplot(gs[0])
w = 0.38
ax.bar(x - w / 2, cPLc, w, color=COL["P_L"], label="P$_L$ (CAMELS-PE)", zorder=3)
ax.bar(x + w / 2, cPIc, w, color=COL["P_I"], label="P$_I$ (IMERG)", zorder=3)
ax.plot(x, cPET, color=COL["PET"], lw=2, marker="o", ms=5, label="PET Penman–Monteith (1998–2016)", zorder=4)
ax.plot(x, 0.5 * cPET, color=COL["PET"], lw=1.4, ls="--", label="0.5 · PET (umbral FAO)", zorder=4)
ax.plot(x, cR, color=COL["R"], lw=2, marker="s", ms=5, label="R (caudal en lámina)", zorder=5)
ax.axvspan(-0.5, 3.5, color="#1b6ca8", alpha=0.06, lw=0)
ax.axvspan(10.5, 11.5, color="#1b6ca8", alpha=0.06, lw=0)
ax.text(1.5, 0.97, "estación húmeda (dic–abr)", transform=ax.get_xaxis_transform(), ha="center", va="top",
        fontsize=9, color="#333")
ax.text(7, 0.97, "estación seca (may–nov)", transform=ax.get_xaxis_transform(), ha="center", va="top",
        fontsize=9, color="#333")
ax.set_ylabel("mm/mes")
ax.set_xlim(-0.5, 11.5)
ax.set_ylim(0, max(cPLc.max(), cPET.max()) * 1.18)
ax.tick_params(axis="x", labelbottom=False)
ax.grid(axis="y", color="#ddd", lw=0.6, zorder=0)
ax.spines[["top", "right"]].set_visible(False)
ax.legend(loc="upper right", bbox_to_anchor=(1, 0.9), fontsize=8.5, frameon=False, ncol=1)
ax.set_title(f"Síntesis del ciclo anual – cuenca El Tambo (P, R: periodo común {INI:%Y-%m} a {FIN:%Y-%m})",
             fontsize=12, loc="left")

# franjas de clasificación
bx = fig.add_subplot(gs[1], sharex=ax)
filas_f = [("P$_L$ vs PET (FAO)", mens["clase_FAO_P_L"]), ("P$_I$ vs PET (FAO)", mens["clase_FAO_P_I"])]
for i, (nombre, clases) in enumerate(filas_f):
    y = 5 - i
    for j, c in enumerate(clases):
        bx.add_patch(plt.Rectangle((j - 0.5, y - 0.42), 1, 0.84, color=CLASE_COL[c], lw=0))
    bx.text(-0.62, y, nombre, ha="right", va="center", fontsize=8.5)
    if i == 0:
        for j, v in enumerate(mens["pct_anios_P_L_ge_media_PET"]):
            bx.text(j, y, f"{v:.0f}%", ha="center", va="center", fontsize=7.5,
                    color="white" if mens["clase_FAO_P_L"][j] == "húmedo" else "#222")
otras = [("Seco Gaussen (P$_L$ < 2T)", mens["seco_Gaussen_P_L"], "#c9a227"),
         ("Aguas altas (Q̄ ≥ media)", mens["aguas_altas_Q_ge_media"], COL["R"]),
         ("R > P$_L$", mens["R_mayor_que_P_L"], "#7f7f7f")]
for i, (nombre, cond, color) in enumerate(otras):
    y = 3 - i
    for j, c in enumerate(cond):
        if c:
            bx.add_patch(plt.Rectangle((j - 0.5, y - 0.42), 1, 0.84, color=color, alpha=0.85, lw=0))
    bx.text(-0.62, y, nombre, ha="right", va="center", fontsize=8.5)
y = 0
bx.text(-0.62, y, "Años con pico en el mes\n(P$_L$ / Q̄)", ha="right", va="center", fontsize=8.5)
for j in range(12):
    bx.text(j, y, f"{n_pico_P[j]} / {n_pico_Q[j]}", ha="center", va="center", fontsize=8,
            color="#222" if n_pico_P[j] + n_pico_Q[j] else "#aaa")
for j in range(13):
    bx.axvline(j - 0.5, color="white", lw=1.5)
bx.set_ylim(-0.6, 5.6)
bx.set_yticks([])
bx.set_xticks(x, MESES)
for s in ("top", "right", "left"):
    bx.spines[s].set_visible(False)
leg = [Patch(color=CLASE_COL[k], label=f"{k}  ({d})") for k, d in
       (("húmedo", "P ≥ PET"), ("intermedio", "P ≥ 0.5·PET"), ("seco", "P < 0.5·PET"))]
bx.legend(handles=leg, title="Clase FAO", loc="upper left",
          bbox_to_anchor=(1.0, 1.02), fontsize=8, title_fontsize=8, frameon=False)
a = pm_L
nota = (f"Aridez UNEP (P/PET, PM 1998–2016): P$_L$ {a.loc['P_L (CAMELS-PE)','AI_P_PET']:.2f} "
        f"({a.loc['P_L (CAMELS-PE)','clase_UNEP']}), PISCOp v3 {a.loc['PISCOp v3.0','AI_P_PET']:.2f}, "
        f"P$_I$ {a.loc['P_I (IMERG)','AI_P_PET']:.2f} ({a.loc['P_I (IMERG)','clase_UNEP']}).  "
        f"Köppen: {kop.loc[0,'clase']} (P$_L$) / {kop.loc[1,'clase']} (P$_I$).  "
        f"% en la fila FAO de P$_L$: años con P$_L$ ≥ 0.5·PET.  Picos: {len(picos)} años hidrológicos "
        f"(Q̄: {n_anios_Q} con dato).")
fig.text(0.125, 0.02, nota, fontsize=8, color="#444", wrap=True)
fig.subplots_adjust(left=0.2, right=0.84, top=0.94, bottom=0.11)
fig.savefig(OUT / "fig_sintesis_ciclo_1_5c.png", dpi=200)
plt.close(fig)

# ------------------------------------------------------------------ resumen en consola
pd.set_option("display.width", 200)
print(f"Latitud del centroide: {LAT:.3f}°")
print("\nAridez y Budyko:\n", aridez.to_string(index=False))
print("\nKöppen:\n", kop.to_string(index=False))
print("\nClasificación mensual:\n", mens.to_string(index=False))
print("\nLáminas estacionales:\n", estac.to_string(index=False))
print("\nCriterios:")
for _, r in crit.iterrows():
    print(f"- {r['rasgo']}: {r['evidencia_principal']} | {r['evidencia_alternativa']} -> {r['lectura_propuesta']} "
          f"[{r['robustez_o_incertidumbre']}]")
