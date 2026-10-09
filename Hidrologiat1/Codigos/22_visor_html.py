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

import numpy as np
import pandas as pd
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


def figura(p, titulo, que, leer_=None, script=None):
    p = Path(p)
    extra = f'<p class="leer"><b>Cómo leerla.</b> {leer_}</p>' if leer_ else ""
    pie = f'<span class="script">{html.escape(script)}</span>' if script else ""
    return f"""
<figure class="fig">
  <button class="zoom" aria-label="Ampliar figura: {html.escape(titulo)}">
    <img loading="lazy" src="{img64(p)}" alt="{html.escape(titulo)}">
  </button>
  <figcaption>
    <h4>{titulo}</h4>
    <p>{que}</p>{extra}
    <p class="meta">{pie}<a href="{rel(p)}">PNG original</a> <span class="nota">(disponible si se clona el repositorio)</span></p>
  </figcaption>
</figure>"""


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
t14 = reg[["id", "variable", "fecha_o_periodo", "anomalia", "decision", "estado"]]
t14.columns = ["Id", "Variable", "Fecha / periodo", "Anomalía", "Decisión", "Estado"]

secciones = [
    seccion("cuenca", "La cuenca y su delimitación",
            "Río Chicama hasta la estación hidrométrica El Tambo (SENAMHI 211003; CAMELS-PE PE_47E0A2A8), en la "
            "vertiente del Pacífico del norte del Perú (−7.574°, −78.707°; 684 m s.n.m.). La cuenca se delimitó "
            "con el modelo digital de elevación Copernicus GLO-30 y se comparó con el polígono oficial de CAMELS-PE.",
            figura(DL / "mapa_cuenca_el_tambo.png", "Mapa de la cuenca",
                   "Relieve (GLO-30), límite de la cuenca, río Chicama, red de drenaje, punto de cierre y estaciones de SENAMHI. "
                   "Área delimitada 2 170 km² (CAMELS-PE: 2 171.6 km²), con cotas entre unos 680 y 4 300 m s.n.m. "
                   "Las flechas indican la <i>hipótesis</i> de fuentes de humedad que se contrastará en el punto 5.",
                   script="06_delimitacion_el_tambo_glo30.py")
            + figura(DL / "comparacion_camels.png", "Delimitación propia frente a CAMELS-PE",
                     "Los dos polígonos casi coinciden: índice de intersección sobre unión (IoU) = 0.994. Se usa la "
                     "delimitación propia para ponderar las mallas (IMERG, PISCOp, ERA5-Land).",
                     script="06_delimitacion_el_tambo_glo30.py")),
    seccion("contexto", "Contexto físico",
            "Rasgos del terreno que ayudan a explicar cómo la cuenca recibe, guarda y entrega el agua.",
            figura(CF / "panel_contexto_fisico.png", "Panel de contexto físico",
                   "(a) pendiente, (b) orientación de laderas, (c) cobertura del suelo, (d) unidades hidrogeológicas "
                   "agrupadas desde la geología de INGEMMET y (e) precipitación media anual en malla.",
                   "La agrupación hidrogeológica es una propuesta del grupo (por porosidad y permeabilidad esperadas) "
                   "que debe contrastarse con INGEMMET y ANA.", script="07b_mapas_contexto_fisico.py")
            + figura(CF / "mapa_e_precipitacion.png", "Precipitación media anual (PISCOp v3.0, fuente de apoyo)",
                     f"Media anual del periodo común en las celdas de 0.1° de PISCOp v3.0, sin interpolar. Aumenta desde la "
                     f"salida hacia las cabeceras del norte y el sureste. CAMELS-PE, la fuente principal de P<sub>L</sub>, solo "
                     f"entrega el promedio de cuenca y no tiene malla; por eso el mapa usa PISCOp v3.0, que da un promedio "
                     f"menor que CAMELS-PE.", script="07b_mapas_contexto_fisico.py")
            + figura(CF / "curva_hipsometrica.png", "Curva hipsométrica",
                     "Distribución del área por elevación. La integral hipsométrica (0.52) indica una cuenca de relieve "
                     "intermedio, con buena parte del área a media ladera.", script="07b_mapas_contexto_fisico.py")),
    seccion("p11", "Series mensuales de precipitación, caudal y temperatura",
            f"Series construidas desde los registros diarios: P<sub>L</sub> = suma diaria; Q̄ = media diaria (m³/s); "
            f"R = caudal en lámina (mm/mes); T = media diaria. Todas las figuras usan el <b>periodo común {PER}</b> "
            f"({N_MESES} meses): los meses con precipitación de CAMELS-PE y de IMERG a la vez.",
            figura(SM / "fig_series_mensuales_1_1.png", "Series mensuales en paneles alineados",
                   "Los cuatro paneles comparten el eje de tiempo. Los meses sin dato quedan como vacíos sin interpolar "
                   "(sombreado gris). Los círculos naranjas marcan meses de caudal con hasta 10 % de días faltantes. "
                   "Desde 2021 la temperatura es una estimación con ERA5-Land (línea discontinua morada), no PISCOt.",
                   "Los máximos de 1998 y 2017 corresponden a eventos de El Niño; el caudal queda alto también en los "
                   "meses siguientes.", script="15_figura_series_1_1.py"), "Punto 1.1"),
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
            + figura(SM / "fig_P_L_vs_IMERG_mensual.png", "P<sub>L</sub> frente a IMERG, mes a mes",
                     "IMERG sigue bien el ritmo de la lluvia, pero sus picos son mucho más bajos, sobre todo en los "
                     "meses más lluviosos.", script="10_imerg_mensual_y_grafica.py")
            + figura(SM / "fig_precipitacion_camels_1998_2025.png", "P<sub>L</sub> de CAMELS-PE",
                     "(a) Acumulado mensual y media móvil de 12 meses; (b) totales anuales (2025 queda fuera por estar "
                     "incompleto).", script="11_grafica_precipitacion_camels.py")
            + figura(SM / "fig_precipitacion_imerg_1998_2025.png", "P<sub>I</sub> de IMERG",
                     "Mismo formato y mismos ejes que la figura anterior, para compararlas a simple vista.",
                     script="16_grafica_precipitacion_imerg.py"), "Punto 1.2"),
    seccion("temp", "Temperatura: PISCOt frente a ERA5-Land",
            f"PISCOt (la temperatura de CAMELS-PE) termina en 2020. Se descargó ERA5-Land para compararlo y para "
            f"estimar 2021 en adelante. En los {TB['n']} meses con ambas fuentes varían igual (r = {TB['r']:.2f}), pero "
            f"ERA5-Land es {abs(TB['sesgo']):.1f} °C más fría de forma casi constante. La estimación desde 2021 "
            f"(ERA5-Land + sesgo de cada mes) tuvo un error medio de {float(val['MAE_C']):.2f} °C al validarla en "
            f"{val['evaluacion']}. Las dos fuentes no se unen en una sola serie.",
            figura(SM / "fig_T_camels_vs_era5land.png", "Temperatura media mensual: CAMELS-PE vs ERA5-Land",
                   "(a) Series superpuestas; (b) ciclo anual medio. Solo unos 0.5 °C del desfase se explican por la "
                   "altitud del relieve del modelo. Otra causa probable es la definición: PISCOt usa (Tmáx + Tmín)/2, "
                   "y ERA5-Land promedia las 24 horas.", script="17_era5land_temperatura.py")),
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
            figura(CC / "fig_disponibilidad_1_4.png", "Disponibilidad año × mes",
                   "Condición de cada mes y variable. Los faltantes de caudal se concentran en septiembre (y agosto–octubre) "
                   "casi todos los años, un patrón sistemático, además de un periodo con más vacíos en 2019–2021.",
                   script="19_control_calidad_1_4.py")
            + tabla(t14, "Registro de anomalías",
                    "Estado: «sin problema» = comprobado, nada que corregir; «incierto» = la duda se declara como "
                    "limitación. El detalle de cada comprobación está en registro_anomalias_1_4.csv."),
            "Punto 1.4"),
    seccion("p15", "Ciclo anual y su variabilidad",
            f"Climatología de 12 meses en el periodo común y rasgos del ciclo. P<sub>L</sub> llega a su máximo en marzo "
            f"({cmes('P_L', 3)['media']:.0f} mm de media) y a su mínimo en julio–agosto; el caudal llega a su máximo en "
            f"abril ({cmes('Q', 4)['media']:.0f} m³/s). Con los criterios propuestos, el régimen es "
            f"<b>{ip.loc['P_L', 'regimen']}</b>.",
            figura(CL / "fig_climatologia_1_5a.png", "Ciclo anual medio con bandas de dispersión",
                   "Línea continua = media; discontinua = mediana; banda oscura = rango intercuartílico; banda clara = "
                   "percentiles 10–90. Las bandas describen cuánto cambian los años entre sí; no son intervalos de "
                   "confianza.", script="20_climatologia_1_5a.py")
            + figura(CL / "fig_anio_mes_1_5a.png", "Mapa año × mes",
                     "Cada celda es un mes de un año, de modo que se ven los años individuales. El caudal va en escala "
                     "logarítmica. En temperatura, la trama marca los valores estimados con ERA5-Land.",
                     script="20_climatologia_1_5a.py")
            + figura(CL / "fig_ciclo_1_5b.png", "Rasgos del ciclo anual",
                     f"(a) Forma y fase: el caudal va {lag_centroide('P_L'):.1f} meses detrás de P<sub>L</sub> según el "
                     f"centroide. (b) El máximo anual de P<sub>L</sub> cae en marzo en el {pct_mar:.0f} % de los años "
                     f"(feb–abr en el {pct_fa:.0f} %), y el del caudal en abril en el {pct_qabr:.0f} %. (c) Los dos "
                     f"subperiodos mantienen el mismo régimen. (d) Variabilidad entre años por mes; en los meses secos "
                     f"el CV es inestable porque la media es casi cero.",
                     "SI = índice de estacionalidad (Walsh & Lawler, 1981); PCI = índice de concentración (Oliver, 1980); "
                     "centroide = Markham (1970).", script="21_ciclo_anual_1_5b.py")
            + figura(CL / "fig_anomalias_1_5b.png", "Anomalías respecto al ciclo anual",
                     f"Diferencia de cada mes respecto a la media de su mes calendario, con la media móvil de 12 meses. "
                     f"Años más húmedos según P<sub>L</sub>: {humedos}; más secos: {secos}.",
                     script="21_ciclo_anual_1_5b.py")
            + tabla(t15, "Indicadores de estacionalidad (periodo común)",
                    "R = caudal en lámina. Mes seco según Bagnouls–Gaussen (P < 2T). Criterios de régimen: unimodal si "
                    "hay un solo pico y el armónico anual domina; bimodal si hay dos picos y el semianual es comparable; "
                    "estacionalidad débil si SI < 0.40."),
            "Punto 1.5"),
]

nav = [("cuenca", "Cuenca"), ("contexto", "Contexto"), ("p11", "1.1 Series"), ("p12", "1.2 IMERG"),
       ("temp", "Temperatura"), ("p13", "1.3 Distribución"), ("p14", "1.4 Calidad"), ("p15", "1.5 Ciclo anual"),
       ("datos", "Datos")]

CSS = """
:root{--bg:#f6f5f2;--card:#ffffff;--ink:#1d1d1b;--ink2:#55534e;--line:#e2e0da;--acc:#1b6ca8;--chip:#e7f0f8;--code:#f0efeb}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){--bg:#151514;--card:#1f1f1d;--ink:#ecebe7;--ink2:#a9a7a1;--line:#34332f;--acc:#6aaee6;--chip:#1c2d3d;--code:#2a2a27}}
:root[data-theme="dark"]{--bg:#151514;--card:#1f1f1d;--ink:#ecebe7;--ink2:#a9a7a1;--line:#34332f;--acc:#6aaee6;--chip:#1c2d3d;--code:#2a2a27}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);font:16px/1.6 system-ui,-apple-system,"Segoe UI",Roboto,Ubuntu,sans-serif}
header.top{background:var(--card);border-bottom:1px solid var(--line)}
.wrap{max-width:1180px;margin:0 auto;padding:0 16px}
header.top .wrap{padding-top:28px;padding-bottom:22px}
h1{font-size:1.75rem;line-height:1.25;margin:0 0 6px}
.sub{color:var(--ink2);margin:0}
.datos-clave{display:flex;flex-wrap:wrap;gap:10px;margin-top:16px}
.dato{background:var(--bg);border:1px solid var(--line);border-radius:8px;padding:8px 12px;font-size:.9rem}
.dato b{display:block;font-size:1.05rem}
nav{position:sticky;top:0;z-index:5;background:var(--card);border-bottom:1px solid var(--line)}
nav .wrap{display:flex;gap:4px;overflow-x:auto;padding-top:8px;padding-bottom:8px}
nav a{white-space:nowrap;color:var(--ink2);text-decoration:none;padding:6px 10px;border-radius:6px;font-size:.92rem}
nav a:hover,nav a:focus-visible{background:var(--chip);color:var(--ink)}
main .wrap{padding-top:8px;padding-bottom:40px}
section{padding-top:36px;scroll-margin-top:56px}
h2{font-size:1.4rem;margin:0 0 8px;display:flex;align-items:center;gap:10px;flex-wrap:wrap}
.chip{font-size:.75rem;font-weight:600;background:var(--chip);color:var(--acc);padding:3px 9px;border-radius:999px}
.intro{color:var(--ink2);max-width:880px;margin:0 0 18px}
.fig{background:var(--card);border:1px solid var(--line);border-radius:12px;margin:0 0 22px;overflow:hidden}
.zoom{all:unset;display:block;cursor:zoom-in;background:#fff;padding:10px}
.zoom:focus-visible{outline:3px solid var(--acc);outline-offset:-3px}
.fig img{display:block;max-width:100%;max-height:82vh;width:auto;height:auto;margin:0 auto}
figcaption{padding:14px 18px 16px;border-top:1px solid var(--line)}
figcaption h4{margin:0 0 6px;font-size:1.05rem}
figcaption p{margin:0 0 8px;max-width:900px}
.leer{color:var(--ink2);font-size:.95rem}
.meta{font-size:.82rem;color:var(--ink2);display:flex;flex-wrap:wrap;gap:10px;align-items:center;margin:0}
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
@media (max-width:600px){h1{font-size:1.35rem}figcaption{padding:12px}.zoom{padding:4px}}
"""

JS = """
const visor=document.getElementById('visor'),vimg=visor.querySelector('img');
document.querySelectorAll('.zoom').forEach(b=>b.addEventListener('click',()=>{const i=b.querySelector('img');vimg.src=i.src;vimg.alt=i.alt;visor.classList.add('abierto');visor.querySelector('button').focus();}));
function cerrar(){visor.classList.remove('abierto');vimg.removeAttribute('src');}
visor.addEventListener('click',cerrar);
document.addEventListener('keydown',e=>{if(e.key==='Escape')cerrar();});
const raiz=document.documentElement,bt=document.querySelector('.tema');
function aplicar(t){if(t)raiz.setAttribute('data-theme',t);else raiz.removeAttribute('data-theme');}
try{aplicar(localStorage.getItem('tema'));}catch(e){}
bt.addEventListener('click',()=>{const oscuro=raiz.getAttribute('data-theme')==='dark'||(!raiz.getAttribute('data-theme')&&matchMedia('(prefers-color-scheme: dark)').matches);const t=oscuro?'light':'dark';aplicar(t);try{localStorage.setItem('tema',t);}catch(e){}});
"""

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
  <p class="sub">Tarea 1 · Hidrología 2026-2S · Universidad Nacional de Colombia, sede Medellín. Visor de los productos
  construidos hasta ahora, con una breve explicación de cada uno.</p>
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
<code>06 → 07a → 13 → 14 → 09 → 10 → 07b → 11 → 12 → 17 → 15 → 16 → 18 → 19 → 20 → 21 → 22</code>.
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
<script>{JS}</script>
</body>
</html>"""

SALIDA.write_text(pagina, encoding="utf-8")
print(f"Visor: {SALIDA} ({SALIDA.stat().st_size / 2**20:.1f} MB, {pagina.count('<figure')} figuras, "
      f"{pagina.count('<table')} tablas)")
