# Estado del proyecto – Tarea 1 Hidrología 2026-2S

**Cuenca:** río Chicama hasta la estación El Tambo (CAMELS-PE `PE_47E0A2A8`, SENAMHI 211003) ·
−7.574°, −78.707°, 684 m s.n.m. · vertiente del Pacífico, norte del Perú.
**Actualizado:** 2026-10-07. Rutas relativas a la raíz del proyecto (`Hidrologia/`).
Entorno: `Hidrologiat1/` (venv, Python 3.14); dependencias en `requirements.txt`.

---

## 1. Datos disponibles

| Conjunto | Ruta | Contenido / periodo | Estado |
|---|---|---|---|
| CAMELS-PE v1.0.1 (ZIP completo) | `Tarea 1/CAMELS_PE/CAMELS-PE_v1.0.1.zip` | DOI 10.5281/zenodo.21195425; MD5 `13f127d381338eee0e35359c08dba199` verificado contra Zenodo | ✔ |
| CAMELS-PE, subconjunto El Tambo | `Tarea 1/El_Tambo/datos/camels_pe/` | Serie diaria `PE_47E0A2A8.csv` 1981–2025 (`prec` PISCOp v2.1, `flow_obs` SENAMHI en mm/d, `tmean` PISCOt v1.2), polígono, salida, atributos, `PROCEDENCIA.json` | ✔ |
| IMERG V07 Final mensual | `Tarea 1/El_Tambo/datos/imerg/manual/*.nc4` | GPM_3IMERGM.07, descarga manual (GES DISC, OPeNDAP). 332 archivos, 1998-01 a 2025-08, sin faltantes ni duplicados. Grupo `Grid`; `precipitation` en mm/h y `gaugeRelativeWeighting`; malla 0.1°, lon −78.95…−78.15, lat −8.15…−7.25 | ✔ descargado · ⏳ sin procesar |
| DEM Copernicus GLO-30 (DSM, 30 m) | `Tarea 1/El_Tambo/datos/dem/` | Tiles originales (`GLO30_tiles/`), mosaico EPSG:4326 y re-proyección UTM 17S a 30 m | ✔ |
| Cobertura ESA WorldCover 2021 (10 m) | `Tarea 1/El_Tambo/datos/cobertura/` | Recorte de la cuenca | ✔ |
| Geología INGEMMET 1:100 000 | `Tarea 1/El_Tambo/datos/geologia/` | Unidades, fallas y cartera minera (GEOCATMIN, ArcGIS REST) | ✔ |
| Suelos SoilGrids WRB (250 m) | `Tarea 1/El_Tambo/datos/suelos/` | Descargado; retirado de los mapas a pedido del grupo | no usado |
| Estaciones SENAMHI | `Tarea 1/El_Tambo/datos/estaciones/senamhi_estaciones_20261006.csv` | Red actual (969 estaciones; mapa web de SENAMHI) | ✔ |
| Límites de países | `Tarea 1/El_Tambo/datos/ne_50m_admin_0_countries.geojson` | Natural Earth 1:50m (mapas de ubicación) | ✔ |
| PISCOp oficial SENAMHI (figshare, mensual) | `Tarea 1/El_Tambo/datos/piscop_senamhi/{v2p1_update,v3p0}/` | v2.1 update 1981–2019 (doi:10.6084/m9.figshare.21127423.v2) y v3.0 1981–2025 (doi:10.6084/m9.figshare.32411886.v1); MD5 verificado; malla 0.1° | ✔ descargado · **no reproduce CAMELS-PE** |
| **PISCOp v3.0 diario** | `…/piscop_senamhi/v3p0/PISCOp_d.nc` (1.5 GB, nacional, **fuera del ZIP**) y `PISCOp_v3_d_El_Tambo_recorte.nc` (2.9 MB, **va en el ZIP**) | 1981-01-01 a 2025-12-31, mm/día, 0.1°; 31 celdas en la cuenca. 0 faltantes, 0 negativos, 0 secuencias constantes; Σ diario = archivo mensual (máx 3.4 mm/mes). Cumple la guía: diario, 28 años comunes con Q (1998–2025), 0 % faltantes | ✔ 2026-10-08 |
| Temperatura ERA5-Land | `Tarea 1/El_Tambo/datos/era5land/` | T2m mensual 1981-01 a 2025-12 (540 meses, sin faltantes), 0.1°, recorte de la cuenca + geopotencial; `PROCEDENCIA.json` | ✔ 2026-10-08 |
| Campos globales SST / ERA5 (punto 5) | — | — | ✗ pendiente |

**DEM:** el proyecto usa **Copernicus GLO-30**. FABDEM v1.2 solo aparece en los atributos topográficos
de CAMELS-PE; no se descargó ni se usa.
**Descartado:** PISCOp del IRI Data Library (decisión del grupo, 2026-10-07). Datos borrados; queda solo
`08_descarga_piscop.py`, fuera del flujo.

### Periodo común de análisis: 1998-01 a 2025-08 (332 meses)

**Decisión del grupo (2026-10-08):** el periodo común es el de los meses con P_L (CAMELS-PE) **y** P_I (IMERG)
a la vez. **Todas las figuras, mapas y diagramas** se construyen sobre él. Se calcula en un solo lugar,
`Codigos/periodo_comun.py` (a partir de los datos, no fijado a mano), y todos los scripts lo importan.
Las tablas del registro completo se conservan solo donde la guía pide distinguirlas (1.3 tabla A; 1.4
resumen de disponibilidad) y para las tendencias de largo plazo (punto 3, si el grupo lo decide). Media anual
en el periodo común = Σ de las 12 medias mensuales (2025 está incompleto).

| Serie | Meses válidos en el periodo común | Observación |
|---|---|---|
| P_L (CAMELS-PE `prec`) | 332/332 (540/540 en 1981–2025) | **0 días faltantes** en todo el registro |
| P_I (IMERG) | 332/332 | Registro 1998-01 a 2025-08; define el fin del periodo común |
| Q / R (CAMELS-PE `flow_obs`) | 305/332 | 3.11 % de días faltantes; 252 meses completos, 53 escalados, 26 excluidos y 1 sin dato |
| T (CAMELS-PE `tmean`) | 276/332 | Sin dato desde 2021-01 (estimación ERA5-Land aparte, 56 meses) |

---

## 2. Scripts (`Hidrologiat1/Codigos/`)

| Script | Función | Estado |
|---|---|---|
| `04_seleccion_cuenca_PE.py` | Tamizado CAMELS-PE: área 100–10 000 km², ≥ 25 años, ≤ 10 % faltantes, celdas IMERG por cuenca | ✔ ejecutado |
| `05_evaluacion_el_tambo.py` | Huecos, ventanas candidatas, balance anual P–R y escalas espaciales de El Tambo | ✔ ejecutado |
| `06_delimitacion_el_tambo_glo30.py` | Delimitación D8 (pysheds) con GLO-30, red de Strahler, río principal, morfometría, estaciones y mapa base | ✔ validado: 2170.0 km² vs 2171.6 (CAMELS-PE), IoU 0.994 |
| `07a_descarga_contexto_fisico.py` | Descarga WorldCover, INGEMMET y SoilGrids | ✔ |
| `07b_mapas_contexto_fisico.py` | Panel de contexto físico (pendiente, orientación, cobertura, grupos hidrogeológicos y (e) precipitación media anual del periodo común 1998-01 a 2025-08 de **PISCOp v3.0** en celdas de 0.1° sin interpolar: 371–1090 mm/año, cuenca 721; CAMELS-PE 841) y ficha | ✔ (e) actualizado 2026-10-08; requiere 14 y 10 |
| `12_mapa_pesos_imerg.py` | Mapa de celdas IMERG con su peso por área y desnivel del terreno dentro de cada celda (`mapa_celdas_imerg_pesos.png`, `relieve_celdas_imerg.csv`) | ✔ 2026-10-08 |
| `utils_mapas.py` | Cuadrícula en grados, barra de escala y flecha de norte | ✔ |
| `09_datos_camels_pe_el_tambo.py` | Verificación MD5 del DOI, extracción del subconjunto, series mensuales y disponibilidad año × mes | ✔ validado (revisión manual de 2010-03 y 1998-03) |
| `11_grafica_precipitacion_camels.py` | Serie de P_L 1981–2025 (mensual con media móvil de 12 meses, y anual) | ✔ |
| `10_imerg_mensual_y_grafica.py` | Lectura de IMERG (manual o automática), promedio ponderado, mm/h → mm/mes, gráfica P_L vs P_I | ✔ ejecutado 2026-10-07 (332 meses; corregida la verificación de malla float32). P_I 45.3 vs P_L 70.3 mm/mes, sesgo −24.9, MAE 29.1, RMSE 48.5, r 0.908 |
| `13_piscop_senamhi_verificacion.py` | Descarga PISCOp oficial (v2.1 update y v3.0, figshare) y compara su promedio de cuenca con `prec` de CAMELS-PE (2 polígonos × 2 métodos de pesos × 3 periodos) → `resultados/precipitacion/verificacion_piscop_senamhi_vs_camels.csv` | ✔ 2026-10-08 |
| `22_visor_html.py` | Genera `visor_resultados.html` (raíz del proyecto): página autocontenida (~6 MB, figuras incrustadas en JPEG, sin dependencias externas) con las 18 figuras, 4 tablas y descripciones; cifras leídas de los CSV de resultados; tema claro/oscuro y ampliación al clic | ✔ 2026-10-08 |
| `21_ciclo_anual_1_5b.py` | Punto 1.5.b: indicadores de estacionalidad (SI, PCI, centroide de Markham, picos, armónicos, meses húmedos y secos de Gaussen), desfases, variabilidad por mes (CV, G1, jackknife), mes del pico por año hidrológico, subperiodos y anomalías por año hidrológico y calendario (`resultados/climatologia/`) | ✔ 2026-10-08 |
| `20_climatologia_1_5a.py` | Punto 1.5.a: climatología de 12 meses en el periodo común (n años, media, mediana, desv. est., Q1, Q3, P10, P90 por mes) de P_L, P_I, Q̄, R y T (PISCOt + estimación desde 2021); ciclo con bandas IQR y P10–P90; mapas año × mes (`resultados/climatologia/`) | ✔ 2026-10-08 |
| `19_control_calidad_1_4.py` | Punto 1.4: mapa y tablas de disponibilidad año × mes, comprobaciones diarias (fechas, códigos, negativos, tmean, flow_sim), secuencias constantes, tramos lineales, picos aislados, casos 1998 y diciembres, revisión manual de IMERG y registro de 17 anomalías (`resultados/control_calidad/`) | ✔ 2026-10-08 |
| `18_estadisticos_1_3.py` | Punto 1.3: tablas de estadísticos (registro completo, meses comunes de las 5 variables y pares P_L–P_I), meses extremos con valores simultáneos, histogramas (Freedman–Diaconis, frecuencia relativa) y diagramas de caja (Tukey 1.5·IQR) | ✔ 2026-10-08 |
| `17_era5land_temperatura.py` | Descarga ERA5-Land mensual T2m 1981–2025 por CDS (doi:10.24381/cds.68d2bb30; credenciales en `~/.cdsapirc`, fuera del ZIP), promedio ponderado por área (30 celdas), comparación con `tmean` CAMELS-PE 1981–2020 y relieve del modelo (`fig_T_camels_vs_era5land.png`, `era5land_t2m_cuenca_mensual.csv`) | ✔ 2026-10-08 |
| `16_grafica_precipitacion_imerg.py` | Punto 1.2: serie de P_I 1998-01 a 2025-08 (mensual con media móvil de 12 meses y totales anuales; 2025 incompleto excluido), mismo formato y ejes que la figura de CAMELS-PE, color #e4572e (`fig_precipitacion_imerg_1998_2025.png`). Media 539 mm/año (1998–2024) | ✔ 2026-10-08 |
| `15_figura_series_1_1.py` | Punto 1.1: P_L, Q̄ (m³/s), R (mm/mes) y T en paneles alineados en el periodo común (1998-01 a 2025-08), vacíos sombreados y meses escalados marcados (`fig_series_mensuales_1_1.png`, `extremos_series_mensuales_1_1.csv`). R > P_L en 98/305 meses (sobre todo may–ago; 1998-03: R 837 vs P 414 mm) | ✔ 2026-10-08 |
| `14_piscop_v3_cuenca.py` | Recorte de PISCOp v3.0 diario, pesos por área (GLO-30), controles de calidad, series diaria y mensual de cuenca (`piscop_v3_cuenca_{diaria,mensual}.csv`) | ✔ 2026-10-08 |
| `08_descarga_piscop.py` | Descarga de PISCOp desde IRI | ✗ descartado (no usar) |
| `01_delimitacion.py`, `03.py` (vacíos), `02_Prueba.py` | Archivos previos (CAMELS-DE) | fuera del flujo |

**Orden de ejecución:** 06 → 07a → 13 → 14 → 09 → 10 → 07b → 11 → 12 → 17 → 15 → 16 → 18 → 19 → 20 → 21 → 22 (visor HTML, siempre al final). El periodo común
(`periodo_comun.py`) necesita 09 y 10, así que 07b, 11 y 15–19 van después; 15 requiere además la extensión de T del 17.

**Productos principales:**
- `Tarea 1/El_Tambo/resultados/delimitacion/`: polígono, red, río principal, `parametros_morfometricos.csv`, `mapa_cuenca_el_tambo.png`.
- `Tarea 1/El_Tambo/resultados/contexto_fisico/`: `panel_contexto_fisico.png`, mapas individuales, `ficha_contexto_fisico.csv`.
- `Tarea 1/El_Tambo/resultados/series_mensuales/`: `series_mensuales_camels_pe*.csv`, `disponibilidad_anio_mes.csv`, `fig_precipitacion_camels_1998_2025.png`.

---

## 3. Decisiones metodológicas

**Fuentes**
- P_L, Q y T se toman solo de CAMELS-PE v1.0.1 (DOI anterior). P_L no es una observación independiente: PISCOp v2.1 combina pluviómetros con CHIRPS, y debe declararse así al compararla con IMERG Final, que también incorpora pluviómetros (GPCC).
- PISCOp público del IRI, promediado sobre la cuenca, daba 20–25 % menos que `prec` de CAMELS-PE (55.5 vs 69.6 mm/mes en 1998–2016; r = 0.98). La causa no está documentada; se descartó y se reporta como incertidumbre de la lluvia de referencia.
- **Verificación con PISCOp oficial de SENAMHI (2026-10-08, script 13):** tampoco reproduce CAMELS-PE. Cuenca El Tambo 1998–2025: v3.0 = 718 vs 841 mm/año (−14.5 %); v2.1 update 1998–2019 = 704 vs 838 (−16 %); r mensual 0.98–0.99. No depende del polígono (GLO-30 vs CAMELS-PE) ni de los pesos (área vs centro de celda). El cociente CAMELS/PISCOp no es constante: 1.06 (dic) a 1.46 (ago); en Puente Coina la diferencia es mayor (−26 % con v3.0). Conclusión: la `prec` de CAMELS-PE no es reproducible con ninguna versión pública de PISCOp (IRI stable/unstable, v2.1 update, v3.0). No se verificó el archivo diario de v3.0 (1.5 GB).

**Delimitación y espacialización**
- Polígono de trabajo: delimitación propia con GLO-30 (UTM 17S, 30 m, D8, umbral de red 5 km²). El cierre se ajustó a la celda con área drenada a ±15 % de la de CAMELS-PE (desplazamiento de 0 m).
- **Promedio de cuenca de productos en malla (IMERG):** media ponderada por el **área de intersección celda–polígono**, calculada en UTM 17S. Las celdas de borde pesan según su fracción dentro de la cuenca y los pesos se renormalizan si una celda no tiene dato. No se usa el píxel de la estación. IMERG: 31 celdas tocan la cuenca, 10 completas (≥ 99 % de la celda dentro), 43.7 % del peso en celdas de borde (≈ 17.8 celdas equivalentes). **Los pesos no consideran la topografía**; el desnivel interno de las celdas (P95−P5) tiene mediana de 1233 m y máximo de 2233 m.
- Mapas en coordenadas geográficas WGS 84 (EPSG:4326); cálculos de área y morfometría en UTM 17S (EPSG:32717).

**Agregación mensual y vacíos**
- P_L,m = Σ diario; solo meses completos (todos lo están).
- R_m = Σ q_d (q ya está en mm/d; no se divide otra vez por el área). Mes completo → suma directa; ≤ 10 % de días faltantes → suma de los días válidos escalada a n_días (`R_flag = escalado`); > 10 % → excluido. No se presentan sumas parciales como acumulados completos.
- Q en m³/s = q[mm/d] · A / 86.4, con A = 2171.585 km² (área CAMELS-PE).
- T mensual = media de días válidos (≤ 10 % faltantes).
- P_I: tasa mm/h × 24 × días del mes.
- **Relleno de vacíos: ninguno.** P_L no tiene vacíos. Se decidió **no rellenar P_L con IMERG corregido**, porque contaminaría la comparación P_I–P_L y los modelos del punto 2. Tampoco se rellena Q: el periodo de análisis empieza en 1998 y los meses incompletos se excluyen. Para T 2021–2025 se usará ERA5-Land evaluado contra PISCOt en 1998–2020, sin unir ambas fuentes en una sola serie.

**Contexto físico**
- Litología agrupada en 6 grupos hidrogeológicos más lagunas (asignación en `unidades_geologicas_cuenca.csv`); pendiente de contrastar con INGEMMET (Serie H) y ANA.
- Hipótesis de humedad: fuente principal amazónica/atlántica por flujo del este en niveles medios y altos (dic–abr); Pacífico secundario, relevante en El Niño. Es una hipótesis a contrastar en el punto 5 (flujo de humedad de ERA5, SST Niño 1+2).

---

- **Decisión del grupo (2026-10-08): P_L principal = `prec` de CAMELS-PE; PISCOp v3.0 = fuente de apoyo** (mapa continuo (e) y contraste de incertidumbre). Razones: coherencia con la base elegida (Q, T, atributos), trazabilidad de la entrega y balance hídrico menos inverosímil con la P más alta (R/P 0.48 vs 0.56; ambas por debajo de lo esperado por Budyko). P v3.0 1998–2025 = 720 mm/año vs 841 CAMELS-PE (r mensual 0.983); cociente anual 0.8–0.9, menor en años muy lluviosos (1983: 0.62; 1998: 0.74). Scripts 09, 10 y 11 siguen con CAMELS-PE. No se agrega comparación de las tres lluvias (decisión del grupo). En el periodo común 1998-01 a 2025-08 (Σ medias mensuales): PISCOp v3.0 721 vs CAMELS-PE 841 mm/año.

## 4. Pendientes inmediatos
1. Interpretar el sesgo de IMERG (P_I/P_L = 0.645; subestima en dic–abr, sobreestima en jun–ago) y declarar que P_L no es una observación independiente.
2. Temperatura: ERA5-Land ya descargado (1981–2025) y comparado (script 17). **Resultado en el periodo común (meses con PISCOt: 1998-01 a 2020-12, 276 meses):** misma variabilidad (r series 0.89, r anomalías 0.90; amplitud del ciclo anual 0.57 vs 0.60 °C, ambos con máximo ago–sep) pero ERA5-Land **2.84 °C más frío**, desfase casi constante en todos los meses. (Con 1981–2020 era −2.92 °C.) Relieve del modelo 2641 m vs 2557 m reales (+84 m) → explica solo ≈ −0.55 °C. Otras causas a evaluar: `tmean` de PISCOt = (Tmax+Tmin)/2 vs media de 24 h de ERA5-Land; sesgo del reanálisis en los Andes; interpolación de PISCOt. Verificar con normales de estaciones SENAMHI cercanas cuál nivel es más realista. Para anomalías, tendencias y correlaciones el desfase constante no pesa; para valores absolutos (PET, aridez) sí. **Extensión de T desde 2021 (2026-10-08, pedido del grupo):** T_est = ERA5-Land + sesgo medio del mes calendario frente a PISCOt en el periodo común 1998–2020 (+2.62 a +2.99 °C); 56 meses estimados, 2021-01 a 2025-08. Validación fuera de muestra (sesgo ajustado en 1998–2010, evaluado en 2011–2020): MAE 0.19 °C, sesgo −0.05, r 0.89; climatología como referencia: MAE 0.40 °C (`validacion_extension_T.csv`). Archivo `temperatura_piscot_extendida_era5land.csv` (columna aparte y bandera de fuente; `T_C` original intacta). Se muestra en la fig. 1.1 como tramo discontinuo morado. En las tendencias del punto 3, analizar PISCOt 1981–2020 por separado y tratar 2021–2025 como estimado (riesgo de salto por cambio de fuente). Enfoque acordado: **no unir fuentes**. PISCOt v1.2 (CAMELS-PE, 1981–2020; no existe versión posterior) para tendencias de largo plazo; ERA5-Land T2m mensual 1998–2025 como serie completa para el periodo común (1.5, 4, 5), evaluada contra PISCOt en 1998–2020 (sesgo por altitud: ERA5-Land ~9 km vs relieve 684–>4000 m). Requiere cuenta CDS y `~/.cdsapirc` (fuera del ZIP).
3. Puntos 1.1 a 1.4, 1.5.a y 1.5.b hechos (scripts 15, 10, 16, 18, 19, 20 y 21). Sigue: 1.5.c (procesos regionales y de la cuenca; síntesis y clasificación hidroclimática) e hipótesis inicial, que son sobre todo interpretación del grupo con bibliografía.
4. Resolver la diferencia PISCOp (IRI, v2.1 update y v3.0 oficiales, todos 15–25 % menores) vs CAMELS-PE con los autores (hllauca@senamhi.gob.pe) o declararla como limitación.
5. Decidir si se borra `08_descarga_piscop.py` antes de armar el ZIP. El ZIP no debe incluir `PISCOp_d.nc` (1.5 GB); basta el recorte de la cuenca.

---

## 5. Notas para el informe

### Punto 1.1 – Series mensuales (fig. `fig_series_mensuales_1_1.png`)
- Criterio de completitud a declarar: la guía da fórmulas "para meses completos"; además de los 256 meses completos de Q, se conservan 53 meses con ≤ 10 % de días faltantes (Q̄ = media de días válidos; R escalada a n_m) y se excluyen 27 (26 excluidos + 1 sin dato). R = Σq_d sin repetir 86.4/A porque `flow_obs` ya está en mm/d.
- Extremos: 03/1998 es el máximo de P_L (414 mm), Q̄ (678 m³/s) y R (837 mm), El Niño 1997–98; segundo pico de Q en 03/2017 (El Niño costero). Mínimos: Q̄ 0.39 m³/s (09/2004), P_L 0.2 mm (07/2012).
- R > P_L en 98/305 meses del periodo común, casi todos en may–ago (estación seca): posible liberación de almacenamiento (la guía advierte que no es error por sí solo). 03/1998: R ≈ 2 × P_L → revisar en 1.4 (almacenamiento, subestimación de P_L en el evento o curva de calibración en crecidas).

### Punto 1.3 – Distribución (figs. `fig_histogramas_1_3.png`, `fig_cajas_1_3.png`; tablas `estadisticos_1_3_*.csv`, `extremos_1_3.csv`)
- Convenciones a declarar: faltantes excluidos sin relleno; percentiles tipo 7 (interpolación lineal, numpy); desv. est. muestral; asimetría G1; histogramas con ancho Freedman–Diaconis redondeado y frecuencia relativa (clases de igual ancho); P_L y P_I con los mismos 332 meses y límites (clases de 20 mm); cajas de Tukey 1.5·IQR. T solo PISCOt (la estimación 2021–2025 no entra).
- Dos muestras: (A) registro completo, solo tabla (P_L 540, P_I 332, Q̄ y R 309, T 480 meses); (B) periodo común 1998-01 a 2025-08 (P_L y P_I 332, Q̄ y R 305, T 276), base de todas las figuras y de los extremos. En (B): P_L media 70.3 / mediana 47.3 mm; P_I 45.3 / 33.7; Q̄ 28.1 / 12.1 m³/s; T 16.22 / 16.18 °C; P_L < 1 mm en 15 meses. Las cifras de asimetría de abajo son del registro completo; en (B) son casi iguales (P_L 1.74, P_I 1.70, Q̄ 7.47, R 7.61, T 0.40).
- Asimetría: P_L 1.83 y P_I 1.70 (asimetría positiva moderada; media ≫ mediana: 66.4 vs 44.3 y 45.3 vs 33.7 mm); Q̄ 7.5 y R 7.7 (muy fuerte: media 27.9 vs mediana 11.1 m³/s, dominada por 1998-03 y 2017-04); T 0.18 (casi simétrica, media = mediana = 16.1 °C, IQR 0.72 °C). Para condición típica usar mediana/IQR en P y Q; media en T.
- Ceros: ningún cero exacto en ninguna variable. P_L < 1 mm en 26 meses (estación seca); P_I nunca < 1 mm (mínimo 1.97) → IMERG no reproduce los meses casi secos; Q̄ < 1 m³/s en 9 meses.
- P_I vs P_L (pares comunes): P_I tiene más meses en 0–80 mm y casi ninguno sobre 200 mm (P95 141 vs 229 mm; máx 265 vs 414) → comprime la distribución en ambos extremos (sobreestima meses secos, subestima húmedos).
- Extremos coherentes: 1998-03 es máximo de P_L, P_I, Q̄ y R (El Niño 1997–98); 2017-03/04 (El Niño costero) segundo máximo de Q̄ y R. 1998-02 (2.º máximo de P_L, 381 mm) no tiene Q̄: mes excluido por faltantes, posible pérdida de registro en la crecida → revisar en 1.4. Mínimos de Q̄ en 2004-09, 2016-12 y 2005-12: dos en diciembre con P_L 68–79 mm → revisar en 1.4 (fin de la recesión con inicio tardío de lluvias, extracciones o problema de datos). T máxima 2015-12/2016-01 (El Niño 2015–16), mínimas 1984 y 1989.

### Punto 1.5.a – Climatología (`resultados/climatologia/`: `fig_climatologia_1_5a.png`, `fig_anio_mes_1_5a.png`, `climatologia_1_5a.csv`)
- Periodo común 1998-01 a 2025-08: ene–ago 28 años, sep–dic 27. P_L y P_I con los mismos meses. Q̄/R con sus meses válidos: 20–28 años por mes (sep solo 20, por los faltantes sistemáticos de estiaje; ver CC-04). T = PISCOt hasta 2020-12 + estimación ERA5-Land desde 2021-01 (4–5 años estimados por mes, marcados con trama en el mapa año × mes). Bandas = dispersión entre años, no intervalos de confianza. Percentiles tipo 7.
- P_L: un solo máximo en marzo (media 213, mediana 204 mm), mínimo jul–ago (~3 mm); 81 % del total anual medio en dic–abr; meseta secundaria oct–nov (~48 mm). P_I: mismo máximo en marzo pero mucho menor (135 mm); la mayor subestimación es febrero (84 vs 159 mm); sobreestima jun–ago (8–11 vs 3–9 mm); 73 % en dic–abr.
- Q̄: máximo en abril (media 90, mediana 70 m³/s), un mes después del máximo de P_L; en marzo media 86 vs mediana 51 (asimetría por 1998 y 2017). Mínimo sep–oct (~3 m³/s), uno o dos meses después del mínimo de lluvia (jul–ago) → recesión. 76 % del caudal anual medio en dic–abr.
- T: ciclo débil (amplitud 0.62 °C); máximo en septiembre (16.7 °C), al final de la estación seca y antes de las lluvias; mínimo en junio (16.1 °C). La dispersión entre años (P10–P90 ≈ 1.5 °C en ene–mar) supera la amplitud del ciclo.
- Mapa año × mes: 1998 y 2017 destacan en P_L y Q̄ (marzo–abril); 2015–16 y 2023–24 (este último estimado) son los más cálidos.

### Punto 1.5.b – Rasgos del ciclo (`fig_ciclo_1_5b.png`, `fig_anomalias_1_5b.png`, `indicadores_estacionalidad_1_5b.csv`, `variabilidad_por_mes_1_5b.csv`, `mes_pico_por_anio_1_5b.csv`, `subperiodos_1_5b.csv`, `anomalias_anio_{hidrologico,calendario}_1_5b.csv`)
- Criterios (propuestos; el grupo puede ajustarlos): SI de Walsh & Lawler (1981), PCI de Oliver (1980), centroide y concentración de Markham (1970); pico = máximo local ≥ media mensual con prominencia ≥ 10 % de la amplitud; UNIMODAL si 1 pico y H1 ≥ 2·H2, BIMODAL si 2 picos y H2 ≥ 0.5·H1, DÉBIL si SI < 0.40; mes húmedo = media ≥ total/12; mes seco = Bagnouls–Gaussen P < 2T. Año hidrológico sep–ago (27 completos; 1997-98 queda fuera del periodo común) y año calendario (27 completos, incluye 1998) para años contrastantes.
- Clasificación: **unimodal** con las tres series. P_L: SI 0.79 (estacional, en el límite con "marcadamente estacional"), PCI 15 (estacionalidad moderada, cerca de "estacional"), pico en marzo, centroide 2.1 (inicio de feb), concentración 0.57, H1 78 % vs H2 15 %; 5 meses húmedos (dic–abr) y 5 secos de Gaussen (may–sep). P_I: SI 0.63, PCI 14 (más uniforme). R: SI 0.89, PCI 19 (más concentrado que la lluvia), pico en abril, 4 meses húmedos (feb–may). T: ciclo débil (0.62 °C) y no unimodal (H1 33 %, H2 38 %): máximo sep, mínimo jun.
- Desfases Q–P: centroide 1.4 meses con P_L (1.5 con P_I); máximo 1 mes (mar → abr); mínimo 2 meses con P_L (jul → sep), 1 con P_I (ago → sep). r entre ciclos máxima con 1 mes de desplazamiento (0.94 P_L; 0.91 P_I). La interpretación del desfase casi no cambia con IMERG.
- ¿El pico medio aparece en la mayoría de los años? P_L: máximo en marzo en 59 % de los años hidrológicos y en feb–abr en 93 %; P_I: marzo 70 %; Q̄ (15 años completos): abril 73 %, feb–abr 100 %.
- Variabilidad: CV de P en meses húmedos 0.3–0.55; jun–sep inestables (media cercana a cero). Q̄: CV máximo en ene (1.26) y mar (1.43) por 1998; jackknife: quitar 1998 cambia la media de marzo de Q̄ en 25 % y la de enero en 22 %; quitar 2017 cambia abril en 15 %. En P_L ningún año mueve la media de un mes húmedo más de 7 %.
- Subperiodos 1998–2011 vs 2012–2025-08: clasificación estable (unimodal, SI 0.80 vs 0.77). El pico de R pasa de marzo a abril, pero el de 1998–2011 lo arrastra 1998 (no es cambio de régimen). P_L: feb menor (173 → 145 mm), oct y dic mayores. T 0.4 °C más alta en el segundo (incluye la estimación 2021–2025). No se atribuye a tendencia: se evalúa en el punto 3.
- Años contrastantes (año calendario, anomalía de P_L): húmedos 1998 (+580 mm), 2023 (+366), 2009 (+182), 2017 (+162); secos 2003 (−293), 2016 (−280), 2005 (−274), 2004 (−267). 2017: R +494 mm, mucho mayor que P_L (+162). 2023: P_L +366 pero P_I −74 (IMERG no capta el evento) y R solo +37 → revisar. 2016 seco y cálido (+0.7 °C, El Niño 2015–16). r de anomalías anuales P_L–P_I 0.67; ene–abr P_L–R 0.75 (n = 25), P_I–R 0.60. En 9 de 27 años P_L y P_I tienen anomalía de signo opuesto.

### Punto 1.4 – Control de calidad (`resultados/control_calidad/`: `fig_disponibilidad_1_4.png`, `registro_anomalias_1_4.csv`)
- Sin problemas: fechas continuas y sin duplicados, NA como único código de faltante, sin negativos, tmin ≤ tmean ≤ tmax, sin tramos interpolados linealmente ni picos aislados en Q, una sola secuencia constante (5 días, 2007-09, estiaje). flow_sim no se usa.
- Hallazgo principal: los faltantes de Q son sistemáticos en estiaje (meses incompletos por mes calendario: sep 26, ago 11, oct 10, jul 8, nov 8; resto ≤ 5) → posible mantenimiento anual; causa no documentada. Cambio de cobertura en 2019–2021 (61, 56 y 24 días faltantes/año).
- 1998-02 excluido (faltan 5 días en plena crecida de El Niño); 1998-03 escalado con R ≈ 2 × P_L (incierto: almacenamiento, P_L subestimada o curva de descarga en crecida).
- Mínimos de diciembre (2005, 2016): datos completos y sin patrones sospechosos; extracciones registradas ANA ≈ 0.002 m³/s (despreciables); hipótesis: fin de recesión con primeras lluvias absorbidas o extracciones no registradas.
- tmean de PISCOt = (tmin+tmax)/2 exacto → parte del desfase con ERA5-Land.
- Revisiones manuales que coinciden: P y R de 2010-03 y 1998-03; IMERG 2010-03 (102.63 mm).
- Registro: 17 entradas (8 sin problema, 5 inciertas, 2 retenidas, 1 excluida, 1 corregida = T 2021–2025 reconstruida). Las decisiones de "incierto" deben revisarse y defenderse por el grupo.

### Punto 1.2 – IMERG (figs. `fig_P_L_vs_IMERG_mensual.png` y `fig_precipitacion_imerg_1998_2025.png`)
- **Procedencia de P_L (lo exige la guía):** P_L = `prec` de CAMELS-PE = PISCOp v2.1 (interpolación de pluviómetros SENAMHI con CHIRP/CHIRPS y climatologías) promediado sobre la cuenca. Es un producto interpolado, **no una observación local independiente**. No es reproducible con las versiones públicas de PISCOp (IRI, v2.1 update, v3.0 dan 15–25 % menos): declararlo como limitación. IMERG Final también incorpora pluviómetros (GPCC): fuentes compartidas posibles.
- **Documentación de IMERG que pide la guía** (ya en la cabecera del script 10 y en la nota de la figura): producto GPM IMERG Final Precipitation L3 1 month 0.1° V07 (GPM_3IMERGM.07, doi:10.5067/GPM/IMERG/3B-MONTH/07); modalidad Final; variable `precipitation` en mm/h (tasa media mensual) → mm/mes = tasa · 24 · n_días; malla 0.1°; 1998-01 a 2025-08 (332 meses, sin faltantes ni duplicados); acceso manual NASA GES DISC (Subset / OPeNDAP, grupo `Grid`), con descarga automática alternativa por OPeNDAP + Earthdata; promedio ponderado por área de intersección celda–cuenca (UTM 17S), 31 celdas, 10 completas, 43.7 % del peso en bordes, ≈ 17.8 celdas equivalentes (mapa `mapa_celdas_imerg_pesos.png`); no se usa el píxel de la estación; los pesos no consideran la topografía (desnivel interno de las celdas: mediana 1233 m).
- **Periodo común P_L–P_I:** 1998-01 a 2025-08, 332 pares mensuales. P_I 45.3 vs P_L 70.3 mm/mes; sesgo −24.9 mm/mes; MAE 29.1; RMSE 48.5; r = 0.908.
- **Cifras para describir la figura de IMERG:** total anual medio 539 mm/año (1998–2024; 2025 incompleto excluido) vs 833 mm/año de P_L en los mismos años (1998–2024; 797 mm/año en 1981–2025); máximo mensual 265 mm (03/1998) vs 414 mm; 03/2017: 202 vs 377 mm; año más seco 2020 (370 mm), más lluvioso 1998 (785 mm). IMERG reproduce los eventos de El Niño 1998 y 2017 con picos mucho menores; media móvil más plana que la de P_L; desde 2019 casi todos los años por debajo de su media → contrastar en el punto 3 (tendencia real vs cambio en los insumos del producto).
- **Balance hídrico como contraste independiente:** con R medido (~402–407 mm/año), R/P = 0.47–0.48 con P_L y 0.72 con P_I; P − R implícito 462 vs 157 mm/año. Ambas P quedan por debajo de lo esperado por Budyko (PET CAMELS-PE 1341 mm/año) → la brecha de IMERG parece venir sobre todo de subestimación de IMERG, sin descartar sesgo de P_L. Hipótesis a sustentar con bibliografía (lluvia orográfica cálida mal detectada por microondas en los Andes).
