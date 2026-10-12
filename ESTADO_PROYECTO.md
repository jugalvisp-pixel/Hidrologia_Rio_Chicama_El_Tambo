# Estado del proyecto – Tarea 1 Hidrología 2026-2S

**Cuenca:** río Chicama hasta la estación El Tambo (CAMELS-PE `PE_47E0A2A8`, SENAMHI 211003) ·
−7.574°, −78.707°, 684 m s.n.m. · vertiente del Pacífico, norte del Perú.
**Actualizado:** 2026-10-12 (puntos 1–4 hechos; falta el 5). Rutas relativas a la raíz del proyecto (`Hidrologia/`).
Entorno: `Hidrologiat1/` (venv, Python 3.14); dependencias en `requirements.txt`.
**Repositorio:** https://github.com/jugalvisp-pixel/Hidrologia_Rio_Chicama_El_Tambo (público) ·
**Visor en línea:** https://jugalvisp-pixel.github.io/Hidrologia_Rio_Chicama_El_Tambo/visor_resultados.html
(ver sección 6).

### Para retomar (2026-10-12)
- **Hecho:** puntos 1 a 4 completos (scripts hasta 33) y visor reorganizado por secciones con ficha "Qué muestra / Cómo se construyó" en cada figura.
- **Siguiente:** (a) más correcciones del visor que pida el grupo, sección por sección; (b) punto 5 **hecho (2026-10-12)**: el grupo debe revisar la interpretación física propuesta en el visor y leer las fuentes; (c) síntesis de los puntos 4 y 5 y cierre de la tarea (guía).
- **Tentativas:** las figuras de inversión térmica (29) y síntesis 1.5.c (23) siguen siendo tentativas; el grupo no las ha aprobado.
- **Entorno:** el venv `Hidrologiat1` tiene `statsmodels` 0.15 y `pymannkendall` 1.4.3 (en `requirements.txt`). Correr los scripts con `Hidrologiat1\Scripts\python.exe` desde `Hidrologiat1/Codigos/`.
- **Visor bloqueado:** si `22_visor_html.py` falla con `OSError [Errno 22]` al escribir, el HTML está abierto en otro programa; cerrarlo y repetir.
- **Commits:** el usuario hace *Commit* y *Push* con GitHub Desktop; recordárselo al terminar.
- **Política del curso:** la IA no interpreta por el grupo. Las referencias citadas de memoria van marcadas para verificar.

---

## 1. Datos disponibles

| Conjunto | Ruta | Contenido / periodo | Estado |
|---|---|---|---|
| CAMELS-PE v1.0.1 (ZIP completo) | `Tarea 1/CAMELS_PE/CAMELS-PE_v1.0.1.zip` | DOI 10.5281/zenodo.21195425; MD5 `13f127d381338eee0e35359c08dba199` verificado contra Zenodo | ✔ |
| CAMELS-PE, subconjunto El Tambo | `Tarea 1/El_Tambo/datos/camels_pe/` | Serie diaria `PE_47E0A2A8.csv` 1981–2025 (`prec` PISCOp v2.1, `flow_obs` SENAMHI en mm/d, `tmean` PISCOt v1.2), polígono, salida, atributos, `PROCEDENCIA.json` | ✔ |
| IMERG V07 Final mensual | `Tarea 1/El_Tambo/datos/imerg/manual/*.nc4` | GPM_3IMERGM.07, descarga manual (GES DISC, OPeNDAP). 332 archivos, 1998-01 a 2025-08, sin faltantes ni duplicados. Grupo `Grid`; `precipitation` en mm/h y `gaugeRelativeWeighting`; malla 0.1°, lon −78.95…−78.15, lat −8.15…−7.25 | ✔ procesado (script 10) |
| DEM Copernicus GLO-30 (DSM, 30 m) | `Tarea 1/El_Tambo/datos/dem/` | Tiles originales (`GLO30_tiles/`), mosaico EPSG:4326 y re-proyección UTM 17S a 30 m | ✔ |
| Cobertura ESA WorldCover 2021 (10 m) | `Tarea 1/El_Tambo/datos/cobertura/` | Recorte de la cuenca | ✔ |
| Geología INGEMMET 1:100 000 | `Tarea 1/El_Tambo/datos/geologia/` | Unidades, fallas y cartera minera (GEOCATMIN, ArcGIS REST) | ✔ |
| Suelos SoilGrids WRB (250 m) | `Tarea 1/El_Tambo/datos/suelos/` | Descargado; retirado de los mapas a pedido del grupo | no usado |
| Estaciones SENAMHI | `Tarea 1/El_Tambo/datos/estaciones/senamhi_estaciones_20261006.csv` | Red actual (969 estaciones; mapa web de SENAMHI) | ✔ |
| Límites de países | `Tarea 1/El_Tambo/datos/ne_50m_admin_0_countries.geojson` | Natural Earth 1:50m (mapas de ubicación) | ✔ |
| PISCOp oficial SENAMHI (figshare, mensual) | `Tarea 1/El_Tambo/datos/piscop_senamhi/{v2p1_update,v3p0}/` | v2.1 update 1981–2019 (doi:10.6084/m9.figshare.21127423.v2) y v3.0 1981–2025 (doi:10.6084/m9.figshare.32411886.v1); MD5 verificado; malla 0.1° | ✔ descargado · **no reproduce CAMELS-PE** |
| **PISCOp v3.0 diario** | `…/piscop_senamhi/v3p0/PISCOp_d.nc` (1.5 GB, nacional, **fuera del ZIP**) y `PISCOp_v3_d_El_Tambo_recorte.nc` (2.9 MB, **va en el ZIP**) | 1981-01-01 a 2025-12-31, mm/día, 0.1°; 31 celdas en la cuenca. 0 faltantes, 0 negativos, 0 secuencias constantes; Σ diario = archivo mensual (máx 3.4 mm/mes). Cumple la guía: diario, 28 años comunes con Q (1998–2025), 0 % faltantes | ✔ 2026-10-08 |
| Temperatura ERA5-Land | `Tarea 1/El_Tambo/datos/era5land/` | T2m mensual 1981-01 a 2025-12 (540 meses, sin faltantes), 0.1°, recorte de la cuenca + geopotencial; `PROCEDENCIA.json` | ✔ 2026-10-08 |
| Índice Niño 1+2 (ERSST v5, CPC/NOAA) | `Tarea 1/El_Tambo/datos/indices_clima/ersst5.nino.mth.91-20.ascii` | Mensual 1950→2026-06, anomalías respecto a 1991–2020; usado en 3.6 | ✔ 2026-10-11 |
| Campos globales SST / ERA5 (punto 5) | — | — | ✗ pendiente |

**DEM:** el proyecto usa **Copernicus GLO-30**. FABDEM v1.2 solo aparece en los atributos topográficos
de CAMELS-PE; no se descargó ni se usa.
**Descartado:** PISCOp del IRI Data Library (decisión del grupo, 2026-10-07). Datos borrados; queda solo
`08_descarga_piscop.py`, fuera del flujo.

### Periodo común de análisis: 1998-01 a 2024-12 (324 meses, 27 años completos)

**Decisiones del grupo:**
- **2026-10-08:** el periodo común es el de los meses con P_L (CAMELS-PE) **y** P_I (IMERG) a la vez. **Todas las figuras, mapas y diagramas** se construyen sobre él.
- **2026-10-11:** solo **años calendario completos**. Antes era 1998-01 a 2025-08 (332 meses); ahora quedan fuera 2025-01..08.

Se calcula en un solo lugar, `Codigos/periodo_comun.py`, a partir de los datos (primer enero a último diciembre comunes, no fijado a mano), y todos los scripts lo importan. Toda la cadena se regeneró el 2026-10-11 con el periodo nuevo.

Las tablas del registro completo se conservan solo donde la guía pide distinguirlas (1.3 tabla A; 1.4 resumen de disponibilidad) y para las tendencias de largo plazo (punto 3, si el grupo lo decide). Con años completos, la media anual (Σ de las 12 medias mensuales) es la media de los totales anuales.

**Cifras con el periodo nuevo:** P_L 833 mm/año, P_I 539 mm/año (65 % de P_L), R 391 mm/año.

| Serie | Meses válidos en el periodo común | Observación |
|---|---|---|
| P_L (CAMELS-PE `prec`) | 324/324 (540/540 en 1981–2025) | **0 días faltantes** en todo el registro |
| P_I (IMERG) | 324/324 | Registro 1998-01 a 2025-08; 2025 queda fuera por incompleto |
| Q / R (CAMELS-PE `flow_obs`) | **324/324 rellenada** (297/324 observada) | Observado: 244 meses completos, 53 escalados, 26 excluidos y 1 sin dato. **Rellenado (2026-10-10):** 244 completos, 76 con días interpolados, 3 escalados y 1 reconstruido con modelo (1998-08). Ver sección 3, "Relleno del caudal" |
| T (CAMELS-PE `tmean`) | 276/324 | Sin dato desde 2021-01 (estimación ERA5-Land aparte, 48 meses en el periodo) |

> Las notas de la sección 5 se escribieron con el periodo anterior (1998-01..2025-08). Las cifras cambian poco, pero **hay que tomar los valores definitivos de los CSV regenerados** al redactar.

---

## 2. Scripts (`Hidrologiat1/Codigos/`)

| Script | Función | Estado |
|---|---|---|
| `04_seleccion_cuenca_PE.py` | Tamizado CAMELS-PE: área 100–10 000 km², ≥ 25 años, ≤ 10 % faltantes, celdas IMERG por cuenca | ✔ ejecutado |
| `05_evaluacion_el_tambo.py` | Huecos, ventanas candidatas, balance anual P–R y escalas espaciales de El Tambo | ✔ ejecutado |
| `06_delimitacion_el_tambo_glo30.py` | Delimitación D8 (pysheds) con GLO-30, red de Strahler, río principal, morfometría, estaciones y mapa base. **Mapa (2026-10-11):** disposición horizontal 18 × 10.5 in: mapa principal a la izquierda; arriba a la derecha "Ubicación en Sudamérica" y "Ubicación en Perú" (con Lima); abajo leyenda, barra de elevación horizontal y notas. Sin IVT ni vientos (decisión del grupo); ya no genera `comparacion_camels.png` (eliminada) | ✔ validado: 2170.0 km² vs 2171.6 (CAMELS-PE), IoU 0.994 |
| `07a_descarga_contexto_fisico.py` | Descarga WorldCover, INGEMMET y SoilGrids | ✔ |
| `07b_mapas_contexto_fisico.py` | Panel de contexto físico **1 × 4** (2026-10-12, antes 2 × 2; 31 × 14.5 in, leyendas debajo de cada mapa, barra de color de (d) en eje insertado para que los 4 mapas tengan el mismo tamaño): (a) orientación, (b) cobertura, (c) **unidades hidrogeológicas**: las 25 unidades de INGEMMET unificadas por comportamiento hídrico en 6 grupos más lagunas. Colores según la convención IAH/UNESCO (azul = poroso; verdes = fisurado o kárstico; marrones = acuitardo o acuífugo). La leyenda va de más a menos permeable e incluye los códigos de cada grupo; los límites se disuelven por grupo. Aparte, `mapa_geologia_unidades_ingemmet.png` muestra las unidades individuales con los colores oficiales de INGEMMET (simbología de GEOCATMIN en `datos/geologia/ingemmet_simbologia_geologia100k.json`), con límites sin los bordes entre hojas (±25 m) y códigos en polígonos ≥ 20 km² (2026-10-10). (d) precipitación media anual del periodo común 1998-01 a 2025-08 de **CHIRPS v2.0**, interpolada de forma bilineal a la malla del DEM, con isoyetas cada 100 mm (píxeles 319–945 mm/año; cuenca 636; CAMELS-PE 841). Ficha. **Pendiente retirada del panel** (decisión del grupo, 2026-10-10; se sigue calculando para la ficha y las zonas planas). Antes (d) era PISCOp v3.0 en celdas de 0.1° sin interpolar | ✔ 2026-10-10; requiere 26 |
| `27_mapas_precipitacion_mensual.py` | 12 mapas de la precipitación media de cada mes calendario sobre la cuenca (CHIRPS 0.05°, periodo común, interpolación bilineal solo para la representación), con escala discretizada de naranja (seco) a azul oscuro (húmedo) en el estilo del ejemplo del curso. Matriz 4 × 3 con una fila por trimestre (DEF, MAM, JJA, SON) y sus tres meses en orden; a la izquierda, el acumulado de P_L de cuenca del trimestre. P_L: DEF 351 mm (42 %), MAM 358 (43 %), JJA 15 (2 %), SON 116 (14 %). Recuadro con el promedio de cuenca de CHIRPS y de P_L (`resultados/climatologia/mapas_precipitacion_mensual_chirps.png`). Rango por mes: marzo 95–216 mm, julio 0–10 mm | ✔ 2026-10-10 |
| `26_chirps_cuenca.py` | CHIRPS v2.0 mensual 0.05° (doi:10.15780/G2RP4Q), recorte del IRI Data Library para el periodo común (332 meses, sin faltantes). **La X del IRI marca el borde oeste de la celda**: se corrige +0.025°, verificado contra el GeoTIFF oficial de UCSB de 2017-03 (diferencia 0.000; sin corregir, ~94 mm). Pesos por área (98 celdas, 71.1 equivalentes) → `chirps_cuenca_mensual.csv`. P media anual 636 mm/año (PISCOp v3 721; CAMELS-PE 841); r mensual con CAMELS-PE 0.915 | ✔ 2026-10-10 |
| `12_mapa_pesos_imerg.py` | Mapa de celdas IMERG con su peso por área y desnivel del terreno dentro de cada celda (`mapa_celdas_imerg_pesos.png`, `relieve_celdas_imerg.csv`) | ✔ 2026-10-08 |
| `utils_mapas.py` | Cuadrícula en grados, barra de escala y flecha de norte | ✔ |
| `09_datos_camels_pe_el_tambo.py` | Verificación MD5 del DOI, extracción del subconjunto, series mensuales y disponibilidad año × mes. **Relleno del caudal** con `relleno_caudal.py`: columnas observadas `R_obs_mm_mes`, `R_flag_obs`, `Q_obs_*`; columnas rellenadas `R_mm_mes`, `Q_mm_d`, `Q_m3_s`, `R_flag` (completo / interpolado_diario / escalado / modelo_mensual), `R_modelo`, `dias_interpolados_Q`. Serie diaria interpolada en `caudal_diario_interpolado.csv` | ✔ validado (revisión manual de 2010-03 y 1998-03); relleno 2026-10-10 |
| `relleno_caudal.py` | Módulo del método de relleno (interpolación de ln q en huecos ≤ 15 días; escalado si ≥ 50 % de días; modelo log-lineal con flow_sim y R del mes anterior) | ✔ 2026-10-10 |
| `28_validacion_relleno_caudal.py` | Validación del relleno: huecos sintéticos (600 por longitud) y validación cruzada dejando un año fuera de 5 modelos mensuales → `validacion_relleno_huecos_diarios.csv`, `validacion_relleno_modelo_mensual.csv`, `resumen_relleno_caudal.csv`, `fig_validacion_relleno_caudal.png` | ✔ 2026-10-10 |
| `11_grafica_precipitacion_camels.py` | Serie de P_L en el periodo común 1998-01 a 2024-12 (mensual con media móvil de 12 meses, y anual) → `fig_precipitacion_camels_1998_2024.png`. Media 833 mm/año (1998–2024) | ✔ 2026-10-08 |
| `10_imerg_mensual_y_grafica.py` | Lectura de IMERG (manual o automática), promedio ponderado, mm/h → mm/mes, gráfica P_L vs P_I | ✔ ejecutado 2026-10-07 (332 meses; corregida la verificación de malla float32). P_I 45.3 vs P_L 70.3 mm/mes, sesgo −24.9, MAE 29.1, RMSE 48.5, r 0.908 |
| `13_piscop_senamhi_verificacion.py` | Descarga PISCOp oficial (v2.1 update y v3.0, figshare) y compara su promedio de cuenca con `prec` de CAMELS-PE (2 polígonos × 2 métodos de pesos × 3 periodos) → `resultados/precipitacion/verificacion_piscop_senamhi_vs_camels.csv` | ✔ 2026-10-08 |
| `22_visor_html.py` | Genera `visor_resultados.html` (raíz): archivo autocontenido (~11 MB, 36 figuras JPEG incrustadas, 7 tablas, sin dependencias externas); cifras leídas de los CSV; tema claro/oscuro; ampliación al clic. **Desde 2026-10-11:** una sección a la vez (pestañas con enlace `#id` y botones anterior/siguiente): Cuenca · Contexto físico · 1.1 Series (incluye temperatura PISCOt vs ERA5-Land) · 1.2 IMERG · 1.3 Distribución · 1.4 Calidad · 1.5a Climatología · 1.5b Ciclo anual · 1.5c Síntesis · 2. Relaciones · 3. Tendencias · 4. Fourier · 5. Clima global (en construcción) · Datos. Ancho máx. 1600 px; en pantallas ≥ 1200 px la ficha va a la derecha de la figura. **Cada figura lleva una ficha "Qué muestra / Cómo se construyó / Cómo leerla"**: el texto de construcción está en el diccionario `COMO` (clave = nombre del PNG) y el script se detiene si falta una entrada. Subtítulo: "Tarea 1 · Hidrología 2026-2S · Universidad Nacional de Colombia - Sede Medellín"; 6 datos clave en una fila. Cuenca: lista de 17 parámetros morfométricos (`parametros_morfometricos.csv`) bajo el mapa | ✔ 2026-10-11 |
| `24_era5_vientos_humedad.py` | Descarga ERA5 mensual por CDS, 1998–2025, 0.5°, dominio 20°N–40°S × 110°O–20°O (Atlántico tropical, Amazonía, Andes, Pacífico oriental con Niño 1+2): u, v, q en 850/700/500/200 hPa; IVT este/norte y tcwv; presión en superficie para enmascarar los niveles bajo el terreno del modelo. Guarda una climatología estacional a 1° (`datos/era5/era5_viento_ivt_climatologia_estacional.nc`, 1.4 MB, versionada); los NetCDF mensuales (~400 MB) quedan fuera de Git | ✔ 2026-10-10 |
| `30_relaciones_2_1.py` | Punto 2.1: dispersión P_I vs P_L (1:1), P_L vs Q̄ (log), P_I vs Q̄, P_L vs R (R = P), coloreada por estación; Pearson y Spearman (originales y anomalías; sensibilidad con solo Q observado completo); sesgo, MAE y RMSE de IMERG por estación; rezagos y lluvia acumulada k = 0..6 (`resultados/relaciones/`) | ✔ 2026-10-11 |
| `31_modelos_2_2_2_3.py` | Puntos 2.2–2.3: modelos P_L ← P_I (lin, prop, factor mensual, potencial) y Q̄ ← P_L / P_I (lineal, lineal con rezagos, log con rezagos, log con intercepto mensual, climatología + anomalías; K = 1–4). Selección por validación cruzada dejando un año fuera **dentro** del ajuste. Evaluación por bloques: A ajuste 1998–2015 / evaluación 2016–2024; B ajuste 2007–2024 / evaluación 1998–2006. Referencias: IMERG sin corregir y climatología del caudal con los años de ajuste. Residuos vs tiempo, valor estimado y mes (`resultados/modelos/`) | ✔ 2026-10-11 |
| `32_tendencias_3.py` | Punto 3 completo: registros y fuentes (3.1), anomalías a y z con climatología fija 1998–2020 (3.2), OLS-HAC (con y sin efectos de mes), MK estacional/Hirsch–Slack, MK Hamed–Rao + Sen con IC de Gilbert, LOWESS con banda por bloques y modelo aditivo, segmentada + Pettitt; todos los meses y 12 subseries; registro completo y periodo común; FDR (BH); sensibilidad (ventanas inicio/fin, sin 1998/2017/2023, ±5 años, fuente: PISCOp v3, CHIRPS, Q solo observado); Niño 1+2 y balance por año hidrológico (3.6) → `resultados/tendencias/` | ✔ 2026-10-11 |
| `33_espectros_4.py` | Punto 4 completo: periodograma FFT (rect.), Hann y Welch (120 m, 50 %, ν calculado); densidad unilateral con Parseval verificado y contraste con `scipy.signal.periodogram`; series centradas, anomalías (ref. 1998–2020) y ambas sin tendencia; registro completo, periodo común y 1998–2020 (T); fondo AR(1) con umbral 95 % puntual y global (Bonferroni sobre ⌊N/2⌋ frecuencias); tabla de picos con resolución y ciclos; fracción de varianza por bandas; estabilidad (ventana, tendencia, extremos P1–P99, ln Q, mitades, periodo común); efecto de los vacíos de Q con Lomb–Scargle → `resultados/espectros/` | ✔ 2026-10-11 |
| `34a_descarga_perfiles.py`, `34b_perfiles_inversion.py`, `35_temporadas_humedad_inversion.py`, `estilo_figuras.py` | **Traídos del repositorio de la compañera** (github.com/anajul1a4/Hidrologia-Chicama-El-Tambo) el 2026-10-12, sin cambios. 34a: radiosondeos IGRA2 PEM00084516 Trujillo (12 UTC, desde 2024-10; `datos/radiosondeo/`) y ERA5 mensual por niveles a las 12 UTC del periodo común (`datos/era5_perfiles/`, CDS). 34b: inversión por sondeo (primer aumento continuo de T ≥ 1 °C bajo 3 500 m), perfiles compuestos por temporada, sondeo típico de la seca y contraste con ERA5 → `resultados/inversion/`. 35: IVT (script 24) + perfil de temperatura por temporada. Resultado: 645 sondeos 2024-10..2026-10, inversión en 91 %; seca +5.0 °C, 814–1188 m, 95 %; húmeda +2.4 °C, 338–717 m, 85 %. **Los sondeos quedan fuera del periodo común** (no hay anteriores); se contrastan con ERA5 12 UTC 1998–2024 (decisión del grupo, 2026-10-12). **Único cambio en 34b:** guarda `perfiles_por_sondeo.nc` (perfiles de cada sondeo en la malla de 25 m) para el 34c; sus figuras no cambian | ✔ 2026-10-12 |
| `34c_contraste_sondeos_era5.py` | (Propio, 2026-10-12, pedido del grupo.) Descarga ERA5 12 UTC 2025-01..2026-09 (`ERA5_perfiles_12UTC_mensual_2025_2026.nc`; 2024-10..12 del archivo del 34a) y separa dos comparaciones: **1** Trujillo (media de los sondeos) vs ERA5 en los mismos meses 2024-10..2026-09 → error de ERA5; **2** ERA5 en esos meses vs ERA5 en el periodo común → ¿fueron típicos? Estabilidad T(1500 m) − T(500 m): may–sep sondeos +1.5 °C, ERA5 mismos meses −0.2, ERA5 periodo común 0.0; dic–abr −2.5 / −3.6 / −3.4. ERA5 es hasta ~2 °C más cálido bajo ~800 m y no forma la inversión, pero sí la capa casi isoterma de may–sep; may–sep 2025–2026 fue 0.5–1 °C más cálido que lo normal con la misma forma del perfil → `fig_contraste_1_trujillo_vs_era5.png`, `fig_contraste_2_era5_periodos.png`, `contraste_sondeos_era5.csv`. Reemplaza en el visor a `fig_perfiles_inversion.png` (34b, que se sigue generando) | ✔ 2026-10-12 |
| `33_espectros_4.py` (revisión 2026-10-12) y `38_onditas_4.py` | **Punto 4 combinado con el visor de la compañera (decisión del usuario):** los cálculos del 33 son idénticos en ambos repositorios. Cambios en el 33: paleta común (`estilo_figuras.py`); misma escala vertical en las dos columnas de cada variable; **eje doble** en las figuras de espectros, fuentes y vacíos (inferior: periodo en años, creciente a la derecha; superior: frecuencia en ciclos/mes, como pide la guía); figura nueva `fig_espectros_porcentaje_4_2.png` (% de la varianza por frecuencia de Fourier, eje lineal en años) con rótulos de precisión honesta: periodo ≈ T con su rango ±½ bin, % y ciclos k, solo picos con k ≥ 3; * = 95 % puntual, ** = también global. Se descartó su par de figuras período/frecuencia separadas (cada una perdía la otra referencia). `38_onditas_4.py` traído sin cambios de su repositorio: Morlet ω0 = 6 (Torrence y Compo 1998) con corrección de sesgo (Liu et al. 2007), AR(1) 95 %, cono de influencia; banda 2–7 años vs Niño 1+2/3.4 (ERSST v5) → `fig_onditas_4_3.png`, `fig_onditas_banda_4_3.png`, `onditas_resumen_4_3.csv`, `onditas_banda_4_3.csv`. r de la componente 2–7 años con Niño 1+2: P_L 0.59, Q 0.46, T ERA5-Land 0.76 (0.86 con Niño 3.4); la banda reúne 12 % / 7 % / 35 % de las anomalías (descriptivo, sin p). En el visor, onditas como subsección «Complemento» al final del punto 4 (el usuario decidió mantenerlas aunque la guía no las pide) | ✔ 2026-10-12 |
| `37a_descarga_campos_globales.py`, `37_correlaciones_5.py` (de la compañera) y `37b_descarga_q_z.py` (propio) | **Punto 5, en evaluación (2026-10-12).** 37a/37 traídos de su repositorio y ejecutados con nuestros datos: resultados **idénticos** a los suyos (máx. |dif| = 0 en resumen_mapas, índices, robustez y todos_meses). Método: correlación entre años (27 eneros… 1998–2024) de anomalías (ref. 1998–2020) cuenca–celda, malla 2°, n_ef por persistencia entre años, FDR Benjamini–Hochberg α = 0.10 por mapa (índices: q < 0.05), robustez por correlación de patrones (Spearman, sin tendencia, sin 1998/2017/2023, IMERG, ℓ = 1 para Q, subperiodos P_L 1981–2002 vs 2003–2024). Campos de ella: SST (ERSST v5), SLP y u200 (ERA5). **Se agregó a su 37 un argumento `alternativo`** (sin argumento hace lo de ella): SST + q 700 hPa (g/kg, enmascarada donde sp < 700 hPa: 0.85 % de celdas-mes) + altura geopotencial Z = z/g0 en 500 hPa (37b descarga q y z en 500/700 y sp, 2°). Salidas en `resultados/correlaciones/alternativo_q_z/`. Hallazgos: Z500 da áreas FDR grandes (P_L dic 40 %, sep 17 %) pero es una sola banda tropical uniforme, no regiones distintas; q700 en ene–mar muestra un dipolo (más humedad en el Pacífico este, menos al norte del ecuador) y el índice «Costa y Pacífico este (q 700)» pasa FDR con P_L en ene, feb y nov (r = 0.77) y con Q en ene. Correlación de patrones mediana con el mapa de SST: SLP −0.27, u200 −0.21, q700 +0.30, Z500 +0.26 (ninguno redundante). **Decisión del usuario (2026-10-12): conjunto final = método de la compañera + SST + presión al nivel del mar + q 700 hPa; índices E y C; FDR q < 0.10 también en índices.** El 37 ahora toma `final` (por defecto, salidas en `resultados/correlaciones/`), `ella` y `alternativo` (en `correlaciones/comparacion_campos/`). E y C reconstruidos (EOF ERSST v5 10°S–10°N, 140°E–80°O, 1981–2025; no hay serie oficial del IGP publicada): r(E, Niño 1+2) = 0.92, r(C, Niño 4) = 0.97, r(E, C) = 0.00 → `indices_E_C_5.csv`. Resultados: E es la asociación más fuerte (Q r 0.73–0.75 ene–mar, FDR ene–may; P_L 0.52–0.64 ene–abr, FDR ene–abr); C con signo opuesto en P_L ene–abr (r −0.52 a −0.27; Spearman −0.59 a −0.29; no pasa FDR); Niño 3.4 ≈ 0 en ene–abr (hipótesis de cancelación E/C por altitud, no demostrada). Visor: sección p5 con 12 figuras, tabla de comparación de los tres conjuntos e **«Interpretación física» (propuesta para revisión del grupo) con referencias verificadas** | ✔ 2026-10-12 |
| `29_era5_inversion_termica.py` | Inversión térmica y subsidencia (ERA5 mensual, 2°S–12°S × 88°O–72°O, 0.25°): T y ω en 23 niveles (1000–200 hPa; descarga en bloques de 7 años por el límite del CDS), presión en superficie, T2m y relieve del modelo. Produce: corte oeste–este a 7–8°S (θ, ω, capas de inversión) en dic–abr y may–sep; ciclo anual de la LTS costera (θ700 − θsup, 83–80°O × 9–6.5°S), de ω500 sobre la cuenca, de la convergencia del IVT y de P_L; correlación de anomalías en dic–abr (`fig_inversion_termica_era5.png`, `ciclo_inversion_convergencia.csv`, `correlaciones_inversion_lluvia.csv`; `datos/era5_inversion/climatologia_inversion.nc` versionada, NetCDF de ~77 MB fuera de Git) | ✔ 2026-10-11 |
| `25_mapa_transporte_humedad.py` | Mapa regional de IVT y viento en 500 hPa para dic–abr y may–sep (`resultados/contexto_fisico/mapa_transporte_humedad_era5.png`). Contrasta el esquema de rutas de humedad del recuadro del mapa de la cuenca (script 06) | ✔ 2026-10-10 |
| `23_sintesis_ciclo_1_5c.py` | Punto 1.5.c (síntesis del ciclo anual; la explicación física con bibliografía queda en pausa): aridez UNEP P/PET con 3 fuentes de P × 2 de PET (Penman–Monteith de CAMELS-PE `pet`, solo hasta 2016-12; Hargreaves–Samani desde Tmin/Tmax de PISCOt hasta 2020-12) y contraste con Budyko; Köppen–Geiger (Peel et al. 2007); clasificación de cada mes (FAO P vs PET, Gaussen, aguas altas, R > P, frecuencia de picos); láminas por estación; tabla de criterios (`criterios_clasificacion_1_5c.csv`, `aridez_budyko_1_5c.csv`, `koppen_1_5c.csv`, `clasificacion_mensual_1_5c.csv`, `laminas_estacionales_1_5c.csv`, `fig_sintesis_ciclo_1_5c.png`) | ✔ 2026-10-10 |
| `21_ciclo_anual_1_5b.py` | Punto 1.5.b: indicadores de estacionalidad (SI, PCI, centroide de Markham, picos, armónicos, meses húmedos y secos de Gaussen), desfases, variabilidad por mes (CV, G1, jackknife), mes del pico por año hidrológico, subperiodos y anomalías por año hidrológico y calendario (`resultados/climatologia/`) | ✔ 2026-10-08 |
| `20_climatologia_1_5a.py` | Punto 1.5.a: climatología de 12 meses en el periodo común (n años, media, mediana, desv. est., Q1, Q3, P10, P90 por mes) de P_L, P_I, Q̄, R y T (PISCOt + estimación desde 2021); ciclo con bandas IQR y P10–P90; mapas año × mes (`resultados/climatologia/`) | ✔ 2026-10-08 |
| `19_control_calidad_1_4.py` | Punto 1.4: mapa y tablas de disponibilidad año × mes, comprobaciones diarias (fechas, códigos, negativos, tmean, flow_sim), secuencias constantes, tramos lineales, picos aislados, casos 1998 y diciembres, revisión manual de IMERG y registro de 17 anomalías (`resultados/control_calidad/`) | ✔ 2026-10-08 |
| `18_estadisticos_1_3.py` | Punto 1.3: tablas de estadísticos (registro completo, meses comunes de las 5 variables y pares P_L–P_I), meses extremos con valores simultáneos, histogramas (Freedman–Diaconis, frecuencia relativa) y diagramas de caja (Tukey 1.5·IQR) | ✔ 2026-10-08 |
| `17_era5land_temperatura.py` | Descarga ERA5-Land mensual T2m 1981–2025 por CDS (doi:10.24381/cds.68d2bb30; credenciales en `~/.cdsapirc`, fuera del ZIP), promedio ponderado por área (30 celdas), comparación con `tmean` CAMELS-PE 1981–2020 y relieve del modelo (`fig_T_camels_vs_era5land.png`, `era5land_t2m_cuenca_mensual.csv`) | ✔ 2026-10-08 |
| `16_grafica_precipitacion_imerg.py` | Punto 1.2: serie de P_I 1998-01 a 2024-12 (mensual con media móvil de 12 meses y totales anuales), mismo formato y ejes que la figura de CAMELS-PE, color #e4572e (`fig_precipitacion_imerg_1998_2024.png`). Media 539 mm/año (1998–2024) | ✔ 2026-10-08 |
| `15_figura_series_1_1.py` | Punto 1.1: P_L, Q̄ (m³/s), R (mm/mes) y T en paneles alineados en el periodo común (1998-01 a 2024-12), meses rellenados marcados por método y meses escalados marcados (`fig_series_mensuales_1_1.png`, `extremos_series_mensuales_1_1.csv`). R > P_L en 98/305 meses (sobre todo may–ago; 1998-03: R 837 vs P 414 mm) | ✔ 2026-10-08 |
| `14_piscop_v3_cuenca.py` | Recorte de PISCOp v3.0 diario, pesos por área (GLO-30), controles de calidad, series diaria y mensual de cuenca (`piscop_v3_cuenca_{diaria,mensual}.csv`) | ✔ 2026-10-08 |
| `08_descarga_piscop.py` | Descarga de PISCOp desde IRI | ✗ descartado (no usar) |
| `01_delimitacion.py`, `03.py` (vacíos), `02_Prueba.py` | Archivos previos (CAMELS-DE) | fuera del flujo |

**Orden de ejecución:** 06 → 07a → 13 → 14 → 09 → 10 → 28 → 26 → 07b → 11 → 12 → 17 → 15 → 16 → 18 → 19 → 20 → 21 → 23 → 24 → 25 → 29 → 27 → 30 → 31 → 32 → 33 → 38 → 34a → 34b → 34c → 35 → 22 (visor HTML, siempre al final). El periodo común
(`periodo_comun.py`) necesita 09 y 10, así que 07b, 11 y 15–19 van después; 15 requiere además la extensión de T del 17.

**Productos principales:**
- `Tarea 1/El_Tambo/resultados/delimitacion/`: polígono, red, río principal, `parametros_morfometricos.csv`, `mapa_cuenca_el_tambo.png`.
- `Tarea 1/El_Tambo/resultados/contexto_fisico/`: `panel_contexto_fisico.png`, mapas individuales, `ficha_contexto_fisico.csv`.
- `Tarea 1/El_Tambo/resultados/series_mensuales/`: `series_mensuales_camels_pe*.csv`, `disponibilidad_anio_mes.csv`, `fig_precipitacion_camels_1998_2024.png` (el CSV `_1998_2025` es la serie completa; cada script la recorta con `periodo_comun.py`).

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
- R_m = Σ q_d (q ya está en mm/d; no se divide otra vez por el área). Criterio observado (columnas `*_obs`): mes completo → suma directa; ≤ 10 % de días faltantes → suma de los días válidos escalada a n_días; > 10 % → excluido. **La serie de análisis es la rellenada** (ver "Relleno del caudal" más abajo). No se presentan sumas parciales como acumulados completos.
- Q en m³/s = q[mm/d] · A / 86.4, con A = 2171.585 km² (área CAMELS-PE).
- T mensual = media de días válidos (≤ 10 % faltantes).
- P_I: tasa mm/h × 24 × días del mes.
- **Relleno de vacíos.** P_L no tiene vacíos. Se decidió **no rellenar P_L con IMERG corregido**, porque contaminaría la comparación P_I–P_L y los modelos del punto 2. Para T 2021–2025 se usa ERA5-Land evaluado contra PISCOt en 1998–2020, sin unir ambas fuentes en una sola serie.
- **Relleno del caudal (decisión del grupo, 2026-10-10: "dado el alcance del trabajo, se aceptan afectaciones en la calidad").** Antes, los meses con > 10 % de días faltantes se excluían. Ahora se aplica, en orden (`relleno_caudal.py`, aplicado en el 09, validado en el 28):
  1. **Huecos diarios de hasta 15 días** con datos a ambos lados: interpolación lineal de ln(q), que sigue la forma de la recesión. Validación con 600 huecos sintéticos por longitud en meses completos: error mediano del total mensual 0.1–1.5 % hasta 7 días y 2–4 % con 10–15 días; P90 ≤ 6 % hasta 7 días y 20 % con 15. El escalado anterior daba 3–8 veces más error.
  2. **Meses con faltantes restantes y datos en ≥ 50 % de los días:** media de los días disponibles × n_días. Error mediano ~10 % con 15 días faltantes.
  3. **Meses con < 50 % de días:** ln R_m = a_mes + b·ln S_m + c·ln R_(m−1), con S = `flow_sim` de CAMELS-PE (PISCO-ARNOVIC v1.1); si falta R_(m−1), solo S. Corrección log-normal +σ²/2. Validación cruzada dejando un año fuera: NSE 0.72, NSE(log) 0.87, error mediano 23 %, sesgo +1.6 %. Mejor que P y P(m−1) solos (NSE 0.58) y que la razón obs/sim (NSE 0.54).

  **Resultado en el periodo común:** 80 meses tocados (262 días interpolados): 53 antes escalados y 23 antes excluidos pasan a interpolado_diario; 3 excluidos quedan escalados (2019-01, 2019-02, 1998-09); 1998-08 (0 días) se reconstruye con el modelo (8.1 mm). Los meses antes excluidos suman 529 mm de R (4.9 % del total). 1998-02 queda en 385 mm y 1998-03 pasa de 837 (escalado) a 852 mm.

  **Implicaciones a declarar:** el relleno mezcla información del modelo y de la lluvia (vía flow_sim) solo en 1998-08. En la relación lluvia–caudal del punto 2 conviene repetir el análisis sin los meses con `R_flag` ≠ completo como prueba de sensibilidad. En 2019–2021 hay más meses rellenados, con mayor incertidumbre. La disponibilidad del 1.4 sigue mostrando lo observado.

**Contexto físico**
- Litología agrupada en 6 grupos hidrogeológicos más lagunas (asignación en `unidades_geologicas_cuenca.csv`); pendiente de contrastar con INGEMMET (Serie H) y ANA.
- Hipótesis de humedad: fuente principal amazónica/atlántica por flujo del este en niveles medios y altos (dic–abr); Pacífico secundario, relevante en El Niño. Es una hipótesis a contrastar en el punto 5 (flujo de humedad de ERA5, SST Niño 1+2).
- **Mapa de la cuenca (2026-10-10, pedido del grupo):**
  - Se quitaron las flechas de humedad del mapa principal.
  - **Recuadro de Sudamérica:** el 2026-10-10 se le agregó el IVT de ERA5 de dic–abr; **el 2026-10-11 el grupo decidió quitarlo**. El mapa de la cuenca vuelve a tener solo ubicación (Sudamérica y Perú con Lima), en disposición horizontal. El transporte de humedad queda en mapa_transporte_humedad_era5.png (script 25). También se eliminó la figura comparacion_camels.png (el IoU 0.994 sigue en los parámetros morfométricos).
  - **Por qué dic–abr:** en esa estación cae el 81 % de la lluvia. En may–sep el IVT sobre la cuenca es mayor, pero no llueve, así que mostrarlo en el recuadro induciría a error. La estación seca y el viento en 500 hPa quedan en la figura del script 25.
  - Un esquema anterior de 4 rutas según bibliografía (Atlántico → Amazonía; cruce de los Andes desde el este por la depresión de Huancabamba; chorro de bajo nivel; Pacífico en El Niño costero) se reemplazó por los datos. Sigue siendo útil para la explicación física de 1.5.c.
  - **Verificadas en la web el 2026-10-12** (existencia, revista, DOI y resumen; el texto completo lo debe leer el grupo): Takahashi et al. 2011 GRL 38 L10704 (índices E y C); Lavado-Casimiro y Espinoza 2014 Rev. Bras. Meteorol. 29(2) 171–182 doi:10.1590/S0102-77862014000200003 (la cita «2012» era errónea); Sulca et al. 2018 IJC 38 420–435 doi:10.1002/joc.5185; Takahashi et al. 2018 BAMS 99(8) S210–S211 (El Niño costero 2017); Garreaud 2018 IJC 38 e1296–e1302 doi:10.1002/joc.5426; Garreaud 2009 Adv. Geosci. 22 3–11 doi:10.5194/adgeo-22-3-2009; Takahashi y Battisti 2007a J. Climate 20 3434–3451 y 2007b 20 5696–5706 (Andes, ZCIT al norte y zona seca); Rau et al. 2017 IJC 37 143–158 doi:10.1002/joc.4693; Espinoza et al. 2020 Front. Earth Sci. 8:64; Wilks 2016 BAMS 97 2263–2273 (FDR; la regla α_FDR = 2α falta confirmarla en el texto); Bretherton et al. 1999 J. Climate 12 1990–2009 (n efectivo). Lista original, **sin verificar**, conservada abajo: Garreaud 2009 (Adv. Geosci.); Garreaud et al. 2009 (Palaeo3); Vera et al. 2006 y Marengo et al. 2004 (J. Climate); Espinoza et al. 2020 (Front. Earth Sci.); Segura et al. 2019 (Clim. Dyn.); Garreaud 2018 y Takahashi et al. 2018 (El Niño costero 2017); Sulca et al. 2018, Lavado-Casimiro et al. 2012 y Rau et al. 2017 (IJC); McPhaden 2006 (en `Bibliografia/`).
- **Inversión térmica y subsidencia (script 29, 2026-10-11).** Hipótesis: el anticiclón del Pacífico Sur genera subsidencia; el aire que desciende se calienta por compresión y, sobre la capa marina fría (afloramiento), forma una inversión que actúa como tapa. Por continuidad de masa, la subsidencia implica divergencia en niveles bajos. Resultados con ERA5 en el periodo común:
  - **Corte vertical a 7–8°S:** en may–sep, inversión extendida desde 88°O hasta la costa a ~850–925 hPa (~1–1.5 km) y descenso sobre el Pacífico y la vertiente occidental. En dic–abr, la inversión queda reducida a parches y hay ascenso sobre los Andes, incluida la cuenca.
  - **Ciclo anual:**
    - LTS costera mínima en feb (17.0 K) y máxima en sep (24.2 K).
    - ω500 sobre la cuenca: descenso en jun–ago (+0.025 a +0.039 Pa/s), ascenso máximo en mar (−0.087).
    - La convergencia del IVT sigue a ω500.
  - **Matiz importante:** en sep–nov la LTS sigue alta (24.2 → 22.6 K) mientras sobre la cuenca ya hay ascenso y convergencia. La inversión costera explica la supresión de lluvia en la estación seca, pero el inicio de las lluvias lo controla el ascenso sobre los Andes (calentamiento continental, monzón, Alta de Bolivia).
  - **Anomalías en dic–abr (n = 139):**
    - r(LTS, P_L) = −0.47 (n efectivo ≈ 108): inversión más débil, más lluvia.
    - r(ω500, P_L) = −0.78.
    - r(convergencia, P_L) = +0.83.
    - Feb–Mar 1998 y Mar 2017 (El Niño costero): LTS 1.2–2.7 K por debajo de lo normal y lluvia muy por encima. Feb 2017 tuvo la inversión débil sin exceso de lluvia, así que la inversión débil no basta por sí sola.
  - **Limitaciones:** ERA5 suaviza los Andes (~30 km). ω y la convergencia vienen del mismo reanálisis, así que su relación con P_L es consistencia física, no prueba de causalidad; la LTS es el indicador más independiente. Referencias a verificar: Klein & Hartmann 1993 (LTS), Wood & Bretherton 2006 (EIS), Garreaud 2009.
- **Contraste con ERA5 (scripts 24–25, climatología 1998-01..2025-08):**
  - **Estación húmeda:** el IVT entra desde el Atlántico, cruza la Amazonía hacia el oeste y gira al sur pegado a los Andes (chorro de bajo nivel).
  - **Sobre la cuenca** (celda de 1° en −7.5°, −78.8°): el IVT va hacia el oeste en ambas estaciones (270° en dic–abr; 264° en may–sep), con 33 y 50 kg m⁻¹ s⁻¹. El agua precipitable es mayor en dic–abr: 28 vs 22 kg/m².
  - **Viento en 500 hPa:** del este todo el año, 3.0 m/s en dic–abr y 5.3 en may–sep.
  - **Frente a la costa:** el IVT se aleja del continente, lo que es coherente con un aporte del Pacífico bajo en años normales.
  - **Lectura:** la magnitud del IVT no explica la estacionalidad (es mayor en la estación seca). Lo que la explica sería la convergencia y el ascenso, más humedad disponible en dic–abr. Para el punto 5 conviene calcular la convergencia del IVT y las anomalías de los años El Niño.
  - **Limitaciones:** ERA5 a 0.25° (promediado a 1°) no resuelve los valles. En 850 hPa el nivel queda bajo el terreno del modelo sobre la cuenca, y esas celdas se enmascaran. Para la cuenca se descartó ERA5-Land, porque no tiene viento en altura ni IVT y su viento de 10 m es el de ERA5 interpolado.

---

- **Decisión del grupo (2026-10-08): P_L principal = `prec` de CAMELS-PE; PISCOp v3.0 = fuente de apoyo** (mapa continuo (e) y contraste de incertidumbre). Razones: coherencia con la base elegida (Q, T, atributos), trazabilidad de la entrega y balance hídrico menos inverosímil con la P más alta (R/P 0.48 vs 0.56; ambas por debajo de lo esperado por Budyko). P v3.0 1998–2025 = 720 mm/año vs 841 CAMELS-PE (r mensual 0.983); cociente anual 0.8–0.9, menor en años muy lluviosos (1983: 0.62; 1998: 0.74). Scripts 09, 10 y 11 siguen con CAMELS-PE. No se agrega comparación de las tres lluvias (decisión del grupo). En el periodo común 1998-01 a 2025-08 (Σ medias mensuales): PISCOp v3.0 721 vs CAMELS-PE 841 mm/año.

## 4. Pendientes inmediatos
0. **Figuras tentativas (aviso del grupo, 2026-10-11):** `fig_inversion_termica_era5.png` (script 29) y `fig_sintesis_ciclo_1_5c.png` (script 23) **no convencen al grupo**: son parciales, no finales. Revisarlas antes de la entrega; por ahora no se modifican.
1. Interpretar el sesgo de IMERG (P_I/P_L = 0.645; subestima en dic–abr, sobreestima en jun–ago) y declarar que P_L no es una observación independiente.
2. Temperatura: ERA5-Land ya descargado (1981–2025) y comparado (script 17). **Resultado en el periodo común (meses con PISCOt: 1998-01 a 2020-12, 276 meses):** misma variabilidad (r series 0.89, r anomalías 0.90; amplitud del ciclo anual 0.57 vs 0.60 °C, ambos con máximo ago–sep) pero ERA5-Land **2.84 °C más frío**, desfase casi constante en todos los meses. (Con 1981–2020 era −2.92 °C.) Relieve del modelo 2641 m vs 2557 m reales (+84 m) → explica solo ≈ −0.55 °C. Otras causas a evaluar: `tmean` de PISCOt = (Tmax+Tmin)/2 vs media de 24 h de ERA5-Land; sesgo del reanálisis en los Andes; interpolación de PISCOt. Verificar con normales de estaciones SENAMHI cercanas cuál nivel es más realista. Para anomalías, tendencias y correlaciones el desfase constante no pesa; para valores absolutos (PET, aridez) sí. **Extensión de T desde 2021 (2026-10-08, pedido del grupo):** T_est = ERA5-Land + sesgo medio del mes calendario frente a PISCOt en el periodo común 1998–2020 (+2.62 a +2.99 °C); 56 meses estimados, 2021-01 a 2025-08. Validación fuera de muestra (sesgo ajustado en 1998–2010, evaluado en 2011–2020): MAE 0.19 °C, sesgo −0.05, r 0.89; climatología como referencia: MAE 0.40 °C (`validacion_extension_T.csv`). Archivo `temperatura_piscot_extendida_era5land.csv` (columna aparte y bandera de fuente; `T_C` original intacta). Se muestra en la fig. 1.1 como tramo discontinuo morado. En las tendencias del punto 3, analizar PISCOt 1981–2020 por separado y tratar 2021–2025 como estimado (riesgo de salto por cambio de fuente). Enfoque acordado: **no unir fuentes**. PISCOt v1.2 (CAMELS-PE, 1981–2020; no existe versión posterior) para tendencias de largo plazo; ERA5-Land T2m mensual 1998–2025 como serie completa para el periodo común (1.5, 4, 5), evaluada contra PISCOt en 1998–2020 (sesgo por altitud: ERA5-Land ~9 km vs relieve 684–>4000 m). Requiere cuenta CDS y `~/.cdsapirc` (fuera del ZIP).
3. Puntos 1.1 a 1.4, 1.5.a, 1.5.b y la evidencia de la síntesis de 1.5.c hechos (scripts 15, 10, 16, 18, 19, 20, 21 y 23). Sigue: el grupo elige y redacta la clasificación (notas 1.5.c en la sección 5); explicación física de 1.5.c con bibliografía (**en pausa por decisión del grupo, 2026-10-10**) e hipótesis inicial (punto 1.4).
4. Resolver la diferencia PISCOp (IRI, v2.1 update y v3.0 oficiales, todos 15–25 % menores) vs CAMELS-PE con los autores (hllauca@senamhi.gob.pe) o declararla como limitación.
5. Decidir si se borra `08_descarga_piscop.py` antes de armar el ZIP. El ZIP no debe incluir `PISCOp_d.nc` (1.5 GB); basta el recorte de la cuenca.
6. Armar el ZIP de entrega (`requirements.txt` regenerado con `pip freeze` el 2026-10-11, incluye statsmodels y pymannkendall).
7. Completar los integrantes en `README.md`.
8. Borrar las carpetas sobrantes `hidrologia-tarea1-el-tambo/` y `Hidrologia_Rio_Chicama_El_Tambo/` (repositorios vacíos creados por error con GitHub Desktop; el `.gitignore` ya las excluye).
9. Puntos siguientes de la guía: 1.5.c e hipótesis inicial (interpretación del grupo), **punto 5 (siguiente)**. Puntos 3 y 4 hechos (scripts 32 y 33): el grupo debe revisar las hipótesis de 3.6 y buscar bibliografía sobre usos del agua y cobertura en el Chicama.
10. **Correcciones del visor en curso (pedido del grupo, 2026-10-11/12):** hechas las 1–3 (secciones individuales, subtítulo, fila de datos clave), el mapa de la cuenca y la ficha de cada figura. **Contexto físico (2026-10-12):** panel en 1 × 4 con la ficha debajo (`figura(..., ancha=True)`, clase CSS `.fig.ancha`); se retiraron del visor `mapa_geologia_unidades_ingemmet.png` y `mapa_d_precipitacion.png` (07b los sigue generando; sus entradas `COMO` se borraron). **Después, la sección se igualó a la del visor de la compañera:** todas sus figuras con `ancha=True` (imagen centrada; ficha debajo con "Qué muestra" y "Cómo se construyó" en dos columnas y "Cómo leerla" en una fila); subsección "Fuentes de humedad e inversión térmica" con el sondeo típico (34b), la figura seca vs húmeda (35) y una subsección "¿Representan los sondeos el periodo común?" con las comparaciones 1 y 2 del 34c. Se retiraron del visor `mapa_transporte_humedad_era5.png` (25), `fig_inversion_termica_era5.png` (29) y `curva_hipsometrica.png`; los scripts siguen. Cifras de la ficha leídas de `inversion_por_sondeo.csv` y de la climatología del 24. **Layout del visor (2026-10-12, tomado del visor de la compañera):** ancho máximo 1200 px; cada figura es una tarjeta (título + chip «imagen»/«interactiva», figura a todo el ancho con botón de ampliar, ficha debajo en dos columnas + «Cómo leerla» en fila completa); subtítulo de sección justificado. **Gráficas interactivas (Plotly 2.35.2 incrustado desde `Codigos/recursos/plotly-2.35.2.min.js`; lógica en `Codigos/visor_graficas.js`, tomada de su `web/graficas.js`; datos `DATOS` en el 22, portados de su `36_pagina_web.py`):** 1.1 series alineadas y temperatura; 1.2 P_L vs IMERG (con CHIRPS y PISCOp ocultas) y totales anuales; 1.4 disponibilidad + serie diaria + panel lateral por variable; 1.5a climatología con bandas; 1.5b anomalías mensuales y anuales; 2 observado vs modelado; 3 anomalías con Sen. Reemplazan a los PNG de 15, 17, 10, 11, 16, 19 (disponibilidad), 20 (climatología) y 21 (anomalías); en 2 y 3 se agregan sin quitar figuras. Su texto «diferencia mayor en 1998, 2017, 2023» no coincide con nuestros datos (son 2023, 1998, 2021): ahora se calcula. Secciones 4 y 5: diferencias con su visor pendientes (decisión del usuario). El grupo puede seguir pidiendo ajustes por sección. Al agregar figuras (p. ej. punto 5), escribir su entrada en `COMO` verificando el script que la genera.
11. **Punto 5 (guía, sección 5):** campos mensuales globales de SST (catálogo NOAA PSL, p. ej. ERSST v5 u OISST) y dos variables atmosféricas justificadas (ERA5 mensual: p. ej. presión al nivel del mar y altura geopotencial de 500 hPa, o viento en 850/200 hPa, o humedad específica). 12 mapas por mes calendario para P_L y Q (y contraste con IMERG) frente a cada campo; Pearson (y Spearman); anomalías con referencia fija compatible (1998–2020 como en los puntos 3–4); escala −1 a 1; máscara de tierra en SST; mínimo de años; sin tendencia vs con tendencia; persistencia (n efectivo o bloques); FDR sobre la familia de celdas; estabilidad por subperiodos, años extremos y fuente; interpretación física con bibliografía. Script sugerido: `37_correlaciones_5.py` (+ `37a_descarga_campos_globales.py`), la numeración del repositorio de la compañera, que **ya tiene el punto 5 hecho** (y 36 página web, 38 onditas, 39 esquema, 40 informe PDF, 41 entrega): revisar con el grupo antes de rehacerlo. Hay que escribir las entradas `COMO` de sus figuras y llenar la sección `p5` del visor.

---

## 5. Notas para el informe

> **Aviso (2026-10-10): relleno del caudal.** Las cifras de Q y R de las notas de abajo se escribieron antes de rellenar el caudal. Las figuras y los CSV ya están regenerados (scripts 15, 18–21 y 23). Hay que actualizar en el texto, entre otras:
> - **R medio anual:** 380 → **390 mm** (Σ de medias mensuales; observado 380).
> - **R > P_L:** 98/305 → **104/332** meses.
> - **R de 1998-03:** 837 → **852 mm**.
> - **R:** SI 0.89 → **0.91**, PCI 19.1 → **18.8**, centroide de R 3.45 → **3.40**.
> - **Balance 1998–2016:** R/P_L 0.46 → **0.48**.
>
> Tomar los valores definitivos de los CSV de `resultados/`.

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

### Punto 1.2 – IMERG (figs. `fig_P_L_vs_IMERG_mensual.png` y `fig_precipitacion_imerg_1998_2024.png`)
- **Procedencia de P_L (lo exige la guía):** P_L = `prec` de CAMELS-PE = PISCOp v2.1 (interpolación de pluviómetros SENAMHI con CHIRP/CHIRPS y climatologías) promediado sobre la cuenca. Es un producto interpolado, **no una observación local independiente**. No es reproducible con las versiones públicas de PISCOp (IRI, v2.1 update, v3.0 dan 15–25 % menos): declararlo como limitación. IMERG Final también incorpora pluviómetros (GPCC): fuentes compartidas posibles.
- **Documentación de IMERG que pide la guía** (ya en la cabecera del script 10 y en la nota de la figura): producto GPM IMERG Final Precipitation L3 1 month 0.1° V07 (GPM_3IMERGM.07, doi:10.5067/GPM/IMERG/3B-MONTH/07); modalidad Final; variable `precipitation` en mm/h (tasa media mensual) → mm/mes = tasa · 24 · n_días; malla 0.1°; 1998-01 a 2025-08 (332 meses, sin faltantes ni duplicados); acceso manual NASA GES DISC (Subset / OPeNDAP, grupo `Grid`), con descarga automática alternativa por OPeNDAP + Earthdata; promedio ponderado por área de intersección celda–cuenca (UTM 17S), 31 celdas, 10 completas, 43.7 % del peso en bordes, ≈ 17.8 celdas equivalentes (mapa `mapa_celdas_imerg_pesos.png`); no se usa el píxel de la estación; los pesos no consideran la topografía (desnivel interno de las celdas: mediana 1233 m).
- **Periodo común P_L–P_I:** 1998-01 a 2025-08, 332 pares mensuales. P_I 45.3 vs P_L 70.3 mm/mes; sesgo −24.9 mm/mes; MAE 29.1; RMSE 48.5; r = 0.908.
- **Cifras para describir la figura de IMERG:** total anual medio 539 mm/año (1998–2024; 2025 incompleto excluido) vs 833 mm/año de P_L en los mismos años (1998–2024; 797 mm/año en 1981–2025); máximo mensual 265 mm (03/1998) vs 414 mm; 03/2017: 202 vs 377 mm; año más seco 2020 (370 mm), más lluvioso 1998 (785 mm). IMERG reproduce los eventos de El Niño 1998 y 2017 con picos mucho menores; media móvil más plana que la de P_L; desde 2019 casi todos los años por debajo de su media → contrastar en el punto 3 (tendencia real vs cambio en los insumos del producto).
- **Balance hídrico como contraste independiente:** con R medido (~402–407 mm/año), R/P = 0.47–0.48 con P_L y 0.72 con P_I; P − R implícito 462 vs 157 mm/año. Ambas P quedan por debajo de lo esperado por Budyko (PET CAMELS-PE 1341 mm/año) → la brecha de IMERG parece venir sobre todo de subestimación de IMERG, sin descartar sesgo de P_L. Hipótesis a sustentar con bibliografía (lluvia orográfica cálida mal detectada por microondas en los Andes).

### Punto 1.5.c – Síntesis del ciclo anual (script 23; `fig_sintesis_ciclo_1_5c.png`, `criterios_clasificacion_1_5c.csv`)
La explicación física con bibliografía está **en pausa** (decisión del grupo, 2026-10-10). Lo que sigue es evidencia cuantitativa y una **propuesta de trabajo**; la clasificación final la elige y redacta el grupo.

**Lo que pide la guía:**
- una clasificación hidroclimática mensual con criterios formulados;
- cada rasgo ligado a evidencia y a procesos;
- alternativas e incertidumbres;
- una hipótesis contrastable con el punto 2 (relaciones) y el punto 3 (tendencias);
- no copiar atributos de CAMELS;
- aridez con un criterio climático.

**Evidencia por rasgo (periodo común 1998-01..2025-08 salvo indicación):**
- **Forma del ciclo:**
  - Unimodal: un solo pico en marzo; el armónico anual explica el 78 % de la varianza del ciclo y el semianual el 15 %. Se mantiene en los dos subperiodos.
  - El máximo anual de P_L cae en feb–abr en el 93 % de los años hidrológicos. Excepciones: 2019-20 y 2023-24, con el máximo en diciembre.
- **Concentración:**
  - P_L: SI 0.79 (estacional, en el límite con «marcadamente estacional con estación seca larga», 0.80); en los subperiodos 0.77–0.80. PCI 15.4, entre las clases moderada y estacional (15–16). El 81 % de la lluvia cae en dic–abr.
  - IMERG es menos concentrada: SI 0.63, PCI 13.7, 73 % en dic–abr.
- **Estaciones:**
  - Húmeda dic–abr (meses ≥ media mensual); seca de Gaussen may–sep (P < 2T); oct–nov de transición.
  - Con P/PET (criterio FAO, PET Penman–Monteith 1998–2016):
    - P_L: solo feb–abr son húmedos (P ≥ PET); ene y dic intermedios; may–nov secos. Mayo es seco en el 84 % de los años (P_L < 0.5·PET) y jun–sep en todos.
    - IMERG: solo marzo es húmedo.
- **Aridez (UNEP P/PET):**

  | Fuente de P | AI (PET Penman–Monteith) | Clase |
  |---|---|---|
  | P_L | 0.62 | subhúmedo seco |
  | PISCOp v3 | 0.53 | subhúmedo seco |
  | IMERG | 0.42 | semiárido |

  Con Hargreaves–Samani los valores son casi iguales (PET 1349 vs 1344 mm/año). **La clase depende de la fuente de P, no de la de PET.** P_L queda cerca del límite de 0.65 (húmedo).
- **Köppen–Geiger (Peel et al. 2007):** Cwb con P_L y con PISCOp v3 (templado, invierno seco, verano templado). BSk con IMERG (MAP 543 < umbral 605 mm). Es una clase de cuenca agregada: T media 16.3 °C sobre 684–4317 m oculta el gradiente entre el valle bajo (más cálido y seco) y la puna.
- **Temperatura:** isotermal. Amplitud del ciclo 0.62 °C (máx sep, mín jun), del orden de la desviación estándar entre años de cada mes (0.36–0.69 °C). La estacionalidad la impone la lluvia, no la energía.
- **Respuesta del caudal:**
  - Unimodal, con máximo en abril (11 de 15 años con dato completo) y aguas altas (Q̄ ≥ media anual) feb–may.
  - Desfase de centroides R − P_L: 1.4 meses (1.5 con IMERG).
  - El caudal es más concentrado que la lluvia (SI de R 0.89 vs 0.79).
  - En la climatología de may–ago, R > P_L.
  - La estación seca (may–nov) aporta el 24 % de R con solo el 19 % de P_L (R/P_L 0.57 vs 0.43 en la húmeda). Es compatible con la liberación de agua almacenada, aunque mayo incluye la recesión de la estación húmeda.
- **Balance de largo plazo:** R/P observado 0.46 (P_L) y 0.69 (IMERG), frente a 0.16 y 0.07 esperados por Budyko con la PET de PISCOeo_pm. El caudal es mucho mayor que el que admite el balance con cualquier P. Explicaciones posibles:
  - P subestimada en las cabeceras, donde hay pocos pluviómetros;
  - PET de Penman–Monteith sobrestimada (no es la ET real);
  - aportes no representados (agua subterránea de fuera de la divisoria, trasvases, retornos de riego o minería);
  - un sesgo de Q.

  El grupo debe discutirlas. Ver también 1.4 (1998-03, R ≈ 2·P).

**Propuesta de clasificación para discutir (no adoptada):** régimen hidroclimático **pluvial tropical andino de vertiente pacífica**:
- **unimodal y estacional**: lluvia dic–abr con máximo en marzo, estación seca may–sep, oct–nov de transición;
- **isotermal**;
- **subhúmedo seco** según P/PET con P_L; semiárido con IMERG;
- caudal unimodal con máximo en abril, desfasado ~1.4 meses de la lluvia, y estiaje sostenido por almacenamiento.

**Alternativas e incertidumbres a declarar:**
1. La clase de aridez y la de Köppen cambian con la fuente de P (C ↔ B).
2. SI y PCI están en los límites entre clases: «estacional» vs «marcadamente estacional».
3. La PET de Penman–Monteith solo cubre 1998–2016.
4. T de cuenca agregada.
5. Balance R/P incompatible con Budyko.
6. Sin nieve ni glaciares relevantes (no hay evidencia de pico de fusión en el ciclo de Q); embalse minero según atributos (`reservoir_use` Mining, volumen sin dato); acreditaciones de agua de ANA: 2 superficiales (0.008 hm³/año) y 17 subterráneas (0.064 hm³/año), en total ≈ 0.03 mm/año sobre la cuenca, despreciable frente a R ≈ 380 mm/año. Solo cuenta los usos formales registrados.

**Hipótesis contrastables (borradores para que el grupo elija):**
- **H1 (punto 2):** la relación P–Q̄ mensual mejora con P desfasada 1 mes (o con P acumulada de 2–3 meses) frente a la simultánea, y es no lineal: umbral de humedecimiento en dic–ene antes de que responda el caudal. En anomalías, la asociación persiste pero es más débil.
- **H2 (punto 2):** IMERG frente a Q̄ tiene una pendiente mayor que P_L frente a Q̄, porque IMERG subestima los picos de la estación húmeda; la dispersión es mayor en mar–abr.
- **H3 (punto 3):** la fase del ciclo es estable (máximo en marzo, sin migración del mes del pico); los cambios entre años se concentran en el volumen de ene–abr. Una tendencia en P_I desde 2019 sin equivalente en P_L ni en R apuntaría a un cambio del producto IMERG, no del clima.
- **H4 (punto 3):** el caudal de estiaje (jul–oct) tiene más persistencia entre años que la lluvia de esos meses (memoria del almacenamiento).

### Punto 2 – Relaciones y modelos (scripts 30 y 31; periodo 1998-01..2024-12, 324 pares)

**2.1 Relaciones:**
- **(a) IMERG vs P_L:** r 0.90, ρ 0.93, pero no concordantes.
  - Sesgo −24.5 mm/mes (−35 %), MAE 28.7, RMSE 48.0.
  - Por estación: DEF −42 %, MAM −39 %, SON −18 %; en JJA IMERG sobreestima (+3.5 mm/mes, +70 % de casi nada).
  - Ejemplo de "correlación alta ≠ concordancia".
- **(b) P_L vs Q̄ del mismo mes:** débil (r 0.57, ρ 0.52), muy dispersa y heterocedástica.
- **(c) P_I vs Q̄:** r 0.50, ρ 0.46.
- **Pearson vs Spearman:** Pearson > Spearman en (b) y (c) por los meses extremos de El Niño (1998-03, 2017-03/04); Spearman > Pearson en (a) por la relación monótona no lineal.
- **(b') En lámina:** R > P_L sobre todo en may–ago (almacenamiento) y en crecidas de El Niño.
- **Rezagos (la lluvia antecede al caudal):**
  - ρ(Q̄_t, P_L,t−1) = 0.81 frente a 0.52 sin rezago.
  - Lluvia acumulada t−3..t: ρ 0.91. Es coherente con el desfase de ~1.4 meses del centroide (1.5.b) y con el almacenamiento.
  - P_I: k = 1 ρ 0.75; acumulada k = 3 ρ 0.87.
- **Anomalías:** la asociación se debilita pero persiste. P_L vs Q̄: ρ 0.37 sin rezago, 0.61 con k = 1, 0.62 acumulada; IMERG vs P_L: r 0.74. Parte de la relación original es el ciclo anual compartido, pero hay señal entre años.
- **Sensibilidad sin meses rellenados** (244 pares): ρ(P_L, Q̄) 0.54 frente a 0.52; conclusiones iguales.

**2.2–2.3 Modelos** (selección por validación cruzada dejando un año fuera dentro del ajuste; evaluación fuera del ajuste en dos particiones por bloques de años):
- **(i) P_L ← P_I:** se elige el **factor por mes calendario** P̂_L = b_j·P_I en ambas particiones.
  - Habilidad frente a IMERG sin corregir: +37 % (A, evaluación 2016–2024; MAE 19.2 vs 30.4, NSE 0.77 vs 0.49) y +40 % (B, evaluación 1998–2006).
  - La regresión lineal mejora menos (+25 %) y da 2 valores negativos.
  - **Relación aprovechable:** corrige sobre todo el sesgo estacional de IMERG.
- **(ii) Q̄ ← P_L:** se elige ln Q̄ = a_mes + Σ b_k·ln(P_L,t−k + 1), con K = 3 (A) o 2 (B).
  - Habilidad frente a la climatología: +41 % (A; MAE 9.5 vs 16.1, NSE 0.52 vs 0.21) y +26 % (B).
  - Sensibilidad solo con meses observados: +44 % (A) y +30 % (B).
  - La regresión lineal simple es **peor que la climatología** (−55 % en A) y da 30 caudales negativos en A.
  - Funciona bien en recesión y estiaje; subestima mucho la crecida de El Niño costero de 2017 (abr-2017: −310 m³/s).
- **(ii) Q̄ ← P_I:** mismo tipo de modelo, pero la habilidad es solo +8 % (A) y +12 % (B). IMERG casi no mejora la climatología del caudal.
- **Estabilidad:** el tipo de modelo elegido es el mismo en ambas particiones; el rezago óptimo varía entre 2 y 3 meses.

**Conclusión preliminar** (la interpretación final es del grupo): existe una relación aprovechable, no lineal (logarítmica), estacional y con memoria de 2–3 meses entre la lluvia local y el caudal. No reproduce las crecidas extremas de El Niño ni conserva masa.

**Fuentes compartidas:** IMERG Final incorpora pluviómetros (GPCC) y P_L incluye información satelital (CHIRPS), así que su concordancia no es del todo independiente.

### Punto 3 – Tendencias (script 32; `resultados/tendencias/`)

**Diseño:**
- Registro completo de cada variable más el periodo común 1998–2024.
- Anomalías con climatología fija 1998–2020 (todas las fuentes comparten la referencia; n de años por mes en `climatologia_referencia_3_2.csv`).
- α = 0.05; FDR de Benjamini–Hochberg sobre 60 pruebas (12 meses × 5 variables, registro completo), por método.
- Verificado numéricamente: restar μ_j no cambia la pendiente OLS ni la de Sen (|dif| < 1e−13); z = a/s_j solo reescala; τ de MK idéntico en X y z (`verificacion_transformaciones_3_3.csv`).
- Implementaciones: statsmodels 0.15 (OLS, HAC, LOWESS), pymannkendall 1.4.3 (Hamed–Rao, estacional, Hirsch–Slack); S de MK y Sen con fechas reales, de cálculo propio y contrastados con pymannkendall.

**Resultados (anomalías, todos los meses, Sen [IC 95 %] por década, p MK Hamed–Rao):**
- **Temperatura:** sube en ambas fuentes y es robusta frente a método, periodo, retiro de años extremos y ENSO.
  - PISCOt 1981–2020: +0.18 [0.11, 0.25] °C (p < 0.001).
  - ERA5-Land 1981–2025: +0.23 [0.21, 0.26] °C (p < 0.001).
  - Por mes: ERA5-Land significativa en los 12 meses tras FDR. El mayor cambio estandarizado es en ago–nov (z ≈ +0.6 a +1.0 por década) y el menor en ene–mar.
  - Casi lineal: LOWESS ≈ OLS; con meseta 2005–2010.
- **P_L:** sin tendencia (+1.2 [−0.3, 3.2] mm/mes en 1981–2025, p 0.12; −0.8 [−3.3, 1.5] en 1998–2024).
  - Marzo da +19 mm/mes/década con OLS (q FDR 0.02), pero no con MK tras FDR (q 0.11). No es robusto.
  - PISCOp v3 coincide (+1.5, p 0.03 en 1981–2025; negativa no significativa en 1998–2024).
- **IMERG:** baja −3.0 [−5.6, −0.4] mm/mes/década (p 0.025; −3.6 en 1998–2024, p 0.003).
  - Desaparece al empezar en 2003 (−1.0, p 0.61). Depende de los años iniciales (1998–2002, húmedos).
  - CHIRPS también es negativa (−3.0, p 0.04). P_L y PISCOp no lo son en el mismo periodo.
  - **Discrepancia entre fuentes:** hipótesis de observación/producto (cambio TRMM→GPM, densidad de estaciones) a contrastar.
- **Caudal:** sin tendencia (Sen −0.6 [−2.7, 1.4] m³/s, p 0.53).
  - El OLS (−8.7) está dominado por la crecida de 1998 (pendiente OLS ≫ Sen).
  - La regresión segmentada prefiere un cambio en 2003 (descenso fuerte 1998–2003 y luego plano), artefacto de empezar con El Niño 1997–98. Pettitt no significativo (p 0.25).
  - Febrero en el periodo común: MK −10 m³/s/década (q 0.016), sensible al inicio de 1998.
- **Saltos:** Pettitt significativo solo en temperatura (cambio en 2002; ERA5-Land p 0.001, PISCOt p 0.046). Sugiere un calentamiento con aceleración más que un salto; la segmentada prefiere la lineal.
- **Ventanas inicio/fin:** el signo de las tendencias de lluvia y caudal cambia con el año inicial (El Niño 1983, 1998, 2017 en los extremos); la temperatura es positiva en casi todas.

**3.6 Posibles causas (borrador de hipótesis para el grupo):**
- **H-T (calentamiento regional):**
  - A favor: dos fuentes independientes en método; persiste sin años El Niño y tras regresar sobre Niño 1+2 (residuo +0.18 y +0.23 °C/década).
  - Faltante: estaciones SENAMHI locales.
  - Alternativa observacional: cambios en las observaciones asimiladas por ERA5 y en la interpolación de PISCOt.
- **H-ENSO (variabilidad, no tendencia):** las anomalías de P y Q se correlacionan con Niño 1+2 (r 0.27–0.28). Los Niños de 1983, 1998, 2017 y 2023 en los extremos del registro cambian el signo de la pendiente: es la explicación más simple de las "tendencias" de lluvia y caudal.
- **H-IMERG (artefacto de producto):** la tendencia negativa de IMERG (y CHIRPS) sin equivalente en P_L/PISCOp ni en el caudal apunta a un cambio de producto o de entradas.
  - Faltante: comparar con pluviómetros independientes y con la serie IMERG por era TRMM/GPM.
- **H-ET/usos del agua:**
  - Por año hidrológico, el residuo del caudal dado P_L baja −27 mm/año/década (p 0.03) y R/P baja −0.023/década (p 0.32), mientras T sube +0.28 °C/década.
  - Compatible con más evapotranspiración o con extracciones aguas arriba.
  - Sensible a la crecida de 1998/99–2000 al inicio del registro.
  - **Falta bibliografía/metadatos:** captaciones, minería en la cabecera (p. ej. Lagunas Norte, desde ~2005 — dato de memoria, NO verificado), cambios de cobertura (MapBiomas Perú tiene series anuales; no descargado), curva de calibración de El Tambo.
- **Nieve/deshielo:** cota máxima 4 317 m (`elev_max` de CAMELS-PE) → sin glaciares relevantes; no aplica (verificar con inventario de glaciares si se cita).
- **Referencias metodológicas** (citadas de memoria, verificar antes de usar):
  - Hirsch et al. 1982; Hirsch & Slack 1984; Hamed & Rao 1998; Sen 1968; Gilbert 1987; Cleveland 1979; Benjamini & Hochberg 1995; Pettitt 1979; Hussain & Mahmud 2019 (pymannkendall).
  - Para tendencias en caudales con persistencia: Khaliq et al. 2009; Yue et al. 2002.

### Punto 4 – Espectros de Fourier (script 33; `resultados/espectros/`)

**Cálculo (4.1):**
- Secuencia mensual cronológica, Δt = 1 mes. Ninguna serie tiene vacíos internos (el script lo comprueba); Q es continuo por el relleno.
- N: P_L 540, IMERG 332, Q 336, T PISCOt 480, T ERA5-Land 540; periodo común 324 (276 para T PISCOt).
- IMERG con N = 332 no tiene la frecuencia anual en un bin exacto (332 no es múltiplo de 12) → algo de fuga. En el periodo común (N = 324) sí lo está.
- Densidad unilateral S_k = 2|X_k|²/N, con Σ S·Δf = varianza (cociente 1.000000 en todos los casos). Igual a `scipy.signal.periodogram`.
- Welch: Hann de 120 meses, 50 % de traslape; ν ≈ 15 (N 540) a ≈ 6 (N 276). Pierde periodos > 10 años. La ventana de Hann atenúa los bordes, así que en Q la crecida de 1998 (al inicio) pesa menos en Welch que en el periodograma.
- **Vacíos de Q:** Lomb–Scargle con solo los 244 meses observados completos del periodo común. La banda de 2–5 meses sube de 30 % a 46 % y bajan las largas. El control (Lomb–Scargle en malla regular = FFT) muestra que la diferencia es del muestreo, no del estimador.

**Resultados (4.2):**
- **Series centradas:** domina la banda anual en la lluvia (P_L 52 %, IMERG 45 % de la varianza), con armónico semianual (11–12 %) por la estación húmeda concentrada.
  - En el caudal la anual pesa solo 28 %, porque las crecidas de 1998 y 2017 reparten varianza en todas las frecuencias. Recortando P1–P99, la fracción anual de las anomalías sube de 10 % a 16 %. Con ln Q no cambia (10 %), pero sube la de > 7 años (8 % → 27 %): en escala logarítmica pesan más los años secos persistentes.
  - En la temperatura el ciclo anual es pequeño (PISCOt 3 %, ERA5-Land 7 %) y dominan 2–7 años (36–47 %) y > 7 años (tendencia + baja frecuencia).
- **Anomalías:**
  - La lluvia es casi ruido blanco (ρ1: P_L 0.32, IMERG 0.06; 2–5 meses = 42–58 % de la varianza).
  - El caudal es más rojo (ρ1 0.51).
  - La temperatura es muy roja (ρ1 0.78–0.83; 40–51 % en 2–7 años).
- **Almacenamiento (P_L vs Q, anomalías 1998–2024):** 2–5 meses: Q 30 % vs P 47 %; 13 m–7 años: Q 29 % vs P 19 %. El caudal filtra las variaciones rápidas y gana persistencia. Sin inferir desfases del espectro: estos vienen del punto 2.
- **Picos interanuales de anomalías (pico [rango ±½ bin], ciclos en el registro):**
  - P_L 60 m [57–64], 9 ciclos.
  - IMERG 111 m (3 ciclos) y 41.5 m.
  - Q 112 m (3 ciclos) y 26 m.
  - T PISCOt 44 m [42–46], 11 ciclos.
  - T ERA5-Land 45 m, 12 ciclos.
  - Varios superan el 95 % puntual frente a AR(1); **ninguno (0 de 144 picos tabulados) supera el umbral global**. No son periodicidades demostradas: hay que hablar de bandas (2–7 años) con pocos ciclos.
- **Estabilidad:** las fracciones por banda son estables frente a la tendencia y al recorte de extremos. El periodo del pico interanual no lo es: P_L pasa de 5.0 a 2.4 años con Hann o en la primera mitad, y T PISCOt de 3.6 a 5.7 años con Hann.
- **Fuentes:** P_L e IMERG tienen la misma forma en las series centradas; en anomalías IMERG es más "blanca" (más varianza a 2–5 meses), coherente con su peor captura de la variabilidad interanual (punto 2). PISCOt y ERA5-Land coinciden.

**Para la interpretación (grupo):** la banda de 2–7 años es compatible con ENSO, pero un pico interanual no basta para atribuirlo (guía 4.2). El punto 3 ya mostró correlación de las anomalías con Niño 1+2 (r ≈ 0.27 P_L y Q, 0.56 T PISCOt), y el punto 5 debe probarlo con mapas.

- **Referencias metodológicas** (de memoria, verificar):
  - Welch 1967; Percival & Walden 1993 (ν de Welch);
  - Gilman, Fuglister & Mitchell 1963 y Torrence & Compo 1998 (fondo AR(1));
  - Lomb 1976 / Scargle 1982;
  - Wilks, *Statistical Methods in the Atmospheric Sciences* (periodograma y χ²).

---

## 6. Repositorio GitHub y visor en línea (2026-10-09/10)

- **Repositorio:** https://github.com/jugalvisp-pixel/Hidrologia_Rio_Chicama_El_Tambo (público, rama `main`). El repositorio Git está en la **raíz del proyecto** (`Hidrologia/.git`). Se maneja con **GitHub Desktop** (Git no está en el PATH; el ejecutable de Git de Desktop está en `%LOCALAPPDATA%\GitHubDesktop\app-*\resources\app\git\cmd\git.exe`).
- **Visor en línea (GitHub Pages, `main` / root, con `.nojekyll`):** https://jugalvisp-pixel.github.io/Hidrologia_Rio_Chicama_El_Tambo/visor_resultados.html. Verificado (HTTP 200) el 2026-10-10. Se actualiza solo 1–2 min después de cada *Push*. No hay `index.html`, así que la URL debe terminar en `/visor_resultados.html`.
- **Qué se sube:** controlado por una lista blanca en `.gitignore` (raíz) más `Hidrologiat1/.gitignore`. El de `Hidrologiat1` es el del venv (`*`), con una excepción para `Codigos/`. Al publicar eran 448 archivos (~71 MB, ninguno > 10 MB):
  - en la raíz: `README.md`, `ESTADO_PROYECTO.md`, `requirements.txt`, `visor_resultados.html`, `.nojekyll`;
  - `Hidrologiat1/Codigos/`;
  - `Tarea 1/El_Tambo/datos/` (recortes) y `Tarea 1/El_Tambo/resultados/`.
- **Excluido a propósito:** venv, DEM, geología, NetCDF nacionales de PISCOp (salvo `PISCOp_v3_d_El_Tambo_recorte.nc`), CAMELS-PE completo y ZIPs, guía y rúbrica del curso, `Bibliografia/`, `Datos IMERG/`, `__pycache__`, credenciales (`.cdsapirc`, `_netrc`). **`DescargaIMERG.py` (raíz) tiene la contraseña de Earthdata en texto plano: nunca subirlo**. Cualquier carpeta o archivo nuevo en la raíz queda excluido hasta agregarlo con `!/<nombre>` en `.gitignore`.
- **Flujo para actualizar:** correr los scripts que cambien → `22_visor_html.py` al final → GitHub Desktop: revisar *Changes* (que no aparezca nada de lo excluido) → *Commit to main* → *Push origin*.
- `README.md` (portada del repositorio): cómo abrir el visor, estructura, reproducción, datos no incluidos y fuentes. Falta completar los integrantes.
