# Datos de la Tarea 1 – cuenca del río Chicama hasta El Tambo

Qué se descargó, qué tratamiento recibió cada conjunto y con qué versión final se trabaja.

- **Cuenca:** río Chicama hasta la estación hidrométrica El Tambo. Es la estación CAMELS-PE `PE_47E0A2A8` (SENAMHI 211003). Está en −7.574°, −78.707°, a 684 m s.n.m., en la vertiente del Pacífico del norte del Perú, y drena 2 170 km².
- **Periodo de análisis (periodo común):** 1998-01 a 2024-12, 324 meses (27 años calendario completos). Son los meses que tienen a la vez la precipitación de CAMELS-PE y la de IMERG, recortados a años completos. Los meses 2025-01 a 2025-08 de IMERG quedan fuera. Se calcula en un solo lugar, `Hidrologiat1/Codigos/periodo_comun.py`, y todas las figuras, mapas y tablas lo usan.
- **Rutas:** todas son relativas a la raíz del repositorio.

---

## 1. Resumen: con qué se trabaja

| Variable | Fuente final | Resolución | Serie de trabajo | Archivo de trabajo |
|---|---|---|---|---|
| **P_L** – precipitación de referencia | CAMELS-PE v1.0.1 (`prec`, PISCOp v2.1 promediado sobre la cuenca) | diaria → mensual | 324/324 meses, sin faltantes | `resultados/series_mensuales/series_mensuales_camels_pe.csv` (`P_L_mm_mes`) |
| **P_I** – precipitación satelital | IMERG V07 Final mensual | 0.1° → promedio de cuenca | 324/324 meses | `resultados/series_mensuales/imerg_cuenca_mensual.csv` |
| **Q̄** – caudal medio mensual | CAMELS-PE (`flow_obs`, SENAMHI) | diario → mensual | 324/324 meses **con relleno** (297 observados) | `series_mensuales_camels_pe.csv` (`Q_m3_s`) |
| **R** – escorrentía en lámina | ídem | diario → mensual | ídem | `series_mensuales_camels_pe.csv` (`R_mm_mes`, bandera `R_flag`) |
| **T** – temperatura media | CAMELS-PE (`tmean`, PISCOt v1.2) hasta 2020-12; desde 2021, estimación con ERA5-Land | diaria → mensual | 276 meses PISCOt + 56 estimados, **en columnas separadas** | `resultados/series_mensuales/temperatura_piscot_extendida_era5land.csv` |
| PET (solo 1.5.c) | CAMELS-PE (`pet`, PISCOeo_pm, Penman–Monteith) | diaria → mensual | 1998-01 a 2016-12 | calculada en el script 23 |
| Campos de apoyo | PISCOp v3.0, CHIRPS v2.0, ERA5 | mallas | solo mapas y contrastes | ver sección 2 |

`Tarea 1/El_Tambo/` se abrevia en lo que sigue: `datos/` = `Tarea 1/El_Tambo/datos/` y `resultados/` = `Tarea 1/El_Tambo/resultados/`.

---

## 2. Conjuntos descargados

### 2.1 CAMELS-PE v1.0.1 — fuente principal de P_L, Q y T

| | |
|---|---|
| Referencia | CAMELS-PE: Catchment Attributes and Meteorology for Large-sample Studies in Peru, v1.0.1. Zenodo, doi:10.5281/zenodo.21195425 |
| Descarga | ZIP completo desde Zenodo (script 09), con MD5 verificado contra la API de Zenodo (`13f127d381338eee0e35359c08dba199`) |
| Extraído | Serie diaria `PE_47E0A2A8.csv` (1981-01-01 a 2025-12-31), polígono y punto de salida de la cuenca, atributos, metadatos de la estación, diccionario de datos y `PROCEDENCIA.json` |
| Variables usadas | `prec` (mm/d; PISCOp v2.1), `flow_obs` (mm/d; SENAMHI), `tmean`, `tmin`, `tmax` (°C; PISCOt v1.2), `pet` (mm/d; PISCOeo_pm v1.0), `flow_sim` (mm/d; PISCO-ARNOVIC v1.1, solo para el relleno) |
| Dónde está | `datos/camels_pe/` (en el repositorio). El ZIP completo (116 MB) no se sube; el script 09 lo descarga |

**Tratamiento:**

1. **Controles** (scripts 09 y 19):
   - índice diario continuo, sin fechas duplicadas;
   - sin negativos en `prec`, `flow_obs` ni `pet`;
   - `NA` como único código de faltante;
   - `tmean` = (`tmin` + `tmax`)/2;
   - revisión manual de 2010-03 y 1998-03, recalculadas desde los datos diarios.
2. **Agregación mensual** (script 09):
   - **P_L:** P_L,m = Σ P_d, solo con meses completos. `prec` no tiene ningún día faltante en 1981–2025.
   - **Escorrentía:** R_m = Σ q_d. El caudal ya viene en mm/d, así que no se repite la conversión por área.
   - **Caudal en m³/s:** Q̄ = q̄ [mm/d] · A / 86.4, con A = 2 171.585 km² (área de CAMELS-PE).
   - **Temperatura:** T_m = media de los días válidos, si faltan ≤ 10 % de los días.
3. **Relleno del caudal.** Decisión del grupo, 2026-10-10: el caudal tiene 3.1 % de días faltantes, concentrados en septiembre y en 2019–2021. Antes, los meses con > 10 % de días faltantes se excluían. Ahora se rellenan con `relleno_caudal.py`, aplicado en el script 09 y validado en el script 28:

   | Paso | Cuándo | Método | Error medido en la validación |
   |---|---|---|---|
   | 1 | Huecos de ≤ 15 días con datos a ambos lados | Interpolación lineal de ln(q), que sigue la forma de la recesión | Mediana del total mensual 0.1–1.5 % (≤ 7 días) y 2–4 % (10–15 días) |
   | 2 | Meses aún incompletos con datos en ≥ 50 % de los días | Media de los días disponibles × días del mes | ~10 % con 15 días faltantes |
   | 3 | Meses con < 50 % de los días | ln R_m = a_mes + b·ln(`flow_sim`_m) + c·ln R_m−1 | Validación cruzada dejando un año fuera: NSE 0.72, error mediano 23 % |

   **Resultado en el periodo común:**

   | Condición del mes | Meses |
   |---|---|
   | Observados completos | 244 |
   | Con días interpolados | 76 |
   | Escalados | 3 (1998-09, 2019-01, 2019-02) |
   | Reconstruidos con el modelo | 1 (1998-08, sin ningún día medido) |

   La tabla conserva lo observado en `R_obs_mm_mes`, `R_flag_obs` y `Q_obs_*`. La bandera `R_flag` indica el origen de cada mes.
4. **`prec` no es una observación independiente.** Es PISCOp v2.1: una interpolación de pluviómetros de SENAMHI con CHIRPS, que los autores de CAMELS-PE promediaron sobre la cuenca. Además, **no se pudo reproducir** con ninguna versión pública de PISCOp (ver 2.3), y eso se declara como limitación.

**Con qué se trabaja:** `resultados/series_mensuales/series_mensuales_camels_pe.csv`, con estas columnas:
- `P_L_mm_mes`
- `R_mm_mes`, `Q_m3_s` y `R_flag` (serie rellenada)
- `T_C`

Además, `caudal_diario_interpolado.csv` tiene la serie diaria con los huecos cortos rellenados, para el punto 4.

### 2.2 IMERG V07 Final mensual — P_I

| | |
|---|---|
| Referencia | GPM IMERG Final Precipitation L3 1 month 0.1° V07 (GPM_3IMERGM.07), NASA GES DISC, doi:10.5067/GPM/IMERG/3B-MONTH/07 |
| Descarga | Manual desde GES DISC (Subset / OPeNDAP, grupo `Grid`), recortada a la cuenca: 332 archivos `.nc4`, 1998-01 a 2025-08. El script 10 permite también la descarga automática (requiere cuenta de NASA Earthdata) |
| Variable | `precipitation`, tasa media mensual en **mm/h** |
| Dónde está | `datos/imerg/manual/` y el recorte unido `datos/imerg/GPM_3IMERGM.07_El_Tambo_recorte.nc` |

**Tratamiento** (scripts 10 y 12):
1. **Conversión a acumulado mensual:** mm/mes = tasa · 24 · días del mes.
2. **Promedio de cuenca ponderado por el área de intersección** entre cada celda de 0.1° y el polígono, calculada en UTM 17S. Intervienen 31 celdas: 10 completas, con el 43.7 % del peso en celdas de borde (≈ 17.8 celdas equivalentes). No se usa el píxel de la estación. Los pesos no consideran la topografía; el desnivel interno de una celda llega a 2 200 m.
3. **Controles:** sin meses faltantes ni duplicados, y revisión manual de 2010-03.

**Con qué se trabaja:** `resultados/series_mensuales/imerg_cuenca_mensual.csv`.

**Resultado frente a P_L:** P_I acumula el 65 % de P_L (539 frente a 833 mm/año), con r = 0.90.

### 2.3 PISCOp (SENAMHI) — verificación y apoyo espacial

| Versión | Referencia | Uso |
|---|---|---|
| PISCOp v2.1 (IRI Data Library) | — | **Descartada** (2026-10-07). Daba 20–25 % menos que `prec` de CAMELS-PE. Datos borrados; queda el script 08, fuera del flujo |
| PISCOp v2.1 update, mensual | figshare, doi:10.6084/m9.figshare.21127423.v2 | Solo verificación (script 13) |
| PISCOp v3.0, mensual y diario | figshare, doi:10.6084/m9.figshare.32411886.v1 | Verificación y campo de apoyo (scripts 13 y 14) |

**Tratamiento:**
1. **Descarga con MD5 verificado** (script 13).
2. **Comparación del promedio de cuenca con `prec` de CAMELS-PE.** Se probaron 2 polígonos × 2 métodos de pesos × 3 periodos. Ninguna versión reproduce CAMELS-PE: v3.0 da 718 frente a 841 mm/año (−14.5 %), con r mensual 0.98.
3. **Recorte diario de v3.0 a la cuenca** (script 14), con 31 celdas, 0 faltantes y 0 negativos, y su serie de cuenca ponderada por área.

**Con qué se trabaja:**
- `datos/piscop_senamhi/v3p0/PISCOp_v3_d_El_Tambo_recorte.nc` (2.9 MB, en el repositorio)
- `resultados/series_mensuales/piscop_v3_cuenca_{diaria,mensual}.csv`

Se usa solo como contraste de magnitud (aridez en 1.5.c). Los NetCDF nacionales (1.5 GB) no se suben.

### 2.4 CHIRPS v2.0 — mapas continuos de precipitación

| | |
|---|---|
| Referencia | Climate Hazards group InfraRed Precipitation with Station data v2.0, mensual, 0.05°. Funk et al. (2015), Scientific Data, doi:10.1038/sdata.2015.66; dataset doi:10.15780/G2RP4Q |
| Descarga | Recorte 79°O–78°O × 8.2°S–7.2°S del IRI Data Library, periodo común (script 26) |
| Dónde está | `datos/chirps/CHIRPS_v2p0_mensual_El_Tambo.nc` (0.5 MB, en el repositorio) |

**Tratamiento** (scripts 26, 07b y 27):
1. **Corrección de la georreferencia.** En el IRI la longitud marca el borde oeste de la celda y no su centro, así que se suma 0.025°. Se verificó contra el GeoTIFF oficial de UCSB (2017-03): la diferencia es 0.000, y sin corregir llegaba a ~94 mm.
2. **Serie de cuenca ponderada por área**, con 98 celdas.
3. **Climatología del periodo común por celda**, interpolada de forma bilineal a la malla del DEM. La interpolación solo sirve para representar el campo de forma continua; el detalle real sigue siendo de 0.05°.

**Con qué se trabaja:** mapas de precipitación media anual (panel de contexto físico, (d)) y de precipitación media de cada mes. CHIRPS da 636 mm/año frente a 841 de P_L, así que sirve para el **patrón espacial, no para la magnitud**.

### 2.5 ERA5-Land — temperatura desde 2021

| | |
|---|---|
| Referencia | ERA5-Land monthly averaged data from 1950 to present, C3S/ECMWF, doi:10.24381/cds.68d2bb30 |
| Descarga | CDS con `cdsapi` (script 17): T2m mensual 1981-01 a 2025-12, 0.1°, recorte de la cuenca, más el geopotencial (relieve del modelo). Requiere cuenta personal en el CDS (`~/.cdsapirc`, fuera del repositorio) |
| Dónde está | `datos/era5land/` |

**Tratamiento** (script 17):
1. **Promedio de cuenca ponderado por área**, con 30 celdas.
2. **Comparación con PISCOt en 1998–2020.** Las dos series varían igual (r 0.89), pero ERA5-Land es 2.84 °C más fría de forma casi constante. La diferencia de relieve del modelo explica solo ≈ 0.55 °C.
3. **Estimación desde 2021:** T_est = ERA5-Land + sesgo medio del mes calendario frente a PISCOt (+2.62 a +2.99 °C). Validación fuera de muestra (ajuste en 1998–2010, evaluación en 2011–2020): MAE 0.19 °C.

**Con qué se trabaja:** `temperatura_piscot_extendida_era5land.csv`. `T_piscot_C` (hasta 2020) y `T_estimada_era5land_C` (2021-01 a 2024-12 en el periodo común) van **en columnas separadas**, sin unir las fuentes. En las figuras, el tramo estimado se dibuja discontinuo y morado. Las tendencias (punto 3) deben tratar ese tramo aparte.

### 2.6 ERA5 — viento y transporte de humedad regional

| | |
|---|---|
| Referencia | ERA5 monthly averaged reanalysis, C3S/ECMWF: niveles de presión (doi:10.24381/cds.6860a573) y nivel único (doi:10.24381/cds.f17050d7) |
| Descarga | CDS (script 24), 1998–2025, malla 0.5°, dominio 20°N–40°S × 110°O–20°O (Atlántico tropical, Amazonía, Andes y Pacífico oriental). Variables: u, v y q en 850, 700, 500 y 200 hPa; flujo de vapor integrado en la vertical (IVT, este y norte); agua precipitable; presión en superficie |
| Dónde está | `datos/era5/`. Los NetCDF mensuales (~230 MB) no se suben; sí la climatología estacional (1.4 MB) |

**Tratamiento** (scripts 24 y 25):
1. **Recorte al periodo común.**
2. **Máscara de niveles bajo el terreno del modelo:** se descartan las celdas con presión en superficie menor que la del nivel.
3. **Climatología estacional a 1°** (dic–abr, may–sep, oct–nov y anual).

**Con qué se trabaja:**
- recuadro de Sudamérica del mapa de la cuenca (IVT de dic–abr);
- mapa regional de IVT y viento en 500 hPa;
- cálculo de la convergencia del IVT sobre la cuenca.

Para la cuenca se descartó ERA5-Land, porque no tiene viento en altura ni IVT.

### 2.6 b ERA5 — temperatura y velocidad vertical (inversión térmica)

| | |
|---|---|
| Referencia | ERA5 monthly averaged reanalysis, C3S/ECMWF (mismos DOI que 2.6) |
| Descarga | CDS (script 29), 1998–2025, malla 0.25°, dominio 2°S–12°S × 88°O–72°O. Temperatura y velocidad vertical ω en 23 niveles (1000–200 hPa), en bloques de 7 años por el límite de tamaño del CDS; presión en superficie y temperatura a 2 m; geopotencial de superficie (relieve del modelo) |
| Dónde está | `datos/era5_inversion/`. Los NetCDF (~77 MB) no se suben; sí la climatología del corte (0.1 MB) |

**Tratamiento** (script 29):
1. **Temperatura potencial**, θ = T·(1000/p)^0.286.
2. **Capas de inversión:** las capas donde T aumenta con la altura.
3. **Estabilidad de la baja troposfera**, LTS = θ(700 hPa) − θ(superficie), frente a la costa (83–80°O × 9–6.5°S).
4. **ω en 500 hPa y convergencia del IVT** sobre las celdas de la cuenca.
5. **Corte vertical promediado en 7–8°S** para dic–abr y may–sep.
6. **Correlación de anomalías mensuales** con P_L en dic–abr.

**Con qué se trabaja:** `fig_inversion_termica_era5.png`, `ciclo_inversion_convergencia.csv` y `correlaciones_inversion_lluvia.csv`, en `resultados/contexto_fisico/`.

### 2.7 Relieve, cobertura, geología y otros

| Conjunto | Referencia / acceso | Tratamiento | Uso | ¿En el repositorio? |
|---|---|---|---|---|
| **Copernicus GLO-30 DSM** (30 m) | ESA / Airbus, tiles COG públicos (script 06) | Mosaico, recorte, reproyección a UTM 17S a 30 m, acondicionamiento hidrológico, direcciones D8 (pysheds) | Delimitación (2 170.0 km²; IoU 0.994 con CAMELS-PE), red de drenaje, morfometría, relieve sombreado | No (175 MB; lo regenera el script 06) |
| **ESA WorldCover 2021** (10 m) | ESA, COG en AWS (script 07a) | Recorte y remuestreo por moda a la malla del DEM | Mapa de cobertura | Sí |
| **Geología INGEMMET 1:100 000** | GEOCATMIN, servicio ArcGIS REST (script 07a) + simbología oficial de la capa (script 07b) | Recorte a la cuenca; 25 unidades agrupadas en 6 grupos hidrogeológicos más lagunas; polígonos disueltos por unidad o grupo (±25 m) para quitar los bordes entre hojas | Mapa hidrogeológico del panel y mapa de unidades con colores oficiales | No los vectores (22 MB; los regenera el 07a); sí la simbología |
| SoilGrids WRB (250 m) | ISRIC, WCS (script 07a) | Recorte | **No se usa** (retirado a pedido del grupo) | Sí |
| Estaciones SENAMHI | Mapa web de estaciones de SENAMHI, copia fechada 2026-10-06 | Filtrado a la zona del mapa | Mapa de la cuenca | Sí |
| Límites de países | Natural Earth 1:50m | — | Mapas de ubicación | Sí |
| **Índice Niño 1+2** (ERSST v5, mensual 1950→) | CPC/NOAA, `https://www.cpc.ncep.noaa.gov/data/indices/ersst5.nino.mth.91-20.ascii` (descarga 2026-10-11; script 32 la repite si falta) | Se usa la columna ANOM de NINO1+2 (anomalías respecto a 1991–2020) | Punto 3.6: ¿la variabilidad ENSO produce tendencias aparentes? | Sí (`datos/indices_clima/`) |

---

## 3. Criterios comunes a todo el análisis

- **Periodo común:** 1998-01 a 2024-12, solo años calendario completos. Cada mes calendario tiene 27 años; la media anual (Σ de las 12 medias mensuales) equivale a la media de los totales anuales.
- **Promedios de cuenca de productos en malla** (IMERG, PISCOp, CHIRPS, ERA5-Land): media ponderada por el área de intersección celda–polígono, calculada en UTM 17S (EPSG:32717) con el polígono propio de GLO-30.
- **Coordenadas:** los mapas van en WGS 84 (EPSG:4326); áreas y morfometría en UTM 17S.
- **Comparaciones entre fuentes:** siempre con los mismos meses válidos (pares comunes).
- **Lo que no se rellena ni se une:**
  - P_L no se rellena con IMERG.
  - PISCOt y ERA5-Land no se unen en una sola serie.
- **Lo que sí se rellena:** solo el caudal, con el método de 2.1 y la bandera `R_flag`. En el punto 2 conviene repetir el análisis solo con los meses `R_flag = completo` como prueba de sensibilidad.
- **Excepción en el punto 3 (tendencias):** cada variable usa además su **registro completo** (P_L 1981–2025, IMERG 1998-01..2025-08, Q 1998–2025, T PISCOt 1981–2020, T ERA5-Land 1981–2025), como pide la guía (3.1). Las anomalías usan una climatología fija 1998–2020, el periodo más largo con las cinco series. Q_obs (solo meses completos) no tiene ningún septiembre completo en 1998–2020, así que sus anomalías de septiembre no se definen.

---

## 4. Lo que no está en el repositorio y cómo obtenerlo

| Archivo | Tamaño | Cómo regenerarlo |
|---|---|---|
| CAMELS-PE completo (ZIP) | 116 MB | Script 09 (Zenodo) |
| DEM GLO-30 | 175 MB | Script 06 |
| Geología INGEMMET (vectores) | 22 MB | Script 07a |
| PISCOp nacionales (v2.1 update y v3.0) | 1.5 GB | Script 13 (figshare) |
| ERA5 mensual regional | 233 MB | Script 24 (cuenta CDS) |

Las credenciales (NASA Earthdata y Copernicus CDS) son personales y nunca se guardan en el repositorio.

Orden de ejecución completo, scripts y decisiones: ver `README.md` y `ESTADO_PROYECTO.md`.
