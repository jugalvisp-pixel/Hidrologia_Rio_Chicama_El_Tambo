# Tarea 1 – Hidrología 2026-2S (UNAL Medellín)

Clasificación hidroclimática mensual de la cuenca del río Chicama hasta la estación El Tambo
(CAMELS-PE `PE_47E0A2A8`, SENAMHI 211003, 2170 km²), Perú.

**Integrantes:** _(completar)_

## Ver los resultados

Abre **`visor_resultados.html`** en cualquier navegador. Es un solo archivo con todas las figuras
incrustadas: funciona sin internet y sin instalar nada.

- Desde GitHub: abre el archivo → botón **Download raw file** → doble clic en el archivo descargado.
- Si el repositorio tiene GitHub Pages activo, se puede ver directamente en línea.

## Estructura

| Ruta | Contenido |
|---|---|
| `visor_resultados.html` | Visor con todos los productos y descripciones breves |
| `Hidrologiat1/Codigos/` | Scripts de Python numerados en orden de trabajo |
| `Tarea 1/El_Tambo/datos/` | Datos ya recortados a la cuenca (CAMELS-PE, IMERG, PISCOp v3, ERA5-Land, cobertura, suelos) |
| `Tarea 1/El_Tambo/resultados/` | Figuras (PNG) y tablas (CSV) generadas por los scripts |
| `ESTADO_PROYECTO.md` | Bitácora: fuentes, decisiones metodológicas y notas para el informe |

Periodo común de análisis: **1998-01 a 2025-08** (meses con precipitación de CAMELS-PE e IMERG).

## Reproducir

```bash
python -m venv .venv
.venv\Scripts\activate          # Windows  (Linux/macOS: source .venv/bin/activate)
pip install -r requirements.txt
cd Hidrologiat1/Codigos
```

Orden de ejecución: `06 → 07a → 13 → 14 → 09 → 10 → 07b → 11 → 12 → 17 → 15 → 16 → 18 → 19 → 20 → 21 → 22`.

No se incluyen en el repositorio, por tamaño, el DEM Copernicus GLO-30, la geología INGEMMET,
los NetCDF nacionales de PISCOp ni el CAMELS-PE completo. Los scripts 06, 07a y 13 los descargan,
y CAMELS-PE v1.0.1 se obtiene de su repositorio oficial. Las descargas de IMERG (script 10) y
ERA5-Land (script 17) requieren cuentas personales gratuitas en NASA Earthdata (`~/_netrc`) y en
Copernicus CDS (`~/.cdsapirc`). Las credenciales nunca se guardan en el repositorio.

## Fuentes

CAMELS-PE v1.0.1 (PISCOp v2.1, PISCOt v1.2, caudales SENAMHI) · GPM IMERG V07 Final mensual ·
PISCOp v3.0 (doi:10.6084/m9.figshare.32411886.v1) · ERA5-Land mensual (doi:10.24381/cds.68d2bb30) ·
Copernicus GLO-30.
