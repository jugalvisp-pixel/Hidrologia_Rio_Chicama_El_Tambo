---
editor_options: 
  markdown: 
    wrap: 72
---

# CAMELS-PE: Catchment Attributes and Meteorology for Large-sample Studies - Peru

## Overview

CAMELS-PE is a large-sample hydrological dataset for Peru that provides
consistent hydroclimatic time series, catchment attributes, and
geospatial information for a set of river basins. The dataset follows
the structure of the CAMELS family of datasets and is designed to
support hydrological modeling, regionalization, and climate–streamflow
analysis.

The dataset includes daily meteorological and hydrological variables,
static catchment attributes, and geospatial representations of
catchments and gauging stations.

------------------------------------------------------------------------

## Dataset Structure

```         
CAMELS-PE/
│
├── 01_metadata/
│   ├── stations.csv
│   └── data_dictionary.csv
│
├── 02_attributes/
│   ├── topographic_attributes.csv
│   ├── climatic_indices.csv
│   ├── hydrological_signatures.csv
│   ├── landcover_attributes.csv
│   ├── geologic_attributes.csv
│   ├── soil_attributes.csv
│   └── human_intervention_attributes.csv
│
├── 03_timeseries/
│   ├── timeseries.csv
│   └── by_catchment/
│       ├── PE_XXXX.csv
│       ├── PE_XXXX.csv
│       └── ...
│
└── 04_geospatial/
    ├── camels_pe_gauges.gpkg
    ├── camels_pe_catchments.gpkg
    └── by_catchment/
        ├── PE_XXXX/
        │   ├── PE_XXXX_outlet.gpkg
        │   └── PE_XXXX_catchment.gpkg
        ├── PE_XXXX/
        │   ├── PE_XXXX_outlet.gpkg
        │   └── PE_XXXX_catchment.gpkg
        └── ...
```

------------------------------------------------------------------------

## Metadata

### stations.csv

Contains information for each gauging station.

### data_dictionary.csv

Description, units, and sources of all data.

------------------------------------------------------------------------

## Catchment Attributes (02_attributes)

Each file contains attributes indexed by `gauge_id`.

### Topographic attributes

Derived from FABDEM v1.2: - Area, elevation statistics, slope

### Climatic indices

Derived from PISCO datasets: - Mean precipitation and PET - Aridity
index - Precipitation seasonality and extremes

### Hydrological signatures

Derived from observed and simulated streamflow: - Runoff ratio, baseflow
index - Flow duration curve metrics - High/low flow statistics

### Land cover

Derived from MapBiomas: - Fractional land cover classes - Dominant land
cover type

### Geology

Derived from GLiM and GLHYMPS: - Lithological classes - Porosity and
permeability

### Soil

Derived from DSOLMap: - Soil type fractions - Dominant soil type

### Human intervention

Derived from ANA: - Water availability accreditations (surface and
groundwater) - Reservoir counts, storage, and use

------------------------------------------------------------------------

## Time Series (03_timeseries)

### timeseries.csv

Daily time series in long format: - `date` - `gauge_id` - `prec`,
`prec_var` - `flow_obs`, `flow_sim` - `pet` - `tmin`, `tmean`, `tmax` -
`srad`, `vprp`

### by_catchment/

Individual CSV files per catchment containing the same variables in wide
format.

### Temporal coverage

1981–2025 (daily)

------------------------------------------------------------------------

## Geospatial Data (04_geospatial)

### camels_pe_gauges.gpkg

Point locations of gauging stations.

### camels_pe_catchments.gpkg

Polygon boundaries of catchments.

### by_catchment/

Optional per-catchment geospatial files: - Outlet point - Catchment
polygon

All geospatial data are provided in WGS84 (EPSG:4326).

------------------------------------------------------------------------

## Data Sources

-   **PISCOp v2.1**: Precipitation
-   **PISCOt v2.1**: Air temperature
-   **PISCOeo_pm v1.0**: Potential evapotranspiration
-   **ERA5-Land**: Solar radiation and vapor pressure
-   **PISCO_ARNOVIC v1.1**: Simulated streamflow
-   **SENAMHI**: Observed streamflow
-   **FABDEM v1.2**: Elevation
-   **MapBiomas (2025)**: Land cover
-   **GLiM / GLHYMPS v2.0**: Geological properties
-   **DSOLMap**: Soil classes
-   **ANA (Peru)**: Human water use and infrastructure

------------------------------------------------------------------------

## Units

-   Precipitation, PET, streamflow: mm d⁻¹
-   Temperature: °C
-   Radiation: MJ m⁻² d⁻¹
-   Vapor pressure: hPa
-   Area: km²

------------------------------------------------------------------------

## Usage

The dataset is designed for: - Hydrological modeling - Regionalization
studies - Climate–streamflow analysis - Machine learning applications

------------------------------------------------------------------------

## Version

CAMELS-PE v1.0\
Period: 1981–2025

------------------------------------------------------------------------

## Citation

If you use this dataset, please cite:

> CAMELS-PE dataset (v1.0), [Year], Zenodo/PANGAEA repository.

(Replace with final DOI after publication)

------------------------------------------------------------------------

## License

This dataset is distributed under the terms of the CC-BY 4.0 license.

------------------------------------------------------------------------

## Contact

For questions or collaborations:

Harold Llauca
[hllauca\@senamhi.gob.pe](mailto:hllauca@senamhi.gob.pe){.email}\
Servicio Nacional de Meteorología e Hidrología del Perú (SENAMHI)
