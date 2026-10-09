"""
Descarga de capas de contexto físico para la cuenca El Tambo (río Chicama):
  - Cobertura del suelo: ESA WorldCover 10 m v200 (2021), lectura por ventana del COG público (AWS).
  - Geología 1:100 000 integrada y fallas: INGEMMET GEOCATMIN (servicio ArcGIS REST).
  - Proyectos mineros (cartera, inversión y exploración): INGEMMET GEOCATMIN.
  - Suelos: SoilGrids 2.0 WRB "MostProbable" (ISRIC, WCS) – descargado con WCS GetCoverage.
Salidas en Tarea 1/El_Tambo/datos/{cobertura,geologia,suelos}/.
"""
from pathlib import Path
import json
import urllib.parse
import urllib.request
import ssl

import geopandas as gpd
import rasterio
from rasterio.windows import from_bounds

ROOT = Path(__file__).resolve().parents[2] / "Tarea 1" / "El_Tambo"
DATOS = ROOT / "datos"
BBOX = (-78.80, -8.05, -78.20, -7.30)   # lon_min, lat_min, lon_max, lat_max
GEOCATMIN = "https://geocatmin.ingemmet.gob.pe/arcgis/rest/services"
# El servidor de GEOCATMIN no envía la cadena de certificados intermedios completa.
CTX = ssl.create_default_context()
CTX.check_hostname = False
CTX.verify_mode = ssl.CERT_NONE


def arcgis_query(layer_url, out):
    if out.exists():
        print("ya existe:", out.name)
        return gpd.read_file(out)
    env = {"xmin": BBOX[0], "ymin": BBOX[1], "xmax": BBOX[2], "ymax": BBOX[3],
           "spatialReference": {"wkid": 4326}}
    q = {"geometry": json.dumps(env), "geometryType": "esriGeometryEnvelope", "inSR": 4326,
         "spatialRel": "esriSpatialRelIntersects", "outFields": "*", "returnGeometry": "true",
         "outSR": 4326, "f": "geojson"}
    url = f"{layer_url}/query?" + urllib.parse.urlencode(q)
    with urllib.request.urlopen(url, context=CTX, timeout=300) as r:
        data = r.read()
    out.write_bytes(data)
    gdf = gpd.read_file(out)
    print(f"{out.name}: {len(gdf)} entidades")
    return gdf


# ---- Cobertura: ESA WorldCover 2021 (tile S09W081 cubre lat -9..-6, lon -81..-78)
wc_out = DATOS / "cobertura" / "ESA_WorldCover_10m_2021_v200_El_Tambo.tif"
wc_out.parent.mkdir(parents=True, exist_ok=True)
if not wc_out.exists():
    url = ("/vsicurl/https://esa-worldcover.s3.eu-central-1.amazonaws.com/v200/2021/map/"
           "ESA_WorldCover_10m_2021_v200_S09W081_Map.tif")
    with rasterio.open(url) as src:
        win = from_bounds(*BBOX, transform=src.transform)
        arr = src.read(1, window=win)
        meta = src.meta.copy()
        meta.update(height=arr.shape[0], width=arr.shape[1],
                    transform=src.window_transform(win), compress="deflate", driver="GTiff")
    with rasterio.open(wc_out, "w", **meta) as dst:
        dst.write(arr, 1)
print("WorldCover:", wc_out.name)

# ---- Geología, fallas y minería (INGEMMET)
geo_dir = DATOS / "geologia"
geo_dir.mkdir(parents=True, exist_ok=True)
geo = arcgis_query(f"{GEOCATMIN}/SERV_GEOLOGIA_100K_INTEGRADA/MapServer/6",
                   geo_dir / "ingemmet_geologia_100k.geojson")
arcgis_query(f"{GEOCATMIN}/SERV_GEOLOGIA_100K_INTEGRADA/MapServer/3",
             geo_dir / "ingemmet_fallas_100k.geojson")
for i, name in [(0, "proyectos_mineros"), (1, "proyectos_inversion"), (2, "proyectos_exploracion")]:
    arcgis_query(f"{GEOCATMIN}/SERV_CARTERA_PROYECTOS_MINEROS/MapServer/{i}",
                 geo_dir / f"ingemmet_{name}.geojson")

# ---- Suelos: SoilGrids WRB (si no existe)
soil = DATOS / "suelos" / "soilgrids_wrb_mostprobable.tif"
soil.parent.mkdir(parents=True, exist_ok=True)
if not soil.exists():
    url = ("https://maps.isric.org/mapserv?map=/map/wrb.map&SERVICE=WCS&VERSION=2.0.1"
           "&REQUEST=GetCoverage&COVERAGEID=MostProbable&FORMAT=image/tiff"
           f"&SUBSET=long({BBOX[0]},{BBOX[2]})&SUBSET=lat({BBOX[1]},{BBOX[3]})"
           "&SUBSETTINGCRS=http://www.opengis.net/def/crs/EPSG/0/4326"
           "&OUTPUTCRS=http://www.opengis.net/def/crs/EPSG/0/4326")
    urllib.request.urlretrieve(url, soil)
print("SoilGrids:", soil.name)

# Unidades geológicas presentes (para definir la agrupación litológica)
cols = [c for c in ["NAME", "UNIDAD", "DESCRIP"] if c in geo.columns]
tab = (geo.to_crs(32717).assign(km2=lambda d: d.area / 1e6)
       .groupby(cols)["km2"].sum().sort_values(ascending=False).round(1))
tab.to_csv(geo_dir / "unidades_geologicas_bbox.csv")
import pandas as pd
pd.set_option("display.width", 250, "display.max_colwidth", 110)
print(tab.to_string())
