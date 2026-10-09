import os
import pandas as pd

# 1. Definir rutas y el ID de la estación de Cottbus
PATH_CAMELS = "./data/camels_de"  # Cambia esta ruta según la ubicación de tus carpetas
GAUGE_ID = "DE400300"

# 2. Buscar información en los atributos topográficos
path_topo = os.path.join(PATH_CAMELS, "attributes", "camels_de_topo.csv")

if os.path.exists(path_topo):
    df_topo = pd.read_csv(path_topo)
    meta_cot = df_topo[df_topo["gauge_id"] == GAUGE_ID]

    print("=== METADATOS DE LA CUENCA DE COTTBUS ===")
    print(f"ID Estación: {GAUGE_ID}")
    print(f"Nombre: {meta_cot['gauge_name'].values[0]}")
    print(f"Río: {meta_cot['river'].values[0]}")
    print(f"Área delimitada MERIT Hydro: {meta_cot['area'].values[0]} km²")
    print(f"Área oficial metadata: {meta_cot['area_metadata'].values[0]} km²")
    print("-" * 40)

# 3. Función para cargar las series de tiempo diarias
def cargar_serie_diaria(directorio_var, gauge_id):
    """Busca el archivo CSV o TXT asociado al gauge_id en una carpeta específica."""
    for root, dirs, files in os.walk(directorio_var):
        for file in files:
            if gauge_id in file and file.endswith((".csv", ".txt")):
                filepath = os.path.join(root, file)
                df = pd.read_csv(
                    filepath, parse_dates=["date"], index_col="date"
                )
                return df
    return None

# 4. Extraer series de tiempo diarias de caudal (streamflow) y precipitación
path_streamflow = os.path.join(PATH_CAMELS, "timeseries", "streamflow")
path_precip = os.path.join(PATH_CAMELS, "timeseries", "precipitation")
path_temp = os.path.join(PATH_CAMELS, "timeseries", "temperature")

df_q = cargar_serie_diaria(path_streamflow, GAUGE_ID)
df_p = cargar_serie_diaria(path_precip, GAUGE_ID)
df_t = cargar_serie_diaria(path_temp, GAUGE_ID)

# 5. Unir todas las variables en un único DataFrame diario
if df_q is not None and df_p is not None:
    # Unir por fecha
    df_cuenca = df_q.join([df_p, df_t], how="inner")

    print("\n=== PRIMEROS REGISTROS DIARIOS DE COTTBUS ===")
    print(df_cuenca.head())

    # Guardar en la carpeta de procesados para no repetir la búsqueda
    os.makedirs("./processed", exist_ok=True)
    df_cuenca.to_csv(f"./processed/serie_diaria_{GAUGE_ID}.csv")
    print(
        f"\n¡Datos extraídos y guardados con éxito en './processed/serie_diaria_{GAUGE_ID}.csv'!"
    )

else:
    print(
        f"No se encontraron los archivos para el código {GAUGE_ID}. Verifica las rutas de entrada."
    )