"""
Datos de la cuenca El Tambo (PE_47E0A2A8) desde CAMELS-PE v1.0.1 (Zenodo) y series mensuales.

Fuente única de P_L, Q y T: CAMELS-PE v1.0.1, DOI 10.5281/zenodo.21195425
(P_L = `prec`, PISCOp v2.1 promediado por los autores sobre la cuenca; Q = `flow_obs`,
SENAMHI, en mm/d; T = `tmean`, PISCOt v1.2). No se usan las mallas de PISCOp del IRI.

1. Descarga el ZIP desde el registro Zenodo si no existe y verifica su MD5 contra la API.
2. Extrae solo los archivos de la cuenca (serie diaria, polígono, salida, metadatos y
   atributos de PE_47E0A2A8) a Tarea 1/El_Tambo/datos/camels_pe/, con un archivo de procedencia.
3. Agrega a escala mensual con criterio de completitud explícito y genera la tabla de
   disponibilidad año × mes (punto 1.4 de la guía).

Criterios de agregación mensual (aplicados a todo el análisis):
  - P_L,m [mm/mes]: suma diaria solo si el mes está completo (0 días faltantes); si no, NaN.
  - Q [mm/d y m3/s] y R_m [mm/mes]: mes válido si faltan ≤ 10 % de los días; Q es la media de
    los días válidos y R_m = Q_media · n_días (no se suman parciales como acumulados completos).
  - T [°C]: media de los días válidos si faltan ≤ 10 % de los días.
  - Q en m3/s = Q[mm/d] · A / 86.4, con A = área CAMELS-PE (la misma usada por los autores para
    expresar el caudal en lámina).
"""
from pathlib import Path
import hashlib
import json
import urllib.request
import zipfile

import numpy as np
import pandas as pd

GID = "PE_47E0A2A8"
RECORD = "21195425"
DOI = "10.5281/zenodo.21195425"
ZIP_NAME = "CAMELS-PE_v1.0.1.zip"
MAX_FALT = 0.10

ROOT = Path(__file__).resolve().parents[2] / "Tarea 1"
ZIP_DIR = ROOT / "CAMELS_PE"
DEST = ROOT / "El_Tambo" / "datos" / "camels_pe"
OUT = ROOT / "El_Tambo" / "resultados" / "series_mensuales"
for p in (ZIP_DIR, DEST, OUT):
    p.mkdir(parents=True, exist_ok=True)


def md5(path, chunk=2**20):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


# ------------------------------------------------------------------ 1. descarga y verificación
rec = json.load(urllib.request.urlopen(f"https://zenodo.org/api/records/{RECORD}", timeout=60))
finfo = next(f for f in rec["files"] if f["key"] == ZIP_NAME)
md5_ref = finfo["checksum"].split(":")[1]
zpath = ZIP_DIR / ZIP_NAME
if not zpath.exists():
    print(f"Descargando {ZIP_NAME} ({finfo['size']/2**20:.0f} MB) desde Zenodo …")
    urllib.request.urlretrieve(finfo["links"]["self"], zpath)
md5_loc = md5(zpath)
if md5_loc != md5_ref:
    raise RuntimeError(f"MD5 no coincide: local {md5_loc} ≠ Zenodo {md5_ref}. Borrar y descargar de nuevo.")
print(f"{ZIP_NAME}: MD5 verificado ({md5_loc}) · registro {RECORD} · versión {rec['metadata']['version']}")

# ------------------------------------------------------------------ 2. extracción del subconjunto
with zipfile.ZipFile(zpath) as z:
    members = [m for m in z.namelist() if GID in m and not m.endswith("/")]
    members += ["CAMELS-PE/01_metadata/data_dictionary.csv", "CAMELS-PE/README.md"]
    for m in members:
        (DEST / Path(m).name).write_bytes(z.read(m))
    st = pd.read_csv(z.open("CAMELS-PE/01_metadata/stations.csv"))
    st[st["gauge_id"] == GID].to_csv(DEST / f"{GID}_metadatos_estacion.csv", index=False)
    attrs = [pd.read_csv(z.open(m)).query("gauge_id == @GID").set_index("gauge_id")
             for m in z.namelist() if m.startswith("CAMELS-PE/02_attributes/") and m.endswith(".csv")]
    pd.concat(attrs, axis=1).T.rename(columns={GID: "valor"}).to_csv(DEST / f"{GID}_atributos.csv")
(DEST / "PROCEDENCIA.json").write_text(json.dumps({
    "dataset": "CAMELS-PE: Catchment Attributes and Meteorology for Large-sample Studies in Peru",
    "version": rec["metadata"]["version"], "doi": DOI, "zenodo_record": RECORD,
    "archivo": ZIP_NAME, "md5": md5_loc, "fecha_verificacion": pd.Timestamp.today().strftime("%Y-%m-%d"),
    "estacion": GID, "archivos_extraidos": sorted(Path(m).name for m in members),
}, indent=2, ensure_ascii=False), encoding="utf-8")
print(f"Extraídos a {DEST}: {sorted(p.name for p in DEST.iterdir())}")

# ------------------------------------------------------------------ 3. series mensuales
d = pd.read_csv(DEST / f"{GID}.csv", parse_dates=["date"]).set_index("date")
area = float(pd.read_csv(DEST / f"{GID}_atributos.csv", index_col=0).loc["area", "valor"])
assert not d.index.duplicated().any(), "fechas duplicadas"
assert (d.index == pd.date_range(d.index[0], d.index[-1], freq="D")).all(), "faltan fechas en el índice"
for v in ["prec", "flow_obs"]:
    assert (d[v].dropna() >= 0).all(), f"valores negativos en {v}"

g = d.resample("MS")
n_dias = g.size()
m = pd.DataFrame(index=n_dias.index)
m.index.name = "fecha"
m["n_dias"] = n_dias
for v, nombre in [("prec", "P"), ("flow_obs", "Q"), ("tmean", "T")]:
    m[f"dias_validos_{nombre}"] = g[v].count()
# Acumulados mensuales a partir de la resolución diaria: P_m = Σ P_d ; R_m = Σ q_d (q en mm/d,
# ya expresado en lámina, por lo que no se repite la conversión 86.4/A).
m["P_L_mm_mes"] = g["prec"].sum(min_count=1).where(m["dias_validos_P"] == m["n_dias"])
m["dias_faltantes_Q"] = m["n_dias"] - m["dias_validos_Q"]
ok_q = (m["dias_faltantes_Q"] / m["n_dias"]) <= MAX_FALT
suma_q = g["flow_obs"].sum(min_count=1)
m["R_suma_dias_validos_mm"] = suma_q.where(m["dias_validos_Q"] > 0)
# Mes completo: R = Σ q_d. Mes con ≤ 10 % de días faltantes: Σ de días válidos escalado a
# n_días (equivale a media de días válidos · n_días). Más faltantes: excluido.
m["R_mm_mes"] = (suma_q * m["n_dias"] / m["dias_validos_Q"]).where(ok_q)
m["R_flag"] = np.select([m["dias_faltantes_Q"] == 0, ok_q, m["dias_validos_Q"] == 0],
                        ["completo", "escalado", "sin_dato"], default="excluido")
m["Q_mm_d"] = g["flow_obs"].mean().where(ok_q)
m["Q_m3_s"] = m["Q_mm_d"] * area / 86.4
ok_t = (1 - m["dias_validos_T"] / m["n_dias"]) <= MAX_FALT
m["T_C"] = g["tmean"].mean().where(ok_t)
m = m.round(4)
cols = ["n_dias", "dias_validos_P", "P_L_mm_mes", "dias_validos_Q", "dias_faltantes_Q",
        "R_suma_dias_validos_mm", "R_mm_mes", "R_flag", "Q_mm_d", "Q_m3_s", "dias_validos_T", "T_C"]
m = m[cols]
m.to_csv(OUT / "series_mensuales_camels_pe.csv")
# Serie del periodo común de análisis (1998-01 .. 2025-12), la que se usa en los puntos 1–5
m.loc["1998":"2025"].to_csv(OUT / "series_mensuales_camels_pe_1998_2025.csv")

# Tabla de disponibilidad año × mes (días válidos) para cada variable
with open(OUT / "disponibilidad_anio_mes.csv", "w", encoding="utf-8") as fh:
    for nombre in ["P", "Q", "T"]:
        t = m[f"dias_validos_{nombre}"].to_frame("v")
        t = t.assign(anio=t.index.year, mes=t.index.month).pivot(index="anio", columns="mes", values="v")
        fh.write(f"# dias validos de {nombre}\n")
        t.to_csv(fh)

# ------------------------------------------------------------------ resumen
print(f"\nÁrea CAMELS-PE usada para Q en m3/s: {area:.3f} km2")
for col in ["P_L_mm_mes", "Q_m3_s", "R_mm_mes", "T_C"]:
    s = m[col]
    print(f"{col:11s} meses válidos {s.notna().sum():3d}/{len(s)} | primer {s.first_valid_index():%Y-%m}"
          f" | último {s.last_valid_index():%Y-%m}")
excl = m.loc["1998":, "Q_mm_d"].isna()
print("\nMeses de Q excluidos desde 1998 (> 10 % de días faltantes):",
      ", ".join(f"{t:%Y-%m}" for t in m.loc["1998":].index[excl]))
print("Meses de T sin dato:", f"{m['T_C'].isna().sum()} ({m.index[m['T_C'].isna()].min():%Y-%m} a "
      f"{m.index[m['T_C'].isna()].max():%Y-%m})")
print("\nMeses de R desde 1998 por condición:",
      m.loc["1998":, "R_flag"].value_counts().to_dict())
# Revisión manual (punto 1.4.2): un mes completo y uno escalado, recalculados desde los diarios
for fecha in ["2010-03", "1998-03"]:
    mes = d.loc[fecha]
    nv, nd = mes["flow_obs"].count(), len(mes)
    print(f"Revisión {fecha}: ΣP = {mes['prec'].sum():.2f} mm (tabla {m.loc[fecha, 'P_L_mm_mes'].iloc[0]})"
          f" | Σq = {mes['flow_obs'].sum():.2f} mm en {nv}/{nd} días → R = {mes['flow_obs'].sum()*nd/nv:.2f}"
          f" mm (tabla {m.loc[fecha, 'R_mm_mes'].iloc[0]}, {m.loc[fecha, 'R_flag'].iloc[0]})")
print(f"\nSalidas: {OUT}")
