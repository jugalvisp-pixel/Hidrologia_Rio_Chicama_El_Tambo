"""
Punto 1.4: control de calidad de las series de la cuenca El Tambo.

1. Disponibilidad año × mes por variable (días válidos y condición del mes) con el mismo criterio
   de completitud del análisis (script 09): P_L solo meses completos; Q̄/R válidos con ≤ 10 % de
   días faltantes (completo o escalado); T con ≤ 10 % de días faltantes; P_I y T de ERA5-Land mensuales.
2. Comprobaciones sobre los datos diarios de CAMELS-PE: fechas continuas y sin duplicados, códigos de
   faltante, negativos, coherencia tmin ≤ tmean ≤ tmax, definición de tmean, variables simuladas.
3. Búsqueda de secuencias constantes, tramos interpolados linealmente, saltos y extremos aislados en
   el caudal diario, cambios de cobertura por año.
4. Revisión manual de un mes de IMERG desde el archivo original (tasa → mm/mes y promedio ponderado).
El registro de anomalías y decisiones se escribe a partir de estas salidas (registro_anomalias_1_4.csv).

Salidas (Tarea 1/El_Tambo/resultados/control_calidad/):
  fig_disponibilidad_1_4.png, disponibilidad_condicion_mes.csv, resumen_disponibilidad_1_4.csv,
  q_secuencias_constantes.csv, q_tramos_lineales.csv, q_picos_aislados.csv
"""
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch

from periodo_comun import periodo_comun, fin_diario

ROOT = Path(__file__).resolve().parents[2] / "Tarea 1" / "El_Tambo"
SM = ROOT / "resultados" / "series_mensuales"
OUT = ROOT / "resultados" / "control_calidad"
OUT.mkdir(parents=True, exist_ok=True)
MAX_FALT = 0.10

d = pd.read_csv(ROOT / "datos" / "camels_pe" / "PE_47E0A2A8.csv", parse_dates=["date"], index_col="date")
m = pd.read_csv(SM / "series_mensuales_camels_pe.csv", parse_dates=["fecha"], index_col="fecha")
imerg = pd.read_csv(SM / "imerg_cuenca_mensual.csv", parse_dates=["fecha"], index_col="fecha")
INI, FIN = periodo_comun()                     # periodo común P_L–P_I: figuras, conteos y registro
FIN_D = fin_diario(FIN)

# ------------------------------------------------------------------ 2. comprobaciones básicas
print("=== Comprobaciones de los datos diarios (CAMELS-PE) ===")
idx = pd.date_range(d.index.min(), d.index.max(), freq="D")
print(f"Fechas: {d.index.min():%Y-%m-%d} a {d.index.max():%Y-%m-%d} | {len(d)} filas, esperadas {len(idx)} | "
      f"duplicadas {d.index.duplicated().sum()} | faltan fechas {len(idx.difference(d.index))}")
for v in ["prec", "flow_obs", "tmin", "tmean", "tmax", "pet"]:
    s = d[v]
    sosp = s[(s <= -99) | (s >= 9999)]
    print(f"  {v:8s} NaN {s.isna().sum():5d} | mín {s.min():9.3f} | máx {s.max():9.3f} | negativos "
          f"{int((s < 0).sum()) if v in ('prec', 'flow_obs', 'pet') else 'n/a':>3} | códigos tipo −999/9999: {len(sosp)}")
t = d[["tmin", "tmean", "tmax"]].dropna()
print(f"  tmin ≤ tmean ≤ tmax: {((t.tmin <= t.tmean) & (t.tmean <= t.tmax)).mean() * 100:.2f} % de los días | "
      f"|tmean − (tmin+tmax)/2| máx = {(t.tmean - (t.tmin + t.tmax) / 2).abs().max():.4f} °C")
fs = d[["flow_obs", "flow_sim"]].dropna()
print(f"  flow_sim (simulado, PISCO_ARNOVIC) NO se usa; r(obs, sim) diario = {fs.corr().iloc[0, 1]:.2f} en {len(fs)} días")

# ------------------------------------------------------------------ 1. disponibilidad año × mes
cond = pd.DataFrame(index=m.index)
cond["P_L"] = np.where(m["dias_validos_P"] == m["n_dias"], "completo",
                       np.where(m["dias_validos_P"] == 0, "sin_dato", "excluido"))
cond["Q"] = m["R_flag"]
falt_T = 1 - m["dias_validos_T"] / m["n_dias"]
cond["T"] = np.select([falt_T == 0, falt_T <= MAX_FALT, m["dias_validos_T"] == 0],
                      ["completo", "escalado", "sin_dato"], "excluido")
cond["P_I"] = np.where(imerg["P_I_mm_mes"].reindex(m.index).notna(), "completo", "sin_dato")
cond.to_csv(OUT / "disponibilidad_condicion_mes.csv")
res = []
dias = {"P_L": "prec", "Q": "flow_obs", "T": "tmean"}
PERIODOS = {"registro completo 1981-01 a 2025-12": (m.index.min(), m.index.max()),
            f"periodo común {INI:%Y-%m} a {FIN:%Y-%m}": (INI, FIN)}
for v in ["P_L", "P_I", "Q", "T"]:
    for per, (a, b) in PERIODOS.items():
        c = cond.loc[a:b, v].value_counts()
        res.append({"variable": v, "periodo": per, **{k: int(c.get(k, 0)) for k in
                    ["completo", "escalado", "excluido", "sin_dato"]},
                    "pct_dias_faltantes": 100 * d.loc[a:fin_diario(b), dias[v]].isna().mean() if v in dias else np.nan})
res = pd.DataFrame(res)
res.round(2).to_csv(OUT / "resumen_disponibilidad_1_4.csv", index=False)
print("\n=== Condición de los meses ===\n", res.round(2).to_string(index=False))

CODES = {"completo": 0, "escalado": 1, "excluido": 2, "sin_dato": 3, "fuera": 4}
cmap = ListedColormap(["#1b6ca8", "#8fbce6", "#e4572e", "#d9d9d9", "#ffffff"])
fig, axs = plt.subplots(1, 4, figsize=(15, 8), sharey=True, gridspec_kw={"wspace": 0.08})
anios = np.arange(INI.year, FIN.year + 1)
cp = cond.loc[f"{INI.year}":f"{FIN.year}"].copy()
cp.loc[(cp.index < INI) | (cp.index > FIN)] = "fuera"                     # meses fuera del periodo común
for ax, (v, titulo) in zip(axs, [("P_L", "P$_L$ (CAMELS-PE, diario)"), ("P_I", "P$_I$ (IMERG, mensual)"),
                                 ("Q", "$\\bar{Q}$ / R (CAMELS-PE, diario)"), ("T", "T (PISCOt, diario)")]):
    g = cp[v].map(CODES).to_frame("c").assign(a=cp.index.year, mm=cp.index.month)
    M = g.pivot(index="a", columns="mm", values="c").reindex(anios)
    ax.imshow(M.values, cmap=cmap, norm=BoundaryNorm(np.arange(-0.5, 5.5), 5), aspect="auto",
              extent=(0.5, 12.5, anios[-1] + 0.5, anios[0] - 0.5), interpolation="nearest")
    for f in cp.index[cp[v] == "fuera"]:
        ax.text(f.month, f.year, "×", ha="center", va="center", fontsize=7, color="#9e9e9e")
    if v in ("Q", "T"):                                        # días válidos en meses no completos
        col = "dias_validos_Q" if v == "Q" else "dias_validos_T"
        for f in cp.index[cp[v].isin(["escalado", "excluido"])]:
            ax.text(f.month, f.year, int(m.loc[f, col]), ha="center", va="center", fontsize=6, color="k")
    ax.set_xticks(range(1, 13), list("EFMAMJJASOND"))
    ax.set_title(titulo, fontsize=10)
axs[0].set_yticks(anios[::2])
axs[0].set_ylabel("Año")
fig.legend(handles=[Patch(fc=c, ec="#9e9e9e" if c == "#ffffff" else c, label=l) for c, l in zip(cmap.colors, [
    "Completo (0 días faltantes)", "Retenido con ≤ 10 % de días faltantes (Q̄ media de días válidos; R escalada)",
    "Excluido (> 10 % de días faltantes)", "Sin dato", "× Fuera del periodo común"])],
           loc="lower center", ncol=3, frameon=False, fontsize=9)
fig.suptitle(f"Disponibilidad año × mes por variable en el periodo común {INI:%Y-%m} a {FIN:%Y-%m} – cuenca El Tambo "
             "(números: días válidos en meses no completos)", fontsize=11.5, y=0.97)
fig.subplots_adjust(top=0.92, bottom=0.1, left=0.05, right=0.99)
fig.savefig(OUT / "fig_disponibilidad_1_4.png", dpi=200)
plt.close(fig)

# ------------------------------------------------------------------ 3. caudal diario: patrones sospechosos
q = d["flow_obs"].loc[INI:FIN_D]
print("\n=== Caudal diario (flow_obs, mm/d) ===")
falt_anio = q.isna().groupby(q.index.year).sum()
print("Días faltantes por año:", {int(k): int(v) for k, v in falt_anio[falt_anio > 0].items()})


def tramos(mask, minimo):
    """Tramos de días consecutivos donde mask es True, de al menos `minimo` días."""
    g = (mask != mask.shift()).cumsum()
    out = []
    for _, s in mask[mask].groupby(g[mask]):
        if len(s) >= minimo:
            out.append((s.index.min(), s.index.max(), len(s)))
    return out


# Secuencias constantes: mismo valor > 0 al menos 5 días seguidos
igual = (q.diff() == 0) & (q > 0)
const = [(a - pd.Timedelta(days=1), b, n + 1, q[a]) for a, b, n in tramos(igual, 4)]
pd.DataFrame(const, columns=["inicio", "fin", "dias", "valor_mm_d"]).to_csv(OUT / "q_secuencias_constantes.csv", index=False)
print(f"Secuencias constantes (≥ 5 días con el mismo valor > 0): {len(const)}",
      [f"{a:%Y-%m-%d}→{b:%Y-%m-%d} ({n} d)" for a, b, n, _ in const[:10]])
# Tramos lineales: segunda diferencia ≈ 0 durante ≥ 7 días con caudal cambiante (posible interpolación)
d2 = q.diff().diff().abs()
lineal = (d2 < 1e-6) & (q.diff().abs() > 1e-6)
lin = [(a - pd.Timedelta(days=2), b, n + 2) for a, b, n in tramos(lineal, 5)]
pd.DataFrame(lin, columns=["inicio", "fin", "dias"]).to_csv(OUT / "q_tramos_lineales.csv", index=False)
print(f"Tramos lineales (≥ 7 días, posible relleno por interpolación): {len(lin)}",
      [f"{a:%Y-%m-%d}→{b:%Y-%m-%d} ({n} d)" for a, b, n in lin[:10]])
# Picos aislados: día > 5× la media de los 3 días anteriores y 3 posteriores, con caudal > P90
vec = pd.concat([q.shift(k) for k in (-3, -2, -1, 1, 2, 3)], axis=1).mean(axis=1)
pico = (q > 5 * vec) & (q > q.quantile(0.9))
picos = pd.DataFrame({"q_mm_d": q[pico], "media_vecinos_mm_d": vec[pico], "P_L_mismo_dia": d["prec"][pico.index[pico]]})
picos.round(3).to_csv(OUT / "q_picos_aislados.csv")
print(f"Picos aislados (> 5× media de ±3 días y > P90): {len(picos)}\n{picos.round(2).to_string() if len(picos) else ''}")
# Saltos de nivel en estiaje: mediana de jun–sep por año (detecta cambios de sección o curva)
estiaje = q[q.index.month.isin([6, 7, 8, 9])].groupby(q.index.year[q.index.month.isin([6, 7, 8, 9])]).median()
print("Mediana de caudal jun–sep por año (mm/d):", estiaje.round(3).to_dict())

# ------------------------------------------------------------------ casos señalados en 1.3
print("\n=== Casos a revisar ===")
for per in ["1998-01", "1998-02", "1998-03"]:
    s = d.loc[per]
    print(f"{per}: P_L Σ {s.prec.sum():.0f} mm | Q días válidos {s.flow_obs.count()}/{len(s)} | "
          f"faltan {', '.join(f'{x:%d}' for x in s.index[s.flow_obs.isna()]) or '-'} | Q máx {s.flow_obs.max():.1f} mm/d")
for per in ["2005-10", "2005-11", "2005-12", "2006-01", "2016-10", "2016-11", "2016-12", "2017-01"]:
    s = d.loc[per]
    print(f"{per}: P_L {s.prec.sum():6.1f} mm | Q̄ {s.flow_obs.mean() * 2171.585 / 86.4:6.2f} m³/s | "
          f"Q mín {s.flow_obs.min() * 2171.585 / 86.4:5.2f} m³/s | días válidos Q {s.flow_obs.count()}/{len(s)}")

# ------------------------------------------------------------------ 4. revisión manual de un mes de IMERG
F = "2010-03"
nc = sorted((ROOT / "datos" / "imerg" / "manual").glob(f"*3IMERG.{F.replace('-', '')}01*"))[0]
ds = xr.open_dataset(nc, group="Grid", decode_times=False)
w = pd.read_csv(ROOT / "datos" / "imerg" / "pesos_imerg_cuenca.csv")
vals = [float(ds["precipitation"].sel(lon=lo, lat=la, method="nearest").squeeze()) for lo, la in zip(w.lon, w.lat)]
tasa = np.average(vals, weights=w.area_en_cuenca_km2)
print(f"\nRevisión manual IMERG {F} ({nc.name}): unidades '{ds['precipitation'].attrs.get('units')}' | "
      f"tasa ponderada {tasa:.5f} mm/h × 24 × 31 = {tasa * 24 * 31:.2f} mm | tabla {imerg.loc[F, 'P_I_mm_mes'].iloc[0]:.2f} mm")

# ------------------------------------------------------------------ registro de anomalías (punto 1.4.4)
# Estado: "sin problema" (comprobado, nada que corregir), "retenido" (se conserva con criterio explícito),
# "excluido", "corregido" o "incierto" (permanece la duda; se declara como limitación).
cq = cond.loc[INI:FIN, "Q"]
sep_falt = cq.isin(["escalado", "excluido"]).groupby(cq.index.month).sum()
n_esc, n_exc = int((cq == "escalado").sum()), int((cq == "excluido").sum())
rp = m.loc[INI:FIN, ["P_L_mm_mes", "R_mm_mes"]].dropna()
n_rp = int((rp["R_mm_mes"] > rp["P_L_mm_mes"]).sum())
n_pl1 = int((m.loc[INI:FIN, "P_L_mm_mes"] < 1).sum())
PER = f"periodo común {INI:%Y-%m} a {FIN:%Y-%m}"
REG = [
    ("CC-01", "Todas (diario)", "1981–2025", "Posibles fechas faltantes, duplicadas o códigos de faltante",
     "Índice diario continuo (16 436 días), 0 duplicados; mínimos/máximos sin valores tipo −999/9999; NA es el único código",
     "Ninguna", "Ninguno", "sin problema"),
    ("CC-02", "P_L, Q, PET", "1981–2025", "Valores negativos", "Conteo de negativos = 0 en prec, flow_obs y pet",
     "Ninguna", "Ninguno", "sin problema"),
    ("CC-03", "Q (flow_sim)", "1981–2025", "La base incluye caudal simulado (PISCO_ARNOVIC) junto al observado",
     f"Se verifica que solo se usa flow_obs; r(obs, sim) diario = {fs.corr().iloc[0, 1]:.2f}",
     "flow_sim no se usa en ningún análisis", "Ninguno", "sin problema"),
    ("CC-04", "Q", "1998–2025, sobre todo septiembre",
     f"Días faltantes concentrados en ago–oct casi todos los años (meses no completos por mes calendario: "
     f"{ {int(k): int(v) for k, v in sep_falt[sep_falt > 0].items()} })",
     "Mapa de disponibilidad año × mes; sin tramos lineales (interpolación) ni picos aislados en los días válidos",
     "Criterio único: ≤ 10 % faltantes → mes retenido (Q̄ media de días válidos, R escalada); > 10 % → excluido; sin relleno",
     f"{n_esc} meses escalados en el {PER} (casi todos de estiaje, caudal casi constante → error de escalado pequeño); "
     "patrón sistemático, causa no documentada (¿mantenimiento anual?)",
     "incierto"),
    ("CC-05", "Q", "2019–2021",
     f"Cambio de cobertura: 61, 56 y 24 días faltantes por año; concentra la mayoría de los {n_exc} meses "
     f"excluidos del {PER}",
     "Conteo anual de faltantes y mapa de disponibilidad", "Meses excluidos quedan vacíos; no se rellenan",
     "Menos pares P–Q en 2019–2021; vacíos a tratar en Fourier (punto 4) y tendencias (punto 3)", "incierto"),
    ("CC-06", "Q / R", "1998-02",
     "Segundo mes más lluvioso de P_L (381 mm) sin caudal: 23/28 días válidos (faltan 10, 11, 23, 26 y 27) durante El Niño 1997–98",
     "Revisión de los días faltantes en los datos diarios", "Excluido (17.9 % de días faltantes > 10 %)",
     "El pico de El Niño 1998 queda subrepresentado en Q; posible pérdida de registro en la crecida", "excluido"),
    ("CC-07", "Q / R", "1998-03",
     "Máximo del registro (Q̄ 678 m³/s, R 837 mm) con 2 días faltantes (24 y 25) y R ≈ 2 × P_L (414 mm)",
     "Revisión manual: Σq = 782.59 mm en 29/31 días → R = 836.56 mm (coincide con la tabla)",
     "Retenido como escalado",
     "R > P puede reflejar almacenamiento de feb, subestimación de P_L en el evento o extrapolación de la curva de descarga en crecida",
     "incierto"),
    ("CC-08", "Q / R vs P_L", f"{n_rp} de {len(rp)} meses del {PER}, sobre todo may–ago",
     "R_m > P_m", "Distribución por mes calendario (estiaje)",
     "Se conservan (la guía advierte que no es error por sí solo)",
     "Indica aporte de almacenamiento (flujo base) en estiaje; se discute en 1.5 y en el balance", "retenido"),
    ("CC-09", "Q", "2005-12 y 2016-12",
     "Mínimos del registro en diciembre (0.54 y 0.45 m³/s) con P_L 68–79 mm; Q̄ baja frente a oct–nov (~0.8 m³/s)",
     "Días completos (31/31 y 30/31), sin secuencias constantes ni tramos lineales; extracciones registradas por ANA "
     "en CAMELS-PE: 0.07 hm³/año (≈ 0.002 m³/s), despreciables; 1 embalse minero sin volumen",
     "Se conservan",
     "Hipótesis: fin de la recesión con primeras lluvias absorbidas por el suelo, o extracciones no registradas al inicio de la siembra",
     "incierto"),
    ("CC-10", "Q", "2007-09-11 a 2007-09-15", "Secuencia de 5 días con el mismo caudal (> 0)",
     "Búsqueda de secuencias constantes ≥ 5 días: única en el registro", "Se conserva (estiaje, caudal casi constante)",
     "Despreciable", "retenido"),
    ("CC-11", "Q", "2003–2005",
     "Medianas de estiaje (jun–sep) muy bajas: 0.047–0.133 mm/d frente a ~0.2–0.5 en otros años",
     "Comparación con P_L anual: 2003–2005 son años secos (547–573 mm)", "Se conservan",
     "Coherente con sequía; no se interpreta como salto de nivel", "sin problema"),
    ("CC-12", "P_L", "1981–2025",
     "P_L es un producto interpolado (PISCOp v2.1), no observación; no reproducible con PISCOp público (15–25 % menor)",
     "Comparación con PISCOp IRI, v2.1 update y v3.0 (scripts 08 y 13)",
     "Se mantiene como P_L principal; PISCOp v3.0 como apoyo",
     "Incertidumbre de la lluvia de referencia; se declara como limitación", "incierto"),
    ("CC-13", "P_L", f"{n_pl1} meses del {PER} (estación seca)", "P_L < 1 mm (ningún cero exacto)",
     "Conteo; meses de jun–ago", "Se conservan: son meses secos reales, no faltantes", "Ninguno", "sin problema"),
    ("CC-14", "P_I", "1998-01 a 2025-08",
     "Registro hasta 2025-08; P_I nunca < 1.97 mm y P95 de 141 mm frente a 229 mm de P_L",
     f"332 meses sin faltantes ni duplicados; unidades mm/hr; revisión manual de {F} desde el archivo original: "
     f"{tasa * 24 * 31:.2f} mm = tabla", "Ninguna (limitación del producto, no error de datos)",
     "IMERG comprime la distribución (sobreestima meses secos, subestima húmedos)", "sin problema"),
    ("CC-15", "T", "1981–2020",
     "tmean de PISCOt = (tmin + tmax)/2 exactamente (diferencia máx 0.0005 °C); tmin ≤ tmean ≤ tmax en 100 % de los días",
     "Cálculo directo sobre los datos diarios", "Se documenta la definición",
     "Explica parte del desfase de −2.9 °C frente a ERA5-Land (media de 24 h)", "sin problema"),
    ("CC-16", "T", "2021–2025",
     "Sin dato en PISCOt (termina en 2020)",
     "ERA5-Land + sesgo mensual 1981–2020; validación fuera de muestra 2011–2020: MAE 0.20 °C",
     "Valores reconstruidos en columna aparte con bandera de fuente; no entran en las estadísticas de 1.3",
     "Tramo marcado en la fig. 1.1; usar con cautela en tendencias (cambio de fuente)", "corregido"),
    ("CC-17", "P_L, R, P_I", "1998-03, 2010-03",
     "Revisión manual de agregación y conversión de unidades",
     "ΣP_d, Σq_d y escalado recalculados desde los diarios (script 09); IMERG 2010-03 desde el archivo original",
     "Ninguna", "Los valores coinciden con las tablas", "sin problema"),
]
reg = pd.DataFrame(REG, columns=["id", "variable", "fecha_o_periodo", "anomalia", "comprobacion", "decision",
                                 "efecto", "estado"])
reg.to_csv(OUT / "registro_anomalias_1_4.csv", index=False, encoding="utf-8-sig")
print(f"\nRegistro de anomalías: {len(reg)} entradas → {OUT / 'registro_anomalias_1_4.csv'}")
print(reg["estado"].value_counts().to_dict())
