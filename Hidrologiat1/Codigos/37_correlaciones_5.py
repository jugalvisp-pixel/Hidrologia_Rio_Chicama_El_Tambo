"""
Punto 5: mapas de correlación de la lluvia y el caudal de la cuenca con el clima global.

Origen: script de la compañera (repositorio anajul1a4/Hidrologia-Chicama-El-Tambo). Cambios en nuestro repositorio
(2026-10-12): conjuntos de campos por argumento (final = SST + presión al nivel del mar + humedad específica en 700 hPa,
decisión del grupo; ella = SST + presión + viento zonal en 200 hPa; alternativo = SST + q 700 + altura geopotencial en
500 hPa; los dos últimos en resultados/correlaciones/comparacion_campos/), índices E y C de Takahashi et al. (2011)
reconstruidos con ERSST v5, y el mismo umbral FDR (q < 0.10) para mapas e índices.

Campos (scripts 37a y 37b, malla 2°): SST (ERSST v5), presión al nivel del mar (ERA5 msl), humedad específica en
700 hPa (ERA5 q, enmascarada bajo el terreno), y para la comparación viento zonal en 200 hPa y altura geopotencial en
500 hPa. Variables de cuenca: P_L (CAMELS-PE) y Q̄ (El Tambo, serie rellenada); P_I (IMERG) para contrastar la
sensibilidad a la fuente de precipitación.

Definición (guía 5.2): para el mes calendario j de la respuesta de la cuenca y rezago ℓ,
    r_j(λ, φ; ℓ) = corr_{t: j(t)=j} [a_X(t), a_Y(λ, φ, t − ℓ)]
correlación entre años (todos los eneros de la cuenca con todos los eneros de la celda, …). ℓ > 0: el campo
antecede a la cuenca; con ℓ = 1, enero se empareja con diciembre del año anterior (se verifica abajo).
  · Anomalías por mes calendario y por celda respecto a una referencia fija 1998–2020 (la del punto 3). Para un mes
    fijo, restar la media no cambia Pearson: la referencia solo importa en el mapa de todos los meses.
  · Periodo base: años comunes 1998–2024 (periodo_comun.py) → 27 pares por mes. P_L también 1981–2024 (44 años)
    para la estabilidad entre subperiodos (1981–2002 frente a 2003–2024).
  · Rezago: ℓ = 0 (base). Para Q̄ además ℓ = 1, fijado de antemano por el desfase lluvia → caudal de 1 mes
    (puntos 1.5.b y 2.1); no se buscó el rezago de mayor correlación.
  · Pearson como referencia; Spearman (correlación de los rangos) en la robustez.
  · Significancia: prueba t con tamaño de muestra efectivo por celda, n_ef = n (1 − ρ_x ρ_y)/(1 + ρ_x ρ_y), con ρ la
    autocorrelación de rezago 1 entre años consecutivos de la subserie de cada mes (Bretherton et al. 1999; producto
    negativo → 0). Mínimo 20 pares por celda.
  · Comparaciones múltiples: familia = todas las celdas de un mapa (un mes, un campo, una variable) en 60°S–60°N;
    Benjamini–Hochberg con α_FDR = 0.10, el doble del α = 0.05 buscado, como recomienda Wilks (2016) para campos
    con correlación espacial. Las celdas que pasan se marcan con puntos.
  · Sin tendencia: se resta a ambas series la recta OLS contra el año (por mes calendario) y se recalcula.
  · Índices: además de los mapas, se correlaciona cada variable con Niño 1+2 y Niño 3.4 (ERSST v5, CPC) y con
    promedios de cajas de cada campo (una prueba por índice; FDR por variable sobre meses × índices).

Salidas (Tarea 1/El_Tambo/resultados/correlaciones/):
  fig_corr_{P_L,Q}_{sst,msl,u200}_5.png   12 mapas por combinación (ℓ = 0, 1998–2024)
  fig_corr_P_I_sst_5.png                  IMERG frente a SST (sensibilidad a la fuente)
  fig_corr_todos_meses_5.png              anomalías de todos los meses juntos
  fig_indices_5_4.png                     correlación por mes con índices (Niño y cajas)
  fig_robustez_5_3.png                    estabilidad de los patrones (correlación de patrones con el mapa base)
  fig_muestra_5_3.png                     pares válidos y n efectivo por celda
  correlaciones_5.nc, resumen_mapas_5.csv, indices_5_4.csv, robustez_5_3.csv
"""
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import xarray as xr
import geopandas as gpd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from scipy import stats
from shapely.affinity import translate

from periodo_comun import periodo_comun
from estilo_figuras import CMAP_CORR, TINTA, TINTA2, GRIS, SIN_DATO, aplicar_estilo

aplicar_estilo()
ROOT = Path(__file__).resolve().parents[2] / "Tarea 1" / "El_Tambo"
SM = ROOT / "resultados" / "series_mensuales"
CG = ROOT / "datos" / "campos_globales"
# Conjuntos de campos (agregado en nuestro repositorio; el método es el de la compañera):
#   final       SST + presión al nivel del mar + humedad específica en 700 hPa (decisión del grupo, 2026-10-12)
#   ella        SST + presión al nivel del mar + viento zonal en 200 hPa (elección original de la compañera)
#   alternativo SST + humedad específica en 700 hPa + altura geopotencial en 500 hPa (contraste)
CONJUNTOS = {"final": ["sst", "msl", "q700"], "ella": ["sst", "msl", "u200"], "alternativo": ["sst", "q700", "z500"]}
CONJUNTO = sys.argv[1] if len(sys.argv) > 1 else "final"
assert CONJUNTO in CONJUNTOS, f"uso: python 37_correlaciones_5.py [{'|'.join(CONJUNTOS)}]"
OUT = ROOT / "resultados" / "correlaciones"
if CONJUNTO != "final":
    OUT = OUT / "comparacion_campos" / "_".join(CONJUNTOS[CONJUNTO])
OUT.mkdir(parents=True, exist_ok=True)
INI, FIN = periodo_comun()
ANIOS_BASE = np.arange(INI.year, FIN.year + 1)                 # 1998–2024
ANIOS_LARGO = np.arange(1981, FIN.year + 1)                    # P_L: 1981–2024
REF = (1998, 2020)
ALFA, ALFA_FDR, N_MIN = 0.05, 0.10, 20   # el mismo α_FDR para mapas e índices (antes: índices con q < 0.05)
LAT_LIM = 60                                                   # dominio de pruebas y mapas: 60°S–60°N
EXTREMOS = [1998, 2017, 2023]                                  # años El Niño en los extremos (puntos 1.3 y 3)
CUENCA = (-7.70, 360 - 78.50)                                  # lat, lon (0–360) aproximados del centro de la cuenca
MESES = ["Ene", "Feb", "Mar", "Abr", "May", "Jun", "Jul", "Ago", "Sep", "Oct", "Nov", "Dic"]
SECOS = {6, 7, 8}                                              # P_L ≈ 0: r inestable

# ------------------------------------------------------------------ series de la cuenca (anomalías por mes)
cam = pd.read_csv(SM / "series_mensuales_camels_pe.csv", index_col="fecha", parse_dates=True)
imerg = pd.read_csv(SM / "imerg_cuenca_mensual.csv", index_col="fecha", parse_dates=True)["P_I_mm_mes"]
CUE = {"P_L": cam["P_L_mm_mes"], "Q": cam["Q_m3_s"], "P_I": imerg}
NOMBRE = {"P_L": "P$_L$ (CAMELS-PE)", "Q": "$\\bar{Q}$ (El Tambo)", "P_I": "P$_I$ (IMERG)"}


def anomalia_serie(s):
    s = s.dropna()
    ref = s[(s.index.year >= REF[0]) & (s.index.year <= REF[1])]
    return s - s.index.month.map(ref.groupby(ref.index.month).mean()).values


AX = {k: anomalia_serie(v) for k, v in CUE.items()}

# ------------------------------------------------------------------ campos globales (anomalías por mes y celda)
def abrir(nc, var, escala=1.0):
    ds = xr.open_dataset(CG / nc)
    ren = {k: v for k, v in {"valid_time": "time", "latitude": "lat", "longitude": "lon"}.items() if k in ds.coords}
    da = ds.rename(ren)[var].squeeze(drop=True) * escala
    da = da.drop_vars([c for c in da.coords if c not in ("time", "lat", "lon")])
    da = da.assign_coords(lon=da.lon % 360, time=da.time.values.astype("datetime64[M]").astype("datetime64[ns]"))
    da = da.sortby("lon").sortby("lat").sel(lat=slice(-LAT_LIM, LAT_LIM), time=slice("1981-01-01", "2025-12-01"))
    da = da.transpose("time", "lat", "lon")
    return da.load()


def abrir_nivel(nc, var, nivel, escala=1.0, mascara=True):
    """Campo de ERA5 en un nivel de presión (archivo del 37b con varios niveles). mascara=True: NaN en los meses y
    celdas donde el nivel queda bajo el terreno del modelo (presión en superficie < nivel)."""
    ds = xr.open_dataset(CG / nc)
    lev = "pressure_level" if "pressure_level" in ds.coords else "level"
    ds = ds.sel({lev: nivel}).drop_vars(lev)
    ds.to_netcdf(tmp := CG / f"_tmp_{var}{nivel}.nc")
    da = abrir(tmp.name, var, escala)
    tmp.unlink()
    if mascara:
        sp = abrir("ERA5_sp_mensual_2grados.nc", "sp", 0.01)                # hPa
        bajo = sp.reindex_like(da, method="nearest") < nivel
        print(f"{var}{nivel}: celdas-mes bajo el terreno enmascaradas: {int(bajo.sum())} "
              f"({100 * float(bajo.mean()):.2f} %)")
        da = da.where(~bajo)
    return da


DEF_CAMPOS = {  # cada campo se abre solo si el conjunto lo usa
    "sst": lambda: dict(da=abrir("ERSST_v5_sst_mensual_1981_2025.nc", "sst"), nombre="SST (ERSST v5)", unid="°C",
                        nivel="superficie del mar"),
    "msl": lambda: dict(da=abrir("ERA5_msl_mensual_2grados.nc", "msl", 0.01), nombre="Presión al nivel del mar (ERA5)",
                        unid="hPa", nivel="nivel del mar"),
    "u200": lambda: dict(da=abrir("ERA5_u200_mensual_2grados.nc", "u"), nombre="Viento zonal u (ERA5)", unid="m/s",
                         nivel="200 hPa (positivo hacia el este)"),
    "q700": lambda: dict(da=abrir_nivel("ERA5_q_z_500_700_mensual_2grados.nc", "q", 700, 1000.0),
                         nombre="Humedad específica q (ERA5)", unid="g/kg", nivel="700 hPa (gris: nivel bajo el terreno)"),
    "z500": lambda: dict(da=abrir_nivel("ERA5_q_z_500_700_mensual_2grados.nc", "z", 500, 1 / 9.80665, mascara=False),
                         nombre="Altura geopotencial Z = z/g₀ (ERA5)", unid="m", nivel="500 hPa"),
}
CAMPOS = {k: DEF_CAMPOS[k]() for k in CONJUNTOS[CONJUNTO]}
for k, c in CAMPOS.items():
    da = c["da"]
    ref = da.sel(time=slice(f"{REF[0]}-01-01", f"{REF[1]}-12-01"))
    c["anom"] = da.groupby("time.month") - ref.groupby("time.month").mean("time")
    c["lat"], c["lon"] = da.lat.values, da.lon.values
    print(f"{k}: {da.sizes['time']} meses ({str(da.time.values[0])[:7]} a {str(da.time.values[-1])[:7]}), "
          f"malla {len(c['lat'])}×{len(c['lon'])}, celdas sin dato: {int(np.isnan(da.values).all(0).sum())}")


def matriz(campo, anios, mes, rezago):
    """Años × celdas del campo en el mes de la cuenca menos el rezago (enero con ℓ = 1 → diciembre anterior)."""
    a = CAMPOS[campo]["anom"]
    fechas = pd.DatetimeIndex([pd.Timestamp(y, mes, 1) - pd.DateOffset(months=rezago) for y in anios])
    sub = a.sel(time=fechas.values)
    assert (pd.DatetimeIndex(sub.time.values) == fechas).all()
    return sub.values.reshape(len(anios), -1)


def serie_cuenca(var, anios, mes):
    s = AX[var]
    return np.array([s.get(pd.Timestamp(y, mes, 1), np.nan) for y in anios])


# verificación del cambio de año con ℓ = 1
f_ = pd.Timestamp(2000, 1, 1) - pd.DateOffset(months=1)
assert (f_.year, f_.month) == (1999, 12)

# ------------------------------------------------------------------ correlación por celda
def r1(m):
    """Autocorrelación de rezago 1 (entre años consecutivos) de cada columna."""
    a, b = m[:-1], m[1:]
    ok = np.isfinite(a) & np.isfinite(b)
    a, b = np.where(ok, a, np.nan), np.where(ok, b, np.nan)
    a = a - np.nanmean(a, 0); b = b - np.nanmean(b, 0)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.nansum(a * b, 0) / np.sqrt(np.nansum(a * a, 0) * np.nansum(b * b, 0))


def pearson(x, Y):
    M = np.isfinite(Y) & np.isfinite(x)[:, None]
    n = M.sum(0)
    X = np.where(M, x[:, None], np.nan)
    Yv = np.where(M, Y, np.nan)
    with np.errstate(invalid="ignore", divide="ignore"):
        xa = X - np.nanmean(X, 0); ya = Yv - np.nanmean(Yv, 0)
        r = np.nansum(xa * ya, 0) / np.sqrt(np.nansum(xa ** 2, 0) * np.nansum(ya ** 2, 0))
    r[(n < N_MIN) | ~np.isfinite(r)] = np.nan
    return r, n


def rangos(m):
    out = np.full_like(m, np.nan, dtype=float)
    for j in range(m.shape[1]):
        ok = np.isfinite(m[:, j])
        if ok.sum():
            out[ok, j] = stats.rankdata(m[ok, j])
    return out


def destendencia(m, anios):
    t = anios - anios.mean()
    out = np.full_like(m, np.nan, dtype=float)
    ok = np.isfinite(m)
    for j in range(m.shape[1]):
        o = ok[:, j]
        if o.sum() >= 3:
            b = np.polyfit(t[o], m[o, j], 1)
            out[o, j] = m[o, j] - np.polyval(b, t[o])
    return out


def mapa(var, campo, anios, mes, rezago=0, metodo="pearson", sin_tend=False):
    x = serie_cuenca(var, anios, mes)
    Y = matriz(campo, anios, mes, rezago)
    if sin_tend:
        x = destendencia(x[:, None], anios)[:, 0]; Y = destendencia(Y, anios)
    if metodo == "spearman":
        x = stats.rankdata(x); Y = rangos(Y)
    r, n = pearson(x, Y)
    rho = np.clip(r1(x[:, None])[0] * r1(Y), 0, None)
    neff = np.clip(n * (1 - rho) / (1 + rho), 3, n).astype(float)
    with np.errstate(invalid="ignore", divide="ignore"):
        t = r * np.sqrt((neff - 2) / (1 - r ** 2))
    p = 2 * stats.t.sf(np.abs(t), neff - 2)
    p[~np.isfinite(r)] = np.nan
    return r, n, neff, p


def fdr(p):
    """Benjamini–Hochberg sobre las celdas válidas; devuelve máscara de celdas que pasan."""
    ok = np.isfinite(p)
    ps = np.sort(p[ok])
    m = len(ps)
    pasa = ps <= ALFA_FDR * np.arange(1, m + 1) / m
    umbral = ps[pasa].max() if pasa.any() else -1
    return ok & (p <= umbral)


def pesos(campo):
    lat = CAMPOS[campo]["lat"]
    return np.repeat(np.cos(np.deg2rad(lat)), len(CAMPOS[campo]["lon"]))


def frac_area(campo, mask, valid):
    w = pesos(campo)
    return float((w * mask).sum() / (w * valid).sum())


def corr_patrones(campo, a, b, lat_max=40):
    """Correlación de patrones centrada y ponderada por área entre dos mapas de r (40°S–40°N)."""
    lat = np.repeat(CAMPOS[campo]["lat"], len(CAMPOS[campo]["lon"]))
    w = pesos(campo)
    ok = np.isfinite(a) & np.isfinite(b) & (np.abs(lat) <= lat_max)
    a, b, w = a[ok], b[ok], w[ok]
    a = a - np.average(a, weights=w); b = b - np.average(b, weights=w)
    return float(np.sum(w * a * b) / np.sqrt(np.sum(w * a * a) * np.sum(w * b * b)))


# ------------------------------------------------------------------ mapas base
BASE = {}            # (var, campo, rezago) → dict de arrays 12 × celdas
filas = []
COMBOS = [(v, c, 0) for v in ("P_L", "Q") for c in CAMPOS] + [("Q", c, 1) for c in CAMPOS] + [("P_I", c, 0) for c in CAMPOS]
for var, campo, rez in COMBOS:
    R, NN, NE, P, S = [], [], [], [], []
    for mes in range(1, 13):
        r, n, neff, p = mapa(var, campo, ANIOS_BASE, mes, rez)
        sig = fdr(p)
        R.append(r); NN.append(n); NE.append(neff); P.append(p); S.append(sig)
        valid = np.isfinite(r)
        lat = np.repeat(CAMPOS[campo]["lat"], len(CAMPOS[campo]["lon"]))
        lon = np.tile(CAMPOS[campo]["lon"], len(CAMPOS[campo]["lat"]))
        imx, imn = np.nanargmax(r), np.nanargmin(r)
        filas.append({"variable": var, "campo": campo, "rezago": rez, "periodo": f"{ANIOS_BASE[0]}–{ANIOS_BASE[-1]}",
                      "mes": mes, "n_pares": int(np.nanmax(n)), "n_ef_mediana": float(np.nanmedian(neff[valid])),
                      "n_ef_min": float(np.nanmin(neff[valid])),
                      "frac_area_p05": frac_area(campo, valid & (p < ALFA), valid),
                      "frac_area_FDR": frac_area(campo, sig, valid),
                      "r_max": r[imx], "lat_rmax": lat[imx], "lon_rmax": lon[imx] - 360 * (lon[imx] > 180),
                      "r_min": r[imn], "lat_rmin": lat[imn], "lon_rmin": lon[imn] - 360 * (lon[imn] > 180)})
    BASE[(var, campo, rez)] = dict(r=np.array(R), n=np.array(NN), neff=np.array(NE), p=np.array(P), sig=np.array(S))
RES = pd.DataFrame(filas)
RES.round(4).to_csv(OUT / "resumen_mapas_5.csv", index=False)
print(RES.groupby(["variable", "campo", "rezago"])[["frac_area_FDR"]].agg(["mean", "max"]).round(3))

# NetCDF con los mapas base (productos del análisis)
dsr = {}
for (var, campo, rez), d in BASE.items():
    sh = (12, len(CAMPOS[campo]["lat"]), len(CAMPOS[campo]["lon"]))
    co = {"mes": np.arange(1, 13), f"lat_{campo}": CAMPOS[campo]["lat"], f"lon_{campo}": CAMPOS[campo]["lon"]}
    dims = ("mes", f"lat_{campo}", f"lon_{campo}")
    for k in ("r", "p", "neff"):
        dsr[f"{k}_{var}_{campo}_l{rez}"] = xr.DataArray(d[k].reshape(sh).astype("float32"), coords=co, dims=dims)
    dsr[f"fdr_{var}_{campo}_l{rez}"] = xr.DataArray(d["sig"].reshape(sh).astype("int8"), coords=co, dims=dims)
xr.Dataset(dsr, attrs={"descripcion": "Correlación entre años de las anomalías mensuales de la cuenca y de cada "
                                      "celda (script 37). fdr = 1: pasa Benjamini–Hochberg α_FDR = 0.10 por mapa.",
                       "periodo": f"{ANIOS_BASE[0]}–{ANIOS_BASE[-1]}", "referencia_anomalias": f"{REF[0]}–{REF[1]}"}
           ).to_netcdf(OUT / "correlaciones_5.nc", encoding={k: {"zlib": True} for k in dsr})

# ------------------------------------------------------------------ figuras de 12 mapas
paises = gpd.read_file(ROOT / "datos" / "ne_50m_admin_0_countries.geojson")
costas = pd.concat([paises.geometry.boundary, paises.geometry.boundary.apply(lambda g: translate(g, xoff=360))])
NIV = np.arange(-1, 1.01, 0.2)
NORM = BoundaryNorm(NIV, CMAP_CORR.N)


def dibujar(ax, campo, r, sig, titulo):
    lat, lon = CAMPOS[campo]["lat"], CAMPOS[campo]["lon"]
    ax.set_facecolor(SIN_DATO)                                       # tierra en SST, celdas inválidas
    R = r.reshape(len(lat), len(lon))
    im = ax.pcolormesh(lon, lat, np.ma.masked_invalid(R), cmap=CMAP_CORR, norm=NORM, shading="nearest")
    S = sig.reshape(len(lat), len(lon))
    jj, ii = np.where(S)
    ax.scatter(lon[ii], lat[jj], s=0.6, c=TINTA, marker=".", lw=0)
    costas.plot(ax=ax, color="#3d3d3d", lw=0.3)
    ax.plot(CUENCA[1], CUENCA[0], marker="*", ms=9, mfc="#ffd400", mec=TINTA, mew=0.7, zorder=5)
    ax.set_xlim(0, 360); ax.set_ylim(-LAT_LIM, LAT_LIM); ax.set_aspect("equal")
    ax.set_xticks(range(0, 361, 60), ["0°", "60°E", "120°E", "180°", "120°O", "60°O", "0°"], fontsize=6.5)
    ax.set_yticks([-60, -30, 0, 30, 60], ["60°S", "30°S", "0°", "30°N", "60°N"], fontsize=6.5)
    ax.set_xlabel(""); ax.set_ylabel("")                             # geopandas pone etiquetas del CRS
    ax.set_title(titulo, fontsize=9, loc="left", pad=2)
    return im


def figura12(var, campo, rez, d, archivo, nota=""):
    c = CAMPOS[campo]
    fig, axs = plt.subplots(4, 3, figsize=(16, 8.4), gridspec_kw={"hspace": 0.32, "wspace": 0.06})
    sub = RES[(RES.variable == var) & (RES.campo == campo) & (RES.rezago == rez)].set_index("mes")
    for mes, ax in zip(range(1, 13), axs.flat):
        f = sub.loc[mes]
        seco = " · mes seco: P≈0, r inestable" if var != "Q" and mes in SECOS else ""
        im = dibujar(ax, campo, d["r"][mes - 1], d["sig"][mes - 1],
                     f"{MESES[mes - 1]} · área FDR {100 * f.frac_area_FDR:.0f} % · n = {int(f.n_pares)}{seco}")
    cb = fig.colorbar(im, ax=axs, orientation="horizontal", fraction=0.025, pad=0.04, ticks=NIV[::2])
    cb.set_label("r de Pearson entre años (anomalías mensuales)")
    rz = "sin rezago (ℓ = 0)" if rez == 0 else f"ℓ = {rez}: el campo antecede {rez} mes a la cuenca"
    fig.suptitle(f"Punto 5 · {NOMBRE[var]} frente a {c['nombre']}, {c['nivel']} · {ANIOS_BASE[0]}–{ANIOS_BASE[-1]} · {rz}\n"
                 f"Cada mapa: correlación de los {len(ANIOS_BASE)} valores de ese mes calendario. Puntos: pasa FDR "
                 f"(Benjamini–Hochberg, α = {ALFA_FDR}, n efectivo por persistencia). Gris: tierra o sin dato. "
                 f"★ cuenca.{nota}", fontsize=10, y=0.98)
    fig.savefig(OUT / archivo, dpi=150, bbox_inches="tight")
    plt.close(fig)


for var, campo, rez in [(v, c, 0) for v in ("P_L", "Q") for c in CAMPOS]:
    figura12(var, campo, rez, BASE[(var, campo, rez)], f"fig_corr_{var}_{campo}_5.png")
figura12("P_I", "sst", 0, BASE[("P_I", "sst", 0)], "fig_corr_P_I_sst_5.png",
         " Sensibilidad a la fuente: comparar con P$_L$ frente a SST.")

# ------------------------------------------------------------------ todos los meses juntos (anomalías)
CORTO = {"sst": "SST", "msl": "SLP", "u200": "u 200 hPa", "q700": "q 700 hPa", "z500": "Z 500 hPa"}
fig, axs = plt.subplots(2, 3, figsize=(16, 5.6))
filas_tm = []
for i, var in enumerate(("P_L", "Q")):
    for j, campo in enumerate(CAMPOS):
        fechas = pd.date_range(f"{ANIOS_BASE[0]}-01-01", f"{ANIOS_BASE[-1]}-12-01", freq="MS")
        x = AX[var].reindex(fechas).values
        Y = CAMPOS[campo]["anom"].sel(time=fechas.values).values.reshape(len(fechas), -1)
        r, n = pearson(x, Y)
        rho = np.clip(r1(x[:, None])[0] * r1(Y), 0, None)
        neff = np.clip(n * (1 - rho) / (1 + rho), 3, n)
        with np.errstate(invalid="ignore", divide="ignore"):
            p = 2 * stats.t.sf(np.abs(r * np.sqrt((neff - 2) / (1 - r ** 2))), neff - 2)
        p[~np.isfinite(r)] = np.nan
        sig = fdr(p)
        valid = np.isfinite(r)
        filas_tm.append({"variable": var, "campo": campo, "n_meses": len(fechas), "n_ef_mediana": float(np.nanmedian(neff[valid])),
                         "frac_area_FDR": frac_area(campo, sig, valid), "r_max": np.nanmax(r), "r_min": np.nanmin(r)})
        im = dibujar(axs[i, j], campo, r, sig, f"{NOMBRE[var]} · {CORTO[campo]} · área FDR "
                                               f"{100 * frac_area(campo, sig, valid):.0f} %")
cb = fig.colorbar(im, ax=axs, orientation="horizontal", fraction=0.04, pad=0.08, ticks=NIV[::2])
cb.set_label("r de Pearson (anomalías respecto a la climatología mensual 1998–2020)")
fig.suptitle(f"Punto 5 · Todos los meses juntos en anomalías ({ANIOS_BASE[0]}-01 a {ANIOS_BASE[-1]}-12, n = 324; ℓ = 0). "
             "Anomalías para que el ciclo anual compartido no domine la asociación. Puntos: FDR α = 0.10.", fontsize=10)
fig.savefig(OUT / "fig_corr_todos_meses_5.png", dpi=150, bbox_inches="tight")
plt.close(fig)
pd.DataFrame(filas_tm).round(4).to_csv(OUT / "todos_meses_5.csv", index=False)

# ------------------------------------------------------------------ índices: Niño y cajas de cada campo
nino = pd.read_csv(ROOT / "datos" / "indices_clima" / "ersst5.nino.mth.91-20.ascii", sep=r"\s+")
nino.index = pd.to_datetime(dict(year=nino.YR, month=nino.MON, day=1))
IND = {"Niño 1+2 (SST)": nino.iloc[:, 3], "Niño 3.4 (SST)": nino.iloc[:, 9]}


def indices_EC():
    """Índices E y C de Takahashi et al. (2011, GRL 38, L10704), reconstruidos con ERSST v5: EOF de las anomalías de SST
    del Pacífico tropical (10°S–10°N, 140°E–80°O, pesos √cos φ), PC1 y PC2 estandarizadas y rotadas 45°:
    E = (PC1 − PC2)/√2 (calentamiento del Pacífico este), C = (PC1 + PC2)/√2 (Pacífico central). No es la serie
    oficial del IGP (no se encontró publicada): EOF calculadas en 1981–2025, no en el registro largo del artículo.
    Signos fijados con Niño 3.4 (PC1 > 0 en El Niño) y para que E siga a Niño 1+2 más que C."""
    a = CAMPOS["sst"]["anom"].sel(lat=slice(-10, 10), lon=slice(140, 280))
    w = np.sqrt(np.cos(np.deg2rad(a.lat.values)))[None, :, None]
    X = (a.values * w).reshape(a.sizes["time"], -1)
    X = X[:, np.isfinite(X).all(0)]
    X = X - X.mean(0)
    U, s, _ = np.linalg.svd(X, full_matrices=False)
    pc = U[:, :2] * s[:2]
    pc = (pc - pc.mean(0)) / pc.std(0)
    t = pd.DatetimeIndex(a.time.values)
    n34 = nino.iloc[:, 9].reindex(t).values
    n12 = nino.iloc[:, 3].reindex(t).values
    ok = np.isfinite(n34)
    if np.corrcoef(pc[ok, 0], n34[ok])[0, 1] < 0:
        pc[:, 0] *= -1
    E, C = (pc[:, 0] - pc[:, 1]) / np.sqrt(2), (pc[:, 0] + pc[:, 1]) / np.sqrt(2)
    if np.corrcoef(E[ok], n12[ok])[0, 1] < np.corrcoef(C[ok], n12[ok])[0, 1]:
        E, C = C, E                                             # cambio de signo de PC2
    var = 100 * s[:2] ** 2 / (s ** 2).sum()
    out = pd.DataFrame({"E": E, "C": C}, index=t)
    n4 = nino.iloc[:, 7].reindex(t)
    print(f"E y C: EOF1 {var[0]:.0f} %, EOF2 {var[1]:.0f} % de la varianza; r(E, Niño 1+2) = "
          f"{out.E.corr(pd.Series(n12, t)):.2f}, r(C, Niño 4) = {out.C.corr(n4):.2f}, r(E, C) = {out.E.corr(out.C):.2f}")
    out.round(4).to_csv(OUT / "indices_E_C_5.csv", index_label="fecha")
    return out


_ec = indices_EC()
IND["E: Pacífico este (Takahashi 2011)"] = _ec.E
IND["C: Pacífico central (Takahashi 2011)"] = _ec.C
CAJAS_TODAS = {  # nombre: (campo, lat_s, lat_n, lon_o, lon_e) en 0–360; cada conjunto usa las de sus campos
    "Atlántico tropical norte (SST)": ("sst", 5, 25, 300, 340),
    "Anticiclón del Pacífico sur (SLP)": ("msl", -35, -20, 260, 280),
    "Pacífico ecuatorial oeste (SLP)": ("msl", -10, 10, 120, 150),
    "Andes norte de Perú (u200)": ("u200", -15, 0, 280, 290),
    "Pacífico ecuatorial central (u200)": ("u200", -5, 5, 180, 220),
    "Amazonía occidental (q 700)": ("q700", -10, 0, 285, 295),      # fuente de humedad al este de los Andes
    "Costa y Pacífico este (q 700)": ("q700", -10, 0, 270, 280),     # humedad sobre el mar frente a la cuenca
    "Pacífico ecuatorial este (Z 500)": ("z500", -5, 5, 240, 280),   # columna cálida de El Niño
    "Pacífico subtropical sur (Z 500)": ("z500", -35, -20, 260, 280),  # anticiclón en altura media
}
CAJAS = {k: v for k, v in CAJAS_TODAS.items() if v[0] in CAMPOS}
for nom, (campo, la0, la1, lo0, lo1) in CAJAS.items():
    a = CAMPOS[campo]["anom"].sel(lat=slice(la0, la1), lon=slice(lo0, lo1))
    w = np.cos(np.deg2rad(a.lat))
    IND[nom] = a.weighted(w).mean(("lat", "lon")).to_series()
filas_i = []
for var in ("P_L", "Q", "P_I"):
    for nom, s in IND.items():
        for mes in range(1, 13):
            x = serie_cuenca(var, ANIOS_BASE, mes)
            y = np.array([s.get(pd.Timestamp(yy, mes, 1), np.nan) for yy in ANIOS_BASE])
            ok = np.isfinite(x) & np.isfinite(y)
            r = stats.pearsonr(x[ok], y[ok])[0]
            rs = stats.spearmanr(x[ok], y[ok])[0]
            rho = max(r1(x[ok, None])[0] * r1(y[ok, None])[0], 0)
            neff = max(min(ok.sum() * (1 - rho) / (1 + rho), ok.sum()), 3)
            p = 2 * stats.t.sf(abs(r) * np.sqrt((neff - 2) / (1 - r ** 2)), neff - 2)
            filas_i.append({"variable": var, "indice": nom, "mes": mes, "n": int(ok.sum()), "n_ef": neff,
                            "pearson_r": r, "spearman_rho": rs, "p": p})
IDX = pd.DataFrame(filas_i)
IDX["q_FDR"] = np.nan
for var, g in IDX.groupby("variable"):                         # familia: meses × índices de cada variable
    IDX.loc[g.index, "q_FDR"] = stats.false_discovery_control(g["p"].values)
IDX.round(4).to_csv(OUT / "indices_5_4.csv", index=False)

fig, axs = plt.subplots(1, 3, figsize=(16, 4.4), sharey=True)
for ax, var in zip(axs, ("P_L", "Q", "P_I")):
    M = IDX[IDX.variable == var].pivot(index="indice", columns="mes", values="pearson_r").reindex(list(IND))
    Qv = IDX[IDX.variable == var].pivot(index="indice", columns="mes", values="q_FDR").reindex(list(IND))
    Pv = IDX[IDX.variable == var].pivot(index="indice", columns="mes", values="p").reindex(list(IND))
    im = ax.imshow(M.values, cmap=CMAP_CORR, norm=NORM, aspect="auto")
    for (i, j), v in np.ndenumerate(M.values):
        marca = "**" if Qv.values[i, j] < ALFA_FDR else ("*" if Pv.values[i, j] < ALFA else "")
        ax.text(j, i, f"{v:+.2f}{marca}", ha="center", va="center", fontsize=6.3,
                color="white" if abs(v) > 0.6 else TINTA)
    ax.set_xticks(range(12), MESES, fontsize=8)
    ax.set_yticks(range(len(IND)), list(IND), fontsize=8)
    ax.set_title(f"{NOMBRE[var]} · {ANIOS_BASE[0]}–{ANIOS_BASE[-1]}", fontsize=10, loc="left")
fig.colorbar(im, ax=axs, fraction=0.015, pad=0.01, ticks=NIV).set_label("r de Pearson")
fig.suptitle("Punto 5.4 · Correlación por mes calendario con índices (ℓ = 0). * p < 0.05 con n efectivo; "
             f"** además q FDR < {ALFA_FDR:.2f}, el mismo umbral de los mapas (familia: meses × índices de cada variable)",
             fontsize=10)
fig.savefig(OUT / "fig_indices_5_4.png", dpi=150, bbox_inches="tight")
plt.close(fig)

# ------------------------------------------------------------------ robustez (5.3)
anios_sin_ext = np.array([y for y in ANIOS_BASE if y not in EXTREMOS])
mitad1, mitad2 = np.arange(1981, 2003), np.arange(2003, FIN.year + 1)
filas_r = []
for var in ("P_L", "Q"):
    for campo in CAMPOS:
        b = BASE[(var, campo, 0)]
        for mes in range(1, 13):
            rb = b["r"][mes - 1]
            variantes = {
                "Spearman": mapa(var, campo, ANIOS_BASE, mes, metodo="spearman")[0],
                "sin tendencia": mapa(var, campo, ANIOS_BASE, mes, sin_tend=True)[0],
                f"sin {', '.join(map(str, EXTREMOS))}": mapa(var, campo, anios_sin_ext, mes)[0],
            }
            if var == "P_L":
                variantes["IMERG en lugar de P_L"] = BASE[("P_I", campo, 0)]["r"][mes - 1]
                variantes["1981–2002 vs 2003–2024"] = None
            else:
                variantes["rezago ℓ = 1"] = BASE[("Q", campo, 1)]["r"][mes - 1]
            for nom, rv in variantes.items():
                if rv is None:
                    ra = mapa(var, campo, mitad1, mes)[0]; rb2 = mapa(var, campo, mitad2, mes)[0]
                    cp = corr_patrones(campo, ra, rb2)
                else:
                    cp = corr_patrones(campo, rb, rv)
                filas_r.append({"variable": var, "campo": campo, "mes": mes, "variante": nom, "corr_patrones": cp})
ROB = pd.DataFrame(filas_r)
ROB.round(3).to_csv(OUT / "robustez_5_3.csv", index=False)

orden_v = list(dict.fromkeys(ROB.variante))
fig, axs = plt.subplots(2, 3, figsize=(17, 6.4), gridspec_kw={"wspace": 0.05})
for i, var in enumerate(("P_L", "Q")):
    for j, campo in enumerate(CAMPOS):
        ax = axs[i, j]
        M = ROB[(ROB.variable == var) & (ROB.campo == campo)].pivot(index="variante", columns="mes", values="corr_patrones")
        M = M.reindex([v for v in orden_v if v in M.index])
        im = ax.imshow(M.values, cmap=CMAP_CORR, norm=NORM, aspect="auto")
        for (a, b_), v in np.ndenumerate(M.values):
            ax.text(b_, a, f"{v:.1f}", ha="center", va="center", fontsize=6.5, color="white" if abs(v) > 0.6 else TINTA)
        ax.set_xticks(range(12), MESES, fontsize=7.5)
        ax.set_yticks(range(len(M)), M.index if j == 0 else [""] * len(M), fontsize=7.5)
        ax.set_title(f"{NOMBRE[var]} · {CAMPOS[campo]['nombre']}", fontsize=9, loc="left")
fig.colorbar(im, ax=axs, fraction=0.015, pad=0.01, ticks=NIV).set_label("correlación de patrones con el mapa base")
fig.suptitle("Punto 5.3 · Estabilidad de los mapas: correlación de patrones (ponderada por área, 40°S–40°N) entre el mapa "
             "base (Pearson, 1998–2024, ℓ = 0) y cada variante.\nP$_L$ 1981–2002 vs 2003–2024 compara los dos subperiodos "
             "del registro largo entre sí. Valores cercanos a 1: el patrón no depende de la elección.", fontsize=10)
fig.savefig(OUT / "fig_robustez_5_3.png", dpi=150, bbox_inches="tight")
plt.close(fig)

# ------------------------------------------------------------------ tamaño de muestra por celda
fig, axs = plt.subplots(1, 2, figsize=(16, 3.6))
d = BASE[("P_L", "sst", 0)]
for ax, (arr, tit, vmax) in zip(axs, ((d["n"][2], "pares válidos por celda (P$_L$–SST, marzo)", len(ANIOS_BASE)),
                                      (d["neff"][2], "n efectivo por celda (P$_L$–SST, marzo)", len(ANIOS_BASE)))):
    lat, lon = CAMPOS["sst"]["lat"], CAMPOS["sst"]["lon"]
    A = np.where(np.isfinite(d["r"][2]), arr, np.nan).reshape(len(lat), len(lon))
    ax.set_facecolor(SIN_DATO)
    im = ax.pcolormesh(lon, lat, np.ma.masked_invalid(A), cmap="Greys", vmin=0, vmax=vmax, shading="nearest")
    costas.plot(ax=ax, color="#3d3d3d", lw=0.3)
    ax.set_xlim(0, 360); ax.set_ylim(-LAT_LIM, LAT_LIM); ax.set_aspect("equal")
    ax.set_xlabel(""); ax.set_ylabel("")
    ax.set_title(tit, fontsize=9.5, loc="left")
    fig.colorbar(im, ax=ax, fraction=0.025, pad=0.01)
fig.suptitle(f"Punto 5.3 · Tamaño de muestra: {len(ANIOS_BASE)} años en todas las celdas válidas (gris = tierra, hielo marino o "
             f"menos de {N_MIN} pares); el n efectivo baja donde la SST persiste de un año a otro.", fontsize=10)
fig.savefig(OUT / "fig_muestra_5_3.png", dpi=150, bbox_inches="tight")
plt.close(fig)

print(f"\nÍndices con q FDR < {ALFA_FDR}:")
print(IDX[IDX.q_FDR < ALFA_FDR][["variable", "indice", "mes", "pearson_r", "n_ef", "q_FDR"]].round(3).to_string(index=False))
print("\nRobustez (mediana de la correlación de patrones por variante):")
print(ROB.groupby(["variable", "campo", "variante"]).corr_patrones.median().round(2).unstack().to_string())
print(f"\nSalidas en {OUT}")
