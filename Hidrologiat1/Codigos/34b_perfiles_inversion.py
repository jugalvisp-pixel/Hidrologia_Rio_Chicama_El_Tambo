"""
Inversión térmica de la costa norte del Perú y su posición respecto a la cuenca El Tambo (apoyo del punto 1.5.c).

Datos (script 34a):
  - Radiosondeos observados de Trujillo, IGRA2 PEM00084516 (−8.04°, −79.06°, 127 m), 12 UTC = 7 a.m. local.
  - ERA5 mensual por niveles de presión, promedio de las 12 UTC, celda más cercana a Trujillo (periodo común).

Método:
  - Cada nivel del sondeo con altura geopotencial (GPH) y temperatura válidas; Td = T − depresión del punto de
    rocío; HR = 100·e(Td)/e_s(T) con la fórmula de Magnus (Bolton, 1980).
  - Inversión: desde la superficie hacia arriba, primer nivel a partir del cual T aumenta con la altura (base) y
    último nivel de ese aumento continuo (tope), por debajo de 3 500 m. Se cuenta como inversión si el
    calentamiento base→tope es ≥ 1 °C. Capa de nubes: tramo con HR ≥ 95 % (interpolada cada 10 m) entre la superficie y el tope.
  - Perfiles compuestos: cada sondeo se interpola linealmente a una malla de 25 m (sin extrapolar) y se calculan
    mediana y rango intercuartílico por temporada (dic–abr húmeda, may–sep seca, oct–nov transición).
  - Sondeo típico: el de la estación seca cuya base, tope e intensidad están más cerca de las medianas de esa
    estación (distancia normalizada por el rango intercuartílico).
  - Parcela: asciende desde la superficie por la adiabática seca (9.8 °C/km) hasta el nivel de condensación por
    ascenso (NCA, iterando T y Td de la parcela) y desde ahí por la adiabática húmeda (pseudoadiabática, Γ_m
    calculada con la razón de mezcla de saturación a la presión del sondeo).

Salidas (Tarea 1/El_Tambo/resultados/inversion/):
  fig_sondeo_tipico_inversion.png, fig_perfiles_inversion.png,
  inversion_por_sondeo.csv, inversion_resumen_mensual.csv, perfiles_compuestos_temporada.csv
"""
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

from periodo_comun import periodo_comun
from estilo_figuras import EST

ROOT = Path(__file__).resolve().parents[2] / "Tarea 1" / "El_Tambo"
TXT = ROOT / "datos" / "radiosondeo" / "PEM00084516-data.txt"
NC = ROOT / "datos" / "era5_perfiles" / "ERA5_perfiles_12UTC_mensual.nc"
MORF = ROOT / "resultados" / "delimitacion" / "parametros_morfometricos.csv"
OUT = ROOT / "resultados" / "inversion"
OUT.mkdir(parents=True, exist_ok=True)
INI, FIN = periodo_comun()

Z_MAX_INV = 3500.0           # m s.n.m., búsqueda de la inversión
UMBRAL_INV = 1.0             # °C de calentamiento base→tope
HR_NUBE = 95.0               # %, capa de nubes (HR interpolada cada 10 m)
MALLA = np.arange(100, 6001, 25.0)
TEMPORADAS = {"dic–abr (húmeda)": [12, 1, 2, 3, 4], "may–sep (seca)": [5, 6, 7, 8, 9], "oct–nov": [10, 11]}
COL_T = {"dic–abr (húmeda)": EST["húmeda"], "may–sep (seca)": EST["seca"], "oct–nov": EST["transición"]}
MESES = ["Ene", "Feb", "Mar", "Abr", "May", "Jun", "Jul", "Ago", "Sep", "Oct", "Nov", "Dic"]
G, CP, RD, LV, EPS = 9.80665, 1004.0, 287.0, 2.501e6, 0.622

morf = pd.read_csv(MORF, index_col=0)["valor"]
Z_CIERRE, Z_MED, Z_MAXC = (float(morf["Cota cierre (m)"]), float(morf["Cota mediana (m s.n.m.)"]),
                           float(morf["Cota máxima (m s.n.m.)"]))


def temporada(mes):
    return next(k for k, v in TEMPORADAS.items() if mes in v)


def es(t):                                                    # hPa, t en °C (Bolton 1980)
    return 6.112 * np.exp(17.67 * t / (t + 243.5))


# ------------------------------------------------------------------ 1. lectura de IGRA2
def num(s, escala=1.0):
    s = s.strip()
    if not s or s in ("-9999", "-8888"):
        return np.nan
    return float(s) / escala


sondeos = []
with TXT.open(encoding="ascii") as f:
    cab, niveles = None, []
    for linea in f:
        if linea.startswith("#"):
            if cab is not None:
                sondeos.append((cab, niveles))
            cab = pd.Timestamp(int(linea[13:17]), int(linea[18:20]), int(linea[21:23]), int(linea[24:26]))
            niveles = []
        else:
            niveles.append({"p_hPa": num(linea[9:15], 100), "z_m": num(linea[16:21]), "T_C": num(linea[22:27], 10),
                            "dpd_C": num(linea[34:39], 10), "superficie": linea[1] == "1"})
    sondeos.append((cab, niveles))

def alturas(d):
    """Altitud de todos los niveles: hipsométrica desde la superficie, ajustada a la GPH reportada.

    En IGRA2 los niveles significativos (los que describen la inversión) suelen venir sin GPH; solo la traen los
    niveles estándar. z_i = z_{i−1} + (R_d·T̄_v/g)·ln(p_{i−1}/p_i), con T_v de T y Td; luego se suma la diferencia
    con la GPH reportada, interpolada en ln p.
    """
    d = d[d.p_hPa.notna()].sort_values("p_hPa", ascending=False).drop_duplicates("p_hPa")
    sup = d[d.superficie]
    if sup.empty:
        return None
    p0 = sup.p_hPa.iloc[0]
    d = d[d.p_hPa <= p0].reset_index(drop=True)                 # fuera los niveles bajo la superficie
    lnp = np.log(d.p_hPa.values)
    ok = d.T_C.notna().values
    if ok.sum() < 5 or not ok[0]:
        return None
    t = np.interp(-lnp, -lnp[ok], d.T_C.values[ok])
    td = d.T_C.values - d.dpd_C.values
    q = np.where(np.isfinite(td), EPS * es(np.where(np.isfinite(td), td, 0)) / d.p_hPa.values, 0)
    tv = (t + 273.15) * (1 + 0.61 * q)
    z = np.empty(len(d))
    z[0] = sup.z_m.iloc[0] if pd.notna(sup.z_m.iloc[0]) else 127.0
    for i in range(1, len(d)):
        z[i] = z[i - 1] + RD * 0.5 * (tv[i] + tv[i - 1]) / G * (lnp[i - 1] - lnp[i])
    rep = d.z_m.notna().to_numpy().copy()
    rep[0] = False
    if rep.any():
        z = z + np.interp(-lnp, -lnp[rep], (d.z_m.values - z)[rep], left=0)
    d["z_m"] = z
    return d[ok]


filas, perfiles, crudos = [], {}, {}
n_sin = 0
for fecha, niv in sondeos:
    d = alturas(pd.DataFrame(niv))
    if d is None:
        n_sin += 1
        continue
    d = d.sort_values("z_m").drop_duplicates("z_m")
    d["Td_C"] = d.T_C - d.dpd_C
    d["HR"] = (100 * es(d.Td_C) / es(d.T_C)).clip(upper=100)
    crudos[fecha] = d
    # inversión
    b = d[d.z_m <= Z_MAX_INV].reset_index(drop=True)
    fila = {"fecha": fecha, "mes": fecha.month, "temporada": temporada(fecha.month), "z_sup_m": b.z_m.iloc[0],
            "T_sup_C": b.T_C.iloc[0], "HR_sup": b.HR.iloc[0], "inversion": False}
    k = next((i for i in range(len(b) - 1) if b.T_C[i + 1] > b.T_C[i]), None)
    if k is not None:
        j = k + 1
        while j + 1 < len(b) and b.T_C[j + 1] >= b.T_C[j]:
            j += 1
        dT = b.T_C[j] - b.T_C[k]
        sobre = d[(d.z_m > b.z_m[j]) & (d.z_m <= b.z_m[j] + 500)]
        bajo = d[d.z_m <= b.z_m[k]]
        zf = np.arange(d.z_m.iloc[0], b.z_m[j] + 1, 10.0)
        nube = zf[np.interp(zf, d.z_m, d.HR) >= HR_NUBE]
        fila.update(inversion=dT >= UMBRAL_INV, base_m=b.z_m[k], tope_m=b.z_m[j], T_base_C=b.T_C[k],
                    T_tope_C=b.T_C[j], intensidad_C=dT, espesor_m=b.z_m[j] - b.z_m[k],
                    gradiente_C_km=1000 * dT / max(b.z_m[j] - b.z_m[k], 1),
                    HR_bajo_base=bajo.HR.mean(), HR_sobre_tope=sobre.HR.mean() if len(sobre) else np.nan,
                    nube_base_m=nube.min() if len(nube) else np.nan,
                    nube_tope_m=nube.max() if len(nube) else np.nan)
    filas.append(fila)
    perfiles[fecha] = (np.interp(MALLA, d.z_m, d.T_C, left=np.nan, right=np.nan),
                       np.interp(MALLA, d.z_m, d.HR, left=np.nan, right=np.nan))

inv = pd.DataFrame(filas).set_index("fecha")
for c in ("base_m", "tope_m", "intensidad_C", "espesor_m"):           # solo cuentan las inversiones ≥ umbral
    inv.loc[~inv.inversion, c] = np.nan
inv.round(2).to_csv(OUT / "inversion_por_sondeo.csv")
print(f"Sondeos sin superficie o con < 5 niveles: {n_sin}")
print(f"Sondeos válidos: {len(inv)} ({inv.index.min():%Y-%m-%d} a {inv.index.max():%Y-%m-%d}); "
      f"con inversión ≥ {UMBRAL_INV} °C: {inv.inversion.sum()} ({100 * inv.inversion.mean():.0f} %)")

# resumen mensual
q = lambda p: (lambda s: s.quantile(p))
res = inv.groupby("mes").agg(
    n=("inversion", "size"), frecuencia_pct=("inversion", lambda s: 100 * s.mean()),
    base_med_m=("base_m", "median"), base_q1_m=("base_m", q(.25)), base_q3_m=("base_m", q(.75)),
    tope_med_m=("tope_m", "median"), tope_q1_m=("tope_m", q(.25)), tope_q3_m=("tope_m", q(.75)),
    intensidad_med_C=("intensidad_C", "median"), HR_bajo_base=("HR_bajo_base", "median"),
    HR_sobre_tope=("HR_sobre_tope", "median")).reindex(range(1, 13))
res.insert(0, "mes_nombre", MESES)
res.round(1).to_csv(OUT / "inversion_resumen_mensual.csv")
print(res.round(0).to_string())
for t in TEMPORADAS:
    s = inv[inv.temporada == t]
    print(f"{t:18s} n={len(s):3d} frecuencia {100 * s.inversion.mean():3.0f} % | base {s.base_m.median():5.0f} m | "
          f"tope {s.tope_m.median():5.0f} m | intensidad {s.intensidad_C.median():.1f} °C | "
          f"HR bajo base {s.HR_bajo_base.median():.0f} % / sobre tope {s.HR_sobre_tope.median():.0f} %")

# ------------------------------------------------------------------ 2. perfiles compuestos por temporada
fechas = np.array(list(perfiles))
TT = np.array([perfiles[f][0] for f in fechas])
HH = np.array([perfiles[f][1] for f in fechas])
# Perfiles de cada sondeo en la malla de 25 m (para el contraste con ERA5 en los mismos meses, script 34c)
xr.Dataset({"T_C": (("fecha", "z_m"), TT), "HR": (("fecha", "z_m"), HH)},
           coords={"fecha": pd.DatetimeIndex(fechas), "z_m": MALLA}).to_netcdf(OUT / "perfiles_por_sondeo.nc")
comp = {}
for t, meses in TEMPORADAS.items():
    m = np.isin([f.month for f in fechas], meses)
    ok = (np.isfinite(TT[m]).sum(0) >= max(10, 0.5 * m.sum()))            # al menos la mitad de los sondeos
    comp[t] = {"n": int(m.sum())}
    for nom, A in (("T", TT[m]), ("HR", HH[m])):
        for p, etq in ((25, "q1"), (50, "med"), (75, "q3")):
            v = np.nanpercentile(A, p, axis=0)
            v[~ok] = np.nan
            comp[t][f"{nom}_{etq}"] = v
pd.DataFrame({"z_m": MALLA, **{f"{t}|{k}": v for t, c in comp.items() for k, v in c.items() if k != "n"}}
             ).round(2).to_csv(OUT / "perfiles_compuestos_temporada.csv", index=False)

# ------------------------------------------------------------------ 3. ERA5 12 UTC en la costa de Trujillo
e5 = None
if NC.exists():
    ds = xr.open_dataset(NC)
    ds = ds.rename({k: v for k, v in {"valid_time": "time", "latitude": "lat", "longitude": "lon",
                                       "pressure_level": "level"}.items() if k in ds.dims or k in ds.coords})
    pt = ds.sel(lat=-8.0, lon=-79.0, method="nearest").sel(time=slice(INI, FIN + pd.offsets.MonthEnd(0) + pd.Timedelta(hours=23)))
    e5 = {}
    for t, meses in TEMPORADAS.items():
        s = pt.sel(time=pt.time.dt.month.isin(meses)).mean("time")
        e5[t] = pd.DataFrame({"z_m": (s["z"] / G).values, "T_C": (s["t"] - 273.15).values, "HR": s["r"].values}
                             ).sort_values("z_m")
    E5_PUNTO = (float(pt.lat), float(pt.lon))
    print(f"ERA5 12 UTC: celda {E5_PUNTO}, {pt.time.size} meses")
else:
    print("Sin ERA5 (ejecutar 34a); las figuras se generan solo con los sondeos.")


# ------------------------------------------------------------------ 4. parcela (adiabática seca → húmeda)
def parcela(d):
    """Perfil de la parcela que asciende desde el primer nivel del sondeo."""
    z0, t0, td0 = d.z_m.iloc[0], d.T_C.iloc[0], d.Td_C.iloc[0]
    pz = lambda z: np.exp(np.interp(z, d.z_m, np.log(d.p_hPa.interpolate())))
    zs = np.arange(z0, 5001, 10.0)
    tp, z_nca = [], None
    t = t0
    for i, z in enumerate(zs):
        if i:
            if z_nca is None:
                t -= 9.8e-3 * 10
            else:
                tk = t + 273.15
                rs = EPS * es(t) / (pz(z) - es(t))
                gm = G * (1 + LV * rs / (RD * tk)) / (CP + LV ** 2 * rs * EPS / (RD * tk ** 2))
                t -= gm * 10
        td = td0 - 1.8e-3 * (z - z0)                           # punto de rocío de la parcela (~1.8 °C/km)
        if z_nca is None and t <= td:
            z_nca = z
        tp.append(t)
    return zs, np.array(tp), z_nca


def referencias(ax, xtxt, ha="right"):
    for z, txt in ((Z_CIERRE, "El Tambo (cierre)"), (Z_MED, "cota mediana de la cuenca"), (Z_MAXC, "cota máxima")):
        ax.axhline(z, color="#52514e", lw=0.8, ls=(0, (1, 2)), zorder=1)
        ax.text(xtxt, z + 40, f"{txt} · {z:,.0f} m".replace(",", " "), fontsize=7.5, color="#52514e",
                ha=ha, va="bottom")


plt.rcParams.update({"font.size": 9.5, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.edgecolor": "#898781", "axes.grid": True, "grid.color": "#e1e0d9", "grid.linewidth": 0.6})

# ------------------------------------------------------------------ 5. figura 1: sondeo típico
seca = inv[(inv.temporada == "may–sep (seca)") & inv.inversion].dropna(subset=["base_m", "tope_m", "nube_base_m"])
ref = {c: (seca[c].median(), seca[c].quantile(.75) - seca[c].quantile(.25)) for c in ("base_m", "tope_m", "intensidad_C")}
dist = sum(((seca[c] - m) / max(iqr, 1e-6)) ** 2 for c, (m, iqr) in ref.items())
ft = dist.idxmin()
d = crudos[ft]
d5 = d[d.z_m <= 5000]
r = inv.loc[ft]
zs, tp, z_nca = parcela(d)

fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 7.2), sharey=True, gridspec_kw={"width_ratios": [1.5, 1]})
for ax in (a1, a2):
    ax.axhspan(r.base_m, r.tope_m, color="#e34948", alpha=0.12, lw=0, zorder=0)
    if pd.notna(r.nube_base_m):
        ax.axhspan(r.nube_base_m, r.nube_tope_m, color="#898781", alpha=0.25, lw=0, zorder=0)
z_dalr = np.array([d.z_m.iloc[0], 5000])
a1.plot(d.T_C.iloc[0] - 9.8e-3 * (z_dalr - z_dalr[0]), z_dalr, color="#52514e", ls="--", lw=1.6,
        label="Adiabática seca (9.8 °C/km)")
a1.plot(tp, zs, color="#898781", ls=":", lw=2.2, label="Parcela: seca hasta el NCA y luego húmeda")
a1.plot(d5.T_C, d5.z_m, color="#0b0b0b", lw=2.2, marker="o", ms=3, label="Temperatura medida")
if z_nca:
    a1.plot(np.interp(z_nca, zs, tp), z_nca, "o", color="#52514e", ms=7, zorder=5)
    a1.annotate(f"NCA {z_nca:.0f} m", (np.interp(z_nca, zs, tp), z_nca), xytext=(-70, -4),
                textcoords="offset points", fontsize=8, color="#52514e")
a1.annotate(f"Inversión: +{r.intensidad_C:.1f} °C en {r.espesor_m:.0f} m", (r.T_tope_C, r.tope_m), xytext=(12, 6),
            textcoords="offset points", fontsize=8.5, color="#b8322f", fontweight="bold")
a1.set(xlabel="Temperatura (°C)", ylabel="Altitud (m s.n.m.)", ylim=(0, 5000),
       xlim=(tp[zs <= 5000].min() - 1, max(d5.T_C.max(), d.T_C.iloc[0]) + 4))
referencias(a1, a1.get_xlim()[0] + 0.3, "left")
a2.plot(d5.HR, d5.z_m, color="#0b0b0b", lw=2.2, marker="o", ms=3, label="Humedad relativa medida")
a2.axvline(HR_NUBE, color="#898781", lw=0.8)
a2.set(xlabel="Humedad relativa (%)", xlim=(0, 102))
referencias(a2, 2, "left")
h1, l1 = a1.get_legend_handles_labels()
h1 += [Patch(color="#e34948", alpha=0.25), Patch(color="#898781", alpha=0.4)]
l1 += ["Capa de inversión térmica", f"Capa de nubes (HR ≥ {HR_NUBE:.0f} %)"]
a1.legend(h1, l1, loc="upper right", fontsize=8, frameon=True, framealpha=0.92)
a1.set_title("(a) Temperatura", loc="left", fontsize=10)
a2.set_title("(b) Humedad relativa", loc="left", fontsize=10)
fig.suptitle(f"Sondeo típico de la estación seca · Trujillo (IGRA2 84516), {ft:%Y-%m-%d} 12 UTC (7 a.m.)",
             fontsize=11.5, fontweight="bold", x=0.02, ha="left")
fig.text(0.02, 0.005, f"Elegido entre {len(seca)} sondeos de may–sep con inversión por estar más cerca de las medianas de "
         f"base ({ref['base_m'][0]:.0f} m),\n"
         f"tope ({ref['tope_m'][0]:.0f} m) e intensidad ({ref['intensidad_C'][0]:.1f} °C). "
         f"Las líneas punteadas horizontales marcan altitudes de la cuenca (GLO-30).", fontsize=7.5, color="#52514e")
fig.tight_layout(rect=(0, 0.035, 1, 0.96))
fig.savefig(OUT / "fig_sondeo_tipico_inversion.png", dpi=200)
plt.close(fig)
print(f"Sondeo típico: {ft:%Y-%m-%d} base {r.base_m:.0f} m, tope {r.tope_m:.0f} m, +{r.intensidad_C:.1f} °C, NCA {z_nca}")

# ------------------------------------------------------------------ 6. figura 2: perfiles compuestos y estadística mensual
fig = plt.figure(figsize=(13, 7.4))
gs = fig.add_gridspec(1, 3, width_ratios=[1.15, 1, 1.25], wspace=0.08)
a1 = fig.add_subplot(gs[0])
a2 = fig.add_subplot(gs[1], sharey=a1)
a3 = fig.add_subplot(gs[2], sharey=a1)
for t in ("may–sep (seca)", "dic–abr (húmeda)", "oct–nov"):
    c, col = comp[t], COL_T[t]
    banda = t != "oct–nov"
    for ax, v in ((a1, "T"), (a2, "HR")):
        if banda:
            ax.fill_betweenx(MALLA, c[f"{v}_q1"], c[f"{v}_q3"], color=col, alpha=0.15, lw=0)
        ax.plot(c[f"{v}_med"], MALLA, color=col, lw=2.2 if banda else 1.4,
                label=f"Sondeos {t} (n = {c['n']})" if v == "T" else None)
        if e5 is not None:
            e = e5[t][e5[t].z_m <= 6000]
            ax.plot(e[v + ("_C" if v == "T" else "")], e.z_m, color=col, lw=1.2, ls="--", marker=".", ms=4,
                    label=f"ERA5 12 UTC {t}" if v == "T" else None)
a1.set(xlabel="Temperatura (°C)", ylabel="Altitud (m s.n.m.)", ylim=(0, 5000), xlim=(-2, 26))
a2.set(xlabel="Humedad relativa (%)", xlim=(0, 102))
plt.setp(a2.get_yticklabels(), visible=False)
for ax in (a1, a2):
    referencias(ax, ax.get_xlim()[1] - 0.3)
a1.legend(loc="lower left", fontsize=7.5, framealpha=0.92)
a1.set_title("(a) Temperatura: mediana y rango intercuartílico", loc="left", fontsize=9.5)
a2.set_title("(b) Humedad relativa", loc="left", fontsize=9.5)

xm = np.arange(1, 13)
for i in xm:
    s = res.loc[i]
    if s.n == 0 or pd.isna(s.base_med_m):
        continue
    a3.bar(i, s.tope_med_m - s.base_med_m, bottom=s.base_med_m, width=0.62, color="#e34948", alpha=0.35, lw=0)
    a3.vlines(i, s.base_q1_m, s.base_q3_m, color="#0b0b0b", lw=1.2)
    a3.vlines(i + 0.12, s.tope_q1_m, s.tope_q3_m, color="#b8322f", lw=1.2)
    a3.text(i, s.tope_q3_m + 60, f"{s.frecuencia_pct:.0f} %\n(n={int(s.n)})", ha="center", va="bottom", fontsize=7,
            color="#52514e")
a3.axhline(Z_CIERRE, color="#52514e", lw=0.8, ls=(0, (1, 2)))
a3.set(xticks=xm, xticklabels=MESES, xlim=(0.4, 12.6), xlabel="Mes")
plt.setp(a3.get_yticklabels(), visible=False)
a3.set_title("(c) Capa de inversión por mes (base → tope)", loc="left", fontsize=9.5)
a3.legend([Patch(color="#e34948", alpha=0.35), Line2D([], [], color="#0b0b0b"), Line2D([], [], color="#b8322f")],
          ["Base a tope (medianas)", "RIC de la base", "RIC del tope"], loc="upper right", fontsize=7.5)
fig.suptitle(f"Perfiles verticales sobre la costa de Trujillo y capa de inversión · sondeos IGRA2 "
             f"{inv.index.min():%Y-%m} a {inv.index.max():%Y-%m} (12 UTC)", fontsize=11.5, fontweight="bold",
             x=0.02, ha="left")
nota = (f"Inversión = primer aumento continuo de T con la altura desde la superficie, por debajo de {Z_MAX_INV:.0f} m y "
        f"≥ {UMBRAL_INV:.0f} °C. En (c), el porcentaje es la frecuencia de sondeos con inversión.\n")
if e5 is not None:
    nota += f"ERA5: promedio mensual de las 12 UTC {INI:%Y-%m} a {FIN:%Y-%m}, celda {E5_PUNTO[0]:.2f}°, {E5_PUNTO[1]:.2f}°."
fig.text(0.02, 0.008, nota, fontsize=7.5, color="#52514e")
fig.subplots_adjust(left=0.06, right=0.99, top=0.9, bottom=0.125)
fig.savefig(OUT / "fig_perfiles_inversion.png", dpi=200)
plt.close(fig)
print("Figuras en", OUT)
