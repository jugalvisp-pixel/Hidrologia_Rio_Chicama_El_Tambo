"""
Estación seca frente a estación húmeda: de dónde llega la humedad y qué la frena sobre la cuenca El Tambo.

Una fila por temporada (may–sep seca arriba, dic–abr húmeda abajo):
  - Izquierda: transporte de vapor integrado en la vertical (IVT) de ERA5, climatología del periodo común
    (color = magnitud, flechas = dirección). Entrada: script 24.
  - Derecha: perfil vertical de temperatura sobre la costa de Trujillo (radiosondeos IGRA2 84516, 12 UTC):
    mediana y rango intercuartílico de la temporada, adiabática seca de referencia desde la superficie y capa de
    inversión (medianas de base y tope de los sondeos con inversión), con las altitudes de la cuenca.
    Entrada: script 34b.
Mismas escalas en las dos filas para comparar directamente.

Salida: Tarea 1/El_Tambo/resultados/inversion/fig_temporadas_humedad_inversion.png
"""
from pathlib import Path

import numpy as np
import pandas as pd
import geopandas as gpd
import xarray as xr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

from periodo_comun import periodo_comun
from estilo_figuras import EST, CMAP_P

ROOT = Path(__file__).resolve().parents[2] / "Tarea 1" / "El_Tambo"
CLIM = ROOT / "datos" / "era5" / "era5_viento_ivt_climatologia_estacional.nc"
INV = ROOT / "resultados" / "inversion"
DL = ROOT / "resultados" / "delimitacion"
INI, FIN = periodo_comun()

TEMP = [("may–sep (seca)", "may-sep", "Estación seca (may–sep)"),
        ("dic–abr (húmeda)", "dic-abr", "Estación húmeda (dic–abr)")]
COL = {"may–sep (seca)": EST["seca"], "dic–abr (húmeda)": EST["húmeda"]}
LON, LAT = (-100, -30), (-35, 15)
PASO = 3

clim = xr.open_dataset(CLIM)
paises = gpd.read_file(ROOT / "datos" / "ne_50m_admin_0_countries.geojson")
cuenca = gpd.read_file(DL / "cuenca_el_tambo_glo30.gpkg").to_crs(4326)
cen = cuenca.geometry.iloc[0].centroid
comp = pd.read_csv(INV / "perfiles_compuestos_temporada.csv")
sond = pd.read_csv(INV / "inversion_por_sondeo.csv", parse_dates=["fecha"])
morf = pd.read_csv(DL / "parametros_morfometricos.csv", index_col=0)["valor"]
Z_REF = [(float(morf["Cota cierre (m)"]), "El Tambo (cierre)"),
         (float(morf["Cota mediana (m s.n.m.)"]), "cota mediana de la cuenca"),
         (float(morf["Cota máxima (m s.n.m.)"]), "cota máxima")]
VMAX = float(clim["ivt"].sel(estacion=["dic-abr", "may-sep"]).quantile(0.99))

plt.rcParams.update({"font.size": 9.5})
fig = plt.figure(figsize=(13, 11))
gs = fig.add_gridspec(2, 2, width_ratios=[1.55, 1], hspace=0.28, wspace=0.26)
for fila, (t, est, titulo) in enumerate(TEMP):
    # --- transporte de humedad
    ax = fig.add_subplot(gs[fila, 0])
    d = clim.sel(estacion=est, lon=slice(*LON), lat=slice(*LAT))
    m = ax.pcolormesh(d.lon, d.lat, d["ivt"], cmap=CMAP_P, vmin=0, vmax=VMAX, shading="nearest")
    paises.boundary.plot(ax=ax, color="#333333", lw=0.4)
    s = slice(None, None, PASO)
    q = ax.quiver(d.lon[s], d.lat[s], d["ivt_e"][s, s], d["ivt_n"][s, s], scale=2500, width=0.0028, color="#1d1d1b")
    ax.add_patch(matplotlib.patches.Rectangle((0.015, 0.015), 0.36, 0.08, transform=ax.transAxes, fc="white",
                                              ec="#999999", lw=0.5, alpha=0.92, zorder=7))
    ax.quiverkey(q, 0.15, 0.055, 250, "250 kg m⁻¹ s⁻¹", labelpos="E", coordinates="axes",
                 fontproperties={"size": 8}, zorder=8)
    cuenca.plot(ax=ax, color="red", zorder=5)
    ax.plot(cen.x, cen.y, "o", ms=10, mfc="none", mec="red", mew=1.8, zorder=6)
    ax.set(xlim=LON, ylim=LAT, xlabel="", ylabel="")
    ax.set_xticks(range(-100, -29, 10))
    ax.set_xticklabels([f"{-x}° O" for x in range(-100, -29, 10)], fontsize=8)
    ax.set_yticks(range(-30, 16, 10))
    ax.set_yticklabels([f"{abs(y)}° {'S' if y < 0 else 'N' if y > 0 else ''}" for y in range(-30, 16, 10)], fontsize=8)
    ax.set_aspect("equal")
    ax.set_title(f"({'a' if fila == 0 else 'c'}) {titulo} · transporte de humedad (IVT, ERA5)", loc="left",
                 fontsize=10.5, fontweight="bold")
    cb = fig.colorbar(m, ax=ax, fraction=0.035, pad=0.015)
    cb.set_label("IVT (kg m⁻¹ s⁻¹)", fontsize=8.5)
    cb.ax.tick_params(labelsize=7.5)

    # --- perfil de temperatura
    a2 = fig.add_subplot(gs[fila, 1])
    z = comp["z_m"].values
    med, q1, q3 = comp[f"{t}|T_med"].values, comp[f"{t}|T_q1"].values, comp[f"{t}|T_q3"].values
    sv = sond[(sond.temporada == t)]
    si = sv[sv.inversion]
    base, tope = si.base_m.median(), si.tope_m.median()
    a2.axhspan(base, tope, color="#e34948", alpha=0.15, lw=0, zorder=0)
    a2.fill_betweenx(z, q1, q3, color=COL[t], alpha=0.18, lw=0)
    a2.plot(med, z, color=COL[t], lw=2.4)
    ok = np.isfinite(med)
    z0, t0 = z[ok][0], med[ok][0]
    zz = np.array([z0, 5000])
    a2.plot(t0 - 9.8e-3 * (zz - z0), zz, color="#52514e", lw=1.2, ls="--")
    for zr, txt in Z_REF:
        a2.axhline(zr, color="#52514e", lw=0.8, ls=(0, (1, 2)))
        a2.text(-1.6, zr + 40, f"{txt} · {zr:,.0f} m".replace(",", " "), fontsize=7.5, color="#52514e")
    a2.text(0.97, 0.70, f"inversión +{si.intensidad_C.median():.1f} °C\n{base:.0f}–{tope:.0f} m · "
            f"{100 * sv.inversion.mean():.0f} % de los días", ha="right", va="center", fontsize=8.5,
            color="#b8322f", fontweight="bold", transform=a2.transAxes,
            bbox=dict(fc="white", ec="#e34948", lw=0.6, alpha=0.9))
    a2.set(xlim=(-2, 26), ylim=(0, 5000), xlabel="Temperatura (°C)", ylabel="Altitud (m s.n.m.)")
    a2.grid(color="#e1e0d9", lw=0.6)
    for sp in ("top", "right"):
        a2.spines[sp].set_visible(False)
    a2.set_title(f"({'b' if fila == 0 else 'd'}) Perfil de temperatura (n = {len(sv)} sondeos)",
                 loc="left", fontsize=10.5, fontweight="bold")
    if fila == 0:
        a2.legend([Line2D([], [], color="#0b0b0b", lw=2.4), Patch(color="#898781", alpha=0.3),
                   Line2D([], [], color="#52514e", lw=1.2, ls="--"), Patch(color="#e34948", alpha=0.3)],
                  ["Mediana de los sondeos (7 a.m.)", "Rango intercuartílico", "Adiabática seca (9.8 °C/km)",
                   "Capa de inversión (base–tope)"], loc="upper right", fontsize=7.8, framealpha=0.92)

fig.suptitle("De dónde llega la humedad y qué la frena: estación seca frente a estación húmeda",
             fontsize=13, fontweight="bold", x=0.02, ha="left", y=0.975)
fig.text(0.02, 0.012,
         f"IVT: flujo de vapor de agua integrado en la vertical, ERA5 mensual {INI:%Y-%m} a {FIN:%Y-%m} (malla de 1°), "
         f"una flecha cada {PASO}°; contorno rojo: cuenca El Tambo. Perfiles: radiosondeos IGRA2 84516 (Trujillo, a ~40 km "
         f"de la salida de la cuenca),\n{sond.fecha.min():%Y-%m} a {sond.fecha.max():%Y-%m}, 12 UTC; inversión = primer "
         "aumento continuo de T con la altura desde la superficie (≥ 1 °C, < 3 500 m). Scripts 24, 34a, 34b y 35.",
         fontsize=7.5, color="#52514e")
fig.subplots_adjust(left=0.05, right=0.985, top=0.92, bottom=0.07)
OUT = INV / "fig_temporadas_humedad_inversion.png"
fig.savefig(OUT, dpi=200)
print("Figura:", OUT)
