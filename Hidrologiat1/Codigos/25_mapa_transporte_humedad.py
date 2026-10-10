"""
Mapas regionales del transporte de humedad y del viento (ERA5, climatología del periodo común) para
contrastar el esquema de fuentes de humedad tomado de la bibliografía (recuadro del mapa del script 06).

Paneles (estación húmeda dic–abr y seca may–sep):
  (a, b) IVT: flujo de vapor de agua integrado en la vertical. Color = magnitud (kg m⁻¹ s⁻¹),
         flechas = dirección hacia donde se mueve la humedad.
  (c, d) Viento en un nivel de presión. Color = rapidez (m/s), flechas = dirección hacia donde sopla.
         Se omiten las celdas donde el nivel queda bajo la superficie del modelo (presión en superficie <
         nivel), porque ahí ERA5 entrega valores extrapolados.
El IVT es la variable que describe el transporte de humedad; el viento de un solo nivel muestra la
circulación a esa altura. La convergencia del flujo (no su magnitud) es la que se relaciona con la lluvia.

Entrada : Tarea 1/El_Tambo/datos/era5/era5_viento_ivt_climatologia_estacional.nc (script 24)
Salida  : Tarea 1/El_Tambo/resultados/contexto_fisico/mapa_transporte_humedad_era5.png
"""
from pathlib import Path

import numpy as np
import xarray as xr
import geopandas as gpd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2] / "Tarea 1" / "El_Tambo"
CLIM = ROOT / "datos" / "era5" / "era5_viento_ivt_climatologia_estacional.nc"
OUT = ROOT / "resultados" / "contexto_fisico"
NIVEL = 500
EXT = (-100, -30, -35, 15)                                 # lon_min, lon_max, lat_min, lat_max
PASO = 3                                                    # una flecha cada 3° (malla de 1°)

clim = xr.open_dataset(CLIM).sel(lon=slice(EXT[0], EXT[1]), lat=slice(EXT[2], EXT[3]))
paises = gpd.read_file(ROOT / "datos" / "ne_50m_admin_0_countries.geojson")
cuenca = gpd.read_file(ROOT / "resultados" / "delimitacion" / "cuenca_el_tambo_glo30_wgs84.geojson")
c = cuenca.geometry.iloc[0].centroid

fig, axs = plt.subplots(2, 2, figsize=(13, 9.6), constrained_layout=True)
paneles = [("ivt", "dic-abr", "(a) IVT – estación húmeda (dic–abr)"),
           ("ivt", "may-sep", "(b) IVT – estación seca (may–sep)"),
           ("vel", "dic-abr", f"(c) Viento en {NIVEL} hPa – estación húmeda (dic–abr)"),
           ("vel", "may-sep", f"(d) Viento en {NIVEL} hPa – estación seca (may–sep)")]
lims = {"ivt": (0, float(clim["ivt"].sel(estacion=["dic-abr", "may-sep"]).quantile(0.99))),
        "vel": (0, float(clim["vel"].sel(nivel=NIVEL, estacion=["dic-abr", "may-sep"]).quantile(0.99)))}
for ax, (var, est, titulo) in zip(axs.flat, paneles):
    d = clim.sel(estacion=est)
    if var == "ivt":
        campo, u, v = d["ivt"], d["ivt_e"], d["ivt_n"]
        cmap, etq, esc = "YlGnBu", "IVT (kg m⁻¹ s⁻¹)", 2500
    else:
        d = d.sel(nivel=NIVEL)
        campo, u, v = d["vel"], d["u"], d["v"]
        cmap, etq, esc = "PuRd", f"Rapidez del viento en {NIVEL} hPa (m/s)", 180
    m = ax.pcolormesh(d.lon, d.lat, campo, cmap=cmap, vmin=lims[var][0], vmax=lims[var][1], shading="nearest")
    if var == "vel":
        bajo = d["bajo_terreno"].where(d["bajo_terreno"] == 1)
        ax.pcolormesh(d.lon, d.lat, bajo, cmap=matplotlib.colors.ListedColormap(["#bdbdbd"]), shading="nearest")
    paises.boundary.plot(ax=ax, color="#333333", lw=0.4)
    s = slice(None, None, PASO)
    q = ax.quiver(d.lon[s], d.lat[s], u[s, s], v[s, s], scale=esc, width=0.0025, color="#1d1d1b")
    ref = 250 if var == "ivt" else 10
    ax.set_xlabel(""); ax.set_ylabel("")
    ax.quiverkey(q, 0.74, 1.03, ref, f"{ref} {'kg m⁻¹ s⁻¹' if var == 'ivt' else 'm/s'}", labelpos="E",
                 coordinates="axes", fontproperties={"size": 8})
    cuenca.plot(ax=ax, color="red", edgecolor="red", lw=1.2, zorder=5)
    ax.plot(c.x, c.y, "o", ms=9, mfc="none", mec="red", mew=1.6, zorder=6)
    ax.set_xlim(EXT[0], EXT[1]); ax.set_ylim(EXT[2], EXT[3]); ax.set_aspect("equal")
    ax.set_title(titulo, fontsize=10.5, loc="left")
    ax.set_xticks(range(-100, -29, 10), [f"{abs(x)}° O" for x in range(-100, -29, 10)], fontsize=8)
    ax.set_yticks(range(-30, 16, 10), [f"{abs(y)}° {'S' if y < 0 else 'N' if y > 0 else ''}"
                                       for y in range(-30, 16, 10)], fontsize=8)
    cb = fig.colorbar(m, ax=ax, shrink=0.85, pad=0.01)
    cb.set_label(etq, fontsize=8.5); cb.ax.tick_params(labelsize=8)
fig.suptitle(f"Transporte de humedad y viento – ERA5, climatología {clim.attrs['periodo']}", fontsize=12.5)
fig.text(0.5, -0.025,
         "Flechas: dirección hacia donde se mueve el flujo (una cada 3°); longitud proporcional a la magnitud. "
         "Gris en (c, d): el nivel queda bajo la superficie del modelo (valores extrapolados, omitidos).\n"
         "Contorno rojo: cuenca El Tambo. Malla de 1° (promedio de 0.5°). Fuente: ERA5 monthly averaged reanalysis "
         "(C3S/ECMWF; doi:10.24381/cds.6860a573 y 10.24381/cds.f17050d7). Script 25_mapa_transporte_humedad.py.",
         ha="center", va="top", fontsize=8, color="#444")
fig.savefig(OUT / "mapa_transporte_humedad_era5.png", dpi=200, bbox_inches="tight")
plt.close(fig)

# Valores cerca de la cuenca (celda de 1° que la contiene y celda al este, sobre los Andes orientales)
for lon in (-78.5, -76.5):
    p = clim.sel(lat=-7.5, lon=lon, method="nearest")
    for est in ("dic-abr", "may-sep"):
        e = p.sel(estacion=est)
        dir_ivt = np.degrees(np.arctan2(float(e.ivt_e), float(e.ivt_n))) % 360
        print(f"lon {float(p.lon):6.1f} {est}: IVT {float(e.ivt):5.0f} hacia {dir_ivt:3.0f}° | "
              f"tcwv {float(e.tcwv):4.1f} kg/m² | sp {float(e.sp):4.0f} hPa")
