"""Utilidades cartográficas comunes para mapas en coordenadas geográficas (EPSG:4326)."""
import numpy as np
from matplotlib.ticker import FuncFormatter, MultipleLocator
from matplotlib_scalebar.scalebar import ScaleBar

M_PER_DEG_Y = 110574.0


def m_por_grado_lon(lat_deg):
    """Metros por grado de longitud a una latitud dada (esfera aproximada)."""
    return 111320.0 * np.cos(np.radians(lat_deg))


def formato_grados(ax, lat0, paso=0.1):
    """Cuadrícula y rótulos en grados con 2 decimales; aspecto corregido por cos(lat)."""
    ax.set_aspect(1 / np.cos(np.radians(lat0)))
    ax.xaxis.set_major_locator(MultipleLocator(paso))
    ax.yaxis.set_major_locator(MultipleLocator(paso))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{abs(v):.2f}° {'O' if v < 0 else 'E'}"))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{abs(v):.2f}° {'S' if v < 0 else 'N'}"))
    ax.set_xlabel("Longitud"); ax.set_ylabel("Latitud")
    ax.tick_params(axis="y", labelrotation=90)
    ax.grid(color="grey", alpha=0.35, linewidth=0.5, linestyle="--")


def escala_y_norte(ax, lat0, fontsize=12):
    """Barra de escala (metros por grado de longitud a lat0) y flecha de norte."""
    ax.add_artist(ScaleBar(m_por_grado_lon(lat0), units="m", location="lower right",
                           box_alpha=0.8, rotation="horizontal-only"))
    ax.annotate("N", xy=(0.93, 0.96), xytext=(0.93, 0.87), xycoords="axes fraction",
                ha="center", va="center", fontsize=fontsize, fontweight="bold",
                arrowprops=dict(facecolor="black", width=4, headwidth=11))
