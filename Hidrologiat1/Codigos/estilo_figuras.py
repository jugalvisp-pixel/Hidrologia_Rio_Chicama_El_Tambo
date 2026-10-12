"""Paleta y escalas comunes a todas las figuras (PNG y página web).

Una variable tiene el mismo color en todo el trabajo y ese color no se usa para otra cosa.
Los colores de variable salen de una paleta categórica validada para daltonismo (los tres primeros
se distinguen en cualquier combinación); métodos, particiones y bandas usan tinta, grises o rampas
ordenadas, nunca un color de variable.

Uso:  from estilo_figuras import VAR, ETQ, EST, CMAP_P, ...
"""
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap, ListedColormap

# ---------------------------------------------------------------- tinta y fondo
TINTA = "#0b0b0b"        # texto, observado, líneas de referencia
TINTA2 = "#52514e"       # texto secundario, segundo método
GRIS = "#898781"         # referencias, apoyo, elementos atenuados
GRIS_CLARO = "#c3c2b7"
GRILLA = "#e1e0d9"
SIN_DATO = "#d9d9d9"

# ---------------------------------------------------------------- variables (identidad fija)
VAR = {
    "P_L": "#2a78d6",     # precipitación de referencia, CAMELS-PE (PISCOp v2.1)
    "P_I": "#eb6834",     # IMERG V07 Final
    "Q": "#1baf7a",       # caudal medio (m³/s)
    "R": "#1baf7a",       # el mismo caudal en lámina (mm/mes): mismo color que Q
    "T": "#c0392b",       # temperatura PISCOt (CAMELS-PE)
    "T_ERA5L": "#7d1f16",  # temperatura ERA5-Land: mismo tono, más oscuro (y línea discontinua)
    "P_v3": "#52514e",    # PISCOp v3.0 (apoyo)
    "P_CHIRPS": "#898781",  # CHIRPS v2.0 (apoyo)
    "Q_obs": "#0e6e4b",   # caudal solo con meses observados completos (sensibilidad)
    "PET": "#3d3d3d",
}
# nombres legibles (matplotlib); el texto nunca va en el color de la serie
ETQ = {
    "P_L": "P$_L$ (CAMELS-PE)",
    "P_I": "P$_I$ (IMERG V07)",
    "Q": "$\\bar{Q}$ (El Tambo)",
    "R": "R ($\\bar{Q}$ en lámina)",
    "T": "T (PISCOt)",
    "T_ERA5L": "T (ERA5-Land)",
    "P_v3": "P (PISCOp v3.0)",
    "P_CHIRPS": "P (CHIRPS v2.0)",
    "Q_obs": "$\\bar{Q}$ solo observado",
}
# estilo de línea por fuente dentro de una misma variable
LS = {"T": "-", "T_ERA5L": "--", "Q": "-", "Q_obs": "--"}

# ---------------------------------------------------------------- estaciones de la cuenca
# Definición única en todo el trabajo: húmeda dic–abr, transición oct–nov, seca may–sep.
ESTACION_MES = {12: "húmeda", 1: "húmeda", 2: "húmeda", 3: "húmeda", 4: "húmeda",
                5: "seca", 6: "seca", 7: "seca", 8: "seca", 9: "seca", 10: "transición", 11: "transición"}
ESTACIONES = ["húmeda", "transición", "seca"]
EST_ETQ = {"húmeda": "húmeda (dic–abr)", "transición": "transición (oct–nov)", "seca": "seca (may–sep)"}
EST = {"húmeda": "#6250d6", "transición": "#898781", "seca": "#c98500"}   # violeta / gris / ocre
EST_MK = {"húmeda": "o", "transición": "s", "seca": "^"}

# ---------------------------------------------------------------- métodos, particiones
METODO = {"OLS": TINTA, "Sen": TINTA2, "LOWESS": TINTA, "ref": GRIS}
METODO_LS = {"OLS": "--", "Sen": ":", "LOWESS": "-"}
METODO_MK = {"OLS": "s", "Sen": "o"}
PARTICION = {"A": TINTA, "B": GRIS}

# ---------------------------------------------------------------- escalas continuas
# Magnitud de agua (lluvia, IVT, pesos): una sola rampa azul, de claro (≈0) a oscuro.
_AZUL = ["#f4f8fd", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b", "#081f40"]
CMAP_P = LinearSegmentedColormap.from_list("agua", _AZUL)
CMAP_P.set_bad(SIN_DATO)
# Clases comunes de precipitación mensual (mm/mes) y anual (mm/año)
CLASES_P_MES = [0, 5, 10, 20, 30, 40, 60, 80, 100, 125, 150, 175, 200, 250]
CLASES_P_ANUAL = list(range(300, 1001, 50))
# Elevación: rampa hipsométrica sin azul (el azul se reserva al agua)
CMAP_ELEV = LinearSegmentedColormap.from_list(
    "elevacion", ["#5d8a4f", "#a7c47f", "#e8e2a8", "#d9b77c", "#b08457", "#8a6a55", "#c9c2bb", "#ffffff"])
# Otra magnitud no hídrica en un mapa (desnivel, etc.): rampa gris-marrón
CMAP_NEUTRA = LinearSegmentedColormap.from_list("neutra", ["#f7f5f2", "#d6cfc6", "#a8998a", "#6f5f52", "#3b3029"])
# Divergentes con centro neutro: agua (seco marrón ↔ húmedo verde azulado), temperatura y correlación (azul ↔ rojo)
CMAP_DIV_AGUA = mpl.colormaps["BrBG"]
CMAP_DIV_T = mpl.colormaps["RdBu_r"]
CMAP_CORR = mpl.colormaps["RdBu_r"]          # punto 5: escala común −1 a 1
# Rampa ordinal (bandas de periodo, de corto a largo): luminosidad monótona, apta para daltonismo
CMAP_ORDINAL = mpl.colormaps["viridis"]


def colores_ordinales(n, lo=0.08, hi=0.92):
    """n colores ordenados de la rampa ordinal."""
    return [CMAP_ORDINAL(lo + (hi - lo) * i / max(n - 1, 1)) for i in range(n)]


# Disponibilidad / control de calidad (estado: siempre con leyenda de texto)
DISP = {"completo": "#256abf", "escalado": "#9ec5f4", "excluido": "#d03b3b", "sin_dato": SIN_DATO, "fuera": "#ffffff"}
CMAP_DISP = ListedColormap([DISP[k] for k in ("completo", "escalado", "excluido", "sin_dato", "fuera")])
# Meses de caudal rellenados (marcadores en tinta, no en un color de variable)
RELLENO = {"interpolado_diario": dict(marker="o", mfc="none", mec=TINTA),
           "escalado": dict(marker="s", mfc="none", mec=TINTA),
           "modelo_mensual": dict(marker="D", mfc=TINTA, mec=TINTA)}
SOMBRA_EXCLUIDO = "#d03b3b"   # alpha bajo
SOMBRA_ESTIMADO = "#898781"   # tramo estimado de T (alpha bajo; neutro para no confundirse con excluidos)
INVERSION = "#d03b3b"         # capa de inversión térmica (alpha bajo)


def aplicar_estilo():
    """Ejes y grilla discretos, comunes a todas las figuras."""
    mpl.rcParams.update({
        "axes.edgecolor": "#3d3d3d", "axes.labelcolor": TINTA, "xtick.color": "#3d3d3d", "ytick.color": "#3d3d3d",
        "grid.color": GRILLA, "grid.linewidth": 0.7, "axes.titlecolor": TINTA, "legend.frameon": False,
    })
