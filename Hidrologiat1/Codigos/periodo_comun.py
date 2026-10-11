"""
Periodo común de análisis: meses en que hay a la vez precipitación de CAMELS-PE (P_L) e IMERG (P_I) válidas,
recortado a AÑOS CALENDARIO COMPLETOS (enero a diciembre). Todas las figuras, mapas y diagramas lo usan.

Decisiones del grupo: periodo común P_L–P_I (2026-10-08); solo años climáticos completos (2026-10-11), para
que cada mes calendario tenga el mismo número de años y los totales anuales no queden truncados. Con los
datos actuales: 1998-01 a 2024-12 (27 años, 324 meses); los meses 2025-01..2025-08 de IMERG quedan fuera.

Se calcula a partir de los datos (no se fija a mano), de modo que si cambia la serie de IMERG el periodo se
actualiza solo. Requiere las salidas de 09 (CAMELS-PE) y 10 (IMERG).
"""
from pathlib import Path

import pandas as pd

SM = Path(__file__).resolve().parents[2] / "Tarea 1" / "El_Tambo" / "resultados" / "series_mensuales"


def periodo_comun():
    """Devuelve (inicio, fin) como primer día del primer y del último mes común (años completos)."""
    pl = pd.read_csv(SM / "series_mensuales_camels_pe.csv", index_col="fecha", parse_dates=True)["P_L_mm_mes"]
    pi = pd.read_csv(SM / "imerg_cuenca_mensual.csv", index_col="fecha", parse_dates=True)["P_I_mm_mes"]
    c = pd.concat([pl, pi], axis=1).dropna()
    if len(c) != len(pd.date_range(c.index.min(), c.index.max(), freq="MS")):
        raise ValueError("El periodo común P_L–P_I tiene meses intermedios sin dato; revisar las series.")
    ini, fin = c.index.min(), c.index.max()
    if ini.month != 1:                                    # primer enero disponible
        ini = pd.Timestamp(f"{ini.year + 1}-01-01")
    if fin.month != 12:                                   # último diciembre disponible
        fin = pd.Timestamp(f"{fin.year - 1}-12-01")
    return ini, fin


def fin_diario(fin):
    """Último día del mes final, para recortar series diarias."""
    return fin + pd.offsets.MonthEnd(0)


if __name__ == "__main__":
    i, f = periodo_comun()
    print(f"Periodo común P_L–P_I (años completos): {i:%Y-%m} a {f:%Y-%m} "
          f"({len(pd.date_range(i, f, freq='MS'))} meses)")
