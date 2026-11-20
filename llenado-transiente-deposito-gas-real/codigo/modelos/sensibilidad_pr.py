"""Sensibilidad de V3 a los parámetros de pared con las tres ecuaciones de estado.

La entrega anterior varió UA_gw, UA_wa y Cw solo con ideal y Van der Waals.
Aquí se repiten exactamente los mismos valores agregando Peng-Robinson,
con el mismo objetivo de 190 bar, la misma conexión (Cd, A) y el mismo suministro.
"""
from dataclasses import replace
from pathlib import Path

import numpy as np

from .base import BASE, csv_file, fig_style, export_fig, plt
from .peng_robinson import MODELS, COLORS, BAR, sim, resumir

ESTUDIOS = (
    ("UA_gas_pared", "Conductancia gas-pared [W/K]", "UA_GW", (0.0, 4.0, 8.0, 16.0)),
    ("UA_pared_amb", "Conductancia pared-ambiente [W/K]", "UA_WA", (1.5, 3.0, 6.0)),
    ("C_pared", "Capacidad térmica de pared [J/K]", "CW", (2000.0, 4000.0, 8000.0)),
)


def parametros_estudio(var, valor):
    p = replace(BASE, P_OBJ=190 * BAR)
    if var == "CW":
        return replace(p, MW=valor / p.CPW)
    return replace(p, **{var: valor})


def ejecutar(datos: Path, figuras: Path) -> list[dict]:
    filas = []
    for tag, xlabel, var, valores in ESTUDIOS:
        for valor in valores:
            p = parametros_estudio(var, valor)
            if p.AREA != BASE.AREA or p.CD != BASE.CD:
                raise RuntimeError("El área y Cd deben conservarse en los barridos.")
            for modelo in MODELS:
                r = resumir(sim("V3", modelo, p))
                filas.append({"parametro": tag, "valor": valor, **r})
        for y, ylabel, slug in (("tf_s", "Tiempo de llenado [s]", "tiempo"),
                                ("Tmax_K", "Temperatura máxima del gas [K]", "Tmax"),
                                ("mf_kg", "Masa final [kg]", "masa"),
                                ("Qgw_kJ", "Calor total gas-pared [kJ]", "calor")):
            fig, ax = plt.subplots(figsize=(8, 4.55))
            for modelo in MODELS:
                rr = [f for f in filas if f["parametro"] == tag and f["modelo"] == modelo]
                ax.plot([f["valor"] for f in rr], [f[y] for f in rr], "o-", lw=2,
                        color=COLORS[modelo], label=modelo)
            fig_style(ax, xlabel, ylabel, f"V3 a 190 bar: sensibilidad a {xlabel.split(' [')[0].lower()}")
            export_fig(fig, figuras / f"sens3_{tag}_{slug}.png")
    csv_file(datos / "sensibilidad_pared_3EOS.csv", filas)
    return filas
