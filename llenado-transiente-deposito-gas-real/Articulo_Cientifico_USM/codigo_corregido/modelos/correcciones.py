"""Experimentos adicionales de control metodológico (no reemplazan la corrida histórica).

1. V1 vs V2 a igual P_obj, EOS, suministro, área y condición inicial.
2. V3: sensibilidad de masa térmica efectiva, sin afirmarla como masa real.
3. Diferenciar verificación computacional de validación física experimental.
"""
from __future__ import annotations
from dataclasses import replace
from pathlib import Path
import csv

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from .base import BASE, relacion_critica
from .peng_robinson import MODELS, sim, resumir, BAR


def ejecutar(datos: Path, graficos: Path):
    datos.mkdir(parents=True, exist_ok=True)
    graficos.mkdir(parents=True, exist_ok=True)
    filas, paredes = [], []
    for objetivo in (100., 190.):
        p = replace(BASE, P_OBJ=objetivo*BAR)
        for modelo in MODELS:
            v1 = resumir(sim('V1', modelo, p))
            v2 = resumir(sim('V2', modelo, p))
            filas.append(dict(eos=modelo, Pobj_bar=objetivo,
                t_V1_s=v1['tf_s'], t_V2_s=v2['tf_s'],
                diferencia_s=v2['tf_s']-v1['tf_s'],
                diferencia_pct=100*(v2['tf_s']-v1['tf_s'])/v1['tf_s'],
                temperatura_V1_K=v1['Tf_K'], temperatura_V2_K=v2['Tf_K'],
                masa_V1_kg=v1['mf_kg'], masa_V2_kg=v2['mf_kg']))
    for masa in (8., 20., 54.):
        p = replace(BASE, P_OBJ=190*BAR, MW=masa)
        for modelo in MODELS:
            r = resumir(sim('V3', modelo, p))
            paredes.append(dict(eos=modelo, masa_efectiva_kg=masa,
                C_pared_J_K=p.CW, t_final_s=r['tf_s'],
                temperatura_final_K=r['Tf_K'], temperatura_pared_final_K=r['Twf_K'],
                masa_final_kg=r['mf_kg'], energia_gas_pared_kJ=r['Qgw_kJ']))

    def guardar(nombre, rows):
        with (datos / nombre).open('w', encoding='utf-8-sig', newline='') as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            w.writeheader(); w.writerows(rows)
    guardar('comparacion_controlada_V1_V2.csv', filas)
    guardar('sensibilidad_masa_pared.csv', paredes)
    for rows, x, y, fname, xlabel, ylabel in [
        (filas,'Pobj_bar','diferencia_s','comparacion_V1_V2.png','Presión objetivo [bar]','Tiempo V2 - V1 [s]'),
        (paredes,'masa_efectiva_kg','temperatura_final_K','sensibilidad_pared_T.png','Masa térmica efectiva [kg]','Temperatura final gas [K]'),
        (paredes,'masa_efectiva_kg','t_final_s','sensibilidad_pared_t.png','Masa térmica efectiva [kg]','Tiempo llenado [s]')]:
        fig, ax = plt.subplots(figsize=(9,5.1))
        for modelo in MODELS:
            r = sorted((z for z in rows if z['eos']==modelo), key=lambda z:z[x])
            ax.plot([z[x] for z in r],[z[y] for z in r],marker='o',label=modelo,linewidth=2)
        ax.set(xlabel=xlabel,ylabel=ylabel)
        ax.grid(alpha=.25);ax.legend();fig.tight_layout();fig.savefig(graficos/fname,dpi=200);plt.close(fig)
    return {'comparacion':filas,'pared':paredes,'presion_critica_bar':BASE.PS*relacion_critica(BASE.GAMMA)/BAR}
