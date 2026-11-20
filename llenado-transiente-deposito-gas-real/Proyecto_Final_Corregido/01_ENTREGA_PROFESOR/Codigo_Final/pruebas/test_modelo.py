"""Pruebas automáticas del modelo. Ejecutar desde Codigo_Final con:  python -m pytest -q"""
import json
import math
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from modelos import base as B
from modelos import peng_robinson as PRM

BAR = 1e5
P190 = replace(B.BASE, P_OBJ=190 * BAR)
P100 = replace(B.BASE, P_OBJ=100 * BAR)


def test_configuracion_coincide():
    cfg = json.loads((Path(__file__).resolve().parents[1] / "configuracion" / "caso_base.json").read_text(encoding="utf-8"))
    p = B.BASE
    pares = [(cfg["deposito"]["V_m3"], p.V), (cfg["deposito"]["P0_Pa"], p.P0), (cfg["deposito"]["T0_K"], p.T0),
             (cfg["suministro"]["PS_Pa"], p.PS), (cfg["suministro"]["TS_K"], p.TS),
             (cfg["caudal"]["M_DOT_V1_kg_s"], p.M_DOT_V1), (cfg["caudal"]["CD"], p.CD),
             (cfg["caudal"]["AREA_m2_calibrada"], p.AREA), (cfg["pared"]["MW_kg"], p.MW),
             (cfg["pared"]["CPW_J_kgK"], p.CPW), (cfg["pared"]["UA_GW_W_K"], p.UA_GW),
             (cfg["pared"]["UA_WA_W_K"], p.UA_WA), (cfg["nitrogeno"]["cp_J_kgK"], p.CP),
             (cfg["nitrogeno"]["Tc_K"], B.TC_N2), (cfg["nitrogeno"]["Pc_Pa"], B.PC_N2),
             (cfg["nitrogeno"]["omega"], PRM.OMEGA), (cfg["numerico"]["rtol"], p.RTOL)]
    for a, b in pares:
        assert a == pytest.approx(b, rel=1e-12)


def test_v1_ideal_coincide_con_solucion_analitica():
    c = B.simular("V1", "Ideal", P100)
    p = c.p
    tf_exacto = (p.P_OBJ - p.P0) * p.V / (p.R * p.GAMMA * p.TS * p.M_DOT_V1)
    assert c.tf == pytest.approx(tf_exacto, abs=1e-8)
    t = np.linspace(0, c.tf, 400)
    m, tg = c.estados(t)[:2]
    exacta = p.GAMMA * p.TS + c.m0 / m * (p.T0 - p.GAMMA * p.TS)
    assert np.max(np.abs(tg - exacta)) < 1e-4


def test_caudal_continuo_en_la_razon_critica():
    p = B.BASE
    pc = p.PS * B.relacion_critica(p.GAMMA)
    izq, der = B.caudal(pc * (1 - 1e-9), p, "V2"), B.caudal(pc * (1 + 1e-9), p, "V2")
    assert izq == pytest.approx(der, rel=1e-6)
    assert B.caudal(0.0, p, "V2") == pytest.approx(0.020, rel=1e-12)
    assert B.caudal(p.PS, p, "V2") == 0.0


@pytest.mark.parametrize("ver", ["V1", "V2", "V3"])
def test_vdw_con_a_b_nulos_es_ideal(ver):
    p = P100 if ver == "V1" else P190
    ideal = B.simular(ver, "Ideal", p)
    vdw0 = B.simular(ver, "Van der Waals", replace(p, A=0.0, B=0.0))
    assert vdw0.tf == pytest.approx(ideal.tf, abs=1e-9)


@pytest.mark.parametrize("modelo", ["Ideal", "Van der Waals", PRM.PR])
def test_v3_sin_conductancias_es_v2(modelo):
    v2 = PRM.sim("V2", modelo, P190)
    v3 = PRM.sim("V3", modelo, replace(P190, UA_GW=0.0, UA_WA=0.0))
    assert v3.tf == pytest.approx(v2.tf, abs=1e-8)


def test_pr_propiedades_y_limite_ideal():
    filas = PRM.verificar_propiedades()      # contiene sus propios assert
    assert len(filas) == 20
    lim = PRM.limite_ideal()
    assert lim[-1]["error_Z"] < 1e-12


def test_raices_en_rama_estable():
    for P in (1 * BAR, 100 * BAR, 200 * BAR, 250 * BAR):
        v = PRM.volumen_pr(P, 300.0)
        P_calc, *_, pv = PRM.propiedades(300.0, v)
        assert v > PRM.B and pv < 0 and P_calc == pytest.approx(P, rel=1e-9)
        vv = B.volumen_vdw(P, 300.0, B.BASE)
        assert vv > B.B_VDW and B.presion_eos("Van der Waals", B.BASE.V / vv, 300.0, B.BASE) == pytest.approx(P, rel=1e-9)


def test_dominio_invalido_se_informa():
    with pytest.raises(ValueError):
        B.presion_eos("Van der Waals", B.BASE.V / (0.5 * B.B_VDW), 300.0, B.BASE)
    with pytest.raises(ValueError):
        PRM.propiedades(100.0, 0.01)          # T < Tc
    with pytest.raises(ValueError):
        B.validar_parametros(replace(B.BASE, P_OBJ=250 * BAR), "V2")   # objetivo sobre el suministro


@pytest.mark.parametrize("modelo", ["Ideal", "Van der Waals"])
def test_balances_cierran_v3(modelo, tmp_path):
    r = B.validacion_balances(B.simular("V3", modelo, P190), tmp_path)
    assert r["masa_rel"] < 1e-9 and r["global_rel"] < 1e-9


def test_balances_cierran_pr_v3():
    r = PRM.balances(PRM.simular_pr("V3", P190))
    assert r["masa_kg_rel"] < 1e-9 and r["global_J_rel"] < 1e-9


def test_formulacion_conservativa_pr():
    fila = PRM.contrastar_integradores(PRM.simular_pr("V3", P190))
    assert fila["error_tf_s"] < 1e-4 and fila["error_T_K"] < 2e-5


@pytest.mark.parametrize("ver,modelo,tf,Tf,mf", [
    ("V1", "Ideal", 198.56, 418.30, 4.027), ("V1", "Van der Waals", 213.36, 387.92, 4.323),
    ("V1", PRM.PR, 208.53, 389.49, 4.227), ("V2", "Ideal", 418.92, 419.09, 7.637),
    ("V2", "Van der Waals", 410.86, 403.67, 7.595), ("V2", PRM.PR, 404.88, 402.67, 7.471),
    ("V3", "Ideal", 456.40, 385.32, 8.307), ("V3", "Van der Waals", 440.75, 379.44, 8.132),
    ("V3", PRM.PR, 433.60, 378.59, 7.985)])
def test_cifras_de_la_presentacion(ver, modelo, tf, Tf, mf):
    """Las cifras que se muestran en la presentación salen de esta misma ejecución."""
    r = PRM.resumir(PRM.sim(ver, modelo, P100 if ver == "V1" else P190))
    assert r["tf_s"] == pytest.approx(tf, abs=0.006)
    assert r["Tf_K"] == pytest.approx(Tf, abs=0.006)
    assert r["mf_kg"] == pytest.approx(mf, abs=0.0006)
    assert math.isfinite(r["Zf"]) and r["Zf"] > 0
