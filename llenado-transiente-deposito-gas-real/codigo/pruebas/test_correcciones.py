"""Pruebas de comparación controlada y sensibilidad de pared."""
from dataclasses import replace
import pytest
from modelos.base import BASE, relacion_critica
from modelos.peng_robinson import BAR, MODELS, sim, resumir

@pytest.mark.parametrize('modelo',MODELS)
def test_v1_v2_coinciden_hasta_100bar(modelo):
    p=replace(BASE,P_OBJ=100*BAR)
    a=resumir(sim('V1',modelo,p));b=resumir(sim('V2',modelo,p))
    assert abs(a['tf_s']-b['tf_s'])<1e-6
    assert abs(a['Tf_K']-b['Tf_K'])<1e-6
    assert abs(a['mf_kg']-b['mf_kg'])<1e-8

@pytest.mark.parametrize('modelo',MODELS)
def test_a_190bar_v2_tarda_mas_pero_mismo_estado(modelo):
    p=replace(BASE,P_OBJ=190*BAR)
    a=resumir(sim('V1',modelo,p));b=resumir(sim('V2',modelo,p))
    assert b['tf_s']>a['tf_s']
    assert abs(a['Tf_K']-b['Tf_K'])<1e-5
    assert abs(a['mf_kg']-b['mf_kg'])<1e-8

@pytest.mark.parametrize('modelo',MODELS)
def test_pared_mas_masiva_enfria_mas(modelo):
    a=resumir(sim('V3',modelo,replace(BASE,MW=8)))
    b=resumir(sim('V3',modelo,replace(BASE,MW=54)))
    assert b['Tf_K']<a['Tf_K']
    assert b['tf_s']>a['tf_s']
