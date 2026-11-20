"""ICF130 · Proyecto computacional · Tema 7
Llenado transiente de un depósito rígido con nitrógeno (gas ideal, Van der Waals y Peng-Robinson).

Uso rápido (desde esta carpeta):
    python main.py                      todo: cálculos, validaciones, figuras y video
    python main.py --sin-video          todo menos el video
    python main.py --solo-video         solo la animación MP4
    python main.py --caso V3 --eos pr   un caso individual, con resumen en pantalla
    python main.py --salida OTRA_CARPETA
"""
from __future__ import annotations

import argparse
import json
import platform
import sys
import time
from dataclasses import asdict, replace
from pathlib import Path

import numpy as np
import scipy
import matplotlib

from modelos import base as B
from modelos import peng_robinson as PRM
from modelos import sensibilidad_pr, figuras_presentacion, correcciones
from modelos.animacion import crear_video

AQUI = Path(__file__).resolve().parent
EOS = {"ideal": "Ideal", "vdw": "Van der Waals", "pr": PRM.PR}


def parametros_version(ver: str):
    return replace(B.BASE, P_OBJ=(100 if ver == "V1" else 190) * B.BAR)


def correr_base(datos: Path, figuras: Path) -> dict:
    """Etapa ideal-Van der Waals de la entrega anterior (mismo orden de cálculo)."""
    datos.mkdir(parents=True, exist_ok=True)
    figuras.mkdir(parents=True, exist_ok=True)
    vmap = {ver: [B.simular(ver, mod, parametros_version(ver)) for mod in B.MODELOS]
            for ver in B.VERSIONES}
    todos = [c for par in vmap.values() for c in par]
    B.csv_file(datos / "resumen.csv", [B.enriquecer_resumen(c) for c in todos])
    balances = [B.validacion_balances(c, datos) for c in todos]
    B.csv_file(datos / "balances.csv", balances)
    B.csv_file(datos / "analitica.csv", [B.validar_ideal(vmap[v][0]) for v in ("V1", "V2")])
    B.figura_balances(balances, figuras)
    B.limites(vmap, datos)
    B.validar_evento(todos, datos)
    B.convergencia(vmap["V3"], datos, figuras)
    B.prueba_rk4_orden(B.BASE, datos, figuras)
    B.barridos(datos, figuras)
    B.figuras_base(vmap, datos, figuras)
    return {"balances": balances}


def casos_nueve() -> dict[str, list]:
    return {ver: [PRM.sim(ver, mod, parametros_version(ver)) for mod in PRM.MODELS]
            for ver in B.VERSIONES}


def caso_individual(ver: str, eos: str, salida: Path) -> None:
    modelo = EOS[eos]
    c = PRM.sim(ver, modelo, parametros_version(ver))
    r = PRM.resumir(c)
    carpeta = salida / "Casos_individuales"
    carpeta.mkdir(parents=True, exist_ok=True)
    PRM.trayectoria(c, carpeta / f"{ver}_{eos}.csv")
    print(f"\n{ver} con {modelo}")
    print(f"  tiempo para llegar a {r['objetivo_bar']:.0f} bar : {r['tf_s']:.2f} s")
    print(f"  temperatura final del gas      : {r['Tf_K']:.2f} K ({r['Tf_K'] - 273.15:.2f} °C)")
    print(f"  temperatura máxima del gas     : {r['Tmax_K']:.2f} K en t = {r['t_Tmax_s']:.1f} s")
    print(f"  masa final                     : {r['mf_kg']:.4f} kg (inicial {r['m0_kg']:.5f} kg)")
    print(f"  factor Z final                 : {r['Zf']:.4f}")
    print(f"  caudal inicial / final         : {r['mdot_0_g_s']:.2f} / {r['mdot_f_g_s']:.2f} g/s")
    if ver == "V3":
        print(f"  pared final                    : {r['Twf_K']:.2f} K")
        print(f"  calor gas->pared / pared->amb  : {r['Qgw_kJ']:.2f} / {r['Qwa_kJ']:.2f} kJ")
    print(f"  serie temporal guardada en     : {carpeta / f'{ver}_{eos}.csv'}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--salida", type=Path, default=AQUI / "resultados")
    modo = ap.add_mutually_exclusive_group()
    modo.add_argument("--sin-video", action="store_true")
    modo.add_argument("--solo-video", action="store_true")
    ap.add_argument("--caso", choices=B.VERSIONES)
    ap.add_argument("--eos", choices=tuple(EOS), default="pr")
    ap.add_argument("--fps", type=int, default=24)
    ap.add_argument("--frames-por-version", type=int, default=120)
    a = ap.parse_args()
    salida = a.salida.resolve()
    if a.caso:
        caso_individual(a.caso, a.eos, salida)
        return
    if a.frames_por_version < 10 or a.fps < 1:
        ap.error("Se necesitan al menos 10 fotogramas por versión y 1 fps.")

    datos, graf = salida / "Datos", salida / "Graficos"
    registro = {"inicio": time.strftime("%Y-%m-%d %H:%M:%S"), "python": sys.version.split()[0],
                "plataforma": platform.platform(), "numpy": np.__version__, "scipy": scipy.__version__,
                "matplotlib": matplotlib.__version__, "CoolProp": PRM.coolprop_version,
                "parametros_base": {k: v for k, v in asdict(B.BASE).items()},
                "constantes": {"R_J_kgK": B.R_N2, "cv_ig_J_kgK": B.CV_N2, "gamma": B.GAMMA,
                               "Tc_K": B.TC_N2, "Pc_Pa": B.PC_N2, "omega": PRM.OMEGA,
                               "a_vdw": B.A_VDW, "b_vdw": B.B_VDW, "a_c_pr": PRM.AC, "b_pr": PRM.B},
                "etapas": {}}
    t0 = time.time()
    if not a.solo_video:
        print("[1/4] Etapa ideal y Van der Waals (casos base, validación, convergencia, barridos)", flush=True)
        correr_base(datos / "base", graf / "base")
        registro["etapas"]["base_s"] = round(time.time() - t0, 1)

        t1 = time.time()
        print("[2/4] Peng-Robinson, contraste con Span y malla de suministro", flush=True)
        PRM.ejecutar(datos / "peng_robinson", graf / "peng_robinson")
        registro["etapas"]["peng_robinson_s"] = round(time.time() - t1, 1)

        t1 = time.time()
        print("[3/4] Sensibilidad de pared con las tres EOS", flush=True)
        sensibilidad_pr.ejecutar(datos / "sensibilidad_3EOS", graf / "sensibilidad_3EOS")
        figuras_presentacion.generar(datos, graf / "presentacion")
        print("[Extra] Comparación controlada y capacidad térmica de pared", flush=True)
        correcciones.ejecutar(datos / "correcciones", graf / "correcciones")
        registro["etapas"]["sensibilidad_y_figuras_s"] = round(time.time() - t1, 1)
    if not a.sin_video:
        t1 = time.time()
        print("[4/4] Animación MP4 (V1-V3, tres EOS)", flush=True)
        crear_video(casos_nueve(), salida / "Animaciones_MP4" / "ICF130_V1_V2_V3_tres_EOS.mp4",
                    frames_por_version=a.frames_por_version, fps=a.fps)
        registro["etapas"]["video_s"] = round(time.time() - t1, 1)
    registro["total_s"] = round(time.time() - t0, 1)
    salida.mkdir(parents=True, exist_ok=True)
    nombre = "registro_ejecucion.json" if not a.solo_video else "registro_video.json"
    (salida / nombre).write_text(json.dumps(registro, indent=2, ensure_ascii=False, default=str) + "\n",
                                 encoding="utf-8")
    print(f"Listo en {registro['total_s']} s. Resultados en {salida}", flush=True)


if __name__ == "__main__":
    main()
