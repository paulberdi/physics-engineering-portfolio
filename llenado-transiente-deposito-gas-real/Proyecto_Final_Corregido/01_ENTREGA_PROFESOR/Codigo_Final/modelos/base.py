"""ICF130. Núcleo del modelo: gas ideal y Van der Waals, etapas V1-V3.

Copia del programa de la entrega anterior (ICF130_PROYECTO_FINAL.py), sin cambios
en la física ni en los métodos numéricos. Lo usan main.py y peng_robinson.py.
La animación antigua (crear_mp4, mapa_ilustrativo) se conserva solo como
antecedente; la animación final está en modelos/animacion.py.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
import shutil
import sys
import scipy
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq, minimize_scalar

# ───────────────────────────────────────────────────────────────────────────────
# 1. CONSTANTES FÍSICAS Y ESCENARIO BASE
# ───────────────────────────────────────────────────────────────────────────────
RU = 8.314462618  # J mol-1 K-1
M_N2 = 28.0134e-3  # kg mol-1
R_N2 = RU / M_N2  # J kg-1 K-1
CP_N2 = 1039.0  # J kg-1 K-1; valor constante de estudio ~300 K
CV_N2 = CP_N2 - R_N2
GAMMA = CP_N2 / CV_N2
TC_N2 = 126.192  # K
PC_N2 = 3.3958e6  # Pa
# Constantes VdW por unidad de masa (a y b para v=V/m):
A_VDW = 27.0 * R_N2**2 * TC_N2**2 / (64.0 * PC_N2)
B_VDW = R_N2 * TC_N2 / (8.0 * PC_N2)
MODELOS = ("Ideal", "Van der Waals")
COLORS = {"Ideal": "#056992", "Van der Waals": "#A34A08"}


def relacion_critica(gamma: float) -> float:
    return (2.0 / (gamma + 1.0)) ** (gamma / (gamma - 1.0))


def factor_estrangulado(r: float, g: float, ts: float) -> float:
    return math.sqrt(g / (r * ts)) * (2.0 / (g + 1.0)) ** ((g + 1.0) / (2.0 * (g - 1.0)))


AREA_REFERENCIA = 0.020 / (0.8 * 200e5 * factor_estrangulado(R_N2, GAMMA, 300.0))


@dataclass(frozen=True)
class Parametros:
    V: float = 0.050  # m3
    P0: float = 1.0e5  # Pa
    T0: float = 300.0  # K
    PS: float = 200e5  # Pa
    TS: float = 300.0  # K
    P_OBJ: float = 190e5  # Pa (V2, V3)
    M_DOT_V1: float = 0.020  # kg/s
    CD: float = 0.80  # coeficiente de descarga
    AREA: float = AREA_REFERENCIA  # m2, FIJA en barridos de PS
    TINF: float = 300.0  # K
    TW0: float = 300.0  # K
    MW: float = 8.0  # kg efectivos de pared
    CPW: float = 500.0  # J/(kg K)
    UA_GW: float = 8.0  # W/K
    UA_WA: float = 3.0  # W/K
    CP: float = CP_N2
    CV: float = CV_N2
    R: float = R_N2
    GAMMA: float = GAMMA
    A: float = A_VDW
    B: float = B_VDW
    RTOL: float = 1e-9
    ATOL: float = 1e-11
    MAX_STEP: float = 0.5
    METHOD: str = "RK45"
    TIME_LIMIT: float = 10000.0

    @property
    def CW(self) -> float:
        return self.MW * self.CPW


BASE = Parametros()


def validar_parametros(p: Parametros, version: str) -> None:
    if version not in ("V1", "V2", "V3"):
        raise ValueError("Versión desconocida")
    if not (0 < p.P0 < p.P_OBJ < p.PS):
        raise ValueError(f"Rango de presiones inválido: P0={p.P0}, Pobj={p.P_OBJ}, PS={p.PS}")
    if min(p.V, p.T0, p.TS, p.CP, p.CV, p.R, p.AREA, p.M_DOT_V1,
           p.CW, p.RTOL, p.ATOL, p.MAX_STEP, p.TIME_LIMIT) <= 0:
        raise ValueError("Algún parámetro positivo es cero/negativo")
    if p.A < 0 or p.B < 0 or p.UA_GW < 0 or p.UA_WA < 0:
        raise ValueError("a, b y las conductancias deben ser no negativas")
    if not (0 < p.CD <= 1 and p.TW0 > 0 and p.TINF > 0):
        raise ValueError("Cd o temperaturas de pared/ambiente inválidas")


def volumen_vdw(P: float, T: float, p: Parametros) -> float:
    """Resuelve la EOS en densidad para estados de gas con T mayor que la crítica."""
    if p.A == 0 and p.B == 0:
        return p.R * T / P
    if p.A and p.B and T <= 8 * p.A / (27 * p.R * p.B):
        raise ValueError("La inversión de Van der Waals requiere la rama supercrítica.")

    def f(rho):
        return p.R * T * rho / (1 - p.B * rho) - p.A * rho**2 - P

    hi = (1 - 1e-10) / p.B if p.B else 2 * P / (p.R * T)
    if f(hi) <= 0:
        raise ValueError("No se acotó la raíz en densidad de Van der Waals.")
    return 1 / brentq(f, 0, hi, xtol=1e-12, rtol=1e-14)


def presion_eos(modelo: str, m: Any, Tg: Any, p: Parametros) -> Any:
    if modelo == "Ideal":
        return m * p.R * Tg / p.V
    v = p.V / m
    if np.any(v <= p.B):
        raise ValueError("Van der Waals fuera de dominio: volumen específico <= b")
    return p.R * Tg / (v - p.B) - p.A / (v**2)


def energia_gas(modelo: str, m: Any, Tg: Any, p: Parametros) -> Any:
    u = m * p.CV * Tg
    return u if modelo == "Ideal" else u - p.A * m**2 / p.V


def entalpia_entrada(modelo: str, p: Parametros) -> float:
    if modelo == "Ideal":
        return p.CP * p.TS
    vs = volumen_vdw(p.PS, p.TS, p)
    return p.CV * p.TS - p.A / vs + p.PS * vs


def masa_inicial(modelo: str, p: Parametros) -> float:
    if modelo == "Ideal":
        return p.P0 * p.V / (p.R * p.T0)
    return p.V / volumen_vdw(p.P0, p.T0, p)


def caudal(Ptan: float, p: Parametros, version: str) -> float:
    """Ecuación de flujo compresible ideal, no una ley constitutiva general."""
    if version == "V1":
        return p.M_DOT_V1
    if Ptan >= p.PS:
        return 0.0
    r = float(max(0.0, Ptan / p.PS))
    crit = relacion_critica(p.GAMMA)
    if r <= crit:
        return p.CD * p.AREA * p.PS * factor_estrangulado(p.R, p.GAMMA, p.TS)
    term = r**(2.0 / p.GAMMA) * -math.expm1((p.GAMMA-1) / p.GAMMA * math.log(r))
    return p.CD * p.AREA * p.PS / math.sqrt(p.R * p.TS) * math.sqrt(
        2.0 * p.GAMMA / (p.GAMMA - 1.0) * term)


@dataclass
class Caso:
    version: str
    modelo: str
    p: Parametros
    m0: float
    h_in: float
    sol: Any
    tf: float

    def estados(self, t: Any) -> np.ndarray:
        tt = np.asarray(t, dtype=float)
        return np.asarray(self.sol.sol(np.clip(tt, 0, self.tf)))

    def muestras(self, n: int = 1201) -> dict[str, np.ndarray]:
        t = np.linspace(0, self.tf, n)
        return self.muestras_en(t)

    def muestras_en(self, t: np.ndarray) -> dict[str, np.ndarray]:
        y = self.estados(t)
        m, tg, tw = y[0], y[1], y[2]
        pres = np.asarray(presion_eos(self.modelo, m, tg, self.p))
        md = np.array([caudal(float(x), self.p, self.version) for x in pres])
        qgw = self.p.UA_GW * (tg - tw) if self.version == "V3" else np.zeros_like(t)
        qwa = self.p.UA_WA * (tw - self.p.TINF) if self.version == "V3" else np.zeros_like(t)
        z = pres * self.p.V / (m * self.p.R * tg)
        return {"t":t, "m":m, "Tg":tg, "Tw":tw, "P":pres, "mdot":md,
                "q_gw":qgw, "q_wa":qwa, "Z":z, "Ein":y[3], "Qgw":y[4], "Qwa":y[5]}


def simular(version: str, modelo: str, p: Parametros, horizonte: float | None = None) -> Caso:
    validar_parametros(p, version)
    if modelo not in MODELOS:
        raise ValueError("Modelo desconocido")
    m0 = masa_inicial(modelo, p)
    hin = entalpia_entrada(modelo, p)
    def rhs(_t: float, y: np.ndarray) -> list[float]:
        m, tg, tw, _ein, _qgw, _qwa = y
        pp = float(presion_eos(modelo, m, tg, p))
        md = caudal(pp, p, version)
        q_gw = p.UA_GW * (tg - tw) if version == "V3" else 0.0
        q_wa = p.UA_WA * (tw - p.TINF) if version == "V3" else 0.0
        # Ugas(VdW)=m*cv*T-a*m²/V => 
        # m*cv*dT/dt = mdot*(hin-cv*T+2*a*m/V)-Q_gw
        corr = 2.0 * p.A * m / p.V if modelo == "Van der Waals" else 0.0
        dtg = (md * (hin - p.CV * tg + corr) - q_gw) / (m * p.CV)
        dtw = (q_gw - q_wa) / p.CW if version == "V3" else 0.0
        return [md, dtg, dtw, md * hin, q_gw, q_wa]
    def event(_t: float, y: np.ndarray) -> float:
        return float(presion_eos(modelo, y[0], y[1], p) - p.P_OBJ)
    event.terminal = True
    event.direction = 1
    y0 = [m0, p.T0, p.TW0, 0., 0., 0.]
    sol = solve_ivp(rhs, (0., p.TIME_LIMIT if horizonte is None else horizonte), y0,
                    dense_output=True, events=event if horizonte is None else None, method=p.METHOD,
                    rtol=p.RTOL, atol=p.ATOL, max_step=p.MAX_STEP)
    if (not sol.success) or (horizonte is None and not len(sol.t_events[0])):
        raise RuntimeError(f"La simulación {version}-{modelo} no alcanzó P objetivo: {sol.message}")
    return Caso(version, modelo, p, m0, hin, sol, float(sol.t_events[0][0] if horizonte is None else sol.t[-1]))

# ───────────────────────────────────────────────────────────────────────────────
# 2. EXPORTACIÓN Y VALIDACIÓN
# ───────────────────────────────────────────────────────────────────────────────
def csv_file(path: Path, filas: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not filas:
        return
    keys = list(dict.fromkeys(k for r in filas for k in r))
    with path.open("w", encoding="utf-8-sig", newline="") as f:
        wr = csv.DictWriter(f, keys)
        wr.writeheader()
        wr.writerows(filas)


def save_samples(c: Caso, path: Path, n: int = 1601) -> None:
    d = c.muestras(n)
    rows = [{"t_s":float(d["t"][i]), "P_bar":float(d["P"][i]/1e5),
             "Tg_K":float(d["Tg"][i]), "Tw_K":float(d["Tw"][i]) if c.version=="V3" else "",
             "masa_kg":float(d["m"][i]), "mdot_kg_s":float(d["mdot"][i]),
             "Qdot_gw_W":float(d["q_gw"][i]), "Qdot_wa_W":float(d["q_wa"][i]),
             "Z":float(d["Z"][i]), "regimen":regimen(float(d["P"][i]),c.p,c.version), "Eing_J":float(d["Ein"][i]),
             "Qgw_acum_J":float(d["Qgw"][i]), "Qwa_acum_J":float(d["Qwa"][i])}
            for i in range(n)]
    csv_file(path, rows)


def resumen_caso(c: Caso) -> dict[str, Any]:
    data = c.muestras(1601)
    critical_pressure = relacion_critica(c.p.GAMMA) * c.p.PS
    t_transition = ""
    if c.version != "V1" and c.p.P0 < critical_pressure < c.p.P_OBJ:
        t_transition = float(brentq(lambda t: float(presion_eos(c.modelo, c.estados(t)[0],
                                                  c.estados(t)[1], c.p) - critical_pressure), 0, c.tf))
    return {"version":c.version,"modelo":c.modelo,"tiempo_final_s":c.tf,
            "P_final_bar":float(data["P"][-1]/1e5),
            "Tgas_final_C":float(data["Tg"][-1]-273.15),
            "Tgas_max_C":float(np.max(data["Tg"])-273.15),
            "Tpared_final_C":float(data["Tw"][-1]-273.15) if c.version=="V3" else "",
            "masa_final_kg":float(data["m"][-1]),
            "caudal_final_kg_s":float(data["mdot"][-1]),
            "Z_final":float(data["Z"][-1]),
            "transicion_subsonico_s":t_transition,
            "mdot_inicial_kg_s":float(data["mdot"][0]),
            "Qgas_pared_acum_kJ":float(data["Qgw"][-1]/1000),
            "Qpared_ambiente_acum_kJ":float(data["Qwa"][-1]/1000)}


def prueba_rk4_orden(p: Parametros, out: Path, figuras: Path) -> list[dict[str, Any]]:
    """Orden experimental de RK4 paso FIJO en la EDO V3 ideal hasta t=180 s.

    Se elige t=180 s, antes del cambio de régimen, de modo que la EDO es suave.
    Se compara contra DOP853 con tolerancias estrictas, a IGUAL tiempo.
    RK4 es una prueba numérica adicional; V3 se ejecuta normalmente con RK45.
    """
    t_end = 180.0
    ref = simular("V3", "Ideal", replace(p, METHOD="DOP853",
                     RTOL=2e-13, ATOL=1e-14, MAX_STEP=0.03))
    y_ref = ref.estados(t_end)
    m0 = masa_inicial("Ideal",p)
    hin = entalpia_entrada("Ideal",p)
    def deriv(y):
        m,tg,tw,_,_,_=y
        P = float(presion_eos("Ideal",m,tg,p))
        md = caudal(P,p,"V3")
        qgw = p.UA_GW*(tg-tw)
        qwa = p.UA_WA*(tw-p.TINF)
        return np.array([md, (md*(hin-p.CV*tg)-qgw)/(m*p.CV),
                         (qgw-qwa)/p.CW, md*hin, qgw, qwa])
    rows=[]
    for ns in (90,180,360,720,1440,2880):
        dt=t_end/ns
        y=np.array([m0,p.T0,p.TW0,0.,0.,0.],dtype=float)
        for _ in range(ns):
            k1=deriv(y)
            k2=deriv(y+dt*k1/2)
            k3=deriv(y+dt*k2/2)
            k4=deriv(y+dt*k3)
            y += dt*(k1+2*k2+2*k3+k4)/6
        rows.append({"n_pasos":ns,"paso_fijo_s":dt,
                     "error_Tgas_K":abs(float(y[1]-y_ref[1])),
                     "error_masa_kg":abs(float(y[0]-y_ref[0]))})
    for i,r in enumerate(rows):
        if i==0:
            r["orden_estimado_T"]=""
        else:
            ep=rows[i-1]["error_Tgas_K"]
            en=r["error_Tgas_K"]
            r["orden_estimado_T"]=(math.log2(ep/en) if min(ep,en)>1e-11
                                    else "no confiable: redondeo")
    csv_file(out/"orden_RK4_paso_fijo_V3.csv",rows)
    fig,ax=plt.subplots(figsize=(8,4.55))
    ax.loglog([r["paso_fijo_s"] for r in rows],
              [r["error_Tgas_K"] for r in rows],"o-",label="RK4 vs DOP853",
              color=COLORS["Ideal"])
    fig_style(ax,"Paso fijo h [s]","Error absoluto Tgas(180 s) [K]",
              "Convergencia de orden 4 en V3-ideal (antes de transición)")
    export_fig(fig,figuras/"orden_RK4_V3.png")
    return rows

# ───────────────────────────────────────────────────────────────────────────────
# 3. GRÁFICOS PARA LATEX — una magnitud por figura
# ───────────────────────────────────────────────────────────────────────────────
def fig_style(ax, xlabel: str, ylabel: str, title: str) -> None:
    ax.set(xlabel=xlabel,ylabel=ylabel,title=title)
    ax.grid(True,color="#9AA6B4",alpha=0.22)
    ax.spines[["top","right"]].set_visible(False)
    ax.legend(loc="best",fontsize=9)


def export_fig(fig, path: Path) -> None:
    path.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(path, dpi=300, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def figura_version(casos: list[Caso], variable: str, ylabel: str, path: Path) -> None:
    fig,ax=plt.subplots(figsize=(8,4.55))
    for c in casos:
        d=c.muestras(1001)
        dmap={"temperatura":d["Tg"]-273.15,"presion":d["P"]/1e5,
              "masa":d["m"],"caudal":d["mdot"],"factor_z":d["Z"],
              "calor":d["q_gw"],"calor_ambiente":d["q_wa"],"pared":d["Tw"]-273.15}
        ax.plot(d["t"],dmap[variable],color=COLORS[c.modelo],lw=2.2,label=c.modelo)
    fig_style(ax,"Tiempo [s]",ylabel,f"{casos[0].version}: {ylabel.split(' [')[0]}")
    export_fig(fig,path)


def comparar_v2_v3(v2: list[Caso],v3: list[Caso],out: Path) -> None:
    fig,ax=plt.subplots(figsize=(8,4.55))
    for vv,li in ((v2,"--"),(v3,"-")):
        for c in vv:
            d=c.muestras(1001)
            ax.plot(d["t"],d["Tg"]-273.15,li,color=COLORS[c.modelo],lw=2,
                    label=f"{c.version} — {c.modelo}")
    fig_style(ax,"Tiempo [s]","Temperatura del gas [°C]","V2 frente a V3: enfriamiento por la pared")
    export_fig(fig,out)


def grafico_barrido(rows: list[dict[str,Any]],xcol: str,ycol: str,
                   ylabel: str,xlabel: str,title: str,out:Path) -> None:
    fig,ax=plt.subplots(figsize=(8,4.55))
    for mod in MODELOS:
        rr=sorted((r for r in rows if r["modelo"]==mod),key=lambda r:r[xcol])
        ax.plot([r[xcol] for r in rr],[r[ycol] for r in rr],"o-",lw=2,
                color=COLORS[mod],label=mod)
    fig_style(ax,xlabel,ylabel,title)
    export_fig(fig,out)


# 4. ANIMACIÓN ANTIGUA (no se usa en la entrega final; ver modelos/animacion.py)
# ───────────────────────────────────────────────────────────────────────────────
VIDEO_COLORS = {"Ideal": "#45C7F2", "Van der Waals": "#FFAB57"}
DARK="#080D15"; BG_AX="#111A27"; MUTED="#ABB8C9"; WHITE="#EFF4FA"


def mapa_ilustrativo(tg_C: float, t0_C: float, avance: float,
                     xs: np.ndarray, ys: np.ndarray) -> np.ndarray:
    """Textura térmica SIMBÓLICA tipo difusión, no resultado de PDE.

    Centra el campo en Tg(t); el patrón es decorativo y no realimenta la EDO.
    """
    xx,yy=np.meshgrid(xs,ys)
    sig=0.10+0.31*math.sqrt(max(float(avance),0.))
    jet=np.exp(-((xx-.085)/sig)**2)*np.exp(-((yy-.53)/.29)**2)
    brillo=0.30*np.sin(5.5*xx-1.9*avance)*np.cos(3.5*yy)
    perfil=jet+0.18*brillo
    perfil-=np.mean(perfil)
    delta=(tg_C-t0_C)
    return np.asarray(tg_C+0.37*delta*perfil,dtype=np.float32)


def crear_mp4(casos_por_version: dict[str,list[Caso]],out:Path,
               frames_por_version:int=64,fps:int=16,dpi:int=120) -> None:
    if shutil.which("ffmpeg") is None:
        raise RuntimeError("ffmpeg no está instalado/en PATH. Instálalo o usa --sin-video.")
    versions=("V1","V2","V3")
    # Precalcula todo para evitar resolver EDO por fotograma.
    DATA={}
    for ver in versions:
        cc=casos_por_version[ver]
        tfin=max(c.tf for c in cc)
        tgrid=np.linspace(0,tfin,frames_por_version)
        d=[]
        for c in cc:
            ts=np.minimum(tgrid,c.tf)
            y=c.estados(ts)
            m,tg,tw=y[:3]
            pp=np.asarray(presion_eos(c.modelo,m,tg,c.p))
            md=np.array([caudal(float(x),c.p,c.version) for x in pp])
            qw=c.p.UA_GW*(tg-tw) if ver=="V3" else np.zeros(len(tgrid))
            qa=c.p.UA_WA*(tw-c.p.TINF) if ver=="V3" else np.zeros(len(tgrid))
            d.append({"qa":qa,"m":m,"Tg":tg,"Tw":tw,"P":pp,"md":md,"qw":qw,"t":tgrid,
                      "tloc":ts,"tf":c.tf,"Z":pp*c.p.V/(m*c.p.R*tg)})
        DATA[ver]=(d,tfin)
    fig=plt.figure(figsize=(16,9),dpi=dpi,facecolor=DARK)
    gs=fig.add_gridspec(3,2,left=.06,right=.92,top=.93,bottom=.085,
                        height_ratios=(.14, .38,.48),hspace=.28,wspace=.13)
    ax_header=fig.add_subplot(gs[0,:]);ax_header.axis("off")
    ax_l=fig.add_subplot(gs[1,0]);ax_r=fig.add_subplot(gs[1,1])
    sub=gs[2,:].subgridspec(2,3,wspace=.23,hspace=.53)
    axs=[fig.add_subplot(sub[i,j]) for i in range(2) for j in range(3)]
    for ax in [ax_l,ax_r]+axs:
        ax.set_facecolor(BG_AX)
        for sp in ax.spines.values():sp.set_color("#2D3B4C")
        ax.tick_params(colors=MUTED,labelsize=9)
        ax.title.set_color(WHITE)
        ax.xaxis.label.set_color(MUTED);ax.yaxis.label.set_color(MUTED)
    ax_header.text(0,1,"Llenado transiente de un depósito · N₂",
                   color=WHITE,fontsize=22,weight="bold",va="top",transform=ax_header.transAxes)
    subtitle=ax_header.text(0,.37,"",color=MUTED,fontsize=12,transform=ax_header.transAxes)
    clock=ax_header.text(1,.90,"",color="#C6E7F4",fontsize=16,weight="bold",ha="right",
                          transform=ax_header.transAxes)
    labels=("GAS IDEAL","VAN DER WAALS")
    for ax,lab,mod in ((ax_l,labels[0],"Ideal"),(ax_r,labels[1],"Van der Waals")):
        ax.set_title(lab,color=VIDEO_COLORS[mod],fontsize=14,weight="bold",pad=8)
        ax.set_xticks([]);ax.set_yticks([])
        ax.set_xlim(0,1);ax.set_ylim(0,1)
    xs=np.linspace(0,1,190);ys=np.linspace(0,1,66)
    zero=np.full((len(ys),len(xs)),BASE.T0-273.15,dtype=np.float32)
    images=[ax.imshow(zero,cmap="inferno",vmin=20,vmax=175,origin="lower",
                      aspect="auto",extent=(0,1,0,1),interpolation="bicubic") for ax in (ax_l,ax_r)]
    # Salida térmica en entrada izquierda: flecha sobre cada depósito (ilustración)
    for ax in (ax_l,ax_r):
        ax.annotate("",xy=(.18,.5),xytext=(-.08,.5),xycoords="axes fraction",
                    arrowprops=dict(arrowstyle="-|>",lw=3,color="#C7EDFF"))
    info=[ax.text(.025,.07,"",color=WHITE,fontsize=11,ha="left",va="bottom",
                  transform=ax.transAxes,fontfamily="monospace",
                  bbox=dict(boxstyle="round,pad=.5",facecolor="#070B12",edgecolor="#4C6074",alpha=.83))
          for ax in (ax_l,ax_r)]
    cax=fig.add_axes((.935,.535,.013,.27))
    cb=fig.colorbar(images[0],cax=cax)
    cb.ax.tick_params(colors=MUTED,labelsize=8)
    cb.set_label("Color ilustrativo [°C]",color=MUTED,fontsize=9)
    titulos = ("Temperatura", "Presión del depósito", "Masa almacenada",
               "Caudal de entrada", "Compresibilidad", "Transferencia térmica")
    unidades = ("Temperatura [°C]", "Presión [bar]", "Masa [kg]",
                "Caudal [g/s]", "Factor Z [-]", "Potencia [W]")
    for ax, title, ylabel in zip(axs, titulos, unidades):
        ax.set_title(title,loc="left",fontsize=10,pad=5,color=WHITE)
        ax.set_xlabel("Tiempo [s]",fontsize=9)
        ax.set_ylabel(ylabel,fontsize=9)
        ax.grid(True,color="#667587",alpha=.20,lw=.65)
    for ax in axs:ax.set_xlim(0,1)
    lines=[]
    for ax in axs:
        pair=[]
        for mod in MODELOS:
            (line,)=ax.plot([],[],color=VIDEO_COLORS[mod],lw=2.1,label=mod)
            pair.append(line)
        lines.append(pair)
    pared_lineas = [axs[0].plot([], [], "--", color=VIDEO_COLORS[mod], lw=1.7)[0]
                    for mod in MODELOS]
    ambiente_lineas = [axs[5].plot([], [], "--", color=VIDEO_COLORS[mod], lw=1.7)[0]
                       for mod in MODELOS]
    fig.text(.06,.023,"Modelo 0D · Mapas tipo difusión ilustrativos (no CFD ni PDE) · Flujo compresible ideal aun para VdW",
             color=MUTED,fontsize=8.8)
    fig.text(.92,.023,"UTFSM  ·  2026-2",ha="right",color=MUTED,fontsize=8.8)

    current_ver=[None]
    def update(i):
        ver=versions[min(i//frames_por_version,2)]
        k=i%frames_por_version
        cases=casos_por_version[ver]
        dlist,tfin=DATA[ver]
        tnow=float(dlist[0]["t"][k])
        if current_ver[0]!=ver:
            current_ver[0]=ver
            descr={"V1":"Caudal constante · depósito adiabático · Pobj = 100 bar",
                   "V2":"Flujo estrangulado → subsónico · adiabático · Pobj = 190 bar",
                   "V3":"Flujo variable + almacenamiento de calor en pared · Pobj = 190 bar"}[ver]
            subtitle.set_text(f"{ver}  /  {descr}")
            for ax in axs:ax.set_xlim(0,tfin)
            t_max=max(float(np.max(d["Tg"]-273.15)) for d in dlist)
            p_max=max(float(np.max(d["P"]/1e5)) for d in dlist)
            md_max=max(float(np.max(d["md"]))*1000 for d in dlist)
            qw_max=max(float(np.max(d["qw"])) for d in dlist)
            m_max=max(float(np.max(d["m"])) for d in dlist)
            z_low=min(.995, min(float(np.min(d["Z"])) for d in dlist) - .006)
            z_high=max(1.005, max(float(np.max(d["Z"])) for d in dlist) + .006)
            for ax, limits in zip(axs, ((18, max(70, t_max+13)), (0,p_max*1.06),
                                       (0,m_max*1.08), (0,max(25,md_max*1.18)),
                                       (z_low,z_high), (0,max(10,qw_max*1.14)))):
                ax.set_ylim(*limits)
            axs[0].set_title("Temperatura · pared en trazos" if ver=="V3" else "Temperatura del gas",
                              loc="left", fontsize=10, color=WHITE, pad=5)
            axs[5].set_title("Potencia térmica · ambiente en trazos" if ver=="V3" else "Intercambio térmico nulo",
                              loc="left", fontsize=10, color=WHITE, pad=5)
        clock.set_text(f"{ver}   |   t = {tnow:6.1f} / {tfin:.1f} s")
        for case,d,im,text_artist in zip(cases,dlist,images,info):
            prog=min(1.,tnow/case.tf)
            tg_c=float(d["Tg"][k]-273.15)
            im.set_data(mapa_ilustrativo(tg_c,BASE.T0-273.15,prog,xs,ys))
            state = "OBJETIVO" if tnow >= case.tf else "ENTRADA"
            fase = regimen(float(d["P"][k]), case.p, ver)
            row=(f"P: {d['P'][k]/1e5:6.1f} bar | Tgas: {tg_c:6.1f} °C\n"
                 f"m: {d['m'][k]:6.3f} kg | Z: {d['Z'][k]:.4f} | {state}"
                 f"\n{fase}  |  Caudal: {d['md'][k]*1000:5.2f} g/s")
            if ver=="V3":row+=f"\nTpared: {d['Tw'][k]-273.15:5.1f} °C | Qgas→pared: {d['qw'][k]:5.1f} W"
            text_artist.set_text(row)
        for col,d in enumerate(dlist):
            sl=slice(0,k+1)
            signals=(d["Tg"][sl]-273.15,d["P"][sl]/1e5,d["m"][sl],
                     d["md"][sl]*1000,d["Z"][sl],d["qw"][sl])
            for ai,sig in enumerate(signals):
                lines[ai][col].set_data(d["tloc"][sl],sig)
            if ver == "V3":
                pared_lineas[col].set_data(d["tloc"][sl], d["Tw"][sl]-273.15)
                ambiente_lineas[col].set_data(d["tloc"][sl], d["qa"][sl])
            else:
                pared_lineas[col].set_data([], [])
                ambiente_lineas[col].set_data([], [])
        return images+info+[clock,subtitle]+[x for pair in lines for x in pair]+pared_lineas+ambiente_lineas

    movie=FuncAnimation(fig,update,frames=frames_por_version*3,interval=1000/fps,
                        blit=False,repeat=False)
    out.parent.mkdir(parents=True,exist_ok=True)
    writer=FFMpegWriter(fps=fps,bitrate=2600,codec="libx264",
                        extra_args=["-pix_fmt","yuv420p","-movflags","+faststart"])
    print(f"[VIDEO] Exportando {out.name} ({frames_por_version*3} fotogramas, {fps} fps)",flush=True)
    movie.save(str(out),writer=writer,dpi=dpi)
    update(frames_por_version*3-1)
    fig.savefig(out.with_suffix(".png"),dpi=130,facecolor=fig.get_facecolor(),bbox_inches="tight")
    plt.close(fig)
    print(f"[VIDEO] Completado: {out}",flush=True)

# ───────────────────────────────────────────────────────────────────────────────

# 5. Verificaciones numéricas y estudios.
VERSIONES = ("V1", "V2", "V3")
BAR = 1e5


def regimen(presion: float, p: Parametros, version: str) -> str:
    if version == "V1":
        return "Caudal constante"
    if presion >= p.PS:
        return "Sin entrada"
    return "Estrangulado" if presion <= p.PS * relacion_critica(p.GAMMA) else "Subsónico"


def etiqueta(c: Caso) -> str:
    return "ideal" if c.modelo == "Ideal" else "vdw"


def presion_t(c: Caso, t):
    m, tg = c.estados(t)[:2]
    return presion_eos(c.modelo, m, tg, c.p)


def pico_temperatura(c: Caso) -> tuple[float, float]:
    t = np.linspace(0, c.tf, 1201)
    temperatures = c.estados(t)[1]
    i = int(np.argmax(temperatures))
    if i in (0, len(t) - 1):
        return float(t[i]), float(temperatures[i])
    opt = minimize_scalar(lambda u: -c.estados(u)[1], bounds=(t[i - 1], t[i + 1]),
                          method="bounded", options={"xatol": 1e-8})
    return float(opt.x), float(-opt.fun)


def validacion_balances(c: Caso, out: Path) -> dict:
    # Cuadratura posterior de las tasas reconstruidas, independiente de los
    # acumuladores de la EDO. Se parte en cada paso y en el cambio de régimen.
    p = c.p
    pc = p.PS * relacion_critica(p.GAMMA)
    cortes = list(c.sol.t)
    if c.version != "V1" and p.P0 < pc < p.P_OBJ:
        cortes.append(brentq(lambda t: presion_t(c, t) - pc, 0, c.tf))
    t = np.unique(cortes)
    nodes, weights = np.polynomial.legendre.leggauss(8)
    mitad = np.diff(t) / 2
    eval_t = (t[:-1, None] + mitad[:, None] * (1 + nodes)).ravel()
    d = c.muestras_en(eval_t)
    tasas = (d["mdot"], d["mdot"] * c.h_in, d["q_gw"], d["q_wa"])
    masa, ein, qgw, qwa = [
        np.r_[0, np.cumsum(mitad * (rate.reshape(-1, len(nodes)) @ weights))]
        for rate in tasas
    ]
    m, tg, tw, ein_ode, qgw_ode, qwa_ode = c.estados(t)
    du = energia_gas(c.modelo, m, tg, p) - energia_gas(c.modelo, c.m0, p.T0, p)
    dw = p.CW * (tw - p.TW0) if c.version == "V3" else np.zeros_like(t)
    rm, rg, rw = m - c.m0 - masa, du - ein + qgw, dw - qgw + qwa
    rt = du + dw - ein + qwa
    nombre = f"residuos_{c.version}_{etiqueta(c)}"
    csv_file(out / f"{nombre}.csv", [
        dict(zip(("t_s", "Rm_kg", "Rg_J", "Rw_J", "Rglobal_J"), row))
        for row in zip(t, rm, rg, rw, rt)
    ])
    result = {"version": c.version, "modelo": c.modelo, "archivo": nombre,
              "discordancia_acumuladores_J": float(np.max(np.abs(
                  np.array([ein, qgw, qwa]) - np.array([ein_ode, qgw_ode, qwa_ode]))))}
    for label, units, error, scale in (
        ("masa", "kg", rm, masa[-1]), ("gas", "J", rg, abs(ein[-1]) + abs(qgw[-1])),
        ("pared", "J", rw, abs(qgw[-1]) + abs(qwa[-1])),
        ("global", "J", rt, abs(ein[-1]) + abs(qwa[-1])),
    ):
        max_abs = float(np.max(np.abs(error)))
        result[f"{label}_abs_{units}"] = max_abs
        result[f"{label}_escala_{units}"] = scale
        result[f"{label}_rel"] = max_abs / scale if scale > 0 else None
    return result


def validar_ideal(c: Caso) -> dict:
    p = c.p
    t = np.linspace(0, c.tf, 2401)
    m, tg = c.estados(t)[:2]
    exacta = p.GAMMA * p.TS + c.m0 / m * (p.T0 - p.GAMMA * p.TS)
    dp_dm = p.R * p.GAMMA * p.TS / p.V
    pc = p.PS * relacion_critica(p.GAMMA)
    cortes = [pc] if c.version == "V2" and p.P0 < pc < p.P_OBJ else None
    tf = quad(lambda P: 1 / (dp_dm * caudal(P, p, c.version)), p.P0, p.P_OBJ,
              points=cortes, epsabs=1e-10, epsrel=1e-12)[0]
    m_exact = c.m0 + (p.P_OBJ - p.P0) / dp_dm
    if c.version == "V1":
        # Solución explícita adicional en todo el intervalo.
        masa = c.m0 + p.M_DOT_V1 * t
        e_m = float(np.max(np.abs(m - masa)))
        e_p = float(np.max(np.abs(presion_t(c, t) - p.P0 - dp_dm * p.M_DOT_V1 * t)))
        exacta = p.GAMMA * p.TS + c.m0 / masa * (p.T0 - p.GAMMA * p.TS)
    else:
        e_m = abs(m[-1] - m_exact)
        e_p = None
    e_T = float(np.max(np.abs(tg - exacta)))
    return {"version": c.version, "tf_referencia_s": tf, "error_tf_s": abs(c.tf - tf),
            "error_tf_rel": abs(c.tf - tf) / tf, "error_T_K": e_T,
            "error_T_rel": e_T / float(np.max(exacta)), "error_m_kg": e_m,
            "error_m_rel": e_m / m_exact, "error_P_Pa": e_p,
            "error_P_rel": e_p / p.P_OBJ if e_p is not None else None}


def comparar(c: Caso, ref: Caso) -> dict:
    """Errores finales en el evento y máximos de trayectoria a igual tiempo."""
    t = np.linspace(0, min(c.tf, ref.tf), 1001)
    y, yr = c.estados(c.tf), ref.estados(ref.tf)
    vals = [
        ("tiempo", "s", abs(c.tf - ref.tf), ref.tf),
        ("temperatura", "K", abs(y[1] - yr[1]), abs(yr[1])),
        ("m", "kg", abs(y[0] - yr[0]), abs(yr[0])),
        ("P_tray", "Pa", np.max(np.abs(presion_t(c, t) - presion_t(ref, t))), ref.p.P_OBJ),
        ("T_tray", "K", np.max(np.abs(c.estados(t)[1] - ref.estados(t)[1])),
         np.max(np.abs(ref.estados(t)[1]))),
        ("m_tray", "kg", np.max(np.abs(c.estados(t)[0] - ref.estados(t)[0])), abs(yr[0])),
    ]
    return {key: float(v) for var, unit, err, scale in vals
            for key, v in ((f"error_{var}_{unit}", err), (f"error_{var}_rel", err / scale))}


def limites(vmap: dict[str, list[Caso]], datos: Path) -> None:
    eos, adia = [], []
    for ver, cases in vmap.items():
        ideal = cases[0]
        for f in (1, 0.1, 0.01, 0.001, 0):
            c = simular(ver, MODELOS[1], replace(ideal.p, A=f * ideal.p.A, B=f * ideal.p.B))
            eos.append({"version": ver, "fraccion_a_b": f, **comparar(c, ideal)})
    for ref in vmap["V2"]:
        c = simular("V3", ref.modelo, replace(ref.p, UA_GW=0, UA_WA=0))
        adia.append({"modelo": ref.modelo, **comparar(c, ref)})
    csv_file(datos / "limite_vdw.csv", eos)
    csv_file(datos / "limite_adiabatico.csv", adia)


def convergencia(v3: list[Caso], datos: Path, figuras: Path) -> None:
    filas = []
    for base in v3:
        p = base.p
        ref = simular("V3", base.modelo, replace(p, METHOD="DOP853", RTOL=2e-13, ATOL=1e-14, MAX_STEP=.0625))
        ref2 = simular("V3", base.modelo, replace(ref.p, MAX_STEP=.125))
        ensayos = [("referencia_refinada", ref2.p, ref2), ("produccion", p, base)]
        ensayos += [("paso", replace(p, RTOL=1e-6, ATOL=1e-8, MAX_STEP=h), None)
                    for h in (16., 8., 4., 2., 1., .5, .25)]
        ensayos += [("tolerancia", replace(p, METHOD=method, MAX_STEP=20, RTOL=r, ATOL=r / 100), None)
                    for method in ("RK45", "DOP853") for r in (1e-5, 1e-7, 1e-9, 1e-11)]
        for group, pp, reuse in ensayos:
            c = reuse or simular("V3", base.modelo, pp)
            filas.append({"estudio": group, "modelo": base.modelo, "metodo": pp.METHOD,
                          "rtol": pp.RTOL, "atol": pp.ATOL, "paso_max_s": pp.MAX_STEP,
                          "pasos": len(c.sol.t) - 1, "nfev": c.sol.nfev, "t_final_s": c.tf,
                          **comparar(c, ref)})
    csv_file(datos / "convergencia.csv", filas)
    for group, x, label in (("paso", "paso_max_s", "Paso máximo [s]"),
                            ("tolerancia", "rtol", "Tolerancia relativa [-]")):
        for y, ylabel, tag in (("error_tiempo_s", "Error de tiempo final [s]", "tiempo"),
                               ("error_temperatura_K", "Error de temperatura final [K]", "temperatura")):
            fig, ax = plt.subplots(figsize=(8, 4.5))
            for model in MODELOS:
                for method in ("RK45", "DOP853"):
                    rows = [r for r in filas if r["estudio"] == group and r["modelo"] == model
                            and r["metodo"] == method]
                    if rows:
                        ax.loglog([r[x] for r in rows], [max(r[y], 1e-14) for r in rows],
                                  "o-" if method == "RK45" else "s--", color=COLORS[model],
                                  label=f"{model} · {method}")
            fig_style(ax, label, ylabel, "V3 frente a la referencia DOP853 refinada")
            export_fig(fig, figuras / f"convergencia_{group}_{tag}.png")


def validar_evento(casos: list[Caso], datos: Path) -> None:
    filas = []
    for c in casos:
        alt = simular(c.version, c.modelo, c.p, horizonte=c.tf + 12)
        for h in (10, 5, 1, .25):
            t = np.arange(0, alt.tf, h)
            q = presion_t(alt, t) - c.p.P_OBJ
            k = int(np.flatnonzero((q[:-1] < 0) & (q[1:] >= 0))[0])
            t_linear = t[k] - q[k] * (t[k + 1] - t[k]) / (q[k + 1] - q[k])
            for metodo, tt in (("umbral", t[k + 1]), ("lineal", t_linear)):
                dt = abs(tt - c.tf)
                dp = abs(float(presion_t(alt, tt)) - c.p.P_OBJ)
                filas.append({"version": c.version, "modelo": c.modelo, "metodo": metodo,
                              "malla_s": h, "tf_evento_s": c.tf, "tf_estimado_s": tt,
                              "error_t_s": dt, "error_t_rel": dt / c.tf,
                              "error_P_Pa": dp, "error_P_rel": dp / c.p.P_OBJ})
    csv_file(datos / "evento.csv", filas)


def enriquecer_resumen(c: Caso) -> dict:
    r = resumen_caso(c)
    t_peak, peak = pico_temperatura(c)
    for name in ("Tgas_final", "Tpared_final"):
        x = r.pop(f"{name}_C")
        r[f"{name}_K"] = x + 273.15 if x != "" else None
    r.pop("Tgas_max_C")
    r.update(t_Tmax_s=t_peak, Tmax_K=peak,
             caudal_medio_kg_s=(r["masa_final_kg"] - c.m0) / c.tf,
             m0_kg=c.m0, h_s_J_kg=c.h_in, area_m2=c.p.AREA)
    return r


def barridos(datos: Path, figuras: Path) -> None:
    casos = (
        ("UA_gas_pared", "Conductancia gas–pared [W/K]", "UA_GW", (0, 4, 8, 16), BASE),
        ("UA_pared_amb", "Conductancia pared–ambiente [W/K]", "UA_WA", (1.5, 3, 6), BASE),
        ("C_pared", "Capacidad térmica de pared [J/K]", "CW", (2000, 4000, 8000), BASE),
        ("P_suministro", "Presión de suministro [bar]", "PS", (150, 200, 250),
         replace(BASE, P_OBJ=140 * BAR)),
        ("T_suministro", "Temperatura de suministro [K]", "TS", (280, 300, 320), BASE),
    )
    tablas, diferencias = [], []
    for tag, xlabel, var, values, base in casos:
        rows = []
        for value in values:
            kw = {"MW": value / base.CPW} if var == "CW" else {var: value * BAR if var == "PS" else value}
            p = replace(base, **kw)
            if p.AREA != BASE.AREA or p.CD != BASE.CD:
                raise RuntimeError("El área y Cd deben conservarse en los barridos.")
            for mod, short in zip(MODELOS, ("ideal", "vdw")):
                c = simular("V3", mod, p)
                rows.append({"parametro": tag, "valor": value, "Pobjetivo_bar": p.P_OBJ / BAR,
                             **enriquecer_resumen(c)})
                slug = str(value).replace(".", "p")
                save_samples(c, datos / f"barrido_{tag}_{slug}_{short}.csv", n=401)
            one, two = rows[-2:]
            indicadores = ("tiempo_final_s", "Tmax_K", "Tgas_final_K", "masa_final_kg", "caudal_medio_kg_s",
                           "mdot_inicial_kg_s", "caudal_final_kg_s",
                           "Qgas_pared_acum_kJ", "Qpared_ambiente_acum_kJ")
            diferencias.append({"parametro": tag, "valor": value,
                                **{f"delta_{k}": two[k] - one[k] for k in indicadores}})
        tablas.extend(rows)
        for y, ylabel, slug in (
            ("tiempo_final_s", "Tiempo de llenado [s]", "tiempo"),
            ("Tmax_K", "Temperatura máxima [K]", "temperatura"),
            ("masa_final_kg", "Masa final [kg]", "masa"),
            ("Qgas_pared_acum_kJ", "Calor total gas–pared [kJ]", "calor"),
        ):
            grafico_barrido(rows, "valor", y, ylabel, xlabel,
                           f"Sensibilidad de V3 a {xlabel.split(' [')[0].lower()}",
                           figuras / f"{tag}_{slug}.png")
    csv_file(datos / "sensibilidad.csv", tablas)
    csv_file(datos / "sensibilidad_diferencias_vdw_menos_ideal.csv", diferencias)


def figuras_base(vmap, datos: Path, figuras: Path) -> None:
    magnitudes = {
        "temperatura": "Temperatura del gas [°C]", "presion": "Presión [bar]",
        "masa": "Masa [kg]", "caudal": "Caudal másico [kg/s]", "factor_z": "Factor Z [-]",
        "pared": "Temperatura de la pared [°C]", "calor": "Potencia gas–pared [W]",
        "calor_ambiente": "Potencia pared–ambiente [W]",
    }
    for ver, cc in vmap.items():
        for c in cc:
            save_samples(c, datos / f"{ver}_{etiqueta(c)}.csv")
        for key, ylabel in magnitudes.items():
            if ver == "V3" or key not in ("pared", "calor", "calor_ambiente"):
                figura_version(cc, key, ylabel, figuras / f"{ver}_{key}.png")
    comparar_v2_v3(vmap["V2"], vmap["V3"], figuras / "V2_V3_comparacion_T.png")
    fig, ax = plt.subplots(figsize=(8, 4.55))
    for c in vmap["V3"]:
        d = c.muestras()
        for key, label, style in (("Tg", "Gas", "-"), ("Tw", "Pared", "--")):
            ax.plot(d["t"], d[key], style, color=COLORS[c.modelo], label=f"{label} · {c.modelo}")
    fig_style(ax, "Tiempo [s]", "Temperatura [K]", "V3: intercambio térmico gas–pared")
    export_fig(fig, figuras / "V3_gas_y_pared.png")
    fig, ax = plt.subplots(figsize=(8, 4.55))
    for ver, style in zip(VERSIONES, (":", "--", "-")):
        for c in vmap[ver]:
            d = c.muestras()
            ax.plot(d["t"], d["P"] / BAR, style, color=COLORS[c.modelo], label=f"{ver} · {c.modelo}")
    fig_style(ax, "Tiempo [s]", "Presión [bar]", "Comparación de V1, V2 y V3")
    export_fig(fig, figuras / "V1_V2_V3_presion.png")


def figura_balances(filas: list[dict], figuras: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 4.55))
    x = np.arange(len(filas))
    for key, label, offset in (("masa_rel", "Masa", -.24), ("gas_rel", "Gas", -.08),
                               ("pared_rel", "Pared", .08), ("global_rel", "Gas + pared", .24)):
        r = [f[key] if f[key] is not None and f[key] > 0 else np.nan for f in filas]
        ax.scatter(x + offset, r, label=label)
    ax.set_xticks(x, [f"{f['version']}·{'I' if f['modelo']=='Ideal' else 'VdW'}" for f in filas])
    ax.set_yscale("log")
    fig_style(ax, "Modelo", "Residuo máximo normalizado [-]", "Cierre integral con cuadratura posterior")
    export_fig(fig, figuras / "validacion_balances.png")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--salida", type=Path, default=Path(__file__).resolve().parent)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--sin-video", action="store_true")
    modes.add_argument("--solo-video", action="store_true")
    parser.add_argument("--frames-por-version", type=int, default=64)
    parser.add_argument("--fps", type=int, default=16)
    args = parser.parse_args()
    if args.frames_por_version < 10 or args.fps < 1:
        parser.error("Se requiere un mínimo de 10 fotogramas por versión y 1 fps.")
    root = args.salida.resolve()
    datos, figuras = root / "datos", root / "figuras"
    vmap = {ver: [simular(ver, mod, replace(BASE, P_OBJ=(100 if ver == "V1" else 190) * BAR))
                  for mod in MODELOS] for ver in VERSIONES}
    all_cases = [c for pair in vmap.values() for c in pair]
    if not args.solo_video:
        for folder in (datos, figuras):
            folder.mkdir(parents=True, exist_ok=True)
        csv_file(datos / "resumen.csv", [enriquecer_resumen(c) for c in all_cases])
        balances = [validacion_balances(c, datos) for c in all_cases]
        csv_file(datos / "balances.csv", balances)
        csv_file(datos / "analitica.csv", [validar_ideal(vmap[ver][0]) for ver in ("V1", "V2")])
        figura_balances(balances, figuras)
        limites(vmap, datos)
        validar_evento(all_cases, datos)
        convergencia(vmap["V3"], datos, figuras)
        prueba_rk4_orden(BASE, datos, figuras)
        barridos(datos, figuras)
        figuras_base(vmap, datos, figuras)
        meta = {"python": sys.version.split()[0], "numpy": np.__version__, "scipy": scipy.__version__,
                "matplotlib": matplotlib.__version__, "R_J_kgK": R_N2, "cp_J_kgK": CP_N2,
                "cv_J_kgK": CV_N2, "gamma": GAMMA, "a_SI": A_VDW, "b_SI": B_VDW,
                "area_m2": BASE.AREA, "CdA_m2": BASE.CD * BASE.AREA,
                "r_crit": relacion_critica(BASE.GAMMA), "rtol": BASE.RTOL, "atol": BASE.ATOL,
                "paso_max_s": BASE.MAX_STEP, "casos_base": len(all_cases)}
        (datos / "metadatos.json").write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")
        print(f"Cálculos, datos y figuras generados en {root}", flush=True)
    if not args.sin_video:
        crear_mp4(vmap, root / "animaciones" / "ICF130_V1_V2_V3_FINAL.mp4",
                   frames_por_version=args.frames_por_version, fps=args.fps)


if __name__ == "__main__":
    main()

