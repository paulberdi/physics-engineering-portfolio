"""ICF130. Ampliación Peng-Robinson y comparación con la EOS de Span.

Copia de ICF130_AMPLIACION.py. Único cambio: importa el núcleo desde modelos/base.py.
"""
from dataclasses import replace
from pathlib import Path
import argparse
import json
import math

import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq, minimize_scalar
from CoolProp import AbstractState, CoolProp as CP, __version__ as coolprop_version
from .base import (
    BASE, MODELOS, Caso, TC_N2, PC_N2, caudal, csv_file, energia_gas,
    entalpia_entrada, fig_style, export_fig, presion_eos, regimen, simular,
    validar_parametros, relacion_critica, volumen_vdw, plt,
)

PR = "Peng–Robinson"
MODELS = (*MODELOS, PR)
TAGS = dict(zip(MODELS, ("ideal", "vdw", "pr")))
COLORS = dict(zip(MODELS, ("#05739e", "#a14e0b", "#007761")))
BAR = 1e5
OMEGA = 0.0372
KAPPA = 0.37464 + 1.54226 * OMEGA - 0.26992 * OMEGA**2
SQ2 = math.sqrt(2)
B = 0.07780 * BASE.R * TC_N2 / PC_N2
AC = 0.45724 * BASE.R**2 * TC_N2**2 / PC_N2


def atraccion(T, f=1.0):
    """a(T), a'(T) y a''(T), en base másica."""
    ac = f * AC
    raiz = np.sqrt(T / TC_N2)
    g = 1 + KAPPA * (1 - raiz)
    dg = -KAPPA / (2 * TC_N2 * raiz)
    ddg = -dg / (2 * T)
    return ac * g**2, 2 * ac * g * dg, 2 * ac * (dg**2 + g * ddg)


def propiedades(T, v, p=BASE, f=1.0):
    """Devuelve P, u, h, cv efectivo, (du/dv)T y (dP/dv)T."""
    b = f * B
    if np.any(np.asarray(T) <= TC_N2) or np.any(np.asarray(v) <= b):
        raise ValueError("Peng–Robinson exige T>Tc y v>b en esta ampliación.")
    a, da, dda = atraccion(T, f)
    D = v**2 + 2 * b * v - b**2
    I = np.log1p(2 * SQ2 * b / (v + (1 - SQ2) * b)) / (2 * SQ2 * b) if b else 1 / v
    P = p.R * T / (v - b) - a / D
    u = p.CV * T + (T * da - a) * I
    cv = p.CV + T * dda * I
    uv = (a - T * da) / D
    pv = -p.R * T / (v - b)**2 + 2 * a * (v + b) / D**2
    return P, u, u + P * v, cv, uv, pv


def volumen_pr(P, T, p=BASE, f=1.0):
    if f == 0:
        return p.R * T / P
    b = f * B
    a = atraccion(T, f)[0]

    def residuo(rho):
        return p.R * T * rho / (1 - b * rho) - a * rho**2 / (
            1 + 2 * b * rho - (b * rho)**2) - P

    rho = brentq(residuo, 0, (1 - 1e-10) / b, xtol=1e-12, rtol=1e-14)
    v = 1 / rho
    if propiedades(T, v, p, f)[-1] >= 0:
        raise ValueError("Raíz mecánicamente inestable.")
    return v


class CasoPR(Caso):
    def muestras_en(self, t):
        m, T, Tw = self.estados(t)[:3]
        P = propiedades(T, self.p.V / m, self.p)[0]
        md = np.array([caudal(x, self.p, self.version) for x in P])
        qgw, qwa = tasas(T, Tw, self.p, self.version)
        return dict(t=t, m=m, Tg=T, Tw=Tw, P=P, mdot=md, q_gw=qgw, q_wa=qwa,
                    Z=P * self.p.V / (m * self.p.R * T))


def tasas(T, Tw, p, ver):
    if ver == "V3":
        return p.UA_GW * (T - Tw), p.UA_WA * (Tw - p.TINF)
    return np.zeros_like(T), np.zeros_like(T)


def simular_pr(ver, p, *, conservativo=False, h_impuesta=None):
    validar_parametros(p, ver)
    v0 = volumen_pr(p.P0, p.T0, p)
    hs = propiedades(p.TS, volumen_pr(p.PS, p.TS, p), p)[2]
    if h_impuesta is not None:
        hs = h_impuesta
    m0 = p.V / v0
    E0 = m0 * propiedades(p.T0, v0, p)[1]

    def temperatura(m, U):
        return brentq(lambda T: m * propiedades(T, p.V / m, p)[1] - U,
                      TC_N2 + 1, 1200, xtol=1e-10)

    def rhs(_t, y):
        m, X, Tw = y
        T = temperatura(m, X) if conservativo else X
        v = p.V / m
        P, u, _, cv, uv, pv = propiedades(T, v, p)
        if cv <= 0 or pv >= 0:
            raise ValueError("Estado fuera de la rama estable.")
        md = caudal(P, p, ver)
        qgw, qwa = tasas(T, Tw, p, ver)
        Udot = md * hs - qgw
        Xdot = Udot if conservativo else (Udot - md * (u - v * uv)) / (m * cv)
        return [md, Xdot, (qgw - qwa) / p.CW]

    def evento(_t, y):
        m, X, _ = y
        T = temperatura(m, X) if conservativo else X
        return propiedades(T, p.V / m, p)[0] - p.P_OBJ

    evento.terminal = True
    evento.direction = 1
    sol = solve_ivp(rhs, (0, p.TIME_LIMIT), [m0, E0 if conservativo else p.T0, p.TW0],
                    events=evento, dense_output=True, method="DOP853" if conservativo else p.METHOD,
                    rtol=2e-12 if conservativo else p.RTOL, atol=p.ATOL, max_step=p.MAX_STEP)
    if not sol.success or not sol.t_events[0].size:
        raise RuntimeError(f"{ver}: no se alcanzó el evento. {sol.message}")
    if conservativo:
        original = sol.sol

        def estados(t):
            y = original(t)
            if y.ndim == 1:
                y[1] = temperatura(y[0], y[1])
            else:
                y[1] = [temperatura(m, U) for m, U in zip(y[0], y[1])]
            return y

        sol.sol = estados
    return CasoPR(ver, PR, p, m0, hs, sol, float(sol.t_events[0][0]))


def sim(ver, modelo, p):
    return simular_pr(ver, p) if modelo == PR else simular(ver, modelo, p)


def resumir(c):
    dat = c.muestras(801)
    i = int(np.argmax(dat["Tg"]))
    t = dat["t"]
    if 0 < i < len(t) - 1:
        peak = minimize_scalar(lambda x: -c.estados(x)[1],
                               bounds=(t[i-1], t[i+1]), method="bounded")
        t_peak, T_peak = peak.x, -peak.fun
    else:
        t_peak, T_peak = t[i], dat["Tg"][i]
    fin = c.estados(c.tf)
    p = c.p
    if c.modelo == PR:
        integrales = cuadratura(c)[1]
        qgw, qwa = integrales[-1, 2:]
    else:
        qgw, qwa = fin[4:6]
    return dict(version=c.version, modelo=c.modelo, Ps_bar=p.PS/BAR, Ts_K=p.TS,
                objetivo_bar=p.P_OBJ/BAR, CdA_m2=p.CD*p.AREA, tf_s=c.tf, Tf_K=fin[1],
                Tmax_K=T_peak, t_Tmax_s=t_peak, mf_kg=fin[0], Zf=dat["Z"][-1],
                Twf_K=fin[2] if c.version == "V3" else None, hs_J_kg=c.h_in,
                m0_kg=c.m0, mdot_0_g_s=1000*dat["mdot"][0],
                mdot_f_g_s=1000*dat["mdot"][-1], mdot_medio_g_s=1000*(fin[0]-c.m0)/c.tf,
                Qgw_kJ=qgw/1000, Qwa_kJ=qwa/1000)


def cuadratura(c):
    t = c.sol.t
    pc = c.p.PS * relacion_critica(c.p.GAMMA)
    if c.version != "V1" and c.p.P0 < pc < c.p.P_OBJ:
        tp = brentq(lambda x: propiedades(c.estados(x)[1], c.p.V/c.estados(x)[0], c.p)[0]-pc,
                    0, c.tf)
        t = np.unique(np.r_[t, tp])
    x, w = np.polynomial.legendre.leggauss(8)
    half = np.diff(t)/2
    d = c.muestras_en((t[:-1,None]+half[:,None]*(1+x)).ravel())
    rates = np.array([d["mdot"], d["mdot"]*c.h_in, d["q_gw"], d["q_wa"]])
    increments = (rates.reshape(4,-1,8) @ w) * half
    return t, np.vstack((np.zeros(4), np.cumsum(increments, axis=1).T))


def balances(c):
    t, I = cuadratura(c)
    m,T,Tw = c.estados(t)[:3]
    dm, Ein, Qg, Qa = I.T
    U = m * propiedades(T,c.p.V/m,c.p)[1]
    du = U-U[0]
    dw = c.p.CW*(Tw-Tw[0])
    errors = [m-c.m0-dm, du-Ein+Qg, dw-Qg+Qa, du+dw-Ein+Qa]
    scales = [dm[-1], abs(Ein[-1])+abs(Qg[-1]), abs(Qg[-1])+abs(Qa[-1]),
              abs(Ein[-1])+abs(Qa[-1])]
    row = dict(version=c.version, Ps_bar=c.p.PS/BAR, Ts_K=c.p.TS, objetivo_bar=c.p.P_OBJ/BAR)
    for label, error, scale in zip(("masa_kg","gas_J","pared_J","global_J"),errors,scales):
        maximum = np.max(np.abs(error))
        row[label] = maximum
        row[label+"_escala"] = scale
        row[label+"_rel"] = maximum/scale if scale else None
    return row


def referencia(estados):
    ref = AbstractState("HEOS", "Nitrogen")
    assert abs(ref.acentric_factor()-OMEGA) < 1e-10
    rows = []
    for origen, T, P in estados:
        ref.update(CP.PT_INPUTS, P, T)
        rho_ref = ref.rhomass()
        R_ref = ref.gas_constant()/ref.molar_mass()
        u_ref = ref.tau()*ref.dalphar_dTau()*R_ref*T
        h_ref = u_ref + P/rho_ref - R_ref*T
        for model in MODELS:
            v = (BASE.R*T/P if model == "Ideal" else volumen_vdw(P,T,BASE)
                 if model == "Van der Waals" else volumen_pr(P,T))
            if model == PR:
                u_r = propiedades(T,v)[1]-BASE.CV*T
                pressure_at_ref = propiedades(T,1/rho_ref)[0]
            else:
                u_r = energia_gas(model,BASE.V/v,T,BASE)/(BASE.V/v)-BASE.CV*T
                pressure_at_ref = presion_eos(model,BASE.V*rho_ref,T,BASE)
            h_r = u_r + P*v - BASE.R*T
            rows.append(dict(origen=origen,T_K=T,P_bar=P/BAR,modelo=model,
                             rho_ref_kg_m3=rho_ref,rho_kg_m3=1/v,
                             error_rho_pct=100*(1/v/rho_ref-1),
                             Z_ref=P/(rho_ref*R_ref*T),Z=P*v/(BASE.R*T),
                             error_P_a_rho_ref_pct=100*(pressure_at_ref/P-1),
                             u_res_ref_J_kg=u_ref,u_res_J_kg=u_r,error_u_J_kg=u_r-u_ref,
                             h_res_ref_J_kg=h_ref,h_res_J_kg=h_r,error_h_J_kg=h_r-h_ref))
    return rows


def verificar_propiedades():
    rows = []
    for T in (280., 300., 340., 420.):
        for Pbar in (1, 50, 150, 200, 250):
            P = Pbar*BAR
            v = volumen_pr(P,T)
            pp,u,h,cv,uv,pv = propiedades(T,v)
            dv, dt = 1e-5*v, 1e-2
            up = (propiedades(T,v+dv)[1]-propiedades(T,v-dv)[1])/(2*dv)
            cpv = (propiedades(T+dt,v)[1]-propiedades(T-dt,v)[1])/(2*dt)
            # Integración directa de la identidad de Maxwell en densidad.
            # (du/dv)_T = (a-Ta')/D, dv = -d rho/rho².
            a, da, _ = atraccion(T)
            residual = quad(lambda rho: -(a-T*da)/(1+2*B*rho-(B*rho)**2),0,1/v)[0]
            row = dict(T_K=T,P_bar=Pbar,error_P_Pa=abs(pp-P),error_h_J_kg=abs(h-u-pp*v),
                       error_u_cuadratura_J_kg=abs(u-BASE.CV*T-residual),
                       error_cv_J_kgK=abs(cpv-cv),error_dudv_rel=abs(up/uv-1),
                       cv_J_kgK=cv,dPdv=pv)
            assert row["error_P_Pa"] < .01 and row["error_h_J_kg"] < 1e-8
            assert row["error_u_cuadratura_J_kg"] < 1e-6 and row["error_cv_J_kgK"] < 1e-5
            assert row["error_dudv_rel"] < 1e-6 and cv > 0 and pv < 0
            rows.append(row)
    return rows


def limite_ideal():
    rows = []
    for factor in (1, .1, .01, .001, 0):
        v = volumen_pr(200*BAR,300,f=factor)
        _,u,h,cv,*_ = propiedades(300,v,f=factor)
        rows.append(dict(fraccion_a_b=factor,error_Z=abs(200*BAR*v/(BASE.R*300)-1),
                         error_u_J_kg=abs(u-BASE.CV*300),error_h_J_kg=abs(h-BASE.CP*300),
                         error_cv_J_kgK=abs(cv-BASE.CV)))
    assert max(v for k,v in rows[-1].items() if k != "fraccion_a_b") < 1e-8
    return rows


def contrastar_integradores(c):
    alt = simular_pr(c.version,replace(c.p,MAX_STEP=.25),conservativo=True)
    t = np.linspace(0,min(c.tf,alt.tf),401)
    ref = alt.estados(t)
    err = np.max(np.abs(c.estados(t)-ref),axis=1)
    row = dict(version=c.version,error_tf_s=abs(c.tf-alt.tf),error_m_kg=err[0],
               error_T_K=err[1],error_Tw_K=err[2],
               error_tf_rel=abs(c.tf-alt.tf)/alt.tf,error_T_rel=err[1]/max(ref[1]))
    assert row["error_tf_s"] < 1e-4 and row["error_T_K"] < 2e-5
    return row


def trayectoria(c, out):
    d = c.muestras(1001)
    cols = dict(t_s=d["t"],P_bar=d["P"]/BAR,T_K=d["Tg"],masa_kg=d["m"],
                Z=d["Z"],mdot_g_s=1000*d["mdot"],qgw_W=d["q_gw"],qwa_W=d["q_wa"])
    if c.version == "V3":
        cols["Tw_K"] = d["Tw"]
    csv_file(out,[dict(zip(cols,x)) for x in zip(*cols.values())])


def graficos_base(base, out):
    for version, cases in base.items():
        for key, ylabel, unit in (("P","Presión [bar]",1/BAR),("Tg","Temperatura [K]",1),
                                  ("m","Masa [kg]",1),("Z","Factor Z [-]",1),
                                  ("mdot","Caudal [g/s]",1000)):
            fig, ax = plt.subplots(figsize=(7.6,4.2))
            for c in cases:
                d = c.muestras(801)
                ax.plot(d["t"],d[key]*unit,color=COLORS[c.modelo],label=c.modelo)
            fig_style(ax,"Tiempo [s]",ylabel,f"{version} · {ylabel.split(' [')[0]}")
            export_fig(fig,out/f"{version}_{key}.png")


def graficos_barrido(delta, ps, ts, out):
    fig, axes = plt.subplots(2,3,figsize=(13,7),layout="constrained")
    for i,(a,b) in enumerate(((PR,"Ideal"),(PR,"Van der Waals"))):
        chosen = [r for r in delta if r["modelo"]==a and r["referencia"]==b]
        for j,(key,title) in enumerate((("tf_s","Tiempo [s]"),("Tf_K","Temperatura final [K]"),
                                       ("mf_kg","Masa final [kg]"))):
            ax = axes[i,j]
            grid = np.array([r["delta_"+key] for r in chosen]).reshape(len(ps),len(ts))
            vmax=max(abs(grid.min()),abs(grid.max()))
            im=ax.imshow(grid,cmap="RdBu_r",vmin=-vmax,vmax=vmax,aspect="auto")
            for k,l in np.ndindex(grid.shape):
                val=grid[k,l]
                ax.text(l,k,f"{val:+.2f}",ha="center",va="center",fontsize=8,
                        color="white" if abs(val)>.62*vmax else "#122535")
            ax.set(xticks=range(len(ts)), xticklabels=ts, yticks=range(len(ps)), yticklabels=ps,
                   xlabel="T suministro [K]",ylabel="P suministro [bar]",title=f"PR − {TAGS[b]} · {title}")
            fig.colorbar(im,ax=ax,shrink=.84)
    export_fig(fig,out/"diferencias_V3_140bar.png")


def grafico_ref(rows, out):
    for key, ylabel in (("error_rho_pct","Error de densidad [%]"),("error_h_J_kg","Error en h residual [J/kg]"),
                        ("error_u_J_kg","Error en u residual [J/kg]")):
        fig, ax=plt.subplots(figsize=(7.6,4.2))
        for name in MODELS:
            r=[x for x in rows if x["origen"]=="malla" and x["T_K"]==300 and x["modelo"]==name]
            ax.plot([x["P_bar"] for x in r],[x[key] for x in r],"o-",color=COLORS[name],label=name)
        ax.axhline(0,color="#666",lw=.7)
        fig_style(ax,"Presión [bar]",ylabel,"Nitrógeno a 300 K · contraste con Span (HEOS)")
        export_fig(fig,out/f"HEOS_{key}.png")


def comparacion_isoterma(casos):
    """Límite T=T0 en V2, con el mismo caudal, vía cuadratura de dt/dP."""
    rows = []
    for c in casos:
        p, model = c.p, c.modelo
        def volumen(P):
            if model == PR:
                return volumen_pr(P,p.T0,p)
            if model == "Van der Waals":
                return volumen_vdw(P,p.T0,p)
            return p.R*p.T0/P

        def dt_dP(P):
            v = volumen(P)
            if model == PR:
                pv = propiedades(p.T0,v,p)[-1]
            elif model == "Van der Waals":
                pv = -p.R*p.T0/(v-p.B)**2 + 2*p.A/v**3
            else:
                pv = -p.R*p.T0/v**2
            return -p.V/(v*v*pv*caudal(P,p,"V2"))

        pc = p.PS*relacion_critica(p.GAMMA)
        tf = quad(dt_dP,p.P0,p.P_OBJ,points=[pc],epsabs=1e-8)[0]
        mf = p.V/volumen(p.P_OBJ)
        def U(m):
            if model == PR:
                return m*propiedades(p.T0,p.V/m,p)[1]
            return energia_gas(model,m,p.T0,p)
        Q = (mf-c.m0)*c.h_in-(U(mf)-U(c.m0))
        if model == "Ideal":
            assert abs(Q-(mf-c.m0)*p.R*p.T0)<1e-8
        rows.append(dict(modelo=model,tf_isotermo_s=tf,tf_adiabatico_s=c.tf,
                         mf_isotermo_kg=mf,mf_adiabatico_kg=c.estados(c.tf)[0],
                         T_isoterma_K=p.T0,Tf_adiabatico_K=c.estados(c.tf)[1],
                         calor_extraido_kJ=Q/1000))
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--salida",type=Path,default=Path(__file__).resolve().parent)
    out=parser.parse_args().salida.resolve()
    ejecutar(out/"datos",out/"figuras")


def ejecutar(dat, fig):
    """Cuerpo original de la ampliación, con carpetas de datos y figuras explícitas."""
    for folder in (dat,fig):
        folder.mkdir(parents=True,exist_ok=True)
    csv_file(dat/"validacion_propiedades.csv",verificar_propiedades())
    csv_file(dat/"limite_PR_ideal.csv",limite_ideal())
    base={ver:[sim(ver,name,replace(BASE,P_OBJ=(100 if ver=="V1" else 190)*BAR))
               for name in MODELS] for ver in ("V1","V2","V3")}
    base_rows=[resumir(c) for cases in base.values() for c in cases]
    csv_file(dat/"casos_base.csv",base_rows)
    for cases in base.values():
        for c in cases:
            trayectoria(c,dat/f"{c.version}_{TAGS[c.modelo]}.csv")
    csv_file(dat/"V2_limite_isotermo.csv",comparacion_isoterma(base["V2"]))
    pr_base=[cases[-1] for cases in base.values()]
    csv_file(dat/"validacion_formulaciones.csv",[contrastar_integradores(c) for c in pr_base])
    v2=pr_base[1]
    off=simular_pr("V3",replace(BASE,UA_GW=0,UA_WA=0))
    t=np.linspace(0,min(off.tf,v2.tf),401)
    err=max(np.abs(off.estados(t)[1]-v2.estados(t)[1]))
    assert abs(off.tf-v2.tf)<1e-8 and err<1e-8
    csv_file(dat/"limite_PR_adiabatico.csv",[dict(error_tf_s=abs(off.tf-v2.tf),error_T_K=err)])
    ps = (150,175,200,225,250)
    ts = (280,290,300,310,320)
    fine, delta, decomp = [], [], []
    pr_cases = pr_base.copy()
    for P in ps:
        for T in ts:
            p=replace(BASE,PS=P*BAR,TS=T,P_OBJ=140*BAR)
            cases=[sim("V3",name,p) for name in MODELS]
            pr_cases.append(cases[-1])
            rows=[resumir(c) for c in cases]
            fine.extend(rows)
            for hi,lo in ((2,0),(2,1),(1,0)):
                x,y=rows[hi],rows[lo]
                deltas={}
                for k in ("tf_s","Tf_K","Tmax_K","mf_kg","Zf","hs_J_kg","Qgw_kJ","Qwa_kJ","mdot_medio_g_s"):
                    value=x[k]-y[k]
                    deltas.update({"delta_"+k:value,"delta_pct_"+k:100*value/y[k]})
                delta.append(dict(Ps_bar=P,Ts_K=T,modelo=x["modelo"],referencia=y["modelo"],**deltas))
            # Cambio secuencial: EOS del depósito a hs fija y luego hs de PR.
            if P in (150,200,250) and T in (280,300,320):
                mid=resumir(simular_pr("V3",p,h_impuesta=cases[0].h_in))
                for k in ("tf_s","Tf_K","mf_kg"):
                    decomp.append(dict(Ps_bar=P,Ts_K=T,magnitud=k,
                                       delta_deposito=mid[k]-rows[0][k],
                                       delta_entalpia=rows[2][k]-mid[k],
                                       delta_total=rows[2][k]-rows[0][k]))
    csv_file(dat/"V3_malla_suministro.csv",fine)
    csv_file(dat/"V3_diferencias.csv",delta)
    csv_file(dat/"V3_descomposicion_secuencial.csv",decomp)
    checks=[balances(c) for c in pr_cases]
    assert max(r["global_J_rel"] for r in checks)<1e-7
    assert max(r["masa_kg_rel"] for r in checks)<1e-7
    csv_file(dat/"validacion_balances_PR.csv",checks)
    fields=[("malla",T,P*BAR) for T in (280,300,320,360,420) for P in (1,50,100,150,190,200,250)]
    fields += [(r["version"]+"_final_"+TAGS[r["modelo"]],r["Tf_K"],r["objetivo_bar"]*BAR)
               for r in base_rows]
    refrows=referencia(fields)
    csv_file(dat/"propiedades_contra_HEOS.csv",refrows)
    graficos_base(base,fig)
    graficos_barrido(delta,ps,ts,fig)
    grafico_ref(refrows,fig)
    (dat/"metadatos.json").write_text(json.dumps(dict(CoolProp=coolprop_version, referencia="HEOS::Nitrogen (Span et al.)",
                                                     Tc_K=TC_N2,Pc_Pa=PC_N2,omega=OMEGA,a_c=AC,b=B,
                                                     CdA_m2=BASE.CD*BASE.AREA),indent=2)+"\n",encoding="utf-8")
    print(f"Verificados {len(base_rows)} casos base, {len(fine)} casos de malla, {len(checks)} balances PR y {len(refrows)} comparaciones HEOS.")
    return base


if __name__=="__main__":
    main()

