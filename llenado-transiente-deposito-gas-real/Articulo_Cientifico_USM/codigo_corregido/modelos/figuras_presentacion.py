"""Figuras en formato de diapositiva (16:9, letra grande) a partir de los CSV generados.

Solo lee resultados ya calculados en esta misma ejecución; no inventa curvas.
"""
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

COL = {"Ideal": "#2563EB", "Van der Waals": "#EA580C", "Peng–Robinson": "#059669"}
ETQ = {"Ideal": "Gas ideal", "Van der Waals": "Van der Waals", "Peng–Robinson": "Peng-Robinson"}
TAG = {"ideal": "Ideal", "vdw": "Van der Waals", "pr": "Peng–Robinson"}
plt.rcParams.update({
    "font.size": 15, "axes.titlesize": 17, "axes.labelsize": 15, "legend.fontsize": 13,
    "xtick.labelsize": 13, "ytick.labelsize": 13, "axes.spines.top": False,
    "axes.spines.right": False, "axes.grid": True, "grid.alpha": 0.25, "figure.dpi": 100,
    "axes.titleweight": "bold", "axes.titlelocation": "left",
})


def leer(path: Path) -> list[dict]:
    with path.open(encoding="utf-8-sig") as f:
        filas = list(csv.DictReader(f))
    for r in filas:
        for k, v in r.items():
            try:
                r[k] = float(v)
            except (TypeError, ValueError):
                pass
    return filas


def serie(datos: Path, ver: str, tag: str) -> dict[str, np.ndarray]:
    filas = leer(datos / "peng_robinson" / f"{ver}_{tag}.csv")
    return {k: np.array([r[k] for r in filas], dtype=float) for k in filas[0]
            if isinstance(filas[0][k], float)}


def guardar(fig, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def generar(datos: Path, out: Path) -> None:
    S = {(v, t): serie(datos, v, t) for v in ("V1", "V2", "V3") for t in TAG}

    # 1. V3: temperatura del gas y de la pared
    fig, ax = plt.subplots(figsize=(10, 5.6))
    for t, m in TAG.items():
        d = S[("V3", t)]
        ax.plot(d["t_s"], d["T_K"] - 273.15, color=COL[m], lw=2.8, label=f"{ETQ[m]} · gas")
        ax.plot(d["t_s"], d["Tw_K"] - 273.15, "--", color=COL[m], lw=2.0, label=f"{ETQ[m]} · pared")
    ax.set(xlabel="Tiempo [s]", ylabel="Temperatura [°C]", title="V3: el gas se calienta y la pared lo sigue con retardo")
    ax.legend(ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.16), frameon=False, fontsize=12.5)
    guardar(fig, out / "f01_V3_temperaturas.png")

    # 2. V3: presión y masa
    fig, axs = plt.subplots(1, 2, figsize=(12.5, 5.0))
    for t, m in TAG.items():
        d = S[("V3", t)]
        axs[0].plot(d["t_s"], d["P_bar"], color=COL[m], lw=2.6, label=ETQ[m])
        axs[1].plot(d["t_s"], d["masa_kg"], color=COL[m], lw=2.6, label=ETQ[m])
    axs[0].axhline(190, color="#DC2626", ls=":", lw=1.4)
    axs[0].set(xlabel="Tiempo [s]", ylabel="Presión [bar]", title="Presión")
    axs[1].set(xlabel="Tiempo [s]", ylabel="Masa [kg]", title="Masa almacenada")
    axs[1].legend(frameon=False, loc="lower right")
    fig.tight_layout()
    guardar(fig, out / "f02_V3_presion_masa.png")

    # 3. Caudal en V2 y cambio de régimen
    fig, ax = plt.subplots(figsize=(10, 5.2))
    for t, m in TAG.items():
        d = S[("V2", t)]
        ax.plot(d["t_s"], d["mdot_g_s"], color=COL[m], lw=2.6, label=ETQ[m])
    d = S[("V2", "ideal")]
    rcrit = 0.528299
    k = int(np.argmax(d["P_bar"] > rcrit * 200))
    ax.axvspan(0, d["t_s"][k], color="#94A3B8", alpha=0.15)
    ax.text(d["t_s"][k] * 0.5, 10.5, "estrangulado\n(Mach 1 en la garganta)", ha="center", fontsize=13)
    ax.text(d["t_s"][k] + (d["t_s"][-1] - d["t_s"][k]) * 0.5, 13.5, "subsónico", ha="center", fontsize=13)
    ax.set(xlabel="Tiempo [s]", ylabel="Caudal [g/s]", ylim=(0, 22.5),
           title=f"V2: el caudal cae cuando P supera r*·Ps ≈ {rcrit * 200:.1f} bar")
    ax.legend(frameon=False, loc="lower left")
    guardar(fig, out / "f03_V2_caudal.png")

    # 4. V1: temperatura y presión
    fig, axs = plt.subplots(1, 2, figsize=(12.5, 5.0))
    for t, m in TAG.items():
        d = S[("V1", t)]
        axs[0].plot(d["t_s"], d["T_K"] - 273.15, color=COL[m], lw=2.6, label=ETQ[m])
        axs[1].plot(d["t_s"], d["P_bar"], color=COL[m], lw=2.6, label=ETQ[m])
    axs[0].set(xlabel="Tiempo [s]", ylabel="Temperatura del gas [°C]", title="Temperatura")
    axs[1].set(xlabel="Tiempo [s]", ylabel="Presión [bar]", title="Presión")
    axs[0].legend(frameon=False, loc="lower right")
    fig.tight_layout()
    guardar(fig, out / "f04_V1_T_P.png")

    # 5. V2 frente a V3
    fig, ax = plt.subplots(figsize=(10, 5.6))
    for t, m in TAG.items():
        ax.plot(S[("V2", t)]["t_s"], S[("V2", t)]["T_K"] - 273.15, ":", color=COL[m], lw=2.4)
        ax.plot(S[("V3", t)]["t_s"], S[("V3", t)]["T_K"] - 273.15, color=COL[m], lw=2.8, label=ETQ[m])
    ax.plot([], [], ":", color="#475569", lw=2.4, label="V2 adiabático")
    ax.plot([], [], color="#475569", lw=2.8, label="V3 con pared")
    ax.set(xlabel="Tiempo [s]", ylabel="Temperatura del gas [°C]",
           title="La pared enfría el gas y alarga el llenado")
    ax.legend(frameon=False, ncol=2, loc="lower right")
    guardar(fig, out / "f05_V2_vs_V3.png")

    # 6. Factor Z en V3
    fig, ax = plt.subplots(figsize=(10, 5.2))
    for t, m in TAG.items():
        d = S[("V3", t)]
        ax.plot(d["P_bar"], d["Z"], color=COL[m], lw=2.6, label=ETQ[m])
    ax.axhline(1, color="#475569", lw=1)
    ax.set(xlabel="Presión del depósito [bar]", ylabel="Z = Pv/(RT) [-]", title="V3: factor de compresibilidad")
    ax.legend(frameon=False, loc="upper left")
    guardar(fig, out / "f06_V3_Z.png")

    # 7. Sensibilidad a UA gas-pared (3 EOS)
    sens = [r for r in leer(datos / "sensibilidad_3EOS" / "sensibilidad_pared_3EOS.csv")
            if r["parametro"] == "UA_gas_pared"]
    fig, axs = plt.subplots(1, 2, figsize=(12.5, 5.0))
    for m in COL:
        rr = [r for r in sens if r["modelo"] == m]
        axs[0].plot([r["valor"] for r in rr], [r["Tmax_K"] - 273.15 for r in rr], "o-", lw=2.4, ms=8, color=COL[m], label=ETQ[m])
        axs[1].plot([r["valor"] for r in rr], [r["tf_s"] for r in rr], "o-", lw=2.4, ms=8, color=COL[m], label=ETQ[m])
    axs[0].set(xlabel="UA gas-pared [W/K]", ylabel="Temperatura máxima [°C]", title="Temperatura máxima")
    axs[1].set(xlabel="UA gas-pared [W/K]", ylabel="Tiempo hasta 190 bar [s]", title="Tiempo de llenado")
    axs[0].legend(frameon=False)
    fig.tight_layout()
    guardar(fig, out / "f07_sensibilidad_UAgw.png")

    # 8. Malla de suministro: PR - ideal (tiempo y temperatura final) a 140 bar
    dif = [r for r in leer(datos / "peng_robinson" / "V3_diferencias.csv")
           if r["modelo"] == "Peng–Robinson" and r["referencia"] == "Ideal"]
    ps = sorted({r["Ps_bar"] for r in dif})
    ts = sorted({r["Ts_K"] for r in dif})
    for clave, titulo, nombre in (("delta_tf_s", "PR − ideal: tiempo hasta 140 bar [s]", "f08_malla_tiempo.png"),
                                  ("delta_Tf_K", "PR − ideal: temperatura final [K]", "f09_malla_temperatura.png")):
        g = np.array([[next(r[clave] for r in dif if r["Ps_bar"] == p and r["Ts_K"] == t) for t in ts] for p in ps])
        vmax = np.max(np.abs(g))
        fig, ax = plt.subplots(figsize=(8.2, 5.8))
        im = ax.imshow(g, cmap="RdBu_r", vmin=-vmax, vmax=vmax, origin="lower", aspect="auto")
        for i, j in np.ndindex(g.shape):
            ax.text(j, i, f"{g[i, j]:+.1f}", ha="center", va="center", fontsize=13,
                    color="white" if abs(g[i, j]) > 0.6 * vmax else "#0F172A")
        ax.set_xticks(range(len(ts)), [f"{t:.0f}" for t in ts])
        ax.set_yticks(range(len(ps)), [f"{p:.0f}" for p in ps])
        ax.set(xlabel="Temperatura de suministro [K]", ylabel="Presión de suministro [bar]", title=titulo)
        ax.grid(False)
        fig.colorbar(im, ax=ax, shrink=0.85)
        guardar(fig, out / nombre)

    # 9. Descomposición secuencial a 200 bar y 300 K
    dec = [r for r in leer(datos / "peng_robinson" / "V3_descomposicion_secuencial.csv")
           if r["Ps_bar"] == 200 and r["Ts_K"] == 300 and r["magnitud"] == "tf_s"][0]
    a, b, c = dec["delta_deposito"], dec["delta_entalpia"], dec["delta_total"]
    fig, ax = plt.subplots(figsize=(7.6, 5.6))
    ax.bar(0, a, color="#2563EB", width=0.6)
    ax.bar(1, b, bottom=a, color="#EA580C", width=0.6)
    ax.bar(2, c, color="#059669", width=0.6)
    for x, y, v in ((0, a, a), (1, a + b, b), (2, c, c)):
        ax.text(x, y + (1.2 if v > 0 else -2.6), f"{v:+.2f} s", ha="center", fontsize=15, weight="bold")
    ax.axhline(0, color="#0F172A", lw=1)
    ax.set_xticks([0, 1, 2], ["solo el\ndepósito a PR", "luego h de\nentrada a PR", "efecto neto\nPR − ideal"])
    ax.set(ylabel="Cambio en el tiempo de llenado [s]", ylim=(min(a, c) - 5, 5),
           title="Ps = 200 bar, Ts = 300 K, objetivo 140 bar")
    ax.grid(axis="x", visible=False)
    guardar(fig, out / "f10_descomposicion.png")

    # 10. Contraste con Span a 300 K
    heos = [r for r in leer(datos / "peng_robinson" / "propiedades_contra_HEOS.csv")
            if r["origen"] == "malla" and r["T_K"] == 300]
    fig, axs = plt.subplots(1, 2, figsize=(12.5, 5.0))
    for m in COL:
        rr = sorted((r for r in heos if r["modelo"] == m), key=lambda r: r["P_bar"])
        axs[0].plot([r["P_bar"] for r in rr], [r["error_rho_pct"] for r in rr], "o-", lw=2.4, ms=7, color=COL[m], label=ETQ[m])
        axs[1].plot([r["P_bar"] for r in rr], [r["error_h_J_kg"] / 1000 for r in rr], "o-", lw=2.4, ms=7, color=COL[m], label=ETQ[m])
    for ax in axs:
        ax.axhline(0, color="#0F172A", lw=1)
    axs[0].set(xlabel="Presión [bar]", ylabel="Error de densidad [%]", title="Densidad frente a Span")
    axs[1].set(xlabel="Presión [bar]", ylabel="Error de h residual [kJ/kg]", title="Entalpía residual frente a Span")
    axs[0].legend(frameon=False)
    fig.tight_layout()
    guardar(fig, out / "f11_contraste_span_300K.png")

    # 11. Orden de RK4
    rk = leer(datos / "base" / "orden_RK4_paso_fijo_V3.csv")
    h = np.array([r["paso_fijo_s"] for r in rk])
    e = np.array([r["error_Tgas_K"] for r in rk])
    fig, ax = plt.subplots(figsize=(9, 5.4))
    ax.loglog(h, e, "o-", color="#2563EB", lw=2.4, ms=9, label="RK4 paso fijo vs DOP853")
    ax.loglog(h, e[0] * (h / h[0]) ** 4, "--", color="#64748B", lw=1.8, label="pendiente 4 (referencia)")
    orden = [r["orden_estimado_T"] for r in rk if isinstance(r["orden_estimado_T"], float)]
    ax.set(xlabel="Paso h [s]", ylabel="Error en Tg(180 s) [K]",
           title=f"Orden observado: {min(orden):.2f} a {max(orden):.2f}")
    ax.legend(frameon=False)
    guardar(fig, out / "f12_orden_rk4.png")

    # 12. Cierre de balances
    bal = leer(datos / "base" / "balances.csv")
    balpr = leer(datos / "peng_robinson" / "validacion_balances_PR.csv")
    etiquetas = [f"{r['version']} {'ideal' if r['modelo'] == 'Ideal' else 'VdW'}" for r in bal]
    valores = [r["global_rel"] for r in bal]
    etiquetas += [f"{r['version']} PR" for r in balpr[:3]]
    valores += [r["global_J_rel"] for r in balpr[:3]]
    colores = ["#2563EB" if "ideal" in x else "#EA580C" if "VdW" in x else "#059669" for x in etiquetas]
    fig, ax = plt.subplots(figsize=(11, 5.0))
    ax.bar(range(len(valores)), valores, color=colores)
    ax.set_yscale("log")
    ax.set_xticks(range(len(valores)), etiquetas, rotation=0, fontsize=12)
    ax.set(ylabel="Residuo energético global relativo [-]",
           title=f"Balance de energía: residuo máximo {max(valores):.1e} (casos base)")
    ax.grid(axis="x", visible=False)
    guardar(fig, out / "f13_balances.png")

    # 13. Límite isotermo
    iso = leer(datos / "peng_robinson" / "V2_limite_isotermo.csv")
    fig, ax = plt.subplots(figsize=(10, 5.2))
    x = np.arange(len(iso))
    ax.bar(x - 0.2, [r["tf_adiabatico_s"] for r in iso], 0.4, color="#F97316", label="adiabático (V2)")
    ax.bar(x + 0.2, [r["tf_isotermo_s"] for r in iso], 0.4, color="#0EA5E9", label="isotermo a 300 K")
    for i, r in enumerate(iso):
        ax.text(i + 0.2, r["tf_isotermo_s"] + 10, f"Q = {r['calor_extraido_kJ']:.0f} kJ", ha="center", fontsize=12)
    ax.set_xticks(x, [ETQ[r["modelo"]] for r in iso])
    ax.set(ylabel="Tiempo hasta 190 bar [s]", ylim=(0, 780), title="Dos extremos: sin calor y con baño perfecto")
    ax.legend(frameon=False, loc="upper center", ncol=2)
    ax.grid(axis="x", visible=False)
    guardar(fig, out / "f14_isotermo.png")
