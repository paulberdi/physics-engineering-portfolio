"""Animación MP4 de V1, V2 y V3 con las tres ecuaciones de estado.

El modelo es 0D: en cada instante hay una sola temperatura del gas Tg(t) y una
sola temperatura de pared Tw(t). Por eso cada estanque se pinta con un color
UNIFORME igual a su Tg calculada, y el anillo exterior con Tw (solo V3, donde
la pared existe en el modelo). No se dibuja ningún campo espacial inventado.
Todas las cifras en pantalla salen de las soluciones de las EDO.
"""
from __future__ import annotations

import shutil
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
from matplotlib.animation import FuncAnimation, FFMpegWriter
from matplotlib.patches import FancyBboxPatch, Rectangle

from .base import BASE, regimen
from .peng_robinson import MODELS, PR, BAR

FONDO = "#0B1220"
PANEL = "#121C2E"
BORDE = "#2B3A52"
TEXTO = "#EEF3FA"
TENUE = "#9FB0C6"
COLORES = {"Ideal": "#4FC3F7", "Van der Waals": "#FFB454", PR: "#5BE3A4"}
NOMBRES = {"Ideal": "GAS IDEAL", "Van der Waals": "VAN DER WAALS", PR: "PENG-ROBINSON"}
DESCRIPCION = {
    "V1": "V1 · caudal fijo 20 g/s · adiabático · objetivo 100 bar",
    "V2": "V2 · caudal compresible (estrangulado → subsónico) · adiabático · objetivo 190 bar",
    "V3": "V3 · caudal compresible + pared con inercia térmica + ambiente · objetivo 190 bar",
}
T_MIN_C, T_MAX_C = 20.0, 150.0
CMAP = mcolors.LinearSegmentedColormap.from_list(
    "termico", plt.get_cmap("inferno")(np.linspace(0.18, 0.97, 256)))
NORM = mcolors.Normalize(T_MIN_C, T_MAX_C)


def preparar(casos: dict[str, list], frames: int) -> dict:
    """Muestrea cada caso en una malla común de tiempo por versión."""
    datos = {}
    for ver, lista in casos.items():
        tfin = max(c.tf for c in lista)
        t = np.linspace(0.0, tfin, frames)
        series = []
        for c in lista:
            d = c.muestras_en(np.minimum(t, c.tf))
            d["t_local"] = np.minimum(t, c.tf)
            d["tf"] = c.tf
            d["p"] = c.p
            series.append(d)
        datos[ver] = (t, tfin, series)
    return datos


def crear_video(casos: dict[str, list], salida: Path, frames_por_version: int = 120,
                pausa: int = 24, fps: int = 24, dpi: int = 120) -> Path:
    if shutil.which("ffmpeg") is None:
        raise RuntimeError("FFmpeg no está instalado o no está en el PATH (ver README).")
    datos = preparar(casos, frames_por_version)
    versiones = ("V1", "V2", "V3")
    por_version = frames_por_version + pausa

    fig = plt.figure(figsize=(16, 9), dpi=dpi, facecolor=FONDO)
    titulo = fig.text(0.04, 0.945, "Llenado transiente de un depósito rígido con N₂",
                      color=TEXTO, fontsize=24, weight="bold", va="center")
    subtitulo = fig.text(0.04, 0.900, "", color=TENUE, fontsize=13.5, va="center")
    reloj = fig.text(0.96, 0.945, "", color="#CFE9FF", fontsize=19, weight="bold",
                     ha="right", va="center", family="monospace")
    fig.text(0.04, 0.022,
             "Modelo 0D: el color del estanque es su temperatura media Tg(t) y el anillo, la temperatura "
             "media de pared Tw(t) (solo V3). No es un campo espacial ni CFD.",
             color=TENUE, fontsize=10.5)
    fig.text(0.96, 0.900, "ICF130 · UTFSM · 2026-2", color=TENUE, fontsize=11, ha="right", va="center")

    # Estanques
    ejes_tanque = []
    xs = (0.035, 0.335, 0.635)
    for x in xs:
        ax = fig.add_axes((x, 0.445, 0.285, 0.405))
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis("off")
        ejes_tanque.append(ax)
    artistas = []
    for ax, modelo in zip(ejes_tanque, MODELS):
        ax.text(0.5, 0.965, NOMBRES[modelo], color=COLORES[modelo], fontsize=15,
                weight="bold", ha="center", va="center")
        pared = FancyBboxPatch((0.12, 0.17), 0.66, 0.66, boxstyle="round,pad=0,rounding_size=0.10",
                               linewidth=0, facecolor="#3A4658")
        gas = FancyBboxPatch((0.155, 0.205), 0.59, 0.59, boxstyle="round,pad=0,rounding_size=0.075",
                             linewidth=0, facecolor=CMAP(NORM(26.85)))
        ax.add_patch(pared)
        ax.add_patch(gas)
        ax.add_patch(Rectangle((0.0, 0.47), 0.125, 0.06, facecolor="#7C8AA0", linewidth=0))
        ax.annotate("", xy=(0.135, 0.50), xytext=(0.0, 0.50),
                    arrowprops=dict(arrowstyle="-|>", color="#DDF3FF", lw=2.4))
        ax.text(0.0, 0.575, "N₂", color="#DDF3FF", fontsize=10.5, ha="left")
        # Manómetro: barra de P/Pobj
        ax.add_patch(Rectangle((0.84, 0.205), 0.05, 0.59, facecolor="#1A2638",
                               edgecolor=BORDE, linewidth=1.2))
        barra = Rectangle((0.84, 0.205), 0.05, 0.0, facecolor=COLORES[modelo], linewidth=0)
        ax.add_patch(barra)
        ax.plot([0.825, 0.905], [0.795, 0.795], color="#FF6B6B", lw=1.6)
        ax.text(0.865, 0.835, "P/Pobj", color=TENUE, fontsize=8.5, ha="center")
        lectura = ax.text(0.45, 0.50, "", color="white", fontsize=12.5, ha="center", va="center",
                          family="monospace", linespacing=1.55,
                          bbox=dict(boxstyle="round,pad=0.45", facecolor="#0B1220", alpha=0.62,
                                    edgecolor="none"))
        estado = ax.text(0.45, 0.085, "", color=TENUE, fontsize=11, ha="center", va="center")
        artistas.append(dict(pared=pared, gas=gas, barra=barra, lectura=lectura, estado=estado))

    cax = fig.add_axes((0.935, 0.50, 0.011, 0.30))
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=NORM, cmap=CMAP), cax=cax)
    cb.set_label("Temperatura media [°C]", color=TENUE, fontsize=10.5)
    cb.ax.tick_params(colors=TENUE, labelsize=9)
    cb.outline.set_edgecolor(BORDE)

    # Curvas
    curvas_def = (("Tg", "Temperatura del gas [°C]", lambda d: d["Tg"] - 273.15),
                  ("P", "Presión [bar]", lambda d: d["P"] / BAR),
                  ("m", "Masa almacenada [kg]", lambda d: d["m"]),
                  ("mdot", "Caudal de entrada [g/s]", lambda d: d["mdot"] * 1000))
    ejes_curva = []
    for i in range(4):
        ax = fig.add_axes((0.045 + i * 0.238, 0.095, 0.195, 0.27))
        ax.set_facecolor(PANEL)
        for s in ax.spines.values():
            s.set_color(BORDE)
        ax.tick_params(colors=TENUE, labelsize=9)
        ax.grid(True, color="#5C6E86", alpha=0.22, lw=0.7)
        ax.set_xlabel("Tiempo [s]", color=TENUE, fontsize=10)
        ax.set_title(curvas_def[i][1], color=TEXTO, fontsize=11.5, loc="left", pad=6)
        ejes_curva.append(ax)
    lineas = [[ax.plot([], [], color=COLORES[m], lw=2.2, label=m)[0] for m in MODELS]
              for ax in ejes_curva]
    lineas_pared = [ejes_curva[0].plot([], [], "--", color=COLORES[m], lw=1.6)[0] for m in MODELS]
    leyenda = ejes_curva[0].legend(loc="lower right", fontsize=8.5, facecolor=PANEL,
                                   edgecolor=BORDE, labelcolor=TEXTO)
    nota_pared = ejes_curva[0].text(0.02, 0.95, "", transform=ejes_curva[0].transAxes,
                                    color=TENUE, fontsize=8.5, va="top")
    linea_obj = ejes_curva[1].axhline(0, color="#FF6B6B", lw=1.0, ls=":", zorder=0)
    actual = {"ver": None}

    def configurar(ver):
        t, tfin, series = datos[ver]
        subtitulo.set_text(DESCRIPCION[ver])
        pobj = series[0]["p"].P_OBJ / BAR
        limites = ((20, max(160, max(np.max(d["Tg"]) for d in series) - 273.15 + 8)),
                   (0, pobj * 1.08), (0, max(np.max(d["m"]) for d in series) * 1.08),
                   (0, 24))
        for ax, lim in zip(ejes_curva, limites):
            ax.set_xlim(0, tfin * 1.02)
            ax.set_ylim(*lim)
        linea_obj.set_ydata([pobj, pobj])
        nota_pared.set_text("trazos: pared Tw" if ver == "V3" else "")

    def actualizar(i):
        ver = versiones[min(i // por_version, 2)]
        k = min(i % por_version, frames_por_version - 1)
        if actual["ver"] != ver:
            actual["ver"] = ver
            configurar(ver)
        t, tfin, series = datos[ver]
        reloj.set_text(f"{ver} │ t = {t[k]:6.1f} s")
        for d, art in zip(series, artistas):
            tg_c = float(d["Tg"][k] - 273.15)
            art["gas"].set_facecolor(CMAP(NORM(tg_c)))
            if ver == "V3":
                art["pared"].set_facecolor(CMAP(NORM(float(d["Tw"][k] - 273.15))))
            else:
                art["pared"].set_facecolor("#3A4658")
            fr = float(np.clip(d["P"][k] / d["p"].P_OBJ, 0, 1))
            art["barra"].set_height(0.59 * fr)
            art["lectura"].set_text(
                f"Tg = {tg_c:6.1f} °C\nP  = {d['P'][k] / BAR:6.1f} bar\n"
                f"m  = {d['m'][k]:6.3f} kg\nZ  = {d['Z'][k]:6.4f}")
            listo = t[k] >= d["tf"] - 1e-9
            if listo:
                art["estado"].set_text(f"objetivo alcanzado en {d['tf']:.1f} s")
                art["estado"].set_color("#FF8A8A")
            else:
                reg = regimen(float(d["P"][k]), d["p"], ver)
                extra = f" · Tw = {d['Tw'][k] - 273.15:5.1f} °C" if ver == "V3" else " · pared sin intercambio"
                art["estado"].set_text(f"{reg}{extra}")
                art["estado"].set_color(TENUE)
        for j, (clave, _, f) in enumerate(curvas_def):
            for d, ln in zip(series, lineas[j]):
                ln.set_data(d["t_local"][:k + 1], f(d)[:k + 1])
        for d, ln in zip(series, lineas_pared):
            if ver == "V3":
                ln.set_data(d["t_local"][:k + 1], d["Tw"][:k + 1] - 273.15)
            else:
                ln.set_data([], [])
        return []

    total = por_version * 3
    anim = FuncAnimation(fig, actualizar, frames=total, interval=1000 / fps, blit=False, repeat=False)
    salida.parent.mkdir(parents=True, exist_ok=True)
    writer = FFMpegWriter(fps=fps, codec="libx264", bitrate=4000,
                          extra_args=["-pix_fmt", "yuv420p", "-preset", "veryfast",
                                      "-movflags", "+faststart"])
    print(f"[VIDEO] {salida.name}: {total} fotogramas a {fps} fps", flush=True)
    anim.save(str(salida), writer=writer, dpi=dpi)
    # Capturas finales de cada versión para el informe y la presentación
    for n, ver in enumerate(versiones):
        actualizar(n * por_version + frames_por_version - 1)
        fig.savefig(salida.with_name(f"{salida.stem}_{ver}_final.png"), dpi=dpi,
                    facecolor=fig.get_facecolor())
    plt.close(fig)
    return salida
