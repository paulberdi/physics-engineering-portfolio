# ICF130 · Proyecto computacional · Tema 7

**Llenado transiente de un depósito rígido con nitrógeno.**
Paul Berdichewsky, Benjamín Arancibia y Andrés Cerón · UTFSM · 2026-2.

El programa resuelve los balances transientes de masa y energía de un depósito de 50 L
que se llena con N₂ desde un suministro a 200 bar y 300 K. Compara tres ecuaciones de
estado (gas ideal, Van der Waals y Peng-Robinson) en tres etapas:

| Etapa | Física | Objetivo |
|---|---|---|
| V1 | caudal fijo de 20 g/s, adiabático | 100 bar |
| V2 | caudal compresible (estrangulado y subsónico), adiabático | 190 bar |
| V3 | V2 + pared con inercia térmica + pérdidas al ambiente | 190 bar |

## Instalación

Requiere Python 3.12 o 3.13 y, solo para el video, FFmpeg.

```bash
python -m venv .venv
# Windows PowerShell
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
# macOS / Linux
.venv/bin/python -m pip install -r requirements.txt
```

FFmpeg:
- Windows: `winget install Gyan.FFmpeg` (luego cerrar y abrir la terminal)
- macOS: `brew install ffmpeg`
- Linux: `sudo apt install ffmpeg`

Comprobar con `ffmpeg -version`. Sin FFmpeg todo funciona salvo el video (`--sin-video`).

## Ejecución

Desde esta carpeta (en Windows usar `.\.venv\Scripts\python.exe` en vez de `python`):

```bash
python main.py                       # todo: cálculos, validaciones, figuras y video (~3 min)
python main.py --sin-video           # todo menos el video (~1 min)
python main.py --solo-video          # solo la animación MP4 (~1,5 min)
python main.py --caso V3 --eos pr    # un caso, con resumen en pantalla (eos: ideal, vdw, pr)
python main.py --salida OTRA_CARPETA # escribir los resultados en otra carpeta
python -m pytest                     # suite de pruebas automáticas
```

No hay que editar el código para cambiar de versión o de ecuación de estado.

## Estructura

```
main.py                      punto de entrada único
modelos/base.py              núcleo ideal y Van der Waals, V1-V3, validaciones y barridos
                             (copia sin cambios físicos del programa de la entrega anterior)
modelos/peng_robinson.py     Peng-Robinson completo, contraste con Span (CoolProp),
                             malla Ps-Ts, límite isotermo y descomposición tanque/entrada
modelos/sensibilidad_pr.py   barridos de pared (UA, Cw) con las tres EOS
modelos/figuras_presentacion.py  figuras 16:9 a partir de los CSV de la misma corrida
modelos/animacion.py         video MP4 V1-V3 con las tres EOS
configuracion/caso_base.json ficha legible de parámetros (una prueba verifica que coincide con el código)
pruebas/test_modelo.py       pruebas automáticas
```

## Salidas (carpeta `resultados/` o la indicada con `--salida`)

```
Datos/base/                 casos ideal-VdW, balances, analítica, límites, convergencia, RK4, evento, barridos
Datos/peng_robinson/        casos con PR, malla Ps-Ts, contraste con Span, isotermo, validaciones PR
Datos/sensibilidad_3EOS/    barridos de pared con las tres EOS
Graficos/...                figuras PNG (300 dpi) y PDF por cada conjunto
Graficos/presentacion/      figuras usadas en las diapositivas
Animaciones_MP4/            video y capturas finales de cada etapa
registro_ejecucion.json     versiones, parámetros, constantes y tiempos de la corrida
```

Los CSV usan UTF-8 con BOM, coma como separador y punto decimal.

## Hipótesis principales

Depósito rígido, una entrada y sin salida, sin trabajo de eje, gas bien mezclado (modelo 0D),
suministro de remanso constante, energía cinética y potencial almacenadas despreciables,
cv ideal constante (742,2 J/(kg K)) como fondo de las tres EOS, conexión modelada con flujo
isentrópico de gas ideal común a las tres EOS (Cd·A fijo), y pared con una sola temperatura (V3).
Detalles, ecuaciones y fuentes en [documentacion/](../documentacion/) y [articulo/](../articulo/).

## Corrección metodológica — 20 noviembre 2026

Se conservan las nueve corridas históricas y se añaden comparaciones controladas en `modelos/correcciones.py`.

- **V1 vs V2 a igual objetivo:** a 100 bar los tres modelos dan resultados indistinguibles porque el flujo V2 sigue estrangulado a 20 g/s. A 190 bar V2 demora aproximadamente 9–11 % más que V1, aislando el efecto del caudal. Se mantiene el mismo estado termodinámico final adiabático.
- **Pared:** estudio adicional de capacidad térmica efectiva para 8, 20 y 54 kg, con `cpw=500 J/(kg K)` y conductancias invariables. **No se afirma que 54 kg sea la masa medida de este depósito:** son escenarios hipotéticos comparables.
- **Resolución EOS:** con variables `(m,T)` y depósito rígido se obtiene `v=V/m` y la presión es explícita en cada RHS; Brent se usa para invertir estados especificados `(P,T)`. La formulación conservativa alternativa de Peng-Robinson invierte la energía interna con Brent.
- **Verificación frente a validación:** los residuos, límites y tests verifican el código; no existe validación del transiente contra mediciones. El caudal se aproxima como ideal para las tres EOS, sin estimar un error porcentual no medido.

Las salidas nuevas están en `Datos/correcciones` y `Graficos/correcciones`; las previas no se sobrescriben. La suite incorpora **9** pruebas adicionales (3 por EOS), además de las 25 originales.

> **Código canónico:** esta es la única copia del programa final en la rama principal. Las nuevas salidas se guardan por defecto en `codigo/resultados/` (ignorada en Git), mientras que los resultados publicados están en [`../resultados/`](../resultados/).
