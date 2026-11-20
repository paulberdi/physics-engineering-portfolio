# Registro de cambios

Respecto del avance (`Avance1+v3.2.pdf`), de la entrega anterior (`ICF130_ENTREGA_FINAL`) y de la ampliación
(`ICF130_AMPLIACION`). Los archivos originales en OneDrive no se modificaron.

## Física
| Cambio | Motivo |
|---|---|
| Ninguno en ideal y Van der Waals | `modelos/base.py` es copia del programa de la entrega anterior; sus 54 CSV se reproducen con diferencias < 1e-6 en las magnitudes físicas |
| Peng-Robinson con u, h y cv residuales (venía de la ampliación) | Se integra sin cambios en el programa único |
| Peng-Robinson agregado a la sensibilidad de pared (UA_gw, UA_wa, Cw) | Antes solo había ideal y VdW; completa la comparación con las tres EOS |
| Se mantiene la conexión con flujo isentrópico de gas ideal para las tres EOS | Hipótesis de comparación declarada: a igual P, igual caudal |

## Visualización
| Cambio | Motivo |
|---|---|
| Animación nueva con tres EOS, estanques de color uniforme = Tg(t) y anillo = Tw(t) | El mapa anterior (`mapa_ilustrativo`) dibujaba un punto caliente en la entrada que no sale de ninguna ecuación y contradice la física (el gas que entra es el más frío). El modelo es 0D |
| Se retiró la referencia a von Kármán | Solo aparecía como referencia [4] del avance v3.2, sin uso en el texto ni en el código; una calle de vórtices no corresponde a este problema |
| Figuras 16:9 para la presentación (`Graficos/presentacion`) | Generadas desde los CSV de la misma corrida |

## Software
| Cambio | Motivo |
|---|---|
| Un solo `main.py` con `--sin-video`, `--solo-video`, `--caso V1/V2/V3 --eos ideal/vdw/pr`, `--salida` | Antes había dos programas separados (`ICF130_PROYECTO_FINAL.py` y `ICF130_AMPLIACION.py`) |
| Paquete `modelos/` | `base.py` (núcleo), `peng_robinson.py` (único cambio: importa desde `.base`; su `main` se separó en `ejecutar()`), `sensibilidad_pr.py`, `figuras_presentacion.py`, `animacion.py` |
| `pruebas/test_modelo.py` con 25 pruebas (pytest) | Las verificaciones existían como `assert` dentro de los programas |
| `configuracion/caso_base.json` | Ficha legible; una prueba verifica que coincide con el código |
| `requirements.txt` unificado (agrega CoolProp y pytest) | El de la entrega anterior no incluía CoolProp |
| `registro_ejecucion.json` | Guarda versiones, parámetros, constantes y tiempos de cada corrida |
