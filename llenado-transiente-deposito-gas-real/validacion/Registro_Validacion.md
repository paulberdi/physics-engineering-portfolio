# Registro de validación

Entorno: Linux, Python 3.13.18, numpy 2.3.3, scipy 1.16.1, matplotlib 3.10.6, CoolProp 7.2.0, pytest 8.4.2, FFmpeg con libx264.
El grupo usó antes Python 3.12.14 con las mismas versiones de librerías.

## 1. Ejecución
| Comprobación | Resultado | Evidencia |
|---|---|---|
| `python main.py --sin-video` sin errores | Sí, 61-63 s | `resultados/registro_ejecucion.json` |
| `python main.py --solo-video` exporta MP4 | Sí, 92 s; 432 fotogramas H.264, 1920×1080, 24 fps, 18 s; decodificado sin errores | `Diagnosticos/video_ffprobe.txt` |
| `python main.py --caso V3 --eos pr` | Sí, imprime resumen y guarda la serie | prueba manual |
| Pruebas automáticas | 25 de 25 pasaban en ~8 s en la ejecución original | `Pruebas/pytest_resultado.txt` |
| Instalación limpia desde `requirements.txt` | Sí; resultados idénticos (diferencia 0) a `02_RESULTADOS` | `Diagnosticos/reproducibilidad_entorno_limpio.txt` |
| Continuidad con el trabajo del grupo | 74 CSV coinciden con los de OneDrive (< 1e-6 en magnitudes físicas) | `Diagnosticos/comparacion_con_originales.txt` |
| Presentación .pptx válida | Pasa la validación de esquema y relaciones; revisión visual de las 21 diapositivas | `Diagnosticos/presentacion_validacion.txt` |
| Documentos LaTeX | Compilan sin errores ni desbordes de línea | `documentacion/latex` |

## 2. Pruebas físicas y numéricas (pytest)
| Prueba | Qué verifica | Tolerancia |
|---|---|---|
| `test_configuracion_coincide` | La ficha JSON dice lo mismo que el código | 1e-12 relativo |
| `test_v1_ideal_coincide_con_solucion_analitica` | tf y T(m) contra la solución exacta | 1e-8 s, 1e-4 K |
| `test_caudal_continuo_en_la_razon_critica` | Empalme de ramas, 20 g/s inicial, caudal nulo a P = Ps | 1e-6 |
| `test_vdw_con_a_b_nulos_es_ideal` (V1, V2, V3) | Límite ideal de Van der Waals | 1e-9 s |
| `test_v3_sin_conductancias_es_v2` (3 EOS) | V3 se reduce a V2 | 1e-8 s |
| `test_pr_propiedades_y_limite_ideal` | P, h = u + Pv, u por cuadratura, cv y u_v por diferencias, en 20 estados; límite ideal | ver código |
| `test_raices_en_rama_estable` | v > b, ∂P/∂v < 0, P recuperada | 1e-9 relativo |
| `test_dominio_invalido_se_informa` | Errores claros fuera del dominio (v ≤ b, T < Tc, objetivo ≥ suministro) | excepción |
| `test_balances_cierran_v3` (ideal, VdW) y `test_balances_cierran_pr_v3` | Masa y energía global | 1e-9 relativo |
| `test_formulacion_conservativa_pr` | Formulación en (m, U) vs en T | 1e-4 s, 2e-5 K |
| `test_cifras_de_la_presentacion` (9 casos) | Cada cifra de la tabla de la presentación | 0,006 s, 0,006 K, 0,0006 kg |

## 3. Resultados de verificación (corrida final)
- Solución exacta V1 ideal: error máximo de T 2,56e-6 K; tf 2,1e-12 s.
- Balance de energía global, residuo relativo máximo: 8,2e-12 (9 casos base); < 1e-10 en los 28 casos PR.
- RK4 de paso fijo: orden observado 4,14 → 4,03 al refinar.
- Formulación conservativa PR: 1,6e-6 K y 3,3e-9 s.
- Rangos físicos: P entre 1 y 190 bar, T entre 300 y 420 K, v > b, cv > 0 y ∂P/∂v < 0 en todas las evaluaciones de PR (el código se detiene si no se cumplen).

## 4. Lo que NO está verificado
- Validación experimental del estanque: no hay mediciones.
- UA, masa de pared y Cd·A: parámetros de estudio, no calibrados.
- Caudal con propiedades reales: no implementado (hipótesis de comparación).
- La reproducción del video en el computador de la sala: probar antes de exponer.
- El enunciado oficial (`PROYECTOS_ICF130.pdf`) no estuvo disponible en esta revisión; los requisitos se tomaron de la documentación del grupo que lo cita.
