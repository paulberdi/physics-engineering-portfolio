# ICF130 · Proyecto computacional final · Tema 7

**Llenado transiente de un depósito con gas real** (gas ideal, Van der Waals y Peng-Robinson).
Paul Berdichewsky · Benjamín Arancibia · Andrés Cerón · UTFSM 2026-2.

## Qué se entrega al profesor
`01_ENTREGA_PROFESOR/` (también en `ICF130_ENTREGA_OFICIAL.zip`):
- `Presentacion_Final.pptx`: presentación editable, 21 diapositivas, con notas del presentador y el video incrustado.
- `Presentacion_Final.pdf`: la misma presentación en PDF (el video aparece como imagen).
- `Codigo_Final/`: programa completo. Ver su `README.md`.

## Carpetas
```
01_ENTREGA_PROFESOR/         presentación y código (entrega oficial)
02_RESULTADOS/
  Datos/                     CSV de todas las simulaciones y pruebas
  Graficos/                  figuras PNG/PDF (base, peng_robinson, sensibilidad_3EOS, presentacion)
  Animaciones_MP4/           video V1-V3 con tres EOS y capturas
  registro_ejecucion.json    versiones, parámetros y tiempos de la corrida
03_DOCUMENTACION_INTERNA/    NO se entrega: material de estudio del grupo
  Informe_Completo.pdf       informe técnico explicativo (14 p.)
  Desarrollo_Matematico.pdf  todas las deducciones paso a paso
  Preguntas_Profesor.pdf     37 preguntas con respuesta oral y explicación
  Guia_Exposicion.pdf        división, tiempos, cifras a memorizar
  Registro_Cambios.md        qué cambió respecto del avance y por qué
  latex/                     fuentes .tex (compilar desde esa carpeta)
04_VALIDACION/               pruebas, diagnósticos y registro de validación
```

## Reproducir todo
```bash
cd 01_ENTREGA_PROFESOR/Codigo_Final
python -m pip install -r requirements.txt
python main.py --salida ../../02_RESULTADOS     # ~3 min con video (requiere FFmpeg)
python -m pytest                                # 25 pruebas
```

## Exposición
Paul: diapositivas 1-5 · Benjamín: 6-12 · Andrés: 13-20. Duración objetivo 18-20 min.
Detalles en `03_DOCUMENTACION_INTERNA/Guia_Exposicion.pdf`.

## Archivos originales
Nada de la carpeta compartida `Proyecto computacional` fue modificado. Este paquete se construyó a partir de
`ICF130_ENTREGA_FINAL`, `ICF130_AMPLIACION` y `ICF130_AVANCE_INTEGRADO`.
