# Proyecto final corregido — ICF130

Entrega técnica correspondiente al llenado transiente de un depósito con gas real (N2). Autores: Paul Berdichewsky, Benjamín Arancibia y Andrés Cerón, UTFSM 2026-2.

- `01_ENTREGA_PROFESOR/`: presentación PDF/PPTX y código Python ejecutable.
- `02_RESULTADOS/`: CSV de correcciones, gráficos, video MP4 y registros de ejecución.
- `03_DOCUMENTACION_INTERNA/`: informe técnico, desarrollo matemático, guía de exposición, preguntas para defensa y fuentes LaTeX.
- `04_VALIDACION/`: pruebas, resultados previos y diagnósticos de consistencia numérica.

**No se modificaron** las versiones V1, V2 y V3 que ya estaban publicadas en este repositorio. Para reproducir:

```bash
cd 01_ENTREGA_PROFESOR/Codigo_Final
python -m pip install -r requirements.txt
python main.py --sin-video
python -m pytest
```

El registro original de 25 pruebas y los experimentos adicionales no constituyen validación experimental. Revisa `03_DOCUMENTACION_INTERNA/Informe_Correcciones.md` para las decisiones de modelado.
