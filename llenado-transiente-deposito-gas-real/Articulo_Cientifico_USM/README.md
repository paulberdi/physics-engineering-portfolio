# Informe científico ICF130 — formato UTFSM

Versión principal: **[Informe_ICF130_Formato_USM.pdf](Informe_ICF130_Formato_USM.pdf)** (21 páginas).

**Fuente para Overleaf:** `main.tex`, `referencias.bib`, `logo-usm.png`, `figuras/` y `tablas/`. Usar **pdfLaTeX**.

```bash
pdflatex -interaction=nonstopmode main.tex
bibtex main
pdflatex -interaction=nonstopmode main.tex
pdflatex -interaction=nonstopmode main.tex
```

`codigo_corregido/` permite reproducir modelos y pruebas. `datos/` contiene los CSV de las figuras y `scripts/` el procesamiento de tablas y figuras. `Version_Anterior/` archiva una presentación editorial previa del manuscrito.

**Precaución:** es un manuscrito académico sin revisión por pares y sin validación experimental de un estanque real. La comparación contra HEOS/Span es una verificación de propiedades. El uso del logo indica afiliación y no implica aval de la universidad.
