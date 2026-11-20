# Artículo científico ICF130 (formato UTFSM)

**Documento principal:** [Informe_ICF130_Formato_USM.pdf](Informe_ICF130_Formato_USM.pdf) (**21 páginas**, según el archivo original).

Aquí están las fuentes independientes del manuscrito:

- **main.tex**, **referencias.bib** y **logo-usm.png**.
- **figuras/** (PDF y PNG) y **tablas/** (fragmentos LaTeX).
- **datos/** (CSV de referencia de las figuras).
- **scripts/** (generación de tablas y figuras).
- **Version_Anterior/** (PDF editorial anterior, conservado como referencia histórica).

La copia anteriormente almacenada en *codigo_corregido/* era idéntica, por SHA de Git, al [código final canónico](../codigo/); se eliminó la duplicación del árbol actual para evitar mantener dos versiones del mismo programa.

## Compilar con pdfLaTeX

~~~bash
pdflatex -interaction=nonstopmode main.tex
bibtex main
pdflatex -interaction=nonstopmode main.tex
pdflatex -interaction=nonstopmode main.tex
~~~

En Overleaf, selecciona **main.tex** y compilador **pdfLaTeX**. Los scripts de figuras leen los CSV de esta misma carpeta, sin depender de rutas del código externo.

**Advertencia:** manuscrito académico sin revisión por pares ni validación experimental del depósito. La comparación HEOS/Span verifica propiedades y no sustituye mediciones físicas. El logo indica afiliación, no aval editorial institucional.
