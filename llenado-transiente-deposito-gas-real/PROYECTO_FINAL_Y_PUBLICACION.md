# Llenado transiente de un depósito con gas real — proyecto computacional ICF130

**Autores:** Paul Berdichewsky, Benjamín Arancibia y Andrés Cerón.  
**Universidad:** Universidad Técnica Federico Santa María (UTFSM), Campus Santiago San Joaquín.  
**Asignatura:** Fundamentos de Termodinámica Aplicada (ICF130), segundo semestre de 2026.

## Estructura y versiones

Este proyecto aborda el llenado transiente de un depósito rígido de 50 L con nitrógeno, mediante balances no estacionarios de masa y energía, tres ecuaciones de estado (ideal, Van der Waals y Peng–Robinson) y una evolución incremental del modelo:

- **Version1, Version2 y Version3:** versiones históricas, conservadas sin cambios.
- **[Proyecto_Final_Corregido](Proyecto_Final_Corregido/):** código Python, presentación editable y PDF, documentación matemática, resultados, animación MP4 y verificaciones; incluye la comparación controlada V1–V2 y sensibilidad de la masa térmica de pared.
- **[Articulo_Cientifico_USM](Articulo_Cientifico_USM/):** informe científico con formato académico UTFSM, PDF final, fuente LaTeX compatible con Overleaf, logotipo institucional, figuras vectoriales y rasterizadas, tablas reproducibles, bibliografía y CSV utilizados.

**Informe para lectura:** [Informe_ICF130_Formato_USM.pdf](Articulo_Cientifico_USM/Informe_ICF130_Formato_USM.pdf)  
**Presentación:** [Presentacion_Final.pdf](Proyecto_Final_Corregido/01_ENTREGA_PROFESOR/Presentacion_Final.pdf)  
**Código:** [main.py](Proyecto_Final_Corregido/01_ENTREGA_PROFESOR/Codigo_Final/main.py)  
**Experimentos corregidos:** [correcciones.py](Proyecto_Final_Corregido/01_ENTREGA_PROFESOR/Codigo_Final/modelos/correcciones.py)

## Modelo físico

Depósito rígido y bien mezclado (0D), una entrada y sin salida, estado de suministro de remanso fijo. El caudal es constante en V1, dependiente de la presión con régimen estrangulado/subsónico en V2, y V3 añade una pared concentrada con intercambios gas–pared–ambiente. La termodinámica real incorpora energía interna consistente con la ecuación de estado, derivada mediante una relación de Maxwell.

**Advertencia metodológica:** la ley de caudal implementada utiliza una aproximación isentrópica de gas ideal común a las tres ecuaciones de estado; los parámetros de pared y conexión no provienen de una calibración experimental. La comparación con HEOS/Span es un contraste de propiedades termodinámicas, no una validación experimental del depósito. El análisis de pared de 8, 20 y 54 kg trata masas *térmicas efectivas* para explorar la sensibilidad.

## Reproducibilidad

Requiere Python 3.12/3.13 y las versiones declaradas en `requirements.txt`. Instalar FFmpeg si se desea reproducir la animación. Desde `Proyecto_Final_Corregido/01_ENTREGA_PROFESOR/Codigo_Final`:

```bash
python -m pip install -r requirements.txt
python main.py --sin-video
python -m pytest
python main.py --caso V3 --eos pr
```

El archivo `04_VALIDACION/Pruebas/pytest_resultado.txt` registra las 25 pruebas originales. Se añadieron 9 pruebas para la comparación controlada y la sensibilidad en `pruebas/test_correcciones.py`. Los registros históricos no sustituyen una ejecución nueva en el entorno del lector.

## Publicación y atribución

Este repositorio documenta un proyecto académico computacional; el manuscrito **no corresponde a una publicación revisada por pares**, ni demuestra validación experimental. El escudo USM identifica la afiliación académica y no implica aval editorial institucional. El código, el manuscrito y la presentación se incluyen para evaluación y reproducibilidad académica. 
