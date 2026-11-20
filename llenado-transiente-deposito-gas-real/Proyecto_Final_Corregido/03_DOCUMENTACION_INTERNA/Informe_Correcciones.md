# Proyecto final corregida — ICF130 (copia independiente)

Se conserva la versión histórica y sus resultados. Esta versión agrega experimentos que separan el efecto del caudal del cambio de presión objetivo, y evalúa la masa térmica efectiva de la pared.

## Resultados centrales comprobados

Para `Ps=200 bar`, `Ts=300 K`, `V=50 L`, misma EOS y mismo objetivo:

| EOS | Objetivo | V1 caudal fijo [s] | V2 caudal compresible [s] | Aumento V2 [s] |
|---|---:|---:|---:|---:|
| Ideal | 100 bar | 198.56 | 198.56 | 0.00 |
| Van der Waals | 100 bar | 213.36 | 213.36 | 0.00 |
| Peng-Robinson | 100 bar | 208.53 | 208.53 | 0.00 |
| Ideal | 190 bar | 379.07 | 418.92 | +39.85 |
| Van der Waals | 190 bar | 376.96 | 410.86 | +33.90 |
| Peng-Robinson | 190 bar | 370.75 | 404.88 | +34.14 |

La coincidencia a 100 bar ocurre porque `P/Ps=0.5 < r_crit≈0.5283`: V2 todavía posee caudal estrangulado. Esto verifica que el nuevo fenómeno actúa solo después del umbral.

## Masa térmica efectiva de pared (V3, Pobjetivo=190 bar)

| EOS | m efectiva | Capacidad [J/K] | t [s] | T gas final [K] |
|---|---:|---:|---:|---:|
| Ideal | 8 kg | 4000 | 456.40 | 385.32 |
| Ideal | 54 kg | 27000 | 469.11 | 376.05 |
| Van der Waals | 8 kg | 4000 | 440.75 | 379.44 |
| Van der Waals | 54 kg | 27000 | 449.80 | 373.34 |
| Peng-Robinson | 8 kg | 4000 | 433.60 | 378.59 |
| Peng-Robinson | 54 kg | 27000 | 442.31 | 372.54 |

La masa de 54 kg es **solo una variación hipotética**, no representa una medición de este equipo. Se mantuvieron constantes las conductancias `UA`, por lo que no se está modelando una modificación geométrica realista.

## Qué se corrigió

1. Se dejó de afirmar que V1, V2 y V3 tienen idéntica presión objetivo: V1 histórica usa 100 bar y V2/V3, 190 bar. La comparación controlada nueva usa objetivos comunes.
2. Se aclara que 8 kg es **masa térmica efectiva asumida**, no necesariamente masa total del cilindro.
3. Se detalla la evaluación explícita de presión para `(m,T)` y la inversión mediante Brent para otros estados, en vez de afirmar que se invierte la EOS en cada paso del método original.
4. Se separa la **verificación** matemática/numérica de la **validación experimental**, que sigue pendiente.
5. Se elimina la estimación no comprobada del error del caudal ideal a 200 bar.
6. Se agregan gráficos, nueve pruebas nuevas, y se aclaran las notas de exposición.

## Reproducibilidad

En `01_ENTREGA_PROFESOR/Codigo_Final`:

```bash
python -m pip install -r requirements.txt
python main.py --sin-video --salida ../../02_RESULTADOS
python -m pytest
```

La comparación de V1-V2 y el barrido de pared se regeneran con `modelos/correcciones.py` automáticamente con `main.py`. El código y los CSV se contrastaron en esta revisión, pero el entorno actual no dispone de CoolProp, por lo que **no se pudo reejecutar la suite completa con la instalación oficial**; la validación de entorno limpio deberá realizarse en el equipo de entrega. Los resultados nuevos se obtuvieron ejecutando los módulos originales de EOS sin las funcionalidades de contraste contra Span.

**No se ha establecido validación experimental del llenado real.**
