# Llenado transiente de un depósito con gas real

Proyecto computacional semestral de **Fundamentos de Termodinámica Aplicada (ICF130)**, Universidad Técnica Federico Santa María, Campus Santiago San Joaquín.

**Autores:** Paul Berdichewsky, Benjamín Arancibia y Andrés Cerón.

## Descripción

Se estudia el llenado transiente de un depósito rígido con nitrógeno desde un suministro a mayor presión. El depósito se modela como un volumen de control abierto y se resuelven conjuntamente los balances transientes de masa y energía, comparando gas ideal y Van der Waals.

El desarrollo se construyó de forma incremental: cada versión introduce una modificación física adicional sobre la anterior para poder identificar su efecto en los resultados.

## Progresión del modelo

### Version1 — Caudal másico constante

- Depósito rígido y modelo de parámetros concentrados (0D).
- Llenado adiabático.
- Caudal másico constante.
- Comparación gas ideal vs. Van der Waals.
- Validación del caso ideal mediante solución analítica y balance energético.
- Incluye el script original, reporte, imagen final y animación MP4.

### Version2 — Caudal compresible variable

- Sustituye el caudal constante por una ley de flujo compresible.
- Distingue régimen estrangulado y subsónico.
- El caudal depende de la presión instantánea del depósito.
- Mantiene la comparación gas ideal vs. Van der Waals.
- Incluye script, reporte, datos CSV, figura y animación MP4.

### Version3 — Pared térmica

- Conserva el modelo de caudal de Version2.
- Introduce una pared con capacidad térmica propia.
- Modela transferencia de calor gas → pared → ambiente.
- Variables dinámicas principales: masa del gas, temperatura del gas y temperatura de la pared.
- Incluye validaciones del balance energético del gas, la pared y el sistema completo.
- Incluye script, reporte, datos CSV, dos figuras y animación MP4.

## Modelo físico base

Los balances generales utilizados son

```text
dm/dt  = mdot
dUg/dt = mdot * h_in + Qdot
```

En Version3 se incorporan además

```text
dUg/dt     = mdot * h_in - Qdot_gw
Cw*dTw/dt  = Qdot_gw - Qdot_wa
Qdot_gw    = UA_gw * (Tg - Tw)
Qdot_wa    = UA_wa * (Tw - Tamb)
```

Para Van der Waals:

```text
P = R*T/(v-b) - a/v^2
```

El modelo es de parámetros concentrados (0D). Las visualizaciones térmicas 2D son ilustrativas y no constituyen una simulación CFD ni una solución espacial de las ecuaciones de transporte.

## Método numérico

El proyecto está implementado en Python y utiliza principalmente:

- `NumPy` para cálculo numérico;
- `SciPy` (`solve_ivp`) para integración temporal;
- `SciPy` (`brentq`) para resolver estados de Van der Waals cuando se conocen presión y temperatura;
- `Matplotlib` para gráficos y animaciones;
- `Pillow` e `imageio-ffmpeg` para las salidas visuales.

Las versiones actuales usan integración RK45 adaptativa y eventos para detener la simulación al alcanzar la presión objetivo. Los tres scripts incluidos fueron comprobados con `python -m py_compile` antes de incorporarlos al repositorio.

## Estructura del proyecto

```text
llenado-transiente-deposito-gas-real/
├── README.md
├── requirements.txt
├── .gitignore
├── Avance1+v3.2.pdf
├── Version1/
│   ├── Avance1.py
│   ├── ideal_vs_vdw_reporte.txt
│   ├── ideal_vs_vdw_sim.mp4
│   └── ideal_vs_vdw_sim_final.png
├── Version2/
│   ├── Avance2.py
│   ├── flujovariable_datos.csv
│   ├── flujovariable_foto_sim.png
│   ├── flujovariable_reporte.txt
│   └── flujovariable_video_sim.mp4
└── Version3/
    ├── Version3.py
    ├── V3_datos.csv
    ├── V3_reporte.txt
    ├── V3_comparacion_V2.png
    ├── V3_resultado_final.png
    └── V3_simulacion.mp4
```

`Avance1+v3.2.pdf` es el informe de avance del proyecto y documenta la progresión V1 → V2 → V3.

## Dependencias

```bash
pip install -r requirements.txt
```

Las animaciones MP4 requieren **FFmpeg**. En este portafolio se conservaron los archivos reales del proyecto —código, datos, reportes, figuras, animaciones e informe— para que el trabajo pueda revisarse de forma completa.
