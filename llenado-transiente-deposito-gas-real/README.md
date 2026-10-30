# Llenado transiente de un depósito con gas real

Proyecto computacional semestral en desarrollo para **Fundamentos de Termodinámica Aplicada (ICF130)**, Universidad Técnica Federico Santa María, Campus Santiago San Joaquín.

**Autores:** Paul Berdichewsky, Benjamín Arancibia y Andrés Cerón.

## Descripción

Se estudia el llenado transiente de un depósito rígido con nitrógeno desde un suministro a mayor presión. El depósito se modela como un volumen de control abierto y se resuelven conjuntamente los balances transientes de masa y energía, comparando gas ideal y Van der Waals.

El desarrollo se construyó de forma incremental: cada versión introduce una modificación física adicional sobre la anterior para poder identificar su efecto en los resultados.

## Progresión del modelo

### V1 — Caudal másico constante

- Depósito rígido y modelo de parámetros concentrados (0D).
- Llenado adiabático.
- Caudal másico constante.
- Comparación gas ideal vs. Van der Waals.
- Validación del caso ideal mediante solución analítica.

### V2 — Caudal compresible variable

- Reemplaza el caudal constante por una ley de flujo compresible.
- Distingue régimen estrangulado y subsónico.
- El caudal depende de la presión instantánea del depósito.
- Mantiene la comparación gas ideal vs. Van der Waals.

### V3 — Pared térmica

- Conserva el modelo de caudal de V2.
- Introduce una pared con capacidad térmica propia.
- Modela transferencia de calor gas → pared → ambiente.
- Variables dinámicas principales: masa del gas, temperatura del gas y temperatura de la pared.
- Incluye validaciones de balance energético del gas, pared y sistema completo.

## Modelo físico base

Los balances generales utilizados son

```text
dm/dt  = mdot
dUg/dt = mdot * h_in + Qdot
```

En V3 se incorporan además

```text
dUg/dt     = mdot * h_in - Qdot_gw
Cw*dTw/dt  = Qdot_gw - Qdot_wa
Qdot_gw    = UA_gw * (Tg - Tw)
Qdot_wa    = UA_wa * (Tw - Tamb)
```

La presión se obtiene con la ecuación de estado correspondiente. Para Van der Waals:

```text
P = R*T/(v-b) - a/v^2
```

## Método numérico

El proyecto está implementado en Python y utiliza principalmente:

- `NumPy` para cálculo numérico;
- `SciPy` (`solve_ivp`) para integración temporal;
- `SciPy` (`brentq`) para resolver estados de Van der Waals cuando se conocen presión y temperatura;
- `Matplotlib` para gráficos y animaciones.

Las versiones actuales usan integración RK45 adaptativa y eventos para detener la simulación al alcanzar la presión objetivo.

## Estado actual del desarrollo

- V1: implementada.
- V2: implementada.
- V3: implementada.
- Validaciones de balances: implementadas.
- Comparación con solución analítica del caso ideal adiabático: implementada.
- Próximos pasos: sensibilidad de parámetros térmicos, rango de validez de Van der Waals y estudio de distintas condiciones de suministro.

El modelo es de parámetros concentrados (0D). Las visualizaciones térmicas 2D usadas durante el desarrollo son ilustrativas y no constituyen una simulación CFD ni una solución espacial de las ecuaciones de transporte.

## Estructura del repositorio

```text
llenado-transiente-deposito-gas-real/
├── README.md
├── requirements.txt
├── .gitignore
├── src/
│   └── README.md
└── reports/
    └── README.md
```

Los scripts V1, V2 y V3 y sus reportes de ejecución existen actualmente en el OneDrive de trabajo del curso. Se mantienen allí mientras se consolidan los nombres y la organización definitiva antes de versionarlos en esta carpeta. Esto evita presentar en el portafolio archivos intermedios como si fueran una versión final estable.

## Dependencias

```bash
pip install -r requirements.txt
```

Las animaciones MP4 requieren **FFmpeg**. Los archivos generados automáticamente —MP4, GIF y otros resultados reproducibles— están excluidos por defecto mediante `.gitignore` para mantener el repositorio liviano.
