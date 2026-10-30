# Llenado transiente de un depósito con gas real

Proyecto computacional en desarrollo para **Fundamentos de Termodinámica Aplicada (ICF130)**, Universidad Técnica Federico Santa María, Campus Santiago San Joaquín.

**Autores:** Paul Berdichewsky, Benjamín Arancibia y Andrés Cerón.

## Descripción

Se estudia el llenado transiente de un depósito rígido con nitrógeno desde un suministro a mayor presión. El modelo trata el depósito como un volumen de control abierto y resuelve conjuntamente los balances transientes de masa y energía, comparando dos ecuaciones de estado:

- gas ideal;
- Van der Waals.

El desarrollo se construyó de forma incremental, de manera que cada versión agrega una modificación física sobre la anterior.

## Progresión del modelo

### V1 — Caudal másico constante

- Depósito rígido y modelo de parámetros concentrados (0D).
- Llenado adiabático.
- Caudal másico constante de referencia.
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

El proyecto está implementado en Python. Se utilizan principalmente:

- `NumPy` para cálculo numérico;
- `SciPy` (`solve_ivp`) para integración temporal;
- `SciPy` (`brentq`) para resolver estados de Van der Waals cuando se conocen presión y temperatura;
- `Matplotlib` para gráficos y animaciones.

Las versiones actuales usan integración RK45 adaptativa y eventos para detener la simulación al alcanzar la presión objetivo.

## Estructura

```text
llenado-transiente-deposito-gas-real/
├── README.md
├── requirements.txt
├── .gitignore
├── src/
│   ├── v1_caudal_constante.py
│   ├── v2_caudal_variable.py
│   └── v3_pared_termica.py
└── reports/
    ├── v1_reporte.txt
    ├── v2_reporte.txt
    └── v3_reporte.txt
```

Los CSV, PNG, GIF y MP4 generados por los scripts se consideran resultados reproducibles y no se incluyen por defecto para mantener el repositorio liviano.

## Instalación

```bash
pip install -r requirements.txt
```

Para exportar animaciones MP4 se necesita **FFmpeg** disponible en el sistema.

## Ejecución

```bash
python src/v1_caudal_constante.py
python src/v2_caudal_variable.py
python src/v3_pared_termica.py
```

Cada script guarda sus resultados en el directorio donde se encuentra el archivo.

## Estado actual

- V1: implementada.
- V2: implementada.
- V3: implementada.
- Validaciones de balances: implementadas.
- Comparación con solución analítica del caso ideal adiabático: implementada.
- Próximos pasos: análisis de sensibilidad de parámetros térmicos, evaluación del rango de validez de Van der Waals y estudio de distintas condiciones del suministro.

## Alcance

El modelo es de parámetros concentrados (0D). Las visualizaciones térmicas 2D incluidas en los scripts son ilustrativas y no constituyen una simulación CFD ni una solución espacial de las ecuaciones de transporte.
