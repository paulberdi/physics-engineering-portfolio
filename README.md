# Paul Berdichewsky
### Portafolio de Ingeniería Civil Física · UTFSM

Soy estudiante de Ingeniería Civil Física en la Universidad Técnica Federico Santa María. Este repositorio reúne proyectos de desarrollo de mayor alcance e informes de laboratorio, separándolos según su naturaleza para que el portafolio sea fácil de recorrer.

## Proyectos de mayor alcance

Estos trabajos tienen desarrollo propio, varias etapas y una estructura reproducible de código, simulación o documentación.

| Proyecto | Área | Contenido |
|---|---|---|
| [Ecualizador analógico de tres bandas](analog-3-band-equalizer/) | Electrónica y señales | Diseño, simulación, montaje, caracterización y selección automática de bandas con MicroPython |
| [Llenado transiente de un depósito con gas real](llenado-transiente-deposito-gas-real/) | Termodinámica computacional | Proyecto semestral ICF130: balances transientes, gas ideal vs. Van der Waals, flujo compresible y pared térmica |

## Informes de laboratorio

La carpeta [`Proyectos/`](Proyectos/) conserva la colección de trabajos experimentales desarrollados en cursos de laboratorio. No se presentan como proyectos semestrales independientes: son informes y experiencias experimentales con sus datos, figuras, código y archivos de apoyo.

- [`Proyectos/informes/`](Proyectos/informes/) — documentos finales en PDF.
- [`Proyectos/laboratorios/`](Proyectos/laboratorios/) — datos, análisis, Tracker, videos, figuras, código y otros materiales asociados a cada experiencia.

## Habilidades representadas

- **Modelado físico y simulación:** formulación de balances, ecuaciones diferenciales y modelos comparativos.
- **Análisis de datos:** Python, regresiones, propagación de incertidumbres y visualización científica.
- **Física experimental:** adquisición de datos, instrumentación y contraste teoría–experimento.
- **Electrónica y señales:** amplificadores operacionales, filtros activos, adquisición ADC y procesamiento digital.
- **Documentación técnica:** informes académicos, figuras reproducibles y organización de resultados.

## Estructura del repositorio

```text
physics-engineering-portfolio/
├── README.md
├── analog-3-band-equalizer/                  # proyecto de electrónica
├── llenado-transiente-deposito-gas-real/    # proyecto semestral de termodinámica
├── Proyectos/                                # colección de informes de laboratorio
│   ├── README.md
│   ├── informes/                             # entregables finales
│   └── laboratorios/                         # datos y material de apoyo
├── inventory.csv                             # inventario de los repositorios importados originalmente
└── VERIFICATION.md                           # procedencia e integridad de la importación original
```

## Autoría y procedencia

Los trabajos académicos pueden ser colaborativos. El ecualizador fue desarrollado por **Paul Berdichewsky y Benjamín Arancibia**. La colección de laboratorios acredita a **Paul Berdichewsky y Rafael Williams**. El proyecto de termodinámica computacional acredita a **Paul Berdichewsky, Benjamín Arancibia y Andrés Cerón**.

`inventory.csv` y `VERIFICATION.md` documentan la procedencia de los materiales que fueron importados desde repositorios previos; no pretenden describir la estructura reorganizada actual del portafolio.
