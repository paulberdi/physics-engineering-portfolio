# Experimental Physics Portfolio

A curated portfolio of experimental physics projects developed at Universidad Técnica Federico Santa María (UTFSM), combining **laboratory measurements, data analysis, physical modeling and reproducible documentation**.

This repository keeps final academic reports in `papers/` and organizes the supporting material for each experiment inside `projects/`. Each documented experiment includes the available raw measurements, processed data, analysis files, figures and final academic deliverable. When the final deliverable is not a paper, such as the viscous-fall poster, it is stored directly with its project.

## Current projects

| Project | Final document | Supporting files |
|---|---|---|
| Análisis de un movimiento en caída libre | [`paper`](./papers/analisis-de-un-movimiento-en-caida-libre.pdf) | [`project`](./projects/analisis-de-un-movimiento-en-caida-libre/) |
| Lanzamiento parabólico: contrastes entre teoría, medición experimental y modelado digital | [`paper`](./papers/lanzamiento-parabolico-contrastes-entre-teoria-medicion-experimental-y-modelado-digital.pdf) | [`project`](./projects/lanzamiento-parabolico/) |
| Movimiento aleatorio en dos dimensiones: análisis del desplazamiento cuadrático medio y clasificación del régimen dinámico | [`paper`](./papers/movimiento-aleatorio-en-dos-dimensiones-analisis-del-desplazamiento-cuadratico-medio-y-clasificacion-del-regimen-dinamico.pdf) | [`project`](./projects/movimiento-aleatorio/) |
| Estudio experimental de la resistencia del aire en el movimiento de un carro sobre riel de aire | [`paper`](./papers/estudio-experimental-de-la-resistencia-del-aire-en-el-movimiento-de-un-carro-sobre-riel-de-aire.pdf) | [`project`](./projects/resistencia-del-aire/) |
| Análisis experimental del comportamiento dinámico de sistemas oscilatorios amortiguados y acoplados | [`paper`](./papers/analisis-experimental-del-comportamiento-dinamico-de-sistemas-oscilatorios-amortiguados-y-acoplados.pdf) | [`project`](./projects/oscilaciones-amortiguadas-y-acopladas/) |
| Caída libre en medios viscosos | [`poster`](./projects/caida-libre-en-medios-viscosos/docs/poster.pdf) | [`project`](./projects/caida-libre-en-medios-viscosos/) |
| Calibración de sensores de temperatura mediante PT100 y sistema manométrico de gas | [`paper`](./papers/calibracion-de-sensores-de-temperatura.pdf) | [`project`](./projects/calibracion-sensores-temperatura/) |

## Repository structure

```text
Proyectos/
├── README.md
├── papers/
│   └── final academic reports
└── projects/
    ├── README.md
    ├── analisis-de-un-movimiento-en-caida-libre/
    ├── lanzamiento-parabolico/
    ├── movimiento-aleatorio/
    ├── resistencia-del-aire/
    ├── caida-libre-en-medios-viscosos/
    ├── oscilaciones-amortiguadas-y-acopladas/
    └── calibracion-sensores-temperatura/
```

## Organization principles

- Final **reports** remain unchanged in `papers/`.
- Final deliverables of another type, such as a **poster**, stay with the corresponding project.
- Raw experimental recordings and instrument-specific files are kept separate from processed results when that distinction is available.
- Analysis code stays with the experiment that generated it.
- Figures and setup material are grouped with their corresponding project.
- Duplicate copies are consolidated when they contain the same data.

This organization is being expanded as the remaining experimental documentation is added.

## Related project

A separate repository documents the design, simulation, implementation and experimental validation of an active three-band analog equalizer:

[`analog-3-band-equalizer`](https://github.com/paulberdi/analog-3-band-equalizer)

## Authors

**Paul Berdichewsky** · **Rafael Williams**  
Civil Engineering in Physics students  
Universidad Técnica Federico Santa María - Santiago, Chile
