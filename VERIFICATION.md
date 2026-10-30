# Verificación de procedencia

Fecha de importación documentada: 2026-10-30.

Este archivo registra la procedencia de los materiales que fueron incorporados originalmente al portafolio desde repositorios previos. No describe la estructura actual del repositorio después de su reorganización.

## Orígenes verificados

- `analog-3-band-equalizer`: commit de origen `06d98b490c5a7b03bad5bb4cb618599c1ab8a880`
- colección experimental `Proyectos`: commit de origen `3fac8ea110a563d8c950b2858d0881d1f7c1d95a`
- archivos verificados en la importación original: 110
- tamaño total de la importación original: 71,701,826 bytes

La comprobación original comparó rutas, tamaños y hashes de blobs Git contra los commits de origen y registró SHA-256 de los archivos en `inventory.csv`.

## Reorganización posterior

Después de la importación, el portafolio se reorganizó para distinguir claramente:

- proyectos de mayor alcance, ubicados como carpetas independientes en la raíz;
- informes de laboratorio, agrupados bajo `Proyectos/informes/`;
- material experimental y de análisis, agrupado bajo `Proyectos/laboratorios/`.

Por esta razón, algunas rutas registradas en `inventory.csv` corresponden a la disposición original y no a la ubicación actual de esos mismos archivos. El inventario se conserva como registro de procedencia, no como índice vigente del árbol del repositorio.
