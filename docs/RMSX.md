# RMSX / Flipbook: asimilación al proyecto

Fuente: <https://github.com/AntunesLab/rmsx>  
Paper: Beruldsen, de Freitas & Antunes, *Scientific Reports* (2026),
DOI [10.1038/s41598-026-39869-7](https://doi.org/10.1038/s41598-026-39869-7).  
Licencia del repo: MIT.

## Qué es

RMSX calcula **RMSF por slices temporales** de una simulación de dinámica
molecular. Combina la lectura de RMSD/RMSF con mapas de movimiento por residuo y
tiempo. Flipbook escribe modelos PDB por slice y los visualiza en Molstar,
ChimeraX o VMD.

Entradas:

- topología PDB/GRO;
- trayectoria DCD/XTC/AMBER;
- R + `Rscript`;
- MDAnalysis;
- visor Molstar/ChimeraX/VMD.

Salidas:

- mapas RMSX;
- RMSD/RMSF;
- PDB por slice con RMSX en B-factor;
- Flipbook HTML/multimodelo.

## Encaje en el proyecto

### Sí encaja

- proteínas B/C priorizadas con estructura;
- análisis de movimiento/allostería de un candidato;
- visualización de modelos múltiples;
- colorear por B-factor/RMSX en el visor 3D;
- enriquecer la ficha de estructura con “dinámica”, no solo “plegado”.

### No encaja

- no es anotación de genoma;
- no reemplaza AlphaFold/ESMFold;
- no funciona sin trayectoria MD;
- no debe correr sobre 647.588 C;
- no demuestra función;
- no cierra brechas;
- no debe presentarse como “estructura experimental”.

## Decisión

No integrar RMSX como dependencia del núcleo CPU. Añadirlo como **extra
opcional de dinámica** para proteínas seleccionadas.

```text
core annotation: CPU, sin R, sin MD
dynamics extra:    PDB/GRO + trayectoria + R/MDAnalysis
```

## Tareas

### D1 — Visor compatible

- soportar PDB multimodelo;
- colorear por B-factor;
- mostrar RMSX/Flipbook como capa de color;
- mantener el aviso: dinámica no es función.

### D2 — Adaptador opcional

- `webapp/dynamics_layer.py`;
- detectar `Rscript`, MDAnalysis y `rmsx`;
- si faltan: `not_configured`, no romper;
- aceptar PDB + trayectoria;
- salida: mapas RMSX, PDB por slice, TSV de residuos;
- no integrarlo al pipeline por defecto.

### D3 — Validación

- una proteína con estructura y trayectoria de prueba;
- comparar RMSX con RMSF global;
- documentar que no hay verdad experimental;
- medir tiempo/RAM.

### D4 — Documentación

- `docs/DYNAMICS.md`;
- requisitos R/MDAnalysis;
- límites de licencia y datos;
- no incluir trayectorias en el repo.

## Criterio de aceptación

- con PDB+trayectoria: produce mapas y Flipbook;
- sin R/MDAnalysis: `not_configured`, sin excepción;
- sin trayectoria: no se ejecuta;
- toda etiqueta: `dynamics_hypothesis`, `closure=no`;
- no se usa para A/B/C ni para brechas.

## Relación con “mejor anotador”

RMSX no mejora la anotación base. Mejora la **interpretación estructural** de un
candidato priorizado. Es un extra de investigación, no una promesa de función.

## Prerrequisitos

- R ≥4.1 para `ggpattern` (mapas enmascarados);
- `Rscript` en PATH;
- MDAnalysis;
- trayectoria MD;
- visor Molstar o ChimeraX/VMD; el visor 3D actual puede ser el punto de
  partida para B-factor/multimodelo.
