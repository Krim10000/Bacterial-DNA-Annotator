# Revisión: ActSeek, JESS/PyJess y EnzyMM

## Resumen

| herramienta | qué hace | licencia | local | estado |
|---|---|---|---|---|
| **EnzyMM/PyJess** | matching de plantillas catalíticas M-CSA | **MIT** | sí, CPU | **integrado** |
| **PyJess/JESS** | motor de template matching 3D | **MIT** | sí, CPU | base de EnzyMM |
| **ActSeek** | busca sitios activos en AFDB usando seed | **no comercial** | sí | **no integrado** |

## EnzyMM/PyJess

- EnzyMM usa PyJess para buscar arreglos geométricos de residuos catalíticos.
- Plantillas derivadas de **M-CSA** (~16 MB).
- No depende de secuencia ni de pliegue global.
- Cubre ~3/4 de mecanismos EC a tercer nivel.
- CLI local: `enzymm -i query.pdb -o results.tsv`.
- Tiempo observado: **~28–40 s por PDB** en un solo hilo; `-j` paraleliza.
- Sin hit = `no_template`; **no** significa ausencia de función.
- Salida: plantilla, M-CSA ID, UniProt, EC, RMSD, residuos mapeados.

**Decisión:** es el motor principal de anotación estructural automática.

## ActSeek

- Publicado en *Bioinformatics* 2025.
- Busca sitios activos similares en AFDB usando una proteína seed y ≥3 residuos.
- Puede usar una base AFDB o una base custom.
- Requiere seed, lista de UniProt y estructura.
- Licencia **non-commercial** de VTT.
- AFDB completo implica descarga/almacenamiento grande.
- Su enfoque es “encontrar proteínas con sitio similar”, no “anotar un PDB contra M-CSA”.

**Decisión:** no se integra por defecto. Se puede usar como búsqueda opcional si:

1. hay seed con residuos catalíticos;
2. el uso es no comercial;
3. se acepta su licencia;
4. se usa una base custom, no AFDB completo.

## JESS

- Algoritmo de template matching 3D.
- PyJess es la interfaz Cython, MIT.
- EnzyMM lo envuelve con plantillas M-CSA.

## Por qué EnzyMM y no ActSeek para el pipeline

| criterio | EnzyMM | ActSeek |
|---|---|---|
| licencia | MIT | no comercial |
| instalación | pip, 16 MB | pip + base AFDB/custom |
| entrada | PDB query | seed + lista |
| salida | plantilla M-CSA/EC | hits AFDB |
| escala | segundos/PDB | base grande |
| integración | directa | requiere seed |
| uso | anotar candidato | minar AFDB |

## Limitaciones comunes

- Match geométrico ≠ función.
- Plantilla M-CSA limitada a enzimas conocidas.
- AFDB es modelo, no estructura experimental.
- Residuos mapeados pueden ser parciales.
- No prueba sustrato, dirección ni cierre.

## Regla

> **EnzyMM dice “esta geometría se parece a una plantilla catalítica”. No dice “esta enzima hace esta reacción”.**

## Validación rápida

- Positivo conocido: **1LYZ** (lisozima, EC 3.2.1.17).
- EnzyMM devolvió 5 matches con M-CSA 203, RMSD 0.45–0.65 Å.
- `2055` y `P77161` AFDB: **sin match** (`no_template`).
- Eso confirma que el motor funciona y que “sin match” es un resultado válido, no un fallo.
