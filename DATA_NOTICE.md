# Aviso sobre los datos

Este repositorio **no redistribuye** las bases de datos grandes que usa el
pipeline. Se descargan en el equipo del usuario y quedan fuera de control de
versiones. Esta nota explica qué se puede compartir, qué no, y cómo citar cada
fuente.

## 1 · Datos que se descargan (no están en el repositorio)

| Fuente | Cómo se obtiene | Condición de uso | Cita mínima |
|---|---|---|---|
| RegulonDB | `python scripts/validation/regulondb_gt.py` | Académico / no comercial | Salgado et al., *Nucleic Acids Research* (2024) |
| KEGG | `python -m webapp.fetch_kegg` | Suscripción; sin redistribución | Kanehisa & Goto (2000) |
| DoriC | `python scripts/oric/02_doric_lookup.py` | Académico; licencia no explícita | Dong et al. (2022) |
| STRING | Descarga desde string-db.org | Acuerdo académico | Szklarczyk et al. (2023) |
| DBTBS | `scripts/validation/fetch_dbtbs_tables.py` | Registro previo; académico | Sierro et al. (2008) |
| BATTER | `scripts/regulation/fetch_batter.py` | Sin licencia declarada | — |
| Pfam | el script de descarga de Pfam | CC0 | Mistry et al. (2021) |
| Rfam | Descarga desde Rfam | CC0 | Kalvari et al. (2021) |
| UniProt / Swiss-Prot | `python -m webapp.setup_data` | CC BY 4.0 | UniProt Consortium (2023) |
| AlphaFold DB | API | CC BY 4.0 | Varadi et al. (2024) |
| PDB | API | CC0 | wwPDB Consortium (2019) |
| M-CSA | Descarga | Académico; consultar términos | Ribeiro et al. (2018) |
| AlphaFill | API | Académico; consultar términos | Hekkelman et al. (2023) |

Ninguna de estas fuentes se incluye en el historial de git. Si necesitas
reproducir un resultado, ejecuta el script de descarga y anota la versión y la
fecha: las bases cambian.

## 2 · Datos que sí se redistribuyen

Solo tablas derivadas pequeñas, generadas por el propio pipeline y citadas:

- `data/reactions/tsv/` — correspondencias Rhea/ChEBI (CC BY 4.0).
- `data/regulondb/gt/` — **no** se versiona; se regenera.
- Tablas de resultados del proyecto (TSV/JSON) — datos derivados del autor,
  bajo los mismos términos que el código, salvo que incorporen datos de una
  fuente con licencia más restrictiva.

Si reutilizas una tabla derivada, cita también la fuente original.

## 3 · Modelos

- **ESM-2**: el código de `fair-esm` es MIT. Los pesos publicados en Hugging
  Face (`facebook/esm2_*`) declaran MIT; verifícalo antes de redistribuirlos.
- **ESMFold**: los pesos tienen condiciones **no comerciales**. El proyecto no
  los incluye; solo llama a la API de Neurosnap o a una instalación local del
  usuario.
- **Modelo de ADN propio**: los checkpoints no se distribuyen en v1.0 porque no
  se pudieron evaluar sin GPU. Se documenta cómo entrenarlo en
  `docs/GPU_STATUS.md`.

## 4 · Uso comercial

El proyecto combina piezas GPL-3.0 (pyrodigal, orcapy, DIAMOND, tRNAscan-SE,
barrnap) y datos de uso académico (RegulonDB, KEGG, DoriC, STRING). **No se
puede usar comercialmente sin revisar cada componente.** La licencia MIT del
código original no concede derechos sobre esas piezas.

## 5 · Privacidad

La aplicación se ejecuta en el equipo del usuario. Los FASTA que sube no se
envían a ningún servidor, salvo que el usuario active explícitamente ESMFold
por Neurosnap o un servicio equivalente; en ese caso la secuencia sale del
equipo y el usuario debe autorizarlo.
