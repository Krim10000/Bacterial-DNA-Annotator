# Licencias de terceros

Este proyecto usa programas, bases de datos y modelos que no son de su autor.
Aquí está cada pieza, su licencia y lo que exige. **La licencia MIT de
`LICENSE` cubre solo el código original**; no reclasifica nada de esta lista.

> Este documento es un resumen de buena fe, no un dictamen legal. Las
> obligaciones exactas están en el texto de cada licencia. Si publicas una
> imagen Docker, un servicio o un binario con estas piezas, revisa los textos
> completos de `LICENSES/`.

## 1 · Herramientas y dependencias

| Componente | Licencia | Uso en el proyecto | Obligación principal |
|---|---|---|---|
| **pyrodigal** | GPL-3.0-or-later | Predicción de CDS (`webapp/annotator.py`) | Copyleft si se distribuye junto con la app |
| **orcapy / ORCA** | GPL-3.0 | Predicción de oriC (`scripts/oric/`) | Copyleft si se redistribuye |
| **DIAMOND** | GPL-3.0 | Búsqueda de homólogos | Copyleft si se redistribuye el binario |
| **tRNAscan-SE** | GPL-3.0 | tRNA | Copyleft si se redistribuye |
| **barrnap** | GPL-3.0-or-later | rRNA | Copyleft si se redistribuye |
| **HMMER** | BSD-3-Clause | Dominios Pfam | Conservar aviso de copyright |
| **Infernal** | BSD-3-Clause | Búsqueda de ARN (Rfam) | Conservar aviso |
| **Foldseek** | GPL-3.0 | Homología estructural | Se descarga aparte; no se versiona |
| **ViennaRNA** | Licencia propia (Universidad de Viena) | Energía de horquillas | Crédito obligatorio; no vender como producto propio |
| **Streamlit** | Apache-2.0 | Interfaz web | Conservar `NOTICE` y aviso |
| **pandas, NumPy, SciPy, scikit-learn, NetworkX, Matplotlib** | BSD-3-Clause | Cálculo y tablas | Conservar aviso |
| **Biopython** | Biopython License (tipo BSD) | FASTA/GFF | Conservar aviso |
| **PyTorch** | Apache-2.0 + BSD/MIT | Modelo de ADN (opcional) | Conservar aviso |
| **fair-esm / ESM-2** | MIT (código) | Puntuación de proteínas | Pesos: verificar términos antes de redistribuir |
| **enzymm / PyJess** | MIT | Filtros estructurales | Conservar aviso |
| **P2Rank** | MIT | Bolsillos | Se invoca por CLI |
| **cobra** | LGPL-2.0+ | **No se usa** (retirado de requirements) | — |
| **vis-network** (`lib/vis-9.1.2`) | Apache-2.0 o MIT | Visualizaciones | Conservar cabecera; se distribuye bajo MIT |
| **tom-select** (`lib/tom-select`) | Apache-2.0 | Selectores | Conservar aviso |
| **sentencepiece, RDKit, OpenMM** | Apache-2.0 / BSD-3 / MIT | Modelo y alostería | Conservar aviso |

### Nota sobre GPL-3.0

pyrodigal, orcapy, DIAMOND, tRNAscan-SE y barrnap son GPL-3.0. El código
original de este repositorio puede licenciarse bajo MIT porque **no incluye**
esas piezas: las instala el usuario o las descarga `install.sh`. Pero si
distribuyes una **obra combinada** (por ejemplo, una imagen Docker con
pyrodigal dentro, o un `.deb` que empaqueta DIAMOND), esa distribución debe
cumplir la GPL-3.0 y ofrecer el código fuente correspondiente.

El `webapp/Dockerfile` instala tRNAscan-SE, barrnap y DIAMOND con `apt`, y
pyrodigal con `pip`. Si publicas esa imagen, trátala como GPL-3.0.

## 2 · Datos y modelos

| Fuente | Licencia / condiciones | ¿Se redistribuye aquí? |
|---|---|---|
| **RegulonDB** | Uso académico, no comercial | No; se descarga con `scripts/validation/regulondb_gt.py` |
| **KEGG** | Suscripción; sin derecho de redistribución | No; el caché local queda fuera de git |
| **DoriC** | Uso académico; sin licencia explícita | No |
| **STRING** | Descarga académica | No |
| **DBTBS** | Registro previo; uso académico | No |
| **BATTER** | Sin licencia declarada | No |
| **Pfam** | CC0 (InterPro) | No (se descarga) |
| **Rfam** | CC0 | No (se descarga) |
| **UniProt / Swiss-Prot** | CC BY 4.0 | No (se descarga); citar |
| **AlphaFold DB** | CC BY 4.0 | No; citar |
| **PDB** | CC0 | No; citar |
| **Rhea** | CC BY 4.0 | Sí, tablas pequeñas; citar |
| **ChEBI** | CC BY 4.0 | Sí, identificadores; citar |
| **M-CSA** | Uso académico; consultar términos | No |
| **AlphaFill** | Uso académico; consultar términos | No |
| **Pesos ESMFold** | Condiciones no comerciales | No |
| **Pesos ESM-2** | MIT declarado en Hugging Face | Se descargan; verificar antes de redistribuir |

Los datos que sí se redistribuyen son tablas derivadas pequeñas y citadas en
`DATA_NOTICE.md`. Las bases grandes se descargan en el equipo del usuario y
quedan fuera de control de versiones (ver `.gitignore`).

## 3 · Textos de licencia

- `LICENSES/GPL-3.0.txt` — GNU General Public License v3.0.
- `LICENSES/Apache-2.0.txt` — Apache License 2.0.
- `NOTICE` — atribuciones de componentes Apache-2.0.

## 4 · Cómo citar

Cada herramienta y base de datos tiene su referencia en la página **Créditos**
de la interfaz (`webapp/credits.py`) y en `CITATION.cff`. Si publicas un
resultado, cita las fuentes que lo hicieron posible, no solo este repositorio.
