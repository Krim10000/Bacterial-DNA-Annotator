# Herramientas que se pueden agregar (y con qué evidencia)

Regla de este módulo: nada entra sin declarar **(a)** de dónde sale, **(b)** si es
caso documentado o sospecha, y **(c)** con qué puntuación y umbral.

## Escala de evidencia

| nivel | qué significa | ejemplo |
|---|---|---|
| **A — documentado experimental** | complejo/interacción con evidencia experimental curada | ComplexPortal `ECO:0000353`, PDB, IntAct |
| **B — documentado inferido** | curado, pero por ortología, parálogos o conocimiento previo | ComplexPortal `ECO:0005544/5546/5547` |
| **C — predicho con score** | modelo que devuelve una métrica calibrada | AlphaFold-Multimer (ipTM), Foldseek (TM) |
| **D — sospecha heurística** | correlación genómica/funcional, sin interfaz | operón, sintenia, coexpresión, STRING |

La distinción A/B no es cosmética: de los **273** complejos usables de K-12,
**133 son `ECO:0000353` (experimental)** y **140 son inferidos** (127 por
conocimiento previo, 12 por parálogos, 1 por ortología). "Curado" no es
sinónimo de "medido".

## Catálogo

| herramienta | qué aporta | forma | nivel | score / umbral | licencia | coste | veredicto |
|---|---|---|---|---|---|---|---|
| **ComplexPortal** (ya usado) | complejos curados + estequiometría | `complextab/<taxid>.tsv`, PSI-MI, API JSON | A/B | estrellas 1-5; `ECO:0000353` vs `554x` | CC BY 4.0 | bajo | **integrado** |
| **EcoCyc** | complejos curados de *E. coli* | flat files BioCyc | A/B | códigos de evidencia | licencia BioCyc (click-through) | medio | solo local; **no redistribuir** |
| **PDB + PDBe API** | ensamblajes biológicos resueltos | REST/JSON (`/pdb/entry/summary`, `/complex/details`) | A | resolución ≤3 Å; método | CC0 | bajo-medio | **recomendado** para homología |
| **PDBePISA** | interfaces y ensamblajes | REST (`/pisa/interfaces`) | A/C | área enterrada, ΔG | libre | medio | complemento |
| **IntAct** | interacciones físicas curadas | PSI-MI / TSV | A | método + ECO | CC BY 4.0 | bajo | **recomendado** (PPI) |
| **BioGRID** | interacciones curadas | TSV / REST | A | método | MIT (libre, académico y comercial) | bajo | **recomendado** (PPI) |
| **ComplexPortal ortología** | transferir complejo a genoma no modelo | API + ortología | B | juicio de curador (`ECO:0005544`) | CC BY 4.0 | medio | clave para los 400 genomas |
| **STRING** (ya usado) | asociación funcional | `protein.links.full` | C/D | combined ≥0.7 alto, 0.4-0.7 medio | CC BY 4.0 | ya está | **integrado**; `database` = circular |
| **AlphaFold-Multimer / ColabFold** | estructura del par | PDB + JSON | C | ipTM ≥0.8 alto, 0.6-0.8 gris, <0.6 fallo; pDockQ ≥0.23 | MIT/Apache | alto (GPU) | para el top de pares |
| **Boltz-2** | co-folding y confianza | JSON (`iptm`, `ptm`, `pae`, `pDockQ`) | C | ipTM ≥0.8 / pDockQ ≥0.23 | MIT | alto (GPU) | **encaja** con `scripts/allostery/04_boltz_modal.py` |
| **Foldseek-Multimer** | buscar el complejo en PDB | binario + DB | B/C | TM-score ≥0.7 fuerte, 0.5 mismo fold | GPL-3.0 | medio | **recomendado**; ya hay binario y DB en `data/db/` |
| **Operon-mapper** | operones con confianza | web / TSV | D | confianza 0-1; >0.5 mismo operón | web gratis | bajo | mejora la feature `same_operon` |
| **UniProt "Subunit"** | anotación curada de subunidades | `.dat` / API | A/B | texto libre (requiere NER) | CC BY 4.0 | medio | ROI bajo automático |
| **Humphreys et al. 2024** | PPI bacterianas por coevolución + estructura | paper / datos | C | score del modelo | por confirmar | alto | vigilar; 19 patógenos |
| **Taghipour et al. 2017** | red PPI de alta confianza de *E. coli* | datos suplementarios | C | confianza de la red | por confirmar | bajo | solo *E. coli* |
| **CoreClust / complexómica** | clustering de AP-MS | software | C | score del clúster | — | alto | requiere datos propios; no |

### No recomendadas para este problema

- **AlphaFold DB**: modelos **monómero**; no dice con quién se combina.
- **KEGG MODULE / MetaCyc pathways**: módulos funcionales, no complejos físicos.
- **CORUM**: complejos de mamífero, no bacterias.
- **MCL / ClusterONE**: algoritmos de agrupamiento, no evidencia. Ya probamos MCL:
  Jaccard 0.43 vs 0.72 del top-2 mutuo (encadena rutas).

## Contrato de puntuación

El score del modelo (`0-1`) y la evidencia son cosas distintas y no se suman:

| evidencia | etiqueta de salida | efecto |
|---|---|---|
| A/B, Jaccard ≥0.8 | `documentado` (+ `CPX-…`, ECO) | no se re-puntúa; se cita |
| A/B, todos los genes en un curado | `subcomplejo` | trozo de algo documentado |
| A/B, ≥50 % de genes en un curado | `parcial` | suele ser operón con dos enzimas |
| C, ipTM ≥0.8 y pDockQ ≥0.23 | `predicción_alta` | evidencia estructural, no experimental |
| C, 0.6 ≤ ipTM < 0.8 | `predicción_media` | zona gris |
| D | `sospecha` | score del modelo + señales declaradas |

Implementado para A/B y D: `annotate_evidence.py` (ver abajo). C queda como
siguiente paso.

## Qué agregar primero

1. **Anotación documentado/sospecha** — hecho en `annotate_evidence.py`. Usa
   ComplexPortal y separa `experimental` de `inferida`. Coste: cero red.
2. **Foldseek-Multimer contra PDB** — convierte un candidato en "hay un
   ensamblaje resuelto con TM-score X". Reutiliza el binario y la DB que ya
   están en `data/db/`. Coste: medio; no requiere GPU.
3. **Boltz-2 para el top-20 de pares** — da ipTM/pDockQ para los candidatos que
   ya pasaron el filtro genómico. Coste: GPU; ya hay un esqueleto Modal en
   `scripts/allostery/04_boltz_modal.py`.

## Trampa de circularidad

ComplexPortal, EcoCyc y el canal `database` de STRING son la misma clase de
evidencia que el gold. Sirven para **anotar** ("esto ya está documentado"), no
para **validar**. La validación de `validate_k12.py` excluye `database` a
propósito; si se usa, R@1 sube de 0.917 a 0.921, poco, pero deja de ser una
prueba independiente.

## Referencias

- ComplexPortal, descarga y webservice: <https://www.ebi.ac.uk/complexportal/download>
- Estrellas y ECO: <https://www.ebi.ac.uk/training/online/courses/complex-portal-quick-tour/searching-and-visualising-data-complex-portal/evidence-and-conclusion-ontology-codes-confidence-scores/>
- BioCyc, licencia de flat files: <https://biocyc.org/download-flatfiles.shtml>
- AlphaFold-Multimer, ipTM: <https://www.ebi.ac.uk/training/online/courses/alphafold/inputs-and-outputs/evaluating-alphafolds-predicted-structures-using-confidence-scores/confidence-scores-in-alphafold-multimer/>
- pDockQ (Bryant et al. 2022): <https://www.nature.com/articles/s41467-022-28865-w>
- Foldseek-Multimer: <https://www.nature.com/articles/s41592-025-02593-7>
- Boltz: <https://github.com/jwohlwend/boltz>
- Operon-mapper: <https://academic.oup.com/bioinformatics/article/34/23/4118/5040321>
- BioGRID (MIT): <https://downloads.thebiogrid.org/>
- Transferencia a complejos homólogos (GO_REF:0000114): <https://geneontology.github.io/GO_REF/0000114>
