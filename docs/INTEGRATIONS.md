# Integraciones: qué asimilar y qué no

## 1. Estado real de oriC: ORCA + DoriC

Sí está implementado, pero **no donde parece**:

| pieza | dónde está | estado |
|---|---|---|
| ORCA (ML) | `scripts/oric/03_orca_predict.py`, `35b_orca_from_gbk.py`, `webapp/oric.py` | implementado para los 400 y para FASTA nuevo |
| DoriC | `scripts/oric/02_doric_lookup.py`, `data/doric/` | lookup de verdad de referencia |
| Marcadores dnaA/dnaN | `scripts/oric/01_validate_markers.py` | implementado |
| GC skew + cajas DnaA + iterones | `webapp/origin_pipeline.py` | implementado |
| Estados explicables | `validated`, `consensus`, `lineage_specific`, `anomalous`, `unresolved` | implementado |
| Vista en la GUI | pestaña **Origen de replicación** | lee tablas precalculadas |

Limitaciones actuales:

- **No es un módulo** de `modules/`; no aparece en `python -m modules --list`.
- El pipeline de FASTA nuevo ya ejecuta ORCA + GC skew + dnaA; **DoriC no**,
  porque un FASTA arbitrario no tiene accession. La pestaña del navegador sigue
  usando datos precalculados para los 400.
- `orcapy` no estaba declarado; ya existe el extra `oric` en `pyproject.toml`.
- El modelo `ORCA_RFC_model_1_4_0.pkl.gz` (37 MB) **no se redistribuye** en el
  repo público; `python -m webapp.setup_data --data-root data --only orca`
  lo descarga desde el repositorio de ORCA. SHA-256 local:
  `8221f68f8e70f26b58b97806ecb0ff3390eb092d912f425c9883331a0eb0cd61`.
- DoriC es una base académica: se consulta localmente, no se redistribuye.

### Qué haría para cerrar esa brecha

1. Añadir extra opcional `oric` con `orcapy==1.0.1`.
2. Documentar descarga/checksum del modelo ORCA.
3. Crear módulo `oric` con `input_kind="fasta"`:
   - ORCA sobre el cromosoma;
   - dnaA/dnaN, GC skew e iterones como apoyo;
   - DoriC solo si el usuario declara un accession conocido;
   - estados y límites iguales a los actuales.
4. Mostrar la pestaña de origen también para genomas subidos.
5. Benchmark contra DoriC con los genomas que sí tienen ground truth.

## 2. Proyectos que conviene asimilar

| proyecto | qué integrar | cómo | licencia | esfuerzo | impacto |
|---|---|---|---|---|---|
| **Bakta** | anotación base y dbxref | importar su GFF3/GenBank, o backend opcional | GPL-3.0 | bajo | alto |
| **Prokka** | anotación base | importar GFF3/GBK | GPL-3.0 | bajo | alto |
| **eggNOG-mapper** | COG/GO/KEGG/EC | importar TSV de función | AGPL-3.0 | bajo | alto |
| **InterProScan** | dominios/familias | importar TSV/JSON | Apache-2.0 | bajo | alto |
| **AlphaFold DB** | estructura por UniProt | lookup por accession antes de ESMFold | servicio | bajo | alto |
| **CarveMe / gapseq** | SBML/FBA | exportar borrador SBML + mapeo Rhea | variable | medio | alto |
| **antiSMASH** | BGCs/metabolitos secundarios | importar JSON/GFF | AGPL-3.0 | medio | medio |
| **Operon-mapper** | operones | benchmark o importar predicciones | no verificado | medio | medio |
| **CheckM / GTDB-Tk** | completitud/contaminación | módulo de calidad | GPL-3.0 | medio | alto |
| **AMRFinderPlus / ResFinder** | resistencia antimicrobiana | módulo/import TSV | GPL-3.0 | medio | medio |
| **PlasmidFinder** | plásmidos | módulo/import | GPL-3.0 | medio | medio |
| **KofamScan** | KO/KEGG | opcional, revisar licencia KEGG | KEGG | medio | medio |
| **Bioconda / Galaxy** | distribución | receta + wrapper | — | medio | alto |

## 3. Regla de integración

Preferir **importar salidas** antes que empaquetar código:

```text
usuario corre Bakta/eggNOG/antiSMASH
  → exporta GFF/TSV/JSON
  → nuestro pipeline lo importa
  → conserva fuente, versión y evidencia
```

Ventajas:

- no se relicencia código GPL/AGPL dentro del nuestro;
- no se duplica mantenimiento;
- el usuario puede usar la versión que ya tiene;
- la trazabilidad queda explícita.

Empaquetar código solo cuando sea imprescindible y con obligaciones de licencia
documentadas.

## 4. Prioridad recomendada

1. **ORCA + DoriC en genomas subidos** (cierra una expectativa ya creada).
2. **Importar Bakta/Prokka** para anotación base.
3. **Importar eggNOG/InterPro** para función y dominios.
4. **AlphaFold DB antes de ESMFold.**
5. **Exportar SBML** para CarveMe/COBRApy.
6. **CheckM/GTDB-Tk** para calidad y taxonomía.
7. **antiSMASH** para BGCs.
8. **Bioconda/Galaxy.**

## 5. Lo que no conviene reimplementar

- PGAP o el envío a NCBI.
- FBA/SBML desde cero: CarveMe/gapseq ya lo hacen.
- antiSMASH: es un proyecto enorme y especializado.
- InterProScan/eggNOG: mejor importar sus resultados.
- Operon-mapper: mejor comparar contra su benchmark que copiar su modelo.

## 6. Riesgos

- **Licencias:** AGPL/GPL y KEGG académico; no redistribuir datos restringidos.
- **Versiones:** cada DB cambia; hay que registrar versión/fecha/checksum.
- **Expectativas:** importar una predicción no la convierte en evidencia nuestra;
  se marca `imported` y se conserva la fuente.
- **ORCA:** solo cromosoma, no plásmidos; modelos y umbrales deben documentarse.

Ver también [`BENCHMARK.md`](BENCHMARK.md) y [`ROADMAP.md`](ROADMAP.md).
