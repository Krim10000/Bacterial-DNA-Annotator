# Orígenes de replicación (oriC): implementación, lecciones y traspaso

**Última revisión:** 2026-11-23. **Estado:** ORCA + GC skew + dnaA funcionan en
los 400 y en un FASTA nuevo; DoriC es opcional y solo con accession conocido.

## 1. Qué hace

```text
FASTA / GenBank
  → ORCA (Random Forest) sobre el cromosoma
  → GC skew, cajas DnaA e iterones como apoyo
  → dnaA/dnaN por anotación funcional
  → DoriC si el usuario da un accession
  → estado + posición + confianza + motivo
```

Estados:

| estado | significado |
|---|---|
| `validated` | ≥2 métodos independientes coinciden (ORCA + DoriC, DoriC + skew, etc.) |
| `consensus` | ORCA + proximidad a dnaA, sin verdad externa |
| `lineage_specific` | linaje sin dnaA (Fusobacteriaceae, endosimbiontes, etc.) |
| `anomalous` | los métodos discrepan; se declara el motivo |
| `unresolved` / `unavailable` | no hay candidato o faltan modelo/dependencias |

## 2. Dónde está

| pieza | archivo |
|---|---|
| ORCA para FASTA nuevo | `webapp/oric.py` |
| Árbol de decisión/estados | `webapp/origin_pipeline.py` |
| Pipeline CLI | `webapp/pipeline.py` (`--accession`, `--no-oric`) |
| GUI | pestaña **🧬 Función y origen** |
| ORCA/DoriC de los 400 | `scripts/oric/` |
| Modelo | `data/models/ORCA_RFC_model_1_4_0.pkl.gz` (descarga con `setup_data --only orca`) |
| DoriC | `data/doric/` (no se redistribuye) |

## 3. Métricas reales de los 400

- 400 genomas con estado: **271 validated**, **77 consensus**,
  **10 lineage_specific**, **42 anomalous**.
- 306 genomas tienen verdad DoriC.
- De los 67 de baja confianza con DoriC, **63 quedaron a ≤50 kb** (mediana 0).
- De las 95 predicciones de baja confianza, **94/95 resueltas** al combinar
  ORCA, dnaA, GC skew y cajas; queda **1 en revisión**.
- 157 predicciones caían en el borde del contig (<20 kb o >len−20 kb);
  93 de las 95 bajas eran genomas circulares.

**Advertencia:** `validated` no es una medida independiente de precisión. Muchos
casos usan DoriC como ancla. La evaluación honesta separa DoriC, consenso y skew.

## 4. Lecciones aprendidas

1. **Inyectar dnaA + dnaN juntos puede empeorar ORCA.** En un caso la
   predicción cayó a 0,02; usando **solo dnaA** subió a 0,42 y el error bajó de
   1,3 Mb a 19 kb. El gate debe probar configuraciones, no asumir que más
   marcadores es mejor.
2. **ORCA espera `accession.version`.** Un FASTA subido no lo trae;
   `BioFile.parse_SeqRecord` falla si el ID no tiene punto. En `webapp/oric.py`
   se añade `.1` cuando falta.
3. **Biopython 1.81 eliminó `Bio.Alphabet`.** El código debe usar `Seq(seq)` y
   `record.annotations["molecule_type"]="DNA"` en versiones nuevas, con
   compatibilidad para 1.76. Ya corregido y cubierto por test.
4. **La circularidad del contig sigue siendo un límite.** ORCA asume linearidad;
   el origen en el borde y el gate de marcadores penalizan la predicción. La
   rotación/duplicación de flancos está planificada, no completa.
5. **DoriC no es una verdad universal.** Es una base curada para genomas con
   accession; no se puede consultar para un FASTA arbitrario. Además, su
   `comparative_support` no es confianza experimental.
6. **`validated` puede ser circular.** Si DoriC se usa como ancla, no es una
   validación independiente. Hay que reportar `validated_by` y `attempted`.
7. **Los linajes sin dnaA existen.** Fusobacteriaceae, endosimbiontes,
   cianobacterias y mollicutes requieren una interpretación distinta; la tabla
   `DNAA_LESS` es un prior, no una rama hardcodeada.
8. **Los iterones son débiles en genomas AT-ricos.** En Fusobacterium dominan
   repeticiones de baja complejidad; se necesita un HMM de RPA/Rep, no contar
   k-mers.
9. **ORCA no modela plásmidos.** Solo se procesa el cromosoma; los plásmidos se
   reportan como `skipped`.
10. **El modelo no se redistribuye.** Se descarga desde ORCA; `orcapy==1.0.1`
    es dependencia opcional y el modelo tiene SHA-256 documentado.

## 5. Cómo se ejecuta

### 400 genomas

```bash
python scripts/oric/01_validate_markers.py
python scripts/oric/02_doric_lookup.py
python scripts/oric/03_orca_predict.py
```

### FASTA nuevo

```bash
python -m webapp.pipeline genoma.fna --outdir results --threads 8
# con DoriC si el genoma tiene accession:
python -m webapp.pipeline genoma.fna --outdir results --accession GCF_000006985.1
```

Salidas:

```text
results/<genoma>.oric.json
results/<genoma>.oric.gff3
results/<genoma>.summary.json   # bloque oric
```

## 6. Validación de transición

`GCF_000006985.1_ASM698v1` procesado como FASTA nuevo:

- ORCA + DoriC + dnaA → `validated`;
- posición **32.396**;
- la tabla precalculada de los 400 también tiene **32.396**;
- cobertura de anotación 87,4 %.

Esto demuestra que el flujo nuevo reproduce el resultado del proyecto para el
núcleo de oriC.

## 7. Límites y trabajo pendiente

- Rotación de contigs circulares en el borde (P0 científico).
- HMM de RPA/Rep para linajes sin dnaA.
- Predecir oriC de plásmidos (ORCA no lo hace; usar PlasmidFinder).
- Revisar los 42 anómalos: 23 skew-ORCA, 16 sin convergencia, 3 solo skew.
- No presentar `validated` como accuracy independiente.
- DoriC: respetar licencia académica y no redistribuir.

## 8. Traspaso

Para retomar la línea:

1. Leer este documento y `docs/PIPELINE_PARITY.md`.
2. Ejecutar `pytest tests/test_oric.py`.
3. Probar el pipeline en demo y en `GCF_000006985.1`.
4. Revisar `data/orca_predictions/ORCA_FIX.md` y `ORCA_QUALITY.md` como
   evidencia histórica (no se publican en el release).
5. No cambiar el estado de un genoma sin declarar `validated_by`, `reason` y
   `attempted`.
