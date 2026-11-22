# Ghost Genes — descubrimiento de genes no anotados

Motor para encontrar marcos abiertos de lectura (ORFs) en regiones que la
anotación no cubre, traducirlos a péptidos —resolviendo notación ambigua— y
priorizarlos con filtros ligeros antes de gastar cómputo en plegamiento.

Diseñado para el proyecto `DNA_TOKEN_BAC` (400 genomas bacterianos), pero el
núcleo no depende de él.

## Flujo

```
FASTA + GFF3
  └─ regiones intergénicas (con padding en los bordes)
       └─ ORFs en 6 marcos (3 por hebra)
            └─ traducción con ambigüedad IUPAC/N
                 ├─ si el codón es inequívoco → 1 cadena
                 ├─ si hay ≤8 combinaciones → varias cadenas
                 └─ si hay más → cadena consenso con X + bandera
       └─ scoring ligero
            ├─ potencial codificante (Markov de codones CDS vs intergénico)
            ├─ Shine-Dalgarno (RBS) y tipo de inicio
            └─ composición proteica (GRAVY, hidrofobicidad, desorden)
       └─ [opcional] DIAMOND
            ├─ vs Swiss-Prot  → función (EC/UniProt)
            └─ vs 400 proteomas → conservación cruzada
       └─ ranking + TSV/JSON
```

## Uso

```bash
# Un genoma
python3 -m ghost.run \
  --genome GCF_000006985.1_ASM698v1 \
  --data-root /ruta/a/dna_poc/data \
  --out ./out \
  --diamond proteomes --db db/proteomes.dmnd \
  --threads 4

# Solo descubrimiento (sin DIAMOND): ~5 s/genoma
python3 -m ghost.run --genome <GCF> --data-root <data> --out ./out

# Validación (oculta K CDS y mide recuperación + control negativo)
python3 -m ghost.validate --genome <GCF> --data-root <data> --out ./out --k 100

# Cruce con brechas metabólicas
python3 -m ghost.gap_link --genome <GCF> --data-root <data> \
  --candidates out/<GCF>_ghost_candidates.tsv \
  --out out/<GCF>_gap_links.tsv
```

Parámetros relevantes: `--min-region-len` (150), `--min-aa` (60),
`--pad` (30), `--alt-starts`, `--max-variants` (8).

## Salidas

- `<genoma>_ghost_candidates.tsv` — una fila por ORF, ordenada por `rank`.
- `<genoma>_ghost_candidates.json` — lo mismo + secuencia `aa` y metadatos.
- `<genoma>_ghost_summary.json` — estadísticas y tiers.
- `<genoma>_ghost_orfs.faa` — péptidos para DIAMOND.
- `<genoma>_gap_links.tsv` — candidatos que podrían cerrar una brecha.

Columnas clave: `codon_llr` (potencial codificante), `rbs`, `ann_frac`
(fracción que invade un gen anotado), `hit`/`self_hit` (conservación),
`rank` (0–1).

## Validación

Se ocultan 100 CDS anotados, se reconstruyen las regiones sin ellos y se mide
si el motor los recupera (identidad ≥90% y solape ≥80%):

| genoma | tamaño | recall | identidad media | FPR umbral 0.75 | FPR umbral 0.55 |
|---|---|---|---|---|---|
| GCF_000006985.1 (Chlorobaculum tepidum) | 2,15 Mb | 0,99 | 0,99 | 0,000 | 0,000 |
| GCF_000007985.2 | 3,81 Mb | 0,99 | 0,99 | 0,000 | 0,000 |

El control negativo baraja dinucleótidos de las regiones (misma composición,
sin estructura de codones): ningún ORF barajado alcanza los umbrales.

Calibración del potencial codificante: CDS reales 0,62 de media (p10 0,28);
ORFs intergénicos aleatorios −0,24.

## Lecciones de implementación

- Escanear **los 3 marcos** por hebra es obligatorio: con solo el marco 0 el
  recall cae a 0,30.
- El **padding** de 30 nt recupera genes cuyo inicio/parada queda recortado
  por genes vecinos solapados (16/17 fallos en la primera versión).
- Un ORF que cae ≥95% dentro de una feature anotada es un fragmento: se
  descarta. `ann_frac` queda como trazabilidad y penaliza el ranking.
- Los hits contra el **mismo genoma** no son conservación: se marcan
  `self_hit` y no suman.

## Límites actuales / siguientes capas

1. **Modelo de proteína ligero (ESM-2 8M)** para separar péptidos plegables de
   basura; requiere `torch` (no instalado en el equipo, sin GPU).
2. **Plegamiento** solo para el top (ESMFold/AlphaFold + Foldseek contra
   `data/db/foldseek_swissprot`).
3. **Sintenía**: exigir que el ORF esté entre los mismos vecinos en varios
   genomas, no solo un hit puntual.
4. **Variantes ambiguas**: hoy se emite el consenso con X; falta escribir las
   cadenas alternativas de los candidatos top.

## Capa ESM-2 (opcional)

El filtro de proteína usa ESM-2 8M (fair-esm) en un venv aislado. Separa
proteínas reales de secuencias barajadas con **AUC 0,987** (media −0,42 vs
−0,53; log-verosimilitud media por residuo).

```bash
export TORCH_HOME=$PWD/.cache/torch   # el caché no puede ir a ~/.cache
.venv/bin/python -m ghost.run \
  --genome <GCF> --data-root <data> --out ./out \
  --esm --esm-top 100 \
  --diamond proteomes --db db/proteomes.dmnd
```

`--esm-top N` puntúa solo los N mejores por ranking (el costo en CPU es
~0,3 s/péptido). Añade `esm_llr` y `esm_score` y ajusta el `rank` un 15%.

## Contexto filogenético

Los hits se anotan con `hit_genome`, `hit_phylum` y `same_phylum`. Un hit a
otro filo se penaliza ×0,25 (los ORFs cortos producen hits composicionales).
Si `len_ratio` (longitud del homólogo / del ORF) ≥1,5 se penaliza ×0,5, y si
≥2,5 se anula: el ORF es un fragmento de dominio, no un gen completo.

## Cierre de brechas dirigido (`ghost.gap_targets`)

Para cada metabolito que el genoma no puede producir, expande sus reacciones
Rhea a los accesos Swiss-Prot que podrían catalizarlas y comprueba:
cuántas ya están anotadas (brecha de mapeo) y si algún ORF fantasma las
recupera (brecha de gen).

```bash
python3 -m ghost.gap_targets --genome <GCF> --data-root <data> \
  --out-dir ./out --db <data>/db/swissprot.dmnd
```

## Capa estructural (`ghost.structure`)

Pliega el top de candidatos con ESMFold y busca homología con Foldseek.

```bash
python3 -m ghost.structure --genome <GCF> --out-dir ./out --top 5 \
  --foldseek-db <data>/db/foldseek_swissprot \
  --foldseek-bin <proyecto>/foldseek/bin/foldseek
```

Estado: **Foldseek validado** (la beta-lactamasa recupera AF-P62593 al 99,6%,
e-value 2·10⁻⁵³). La API pública de ESMFold devuelve 504 para secuencias de
65–120 aa en el momento de escribir esto; el código reintenta y degrada sin
romper el pipeline. Alternativa pendiente: ESMFold local (requiere compilar
`openfold`, sin GPU).

## Resultado del piloto (10 genomas)

3.118 ORFs, 178 con homólogo, 1 tier alto y 30 medios. El mejor candidato:

- `GCF_000008625.1:1348671-1349030(+)`, 120 aa, rank 0,833, ESM −0,27, RBS 0,95.
- 79,2% de identidad, 88,3% de cobertura, homólogo de igual longitud (0,93)
  en un genoma del mismo filo; función del homólogo: **tioredoxina 1 (trxA)**.
- `ann_frac` 0,11: solo roza un gen vecino, no es un fragmento.

## Priorización de brechas: qué es un objetivo real

Tres intentos, en orden de utilidad:

1. **Rareza global** (`ghost.gap_demand`): cuenta en cuántos genomas falta el
   metabolito. Filtra dITP (379/400), pero deja pasar alcaloides vegetales.
2. **Ancestro de biomasa** (`ghost.biomass`): cierra la red Rhea desde los
   precursores de biomasa. Resultó demasiado permisivo (13.074 ancestros, 43%
   de Rhea): toda brecha queda "en ruta".
3. **Conservación de la reacción por filo** (`ghost.gap_conservation`): para
   cada brecha, ¿algún genoma del mismo filo tiene alguna reacción que la
   produzca?
   - Sí → **pérdida de gen** (objetivo real).
   - No → el metabolito no es del filo (**artefacto**, se descarta).

Resultados:

| genoma | filo | brechas | pérdida de gen | artefacto |
|---|---|---|---|---|
| GCF_000006985.1 | Chlorobiota | 43 | **11** | 32 |
| GCF_000008625.1 | Aquificota | 31 | **6** | 25 |

En *C. tepidum* sobreviven objetivos coherentes: geranilgeranil difosfato
(carotenoides), prefrenato (aminoácidos aromáticos), 6-hidroximetil-7,8-
dihidropterina (folato), S-metil-5'-tioadenosina (ciclo MTA), dUTP,
adenilil-molibdopterina. Se descartan hexadecanoil-CoA, octadecanoil-CoA y
demás acyl-CoA (el filo no usa esa química) y los alcaloides vegetales.

Combinado con `gap_targets`:

```bash
python3 -m ghost.gap_conservation --genome <GCF> --data-root <data> --out-dir ./out
python3 -m ghost.gap_targets --genome <GCF> --data-root <data> --out-dir ./out \
  --db <data>/db/swissprot.dmnd \
  --conservation-tsv out/<GCF>_gap_conservation.tsv
```

## Interfaz gráfica (GitHub-friendly)

```bash
pip install -r webapp/requirements.txt
streamlit run webapp/app.py
```

Sube uno o varios FASTA, marca CDS / tRNA / rRNA / ORFs no anotados / ESM-2 y
descarga GFF3, proteínas y TSV de candidatos. Incluye un genoma de ejemplo de
150 kb (`webapp/samples/demo.fna`), Dockerfile y CI. Ver `webapp/README.md`.

![Interfaz](webapp/screenshot_results.png)

## Uso instalado en el proyecto

```bash
cd ~/Documentos/GitHub/DNA_TOKEN_BAC/dna_poc

bash scripts/ghost/run_one.sh GCF_000006985.1_ASM698v1   # pipeline completo
bash scripts/ghost/run_batch.sh 20                       # lote de prueba
bash scripts/ghost/run_batch.sh                          # 400 (reanudable)
```

Orquestador explícito (candidatos → brechas → reporte) y manifiesto:

```bash
PYTHONPATH=scripts/ghost .venv/bin/python -m ghost.analyze \
  --genome <GCF> --data-root data --out out
```

Pruebas:

```bash
cd scripts/ghost
PYTHONPATH=. ../../.venv/bin/python -m unittest discover -s tests -v
```

## Deduplicación global (95%)

Agrupa los candidatos de varios genomas y deja un representante por cluster con
el número de genomas donde aparece. Es lo que permite correr ESM-2 y Pfam una
vez por cluster y no una vez por ORF.

```bash
PYTHONPATH=scripts/ghost .venv/bin/python -m ghost.dedup \
  --inputs "out/*_ghost_candidates.json" --out out/global --identity 95
```

Salidas: `dedup_clusters.tsv`, `dedup_members.tsv`, `representatives.faa`,
`dedup_summary.json`.

Umbral elegido empíricamente (ver `scripts/ghost/DEDUP.md`): 95%
(precisión 0,964, recall 1,0).

## Tier 2: proteínas ya anotadas como "hypothetical"

Extrae las proteínas cuyo mejor hit es hypothetical/uncharacterized/DUF y las
trata como candidatas (no hay que releer ADN):

```bash
PYTHONPATH=scripts/ghost .venv/bin/python -m ghost.tiers \
  --data-root data --out out/tier2
```

Da ~53.650 candidatas. Al deduplicarlas al 95% aparecen **1.164 clusters
conservados entre genomas**; los hay con anotaciones tan sugerentes como
"uncharacterized acyl-CoA thioester hydrolase" o "uncharacterized
metalloprotease". Ver `out/tier2/CONSERVED.md`.

## Tier 3: modelos alternativos de genes imprecisos

Selecciona CDS con evidencia débil (baja confianza, parcial, corto, inicio
raro o varios inicios en fase) y re-enumera la ventana ±300 nt en 6 marcos.
Conserva extensiones 5', truncamientos y cambios de marco que no son idénticos
al gen anotado.

```bash
PYTHONPATH=scripts/ghost .venv/bin/python -m ghost.tier3 \
  --data-root data --out out/tier3
```

## Resultados sobre 310 genomas de vida libre

- 465.324 ORFs, 456.952 clusters (95%), 2.139 multi-genoma.
- Brechas: 9.504 pérdidas de gen en 828 metabolitos distintos.
- Objetivos más frecuentes: precursor de riboflavina, molibdopterina,
  timidina/citidina (salvamento), metilglioxal, MTA, folato, hemo.
- Detalle en `out_fl/GAPS.md`.

## Refinamiento de brechas y cierre

El modelo crudo mezcla transporte, iones y reacciones no enzimáticas. El
refinador exige: productor enzimático no-transporte, metabolito consumido en la
red del genoma, sin iones/ROS/currency, y conservación por filo.

```bash
PYTHONPATH=scripts/ghost .venv/bin/python -m ghost.gap_refine \
  --data-root data --out out_fl/refined

PYTHONPATH=scripts/ghost .venv/bin/python -m ghost.gap_closure \
  --data-root data --refined out_fl/refined/refined_gaps.tsv \
  --candidates-dir out_fl --out out_fl/closures.tsv
```

Resultado: 8.655 brechas refinadas y **1 cierre con evidencia de dominio**
(glioxalasa → (R)-S-lactoilglutationa). Ver `out_fl/CLOSURES.md`.

## Desiertos de anotación

Busca tramos largos sin ninguna señal (CDS, tRNA/rRNA, oriC, transposones,
defensa, ARN) y los cruza con los ORFs fantasma:

```bash
PYTHONPATH=scripts/ghost .venv/bin/python -m ghost.empty_regions \
  --data-root data --out-fl out_fl --out out_fl/empty_regions --min-len 3000
```

En los 310 genomas: 191 regiones ≥3 kb en 93 genomas; 4 contienen ORFs
conservados (candidatos fuertes a genes no anotados). Ver
`out_fl/EMPTY_REGIONS.md`.

## Etiquetado con procedencia

Toda fila de salida lleva `label`, `label_source`, `label_evidence`,
`label_confidence`, `label_interpretation`, `label_pipeline` y `label_date`.
Ninguna etiqueta existe sin fuente; ver `LABELS.md`.
