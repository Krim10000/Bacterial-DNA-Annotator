# Benchmark held-out (Bakta / Prokka / DFAST)

**Fecha:** 2026-11-26 (UTC)  
**Salida principal:** [`summary.md`](summary.md), [`summary.json`](summary.json)  
**Runner:** `benchmarks/external/external_benchmark.py` (sin cambios)  
**Agregador nuevo:** `scripts/benchmark/aggregate_heldout.py`  
**Imagen:** `dna-token-bac/benchmark-external:cpu`
`sha256:2fdb291dc37749aebde733a5dc32e9ca4d4460b3456dcc6dbb4b1d7732d6307b`  
**Herramientas:** Bakta 1.12.2, Prokka 1.15.6, DFAST 1.6.2  
**Hilos:** 4 por genoma · **timeout:** 5400 s por herramienta · **umbral:** 0.5
**Bases montadas:** `data/db/bakta_db` → `/data/bakta_db:ro`;
`data/db/dfast_db` → `/data/dfast_db:ro`.

## Selección

La selección se fijó antes de ejecutar y excluye los cuatro genomas usados para
tuning (`cglutamicum`, `mtb_h37rv`, `paeruginosa`, `pputida`) y el gold de
60 kb. Se eligieron genomas completos de RefSeq con FASTA y GFF locales:

| id | especie | accession | contigs | CDS/tRNA/rRNA de referencia |
|---|---|---|---:|---:|
| `saureus` | *Staphylococcus aureus* | GCF_000013425.1 | 1 | 2767 / 59 / 16 |
| `cjejuni` | *Campylobacter jejuni* | GCF_000009085.1 | 1 | 1572 / 43 / 9 |
| `spneumoniae` | *Streptococcus pneumoniae* | GCF_001457635.1 | 1 | 1931 / 58 / 12 |

Se intentó además `btheta` (*Bacteroides thetaiotaomicron*, GCF_014131755.1,
2 contigs: cromosoma + plásmido). Bakta y DFAST renombran los contigs a
`contig_1`/`contig_2` y `sequence1`/`sequence2`; el runner no empareja seqids
genéricos cuando hay más de un contig, así que esos dos quedaron
`not_comparable`. Prokka sí fue comparable (F1 0.9639). `btheta` se conserva en
`btheta.json`/`btheta.md` como limitación metodológica y **no entra en el
agregado**. No se re-ejecutó con `bakta --keep-contig-headers` ni
`dfast --use_original_name t`; queda como trabajo futuro.

## Comando ejecutado (plantilla por genoma)

```bash
docker run --rm \
  -v /home/camilo/Documentos/GitHub/DNA_TOKEN_BAC/dna_poc:/work \
  -v /home/camilo/Documentos/GitHub/DNA_TOKEN_BAC/dna_poc/data/db/bakta_db:/data/bakta_db:ro \
  -v /home/camilo/Documentos/GitHub/DNA_TOKEN_BAC/dna_poc/data/db/dfast_db:/data/dfast_db:ro \
  dna-token-bac/benchmark-external:cpu \
  --genome /work/data/reference_genomes/<id>/<id>_genomic.fna.gz \
  --reference /work/data/reference_genomes/<id>/<id>_genomic.gff.gz \
  --tools all --run --threads 4 --timeout 5400 \
  --bakta-db /data/bakta_db \
  --extra-arg dfast=--dbroot --extra-arg dfast=/data/dfast_db \
  --save-gffs /work/reports/external_benchmark/heldout/<id>/evidence \
  --json /work/reports/external_benchmark/heldout/<id>.json \
  --md /work/reports/external_benchmark/heldout/<id>.md
```

Los cuatro contenedores terminaron con exit code 0. Los tres JSON del set
agregado (`saureus`, `cjejuni`, `spneumoniae`) tienen `status=ok` en Bakta,
Prokka y DFAST. El detalle de la corrida está en `run.log` y
`run_spneumoniae.log`.

## Agregado

```bash
python3 scripts/benchmark/aggregate_heldout.py \
  --input-dir reports/external_benchmark/heldout \
  --run-log reports/external_benchmark/heldout/run.log \
  --exclude-genome btheta \
  --exclude-reason "Genoma multi-contig: Bakta renombró los contigs a contig_1/contig_2 y DFAST a sequence1/sequence2; el runner no empareja seqids genéricos con más de un contig, por lo que ambos quedaron not_comparable. Se conserva como limitación y no entra en el agregado."
```

## Lectura corta

- F1 micro: Bakta 0.9447, Prokka 0.9449, DFAST 0.9449 (diferencia ≤ 0.02 pp).
- No hay ganador: las diferencias son muy pequeñas y no hay réplicas ni gold
  curado a mano.
- DFAST es el más rápido (852 s sumados vs 1452 s Prokka y 1965 s Bakta), pero
  el tiempo no se midió con réplicas.
- `saureus` es el genoma con más FN (202) y menor recall (≈0.929) en las tres
  herramientas; `spneumoniae` es el de más FP (245) y menor F1 (0.9207).
  Son diferencias de anotación de referencia, no de una herramienta frente a
  otra: las tres predicen casi el mismo conjunto de CDS.

## Límites

La referencia es el gold de esta comparación y hereda sus criterios y errores;
el emparejamiento es 1:1 voraz al 0.5 de solapamiento recíproco; no se evalúa
función, dominios ni rutas; los tiempos son de una sola corrida secuencial por
genoma; y no se declara superioridad entre herramientas. Ver la sección de
limitaciones de `summary.md`.
