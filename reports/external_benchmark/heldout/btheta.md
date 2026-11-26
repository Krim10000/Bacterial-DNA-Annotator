# Benchmark externo CPU: Bakta, Prokka y DFAST

**Generado:** 2026-11-26T02:03:09Z

**Esquema:** `external-benchmark-run/v1`

**Alcance:** Comparación local de CDS/tRNA/rRNA de Bakta, Prokka y DFAST contra una referencia GFF3/GBFF. CPU, sin descargas de datos.

**Comando reproducible:**

```bash
/opt/conda/bin/python /opt/benchmark/external/external_benchmark.py --genome /work/data/reference_genomes/btheta/btheta_genomic.fna.gz --reference /work/data/reference_genomes/btheta/btheta_genomic.gff.gz --tools all --run --threads 4 --timeout 5400 --bakta-db /data/bakta_db --extra-arg dfast=--dbroot --extra-arg dfast=/data/dfast_db --save-gffs /work/reports/external_benchmark/heldout/btheta/evidence --json /work/reports/external_benchmark/heldout/btheta.json --md /work/reports/external_benchmark/heldout/btheta.md
```

Este informe no declara superioridad entre herramientas. Las métricas comparan features contra la referencia local; no evalúan función, dominios ni rutas.

## Entorno

- Python: `3.11.18` (`/opt/conda/bin/python`)
- Plataforma: `Linux 6.8.0-149-generic x86_64`
- CPU lógicas: `8`
- Umbral de solapamiento recíproco: **0.50**
- `--run`: **sí**
- Pseudogenes excluidos: **sí**
- Normalización de seqid: **sí**

### Bakta DB

- Estado: `ok` (válida: **sí**).
- Ruta: `/data/bakta_db` (origen: `argument`).
- Versión: `0` (fuente: `version.json`).
- Marcadores: `bakta.db, version.json`.

| herramienta | binario | disponible | versión | PATH / razón |
|---|---|---|---|---|
| Bakta | `bakta` | sí | 1.12.2 | /usr/local/bin/bakta |
| Prokka | `prokka` | sí | 1.15.6 | /usr/local/bin/prokka |
| DFAST | `dfast` | sí | 1.6.2 | /usr/local/bin/dfast |

- Herramientas disponibles: **bakta, prokka, dfast**.
- Herramientas ausentes: **ninguna**.

## Entradas

- FASTA: `data/reference_genomes/btheta/btheta_genomic.fna.gz` (SHA256 `8d3fe1f9dbe7673e…`, 2 contigs, 6304193 bp).
- Referencia: `data/reference_genomes/btheta/btheta_genomic.gff.gz` (formato `gff3`, SHA256 `a451aa59f5b03324…`).
- Conteos de referencia: CDS **4710**, tRNA **70**, rRNA **15**.
- CDS pseudo excluidos en referencia: **91**.
- SHA256 FASTA: `8d3fe1f9dbe7673e82c9095b19a6f5dc775ea6a58bea0d4ec7556a20bd32bfff`.
- SHA256 referencia: `a451aa59f5b03324fabcdbfa0fe3813c116a4ea53f306da17c599754c31ebc7d`.

## Validación

- Inicio: `2026-11-26T01:06:31Z`; fin: `2026-11-26T02:03:09Z`; total: **3398.31 s**.
- JSON declarado: `reports/external_benchmark/heldout/btheta.json`.
- Markdown declarado: `reports/external_benchmark/heldout/btheta.md`.

| check | estado | detalle |
|---|---|---|
| genome_path | ok | FASTA existe y es legible |
| genome_sha256 | ok | SHA256 8d3fe1f9dbe7673e82c9095b19a6f5dc775ea6a58bea0d4ec7556a20bd32bfff |
| reference_path | ok | referencia gff3 con 4795 features |
| reference_sha256 | ok | SHA256 a451aa59f5b03324fabcdbfa0fe3813c116a4ea53f306da17c599754c31ebc7d |
| reference_features | ok | CDS=4710 tRNA=70 rRNA=15 pseudo_excluidos=91 |
| tool_bakta | ok | binario ejecutable: /usr/local/bin/bakta |
| tool_prokka | ok | binario ejecutable: /usr/local/bin/prokka |
| tool_dfast | ok | binario ejecutable: /usr/local/bin/dfast |
| bakta_db | ok | estado=ok; /data/bakta_db |
| result_bakta | warning | estado=not_comparable; No hay contigs en común entre referencia y herramienta tras normalizar seqid. |
| result_prokka | ok | estado=ok; sin incidencias |
| result_dfast | warning | estado=not_comparable; No hay contigs en común entre referencia y herramienta tras normalizar seqid. |

## Resultados

| herramienta | fuente | estado | F1 CDS | F1 tRNA | F1 rRNA | precisión global | recall global | F1 global | segundos |
|---|---|---|---|---|---|---|---|---|---|
| Bakta | executed | no comparable (not_comparable) | — | — | — | — | — | — | 1397.9 |
| Prokka | executed | ok | 0.9633 | 1.0000 | 1.0000 | 0.9433 | 0.9854 | 0.9639 | 1457.7 |
| DFAST | executed | no comparable (not_comparable) | — | — | — | — | — | — | 504.4 |

La precisión/recall/F1 se calculan contra la referencia local. Si una herramienta no está instalada, su estado es `not_available` y no hay métricas que mostrar.

## Bakta

- Estado: **no comparable (not_comparable)**.
- Fuente: `executed`.
- Razón: No hay contigs en común entre referencia y herramienta tras normalizar seqid..
- Bakta DB: `ok` (válida: **sí**, ruta `/data/bakta_db`, versión `0`).
- Versión: `bakta 1.12.2` (estado `ok`).
- Ejecución: **ok**, 1397.91 s, RSS máx. 2676.7 MiB (`2026-11-26T01:06:39Z` → `2026-11-26T01:29:57Z`).
- Salida: `/tmp/external_benchmark_bakta_52f3dpyu/out/benchmark.gff3` (7438927 bytes, SHA256 `e5fce99e2d193139…`).
- Comando: `/usr/local/bin/bakta --output /tmp/external_benchmark_bakta_52f3dpyu/out --prefix benchmark --threads 4 --force --skip-plot --db /data/bakta_db /tmp/benchmark_genome_bdz3sme8.fna`.
- Anotación comparada: `reports/external_benchmark/heldout/btheta/evidence/bakta/benchmark.gff3` (formato `gff3`, SHA256 `e5fce99e2d193139…`).

### Features (Bakta)

| tipo | ref | query | TP | FP | FN | precisión | recall | F1 | exactos | borde distinto | solape medio |
|---|---|---|---|---|---|---|---|---|---|---|---|
| — | — | — | — | — | — | — | — | — | — | — | — |

### Atributos (Bakta)

| tipo | pares | product en ambos | product literal | product tokens | gene literal | locus literal |
|---|---|---|---|---|---|---|
| — | — | — | — | — | — | — |

### Hipotéticas (Bakta)

| tipo | ref CDS | ref con product | ref hyp | % ref | query CDS | query con product | query hyp | % query |
|---|---|---|---|---|---|---|---|---|
| — | — | — | — | — | — | — | — | — |

#### Diferencias

- Falsos negativos (referencia sin pareja): **0**.
- Falsos positivos (herramienta sin pareja): **0**.
- Pares solo con hebra distinta: **0**.
- Pares emparejados con algún borde distinto: **0**.

- Aviso: No hay contigs en común entre referencia y herramienta tras normalizar seqid.

## Prokka

- Estado: **ok**.
- Fuente: `executed`.
- Versión: `prokka 1.15.6` (estado `ok`).
- Ejecución: **ok**, 1457.68 s, RSS máx. 3573.0 MiB (`2026-11-26T01:29:58Z` → `2026-11-26T01:54:16Z`).
- Salida: `/tmp/external_benchmark_prokka_7lpxfyyq/out/benchmark.gff` (7456097 bytes, SHA256 `b7497b001084ce69…`).
- Comando: `/usr/local/bin/prokka --outdir /tmp/external_benchmark_prokka_7lpxfyyq/out --prefix benchmark --cpus 4 --force --quiet /tmp/benchmark_genome_bdz3sme8.fna`.
- Anotación comparada: `reports/external_benchmark/heldout/btheta/evidence/prokka/benchmark.gff` (formato `gff3`, SHA256 `b7497b001084ce69…`).

### Features (Prokka)

| tipo | ref | query | TP | FP | FN | precisión | recall | F1 | exactos | borde distinto | solape medio |
|---|---|---|---|---|---|---|---|---|---|---|---|
| CDS | 4710 | 4924 | 4640 | 284 | 70 | 0.9423 | 0.9851 | 0.9633 | 4168 | 472 | 0.9947 |
| tRNA | 70 | 70 | 70 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 7 | 63 | 0.9882 |
| rRNA | 15 | 15 | 15 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 15 | 0 | 1.0000 |

### Atributos (Prokka)

| tipo | pares | product en ambos | product literal | product tokens | gene literal | locus literal |
|---|---|---|---|---|---|---|
| CDS | 4640 | 4640 | 0.1698 | 0.3280 | 0.6538 | 0.0000 |
| tRNA | 70 | 70 | 0.0000 | 0.6607 | — | 0.0000 |
| rRNA | 15 | 15 | 1.0000 | 1.0000 | — | 0.0000 |

### Hipotéticas (Prokka)

| tipo | ref CDS | ref con product | ref hyp | % ref | query CDS | query con product | query hyp | % query |
|---|---|---|---|---|---|---|---|---|
| CDS | 4710 | 4710 | 965 | 20.5 | 4924 | 4924 | 2852 | 57.9 |
| tRNA | 0 | 0 | 0 | — | 0 | 0 | 0 | — |
| rRNA | 0 | 0 | 0 | — | 0 | 0 | 0 | — |

#### Diferencias

- Falsos negativos (referencia sin pareja): **70**.
- Falsos positivos (herramienta sin pareja): **284**.
- Pares solo con hebra distinta: **3**.
- Pares emparejados con algún borde distinto: **535**.

Las listas de diferencias están truncadas a 100 elementos por lista; los conteos son completos.

> Primeras diferencias no emparejadas

| tipo | feature | seqid | coordenadas | hebra | product/locus |
|---|---|---|---|---|---|
| FN | CDS | NZ_CP040529.1 | 3328..3678 | - | hypothetical protein |
| FN | CDS | NZ_CP040529.1 | 21605..22045 | + | M23 family metallopeptidase |
| FN | CDS | NZ_CP040529.1 | 32668..33498 | + | virulence RhuM family protein |
| FN | CDS | NZ_CP040530.1 | 180752..181282 | - | hypothetical protein |
| FN | CDS | NZ_CP040530.1 | 181294..181707 | - | hypothetical protein |
| FN | CDS | NZ_CP040530.1 | 190977..191198 | + | hypothetical protein |
| FN | CDS | NZ_CP040530.1 | 266185..266343 | + | hypothetical protein |
| FN | CDS | NZ_CP040530.1 | 333100..333180 | - | K(+)-transporting ATPase subunit F |
| FN | CDS | NZ_CP040530.1 | 369162..369365 | - | helix-turn-helix domain-containing protein |
| FN | CDS | NZ_CP040530.1 | 373403..373609 | + | hypothetical protein |
| FN | CDS | NZ_CP040530.1 | 390556..390777 | + | hypothetical protein |
| FN | CDS | NZ_CP040530.1 | 548791..549111 | - | DUF4133 domain-containing protein |
| FN | CDS | NZ_CP040530.1 | 552900..553325 | + | hypothetical protein |
| FN | CDS | NZ_CP040530.1 | 591156..591278 | - | hypothetical protein |
| FN | CDS | NZ_CP040530.1 | 591406..591489 | - | hypothetical protein |
| FN | CDS | NZ_CP040530.1 | 600631..600756 | + | hypothetical protein |
| FN | CDS | NZ_CP040530.1 | 602638..602820 | + | hypothetical protein |
| FN | CDS | NZ_CP040530.1 | 638954..639184 | - | hypothetical protein |
| FN | CDS | NZ_CP040530.1 | 639243..639575 | + | hypothetical protein |
| FN | CDS | NZ_CP040530.1 | 643992..644138 | + | hypothetical protein |
| FN | CDS | NZ_CP040530.1 | 696334..696525 | - | preprotein translocase subunit SecE |
| FN | CDS | NZ_CP040530.1 | 702631..702765 | + | hypothetical protein |
| FN | CDS | NZ_CP040530.1 | 987544..987675 | - | hypothetical protein |
| FN | CDS | NZ_CP040530.1 | 1055923..1056273 | + | helix-turn-helix domain-containing protein |
| FN | CDS | NZ_CP040530.1 | 1086757..1086996 | - | helix-turn-helix domain-containing protein |
| FN | CDS | NZ_CP040530.1 | 1273652..1273852 | + | hypothetical protein |
| FN | CDS | NZ_CP040530.1 | 1404589..1404750 | + | hypothetical protein |
| FN | CDS | NZ_CP040530.1 | 1720689..1720889 | - | hypothetical protein |
| FN | CDS | NZ_CP040530.1 | 1842133..1843065 | + | family 43 glycosylhydrolase |
| FN | CDS | NZ_CP040530.1 | 2251644..2251802 | - | hypothetical protein |
| FN | CDS | NZ_CP040530.1 | 2541650..2541865 | + | tyrosine-type recombinase/integrase |
| FN | CDS | NZ_CP040530.1 | 2544102..2544842 | + | hypothetical protein |
| FN | CDS | NZ_CP040530.1 | 2547191..2547463 | + | hypothetical protein |
| FN | CDS | NZ_CP040530.1 | 2547736..2548158 | - | hypothetical protein |
| FN | CDS | NZ_CP040530.1 | 3264364..3264744 | + | restriction endonuclease subunit S |
| FN | CDS | NZ_CP040530.1 | 3285351..3285902 | - | DUF6452 family protein |
| FN | CDS | NZ_CP040530.1 | 3524165..3524413 | + | hypothetical protein |
| FN | CDS | NZ_CP040530.1 | 3558889..3559179 | + | hypothetical protein |
| FN | CDS | NZ_CP040530.1 | 3560021..3560425 | + | conjugative transposon protein TraM |
| FN | CDS | NZ_CP040530.1 | 3582899..3583072 | + | hypothetical protein |
| FN | CDS | NZ_CP040530.1 | 3596983..3597138 | + | helix-turn-helix domain-containing protein |
| FN | CDS | NZ_CP040530.1 | 3637864..3638010 | - | hypothetical protein |
| FN | CDS | NZ_CP040530.1 | 3644559..3644888 | + | hypothetical protein |
| FN | CDS | NZ_CP040530.1 | 3647162..3647332 | + | hypothetical protein |
| FN | CDS | NZ_CP040530.1 | 3654686..3654898 | + | hypothetical protein |
| FN | CDS | NZ_CP040530.1 | 3759135..3759347 | - | hypothetical protein |
| FN | CDS | NZ_CP040530.1 | 3912048..3912269 | + | hypothetical protein |
| FN | CDS | NZ_CP040530.1 | 4123154..4123903 | - | glycerophosphodiester phosphodiesterase |
| FN | CDS | NZ_CP040530.1 | 4198195..4198551 | + | hypothetical protein |
| FN | CDS | NZ_CP040530.1 | 4215749..4215880 | + | hypothetical protein |
| FN | CDS | NZ_CP040530.1 | 4332681..4332986 | + | sugar phosphate nucleotidyltransferase |
| FN | CDS | NZ_CP040530.1 | 4613090..4613248 | - | hypothetical protein |
| FN | CDS | NZ_CP040530.1 | 4937610..4937759 | + | hypothetical protein |
| FN | CDS | NZ_CP040530.1 | 4961721..4962092 | - | helix-turn-helix domain-containing protein |
| FN | CDS | NZ_CP040530.1 | 4975432..4975641 | - | hypothetical protein |
| FN | CDS | NZ_CP040530.1 | 4996863..4997042 | - | hypothetical protein |
| FN | CDS | NZ_CP040530.1 | 4998355..4998606 | + | hypothetical protein |
| FN | CDS | NZ_CP040530.1 | 5000256..5001503 | - | site-specific integrase |
| FN | CDS | NZ_CP040530.1 | 5002898..5003296 | - | hypothetical protein |
| FN | CDS | NZ_CP040530.1 | 5081127..5081489 | + | four helix bundle protein |
| FN | CDS | NZ_CP040530.1 | 5097989..5098156 | - | hypothetical protein |
| FN | CDS | NZ_CP040530.1 | 5508537..5508689 | + | hypothetical protein |
| FN | CDS | NZ_CP040530.1 | 5533363..5533485 | + | hypothetical protein |
| FN | CDS | NZ_CP040530.1 | 5803103..5804011 | - | family 43 glycosylhydrolase |
| FN | CDS | NZ_CP040530.1 | 5993387..5993563 | - | hypothetical protein |
| FN | CDS | NZ_CP040530.1 | 5996709..5996912 | - | helix-turn-helix domain-containing protein |
| FN | CDS | NZ_CP040530.1 | 5997211..5997378 | + | hypothetical protein |
| FN | CDS | NZ_CP040530.1 | 6021583..6021837 | + | hypothetical protein |
| FN | CDS | NZ_CP040530.1 | 6114767..6114922 | - | hypothetical protein |
| FN | CDS | NZ_CP040530.1 | 6157333..6157707 | - | DUF4374 domain-containing protein |
| FP | CDS | NZ_CP040529.1 | 4..372 | - | hypothetical protein |
| FP | CDS | NZ_CP040529.1 | 3779..3979 | + | hypothetical protein |
| FP | CDS | NZ_CP040529.1 | 6967..7083 | + | hypothetical protein |
| FP | CDS | NZ_CP040529.1 | 9454..9639 | + | hypothetical protein |
| FP | CDS | NZ_CP040529.1 | 21419..22003 | - | hypothetical protein |
| FP | CDS | NZ_CP040529.1 | 27028..27258 | + | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 38346..38651 | - | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 50526..51362 | - | Beta-galactosidase |
| FP | CDS | NZ_CP040530.1 | 51405..54110 | - | Beta-galactosidase |
| FP | CDS | NZ_CP040530.1 | 78537..78671 | + | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 168030..169049 | - | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 173351..173953 | - | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 176238..178265 | - | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 200197..200391 | - | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 200401..201399 | - | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 220908..221666 | + | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 222449..222553 | + | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 228517..228753 | + | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 234175..234291 | - | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 242088..242186 | - | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 350389..352869 | + | Coenzyme A disulfide reductase |
| FP | CDS | NZ_CP040530.1 | 369366..369485 | - | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 373745..374002 | + | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 374102..374620 | - | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 374794..375282 | - | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 388479..388643 | - | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 389619..390533 | + | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 391818..391937 | + | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 392603..392896 | - | Tyrosine recombinase XerC |
| FP | CDS | NZ_CP040530.1 | 393385..393594 | + | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 408662..409438 | - | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 433886..434002 | + | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 480796..480930 | - | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 481078..481167 | + | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 491024..491179 | - | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 524656..524823 | - | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 527966..528460 | + | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 528454..529116 | + | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 536793..536894 | - | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 536929..537243 | - | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 548791..549528 | - | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 580366..580497 | - | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 593083..593295 | - | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 593455..593655 | - | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 595979..596230 | - | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 633780..633890 | + | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 639015..639419 | - | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 696359..696478 | + | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 702908..703093 | + | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 728013..728105 | - | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 728715..728891 | - | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 740587..740691 | - | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 800173..800268 | + | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 800258..800353 | + | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 932025..932174 | + | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 932116..932259 | + | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 993232..993327 | + | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 1008529..1008642 | + | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 1028455..1028664 | + | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 1032563..1032814 | + | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 1067977..1068186 | - | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 1085119..1085226 | - | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 1103383..1103511 | + | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 1154256..1154633 | - | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 1154649..1154762 | - | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 1188090..1188215 | - | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 1202949..1203068 | - | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 1259353..1259565 | + | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 1269579..1269680 | + | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 1286301..1286504 | - | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 1307930..1308289 | + | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 1309301..1309576 | + | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 1316583..1316735 | - | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 1320430..1320756 | - | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 1323368..1323472 | + | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 1338299..1338472 | - | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 1356890..1357063 | + | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 1400577..1400957 | + | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 1403805..1404038 | + | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 1423653..1423790 | - | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 1464543..1464695 | + | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 1476658..1476813 | + | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 1537978..1538139 | - | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 1544705..1544830 | - | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 1695162..1695467 | - | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 1699157..1699264 | - | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 1732663..1732779 | - | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 1738561..1738692 | + | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 1744383..1744532 | - | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 1756281..1756415 | - | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 1782711..1783235 | - | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 1783467..1783742 | - | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 1818974..1819489 | + | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 1841152..1843065 | + | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 1875206..1875304 | + | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 1876850..1877080 | + | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 1877993..1878244 | - | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 1879592..1880104 | + | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 1880114..1881682 | + | hypothetical protein |
| FP | CDS | NZ_CP040530.1 | 1882537..1882797 | + | hypothetical protein |

> Primeros pares con coordenadas no exactas

| tipo | seqid | ref | Δstart | Δend | solape |
|---|---|---|---|---|---|
| CDS | NZ_CP040529.1 | 5801..6358 | +78 | +0 | 0.8602 |
| CDS | NZ_CP040529.1 | 8101..9141 | +0 | +54 | 0.9507 |
| CDS | NZ_CP040529.1 | 9925..10521 | -24 | +0 | 0.9614 |
| CDS | NZ_CP040529.1 | 17590..18072 | -87 | +0 | 0.8474 |
| CDS | NZ_CP040529.1 | 19634..20401 | -123 | +0 | 0.8620 |
| CDS | NZ_CP040529.1 | 22291..22959 | -39 | +0 | 0.9449 |
| CDS | NZ_CP040529.1 | 27464..28756 | -129 | +0 | 0.9093 |
| CDS | NZ_CP040530.1 | 22045..23691 | +0 | +39 | 0.9769 |
| CDS | NZ_CP040530.1 | 26551..28098 | +0 | -36 | 0.9767 |
| CDS | NZ_CP040530.1 | 46740..47969 | -3 | +0 | 0.9976 |
| CDS | NZ_CP040530.1 | 57301..58272 | +0 | -21 | 0.9784 |
| CDS | NZ_CP040530.1 | 63633..64601 | +0 | -3 | 0.9969 |
| CDS | NZ_CP040530.1 | 68112..70001 | +0 | +15 | 0.9921 |
| CDS | NZ_CP040530.1 | 86707..87669 | +0 | +45 | 0.9554 |
| CDS | NZ_CP040530.1 | 92423..92833 | +0 | -93 | 0.7737 |
| CDS | NZ_CP040530.1 | 105843..106283 | +0 | +3 | 0.9932 |
| CDS | NZ_CP040530.1 | 116529..117488 | +0 | +96 | 0.9091 |
| CDS | NZ_CP040530.1 | 138225..141104 | +0 | +21 | 0.9928 |
| CDS | NZ_CP040530.1 | 144334..147288 | +0 | +21 | 0.9929 |
| CDS | NZ_CP040530.1 | 149804..152899 | +15 | +0 | 0.9952 |
| CDS | NZ_CP040530.1 | 158404..159003 | +72 | +0 | 0.8800 |
| CDS | NZ_CP040530.1 | 180163..180858 | +0 | +63 | 0.9170 |
| CDS | NZ_CP040530.1 | 183415..185427 | +3 | +0 | 0.9985 |
| CDS | NZ_CP040530.1 | 189242..189682 | +0 | +15 | 0.9671 |
| CDS | NZ_CP040530.1 | 226391..226759 | -12 | +0 | 0.9685 |
| CDS | NZ_CP040530.1 | 229283..229582 | -66 | +0 | 0.8197 |
| CDS | NZ_CP040530.1 | 248267..250099 | -6 | +0 | 0.9967 |
| CDS | NZ_CP040530.1 | 265460..265885 | -111 | +0 | 0.7933 |
| CDS | NZ_CP040530.1 | 294960..295526 | +6 | +0 | 0.9894 |
| CDS | NZ_CP040530.1 | 310727..313144 | +18 | +0 | 0.9926 |
| CDS | NZ_CP040530.1 | 313566..314054 | +0 | +3 | 0.9939 |
| CDS | NZ_CP040530.1 | 320865..322763 | +0 | -48 | 0.9747 |
| CDS | NZ_CP040530.1 | 326643..327770 | +0 | +6 | 0.9947 |
| CDS | NZ_CP040530.1 | 329317..331350 | +0 | -90 | 0.9558 |
| CDS | NZ_CP040530.1 | 359152..359955 | +0 | +84 | 0.9054 |
| CDS | NZ_CP040530.1 | 379324..381186 | +0 | +105 | 0.9466 |
| CDS | NZ_CP040530.1 | 400982..402193 | +0 | -27 | 0.9777 |
| CDS | NZ_CP040530.1 | 406919..407251 | -57 | +0 | 0.8538 |
| CDS | NZ_CP040530.1 | 438610..440562 | +0 | +6 | 0.9969 |
| CDS | NZ_CP040530.1 | 442778..443470 | -48 | +0 | 0.9352 |
| CDS | NZ_CP040530.1 | 445944..446453 | +0 | +6 | 0.9884 |
| CDS | NZ_CP040530.1 | 457108..458412 | -54 | +0 | 0.9603 |
| CDS | NZ_CP040530.1 | 462980..464071 | +21 | +0 | 0.9808 |
| CDS | NZ_CP040530.1 | 475297..476472 | -174 | +0 | 0.8711 |
| CDS | NZ_CP040530.1 | 479891..480853 | +27 | +0 | 0.9720 |
| CDS | NZ_CP040530.1 | 492198..493082 | +225 | +0 | 0.7458 |
| CDS | NZ_CP040530.1 | 538569..539486 | +0 | +12 | 0.9871 |
| CDS | NZ_CP040530.1 | 539546..540844 | +0 | -6 | 0.9954 |
| CDS | NZ_CP040530.1 | 567627..568946 | +0 | -9 | 0.9932 |
| CDS | NZ_CP040530.1 | 588270..589427 | +27 | +0 | 0.9767 |
| CDS | NZ_CP040530.1 | 593736..594017 | +0 | -117 | 0.5851 |
| CDS | NZ_CP040530.1 | 594081..595034 | +0 | -120 | 0.8742 |
| CDS | NZ_CP040530.1 | 602742..603785 | -63 | +0 | 0.9431 |
| CDS | NZ_CP040530.1 | 623156..624100 | +0 | -3 | 0.9968 |
| CDS | NZ_CP040530.1 | 645750..646691 | +0 | -15 | 0.9841 |
| CDS | NZ_CP040530.1 | 646645..648039 | +0 | +21 | 0.9852 |
| CDS | NZ_CP040530.1 | 661723..662370 | +0 | +21 | 0.9686 |
| CDS | NZ_CP040530.1 | 676324..676593 | +0 | -12 | 0.9556 |
| CDS | NZ_CP040530.1 | 724258..724881 | +0 | +15 | 0.9765 |
| CDS | NZ_CP040530.1 | 728144..728677 | +0 | +12 | 0.9780 |
| CDS | NZ_CP040530.1 | 741684..743792 | +15 | +0 | 0.9929 |
| CDS | NZ_CP040530.1 | 768706..770031 | +0 | +60 | 0.9567 |
| CDS | NZ_CP040530.1 | 778776..779579 | +0 | +57 | 0.9338 |
| CDS | NZ_CP040530.1 | 781664..782587 | -42 | +0 | 0.9565 |
| CDS | NZ_CP040530.1 | 788885..790492 | -12 | +0 | 0.9926 |
| CDS | NZ_CP040530.1 | 818093..819214 | +3 | +0 | 0.9973 |
| CDS | NZ_CP040530.1 | 820923..824840 | +0 | +12 | 0.9969 |
| CDS | NZ_CP040530.1 | 827623..829665 | +30 | +0 | 0.9853 |
| CDS | NZ_CP040530.1 | 834861..836948 | +3 | +0 | 0.9986 |
| CDS | NZ_CP040530.1 | 840595..842445 | +12 | +0 | 0.9935 |
| CDS | NZ_CP040530.1 | 866356..869271 | +0 | +48 | 0.9838 |
| CDS | NZ_CP040530.1 | 883585..885966 | +0 | -18 | 0.9924 |
| CDS | NZ_CP040530.1 | 894511..895371 | +0 | +15 | 0.9829 |
| CDS | NZ_CP040530.1 | 895379..896320 | +0 | -12 | 0.9873 |
| CDS | NZ_CP040530.1 | 901702..902793 | +0 | +171 | 0.8646 |
| CDS | NZ_CP040530.1 | 918352..921408 | +0 | +18 | 0.9941 |
| CDS | NZ_CP040530.1 | 922515..923783 | +0 | +90 | 0.9338 |
| CDS | NZ_CP040530.1 | 979364..980482 | +0 | +6 | 0.9947 |
| CDS | NZ_CP040530.1 | 1011305..1012570 | -21 | +0 | 0.9837 |
| CDS | NZ_CP040530.1 | 1048808..1050238 | +0 | +24 | 0.9835 |
| CDS | NZ_CP040530.1 | 1050609..1054649 | +0 | -42 | 0.9896 |
| CDS | NZ_CP040530.1 | 1075003..1076388 | -63 | +0 | 0.9565 |
| CDS | NZ_CP040530.1 | 1076919..1077737 | +3 | +0 | 0.9963 |
| CDS | NZ_CP040530.1 | 1081691..1084303 | +66 | +0 | 0.9747 |
| CDS | NZ_CP040530.1 | 1099271..1100461 | +3 | +0 | 0.9975 |
| CDS | NZ_CP040530.1 | 1103815..1104390 | +0 | +36 | 0.9412 |
| CDS | NZ_CP040530.1 | 1105551..1108928 | -9 | +0 | 0.9973 |
| CDS | NZ_CP040530.1 | 1118985..1120835 | +0 | +60 | 0.9686 |
| CDS | NZ_CP040530.1 | 1174673..1174840 | +48 | +0 | 0.7143 |
| CDS | NZ_CP040530.1 | 1174915..1175790 | -75 | +0 | 0.9211 |
| CDS | NZ_CP040530.1 | 1198460..1199743 | +0 | -24 | 0.9813 |
| CDS | NZ_CP040530.1 | 1212633..1213436 | +0 | +66 | 0.9241 |
| CDS | NZ_CP040530.1 | 1220510..1222024 | +0 | +3 | 0.9980 |
| CDS | NZ_CP040530.1 | 1235209..1237206 | +0 | +132 | 0.9380 |
| CDS | NZ_CP040530.1 | 1237462..1241493 | -78 | +0 | 0.9810 |
| CDS | NZ_CP040530.1 | 1243566..1244495 | +0 | -54 | 0.9419 |
| CDS | NZ_CP040530.1 | 1247751..1250915 | +0 | +21 | 0.9934 |
| CDS | NZ_CP040530.1 | 1254918..1256432 | +0 | +105 | 0.9352 |
| CDS | NZ_CP040530.1 | 1279158..1281746 | +0 | +30 | 0.9885 |
| CDS | NZ_CP040530.1 | 1313990..1315177 | +0 | -21 | 0.9823 |

- Aviso: 1 feature(s) GFF multi-segmento fusionados por ID; se comparan por su intervalo envolvente.

## DFAST

- Estado: **no comparable (not_comparable)**.
- Fuente: `executed`.
- Razón: No hay contigs en común entre referencia y herramienta tras normalizar seqid..
- Versión: `DFAST ver. 1.6.2` (estado `ok`).
- Ejecución: **ok**, 504.35 s, RSS máx. 3573.0 MiB (`2026-11-26T01:54:44Z` → `2026-11-26T02:03:08Z`).
- Salida: `/tmp/external_benchmark_dfast_dmybgdsd/out/genome.gff` (10599890 bytes, SHA256 `2a45426296fa0561…`).
- Comando: `/usr/local/bin/dfast --genome /tmp/benchmark_genome_bdz3sme8.fna --out /tmp/external_benchmark_dfast_dmybgdsd/out --cpu 4 --force --dbroot /data/dfast_db`.
- Anotación comparada: `reports/external_benchmark/heldout/btheta/evidence/dfast/genome.gff` (formato `gff3`, SHA256 `2a45426296fa0561…`).

### Features (DFAST)

| tipo | ref | query | TP | FP | FN | precisión | recall | F1 | exactos | borde distinto | solape medio |
|---|---|---|---|---|---|---|---|---|---|---|---|
| — | — | — | — | — | — | — | — | — | — | — | — |

### Atributos (DFAST)

| tipo | pares | product en ambos | product literal | product tokens | gene literal | locus literal |
|---|---|---|---|---|---|---|
| — | — | — | — | — | — | — |

### Hipotéticas (DFAST)

| tipo | ref CDS | ref con product | ref hyp | % ref | query CDS | query con product | query hyp | % query |
|---|---|---|---|---|---|---|---|---|
| — | — | — | — | — | — | — | — | — |

#### Diferencias

- Falsos negativos (referencia sin pareja): **0**.
- Falsos positivos (herramienta sin pareja): **0**.
- Pares solo con hebra distinta: **0**.
- Pares emparejados con algún borde distinto: **0**.

- Aviso: No hay contigs en común entre referencia y herramienta tras normalizar seqid.

## Herramientas no disponibles

- Ninguna herramienta seleccionada quedó `not_available`.

## Limitaciones

- La referencia es el gold de esta comparación; sus criterios y errores se heredan.
- El emparejamiento es 1:1 voraz por solapamiento recíproco y no un óptimo bipartito.
- La coincidencia literal de 'product' no mide función; los nombres dependen de cada base.
- No se comparan ncRNA, tmRNA, pseudogenes (si se excluyen), dominios ni rutas.
- Una herramienta ausente no es un cero: queda como not_available sin métricas.
- La Bakta DB se detecta localmente; no se descarga y su ausencia deja Bakta en not_available.
- No se declara superioridad entre herramientas sin un gold independiente y ejecuciones pareadas.
- Las ubicaciones compuestas de GenBank se comparan por su intervalo envolvente.

## Reproducción

Comando exacto con `--run` (Bakta requiere una Bakta DB local; no se descarga):

```bash
/opt/conda/bin/python benchmarks/external/external_benchmark.py --genome data/reference_genomes/btheta/btheta_genomic.fna.gz --reference data/reference_genomes/btheta/btheta_genomic.gff.gz --tools bakta,prokka,dfast --run --threads 1 --json reports/external_benchmark/benchmark.json --md reports/external_benchmark/benchmark.md --bakta-db /data/bakta_db
```

Para comparar salidas producidas en otro host, usar `--query-gff TOOL=RUTA` sin `--run`; el origen queda como `precomputed`.
