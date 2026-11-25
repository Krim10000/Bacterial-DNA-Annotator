# Benchmark externo CPU: Bakta, Prokka y DFAST

**Generado:** 2026-11-25T23:34:26Z

**Esquema:** `external-benchmark-run/v1`

**Alcance:** Comparación local de CDS/tRNA/rRNA de Bakta, Prokka y DFAST contra una referencia GFF3/GBFF. CPU, sin descargas de datos.

**Comando reproducible:**

```bash
/opt/conda/bin/python /opt/benchmark/external/external_benchmark.py --genome /work/data/reference_genomes/paeruginosa/paeruginosa_genomic.fna.gz --reference /work/data/reference_genomes/paeruginosa/paeruginosa_genomic.gff.gz --tools all --run --threads 4 --bakta-db /data/bakta_db --extra-arg dfast=--dbroot --extra-arg dfast=/data/dfast_db --save-gffs /work/reports/external_benchmark/multi_v2/paeruginosa/evidence --json /work/reports/external_benchmark/multi_v2/paeruginosa.json --md /work/reports/external_benchmark/multi_v2/paeruginosa.md
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

- FASTA: `data/reference_genomes/paeruginosa/paeruginosa_genomic.fna.gz` (SHA256 `9433cb94b3f21fa4…`, 1 contig, 6264404 bp).
- Referencia: `data/reference_genomes/paeruginosa/paeruginosa_genomic.gff.gz` (formato `gff3`, SHA256 `e6b8bc796af95342…`).
- Conteos de referencia: CDS **5572**, tRNA **63**, rRNA **13**.
- CDS pseudo excluidos en referencia: **0**.
- SHA256 FASTA: `9433cb94b3f21fa40e94d0d9454996fcca39da1d49074b9b22b6a459e63ca510`.
- SHA256 referencia: `e6b8bc796af95342acd5c170ecf5003e75135c894b1c22f8c4b46c51c4bdf7f0`.

## Validación

- Inicio: `2026-11-25T22:53:02Z`; fin: `2026-11-25T23:34:26Z`; total: **2484.34 s**.
- JSON declarado: `reports/external_benchmark/multi_v2/paeruginosa.json`.
- Markdown declarado: `reports/external_benchmark/multi_v2/paeruginosa.md`.

| check | estado | detalle |
|---|---|---|
| genome_path | ok | FASTA existe y es legible |
| genome_sha256 | ok | SHA256 9433cb94b3f21fa40e94d0d9454996fcca39da1d49074b9b22b6a459e63ca510 |
| reference_path | ok | referencia gff3 con 5648 features |
| reference_sha256 | ok | SHA256 e6b8bc796af95342acd5c170ecf5003e75135c894b1c22f8c4b46c51c4bdf7f0 |
| reference_features | ok | CDS=5572 tRNA=63 rRNA=13 pseudo_excluidos=0 |
| tool_bakta | ok | binario ejecutable: /usr/local/bin/bakta |
| tool_prokka | ok | binario ejecutable: /usr/local/bin/prokka |
| tool_dfast | ok | binario ejecutable: /usr/local/bin/dfast |
| bakta_db | ok | estado=ok; /data/bakta_db |
| result_bakta | ok | estado=ok; sin incidencias |
| result_prokka | ok | estado=ok; sin incidencias |
| result_dfast | ok | estado=ok; sin incidencias |

## Resultados

| herramienta | fuente | estado | F1 CDS | F1 tRNA | F1 rRNA | precisión global | recall global | F1 global | segundos |
|---|---|---|---|---|---|---|---|---|---|
| Bakta | executed | ok | 0.9847 | 1.0000 | 0.9600 | 0.9757 | 0.9942 | 0.9848 | 963.7 |
| Prokka | executed | ok | 0.9837 | 0.9265 | 0.9600 | 0.9738 | 0.9924 | 0.9830 | 1055.7 |
| DFAST | executed | ok | 0.9845 | 0.9474 | 0.9600 | 0.9745 | 0.9938 | 0.9840 | 404.6 |

La precisión/recall/F1 se calculan contra la referencia local. Si una herramienta no está instalada, su estado es `not_available` y no hay métricas que mostrar.

## Bakta

- Estado: **ok**.
- Fuente: `executed`.
- Bakta DB: `ok` (válida: **sí**, ruta `/data/bakta_db`, versión `0`).
- Versión: `bakta 1.12.2` (estado `ok`).
- Ejecución: **ok**, 963.71 s, RSS máx. 2676.2 MiB (`2026-11-25T22:53:05Z` → `2026-11-25T23:09:09Z`).
- Salida: `/tmp/external_benchmark_bakta_08kxy81f/out/benchmark.gff3` (7617253 bytes, SHA256 `d9b2c702fa6cabb7…`).
- Comando: `/usr/local/bin/bakta --output /tmp/external_benchmark_bakta_08kxy81f/out --prefix benchmark --threads 4 --force --skip-plot --db /data/bakta_db /tmp/benchmark_genome_ucxxuz21.fna`.
- Anotación comparada: `reports/external_benchmark/multi_v2/paeruginosa/evidence/bakta/benchmark.gff3` (formato `gff3`, SHA256 `d9b2c702fa6cabb7…`).

### Features (Bakta)

| tipo | ref | query | TP | FP | FN | precisión | recall | F1 | exactos | borde distinto | solape medio |
|---|---|---|---|---|---|---|---|---|---|---|---|
| CDS | 5572 | 5680 | 5540 | 140 | 32 | 0.9754 | 0.9943 | 0.9847 | 5096 | 444 | 0.9932 |
| tRNA | 63 | 63 | 63 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 16 | 47 | 0.9899 |
| rRNA | 13 | 12 | 12 | 0 | 1 | 1.0000 | 0.9231 | 0.9600 | 2 | 10 | 0.9886 |

### Atributos (Bakta)

| tipo | pares | product en ambos | product literal | product tokens | gene literal | locus literal |
|---|---|---|---|---|---|---|
| CDS | 5540 | 5540 | 0.1957 | 0.3686 | 0.8734 | 0.0000 |
| tRNA | 63 | 63 | 0.0000 | 0.6468 | — | 0.0000 |
| rRNA | 12 | 12 | 1.0000 | 1.0000 | — | 0.0000 |

### Hipotéticas (Bakta)

| tipo | ref CDS | ref con product | ref hyp | % ref | query CDS | query con product | query hyp | % query |
|---|---|---|---|---|---|---|---|---|
| CDS | 5572 | 5572 | 2255 | 40.5 | 5680 | 5680 | 683 | 12.0 |
| tRNA | 0 | 0 | 0 | — | 0 | 0 | 0 | — |
| rRNA | 0 | 0 | 0 | — | 0 | 0 | 0 | — |

#### Diferencias

- Falsos negativos (referencia sin pareja): **33**.
- Falsos positivos (herramienta sin pareja): **140**.
- Pares solo con hebra distinta: **5**.
- Pares emparejados con algún borde distinto: **501**.

Las listas de diferencias están truncadas a 100 elementos por lista; los conteos son completos.

> Primeras diferencias no emparejadas

| tipo | feature | seqid | coordenadas | hebra | product/locus |
|---|---|---|---|---|---|
| FN | CDS | NC_002516.2 | 153696..153836 | + | hypothetical protein |
| FN | CDS | NC_002516.2 | 496362..496478 | - | hypothetical protein |
| FN | CDS | NC_002516.2 | 680893..681351 | + | hypothetical protein |
| FN | CDS | NC_002516.2 | 789654..790031 | + | hypothetical protein |
| FN | CDS | NC_002516.2 | 1065970..1066296 | + | colicin immunity protein |
| FN | CDS | NC_002516.2 | 1245665..1245928 | + | pyocin-S2 immunity protein |
| FN | CDS | NC_002516.2 | 1444720..1445271 | - | hypothetical protein |
| FN | CDS | NC_002516.2 | 1462696..1463403 | + | hypothetical protein |
| FN | CDS | NC_002516.2 | 1498338..1498829 | + | hypothetical protein |
| FN | CDS | NC_002516.2 | 1775849..1776034 | + | potassium-transporting ATPase subunit F |
| FN | CDS | NC_002516.2 | 1844144..1845253 | - | translocation protein in type III secretion |
| FN | CDS | NC_002516.2 | 2059569..2060552 | + | hypothetical protein |
| FN | CDS | NC_002516.2 | 2448033..2448533 | - | hypothetical protein |
| FN | CDS | NC_002516.2 | 2760872..2761351 | - | hypothetical protein |
| FN | CDS | NC_002516.2 | 2803345..2803623 | + | hypothetical protein |
| FN | CDS | NC_002516.2 | 2904566..2904907 | - | hypothetical protein |
| FN | CDS | NC_002516.2 | 3017818..3018804 | - | hypothetical protein |
| FN | CDS | NC_002516.2 | 3206915..3207166 | + | outer membrane lipoprotein OprI |
| FN | CDS | NC_002516.2 | 3290694..3291173 | - | hypothetical protein |
| FN | CDS | NC_002516.2 | 3414401..3414613 | + | ribosome modulation factor |
| FN | CDS | NC_002516.2 | 3527670..3528212 | - | hypothetical protein |
| FN | CDS | NC_002516.2 | 3535971..3536579 | - | imidazole glycerol phosphate synthase subunit HisH |
| FN | CDS | NC_002516.2 | 3684163..3684717 | - | hypothetical protein |
| FN | CDS | NC_002516.2 | 3684905..3685762 | - | hypothetical protein |
| FN | CDS | NC_002516.2 | 4008938..4009051 | - | hypothetical protein |
| FN | CDS | NC_002516.2 | 4165719..4165880 | - | MexR antirepressor ArmR |
| FN | CDS | NC_002516.2 | 4327395..4327559 | + | hypothetical protein |
| FN | CDS | NC_002516.2 | 4783771..4784139 | - | preprotein translocase subunit SecE |
| FN | CDS | NC_002516.2 | 4814213..4814899 | + | hypothetical protein |
| FN | CDS | NC_002516.2 | 5069082..5069531 | - | type 4 fimbrial protein PilA |
| FN | CDS | NC_002516.2 | 5206720..5207007 | - | hypothetical protein |
| FN | CDS | NC_002516.2 | 5332354..5332743 | - | preprotein translocase subunit SecG |
| FN | rRNA | NC_002516.2 | 1668840..1669799 | + | 4.5S ribosomal RNA |
| FP | CDS | NC_002516.2 | 56169..56447 | + | Transposase |
| FP | CDS | NC_002516.2 | 65541..65678 | + | Lipoprotein |
| FP | CDS | NC_002516.2 | 68883..69053 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 69056..69190 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 69272..69526 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 153473..154357 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 352227..352316 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 375511..375978 | + | DUF2269 domain-containing protein |
| FP | CDS | NC_002516.2 | 413364..413657 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 433933..434103 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 527176..527361 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 588695..588931 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 628213..628338 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 674046..674216 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 699125..699265 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 740138..740461 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 788927..789142 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 789144..789356 | + | DUF2188 domain-containing protein |
| FP | CDS | NC_002516.2 | 794268..794369 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 805596..805784 | + | Lipoprotein |
| FP | CDS | NC_002516.2 | 816549..817025 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 828855..829307 | + | Toxin CptA |
| FP | CDS | NC_002516.2 | 944467..944640 | - | DUF1427 domain-containing protein |
| FP | CDS | NC_002516.2 | 1029541..1029732 | - | Thiolase |
| FP | CDS | NC_002516.2 | 1064340..1064558 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 1072125..1072766 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 1099361..1099672 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 1116213..1116500 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 1215284..1215727 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 1215856..1215987 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 1340734..1340889 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 1376094..1376330 | - | (2Fe-2S)-binding protein |
| FP | CDS | NC_002516.2 | 1444777..1445508 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 1463065..1463403 | + | Secreted protein |
| FP | CDS | NC_002516.2 | 1479788..1480048 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 1498346..1499005 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 1697188..1697919 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 1697865..1698344 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 1775945..1776034 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 1835118..1835309 | + | Cardiolipin synthase N-terminal domain-containing protein |
| FP | CDS | NC_002516.2 | 1844144..1844542 | - | Type III secretion system needle length determinant |
| FP | CDS | NC_002516.2 | 1878447..1878548 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 1948389..1948505 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 1997076..1997192 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 2017863..2017970 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 2059713..2060000 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 2060169..2060552 | + | Lipoprotein |
| FP | CDS | NC_002516.2 | 2080476..2080664 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 2094139..2094285 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 2117030..2118496 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 2119751..2120179 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 2165064..2165204 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 2207951..2208169 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 2238860..2239420 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 2281755..2281964 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 2351431..2351880 | - | YetF C-terminal domain-containing protein |
| FP | CDS | NC_002516.2 | 2355684..2355806 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 2361293..2361457 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 2364180..2364314 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 2400441..2400671 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 2406961..2407131 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 2412082..2412348 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 2472077..2472370 | - | Inner centromere protein |
| FP | CDS | NC_002516.2 | 2473043..2473171 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 2562140..2562283 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 2725768..2725977 | - | Protein AaeX |
| FP | CDS | NC_002516.2 | 2752866..2753465 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 2753519..2754445 | - | enterochelin esterase domain-containing protein |
| FP | CDS | NC_002516.2 | 2803504..2803623 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 2845916..2846170 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 2881520..2881711 | + | DUF2188 domain-containing protein |
| FP | CDS | NC_002516.2 | 2893100..2893246 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 2901752..2901961 | - | Transmembrane protein |
| FP | CDS | NC_002516.2 | 2904566..2904712 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 3010833..3011099 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 3016675..3016884 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 3017818..3018261 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 3031561..3031752 | - | Immune protein Tsi5 |
| FP | CDS | NC_002516.2 | 3079850..3083482 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 3083485..3086145 | + | AAA domain-containing protein |
| FP | CDS | NC_002516.2 | 3087813..3088037 | + | Toxic protein SymE |
| FP | CDS | NC_002516.2 | 3089644..3090249 | - | Restriction endonuclease type IV Mrr domain-containing protein |
| FP | CDS | NC_002516.2 | 3108663..3108752 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 3110522..3110743 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 3112422..3112655 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 3113767..3114084 | - | Secreted protein |
| FP | CDS | NC_002516.2 | 3123059..3123268 | - | Cation-transporting P-type ATPase C-terminal domain-containing protein |
| FP | CDS | NC_002516.2 | 3125429..3125524 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 3129158..3129424 | - | Plasmid replication protein RepB |
| FP | CDS | NC_002516.2 | 3147641..3147733 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 3150886..3152202 | - | exo-alpha-sialidase |
| FP | CDS | NC_002516.2 | 3218444..3218572 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 3290694..3290825 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 3290847..3291119 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 3306112..3306294 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 3347768..3348076 | + | Phosphodiesterase |
| FP | CDS | NC_002516.2 | 3550746..3552986 | - | 3-phosphoshikimate 1-carboxyvinyltransferase |
| FP | CDS | NC_002516.2 | 3588444..3588656 | - | Fe-S protein |
| FP | CDS | NC_002516.2 | 3607065..3608120 | - | Glycosyltransferase family 1 protein |
| FP | CDS | NC_002516.2 | 3647387..3647746 | - | GGDEF domain-containing protein |

> Primeros pares con coordenadas no exactas

| tipo | seqid | ref | Δstart | Δend | solape |
|---|---|---|---|---|---|
| CDS | NC_002516.2 | 10434..12488 | +0 | +3 | 0.9985 |
| CDS | NC_002516.2 | 29896..30498 | +219 | +0 | 0.6368 |
| CDS | NC_002516.2 | 34624..35538 | +15 | +0 | 0.9836 |
| CDS | NC_002516.2 | 40589..40816 | +6 | +0 | 0.9737 |
| CDS | NC_002516.2 | 57212..58594 | +0 | -129 | 0.9067 |
| CDS | NC_002516.2 | 61879..62388 | -156 | +0 | 0.7658 |
| CDS | NC_002516.2 | 73468..73923 | +30 | +0 | 0.9342 |
| CDS | NC_002516.2 | 91521..94826 | +0 | +225 | 0.9363 |
| CDS | NC_002516.2 | 104006..104851 | +60 | +0 | 0.9291 |
| CDS | NC_002516.2 | 110976..112907 | -9 | +0 | 0.9954 |
| CDS | NC_002516.2 | 117931..119130 | +87 | +0 | 0.9275 |
| CDS | NC_002516.2 | 127378..128502 | +24 | +0 | 0.9787 |
| CDS | NC_002516.2 | 138818..140167 | +15 | +0 | 0.9889 |
| CDS | NC_002516.2 | 151936..153306 | +0 | -72 | 0.9475 |
| CDS | NC_002516.2 | 160552..161448 | +0 | +177 | 0.8352 |
| CDS | NC_002516.2 | 169361..169906 | -87 | +0 | 0.8626 |
| CDS | NC_002516.2 | 173811..174638 | -90 | +0 | 0.9020 |
| CDS | NC_002516.2 | 177307..178458 | +48 | +0 | 0.9583 |
| CDS | NC_002516.2 | 192965..193672 | +0 | +120 | 0.8551 |
| CDS | NC_002516.2 | 283818..284105 | +0 | -33 | 0.8854 |
| CDS | NC_002516.2 | 284197..284673 | +0 | -12 | 0.9748 |
| CDS | NC_002516.2 | 286338..287024 | +0 | +15 | 0.9786 |
| CDS | NC_002516.2 | 342292..343173 | -30 | +0 | 0.9671 |
| CDS | NC_002516.2 | 348440..349051 | +0 | -36 | 0.9412 |
| CDS | NC_002516.2 | 349050..350069 | +78 | +0 | 0.9235 |
| CDS | NC_002516.2 | 350890..351588 | +69 | +0 | 0.9013 |
| CDS | NC_002516.2 | 358931..359920 | +57 | +0 | 0.9424 |
| CDS | NC_002516.2 | 367457..368422 | +36 | +0 | 0.9627 |
| CDS | NC_002516.2 | 368462..370405 | +0 | +33 | 0.9833 |
| CDS | NC_002516.2 | 380903..382033 | +93 | +0 | 0.9178 |
| CDS | NC_002516.2 | 400248..401072 | -54 | +0 | 0.9386 |
| CDS | NC_002516.2 | 416009..417406 | +0 | -27 | 0.9807 |
| CDS | NC_002516.2 | 424670..425341 | -63 | +0 | 0.9143 |
| CDS | NC_002516.2 | 447391..448350 | +0 | -66 | 0.9313 |
| CDS | NC_002516.2 | 464568..465983 | +0 | -42 | 0.9703 |
| CDS | NC_002516.2 | 482111..482665 | +0 | -15 | 0.9730 |
| CDS | NC_002516.2 | 484404..484838 | +18 | +0 | 0.9586 |
| CDS | NC_002516.2 | 492080..493357 | +0 | -3 | 0.9977 |
| CDS | NC_002516.2 | 505047..505586 | +135 | +0 | 0.7500 |
| CDS | NC_002516.2 | 508959..509753 | +3 | +0 | 0.9962 |
| CDS | NC_002516.2 | 518083..520635 | -246 | +0 | 0.9121 |
| CDS | NC_002516.2 | 535085..535489 | +3 | +0 | 0.9926 |
| CDS | NC_002516.2 | 538217..539143 | +0 | +39 | 0.9596 |
| CDS | NC_002516.2 | 550867..551793 | +0 | -9 | 0.9903 |
| CDS | NC_002516.2 | 556275..557276 | +0 | +81 | 0.9252 |
| CDS | NC_002516.2 | 558377..559078 | +0 | -33 | 0.9530 |
| CDS | NC_002516.2 | 581133..581660 | -45 | +0 | 0.9215 |
| CDS | NC_002516.2 | 590896..591402 | +0 | -243 | 0.5207 |
| CDS | NC_002516.2 | 596970..597479 | +0 | +72 | 0.8763 |
| CDS | NC_002516.2 | 614952..615608 | +0 | -231 | 0.6484 |
| CDS | NC_002516.2 | 624994..627765 | +0 | -45 | 0.9838 |
| CDS | NC_002516.2 | 641073..641426 | -63 | +0 | 0.8489 |
| CDS | NC_002516.2 | 675390..675839 | +96 | +0 | 0.7867 |
| CDS | NC_002516.2 | 716731..717153 | +9 | +0 | 0.9787 |
| CDS | NC_002516.2 | 718729..720072 | +0 | +96 | 0.9333 |
| CDS | NC_002516.2 | 727608..730703 | +0 | -21 | 0.9932 |
| CDS | NC_002516.2 | 740203..741348 | +261 | +0 | 0.7723 |
| CDS | NC_002516.2 | 770156..770818 | -15 | +0 | 0.9779 |
| CDS | NC_002516.2 | 784698..785174 | +54 | +0 | 0.8868 |
| CDS | NC_002516.2 | 785969..786925 | +105 | +0 | 0.8903 |
| CDS | NC_002516.2 | 791370..792632 | +186 | +0 | 0.8527 |
| CDS | NC_002516.2 | 801967..802239 | +0 | -60 | 0.7802 |
| CDS | NC_002516.2 | 804381..804560 | +6 | +0 | 0.9667 |
| CDS | NC_002516.2 | 842878..843618 | +0 | +18 | 0.9763 |
| CDS | NC_002516.2 | 856942..857871 | +0 | -12 | 0.9871 |
| CDS | NC_002516.2 | 860175..863300 | +0 | -132 | 0.9578 |
| CDS | NC_002516.2 | 863455..863664 | +6 | +0 | 0.9714 |
| CDS | NC_002516.2 | 875953..876603 | +39 | +0 | 0.9401 |
| CDS | NC_002516.2 | 876617..878608 | +0 | +330 | 0.8579 |
| CDS | NC_002516.2 | 879027..879560 | +237 | +0 | 0.5562 |
| CDS | NC_002516.2 | 899138..899656 | +0 | -153 | 0.7052 |
| CDS | NC_002516.2 | 925007..927214 | +0 | -222 | 0.8995 |
| CDS | NC_002516.2 | 935427..935975 | +45 | +0 | 0.9180 |
| CDS | NC_002516.2 | 937437..938039 | -36 | +0 | 0.9437 |
| CDS | NC_002516.2 | 977420..977743 | +0 | -63 | 0.8056 |
| CDS | NC_002516.2 | 981422..982888 | +3 | +0 | 0.9980 |
| CDS | NC_002516.2 | 1016657..1017994 | -3 | +0 | 0.9978 |
| CDS | NC_002516.2 | 1039044..1039229 | -9 | +0 | 0.9538 |
| CDS | NC_002516.2 | 1045832..1046461 | +201 | +0 | 0.6810 |
| CDS | NC_002516.2 | 1058812..1059606 | +147 | +0 | 0.8151 |
| CDS | NC_002516.2 | 1061207..1062034 | +0 | -90 | 0.8913 |
| CDS | NC_002516.2 | 1071239..1071877 | +0 | +66 | 0.9064 |
| CDS | NC_002516.2 | 1083972..1085564 | +21 | +0 | 0.9868 |
| CDS | NC_002516.2 | 1097296..1098288 | +93 | +0 | 0.9063 |
| CDS | NC_002516.2 | 1105888..1107000 | +48 | +0 | 0.9569 |
| CDS | NC_002516.2 | 1113129..1114718 | +96 | +0 | 0.9396 |
| CDS | NC_002516.2 | 1118158..1119519 | -114 | +0 | 0.9228 |
| CDS | NC_002516.2 | 1131674..1133947 | +24 | +0 | 0.9894 |
| CDS | NC_002516.2 | 1134024..1135202 | +33 | +0 | 0.9720 |
| CDS | NC_002516.2 | 1151435..1151908 | +75 | +0 | 0.8418 |
| CDS | NC_002516.2 | 1151919..1152575 | +0 | -57 | 0.9132 |
| CDS | NC_002516.2 | 1153637..1155547 | +69 | +0 | 0.9639 |
| CDS | NC_002516.2 | 1211888..1212571 | +0 | -129 | 0.8114 |
| CDS | NC_002516.2 | 1224641..1225330 | +18 | +0 | 0.9739 |
| CDS | NC_002516.2 | 1225757..1226278 | +0 | -108 | 0.7931 |
| CDS | NC_002516.2 | 1227474..1228205 | -57 | +0 | 0.9278 |
| CDS | NC_002516.2 | 1228212..1229198 | +0 | +30 | 0.9705 |
| CDS | NC_002516.2 | 1256094..1257452 | -258 | +0 | 0.8404 |
| CDS | NC_002516.2 | 1299558..1300175 | +6 | +0 | 0.9903 |
| CDS | NC_002516.2 | 1315916..1317202 | +0 | +192 | 0.8702 |

- Aviso: 1 feature(s) GFF multi-segmento fusionados por ID; se comparan por su intervalo envolvente.

## Prokka

- Estado: **ok**.
- Fuente: `executed`.
- Versión: `prokka 1.15.6` (estado `ok`).
- Ejecución: **ok**, 1055.68 s, RSS máx. 3217.3 MiB (`2026-11-25T23:09:26Z` → `2026-11-25T23:27:02Z`).
- Salida: `/tmp/external_benchmark_prokka_absn2p97/out/benchmark.gff` (7691018 bytes, SHA256 `6fc7728a704f6d6a…`).
- Comando: `/usr/local/bin/prokka --outdir /tmp/external_benchmark_prokka_absn2p97/out --prefix benchmark --cpus 4 --force --quiet /tmp/benchmark_genome_ucxxuz21.fna`.
- Anotación comparada: `reports/external_benchmark/multi_v2/paeruginosa/evidence/prokka/benchmark.gff` (formato `gff3`, SHA256 `6fc7728a704f6d6a…`).

### Features (Prokka)

| tipo | ref | query | TP | FP | FN | precisión | recall | F1 | exactos | borde distinto | solape medio |
|---|---|---|---|---|---|---|---|---|---|---|---|
| CDS | 5572 | 5671 | 5530 | 141 | 42 | 0.9751 | 0.9925 | 0.9837 | 5086 | 444 | 0.9931 |
| tRNA | 63 | 73 | 63 | 10 | 0 | 0.8630 | 1.0000 | 0.9265 | 16 | 47 | 0.9899 |
| rRNA | 13 | 12 | 12 | 0 | 1 | 1.0000 | 0.9231 | 0.9600 | 2 | 10 | 0.9886 |

### Atributos (Prokka)

| tipo | pares | product en ambos | product literal | product tokens | gene literal | locus literal |
|---|---|---|---|---|---|---|
| CDS | 5530 | 5530 | 0.3886 | 0.5211 | 0.6004 | 0.0000 |
| tRNA | 63 | 63 | 0.0000 | 0.6667 | — | 0.0000 |
| rRNA | 12 | 12 | 1.0000 | 1.0000 | — | 0.0000 |

### Hipotéticas (Prokka)

| tipo | ref CDS | ref con product | ref hyp | % ref | query CDS | query con product | query hyp | % query |
|---|---|---|---|---|---|---|---|---|
| CDS | 5572 | 5572 | 2255 | 40.5 | 5671 | 5671 | 1994 | 35.2 |
| tRNA | 0 | 0 | 0 | — | 0 | 0 | 0 | — |
| rRNA | 0 | 0 | 0 | — | 0 | 0 | 0 | — |

#### Diferencias

- Falsos negativos (referencia sin pareja): **43**.
- Falsos positivos (herramienta sin pareja): **151**.
- Pares solo con hebra distinta: **5**.
- Pares emparejados con algún borde distinto: **501**.

Las listas de diferencias están truncadas a 100 elementos por lista; los conteos son completos.

> Primeras diferencias no emparejadas

| tipo | feature | seqid | coordenadas | hebra | product/locus |
|---|---|---|---|---|---|
| FN | CDS | NC_002516.2 | 153696..153836 | + | hypothetical protein |
| FN | CDS | NC_002516.2 | 496362..496478 | - | hypothetical protein |
| FN | CDS | NC_002516.2 | 680893..681351 | + | hypothetical protein |
| FN | CDS | NC_002516.2 | 789654..790031 | + | hypothetical protein |
| FN | CDS | NC_002516.2 | 1065970..1066296 | + | colicin immunity protein |
| FN | CDS | NC_002516.2 | 1245665..1245928 | + | pyocin-S2 immunity protein |
| FN | CDS | NC_002516.2 | 1399231..1400505 | - | hypothetical protein |
| FN | CDS | NC_002516.2 | 1444720..1445271 | - | hypothetical protein |
| FN | CDS | NC_002516.2 | 1462696..1463403 | + | hypothetical protein |
| FN | CDS | NC_002516.2 | 1498338..1498829 | + | hypothetical protein |
| FN | CDS | NC_002516.2 | 1775849..1776034 | + | potassium-transporting ATPase subunit F |
| FN | CDS | NC_002516.2 | 1844144..1845253 | - | translocation protein in type III secretion |
| FN | CDS | NC_002516.2 | 1923044..1925770 | - | assimilatory nitrate reductase |
| FN | CDS | NC_002516.2 | 2059569..2060552 | + | hypothetical protein |
| FN | CDS | NC_002516.2 | 2099452..2103297 | + | cobaltochelatase subunit CobN |
| FN | CDS | NC_002516.2 | 2448033..2448533 | - | hypothetical protein |
| FN | CDS | NC_002516.2 | 2760872..2761351 | - | hypothetical protein |
| FN | CDS | NC_002516.2 | 2803345..2803623 | + | hypothetical protein |
| FN | CDS | NC_002516.2 | 2904566..2904907 | - | hypothetical protein |
| FN | CDS | NC_002516.2 | 3017818..3018804 | - | hypothetical protein |
| FN | CDS | NC_002516.2 | 3057902..3058921 | - | transcriptional regulator |
| FN | CDS | NC_002516.2 | 3195900..3196691 | + | enoyl-CoA hydratase |
| FN | CDS | NC_002516.2 | 3206915..3207166 | + | outer membrane lipoprotein OprI |
| FN | CDS | NC_002516.2 | 3290694..3291173 | - | hypothetical protein |
| FN | CDS | NC_002516.2 | 3392268..3393485 | - | molybdenum cofactor biosynthesis protein A |
| FN | CDS | NC_002516.2 | 3414401..3414613 | + | ribosome modulation factor |
| FN | CDS | NC_002516.2 | 3527670..3528212 | - | hypothetical protein |
| FN | CDS | NC_002516.2 | 3535971..3536579 | - | imidazole glycerol phosphate synthase subunit HisH |
| FN | CDS | NC_002516.2 | 3684163..3684717 | - | hypothetical protein |
| FN | CDS | NC_002516.2 | 3684905..3685762 | - | hypothetical protein |
| FN | CDS | NC_002516.2 | 3764476..3765087 | - | hypothetical protein |
| FN | CDS | NC_002516.2 | 4008938..4009051 | - | hypothetical protein |
| FN | CDS | NC_002516.2 | 4062426..4062899 | - | 2-C-methyl-D-erythritol 2,4-cyclodiphosphate synthase |
| FN | CDS | NC_002516.2 | 4165719..4165880 | - | MexR antirepressor ArmR |
| FN | CDS | NC_002516.2 | 4327395..4327559 | + | hypothetical protein |
| FN | CDS | NC_002516.2 | 4343540..4344835 | - | nitrite extrusion protein 1 |
| FN | CDS | NC_002516.2 | 4783771..4784139 | - | preprotein translocase subunit SecE |
| FN | CDS | NC_002516.2 | 4814213..4814899 | + | hypothetical protein |
| FN | CDS | NC_002516.2 | 5069082..5069531 | - | type 4 fimbrial protein PilA |
| FN | CDS | NC_002516.2 | 5206720..5207007 | - | hypothetical protein |
| FN | CDS | NC_002516.2 | 5332354..5332743 | - | preprotein translocase subunit SecG |
| FN | CDS | NC_002516.2 | 6204444..6206210 | - | sensor histidine kinase MifS |
| FN | rRNA | NC_002516.2 | 1668840..1669799 | + | 4.5S ribosomal RNA |
| FP | CDS | NC_002516.2 | 56169..56447 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 65541..65678 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 68883..69053 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 69056..69190 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 69272..69526 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 153473..154357 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 352227..352316 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 375511..375978 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 413364..413657 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 433933..434103 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 527176..527361 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 588695..588931 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 628213..628338 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 674046..674216 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 699125..699265 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 740138..740461 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 788927..789142 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 789144..789356 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 794268..794369 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 805596..805784 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 816549..817025 | + | HTH-type transcriptional activator RhaR |
| FP | CDS | NC_002516.2 | 828855..829307 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 944467..944640 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 1029541..1029732 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 1064340..1064558 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 1072125..1072766 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 1099361..1099672 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 1116213..1116500 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 1215284..1215727 | + | Peptide deformylase |
| FP | CDS | NC_002516.2 | 1215856..1215987 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 1340734..1340889 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 1376094..1376330 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 1444777..1445508 | + | Trans-aconitate 2-methyltransferase |
| FP | CDS | NC_002516.2 | 1463065..1463403 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 1479788..1480048 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 1498346..1499005 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 1697188..1697919 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 1697865..1698344 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 1775945..1776034 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 1835118..1835309 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 1844144..1844542 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 1878447..1878548 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 1948389..1948505 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 1997076..1997192 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 2017863..2017970 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 2059713..2060000 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 2060169..2060552 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 2080476..2080664 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 2094139..2094285 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 2117030..2118496 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 2119751..2120179 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 2165064..2165204 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 2207951..2208169 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 2238860..2239420 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 2281755..2281964 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 2351431..2351880 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 2355684..2355806 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 2361293..2361457 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 2364180..2364314 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 2400441..2400671 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 2406961..2407131 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 2412082..2412348 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 2472077..2472370 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 2473043..2473171 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 2562140..2562283 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 2725768..2725977 | - | Protein AaeX |
| FP | CDS | NC_002516.2 | 2752866..2753465 | - | Iron(III) enterobactin esterase |
| FP | CDS | NC_002516.2 | 2753519..2754445 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 2803504..2803623 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 2845916..2846170 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 2881520..2881711 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 2893100..2893246 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 2901752..2901961 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 2904566..2904712 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 3010833..3011099 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 3016675..3016884 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 3017818..3018261 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 3031561..3031752 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 3079850..3083482 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 3083485..3086145 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 3087813..3088037 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 3089644..3090249 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 3108663..3108752 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 3110522..3110743 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 3112422..3112655 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 3113767..3114084 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 3123059..3123268 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 3125429..3125524 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 3129158..3129424 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 3147641..3147733 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 3150886..3152202 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 3218444..3218572 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 3290694..3290825 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 3290847..3291119 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 3306112..3306294 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 3347768..3348076 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 3550746..3552986 | - | 3-phosphoshikimate 1-carboxyvinyltransferase |
| FP | CDS | NC_002516.2 | 3588444..3588656 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 3607065..3608120 | - | GDP-mannose-dependent alpha-mannosyltransferase |
| FP | CDS | NC_002516.2 | 3647387..3647746 | - | hypothetical protein |

> Primeros pares con coordenadas no exactas

| tipo | seqid | ref | Δstart | Δend | solape |
|---|---|---|---|---|---|
| CDS | NC_002516.2 | 10434..12488 | +0 | +3 | 0.9985 |
| CDS | NC_002516.2 | 29896..30498 | +219 | +0 | 0.6368 |
| CDS | NC_002516.2 | 34624..35538 | +15 | +0 | 0.9836 |
| CDS | NC_002516.2 | 40589..40816 | +6 | +0 | 0.9737 |
| CDS | NC_002516.2 | 57212..58594 | +0 | -129 | 0.9067 |
| CDS | NC_002516.2 | 61879..62388 | -156 | +0 | 0.7658 |
| CDS | NC_002516.2 | 73468..73923 | +30 | +0 | 0.9342 |
| CDS | NC_002516.2 | 91521..94826 | +0 | +225 | 0.9363 |
| CDS | NC_002516.2 | 104006..104851 | +60 | +0 | 0.9291 |
| CDS | NC_002516.2 | 110976..112907 | -9 | +0 | 0.9954 |
| CDS | NC_002516.2 | 117931..119130 | +87 | +0 | 0.9275 |
| CDS | NC_002516.2 | 127378..128502 | +24 | +0 | 0.9787 |
| CDS | NC_002516.2 | 138818..140167 | +15 | +0 | 0.9889 |
| CDS | NC_002516.2 | 151936..153306 | +0 | -72 | 0.9475 |
| CDS | NC_002516.2 | 160552..161448 | +0 | +177 | 0.8352 |
| CDS | NC_002516.2 | 169361..169906 | -87 | +0 | 0.8626 |
| CDS | NC_002516.2 | 173811..174638 | -90 | +0 | 0.9020 |
| CDS | NC_002516.2 | 177307..178458 | +48 | +0 | 0.9583 |
| CDS | NC_002516.2 | 192965..193672 | +0 | +120 | 0.8551 |
| CDS | NC_002516.2 | 283818..284105 | +0 | -33 | 0.8854 |
| CDS | NC_002516.2 | 284197..284673 | +0 | -12 | 0.9748 |
| CDS | NC_002516.2 | 286338..287024 | +0 | +15 | 0.9786 |
| CDS | NC_002516.2 | 342292..343173 | -30 | +0 | 0.9671 |
| CDS | NC_002516.2 | 348440..349051 | +0 | -36 | 0.9412 |
| CDS | NC_002516.2 | 349050..350069 | +78 | +0 | 0.9235 |
| CDS | NC_002516.2 | 350890..351588 | +69 | +0 | 0.9013 |
| CDS | NC_002516.2 | 358931..359920 | +57 | +0 | 0.9424 |
| CDS | NC_002516.2 | 367457..368422 | +36 | +0 | 0.9627 |
| CDS | NC_002516.2 | 368462..370405 | +0 | +33 | 0.9833 |
| CDS | NC_002516.2 | 380903..382033 | +93 | +0 | 0.9178 |
| CDS | NC_002516.2 | 400248..401072 | -54 | +0 | 0.9386 |
| CDS | NC_002516.2 | 416009..417406 | +0 | -27 | 0.9807 |
| CDS | NC_002516.2 | 424670..425341 | -63 | +0 | 0.9143 |
| CDS | NC_002516.2 | 447391..448350 | +0 | -66 | 0.9313 |
| CDS | NC_002516.2 | 464568..465983 | +0 | -42 | 0.9703 |
| CDS | NC_002516.2 | 482111..482665 | +0 | -15 | 0.9730 |
| CDS | NC_002516.2 | 484404..484838 | +18 | +0 | 0.9586 |
| CDS | NC_002516.2 | 492080..493357 | +0 | -3 | 0.9977 |
| CDS | NC_002516.2 | 505047..505586 | +135 | +0 | 0.7500 |
| CDS | NC_002516.2 | 508959..509753 | +3 | +0 | 0.9962 |
| CDS | NC_002516.2 | 518083..520635 | -246 | +0 | 0.9121 |
| CDS | NC_002516.2 | 535085..535489 | +3 | +0 | 0.9926 |
| CDS | NC_002516.2 | 538217..539143 | +0 | +39 | 0.9596 |
| CDS | NC_002516.2 | 550867..551793 | +0 | -9 | 0.9903 |
| CDS | NC_002516.2 | 556275..557276 | +0 | +81 | 0.9252 |
| CDS | NC_002516.2 | 558377..559078 | +0 | -33 | 0.9530 |
| CDS | NC_002516.2 | 581133..581660 | -45 | +0 | 0.9215 |
| CDS | NC_002516.2 | 590896..591402 | +0 | -243 | 0.5207 |
| CDS | NC_002516.2 | 596970..597479 | +0 | +72 | 0.8763 |
| CDS | NC_002516.2 | 614952..615608 | +0 | -231 | 0.6484 |
| CDS | NC_002516.2 | 624994..627765 | +0 | -45 | 0.9838 |
| CDS | NC_002516.2 | 641073..641426 | -63 | +0 | 0.8489 |
| CDS | NC_002516.2 | 675390..675839 | +96 | +0 | 0.7867 |
| CDS | NC_002516.2 | 716731..717153 | +9 | +0 | 0.9787 |
| CDS | NC_002516.2 | 718729..720072 | +0 | +96 | 0.9333 |
| CDS | NC_002516.2 | 727608..730703 | +0 | -21 | 0.9932 |
| CDS | NC_002516.2 | 740203..741348 | +261 | +0 | 0.7723 |
| CDS | NC_002516.2 | 770156..770818 | -15 | +0 | 0.9779 |
| CDS | NC_002516.2 | 784698..785174 | +54 | +0 | 0.8868 |
| CDS | NC_002516.2 | 785969..786925 | +105 | +0 | 0.8903 |
| CDS | NC_002516.2 | 791370..792632 | +186 | +0 | 0.8527 |
| CDS | NC_002516.2 | 801967..802239 | +0 | -60 | 0.7802 |
| CDS | NC_002516.2 | 804381..804560 | +6 | +0 | 0.9667 |
| CDS | NC_002516.2 | 842878..843618 | +0 | +18 | 0.9763 |
| CDS | NC_002516.2 | 856942..857871 | +0 | -12 | 0.9871 |
| CDS | NC_002516.2 | 860175..863300 | +0 | -132 | 0.9578 |
| CDS | NC_002516.2 | 863455..863664 | +6 | +0 | 0.9714 |
| CDS | NC_002516.2 | 875953..876603 | +39 | +0 | 0.9401 |
| CDS | NC_002516.2 | 876617..878608 | +0 | +330 | 0.8579 |
| CDS | NC_002516.2 | 879027..879560 | +237 | +0 | 0.5562 |
| CDS | NC_002516.2 | 899138..899656 | +0 | -153 | 0.7052 |
| CDS | NC_002516.2 | 925007..927214 | +0 | -222 | 0.8995 |
| CDS | NC_002516.2 | 935427..935975 | +45 | +0 | 0.9180 |
| CDS | NC_002516.2 | 937437..938039 | -36 | +0 | 0.9437 |
| CDS | NC_002516.2 | 977420..977743 | +0 | -63 | 0.8056 |
| CDS | NC_002516.2 | 981422..982888 | +3 | +0 | 0.9980 |
| CDS | NC_002516.2 | 1016657..1017994 | -3 | +0 | 0.9978 |
| CDS | NC_002516.2 | 1039044..1039229 | -9 | +0 | 0.9538 |
| CDS | NC_002516.2 | 1045832..1046461 | +201 | +0 | 0.6810 |
| CDS | NC_002516.2 | 1058812..1059606 | +147 | +0 | 0.8151 |
| CDS | NC_002516.2 | 1061207..1062034 | +0 | -90 | 0.8913 |
| CDS | NC_002516.2 | 1071239..1071877 | +0 | +66 | 0.9064 |
| CDS | NC_002516.2 | 1083972..1085564 | +21 | +0 | 0.9868 |
| CDS | NC_002516.2 | 1097296..1098288 | +93 | +0 | 0.9063 |
| CDS | NC_002516.2 | 1105888..1107000 | +48 | +0 | 0.9569 |
| CDS | NC_002516.2 | 1113129..1114718 | +96 | +0 | 0.9396 |
| CDS | NC_002516.2 | 1118158..1119519 | -114 | +0 | 0.9228 |
| CDS | NC_002516.2 | 1131674..1133947 | +24 | +0 | 0.9894 |
| CDS | NC_002516.2 | 1134024..1135202 | +33 | +0 | 0.9720 |
| CDS | NC_002516.2 | 1151435..1151908 | +75 | +0 | 0.8418 |
| CDS | NC_002516.2 | 1151919..1152575 | +0 | -57 | 0.9132 |
| CDS | NC_002516.2 | 1153637..1155547 | +69 | +0 | 0.9639 |
| CDS | NC_002516.2 | 1211888..1212571 | +0 | -129 | 0.8114 |
| CDS | NC_002516.2 | 1224641..1225330 | +18 | +0 | 0.9739 |
| CDS | NC_002516.2 | 1225757..1226278 | +0 | -108 | 0.7931 |
| CDS | NC_002516.2 | 1227474..1228205 | -57 | +0 | 0.9278 |
| CDS | NC_002516.2 | 1228212..1229198 | +0 | +30 | 0.9705 |
| CDS | NC_002516.2 | 1256094..1257452 | -258 | +0 | 0.8404 |
| CDS | NC_002516.2 | 1299558..1300175 | +6 | +0 | 0.9903 |
| CDS | NC_002516.2 | 1315916..1317202 | +0 | +192 | 0.8702 |

- Aviso: 1 feature(s) GFF multi-segmento fusionados por ID; se comparan por su intervalo envolvente.

## DFAST

- Estado: **ok**.
- Fuente: `executed`.
- Versión: `DFAST ver. 1.6.2` (estado `ok`).
- Ejecución: **ok**, 404.64 s, RSS máx. 3217.3 MiB (`2026-11-25T23:27:22Z` → `2026-11-25T23:34:07Z`).
- Salida: `/tmp/external_benchmark_dfast_qbdgcebo/out/genome.gff` (11122636 bytes, SHA256 `ad7f40527096db67…`).
- Comando: `/usr/local/bin/dfast --genome /tmp/benchmark_genome_ucxxuz21.fna --out /tmp/external_benchmark_dfast_qbdgcebo/out --cpu 4 --force --dbroot /data/dfast_db`.
- Anotación comparada: `reports/external_benchmark/multi_v2/paeruginosa/evidence/dfast/genome.gff` (formato `gff3`, SHA256 `ad7f40527096db67…`).

### Features (DFAST)

| tipo | ref | query | TP | FP | FN | precisión | recall | F1 | exactos | borde distinto | solape medio |
|---|---|---|---|---|---|---|---|---|---|---|---|
| CDS | 5572 | 5678 | 5538 | 140 | 34 | 0.9753 | 0.9939 | 0.9845 | 5095 | 443 | 0.9932 |
| tRNA | 63 | 70 | 63 | 7 | 0 | 0.9000 | 1.0000 | 0.9474 | 16 | 47 | 0.9899 |
| rRNA | 13 | 12 | 12 | 0 | 1 | 1.0000 | 0.9231 | 0.9600 | 0 | 12 | 0.9763 |

### Atributos (DFAST)

| tipo | pares | product en ambos | product literal | product tokens | gene literal | locus literal |
|---|---|---|---|---|---|---|
| CDS | 5538 | 5538 | 0.2250 | 0.4181 | 0.8489 | 0.0000 |
| tRNA | 63 | 63 | 1.0000 | 1.0000 | — | 0.0000 |
| rRNA | 12 | 12 | 1.0000 | 1.0000 | — | 0.0000 |

### Hipotéticas (DFAST)

| tipo | ref CDS | ref con product | ref hyp | % ref | query CDS | query con product | query hyp | % query |
|---|---|---|---|---|---|---|---|---|
| CDS | 5572 | 5572 | 2255 | 40.5 | 5678 | 5678 | 791 | 13.9 |
| tRNA | 0 | 0 | 0 | — | 0 | 0 | 0 | — |
| rRNA | 0 | 0 | 0 | — | 0 | 0 | 0 | — |

#### Diferencias

- Falsos negativos (referencia sin pareja): **35**.
- Falsos positivos (herramienta sin pareja): **147**.
- Pares solo con hebra distinta: **5**.
- Pares emparejados con algún borde distinto: **502**.

Las listas de diferencias están truncadas a 100 elementos por lista; los conteos son completos.

> Primeras diferencias no emparejadas

| tipo | feature | seqid | coordenadas | hebra | product/locus |
|---|---|---|---|---|---|
| FN | CDS | NC_002516.2 | 153696..153836 | + | hypothetical protein |
| FN | CDS | NC_002516.2 | 496362..496478 | - | hypothetical protein |
| FN | CDS | NC_002516.2 | 680893..681351 | + | hypothetical protein |
| FN | CDS | NC_002516.2 | 789654..790031 | + | hypothetical protein |
| FN | CDS | NC_002516.2 | 1065970..1066296 | + | colicin immunity protein |
| FN | CDS | NC_002516.2 | 1245665..1245928 | + | pyocin-S2 immunity protein |
| FN | CDS | NC_002516.2 | 1444720..1445271 | - | hypothetical protein |
| FN | CDS | NC_002516.2 | 1462696..1463403 | + | hypothetical protein |
| FN | CDS | NC_002516.2 | 1498338..1498829 | + | hypothetical protein |
| FN | CDS | NC_002516.2 | 1775849..1776034 | + | potassium-transporting ATPase subunit F |
| FN | CDS | NC_002516.2 | 1844144..1845253 | - | translocation protein in type III secretion |
| FN | CDS | NC_002516.2 | 2059569..2060552 | + | hypothetical protein |
| FN | CDS | NC_002516.2 | 2448033..2448533 | - | hypothetical protein |
| FN | CDS | NC_002516.2 | 2760872..2761351 | - | hypothetical protein |
| FN | CDS | NC_002516.2 | 2803345..2803623 | + | hypothetical protein |
| FN | CDS | NC_002516.2 | 2904566..2904907 | - | hypothetical protein |
| FN | CDS | NC_002516.2 | 3017818..3018804 | - | hypothetical protein |
| FN | CDS | NC_002516.2 | 3195900..3196691 | + | enoyl-CoA hydratase |
| FN | CDS | NC_002516.2 | 3206915..3207166 | + | outer membrane lipoprotein OprI |
| FN | CDS | NC_002516.2 | 3290694..3291173 | - | hypothetical protein |
| FN | CDS | NC_002516.2 | 3414401..3414613 | + | ribosome modulation factor |
| FN | CDS | NC_002516.2 | 3527670..3528212 | - | hypothetical protein |
| FN | CDS | NC_002516.2 | 3535971..3536579 | - | imidazole glycerol phosphate synthase subunit HisH |
| FN | CDS | NC_002516.2 | 3684163..3684717 | - | hypothetical protein |
| FN | CDS | NC_002516.2 | 3684905..3685762 | - | hypothetical protein |
| FN | CDS | NC_002516.2 | 4008938..4009051 | - | hypothetical protein |
| FN | CDS | NC_002516.2 | 4062426..4062899 | - | 2-C-methyl-D-erythritol 2,4-cyclodiphosphate synthase |
| FN | CDS | NC_002516.2 | 4165719..4165880 | - | MexR antirepressor ArmR |
| FN | CDS | NC_002516.2 | 4327395..4327559 | + | hypothetical protein |
| FN | CDS | NC_002516.2 | 4783771..4784139 | - | preprotein translocase subunit SecE |
| FN | CDS | NC_002516.2 | 4814213..4814899 | + | hypothetical protein |
| FN | CDS | NC_002516.2 | 5069082..5069531 | - | type 4 fimbrial protein PilA |
| FN | CDS | NC_002516.2 | 5206720..5207007 | - | hypothetical protein |
| FN | CDS | NC_002516.2 | 5332354..5332743 | - | preprotein translocase subunit SecG |
| FN | rRNA | NC_002516.2 | 1668840..1669799 | + | 4.5S ribosomal RNA |
| FP | CDS | NC_002516.2 | 56169..56447 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 65541..65678 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 68883..69053 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 69056..69190 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 69272..69526 | + | type II toxin-antitoxin system death-on-curing family toxin |
| FP | CDS | NC_002516.2 | 153473..154357 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 352227..352316 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 375511..375978 | + | DUF2269 domain-containing protein |
| FP | CDS | NC_002516.2 | 413364..413657 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 433933..434103 | - | DUF167 family protein |
| FP | CDS | NC_002516.2 | 527176..527361 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 588695..588931 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 628213..628338 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 674046..674216 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 699125..699265 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 740138..740461 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 788927..789142 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 789144..789356 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 794268..794369 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 805596..805784 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 816549..817025 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 828855..829307 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 944467..944640 | - | XapX domain-containing protein |
| FP | CDS | NC_002516.2 | 1029541..1029732 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 1064340..1064558 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 1072125..1072766 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 1099361..1099672 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 1116213..1116500 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 1215284..1215727 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 1215856..1215987 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 1340734..1340889 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 1376094..1376330 | - | 2Fe-2S iron-sulfur cluster-binding protein |
| FP | CDS | NC_002516.2 | 1444777..1445508 | + | class I SAM-dependent methyltransferase |
| FP | CDS | NC_002516.2 | 1463065..1463403 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 1479788..1480048 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 1498346..1499005 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 1697188..1697919 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 1697865..1698344 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 1775945..1776034 | + | K(+)-transporting ATPase subunit F |
| FP | CDS | NC_002516.2 | 1835118..1835309 | + | PLDc N-terminal domain-containing protein |
| FP | CDS | NC_002516.2 | 1844144..1844542 | - | type III secretion system needle length determinant |
| FP | CDS | NC_002516.2 | 1878447..1878548 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 1948389..1948505 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 1997076..1997192 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 2017863..2017970 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 2059713..2060000 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 2060169..2060552 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 2080476..2080664 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 2094139..2094285 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 2117030..2118496 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 2119751..2120179 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 2165064..2165204 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 2207951..2208169 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 2238860..2239420 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 2281755..2281964 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 2351431..2351880 | - | DUF421 domain-containing protein |
| FP | CDS | NC_002516.2 | 2355684..2355806 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 2361293..2361457 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 2364180..2364314 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 2400441..2400671 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 2406961..2407131 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 2412082..2412348 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 2472077..2472370 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 2473043..2473171 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 2562140..2562283 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 2725768..2725977 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 2752866..2753465 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 2753519..2754445 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 2803504..2803623 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 2845916..2846170 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 2881520..2881711 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 2893100..2893246 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 2901752..2901961 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 2904566..2904712 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 3010833..3011099 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 3016675..3016884 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 3017818..3018261 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 3031561..3031752 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 3079850..3083482 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 3083485..3086145 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 3087813..3088037 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 3089644..3090249 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 3108663..3108752 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 3110522..3110743 | - | ParD-like family protein |
| FP | CDS | NC_002516.2 | 3112422..3112655 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 3113767..3114084 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 3123059..3123268 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 3125429..3125524 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 3129158..3129424 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 3147641..3147733 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 3150886..3152202 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 3218444..3218572 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 3290694..3290825 | - | hypothetical protein |
| FP | CDS | NC_002516.2 | 3290847..3291119 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 3306112..3306294 | + | CbtB-domain containing protein |
| FP | CDS | NC_002516.2 | 3347768..3348076 | + | hypothetical protein |
| FP | CDS | NC_002516.2 | 3550746..3552986 | - | bifunctional prephenate dehydrogenase/3-phosphoshikimate 1-carboxyvinyltransferase |
| FP | CDS | NC_002516.2 | 3588444..3588656 | - | DUF1289 domain-containing protein |
| FP | CDS | NC_002516.2 | 3607065..3608120 | - | glycosyltransferase family 1 protein |
| FP | CDS | NC_002516.2 | 3647387..3647746 | - | hypothetical protein |

> Primeros pares con coordenadas no exactas

| tipo | seqid | ref | Δstart | Δend | solape |
|---|---|---|---|---|---|
| CDS | NC_002516.2 | 10434..12488 | +0 | +3 | 0.9985 |
| CDS | NC_002516.2 | 29896..30498 | +219 | +0 | 0.6368 |
| CDS | NC_002516.2 | 34624..35538 | +15 | +0 | 0.9836 |
| CDS | NC_002516.2 | 40589..40816 | +6 | +0 | 0.9737 |
| CDS | NC_002516.2 | 57212..58594 | +0 | -129 | 0.9067 |
| CDS | NC_002516.2 | 61879..62388 | -156 | +0 | 0.7658 |
| CDS | NC_002516.2 | 73468..73923 | +30 | +0 | 0.9342 |
| CDS | NC_002516.2 | 91521..94826 | +0 | +225 | 0.9363 |
| CDS | NC_002516.2 | 104006..104851 | +60 | +0 | 0.9291 |
| CDS | NC_002516.2 | 110976..112907 | -9 | +0 | 0.9954 |
| CDS | NC_002516.2 | 117931..119130 | +87 | +0 | 0.9275 |
| CDS | NC_002516.2 | 127378..128502 | +24 | +0 | 0.9787 |
| CDS | NC_002516.2 | 138818..140167 | +15 | +0 | 0.9889 |
| CDS | NC_002516.2 | 151936..153306 | +0 | -72 | 0.9475 |
| CDS | NC_002516.2 | 160552..161448 | +0 | +177 | 0.8352 |
| CDS | NC_002516.2 | 169361..169906 | -87 | +0 | 0.8626 |
| CDS | NC_002516.2 | 173811..174638 | -90 | +0 | 0.9020 |
| CDS | NC_002516.2 | 177307..178458 | +48 | +0 | 0.9583 |
| CDS | NC_002516.2 | 192965..193672 | +0 | +120 | 0.8551 |
| CDS | NC_002516.2 | 283818..284105 | +0 | -33 | 0.8854 |
| CDS | NC_002516.2 | 284197..284673 | +0 | -12 | 0.9748 |
| CDS | NC_002516.2 | 286338..287024 | +0 | +15 | 0.9786 |
| CDS | NC_002516.2 | 342292..343173 | -30 | +0 | 0.9671 |
| CDS | NC_002516.2 | 348440..349051 | +0 | -36 | 0.9412 |
| CDS | NC_002516.2 | 349050..350069 | +78 | +0 | 0.9235 |
| CDS | NC_002516.2 | 350890..351588 | +69 | +0 | 0.9013 |
| CDS | NC_002516.2 | 358931..359920 | +57 | +0 | 0.9424 |
| CDS | NC_002516.2 | 367457..368422 | +36 | +0 | 0.9627 |
| CDS | NC_002516.2 | 368462..370405 | +0 | +33 | 0.9833 |
| CDS | NC_002516.2 | 380903..382033 | +93 | +0 | 0.9178 |
| CDS | NC_002516.2 | 400248..401072 | -54 | +0 | 0.9386 |
| CDS | NC_002516.2 | 416009..417406 | +0 | -27 | 0.9807 |
| CDS | NC_002516.2 | 424670..425341 | -63 | +0 | 0.9143 |
| CDS | NC_002516.2 | 447391..448350 | +0 | -66 | 0.9313 |
| CDS | NC_002516.2 | 464568..465983 | +0 | -42 | 0.9703 |
| CDS | NC_002516.2 | 482111..482665 | +0 | -15 | 0.9730 |
| CDS | NC_002516.2 | 484404..484838 | +18 | +0 | 0.9586 |
| CDS | NC_002516.2 | 492080..493357 | +0 | -3 | 0.9977 |
| CDS | NC_002516.2 | 505047..505586 | +135 | +0 | 0.7500 |
| CDS | NC_002516.2 | 508959..509753 | +3 | +0 | 0.9962 |
| CDS | NC_002516.2 | 518083..520635 | -246 | +0 | 0.9121 |
| CDS | NC_002516.2 | 535085..535489 | +3 | +0 | 0.9926 |
| CDS | NC_002516.2 | 538217..539143 | +0 | +39 | 0.9596 |
| CDS | NC_002516.2 | 550867..551793 | +0 | -9 | 0.9903 |
| CDS | NC_002516.2 | 556275..557276 | +0 | +81 | 0.9252 |
| CDS | NC_002516.2 | 558377..559078 | +0 | -33 | 0.9530 |
| CDS | NC_002516.2 | 581133..581660 | -45 | +0 | 0.9215 |
| CDS | NC_002516.2 | 590896..591402 | +0 | -243 | 0.5207 |
| CDS | NC_002516.2 | 596970..597479 | +0 | +72 | 0.8763 |
| CDS | NC_002516.2 | 614952..615608 | +0 | -231 | 0.6484 |
| CDS | NC_002516.2 | 624994..627765 | +0 | -45 | 0.9838 |
| CDS | NC_002516.2 | 641073..641426 | -63 | +0 | 0.8489 |
| CDS | NC_002516.2 | 675390..675839 | +96 | +0 | 0.7867 |
| CDS | NC_002516.2 | 716731..717153 | +9 | +0 | 0.9787 |
| CDS | NC_002516.2 | 718729..720072 | +0 | +96 | 0.9333 |
| CDS | NC_002516.2 | 727608..730703 | +0 | -21 | 0.9932 |
| CDS | NC_002516.2 | 740203..741348 | +261 | +0 | 0.7723 |
| CDS | NC_002516.2 | 770156..770818 | -15 | +0 | 0.9779 |
| CDS | NC_002516.2 | 784698..785174 | +54 | +0 | 0.8868 |
| CDS | NC_002516.2 | 785969..786925 | +105 | +0 | 0.8903 |
| CDS | NC_002516.2 | 791370..792632 | +186 | +0 | 0.8527 |
| CDS | NC_002516.2 | 801967..802239 | +0 | -60 | 0.7802 |
| CDS | NC_002516.2 | 804381..804560 | +6 | +0 | 0.9667 |
| CDS | NC_002516.2 | 842878..843618 | +0 | +18 | 0.9763 |
| CDS | NC_002516.2 | 856942..857871 | +0 | -12 | 0.9871 |
| CDS | NC_002516.2 | 860175..863300 | +0 | -132 | 0.9578 |
| CDS | NC_002516.2 | 863455..863664 | +6 | +0 | 0.9714 |
| CDS | NC_002516.2 | 875953..876603 | +39 | +0 | 0.9401 |
| CDS | NC_002516.2 | 876617..878608 | +0 | +330 | 0.8579 |
| CDS | NC_002516.2 | 879027..879560 | +237 | +0 | 0.5562 |
| CDS | NC_002516.2 | 899138..899656 | +0 | -153 | 0.7052 |
| CDS | NC_002516.2 | 925007..927214 | +0 | -222 | 0.8995 |
| CDS | NC_002516.2 | 935427..935975 | +45 | +0 | 0.9180 |
| CDS | NC_002516.2 | 937437..938039 | -36 | +0 | 0.9437 |
| CDS | NC_002516.2 | 977420..977743 | +0 | -63 | 0.8056 |
| CDS | NC_002516.2 | 981422..982888 | +3 | +0 | 0.9980 |
| CDS | NC_002516.2 | 1016657..1017994 | -3 | +0 | 0.9978 |
| CDS | NC_002516.2 | 1039044..1039229 | -9 | +0 | 0.9538 |
| CDS | NC_002516.2 | 1045832..1046461 | +201 | +0 | 0.6810 |
| CDS | NC_002516.2 | 1058812..1059606 | +147 | +0 | 0.8151 |
| CDS | NC_002516.2 | 1061207..1062034 | +0 | -90 | 0.8913 |
| CDS | NC_002516.2 | 1071239..1071877 | +0 | +66 | 0.9064 |
| CDS | NC_002516.2 | 1083972..1085564 | +21 | +0 | 0.9868 |
| CDS | NC_002516.2 | 1097296..1098288 | +93 | +0 | 0.9063 |
| CDS | NC_002516.2 | 1105888..1107000 | +48 | +0 | 0.9569 |
| CDS | NC_002516.2 | 1113129..1114718 | +96 | +0 | 0.9396 |
| CDS | NC_002516.2 | 1118158..1119519 | -114 | +0 | 0.9228 |
| CDS | NC_002516.2 | 1131674..1133947 | +24 | +0 | 0.9894 |
| CDS | NC_002516.2 | 1134024..1135202 | +33 | +0 | 0.9720 |
| CDS | NC_002516.2 | 1151435..1151908 | +75 | +0 | 0.8418 |
| CDS | NC_002516.2 | 1151919..1152575 | +0 | -57 | 0.9132 |
| CDS | NC_002516.2 | 1153637..1155547 | +69 | +0 | 0.9639 |
| CDS | NC_002516.2 | 1211888..1212571 | +0 | -129 | 0.8114 |
| CDS | NC_002516.2 | 1224641..1225330 | +18 | +0 | 0.9739 |
| CDS | NC_002516.2 | 1225757..1226278 | +0 | -108 | 0.7931 |
| CDS | NC_002516.2 | 1227474..1228205 | -57 | +0 | 0.9278 |
| CDS | NC_002516.2 | 1228212..1229198 | +0 | +30 | 0.9705 |
| CDS | NC_002516.2 | 1256094..1257452 | -258 | +0 | 0.8404 |
| CDS | NC_002516.2 | 1299558..1300175 | +6 | +0 | 0.9903 |
| CDS | NC_002516.2 | 1315916..1317202 | +0 | +192 | 0.8702 |

- Aviso: 1 feature(s) GFF multi-segmento fusionados por ID; se comparan por su intervalo envolvente.

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
/opt/conda/bin/python benchmarks/external/external_benchmark.py --genome data/reference_genomes/paeruginosa/paeruginosa_genomic.fna.gz --reference data/reference_genomes/paeruginosa/paeruginosa_genomic.gff.gz --tools bakta,prokka,dfast --run --threads 1 --json reports/external_benchmark/benchmark.json --md reports/external_benchmark/benchmark.md --bakta-db /data/bakta_db
```

Para comparar salidas producidas en otro host, usar `--query-gff TOOL=RUTA` sin `--run`; el origen queda como `precomputed`.
