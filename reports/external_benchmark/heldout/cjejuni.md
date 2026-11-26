# Benchmark externo CPU: Bakta, Prokka y DFAST

**Generado:** 2026-11-26T02:22:31Z

**Esquema:** `external-benchmark-run/v1`

**Alcance:** Comparación local de CDS/tRNA/rRNA de Bakta, Prokka y DFAST contra una referencia GFF3/GBFF. CPU, sin descargas de datos.

**Comando reproducible:**

```bash
/opt/conda/bin/python /opt/benchmark/external/external_benchmark.py --genome /work/data/reference_genomes/cjejuni/cjejuni_genomic.fna.gz --reference /work/data/reference_genomes/cjejuni/cjejuni_genomic.gff.gz --tools all --run --threads 4 --timeout 5400 --bakta-db /data/bakta_db --extra-arg dfast=--dbroot --extra-arg dfast=/data/dfast_db --save-gffs /work/reports/external_benchmark/heldout/cjejuni/evidence --json /work/reports/external_benchmark/heldout/cjejuni.json --md /work/reports/external_benchmark/heldout/cjejuni.md
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

- FASTA: `data/reference_genomes/cjejuni/cjejuni_genomic.fna.gz` (SHA256 `1a0c7e733310de29…`, 1 contig, 1641481 bp).
- Referencia: `data/reference_genomes/cjejuni/cjejuni_genomic.gff.gz` (formato `gff3`, SHA256 `592fae15135b0ca6…`).
- Conteos de referencia: CDS **1572**, tRNA **43**, rRNA **9**.
- CDS pseudo excluidos en referencia: **0**.
- SHA256 FASTA: `1a0c7e733310de29a6c715cf5e6825727f9ba3c9718b504cd5f5af03591599d9`.
- SHA256 referencia: `592fae15135b0ca61d29b40c06f21fcc0f16bb1e37dbb85da78ecb6760b8a7e1`.

## Validación

- Inicio: `2026-11-26T02:03:11Z`; fin: `2026-11-26T02:22:31Z`; total: **1160.43 s**.
- JSON declarado: `reports/external_benchmark/heldout/cjejuni.json`.
- Markdown declarado: `reports/external_benchmark/heldout/cjejuni.md`.

| check | estado | detalle |
|---|---|---|
| genome_path | ok | FASTA existe y es legible |
| genome_sha256 | ok | SHA256 1a0c7e733310de29a6c715cf5e6825727f9ba3c9718b504cd5f5af03591599d9 |
| reference_path | ok | referencia gff3 con 1624 features |
| reference_sha256 | ok | SHA256 592fae15135b0ca61d29b40c06f21fcc0f16bb1e37dbb85da78ecb6760b8a7e1 |
| reference_features | ok | CDS=1572 tRNA=43 rRNA=9 pseudo_excluidos=0 |
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
| Bakta | executed | ok | 0.9638 | 0.9885 | 1.0000 | 0.9404 | 0.9901 | 0.9646 | 545.8 |
| Prokka | executed | ok | 0.9641 | 0.9885 | 1.0000 | 0.9404 | 0.9908 | 0.9649 | 375.6 |
| DFAST | executed | ok | 0.9641 | 0.9885 | 1.0000 | 0.9409 | 0.9901 | 0.9649 | 219.9 |

La precisión/recall/F1 se calculan contra la referencia local. Si una herramienta no está instalada, su estado es `not_available` y no hay métricas que mostrar.

## Bakta

- Estado: **ok**.
- Fuente: `executed`.
- Bakta DB: `ok` (válida: **sí**, ruta `/data/bakta_db`, versión `0`).
- Versión: `bakta 1.12.2` (estado `ok`).
- Ejecución: **ok**, 545.76 s, RSS máx. 2652.0 MiB (`2026-11-26T02:03:17Z` → `2026-11-26T02:12:23Z`).
- Salida: `/tmp/external_benchmark_bakta_kxs2xxf1/out/benchmark.gff3` (2033652 bytes, SHA256 `a6305aee46e35808…`).
- Comando: `/usr/local/bin/bakta --output /tmp/external_benchmark_bakta_kxs2xxf1/out --prefix benchmark --threads 4 --force --skip-plot --db /data/bakta_db /tmp/benchmark_genome_0c5cy8jf.fna`.
- Anotación comparada: `reports/external_benchmark/heldout/cjejuni/evidence/bakta/benchmark.gff3` (formato `gff3`, SHA256 `a6305aee46e35808…`).

### Features (Bakta)

| tipo | ref | query | TP | FP | FN | precisión | recall | F1 | exactos | borde distinto | solape medio |
|---|---|---|---|---|---|---|---|---|---|---|---|
| CDS | 1572 | 1657 | 1556 | 101 | 16 | 0.9390 | 0.9898 | 0.9638 | 1445 | 111 | 0.9965 |
| tRNA | 43 | 44 | 43 | 1 | 0 | 0.9773 | 1.0000 | 0.9885 | 43 | 0 | 1.0000 |
| rRNA | 9 | 9 | 9 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 0 | 9 | 0.9519 |

### Atributos (Bakta)

| tipo | pares | product en ambos | product literal | product tokens | gene literal | locus literal |
|---|---|---|---|---|---|---|
| CDS | 1556 | 1556 | 0.3425 | 0.5264 | 0.8589 | 0.0000 |
| tRNA | 43 | 43 | 0.0000 | 0.6473 | 0.0000 | 0.0000 |
| rRNA | 9 | 9 | 1.0000 | 1.0000 | — | 0.0000 |

### Hipotéticas (Bakta)

| tipo | ref CDS | ref con product | ref hyp | % ref | query CDS | query con product | query hyp | % query |
|---|---|---|---|---|---|---|---|---|
| CDS | 1572 | 1572 | 192 | 12.2 | 1657 | 1657 | 223 | 13.5 |
| tRNA | 0 | 0 | 0 | — | 0 | 0 | 0 | — |
| rRNA | 0 | 0 | 0 | — | 0 | 0 | 0 | — |

#### Diferencias

- Falsos negativos (referencia sin pareja): **16**.
- Falsos positivos (herramienta sin pareja): **102**.
- Pares solo con hebra distinta: **0**.
- Pares emparejados con algún borde distinto: **120**.

Las listas de diferencias están truncadas a 100 elementos por lista; los conteos son completos.

> Primeras diferencias no emparejadas

| tipo | feature | seqid | coordenadas | hebra | product/locus |
|---|---|---|---|---|---|
| FN | CDS | NC_002163.1 | 167050..167794 | + | hypothetical protein |
| FN | CDS | NC_002163.1 | 315108..315233 | + | hypothetical protein |
| FN | CDS | NC_002163.1 | 383559..383657 | + | hypothetical protein |
| FN | CDS | NC_002163.1 | 460364..460513 | + | exporting protein |
| FN | CDS | NC_002163.1 | 627118..627306 | + | periplasmic protein |
| FN | CDS | NC_002163.1 | 702810..702917 | + | hypothetical protein |
| FN | CDS | NC_002163.1 | 703027..703119 | + | hypothetical protein |
| FN | CDS | NC_002163.1 | 767808..768032 | + | hypothetical protein |
| FN | CDS | NC_002163.1 | 815432..815524 | - | hypothetical protein |
| FN | CDS | NC_002163.1 | 908693..908833 | + | hypothetical protein |
| FN | CDS | NC_002163.1 | 920967..921203 | - | hypothetical protein |
| FN | CDS | NC_002163.1 | 1091407..1091529 | - | hypothetical protein |
| FN | CDS | NC_002163.1 | 1133172..1133363 | - | integral membrane protein |
| FN | CDS | NC_002163.1 | 1136955..1137581 | + | 5-formyltetrahydrofolate cyclo-ligase |
| FN | CDS | NC_002163.1 | 1522347..1522460 | + | 50S ribosomal protein L36 |
| FN | CDS | NC_002163.1 | 1640904..1641386 | - | Holliday junction resolvase |
| FP | CDS | NC_002163.1 | 49008..50156 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 67734..67958 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 67897..68136 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 68337..68504 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 71155..71562 | - | hypothetical protein |
| FP | CDS | NC_002163.1 | 83292..83642 | - | Cation-binding protein |
| FP | CDS | NC_002163.1 | 83618..83914 | - | hypothetical protein |
| FP | CDS | NC_002163.1 | 167104..167322 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 167432..167794 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 181866..182075 | - | Protein YjdM C-terminal domain-containing protein |
| FP | CDS | NC_002163.1 | 197123..197464 | - | DUF4878 domain-containing protein |
| FP | CDS | NC_002163.1 | 206063..206383 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 206897..207010 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 207073..207189 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 207128..207313 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 207392..207661 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 208012..208215 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 208487..208669 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 209067..209198 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 226807..227934 | - | MCP-domain signal transduction protein |
| FP | CDS | NC_002163.1 | 227910..228464 | - | hypothetical protein |
| FP | CDS | NC_002163.1 | 267638..268597 | - | MFS transporter |
| FP | CDS | NC_002163.1 | 268588..268809 | - | hypothetical protein |
| FP | CDS | NC_002163.1 | 412341..412928 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 412937..413476 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 413504..414433 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 420548..420892 | - | hypothetical protein |
| FP | CDS | NC_002163.1 | 449602..449865 | + | Altronate hydrolase N-terminus |
| FP | CDS | NC_002163.1 | 449868..451034 | + | UxaA family hydrolase |
| FP | CDS | NC_002163.1 | 451046..452278 | + | Major facilitator superfamily (MFS) profile domain-containing protein |
| FP | CDS | NC_002163.1 | 455461..455682 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 455718..456899 | + | Lactaldehyde dehydrogenase |
| FP | CDS | NC_002163.1 | 465495..466148 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 485719..486201 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 486213..486488 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 486642..487343 | + | Na/Pi cotransporter family protein |
| FP | CDS | NC_002163.1 | 527885..528163 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 529974..530255 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 530261..531130 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 577027..577638 | + | DUF2920 domain-containing protein |
| FP | CDS | NC_002163.1 | 577608..578255 | + | DUF2920 domain-containing protein |
| FP | CDS | NC_002163.1 | 587868..588383 | + | Autotransporter |
| FP | CDS | NC_002163.1 | 615011..615277 | - | hypothetical protein |
| FP | CDS | NC_002163.1 | 615765..615965 | - | hypothetical protein |
| FP | CDS | NC_002163.1 | 615977..616192 | - | hypothetical protein |
| FP | CDS | NC_002163.1 | 627316..627576 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 627610..627831 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 627909..628226 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 628193..629062 | + | Potassium-transporting ATPase A chain |
| FP | CDS | NC_002163.1 | 631062..631472 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 631667..633487 | + | sensor histidine kinase KdpD |
| FP | CDS | NC_002163.1 | 691596..692672 | + | Haemagglutination activity domain containing protein |
| FP | CDS | NC_002163.1 | 692878..693051 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 693048..693212 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 693191..694117 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 694230..694532 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 694772..694939 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 695197..695514 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 703566..703811 | + | Transposase |
| FP | CDS | NC_002163.1 | 703986..704369 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 765304..766059 | + | Flagellar biosynthesis protein FlgG |
| FP | CDS | NC_002163.1 | 766075..766419 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 766370..766597 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 812155..812391 | + | Aryl-sulfate sulfotransferase |
| FP | CDS | NC_002163.1 | 812961..813326 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 813269..813448 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 813421..813591 | + | Arylsulfate sulfotransferase |
| FP | CDS | NC_002163.1 | 814351..814563 | - | Cytochrome C |
| FP | CDS | NC_002163.1 | 814539..815099 | - | c-type cytochrome |
| FP | CDS | NC_002163.1 | 815050..815337 | - | hypothetical protein |
| FP | CDS | NC_002163.1 | 890799..891473 | - | Methyl-accepting chemotaxis protein |
| FP | CDS | NC_002163.1 | 891663..892379 | - | histidine kinase |
| FP | CDS | NC_002163.1 | 907066..907308 | + | Hemagglutination domain protein |
| FP | CDS | NC_002163.1 | 919988..920953 | - | Uncharacterized protein Cj0987c |
| FP | CDS | NC_002163.1 | 1001218..1001577 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 1001574..1001837 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 1079724..1080041 | - | Glycosyltransferase family 92 protein |
| FP | CDS | NC_002163.1 | 1228282..1228641 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 1228638..1229075 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 1250852..1251394 | + | Transferase |
| FP | CDS | NC_002163.1 | 1251467..1251787 | + | 3-hydroxyacyl-CoA dehydrogenase NAD binding domain-containing protein |
| FP | CDS | NC_002163.1 | 1251775..1252161 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 1253417..1253686 | + | Class I SAM-dependent methyltransferase |
| FP | CDS | NC_002163.1 | 1263507..1263689 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 1325668..1326147 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 1326150..1326605 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 1326571..1327107 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 1327120..1327257 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 1327236..1328282 | + | Cystathionine beta-lyase |
| FP | CDS | NC_002163.1 | 1329668..1330444 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 1330451..1331008 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 1403970..1404182 | - | hypothetical protein |
| FP | CDS | NC_002163.1 | 1404311..1405150 | - | hypothetical protein |
| FP | CDS | NC_002163.1 | 1448603..1449145 | - | 4Fe-4S Mo/W bis-MGD-type domain-containing protein |
| FP | CDS | NC_002163.1 | 1460423..1460707 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 1460737..1461483 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 1490206..1490379 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 1490511..1491152 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 1491162..1491341 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 1597008..1597658 | + | Autotransporter |

> Primeros pares con coordenadas no exactas

| tipo | seqid | ref | Δstart | Δend | solape |
|---|---|---|---|---|---|
| CDS | NC_002163.1 | 6703..8010 | +6 | +0 | 0.9954 |
| CDS | NC_002163.1 | 27402..28304 | +0 | -57 | 0.9369 |
| CDS | NC_002163.1 | 46424..50156 | +0 | -1129 | 0.6976 |
| CDS | NC_002163.1 | 59493..61289 | +3 | +0 | 0.9983 |
| CDS | NC_002163.1 | 70274..71098 | +0 | -36 | 0.9564 |
| CDS | NC_002163.1 | 121680..122369 | +3 | +0 | 0.9957 |
| CDS | NC_002163.1 | 127776..128237 | +0 | -12 | 0.9740 |
| CDS | NC_002163.1 | 136612..136869 | +60 | +0 | 0.7674 |
| CDS | NC_002163.1 | 148819..150600 | +12 | +0 | 0.9933 |
| CDS | NC_002163.1 | 169054..169962 | +0 | -3 | 0.9967 |
| CDS | NC_002163.1 | 172563..172694 | +0 | -6 | 0.9545 |
| CDS | NC_002163.1 | 177995..179200 | +6 | +0 | 0.9950 |
| CDS | NC_002163.1 | 196534..197139 | +0 | +33 | 0.9484 |
| CDS | NC_002163.1 | 199236..201233 | +3 | +0 | 0.9985 |
| CDS | NC_002163.1 | 210767..211612 | +6 | +0 | 0.9929 |
| CDS | NC_002163.1 | 211616..212803 | +6 | +0 | 0.9949 |
| CDS | NC_002163.1 | 214057..215430 | +0 | -3 | 0.9978 |
| CDS | NC_002163.1 | 217475..218038 | +0 | -6 | 0.9894 |
| CDS | NC_002163.1 | 229831..230310 | +3 | +0 | 0.9938 |
| CDS | NC_002163.1 | 232544..233704 | +24 | +0 | 0.9793 |
| CDS | NC_002163.1 | 251784..253037 | +30 | +0 | 0.9761 |
| CDS | NC_002163.1 | 269830..270489 | +3 | +0 | 0.9955 |
| CDS | NC_002163.1 | 339071..339283 | +15 | +0 | 0.9296 |
| CDS | NC_002163.1 | 390162..390497 | +3 | +0 | 0.9911 |
| CDS | NC_002163.1 | 404585..404995 | +3 | +0 | 0.9927 |
| CDS | NC_002163.1 | 409015..409683 | +0 | -3 | 0.9955 |
| CDS | NC_002163.1 | 419857..420426 | +0 | -51 | 0.9105 |
| CDS | NC_002163.1 | 421850..422479 | +0 | +21 | 0.9677 |
| CDS | NC_002163.1 | 423761..424027 | +0 | -24 | 0.9101 |
| CDS | NC_002163.1 | 438992..443128 | +9 | +0 | 0.9978 |
| CDS | NC_002163.1 | 477755..479089 | +15 | +0 | 0.9888 |
| CDS | NC_002163.1 | 490474..491475 | +0 | -222 | 0.7784 |
| CDS | NC_002163.1 | 505212..506921 | +3 | +0 | 0.9982 |
| CDS | NC_002163.1 | 517060..518352 | +6 | +0 | 0.9954 |
| CDS | NC_002163.1 | 526556..527158 | +42 | +0 | 0.9303 |
| CDS | NC_002163.1 | 544189..545331 | +84 | +0 | 0.9265 |
| CDS | NC_002163.1 | 549825..550535 | +3 | +0 | 0.9958 |
| CDS | NC_002163.1 | 550525..550746 | +0 | -6 | 0.9730 |
| CDS | NC_002163.1 | 587868..591303 | +499 | +0 | 0.8548 |
| CDS | NC_002163.1 | 597828..598667 | -39 | +0 | 0.9556 |
| CDS | NC_002163.1 | 622768..623013 | -39 | +0 | 0.8632 |
| CDS | NC_002163.1 | 623410..624138 | -24 | +0 | 0.9681 |
| CDS | NC_002163.1 | 641890..643395 | +3 | +0 | 0.9980 |
| CDS | NC_002163.1 | 649946..650884 | +0 | -6 | 0.9936 |
| CDS | NC_002163.1 | 652530..653918 | +6 | +0 | 0.9957 |
| CDS | NC_002163.1 | 676227..676976 | +0 | -3 | 0.9960 |
| CDS | NC_002163.1 | 683848..684657 | +3 | +0 | 0.9963 |
| CDS | NC_002163.1 | 689656..691584 | +24 | +0 | 0.9876 |
| CDS | NC_002163.1 | 707840..708634 | +3 | +0 | 0.9962 |
| CDS | NC_002163.1 | 714777..716612 | +0 | -3 | 0.9984 |
| CDS | NC_002163.1 | 735708..736490 | +3 | +0 | 0.9962 |
| CDS | NC_002163.1 | 736487..737011 | +3 | +0 | 0.9943 |
| CDS | NC_002163.1 | 738526..738774 | +3 | +0 | 0.9880 |
| CDS | NC_002163.1 | 741167..742441 | +0 | -6 | 0.9953 |
| CDS | NC_002163.1 | 774427..775152 | +3 | +0 | 0.9959 |
| CDS | NC_002163.1 | 797653..798840 | +0 | -27 | 0.9773 |
| CDS | NC_002163.1 | 811319..812119 | +15 | +0 | 0.9813 |
| CDS | NC_002163.1 | 829486..830148 | +0 | +9 | 0.9866 |
| CDS | NC_002163.1 | 840501..840680 | +0 | -33 | 0.8167 |
| CDS | NC_002163.1 | 900703..901044 | +0 | +30 | 0.9194 |
| CDS | NC_002163.1 | 904409..906697 | +297 | +0 | 0.8702 |
| CDS | NC_002163.1 | 910711..911601 | +15 | +0 | 0.9832 |
| CDS | NC_002163.1 | 921397..922158 | +0 | -3 | 0.9961 |
| CDS | NC_002163.1 | 956627..957076 | +0 | -132 | 0.7067 |
| CDS | NC_002163.1 | 987019..989229 | +0 | -3 | 0.9986 |
| CDS | NC_002163.1 | 1011696..1012148 | +66 | +0 | 0.8543 |
| CDS | NC_002163.1 | 1018565..1018960 | +0 | -108 | 0.7273 |
| CDS | NC_002163.1 | 1019786..1020970 | +0 | -12 | 0.9899 |
| CDS | NC_002163.1 | 1027809..1028081 | +0 | -42 | 0.8462 |
| CDS | NC_002163.1 | 1044593..1045219 | +0 | +12 | 0.9812 |
| CDS | NC_002163.1 | 1055060..1055713 | +0 | -45 | 0.9312 |
| CDS | NC_002163.1 | 1062131..1063210 | +0 | -3 | 0.9972 |
| CDS | NC_002163.1 | 1079195..1080041 | +0 | -319 | 0.6234 |
| CDS | NC_002163.1 | 1085907..1088264 | +0 | -6 | 0.9975 |
| CDS | NC_002163.1 | 1104012..1104635 | +0 | -6 | 0.9904 |
| CDS | NC_002163.1 | 1125011..1125907 | +0 | -3 | 0.9967 |
| CDS | NC_002163.1 | 1129228..1130016 | +18 | +0 | 0.9772 |
| CDS | NC_002163.1 | 1194260..1194796 | +0 | -33 | 0.9385 |
| CDS | NC_002163.1 | 1211677..1212912 | +0 | -3 | 0.9976 |
| CDS | NC_002163.1 | 1248624..1249595 | -21 | +0 | 0.9789 |
| CDS | NC_002163.1 | 1253417..1254092 | +316 | +0 | 0.5325 |
| CDS | NC_002163.1 | 1259677..1261638 | +9 | +0 | 0.9954 |
| CDS | NC_002163.1 | 1261635..1263494 | -21 | +0 | 0.9888 |
| CDS | NC_002163.1 | 1263507..1265454 | +166 | +0 | 0.9148 |
| CDS | NC_002163.1 | 1265451..1267337 | -21 | +0 | 0.9890 |
| CDS | NC_002163.1 | 1290427..1290942 | +0 | -24 | 0.9535 |
| CDS | NC_002163.1 | 1293659..1294111 | +0 | -42 | 0.9073 |
| CDS | NC_002163.1 | 1339759..1341147 | +0 | -18 | 0.9870 |
| CDS | NC_002163.1 | 1370667..1372004 | +0 | -3 | 0.9978 |
| CDS | NC_002163.1 | 1372695..1373867 | +0 | +57 | 0.9537 |
| CDS | NC_002163.1 | 1375007..1377337 | +0 | -3 | 0.9987 |
| CDS | NC_002163.1 | 1397726..1398409 | -3 | +0 | 0.9956 |
| CDS | NC_002163.1 | 1423309..1423572 | +0 | +3 | 0.9888 |
| CDS | NC_002163.1 | 1434062..1435549 | +0 | -18 | 0.9879 |
| CDS | NC_002163.1 | 1446341..1449145 | +0 | -564 | 0.7989 |
| CDS | NC_002163.1 | 1449311..1450024 | +0 | -12 | 0.9832 |
| CDS | NC_002163.1 | 1471517..1472236 | +0 | -33 | 0.9542 |
| CDS | NC_002163.1 | 1477635..1478084 | -30 | +0 | 0.9375 |
| CDS | NC_002163.1 | 1485626..1487092 | +0 | -15 | 0.9898 |
| CDS | NC_002163.1 | 1492003..1493991 | +33 | +0 | 0.9834 |

- Aviso: 7 feature(s) GFF multi-segmento fusionados por ID; se comparan por su intervalo envolvente.

## Prokka

- Estado: **ok**.
- Fuente: `executed`.
- Versión: `prokka 1.15.6` (estado `ok`).
- Ejecución: **ok**, 375.60 s, RSS máx. 3047.9 MiB (`2026-11-26T02:12:27Z` → `2026-11-26T02:18:43Z`).
- Salida: `/tmp/external_benchmark_prokka_3_stn7jf/out/benchmark.gff` (2059506 bytes, SHA256 `855a89df50573183…`).
- Comando: `/usr/local/bin/prokka --outdir /tmp/external_benchmark_prokka_3_stn7jf/out --prefix benchmark --cpus 4 --force --quiet /tmp/benchmark_genome_0c5cy8jf.fna`.
- Anotación comparada: `reports/external_benchmark/heldout/cjejuni/evidence/prokka/benchmark.gff` (formato `gff3`, SHA256 `855a89df50573183…`).

### Features (Prokka)

| tipo | ref | query | TP | FP | FN | precisión | recall | F1 | exactos | borde distinto | solape medio |
|---|---|---|---|---|---|---|---|---|---|---|---|
| CDS | 1572 | 1658 | 1557 | 101 | 15 | 0.9391 | 0.9905 | 0.9641 | 1445 | 112 | 0.9964 |
| tRNA | 43 | 44 | 43 | 1 | 0 | 0.9773 | 1.0000 | 0.9885 | 42 | 1 | 0.9991 |
| rRNA | 9 | 9 | 9 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 0 | 9 | 0.9519 |

### Atributos (Prokka)

| tipo | pares | product en ambos | product literal | product tokens | gene literal | locus literal |
|---|---|---|---|---|---|---|
| CDS | 1557 | 1557 | 0.2897 | 0.4720 | 0.7320 | 0.0000 |
| tRNA | 43 | 43 | 0.0000 | 0.6667 | — | 0.0000 |
| rRNA | 9 | 9 | 1.0000 | 1.0000 | — | 0.0000 |

### Hipotéticas (Prokka)

| tipo | ref CDS | ref con product | ref hyp | % ref | query CDS | query con product | query hyp | % query |
|---|---|---|---|---|---|---|---|---|
| CDS | 1572 | 1572 | 192 | 12.2 | 1658 | 1658 | 607 | 36.6 |
| tRNA | 0 | 0 | 0 | — | 0 | 0 | 0 | — |
| rRNA | 0 | 0 | 0 | — | 0 | 0 | 0 | — |

#### Diferencias

- Falsos negativos (referencia sin pareja): **15**.
- Falsos positivos (herramienta sin pareja): **102**.
- Pares solo con hebra distinta: **0**.
- Pares emparejados con algún borde distinto: **122**.

Las listas de diferencias están truncadas a 100 elementos por lista; los conteos son completos.

> Primeras diferencias no emparejadas

| tipo | feature | seqid | coordenadas | hebra | product/locus |
|---|---|---|---|---|---|
| FN | CDS | NC_002163.1 | 167050..167794 | + | hypothetical protein |
| FN | CDS | NC_002163.1 | 315108..315233 | + | hypothetical protein |
| FN | CDS | NC_002163.1 | 383559..383657 | + | hypothetical protein |
| FN | CDS | NC_002163.1 | 460364..460513 | + | exporting protein |
| FN | CDS | NC_002163.1 | 627118..627306 | + | periplasmic protein |
| FN | CDS | NC_002163.1 | 702810..702917 | + | hypothetical protein |
| FN | CDS | NC_002163.1 | 703027..703119 | + | hypothetical protein |
| FN | CDS | NC_002163.1 | 767808..768032 | + | hypothetical protein |
| FN | CDS | NC_002163.1 | 815432..815524 | - | hypothetical protein |
| FN | CDS | NC_002163.1 | 908693..908833 | + | hypothetical protein |
| FN | CDS | NC_002163.1 | 920967..921203 | - | hypothetical protein |
| FN | CDS | NC_002163.1 | 1091407..1091529 | - | hypothetical protein |
| FN | CDS | NC_002163.1 | 1133172..1133363 | - | integral membrane protein |
| FN | CDS | NC_002163.1 | 1136955..1137581 | + | 5-formyltetrahydrofolate cyclo-ligase |
| FN | CDS | NC_002163.1 | 1522347..1522460 | + | 50S ribosomal protein L36 |
| FP | CDS | NC_002163.1 | 49008..50156 | + | Type II methyltransferase M.PaeR7I |
| FP | CDS | NC_002163.1 | 67734..67958 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 67897..68136 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 68337..68504 | + | Putative malate transporter YflS |
| FP | CDS | NC_002163.1 | 71155..71562 | - | hypothetical protein |
| FP | CDS | NC_002163.1 | 83292..83642 | - | hypothetical protein |
| FP | CDS | NC_002163.1 | 83618..83914 | - | Bacteriohemerythrin |
| FP | CDS | NC_002163.1 | 167104..167322 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 167432..167794 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 181866..182075 | - | hypothetical protein |
| FP | CDS | NC_002163.1 | 197123..197464 | - | hypothetical protein |
| FP | CDS | NC_002163.1 | 206063..206383 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 206897..207010 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 207073..207189 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 207128..207313 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 207392..207661 | + | Temperature-sensitive hemagglutinin tsh autotransporter |
| FP | CDS | NC_002163.1 | 208012..208215 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 208487..208669 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 209067..209198 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 226807..227934 | - | hypothetical protein |
| FP | CDS | NC_002163.1 | 227910..228464 | - | hypothetical protein |
| FP | CDS | NC_002163.1 | 267638..268597 | - | Phosphoglycerate transporter protein |
| FP | CDS | NC_002163.1 | 268588..268809 | - | Phosphoglycerate transporter protein |
| FP | CDS | NC_002163.1 | 412341..412928 | + | Vitamin B12 transporter BtuB |
| FP | CDS | NC_002163.1 | 412937..413476 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 413504..414433 | + | Ferrienterobactin receptor |
| FP | CDS | NC_002163.1 | 420548..420892 | - | hypothetical protein |
| FP | CDS | NC_002163.1 | 449602..449865 | + | Altronate dehydratase |
| FP | CDS | NC_002163.1 | 449868..451034 | + | Altronate dehydratase |
| FP | CDS | NC_002163.1 | 451046..452278 | + | Putative tartrate transporter |
| FP | CDS | NC_002163.1 | 455461..455682 | + | Lactaldehyde dehydrogenase |
| FP | CDS | NC_002163.1 | 455718..456899 | + | Lactaldehyde dehydrogenase |
| FP | CDS | NC_002163.1 | 465495..466148 | + | Ammonium transporter |
| FP | CDS | NC_002163.1 | 485719..486201 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 486213..486488 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 486642..487343 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 527885..528163 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 529974..530255 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 530261..531130 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 577027..577638 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 577608..578255 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 587868..588383 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 615011..615277 | - | Peptide transporter YePEPT |
| FP | CDS | NC_002163.1 | 615765..615965 | - | hypothetical protein |
| FP | CDS | NC_002163.1 | 615977..616192 | - | Peptide transporter YePEPT |
| FP | CDS | NC_002163.1 | 627316..627576 | + | Potassium-transporting ATPase potassium-binding subunit |
| FP | CDS | NC_002163.1 | 627610..627831 | + | Potassium-transporting ATPase potassium-binding subunit 2 |
| FP | CDS | NC_002163.1 | 627909..628226 | + | Potassium-transporting ATPase potassium-binding subunit 2 |
| FP | CDS | NC_002163.1 | 628193..629062 | + | Potassium-transporting ATPase potassium-binding subunit |
| FP | CDS | NC_002163.1 | 631062..631472 | + | Potassium-transporting ATPase KdpC subunit |
| FP | CDS | NC_002163.1 | 631667..633487 | + | Sensor protein KdpD |
| FP | CDS | NC_002163.1 | 691596..692672 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 692878..693051 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 693048..693212 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 693191..694117 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 694230..694532 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 694772..694939 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 695197..695514 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 703566..703811 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 703986..704369 | + | IS200/IS605 family transposase ISCco1 |
| FP | CDS | NC_002163.1 | 765304..766059 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 766075..766419 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 766370..766597 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 812155..812391 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 812961..813326 | + | Arylsulfate sulfotransferase AssT |
| FP | CDS | NC_002163.1 | 813269..813448 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 813421..813591 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 814351..814563 | - | hypothetical protein |
| FP | CDS | NC_002163.1 | 814539..815099 | - | Thiosulfate dehydrogenase |
| FP | CDS | NC_002163.1 | 815050..815337 | - | hypothetical protein |
| FP | CDS | NC_002163.1 | 890799..891473 | - | Methyl-accepting chemotaxis protein 4 |
| FP | CDS | NC_002163.1 | 891663..892379 | - | hypothetical protein |
| FP | CDS | NC_002163.1 | 907066..907308 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 919988..920953 | - | hypothetical protein |
| FP | CDS | NC_002163.1 | 1001218..1001577 | + | Putative NAD(P)H nitroreductase YfkO |
| FP | CDS | NC_002163.1 | 1001574..1001837 | + | Putative NAD(P)H nitroreductase |
| FP | CDS | NC_002163.1 | 1079724..1080041 | - | hypothetical protein |
| FP | CDS | NC_002163.1 | 1228282..1228641 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 1228638..1229075 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 1250852..1251394 | + | UDP-N-acetylbacillosamine N-acetyltransferase |
| FP | CDS | NC_002163.1 | 1251467..1251787 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 1251775..1252161 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 1253417..1253686 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 1263507..1263689 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 1325668..1326147 | + | Putative cryptic C4-dicarboxylate transporter DcuD |
| FP | CDS | NC_002163.1 | 1326150..1326605 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 1326571..1327107 | + | Putative cryptic C4-dicarboxylate transporter DcuD |
| FP | CDS | NC_002163.1 | 1327120..1327257 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 1327236..1328282 | + | Cystathionine beta-lyase |
| FP | CDS | NC_002163.1 | 1329668..1330444 | + | putative protein YxeQ |
| FP | CDS | NC_002163.1 | 1330451..1331008 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 1403970..1404182 | - | hypothetical protein |
| FP | CDS | NC_002163.1 | 1404311..1405150 | - | hypothetical protein |
| FP | CDS | NC_002163.1 | 1448603..1449145 | - | Formate dehydrogenase 2 subunit alpha (cytochrome c-553) |
| FP | CDS | NC_002163.1 | 1460423..1460707 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 1460737..1461483 | + | Putative cryptic C4-dicarboxylate transporter DcuD |
| FP | CDS | NC_002163.1 | 1490206..1490379 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 1490511..1491152 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 1491162..1491341 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 1597008..1597658 | + | hypothetical protein |

> Primeros pares con coordenadas no exactas

| tipo | seqid | ref | Δstart | Δend | solape |
|---|---|---|---|---|---|
| CDS | NC_002163.1 | 6703..8010 | +6 | +0 | 0.9954 |
| CDS | NC_002163.1 | 27402..28304 | +0 | -57 | 0.9369 |
| CDS | NC_002163.1 | 46424..50156 | +0 | -1129 | 0.6976 |
| CDS | NC_002163.1 | 59493..61289 | +3 | +0 | 0.9983 |
| CDS | NC_002163.1 | 70274..71098 | +0 | -36 | 0.9564 |
| CDS | NC_002163.1 | 121680..122369 | +3 | +0 | 0.9957 |
| CDS | NC_002163.1 | 127776..128237 | +0 | -12 | 0.9740 |
| CDS | NC_002163.1 | 136612..136869 | +60 | +0 | 0.7674 |
| CDS | NC_002163.1 | 148819..150600 | +12 | +0 | 0.9933 |
| CDS | NC_002163.1 | 169054..169962 | +0 | -3 | 0.9967 |
| CDS | NC_002163.1 | 172563..172694 | +0 | -6 | 0.9545 |
| CDS | NC_002163.1 | 177995..179200 | +6 | +0 | 0.9950 |
| CDS | NC_002163.1 | 196534..197139 | +0 | +33 | 0.9484 |
| CDS | NC_002163.1 | 199236..201233 | +3 | +0 | 0.9985 |
| CDS | NC_002163.1 | 210767..211612 | +6 | +0 | 0.9929 |
| CDS | NC_002163.1 | 211616..212803 | +6 | +0 | 0.9949 |
| CDS | NC_002163.1 | 214057..215430 | +0 | -3 | 0.9978 |
| CDS | NC_002163.1 | 217475..218038 | +0 | -6 | 0.9894 |
| CDS | NC_002163.1 | 229831..230310 | +3 | +0 | 0.9938 |
| CDS | NC_002163.1 | 232544..233704 | +24 | +0 | 0.9793 |
| CDS | NC_002163.1 | 251784..253037 | +30 | +0 | 0.9761 |
| CDS | NC_002163.1 | 269830..270489 | +3 | +0 | 0.9955 |
| CDS | NC_002163.1 | 339071..339283 | +15 | +0 | 0.9296 |
| CDS | NC_002163.1 | 390162..390497 | +3 | +0 | 0.9911 |
| CDS | NC_002163.1 | 404585..404995 | +3 | +0 | 0.9927 |
| CDS | NC_002163.1 | 409015..409683 | +0 | -3 | 0.9955 |
| CDS | NC_002163.1 | 419857..420426 | +0 | -51 | 0.9105 |
| CDS | NC_002163.1 | 421850..422479 | +0 | +21 | 0.9677 |
| CDS | NC_002163.1 | 423761..424027 | +0 | -24 | 0.9101 |
| CDS | NC_002163.1 | 438992..443128 | +9 | +0 | 0.9978 |
| CDS | NC_002163.1 | 477755..479089 | +15 | +0 | 0.9888 |
| CDS | NC_002163.1 | 490474..491475 | +0 | -222 | 0.7784 |
| CDS | NC_002163.1 | 505212..506921 | +3 | +0 | 0.9982 |
| CDS | NC_002163.1 | 517060..518352 | +6 | +0 | 0.9954 |
| CDS | NC_002163.1 | 526556..527158 | +42 | +0 | 0.9303 |
| CDS | NC_002163.1 | 544189..545331 | +84 | +0 | 0.9265 |
| CDS | NC_002163.1 | 549825..550535 | +3 | +0 | 0.9958 |
| CDS | NC_002163.1 | 550525..550746 | +0 | -6 | 0.9730 |
| CDS | NC_002163.1 | 587868..591303 | +499 | +0 | 0.8548 |
| CDS | NC_002163.1 | 597828..598667 | -39 | +0 | 0.9556 |
| CDS | NC_002163.1 | 622768..623013 | -39 | +0 | 0.8632 |
| CDS | NC_002163.1 | 623410..624138 | -24 | +0 | 0.9681 |
| CDS | NC_002163.1 | 641890..643395 | +3 | +0 | 0.9980 |
| CDS | NC_002163.1 | 649946..650884 | +0 | -6 | 0.9936 |
| CDS | NC_002163.1 | 652530..653918 | +6 | +0 | 0.9957 |
| CDS | NC_002163.1 | 676227..676976 | +0 | -3 | 0.9960 |
| CDS | NC_002163.1 | 683848..684657 | +3 | +0 | 0.9963 |
| CDS | NC_002163.1 | 689656..691584 | +24 | +0 | 0.9876 |
| CDS | NC_002163.1 | 707840..708634 | +3 | +0 | 0.9962 |
| CDS | NC_002163.1 | 714777..716612 | +0 | -3 | 0.9984 |
| CDS | NC_002163.1 | 735708..736490 | +3 | +0 | 0.9962 |
| CDS | NC_002163.1 | 736487..737011 | +3 | +0 | 0.9943 |
| CDS | NC_002163.1 | 738526..738774 | +3 | +0 | 0.9880 |
| CDS | NC_002163.1 | 741167..742441 | +0 | -6 | 0.9953 |
| CDS | NC_002163.1 | 774427..775152 | +3 | +0 | 0.9959 |
| CDS | NC_002163.1 | 797653..798840 | +0 | -27 | 0.9773 |
| CDS | NC_002163.1 | 811319..812119 | +15 | +0 | 0.9813 |
| CDS | NC_002163.1 | 829486..830148 | +0 | +9 | 0.9866 |
| CDS | NC_002163.1 | 840501..840680 | +0 | -33 | 0.8167 |
| CDS | NC_002163.1 | 900703..901044 | +0 | +30 | 0.9194 |
| CDS | NC_002163.1 | 904409..906697 | +297 | +0 | 0.8702 |
| CDS | NC_002163.1 | 910711..911601 | +15 | +0 | 0.9832 |
| CDS | NC_002163.1 | 921397..922158 | +0 | -3 | 0.9961 |
| CDS | NC_002163.1 | 956627..957076 | +0 | -132 | 0.7067 |
| CDS | NC_002163.1 | 987019..989229 | +0 | -3 | 0.9986 |
| CDS | NC_002163.1 | 1011696..1012148 | +66 | +0 | 0.8543 |
| CDS | NC_002163.1 | 1018565..1018960 | +0 | -108 | 0.7273 |
| CDS | NC_002163.1 | 1019786..1020970 | +0 | -12 | 0.9899 |
| CDS | NC_002163.1 | 1027809..1028081 | +0 | -42 | 0.8462 |
| CDS | NC_002163.1 | 1044593..1045219 | +0 | +12 | 0.9812 |
| CDS | NC_002163.1 | 1055060..1055713 | +0 | -45 | 0.9312 |
| CDS | NC_002163.1 | 1062131..1063210 | +0 | -3 | 0.9972 |
| CDS | NC_002163.1 | 1079195..1080041 | +0 | -319 | 0.6234 |
| CDS | NC_002163.1 | 1085907..1088264 | +0 | -6 | 0.9975 |
| CDS | NC_002163.1 | 1104012..1104635 | +0 | -6 | 0.9904 |
| CDS | NC_002163.1 | 1125011..1125907 | +0 | -3 | 0.9967 |
| CDS | NC_002163.1 | 1129228..1130016 | +18 | +0 | 0.9772 |
| CDS | NC_002163.1 | 1194260..1194796 | +0 | -33 | 0.9385 |
| CDS | NC_002163.1 | 1211677..1212912 | +0 | -3 | 0.9976 |
| CDS | NC_002163.1 | 1248624..1249595 | -21 | +0 | 0.9789 |
| CDS | NC_002163.1 | 1253417..1254092 | +316 | +0 | 0.5325 |
| CDS | NC_002163.1 | 1259677..1261638 | +9 | +0 | 0.9954 |
| CDS | NC_002163.1 | 1261635..1263494 | -21 | +0 | 0.9888 |
| CDS | NC_002163.1 | 1263507..1265454 | +166 | +0 | 0.9148 |
| CDS | NC_002163.1 | 1265451..1267337 | -21 | +0 | 0.9890 |
| CDS | NC_002163.1 | 1290427..1290942 | +0 | -24 | 0.9535 |
| CDS | NC_002163.1 | 1293659..1294111 | +0 | -42 | 0.9073 |
| CDS | NC_002163.1 | 1339759..1341147 | +0 | -18 | 0.9870 |
| CDS | NC_002163.1 | 1370667..1372004 | +0 | -3 | 0.9978 |
| CDS | NC_002163.1 | 1372695..1373867 | +0 | +57 | 0.9537 |
| CDS | NC_002163.1 | 1375007..1377337 | +0 | -3 | 0.9987 |
| CDS | NC_002163.1 | 1397726..1398409 | -3 | +0 | 0.9956 |
| CDS | NC_002163.1 | 1423309..1423572 | +0 | +3 | 0.9888 |
| CDS | NC_002163.1 | 1434062..1435549 | +0 | -18 | 0.9879 |
| CDS | NC_002163.1 | 1446341..1449145 | +0 | -564 | 0.7989 |
| CDS | NC_002163.1 | 1449311..1450024 | +0 | -12 | 0.9832 |
| CDS | NC_002163.1 | 1471517..1472236 | +0 | -33 | 0.9542 |
| CDS | NC_002163.1 | 1477635..1478084 | -30 | +0 | 0.9375 |
| CDS | NC_002163.1 | 1485626..1487092 | +0 | -15 | 0.9898 |
| CDS | NC_002163.1 | 1492003..1493991 | +33 | +0 | 0.9834 |

- Aviso: 7 feature(s) GFF multi-segmento fusionados por ID; se comparan por su intervalo envolvente.

## DFAST

- Estado: **ok**.
- Fuente: `executed`.
- Versión: `DFAST ver. 1.6.2` (estado `ok`).
- Ejecución: **ok**, 219.91 s, RSS máx. 3047.9 MiB (`2026-11-26T02:18:48Z` → `2026-11-26T02:22:28Z`).
- Salida: `/tmp/external_benchmark_dfast_ukzj4b4t/out/genome.gff` (3011966 bytes, SHA256 `e59fcf123a62c5ba…`).
- Comando: `/usr/local/bin/dfast --genome /tmp/benchmark_genome_0c5cy8jf.fna --out /tmp/external_benchmark_dfast_ukzj4b4t/out --cpu 4 --force --dbroot /data/dfast_db`.
- Anotación comparada: `reports/external_benchmark/heldout/cjejuni/evidence/dfast/genome.gff` (formato `gff3`, SHA256 `e59fcf123a62c5ba…`).

### Features (DFAST)

| tipo | ref | query | TP | FP | FN | precisión | recall | F1 | exactos | borde distinto | solape medio |
|---|---|---|---|---|---|---|---|---|---|---|---|
| CDS | 1572 | 1656 | 1556 | 100 | 16 | 0.9396 | 0.9898 | 0.9641 | 1446 | 110 | 0.9966 |
| tRNA | 43 | 44 | 43 | 1 | 0 | 0.9773 | 1.0000 | 0.9885 | 42 | 1 | 0.9991 |
| rRNA | 9 | 9 | 9 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 0 | 9 | 0.9376 |

### Atributos (DFAST)

| tipo | pares | product en ambos | product literal | product tokens | gene literal | locus literal |
|---|---|---|---|---|---|---|
| CDS | 1556 | 1556 | 0.2898 | 0.5161 | 0.9235 | 0.0000 |
| tRNA | 43 | 43 | 1.0000 | 1.0000 | — | 0.0000 |
| rRNA | 9 | 9 | 1.0000 | 1.0000 | — | 0.0000 |

### Hipotéticas (DFAST)

| tipo | ref CDS | ref con product | ref hyp | % ref | query CDS | query con product | query hyp | % query |
|---|---|---|---|---|---|---|---|---|
| CDS | 1572 | 1572 | 192 | 12.2 | 1656 | 1656 | 247 | 14.9 |
| tRNA | 0 | 0 | 0 | — | 0 | 0 | 0 | — |
| rRNA | 0 | 0 | 0 | — | 0 | 0 | 0 | — |

#### Diferencias

- Falsos negativos (referencia sin pareja): **16**.
- Falsos positivos (herramienta sin pareja): **101**.
- Pares solo con hebra distinta: **0**.
- Pares emparejados con algún borde distinto: **120**.

Las listas de diferencias están truncadas a 100 elementos por lista; los conteos son completos.

> Primeras diferencias no emparejadas

| tipo | feature | seqid | coordenadas | hebra | product/locus |
|---|---|---|---|---|---|
| FN | CDS | NC_002163.1 | 167050..167794 | + | hypothetical protein |
| FN | CDS | NC_002163.1 | 315108..315233 | + | hypothetical protein |
| FN | CDS | NC_002163.1 | 383559..383657 | + | hypothetical protein |
| FN | CDS | NC_002163.1 | 460364..460513 | + | exporting protein |
| FN | CDS | NC_002163.1 | 627118..627306 | + | periplasmic protein |
| FN | CDS | NC_002163.1 | 702810..702917 | + | hypothetical protein |
| FN | CDS | NC_002163.1 | 703027..703119 | + | hypothetical protein |
| FN | CDS | NC_002163.1 | 767808..768032 | + | hypothetical protein |
| FN | CDS | NC_002163.1 | 815432..815524 | - | hypothetical protein |
| FN | CDS | NC_002163.1 | 908693..908833 | + | hypothetical protein |
| FN | CDS | NC_002163.1 | 920967..921203 | - | hypothetical protein |
| FN | CDS | NC_002163.1 | 1091407..1091529 | - | hypothetical protein |
| FN | CDS | NC_002163.1 | 1133172..1133363 | - | integral membrane protein |
| FN | CDS | NC_002163.1 | 1136955..1137581 | + | 5-formyltetrahydrofolate cyclo-ligase |
| FN | CDS | NC_002163.1 | 1522347..1522460 | + | 50S ribosomal protein L36 |
| FN | CDS | NC_002163.1 | 1640904..1641386 | - | Holliday junction resolvase |
| FP | CDS | NC_002163.1 | 49008..50156 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 67734..67958 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 67897..68136 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 68337..68504 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 71155..71562 | - | hypothetical protein |
| FP | CDS | NC_002163.1 | 83292..83642 | - | hypothetical protein |
| FP | CDS | NC_002163.1 | 83618..83914 | - | hemerythrin domain-containing protein |
| FP | CDS | NC_002163.1 | 167104..167322 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 167432..167794 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 181866..182075 | - | alkylphosphonate utilization protein |
| FP | CDS | NC_002163.1 | 197123..197464 | - | hypothetical protein |
| FP | CDS | NC_002163.1 | 206063..206383 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 206897..207010 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 207073..207189 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 207128..207313 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 207392..207661 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 208012..208215 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 208487..208669 | + | autotransporter outer membrane beta-barrel domain-containing protein |
| FP | CDS | NC_002163.1 | 209067..209198 | + | autotransporter outer membrane beta-barrel domain-containing protein |
| FP | CDS | NC_002163.1 | 226807..227934 | - | hypothetical protein |
| FP | CDS | NC_002163.1 | 227910..228464 | - | hypothetical protein |
| FP | CDS | NC_002163.1 | 267638..268597 | - | hypothetical protein |
| FP | CDS | NC_002163.1 | 268588..268809 | - | hypothetical protein |
| FP | CDS | NC_002163.1 | 412341..412928 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 412937..413476 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 413504..414433 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 420548..420892 | - | hypothetical protein |
| FP | CDS | NC_002163.1 | 449602..449865 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 449868..451034 | + | UxaA family hydrolase |
| FP | CDS | NC_002163.1 | 451046..452278 | + | MFS transporter |
| FP | CDS | NC_002163.1 | 455461..455682 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 455718..456899 | + | aldehyde dehydrogenase |
| FP | CDS | NC_002163.1 | 465495..466148 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 485719..486201 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 486213..486488 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 486642..487343 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 527885..528163 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 529974..530255 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 530261..531130 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 577027..577638 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 577608..578255 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 587868..588383 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 615011..615277 | - | hypothetical protein |
| FP | CDS | NC_002163.1 | 615765..615965 | - | hypothetical protein |
| FP | CDS | NC_002163.1 | 615977..616192 | - | hypothetical protein |
| FP | CDS | NC_002163.1 | 627316..627576 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 627610..627831 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 627909..628226 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 628193..629062 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 631062..631472 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 631667..633487 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 691596..692672 | + | filamentous hemagglutinin N-terminal domain-containing protein |
| FP | CDS | NC_002163.1 | 692878..693051 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 693048..693212 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 693191..694117 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 694230..694532 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 694772..694939 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 695197..695514 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 703566..703811 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 703986..704369 | + | IS200/IS605 family accessory protein TnpB-related protein |
| FP | CDS | NC_002163.1 | 765304..766059 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 766075..766419 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 766370..766597 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 812155..812391 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 812961..813326 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 813269..813448 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 813421..813591 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 814351..814563 | - | hypothetical protein |
| FP | CDS | NC_002163.1 | 814539..815099 | - | hypothetical protein |
| FP | CDS | NC_002163.1 | 815050..815337 | - | hypothetical protein |
| FP | CDS | NC_002163.1 | 890799..891473 | - | hypothetical protein |
| FP | CDS | NC_002163.1 | 891663..892379 | - | hypothetical protein |
| FP | CDS | NC_002163.1 | 907066..907308 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 919988..920953 | - | hypothetical protein |
| FP | CDS | NC_002163.1 | 1001218..1001577 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 1001574..1001837 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 1079724..1080041 | - | hypothetical protein |
| FP | CDS | NC_002163.1 | 1228282..1228641 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 1228638..1229075 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 1250852..1251394 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 1251467..1251787 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 1251775..1252161 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 1253417..1253686 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 1263507..1263689 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 1325668..1326147 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 1326150..1326605 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 1326571..1327107 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 1327120..1327257 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 1327236..1328282 | + | cystathionine beta-lyase |
| FP | CDS | NC_002163.1 | 1329668..1330444 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 1330451..1331008 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 1403970..1404182 | - | hypothetical protein |
| FP | CDS | NC_002163.1 | 1404311..1405150 | - | hypothetical protein |
| FP | CDS | NC_002163.1 | 1460423..1460707 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 1460737..1461483 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 1490206..1490379 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 1490511..1491152 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 1491162..1491341 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 1597008..1597658 | + | hypothetical protein |
| FP | CDS | NC_002163.1 | 1633895..1634119 | - | hypothetical protein |

> Primeros pares con coordenadas no exactas

| tipo | seqid | ref | Δstart | Δend | solape |
|---|---|---|---|---|---|
| CDS | NC_002163.1 | 6703..8010 | +6 | +0 | 0.9954 |
| CDS | NC_002163.1 | 27402..28304 | +0 | -57 | 0.9369 |
| CDS | NC_002163.1 | 46424..50156 | +0 | -1129 | 0.6976 |
| CDS | NC_002163.1 | 59493..61289 | +3 | +0 | 0.9983 |
| CDS | NC_002163.1 | 70274..71098 | +0 | -36 | 0.9564 |
| CDS | NC_002163.1 | 121680..122369 | +3 | +0 | 0.9957 |
| CDS | NC_002163.1 | 127776..128237 | +0 | -12 | 0.9740 |
| CDS | NC_002163.1 | 136612..136869 | +60 | +0 | 0.7674 |
| CDS | NC_002163.1 | 148819..150600 | +12 | +0 | 0.9933 |
| CDS | NC_002163.1 | 169054..169962 | +0 | -3 | 0.9967 |
| CDS | NC_002163.1 | 172563..172694 | +0 | -6 | 0.9545 |
| CDS | NC_002163.1 | 177995..179200 | +6 | +0 | 0.9950 |
| CDS | NC_002163.1 | 196534..197139 | +0 | +33 | 0.9484 |
| CDS | NC_002163.1 | 199236..201233 | +3 | +0 | 0.9985 |
| CDS | NC_002163.1 | 210767..211612 | +6 | +0 | 0.9929 |
| CDS | NC_002163.1 | 211616..212803 | +6 | +0 | 0.9949 |
| CDS | NC_002163.1 | 214057..215430 | +0 | -3 | 0.9978 |
| CDS | NC_002163.1 | 217475..218038 | +0 | -6 | 0.9894 |
| CDS | NC_002163.1 | 229831..230310 | +3 | +0 | 0.9938 |
| CDS | NC_002163.1 | 232544..233704 | +24 | +0 | 0.9793 |
| CDS | NC_002163.1 | 251784..253037 | +30 | +0 | 0.9761 |
| CDS | NC_002163.1 | 269830..270489 | +3 | +0 | 0.9955 |
| CDS | NC_002163.1 | 339071..339283 | +15 | +0 | 0.9296 |
| CDS | NC_002163.1 | 390162..390497 | +3 | +0 | 0.9911 |
| CDS | NC_002163.1 | 404585..404995 | +3 | +0 | 0.9927 |
| CDS | NC_002163.1 | 409015..409683 | +0 | -3 | 0.9955 |
| CDS | NC_002163.1 | 419857..420426 | +0 | -51 | 0.9105 |
| CDS | NC_002163.1 | 421850..422479 | +0 | +21 | 0.9677 |
| CDS | NC_002163.1 | 423761..424027 | +0 | -24 | 0.9101 |
| CDS | NC_002163.1 | 438992..443128 | +9 | +0 | 0.9978 |
| CDS | NC_002163.1 | 477755..479089 | +15 | +0 | 0.9888 |
| CDS | NC_002163.1 | 490474..491475 | +0 | -222 | 0.7784 |
| CDS | NC_002163.1 | 505212..506921 | +3 | +0 | 0.9982 |
| CDS | NC_002163.1 | 517060..518352 | +6 | +0 | 0.9954 |
| CDS | NC_002163.1 | 526556..527158 | +42 | +0 | 0.9303 |
| CDS | NC_002163.1 | 544189..545331 | +84 | +0 | 0.9265 |
| CDS | NC_002163.1 | 549825..550535 | +3 | +0 | 0.9958 |
| CDS | NC_002163.1 | 550525..550746 | +0 | -6 | 0.9730 |
| CDS | NC_002163.1 | 587868..591303 | +499 | +0 | 0.8548 |
| CDS | NC_002163.1 | 597828..598667 | -39 | +0 | 0.9556 |
| CDS | NC_002163.1 | 622768..623013 | -39 | +0 | 0.8632 |
| CDS | NC_002163.1 | 623410..624138 | -24 | +0 | 0.9681 |
| CDS | NC_002163.1 | 641890..643395 | +3 | +0 | 0.9980 |
| CDS | NC_002163.1 | 649946..650884 | +0 | -6 | 0.9936 |
| CDS | NC_002163.1 | 652530..653918 | +6 | +0 | 0.9957 |
| CDS | NC_002163.1 | 676227..676976 | +0 | -3 | 0.9960 |
| CDS | NC_002163.1 | 683848..684657 | +3 | +0 | 0.9963 |
| CDS | NC_002163.1 | 689656..691584 | +24 | +0 | 0.9876 |
| CDS | NC_002163.1 | 707840..708634 | +3 | +0 | 0.9962 |
| CDS | NC_002163.1 | 714777..716612 | +0 | -3 | 0.9984 |
| CDS | NC_002163.1 | 735708..736490 | +3 | +0 | 0.9962 |
| CDS | NC_002163.1 | 736487..737011 | +3 | +0 | 0.9943 |
| CDS | NC_002163.1 | 738526..738774 | +3 | +0 | 0.9880 |
| CDS | NC_002163.1 | 741167..742441 | +0 | -6 | 0.9953 |
| CDS | NC_002163.1 | 774427..775152 | +3 | +0 | 0.9959 |
| CDS | NC_002163.1 | 797653..798840 | +0 | -27 | 0.9773 |
| CDS | NC_002163.1 | 811319..812119 | +15 | +0 | 0.9813 |
| CDS | NC_002163.1 | 829486..830148 | +0 | +9 | 0.9866 |
| CDS | NC_002163.1 | 840501..840680 | +0 | -33 | 0.8167 |
| CDS | NC_002163.1 | 900703..901044 | +0 | +30 | 0.9194 |
| CDS | NC_002163.1 | 904409..906697 | +297 | +0 | 0.8702 |
| CDS | NC_002163.1 | 910711..911601 | +15 | +0 | 0.9832 |
| CDS | NC_002163.1 | 921397..922158 | +0 | -3 | 0.9961 |
| CDS | NC_002163.1 | 956627..957076 | +0 | -132 | 0.7067 |
| CDS | NC_002163.1 | 987019..989229 | +0 | -3 | 0.9986 |
| CDS | NC_002163.1 | 1011696..1012148 | +66 | +0 | 0.8543 |
| CDS | NC_002163.1 | 1018565..1018960 | +0 | -108 | 0.7273 |
| CDS | NC_002163.1 | 1019786..1020970 | +0 | -12 | 0.9899 |
| CDS | NC_002163.1 | 1027809..1028081 | +0 | -42 | 0.8462 |
| CDS | NC_002163.1 | 1044593..1045219 | +0 | +12 | 0.9812 |
| CDS | NC_002163.1 | 1055060..1055713 | +0 | -45 | 0.9312 |
| CDS | NC_002163.1 | 1062131..1063210 | +0 | -3 | 0.9972 |
| CDS | NC_002163.1 | 1079195..1080041 | +0 | -319 | 0.6234 |
| CDS | NC_002163.1 | 1085907..1088264 | +0 | -6 | 0.9975 |
| CDS | NC_002163.1 | 1104012..1104635 | +0 | -6 | 0.9904 |
| CDS | NC_002163.1 | 1125011..1125907 | +0 | -3 | 0.9967 |
| CDS | NC_002163.1 | 1129228..1130016 | +18 | +0 | 0.9772 |
| CDS | NC_002163.1 | 1194260..1194796 | +0 | -33 | 0.9385 |
| CDS | NC_002163.1 | 1211677..1212912 | +0 | -3 | 0.9976 |
| CDS | NC_002163.1 | 1248624..1249595 | -21 | +0 | 0.9789 |
| CDS | NC_002163.1 | 1253417..1254092 | +316 | +0 | 0.5325 |
| CDS | NC_002163.1 | 1259677..1261638 | +9 | +0 | 0.9954 |
| CDS | NC_002163.1 | 1261635..1263494 | -21 | +0 | 0.9888 |
| CDS | NC_002163.1 | 1263507..1265454 | +166 | +0 | 0.9148 |
| CDS | NC_002163.1 | 1265451..1267337 | -21 | +0 | 0.9890 |
| CDS | NC_002163.1 | 1290427..1290942 | +0 | -24 | 0.9535 |
| CDS | NC_002163.1 | 1293659..1294111 | +0 | -42 | 0.9073 |
| CDS | NC_002163.1 | 1339759..1341147 | +0 | -18 | 0.9870 |
| CDS | NC_002163.1 | 1370667..1372004 | +0 | -3 | 0.9978 |
| CDS | NC_002163.1 | 1372695..1373867 | +0 | +57 | 0.9537 |
| CDS | NC_002163.1 | 1375007..1377337 | +0 | -3 | 0.9987 |
| CDS | NC_002163.1 | 1397726..1398409 | -3 | +0 | 0.9956 |
| CDS | NC_002163.1 | 1423309..1423572 | +0 | +3 | 0.9888 |
| CDS | NC_002163.1 | 1434062..1435549 | +0 | -18 | 0.9879 |
| CDS | NC_002163.1 | 1449311..1450024 | +0 | -12 | 0.9832 |
| CDS | NC_002163.1 | 1471517..1472236 | +0 | -33 | 0.9542 |
| CDS | NC_002163.1 | 1477635..1478084 | -30 | +0 | 0.9375 |
| CDS | NC_002163.1 | 1485626..1487092 | +0 | -15 | 0.9898 |
| CDS | NC_002163.1 | 1492003..1493991 | +33 | +0 | 0.9834 |
| CDS | NC_002163.1 | 1511274..1512068 | +0 | -6 | 0.9925 |

- Aviso: 7 feature(s) GFF multi-segmento fusionados por ID; se comparan por su intervalo envolvente.

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
/opt/conda/bin/python benchmarks/external/external_benchmark.py --genome data/reference_genomes/cjejuni/cjejuni_genomic.fna.gz --reference data/reference_genomes/cjejuni/cjejuni_genomic.gff.gz --tools bakta,prokka,dfast --run --threads 1 --json reports/external_benchmark/benchmark.json --md reports/external_benchmark/benchmark.md --bakta-db /data/bakta_db
```

Para comparar salidas producidas en otro host, usar `--query-gff TOOL=RUTA` sin `--run`; el origen queda como `precomputed`.
