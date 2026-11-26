# Benchmark externo CPU: Bakta, Prokka y DFAST

**Generado:** 2026-11-26T01:06:28Z

**Esquema:** `external-benchmark-run/v1`

**Alcance:** Comparación local de CDS/tRNA/rRNA de Bakta, Prokka y DFAST contra una referencia GFF3/GBFF. CPU, sin descargas de datos.

**Comando reproducible:**

```bash
/opt/conda/bin/python /opt/benchmark/external/external_benchmark.py --genome /work/data/reference_genomes/saureus/saureus_genomic.fna.gz --reference /work/data/reference_genomes/saureus/saureus_genomic.gff.gz --tools all --run --threads 4 --timeout 5400 --bakta-db /data/bakta_db --extra-arg dfast=--dbroot --extra-arg dfast=/data/dfast_db --save-gffs /work/reports/external_benchmark/heldout/saureus/evidence --json /work/reports/external_benchmark/heldout/saureus.json --md /work/reports/external_benchmark/heldout/saureus.md
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

- FASTA: `data/reference_genomes/saureus/saureus_genomic.fna.gz` (SHA256 `f30fbb37ecacc1c6…`, 1 contig, 2821361 bp).
- Referencia: `data/reference_genomes/saureus/saureus_genomic.gff.gz` (formato `gff3`, SHA256 `b39f4aaab06bcdfe…`).
- Conteos de referencia: CDS **2767**, tRNA **59**, rRNA **16**.
- CDS pseudo excluidos en referencia: **0**.
- SHA256 FASTA: `f30fbb37ecacc1c6baca951ada10065a39a29e019c81f664a5d1fca3d8818adb`.
- SHA256 referencia: `b39f4aaab06bcdfecd9080e8605039fc76e18f8f0ecb3043cba02c38b69f1571`.

## Validación

- Inicio: `2026-11-26T00:34:43Z`; fin: `2026-11-26T01:06:28Z`; total: **1905.07 s**.
- JSON declarado: `reports/external_benchmark/heldout/saureus.json`.
- Markdown declarado: `reports/external_benchmark/heldout/saureus.md`.

| check | estado | detalle |
|---|---|---|
| genome_path | ok | FASTA existe y es legible |
| genome_sha256 | ok | SHA256 f30fbb37ecacc1c6baca951ada10065a39a29e019c81f664a5d1fca3d8818adb |
| reference_path | ok | referencia gff3 con 2842 features |
| reference_sha256 | ok | SHA256 b39f4aaab06bcdfecd9080e8605039fc76e18f8f0ecb3043cba02c38b69f1571 |
| reference_features | ok | CDS=2767 tRNA=59 rRNA=16 pseudo_excluidos=0 |
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
| Bakta | executed | ok | 0.9496 | 0.9833 | 1.0000 | 0.9735 | 0.9289 | 0.9507 | 862.0 |
| Prokka | executed | ok | 0.9502 | 0.9752 | 1.0000 | 0.9745 | 0.9286 | 0.9510 | 656.0 |
| DFAST | executed | ok | 0.9502 | 0.9752 | 1.0000 | 0.9745 | 0.9286 | 0.9510 | 348.5 |

La precisión/recall/F1 se calculan contra la referencia local. Si una herramienta no está instalada, su estado es `not_available` y no hay métricas que mostrar.

## Bakta

- Estado: **ok**.
- Fuente: `executed`.
- Bakta DB: `ok` (válida: **sí**, ruta `/data/bakta_db`, versión `0`).
- Versión: `bakta 1.12.2` (estado `ok`).
- Ejecución: **ok**, 862.05 s, RSS máx. 2980.4 MiB (`2026-11-26T00:34:51Z` → `2026-11-26T00:49:13Z`).
- Salida: `/tmp/external_benchmark_bakta_103880wb/out/benchmark.gff3` (3460782 bytes, SHA256 `04e5ed20cae06b3a…`).
- Comando: `/usr/local/bin/bakta --output /tmp/external_benchmark_bakta_103880wb/out --prefix benchmark --threads 4 --force --skip-plot --db /data/bakta_db /tmp/benchmark_genome_o09cjdbt.fna`.
- Anotación comparada: `reports/external_benchmark/heldout/saureus/evidence/bakta/benchmark.gff3` (formato `gff3`, SHA256 `04e5ed20cae06b3a…`).

### Features (Bakta)

| tipo | ref | query | TP | FP | FN | precisión | recall | F1 | exactos | borde distinto | solape medio |
|---|---|---|---|---|---|---|---|---|---|---|---|
| CDS | 2767 | 2635 | 2565 | 70 | 202 | 0.9734 | 0.9270 | 0.9496 | 2274 | 291 | 0.9929 |
| tRNA | 59 | 61 | 59 | 2 | 0 | 0.9672 | 1.0000 | 0.9833 | 56 | 3 | 0.9981 |
| rRNA | 16 | 16 | 16 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 11 | 5 | 0.9994 |

### Atributos (Bakta)

| tipo | pares | product en ambos | product literal | product tokens | gene literal | locus literal |
|---|---|---|---|---|---|---|
| CDS | 2565 | 2565 | 0.1840 | 0.3537 | 0.8611 | 0.0000 |
| tRNA | 59 | 59 | 0.0000 | 0.6455 | — | 0.0000 |
| rRNA | 16 | 16 | 1.0000 | 1.0000 | — | 0.0000 |

### Hipotéticas (Bakta)

| tipo | ref CDS | ref con product | ref hyp | % ref | query CDS | query con product | query hyp | % query |
|---|---|---|---|---|---|---|---|---|
| CDS | 2767 | 2767 | 1496 | 54.1 | 2635 | 2635 | 385 | 14.6 |
| tRNA | 0 | 0 | 0 | — | 0 | 0 | 0 | — |
| rRNA | 0 | 0 | 0 | — | 0 | 0 | 0 | — |

#### Diferencias

- Falsos negativos (referencia sin pareja): **202**.
- Falsos positivos (herramienta sin pareja): **72**.
- Pares solo con hebra distinta: **1**.
- Pares emparejados con algún borde distinto: **299**.

Las listas de diferencias están truncadas a 100 elementos por lista; los conteos son completos.

> Primeras diferencias no emparejadas

| tipo | feature | seqid | coordenadas | hebra | product/locus |
|---|---|---|---|---|---|
| FN | CDS | NC_007795.1 | 33555..34034 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 34375..34509 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 34473..34565 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 37019..37120 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 42281..42382 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 43976..44179 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 93378..93500 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 113507..113710 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 114927..115214 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 166583..166684 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 213824..213931 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 227920..228027 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 229502..229645 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 232547..232645 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 259676..259810 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 263531..263710 | + | NADH dehydrogenase subunit |
| FN | CDS | NC_007795.1 | 302336..302476 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 318835..318927 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 336647..337054 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 363376..363534 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 366986..367195 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 374638..374739 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 391435..391572 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 399186..399287 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 412193..412321 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 418467..418883 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 444428..444529 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 448172..448270 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 469421..469552 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 472097..472207 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 489935..490027 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 542832..543005 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 547608..547865 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 584859..584960 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 589698..589814 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 639816..639908 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 666157..666273 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 668104..668199 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 707687..707893 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 721050..721148 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 722569..722670 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 730469..730573 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 730850..731005 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 752513..752617 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 754050..754157 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 789868..790014 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 791370..791525 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 811872..811985 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 823343..823480 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 826353..826451 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 833949..834059 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 860487..860678 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 942257..942364 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 942903..943004 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 952891..952995 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 963913..964023 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 976868..976990 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 993974..994093 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 1003967..1004122 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 1006242..1006364 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 1025649..1025816 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 1035362..1035565 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 1048947..1049240 | + | iron compound ABC transporter permease |
| FN | CDS | NC_007795.1 | 1050415..1050522 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 1062304..1062399 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 1070176..1070343 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 1071183..1071275 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 1074500..1074745 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 1075291..1075509 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 1151425..1151523 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 1180529..1180621 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 1181622..1181735 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 1244641..1244838 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 1246952..1247062 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 1248246..1248419 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 1249705..1249827 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 1268125..1268229 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 1273637..1273741 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 1279048..1279146 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 1296727..1296825 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 1299635..1299802 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 1325145..1325330 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 1331232..1331339 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 1362554..1362685 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 1373947..1374168 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 1445100..1445414 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 1450026..1450217 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 1450171..1450308 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 1450334..1450450 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 1453447..1453563 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 1458909..1459028 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 1490813..1490917 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 1503842..1504015 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 1528896..1529150 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 1532809..1532901 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 1556581..1556727 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 1562274..1562477 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 1597626..1597754 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 1638847..1638975 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 1661240..1661359 | - | hypothetical protein |
| FP | CDS | NC_007795.1 | 113771..113869 | - | hypothetical protein |
| FP | CDS | NC_007795.1 | 114553..114843 | + | hypothetical protein |
| FP | CDS | NC_007795.1 | 115867..116088 | + | Transposase IS200-like domain-containing protein |
| FP | CDS | NC_007795.1 | 264380..264784 | + | Truncated transposase |
| FP | CDS | NC_007795.1 | 412583..412645 | - | Phenol-soluble modulin alpha 4 peptide |
| FP | CDS | NC_007795.1 | 412709..412777 | - | Phenol-soluble modulin alpha 3 peptide |
| FP | CDS | NC_007795.1 | 412829..412894 | - | Phenol-soluble modulin alpha 2 peptide |
| FP | CDS | NC_007795.1 | 412926..412991 | - | Phenol-soluble modulin alpha 1 peptide |
| FP | CDS | NC_007795.1 | 591277..591744 | + | DUF443 domain-containing protein |
| FP | CDS | NC_007795.1 | 613988..614191 | - | DUF1659 domain-containing protein |
| FP | CDS | NC_007795.1 | 620050..620352 | + | Na+/H+ antiporter Mnh2 subunit F |
| FP | CDS | NC_007795.1 | 796106..797515 | + | secreted von Willebrand factor-binding protein precursor |
| FP | CDS | NC_007795.1 | 889974..890408 | + | Cell surface like-protein Map-w |
| FP | CDS | NC_007795.1 | 947205..947912 | + | hypothetical protein |
| FP | CDS | NC_007795.1 | 947905..948687 | + | hypothetical protein |
| FP | CDS | NC_007795.1 | 1041899..1042030 | + | hypothetical protein |
| FP | CDS | NC_007795.1 | 1111031..1111225 | - | hypothetical protein |
| FP | CDS | NC_007795.1 | 1180820..1180921 | + | hypothetical protein |
| FP | CDS | NC_007795.1 | 1180899..1181588 | + | hypothetical protein |
| FP | CDS | NC_007795.1 | 1245423..1245512 | + | hypothetical protein |
| FP | CDS | NC_007795.1 | 1250775..1251065 | + | hypothetical protein |
| FP | CDS | NC_007795.1 | 1280956..1283799 | + | Nuclease SbcCD subunit C |
| FP | CDS | NC_007795.1 | 1349111..1349503 | - | Transposase DDE domain protein |
| FP | CDS | NC_007795.1 | 1418908..1419024 | - | hypothetical protein |
| FP | CDS | NC_007795.1 | 1495120..1495428 | - | hypothetical protein |
| FP | CDS | NC_007795.1 | 1562458..1562874 | - | hypothetical protein |
| FP | CDS | NC_007795.1 | 1562911..1563921 | - | Membrane peptidase, contains TPR repeat domain |
| FP | CDS | NC_007795.1 | 1631841..1632614 | - | tRNA-uridine 2-sulfurtransferase |
| FP | CDS | NC_007795.1 | 1632623..1632958 | - | hypothetical protein |
| FP | CDS | NC_007795.1 | 1636131..1636235 | - | hypothetical protein |
| FP | CDS | NC_007795.1 | 1636204..1636554 | - | hypothetical protein |
| FP | CDS | NC_007795.1 | 1705674..1706480 | - | IS5/IS1182 family transposase |
| FP | CDS | NC_007795.1 | 1706498..1707325 | - | Truncated transposase, putative |
| FP | CDS | NC_007795.1 | 1795368..1795592 | - | hypothetical protein |
| FP | CDS | NC_007795.1 | 1795589..1796080 | - | Transposase InsH N-terminal domain-containing protein |
| FP | CDS | NC_007795.1 | 1813964..1814788 | + | Truncated transposase, putative |
| FP | CDS | NC_007795.1 | 1814821..1815429 | + | Transposase IS4-like domain-containing protein |
| FP | CDS | NC_007795.1 | 1897329..1898531 | + | ISL3 family IS1181 transposase |
| FP | CDS | NC_007795.1 | 1960237..1960566 | - | hypothetical protein |
| FP | CDS | NC_007795.1 | 2031738..2031938 | + | Phospholipase C |
| FP | CDS | NC_007795.1 | 2059099..2059299 | - | Phage protein |
| FP | CDS | NC_007795.1 | 2065119..2065268 | - | Phage protein |
| FP | CDS | NC_007795.1 | 2069455..2069616 | - | hypothetical protein |
| FP | CDS | NC_007795.1 | 2069848..2070204 | - | hypothetical protein |
| FP | CDS | NC_007795.1 | 2074631..2075455 | + | Phospholipase C |
| FP | CDS | NC_007795.1 | 2133600..2134388 | - | hypothetical protein |
| FP | CDS | NC_007795.1 | 2134363..2134590 | - | hypothetical protein |
| FP | CDS | NC_007795.1 | 2212638..2212841 | - | hypothetical protein |
| FP | CDS | NC_007795.1 | 2237284..2238486 | + | ISL3 family IS1181 transposase |
| FP | CDS | NC_007795.1 | 2244932..2245879 | - | Iron-sulfur cluster carrier protein |
| FP | CDS | NC_007795.1 | 2264437..2265237 | + | Transposase |
| FP | CDS | NC_007795.1 | 2265255..2265554 | + | hypothetical protein |
| FP | CDS | NC_007795.1 | 2265717..2266061 | + | Truncated transposase |
| FP | CDS | NC_007795.1 | 2289534..2289959 | + | hypothetical protein |
| FP | CDS | NC_007795.1 | 2349715..2349897 | - | hypothetical protein |
| FP | CDS | NC_007795.1 | 2350014..2350880 | - | hypothetical protein |
| FP | CDS | NC_007795.1 | 2390182..2391495 | - | Na+/H+ antiporter NhaC-like C-terminal domain-containing protein |
| FP | CDS | NC_007795.1 | 2420344..2420622 | - | hypothetical protein |
| FP | CDS | NC_007795.1 | 2420613..2421551 | - | hypothetical protein |
| FP | CDS | NC_007795.1 | 2423442..2423579 | - | hypothetical protein |
| FP | CDS | NC_007795.1 | 2423661..2423897 | - | hypothetical protein |
| FP | CDS | NC_007795.1 | 2425635..2425991 | + | Truncated transposase |
| FP | CDS | NC_007795.1 | 2457007..2457645 | + | hypothetical protein |
| FP | CDS | NC_007795.1 | 2457715..2457921 | + | hypothetical protein |
| FP | CDS | NC_007795.1 | 2558843..2559955 | - | hypothetical protein |
| FP | CDS | NC_007795.1 | 2591200..2591868 | - | hypothetical protein |
| FP | CDS | NC_007795.1 | 2591843..2593024 | - | hypothetical protein |
| FP | CDS | NC_007795.1 | 2671337..2671855 | - | hypothetical protein |
| FP | CDS | NC_007795.1 | 2671852..2672676 | - | hypothetical protein |
| FP | CDS | NC_007795.1 | 2721168..2721272 | - | Uncharacterized protein SA2420-1 |
| FP | tRNA | NC_007795.1 | 1899027..1899101 | - | tRNA-Gly(tcc) |
| FP | tRNA | NC_007795.1 | 1899115..1899189 | - | tRNA-Gly(tcc) |

> Primeros pares con coordenadas no exactas

| tipo | seqid | ref | Δstart | Δend | solape |
|---|---|---|---|---|---|
| CDS | NC_007795.1 | 3670..3915 | +9 | +0 | 0.9634 |
| CDS | NC_007795.1 | 9755..10456 | +0 | +111 | 0.8635 |
| CDS | NC_007795.1 | 17365..18303 | +15 | +0 | 0.9840 |
| CDS | NC_007795.1 | 27515..28798 | -51 | +0 | 0.9618 |
| CDS | NC_007795.1 | 42681..43613 | +0 | +54 | 0.9453 |
| CDS | NC_007795.1 | 48313..49056 | +9 | +0 | 0.9879 |
| CDS | NC_007795.1 | 56284..57054 | -33 | +0 | 0.9590 |
| CDS | NC_007795.1 | 58855..61053 | -39 | +0 | 0.9826 |
| CDS | NC_007795.1 | 64265..65926 | +0 | +6 | 0.9964 |
| CDS | NC_007795.1 | 76521..77435 | +0 | +84 | 0.9159 |
| CDS | NC_007795.1 | 79885..80730 | -135 | +0 | 0.8624 |
| CDS | NC_007795.1 | 83505..84761 | +84 | +0 | 0.9332 |
| CDS | NC_007795.1 | 86396..88246 | +72 | +0 | 0.9611 |
| CDS | NC_007795.1 | 88221..88997 | -3 | +0 | 0.9962 |
| CDS | NC_007795.1 | 93570..94550 | +9 | +0 | 0.9908 |
| CDS | NC_007795.1 | 101370..102125 | +0 | -15 | 0.9802 |
| CDS | NC_007795.1 | 103075..104427 | -6 | +0 | 0.9956 |
| CDS | NC_007795.1 | 144157..145188 | -3 | +0 | 0.9971 |
| CDS | NC_007795.1 | 157205..157684 | +0 | +15 | 0.9697 |
| CDS | NC_007795.1 | 159987..161018 | +0 | +3 | 0.9971 |
| CDS | NC_007795.1 | 164125..164622 | +0 | +60 | 0.8925 |
| CDS | NC_007795.1 | 166663..168660 | +0 | +48 | 0.9765 |
| CDS | NC_007795.1 | 169309..170349 | -15 | +0 | 0.9858 |
| CDS | NC_007795.1 | 172714..173589 | -3 | +0 | 0.9966 |
| CDS | NC_007795.1 | 181951..183429 | +45 | +0 | 0.9696 |
| CDS | NC_007795.1 | 199583..200155 | +0 | -78 | 0.8639 |
| CDS | NC_007795.1 | 220761..222338 | +0 | -15 | 0.9905 |
| CDS | NC_007795.1 | 222779..223921 | -87 | +0 | 0.9293 |
| CDS | NC_007795.1 | 224235..225647 | +0 | +63 | 0.9573 |
| CDS | NC_007795.1 | 228482..229291 | -144 | +0 | 0.8491 |
| CDS | NC_007795.1 | 239969..240643 | -42 | +0 | 0.9414 |
| CDS | NC_007795.1 | 241836..243377 | -153 | +0 | 0.9097 |
| CDS | NC_007795.1 | 251527..253266 | -15 | +0 | 0.9915 |
| CDS | NC_007795.1 | 254122..254565 | +6 | +0 | 0.9865 |
| CDS | NC_007795.1 | 274790..275683 | +0 | +9 | 0.9900 |
| CDS | NC_007795.1 | 279728..280009 | +39 | +0 | 0.8617 |
| CDS | NC_007795.1 | 286569..287243 | -3 | +0 | 0.9956 |
| CDS | NC_007795.1 | 287243..287560 | +6 | +0 | 0.9811 |
| CDS | NC_007795.1 | 287570..289414 | -9 | +0 | 0.9951 |
| CDS | NC_007795.1 | 297408..298685 | +0 | +30 | 0.9771 |
| CDS | NC_007795.1 | 302466..303416 | -81 | +0 | 0.9215 |
| CDS | NC_007795.1 | 312595..313908 | +0 | +48 | 0.9648 |
| CDS | NC_007795.1 | 339971..341653 | -30 | +0 | 0.9825 |
| CDS | NC_007795.1 | 341746..342402 | +0 | -51 | 0.9224 |
| CDS | NC_007795.1 | 342464..342688 | +0 | -9 | 0.9600 |
| CDS | NC_007795.1 | 362010..362948 | +0 | +36 | 0.9631 |
| CDS | NC_007795.1 | 362987..363298 | +0 | -21 | 0.9327 |
| CDS | NC_007795.1 | 368020..368403 | +0 | +3 | 0.9922 |
| CDS | NC_007795.1 | 369357..369650 | +0 | -36 | 0.8776 |
| CDS | NC_007795.1 | 372329..373084 | -51 | +0 | 0.9368 |
| CDS | NC_007795.1 | 399449..400921 | -84 | +0 | 0.9461 |
| CDS | NC_007795.1 | 407345..408082 | -78 | +0 | 0.9044 |
| CDS | NC_007795.1 | 410343..410603 | +0 | -66 | 0.7471 |
| CDS | NC_007795.1 | 431034..432152 | +0 | -6 | 0.9946 |
| CDS | NC_007795.1 | 445284..446957 | -24 | +0 | 0.9859 |
| CDS | NC_007795.1 | 457344..458183 | -87 | +0 | 0.9061 |
| CDS | NC_007795.1 | 503169..503630 | -9 | +0 | 0.9809 |
| CDS | NC_007795.1 | 503727..504215 | -78 | +0 | 0.8624 |
| CDS | NC_007795.1 | 504205..505212 | -3 | +0 | 0.9970 |
| CDS | NC_007795.1 | 509574..510629 | -18 | +0 | 0.9832 |
| CDS | NC_007795.1 | 513063..513710 | -33 | +0 | 0.9515 |
| CDS | NC_007795.1 | 525875..529471 | -27 | +0 | 0.9925 |
| CDS | NC_007795.1 | 574840..575490 | +0 | +114 | 0.8510 |
| CDS | NC_007795.1 | 578134..578898 | -72 | +0 | 0.9140 |
| CDS | NC_007795.1 | 584361..584777 | +0 | +24 | 0.9456 |
| CDS | NC_007795.1 | 586908..587465 | -69 | +0 | 0.8900 |
| CDS | NC_007795.1 | 588392..588937 | -72 | +0 | 0.8835 |
| CDS | NC_007795.1 | 590589..591230 | +9 | +0 | 0.9860 |
| CDS | NC_007795.1 | 593183..593752 | -72 | +0 | 0.8879 |
| CDS | NC_007795.1 | 594375..594635 | +15 | +0 | 0.9425 |
| CDS | NC_007795.1 | 599248..599769 | +15 | +0 | 0.9713 |
| CDS | NC_007795.1 | 601586..601969 | -33 | +0 | 0.9209 |
| CDS | NC_007795.1 | 605090..605917 | -60 | +0 | 0.9324 |
| CDS | NC_007795.1 | 606065..607015 | -36 | +0 | 0.9635 |
| CDS | NC_007795.1 | 612623..613516 | +0 | +36 | 0.9613 |
| CDS | NC_007795.1 | 620306..620764 | +21 | +0 | 0.9542 |
| CDS | NC_007795.1 | 621143..623146 | -39 | +0 | 0.9809 |
| CDS | NC_007795.1 | 685710..686057 | +0 | -63 | 0.8190 |
| CDS | NC_007795.1 | 712364..714247 | +6 | +0 | 0.9968 |
| CDS | NC_007795.1 | 725114..725848 | +0 | +132 | 0.8478 |
| CDS | NC_007795.1 | 727082..729247 | +60 | +0 | 0.9723 |
| CDS | NC_007795.1 | 731044..732015 | -42 | +0 | 0.9586 |
| CDS | NC_007795.1 | 747288..748226 | -144 | +0 | 0.8670 |
| CDS | NC_007795.1 | 753021..753896 | -117 | +0 | 0.8822 |
| CDS | NC_007795.1 | 757511..758527 | +36 | +0 | 0.9646 |
| CDS | NC_007795.1 | 773243..774145 | +0 | -33 | 0.9635 |
| CDS | NC_007795.1 | 804682..805242 | -21 | +0 | 0.9639 |
| CDS | NC_007795.1 | 812195..812581 | -3 | +0 | 0.9923 |
| CDS | NC_007795.1 | 819904..821154 | +9 | +0 | 0.9928 |
| CDS | NC_007795.1 | 832811..833695 | -75 | +0 | 0.9219 |
| CDS | NC_007795.1 | 850790..851119 | +0 | +12 | 0.9649 |
| CDS | NC_007795.1 | 855577..855921 | -33 | +0 | 0.9127 |
| CDS | NC_007795.1 | 867744..868211 | -108 | +0 | 0.8125 |
| CDS | NC_007795.1 | 879806..880072 | -42 | +0 | 0.8641 |
| CDS | NC_007795.1 | 912569..912886 | +0 | +48 | 0.8689 |
| CDS | NC_007795.1 | 925707..926915 | +0 | -18 | 0.9851 |
| CDS | NC_007795.1 | 929948..930208 | +9 | +0 | 0.9655 |
| CDS | NC_007795.1 | 933094..935418 | +15 | +0 | 0.9935 |
| CDS | NC_007795.1 | 936932..938443 | -3 | +0 | 0.9980 |
| CDS | NC_007795.1 | 938930..939460 | +0 | +39 | 0.9316 |


## Prokka

- Estado: **ok**.
- Fuente: `executed`.
- Versión: `prokka 1.15.6` (estado `ok`).
- Ejecución: **ok**, 655.99 s, RSS máx. 3069.6 MiB (`2026-11-26T00:49:23Z` → `2026-11-26T01:00:19Z`).
- Salida: `/tmp/external_benchmark_prokka_vg2a8vuq/out/benchmark.gff` (3492972 bytes, SHA256 `eff69adac55d031c…`).
- Comando: `/usr/local/bin/prokka --outdir /tmp/external_benchmark_prokka_vg2a8vuq/out --prefix benchmark --cpus 4 --force --quiet /tmp/benchmark_genome_o09cjdbt.fna`.
- Anotación comparada: `reports/external_benchmark/heldout/saureus/evidence/prokka/benchmark.gff` (formato `gff3`, SHA256 `eff69adac55d031c…`).

### Features (Prokka)

| tipo | ref | query | TP | FP | FN | precisión | recall | F1 | exactos | borde distinto | solape medio |
|---|---|---|---|---|---|---|---|---|---|---|---|
| CDS | 2767 | 2630 | 2564 | 66 | 203 | 0.9749 | 0.9266 | 0.9502 | 2273 | 291 | 0.9929 |
| tRNA | 59 | 62 | 59 | 3 | 0 | 0.9516 | 1.0000 | 0.9752 | 43 | 16 | 0.9934 |
| rRNA | 16 | 16 | 16 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 11 | 5 | 0.9994 |

### Atributos (Prokka)

| tipo | pares | product en ambos | product literal | product tokens | gene literal | locus literal |
|---|---|---|---|---|---|---|
| CDS | 2564 | 2564 | 0.3783 | 0.5075 | 0.8303 | 0.0000 |
| tRNA | 59 | 59 | 0.0000 | 0.6667 | — | 0.0000 |
| rRNA | 16 | 16 | 1.0000 | 1.0000 | — | 0.0000 |

### Hipotéticas (Prokka)

| tipo | ref CDS | ref con product | ref hyp | % ref | query CDS | query con product | query hyp | % query |
|---|---|---|---|---|---|---|---|---|
| CDS | 2767 | 2767 | 1496 | 54.1 | 2630 | 2630 | 863 | 32.8 |
| tRNA | 0 | 0 | 0 | — | 0 | 0 | 0 | — |
| rRNA | 0 | 0 | 0 | — | 0 | 0 | 0 | — |

#### Diferencias

- Falsos negativos (referencia sin pareja): **203**.
- Falsos positivos (herramienta sin pareja): **69**.
- Pares solo con hebra distinta: **1**.
- Pares emparejados con algún borde distinto: **312**.

Las listas de diferencias están truncadas a 100 elementos por lista; los conteos son completos.

> Primeras diferencias no emparejadas

| tipo | feature | seqid | coordenadas | hebra | product/locus |
|---|---|---|---|---|---|
| FN | CDS | NC_007795.1 | 33555..34034 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 34375..34509 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 34473..34565 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 37019..37120 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 42281..42382 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 43976..44179 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 93378..93500 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 113507..113710 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 114927..115214 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 166583..166684 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 213824..213931 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 227920..228027 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 229502..229645 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 232547..232645 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 259676..259810 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 263531..263710 | + | NADH dehydrogenase subunit |
| FN | CDS | NC_007795.1 | 302336..302476 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 318835..318927 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 336647..337054 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 363376..363534 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 366986..367195 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 374638..374739 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 391435..391572 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 399186..399287 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 412193..412321 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 418467..418883 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 444428..444529 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 448172..448270 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 469421..469552 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 472097..472207 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 489935..490027 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 542832..543005 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 547608..547865 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 584859..584960 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 589698..589814 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 639816..639908 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 666157..666273 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 668104..668199 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 707687..707893 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 721050..721148 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 722569..722670 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 730469..730573 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 730850..731005 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 752513..752617 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 754050..754157 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 789868..790014 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 791370..791525 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 811872..811985 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 823343..823480 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 826353..826451 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 833949..834059 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 860487..860678 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 942257..942364 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 942903..943004 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 952891..952995 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 963913..964023 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 976868..976990 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 993974..994093 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 1003967..1004122 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 1006242..1006364 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 1025649..1025816 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 1035362..1035565 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 1048947..1049240 | + | iron compound ABC transporter permease |
| FN | CDS | NC_007795.1 | 1050415..1050522 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 1062304..1062399 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 1070176..1070343 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 1071183..1071275 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 1074500..1074745 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 1075291..1075509 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 1151425..1151523 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 1180529..1180621 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 1181622..1181735 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 1244641..1244838 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 1246952..1247062 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 1248246..1248419 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 1249705..1249827 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 1268125..1268229 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 1273637..1273741 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 1279048..1279146 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 1296727..1296825 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 1299635..1299802 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 1325145..1325330 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 1331232..1331339 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 1362554..1362685 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 1373947..1374168 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 1445100..1445414 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 1450026..1450217 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 1450171..1450308 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 1450334..1450450 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 1453447..1453563 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 1458909..1459028 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 1490813..1490917 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 1503842..1504015 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 1528896..1529150 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 1532809..1532901 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 1556581..1556727 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 1562274..1562477 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 1597626..1597754 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 1638847..1638975 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 1661240..1661359 | - | hypothetical protein |
| FP | CDS | NC_007795.1 | 113771..113869 | - | hypothetical protein |
| FP | CDS | NC_007795.1 | 114553..114843 | + | Replication initiation protein |
| FP | CDS | NC_007795.1 | 115867..116088 | + | IS200/IS605 family transposase ISSep3 |
| FP | CDS | NC_007795.1 | 264380..264784 | + | hypothetical protein |
| FP | CDS | NC_007795.1 | 591277..591744 | + | hypothetical protein |
| FP | CDS | NC_007795.1 | 613988..614191 | - | hypothetical protein |
| FP | CDS | NC_007795.1 | 620050..620352 | + | Na(+)/H(+) antiporter subunit F |
| FP | CDS | NC_007795.1 | 796106..797515 | + | hypothetical protein |
| FP | CDS | NC_007795.1 | 889974..890408 | + | Protein map |
| FP | CDS | NC_007795.1 | 947205..947912 | + | Poly(ribitol-phosphate) alpha-N-acetylglucosaminyltransferase |
| FP | CDS | NC_007795.1 | 947905..948687 | + | Poly(ribitol-phosphate) alpha-N-acetylglucosaminyltransferase |
| FP | CDS | NC_007795.1 | 1041899..1042030 | + | Large ribosomal subunit protein bL32 |
| FP | CDS | NC_007795.1 | 1111031..1111225 | - | hypothetical protein |
| FP | CDS | NC_007795.1 | 1180820..1180921 | + | Small ribosomal subunit protein uS2 |
| FP | CDS | NC_007795.1 | 1180899..1181588 | + | Small ribosomal subunit protein uS2 |
| FP | CDS | NC_007795.1 | 1245423..1245512 | + | hypothetical protein |
| FP | CDS | NC_007795.1 | 1250775..1251065 | + | hypothetical protein |
| FP | CDS | NC_007795.1 | 1280956..1283799 | + | Nuclease SbcCD subunit C |
| FP | CDS | NC_007795.1 | 1349111..1349503 | - | hypothetical protein |
| FP | CDS | NC_007795.1 | 1418908..1419024 | - | hypothetical protein |
| FP | CDS | NC_007795.1 | 1495120..1495428 | - | hypothetical protein |
| FP | CDS | NC_007795.1 | 1562458..1562874 | - | hypothetical protein |
| FP | CDS | NC_007795.1 | 1562911..1563921 | - | Rhomboid protease GluP |
| FP | CDS | NC_007795.1 | 1631841..1632614 | - | tRNA-specific 2-thiouridylase MnmA |
| FP | CDS | NC_007795.1 | 1632623..1632958 | - | tRNA-specific 2-thiouridylase MnmA |
| FP | CDS | NC_007795.1 | 1636131..1636235 | - | hypothetical protein |
| FP | CDS | NC_007795.1 | 1636204..1636554 | - | HTH-type transcriptional regulator CymR |
| FP | CDS | NC_007795.1 | 1705674..1706480 | - | IS1182 family transposase ISSau3 |
| FP | CDS | NC_007795.1 | 1706498..1707325 | - | IS1182 family transposase ISSau3 |
| FP | CDS | NC_007795.1 | 1795368..1795592 | - | hypothetical protein |
| FP | CDS | NC_007795.1 | 1795589..1796080 | - | IS1182 family transposase ISSau3 |
| FP | CDS | NC_007795.1 | 1813964..1814788 | + | IS1182 family transposase ISSau3 |
| FP | CDS | NC_007795.1 | 1814821..1815429 | + | IS1182 family transposase ISSau3 |
| FP | CDS | NC_007795.1 | 1897329..1898531 | + | ISL3 family transposase IS1181 |
| FP | CDS | NC_007795.1 | 1960237..1960566 | - | hypothetical protein |
| FP | CDS | NC_007795.1 | 2031738..2031938 | + | Phospholipase C |
| FP | CDS | NC_007795.1 | 2059099..2059299 | - | hypothetical protein |
| FP | CDS | NC_007795.1 | 2065119..2065268 | - | hypothetical protein |
| FP | CDS | NC_007795.1 | 2069455..2069616 | - | hypothetical protein |
| FP | CDS | NC_007795.1 | 2069848..2070204 | - | hypothetical protein |
| FP | CDS | NC_007795.1 | 2074631..2075455 | + | Phospholipase C |
| FP | CDS | NC_007795.1 | 2133600..2134388 | - | Phosphoserine phosphatase RsbU |
| FP | CDS | NC_007795.1 | 2134363..2134590 | - | hypothetical protein |
| FP | CDS | NC_007795.1 | 2212638..2212841 | - | hypothetical protein |
| FP | CDS | NC_007795.1 | 2237284..2238486 | + | ISL3 family transposase IS1181 |
| FP | CDS | NC_007795.1 | 2244932..2245879 | - | Iron-sulfur cluster carrier protein |
| FP | CDS | NC_007795.1 | 2264437..2265237 | + | hypothetical protein |
| FP | CDS | NC_007795.1 | 2265255..2265554 | + | IS1182 family transposase ISSau3 |
| FP | CDS | NC_007795.1 | 2265717..2266061 | + | IS1182 family transposase ISSau3 |
| FP | CDS | NC_007795.1 | 2289534..2289959 | + | 65 kDa membrane protein |
| FP | CDS | NC_007795.1 | 2349715..2349897 | - | Putative acyl-CoA dehydrogenase YdbM |
| FP | CDS | NC_007795.1 | 2350014..2350880 | - | Putative acyl-CoA dehydrogenase YdbM |
| FP | CDS | NC_007795.1 | 2390182..2391495 | - | Malate-2H(+)/Na(+)-lactate antiporter |
| FP | CDS | NC_007795.1 | 2420344..2420622 | - | hypothetical protein |
| FP | CDS | NC_007795.1 | 2420613..2421551 | - | Bicyclomycin resistance protein |
| FP | CDS | NC_007795.1 | 2423442..2423579 | - | HTH-type transcriptional regulator TcaR |
| FP | CDS | NC_007795.1 | 2423661..2423897 | - | HTH-type transcriptional regulator TcaR |
| FP | CDS | NC_007795.1 | 2425635..2425991 | + | hypothetical protein |
| FP | CDS | NC_007795.1 | 2457007..2457645 | + | hypothetical protein |
| FP | CDS | NC_007795.1 | 2457715..2457921 | + | hypothetical protein |
| FP | CDS | NC_007795.1 | 2558843..2559955 | - | hypothetical protein |
| FP | CDS | NC_007795.1 | 2591200..2591868 | - | hypothetical protein |
| FP | CDS | NC_007795.1 | 2591843..2593024 | - | hypothetical protein |
| FP | CDS | NC_007795.1 | 2671337..2671855 | - | hypothetical protein |
| FP | CDS | NC_007795.1 | 2671852..2672676 | - | hypothetical protein |
| FP | CDS | NC_007795.1 | 2721168..2721272 | - | putative protein.1 |
| FP | tRNA | NC_007795.1 | 1899027..1899101 | - | tRNA-Gly(tcc) |
| FP | tRNA | NC_007795.1 | 1899115..1899189 | - | tRNA-Gly(tcc) |
| FP | tRNA | NC_007795.1 | 2340343..2340435 | - | tRNA-Gln(ttg) |

> Primeros pares con coordenadas no exactas

| tipo | seqid | ref | Δstart | Δend | solape |
|---|---|---|---|---|---|
| CDS | NC_007795.1 | 3670..3915 | +9 | +0 | 0.9634 |
| CDS | NC_007795.1 | 9755..10456 | +0 | +111 | 0.8635 |
| CDS | NC_007795.1 | 17365..18303 | +15 | +0 | 0.9840 |
| CDS | NC_007795.1 | 27515..28798 | -51 | +0 | 0.9618 |
| CDS | NC_007795.1 | 42681..43613 | +0 | +54 | 0.9453 |
| CDS | NC_007795.1 | 48313..49056 | +9 | +0 | 0.9879 |
| CDS | NC_007795.1 | 56284..57054 | -33 | +0 | 0.9590 |
| CDS | NC_007795.1 | 58855..61053 | -39 | +0 | 0.9826 |
| CDS | NC_007795.1 | 64265..65926 | +0 | +6 | 0.9964 |
| CDS | NC_007795.1 | 76521..77435 | +0 | +84 | 0.9159 |
| CDS | NC_007795.1 | 79885..80730 | -135 | +0 | 0.8624 |
| CDS | NC_007795.1 | 83505..84761 | +84 | +0 | 0.9332 |
| CDS | NC_007795.1 | 86396..88246 | +72 | +0 | 0.9611 |
| CDS | NC_007795.1 | 88221..88997 | -3 | +0 | 0.9962 |
| CDS | NC_007795.1 | 93570..94550 | +9 | +0 | 0.9908 |
| CDS | NC_007795.1 | 101370..102125 | +0 | -15 | 0.9802 |
| CDS | NC_007795.1 | 103075..104427 | -6 | +0 | 0.9956 |
| CDS | NC_007795.1 | 144157..145188 | -3 | +0 | 0.9971 |
| CDS | NC_007795.1 | 157205..157684 | +0 | +15 | 0.9697 |
| CDS | NC_007795.1 | 159987..161018 | +0 | +3 | 0.9971 |
| CDS | NC_007795.1 | 164125..164622 | +0 | +60 | 0.8925 |
| CDS | NC_007795.1 | 166663..168660 | +0 | +48 | 0.9765 |
| CDS | NC_007795.1 | 169309..170349 | -15 | +0 | 0.9858 |
| CDS | NC_007795.1 | 172714..173589 | -3 | +0 | 0.9966 |
| CDS | NC_007795.1 | 181951..183429 | +45 | +0 | 0.9696 |
| CDS | NC_007795.1 | 199583..200155 | +0 | -78 | 0.8639 |
| CDS | NC_007795.1 | 220761..222338 | +0 | -15 | 0.9905 |
| CDS | NC_007795.1 | 222779..223921 | -87 | +0 | 0.9293 |
| CDS | NC_007795.1 | 224235..225647 | +0 | +63 | 0.9573 |
| CDS | NC_007795.1 | 228482..229291 | -144 | +0 | 0.8491 |
| CDS | NC_007795.1 | 239969..240643 | -42 | +0 | 0.9414 |
| CDS | NC_007795.1 | 241836..243377 | -153 | +0 | 0.9097 |
| CDS | NC_007795.1 | 251527..253266 | -15 | +0 | 0.9915 |
| CDS | NC_007795.1 | 254122..254565 | +6 | +0 | 0.9865 |
| CDS | NC_007795.1 | 274790..275683 | +0 | +9 | 0.9900 |
| CDS | NC_007795.1 | 279728..280009 | +39 | +0 | 0.8617 |
| CDS | NC_007795.1 | 286569..287243 | -3 | +0 | 0.9956 |
| CDS | NC_007795.1 | 287243..287560 | +6 | +0 | 0.9811 |
| CDS | NC_007795.1 | 287570..289414 | -9 | +0 | 0.9951 |
| CDS | NC_007795.1 | 297408..298685 | +0 | +30 | 0.9771 |
| CDS | NC_007795.1 | 302466..303416 | -81 | +0 | 0.9215 |
| CDS | NC_007795.1 | 312595..313908 | +0 | +48 | 0.9648 |
| CDS | NC_007795.1 | 339971..341653 | -30 | +0 | 0.9825 |
| CDS | NC_007795.1 | 341746..342402 | +0 | -51 | 0.9224 |
| CDS | NC_007795.1 | 342464..342688 | +0 | -9 | 0.9600 |
| CDS | NC_007795.1 | 362010..362948 | +0 | +36 | 0.9631 |
| CDS | NC_007795.1 | 362987..363298 | +0 | -21 | 0.9327 |
| CDS | NC_007795.1 | 368020..368403 | +0 | +3 | 0.9922 |
| CDS | NC_007795.1 | 369357..369650 | +0 | -36 | 0.8776 |
| CDS | NC_007795.1 | 372329..373084 | -51 | +0 | 0.9368 |
| CDS | NC_007795.1 | 399449..400921 | -84 | +0 | 0.9461 |
| CDS | NC_007795.1 | 407345..408082 | -78 | +0 | 0.9044 |
| CDS | NC_007795.1 | 410343..410603 | +0 | -66 | 0.7471 |
| CDS | NC_007795.1 | 431034..432152 | +0 | -6 | 0.9946 |
| CDS | NC_007795.1 | 445284..446957 | -24 | +0 | 0.9859 |
| CDS | NC_007795.1 | 457344..458183 | -87 | +0 | 0.9061 |
| CDS | NC_007795.1 | 503169..503630 | -9 | +0 | 0.9809 |
| CDS | NC_007795.1 | 503727..504215 | -78 | +0 | 0.8624 |
| CDS | NC_007795.1 | 504205..505212 | -3 | +0 | 0.9970 |
| CDS | NC_007795.1 | 509574..510629 | -18 | +0 | 0.9832 |
| CDS | NC_007795.1 | 513063..513710 | -33 | +0 | 0.9515 |
| CDS | NC_007795.1 | 525875..529471 | -27 | +0 | 0.9925 |
| CDS | NC_007795.1 | 574840..575490 | +0 | +114 | 0.8510 |
| CDS | NC_007795.1 | 578134..578898 | -72 | +0 | 0.9140 |
| CDS | NC_007795.1 | 584361..584777 | +0 | +24 | 0.9456 |
| CDS | NC_007795.1 | 586908..587465 | -69 | +0 | 0.8900 |
| CDS | NC_007795.1 | 588392..588937 | -72 | +0 | 0.8835 |
| CDS | NC_007795.1 | 590589..591230 | +9 | +0 | 0.9860 |
| CDS | NC_007795.1 | 593183..593752 | -72 | +0 | 0.8879 |
| CDS | NC_007795.1 | 594375..594635 | +15 | +0 | 0.9425 |
| CDS | NC_007795.1 | 599248..599769 | +15 | +0 | 0.9713 |
| CDS | NC_007795.1 | 601586..601969 | -33 | +0 | 0.9209 |
| CDS | NC_007795.1 | 605090..605917 | -60 | +0 | 0.9324 |
| CDS | NC_007795.1 | 606065..607015 | -36 | +0 | 0.9635 |
| CDS | NC_007795.1 | 612623..613516 | +0 | +36 | 0.9613 |
| CDS | NC_007795.1 | 620306..620764 | +21 | +0 | 0.9542 |
| CDS | NC_007795.1 | 621143..623146 | -39 | +0 | 0.9809 |
| CDS | NC_007795.1 | 685710..686057 | +0 | -63 | 0.8190 |
| CDS | NC_007795.1 | 712364..714247 | +6 | +0 | 0.9968 |
| CDS | NC_007795.1 | 725114..725848 | +0 | +132 | 0.8478 |
| CDS | NC_007795.1 | 727082..729247 | +60 | +0 | 0.9723 |
| CDS | NC_007795.1 | 731044..732015 | -42 | +0 | 0.9586 |
| CDS | NC_007795.1 | 747288..748226 | -144 | +0 | 0.8670 |
| CDS | NC_007795.1 | 753021..753896 | -117 | +0 | 0.8822 |
| CDS | NC_007795.1 | 757511..758527 | +36 | +0 | 0.9646 |
| CDS | NC_007795.1 | 773243..774145 | +0 | -33 | 0.9635 |
| CDS | NC_007795.1 | 804682..805242 | -21 | +0 | 0.9639 |
| CDS | NC_007795.1 | 812195..812581 | -3 | +0 | 0.9923 |
| CDS | NC_007795.1 | 819904..821154 | +9 | +0 | 0.9928 |
| CDS | NC_007795.1 | 832811..833695 | -75 | +0 | 0.9219 |
| CDS | NC_007795.1 | 850790..851119 | +0 | +12 | 0.9649 |
| CDS | NC_007795.1 | 855577..855921 | -33 | +0 | 0.9127 |
| CDS | NC_007795.1 | 867744..868211 | -108 | +0 | 0.8125 |
| CDS | NC_007795.1 | 879806..880072 | -42 | +0 | 0.8641 |
| CDS | NC_007795.1 | 912569..912886 | +0 | +48 | 0.8689 |
| CDS | NC_007795.1 | 925707..926915 | +0 | -18 | 0.9851 |
| CDS | NC_007795.1 | 929948..930208 | +9 | +0 | 0.9655 |
| CDS | NC_007795.1 | 933094..935418 | +15 | +0 | 0.9935 |
| CDS | NC_007795.1 | 936932..938443 | -3 | +0 | 0.9980 |
| CDS | NC_007795.1 | 938930..939460 | +0 | +39 | 0.9316 |


## DFAST

- Estado: **ok**.
- Fuente: `executed`.
- Versión: `DFAST ver. 1.6.2` (estado `ok`).
- Ejecución: **ok**, 348.46 s, RSS máx. 3069.6 MiB (`2026-11-26T01:00:30Z` → `2026-11-26T01:06:19Z`).
- Salida: `/tmp/external_benchmark_dfast_e9igq4wx/out/genome.gff` (4956296 bytes, SHA256 `a004c09c62e70480…`).
- Comando: `/usr/local/bin/dfast --genome /tmp/benchmark_genome_o09cjdbt.fna --out /tmp/external_benchmark_dfast_e9igq4wx/out --cpu 4 --force --dbroot /data/dfast_db`.
- Anotación comparada: `reports/external_benchmark/heldout/saureus/evidence/dfast/genome.gff` (formato `gff3`, SHA256 `a004c09c62e70480…`).

### Features (DFAST)

| tipo | ref | query | TP | FP | FN | precisión | recall | F1 | exactos | borde distinto | solape medio |
|---|---|---|---|---|---|---|---|---|---|---|---|
| CDS | 2767 | 2630 | 2564 | 66 | 203 | 0.9749 | 0.9266 | 0.9502 | 2273 | 291 | 0.9929 |
| tRNA | 59 | 62 | 59 | 3 | 0 | 0.9516 | 1.0000 | 0.9752 | 43 | 16 | 0.9934 |
| rRNA | 16 | 16 | 16 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 0 | 16 | 0.9916 |

### Atributos (DFAST)

| tipo | pares | product en ambos | product literal | product tokens | gene literal | locus literal |
|---|---|---|---|---|---|---|
| CDS | 2564 | 2564 | 0.2539 | 0.4093 | 0.9301 | 0.0000 |
| tRNA | 59 | 59 | 1.0000 | 1.0000 | — | 0.0000 |
| rRNA | 16 | 16 | 1.0000 | 1.0000 | — | 0.0000 |

### Hipotéticas (DFAST)

| tipo | ref CDS | ref con product | ref hyp | % ref | query CDS | query con product | query hyp | % query |
|---|---|---|---|---|---|---|---|---|
| CDS | 2767 | 2767 | 1496 | 54.1 | 2630 | 2630 | 454 | 17.3 |
| tRNA | 0 | 0 | 0 | — | 0 | 0 | 0 | — |
| rRNA | 0 | 0 | 0 | — | 0 | 0 | 0 | — |

#### Diferencias

- Falsos negativos (referencia sin pareja): **203**.
- Falsos positivos (herramienta sin pareja): **69**.
- Pares solo con hebra distinta: **1**.
- Pares emparejados con algún borde distinto: **323**.

Las listas de diferencias están truncadas a 100 elementos por lista; los conteos son completos.

> Primeras diferencias no emparejadas

| tipo | feature | seqid | coordenadas | hebra | product/locus |
|---|---|---|---|---|---|
| FN | CDS | NC_007795.1 | 33555..34034 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 34375..34509 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 34473..34565 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 37019..37120 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 42281..42382 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 43976..44179 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 93378..93500 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 113507..113710 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 114927..115214 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 166583..166684 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 213824..213931 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 227920..228027 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 229502..229645 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 232547..232645 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 259676..259810 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 263531..263710 | + | NADH dehydrogenase subunit |
| FN | CDS | NC_007795.1 | 302336..302476 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 318835..318927 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 336647..337054 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 363376..363534 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 366986..367195 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 374638..374739 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 391435..391572 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 399186..399287 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 412193..412321 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 418467..418883 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 444428..444529 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 448172..448270 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 469421..469552 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 472097..472207 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 489935..490027 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 542832..543005 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 547608..547865 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 584859..584960 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 589698..589814 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 639816..639908 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 666157..666273 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 668104..668199 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 707687..707893 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 721050..721148 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 722569..722670 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 730469..730573 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 730850..731005 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 752513..752617 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 754050..754157 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 789868..790014 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 791370..791525 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 811872..811985 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 823343..823480 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 826353..826451 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 833949..834059 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 860487..860678 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 942257..942364 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 942903..943004 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 952891..952995 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 963913..964023 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 976868..976990 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 993974..994093 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 1003967..1004122 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 1006242..1006364 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 1025649..1025816 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 1035362..1035565 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 1048947..1049240 | + | iron compound ABC transporter permease |
| FN | CDS | NC_007795.1 | 1050415..1050522 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 1062304..1062399 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 1070176..1070343 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 1071183..1071275 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 1074500..1074745 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 1075291..1075509 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 1151425..1151523 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 1180529..1180621 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 1181622..1181735 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 1244641..1244838 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 1246952..1247062 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 1248246..1248419 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 1249705..1249827 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 1268125..1268229 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 1273637..1273741 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 1279048..1279146 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 1296727..1296825 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 1299635..1299802 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 1325145..1325330 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 1331232..1331339 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 1362554..1362685 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 1373947..1374168 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 1445100..1445414 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 1450026..1450217 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 1450171..1450308 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 1450334..1450450 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 1453447..1453563 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 1458909..1459028 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 1490813..1490917 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 1503842..1504015 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 1528896..1529150 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 1532809..1532901 | + | hypothetical protein |
| FN | CDS | NC_007795.1 | 1556581..1556727 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 1562274..1562477 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 1597626..1597754 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 1638847..1638975 | - | hypothetical protein |
| FN | CDS | NC_007795.1 | 1661240..1661359 | - | hypothetical protein |
| FP | CDS | NC_007795.1 | 113771..113869 | - | hypothetical protein |
| FP | CDS | NC_007795.1 | 114553..114843 | + | hypothetical protein |
| FP | CDS | NC_007795.1 | 115867..116088 | + | hypothetical protein |
| FP | CDS | NC_007795.1 | 264380..264784 | + | hypothetical protein |
| FP | CDS | NC_007795.1 | 591277..591744 | + | hypothetical protein |
| FP | CDS | NC_007795.1 | 613988..614191 | - | DUF1659 domain-containing protein |
| FP | CDS | NC_007795.1 | 620050..620352 | + | Na+/H+ antiporter Mnh2 subunit F |
| FP | CDS | NC_007795.1 | 796106..797515 | + | von Willebrand factor binding protein Vwb |
| FP | CDS | NC_007795.1 | 889974..890408 | + | MAP domain-containing protein |
| FP | CDS | NC_007795.1 | 947205..947912 | + | hypothetical protein |
| FP | CDS | NC_007795.1 | 947905..948687 | + | hypothetical protein |
| FP | CDS | NC_007795.1 | 1041899..1042030 | + | 50S ribosomal protein L32 |
| FP | CDS | NC_007795.1 | 1111031..1111225 | - | hypothetical protein |
| FP | CDS | NC_007795.1 | 1180820..1180921 | + | hypothetical protein |
| FP | CDS | NC_007795.1 | 1180899..1181588 | + | 30S ribosomal protein S2 |
| FP | CDS | NC_007795.1 | 1245423..1245512 | + | hypothetical protein |
| FP | CDS | NC_007795.1 | 1250775..1251065 | + | hypothetical protein |
| FP | CDS | NC_007795.1 | 1280956..1283799 | + | SMC family ATPase |
| FP | CDS | NC_007795.1 | 1349111..1349503 | - | hypothetical protein |
| FP | CDS | NC_007795.1 | 1418908..1419024 | - | hypothetical protein |
| FP | CDS | NC_007795.1 | 1495120..1495428 | - | hypothetical protein |
| FP | CDS | NC_007795.1 | 1562458..1562874 | - | hypothetical protein |
| FP | CDS | NC_007795.1 | 1562911..1563921 | - | hypothetical protein |
| FP | CDS | NC_007795.1 | 1631841..1632614 | - | tRNA 2-thiouridine(34) synthase MnmA |
| FP | CDS | NC_007795.1 | 1632623..1632958 | - | hypothetical protein |
| FP | CDS | NC_007795.1 | 1636131..1636235 | - | hypothetical protein |
| FP | CDS | NC_007795.1 | 1636204..1636554 | - | Rrf2 family transcriptional regulator |
| FP | CDS | NC_007795.1 | 1705674..1706480 | - | hypothetical protein |
| FP | CDS | NC_007795.1 | 1706498..1707325 | - | hypothetical protein |
| FP | CDS | NC_007795.1 | 1795368..1795592 | - | hypothetical protein |
| FP | CDS | NC_007795.1 | 1795589..1796080 | - | hypothetical protein |
| FP | CDS | NC_007795.1 | 1813964..1814788 | + | hypothetical protein |
| FP | CDS | NC_007795.1 | 1814821..1815429 | + | hypothetical protein |
| FP | CDS | NC_007795.1 | 1897329..1898531 | + | ISL3-like element IS1181 family transposase |
| FP | CDS | NC_007795.1 | 1960237..1960566 | - | hypothetical protein |
| FP | CDS | NC_007795.1 | 2031738..2031938 | + | hypothetical protein |
| FP | CDS | NC_007795.1 | 2059099..2059299 | - | DUF1514 family protein |
| FP | CDS | NC_007795.1 | 2065119..2065268 | - | hypothetical protein |
| FP | CDS | NC_007795.1 | 2069455..2069616 | - | hypothetical protein |
| FP | CDS | NC_007795.1 | 2069848..2070204 | - | hypothetical protein |
| FP | CDS | NC_007795.1 | 2074631..2075455 | + | sphingomyelin phosphodiesterase |
| FP | CDS | NC_007795.1 | 2133600..2134388 | - | PP2C family protein-serine/threonine phosphatase |
| FP | CDS | NC_007795.1 | 2134363..2134590 | - | hypothetical protein |
| FP | CDS | NC_007795.1 | 2212638..2212841 | - | hypothetical protein |
| FP | CDS | NC_007795.1 | 2237284..2238486 | + | ISL3-like element IS1181 family transposase |
| FP | CDS | NC_007795.1 | 2244932..2245879 | - | Mrp/NBP35 family ATP-binding protein |
| FP | CDS | NC_007795.1 | 2264437..2265237 | + | hypothetical protein |
| FP | CDS | NC_007795.1 | 2265255..2265554 | + | hypothetical protein |
| FP | CDS | NC_007795.1 | 2265717..2266061 | + | hypothetical protein |
| FP | CDS | NC_007795.1 | 2289534..2289959 | + | MAP domain-containing protein |
| FP | CDS | NC_007795.1 | 2349715..2349897 | - | hypothetical protein |
| FP | CDS | NC_007795.1 | 2350014..2350880 | - | acyl-CoA dehydrogenase family protein |
| FP | CDS | NC_007795.1 | 2390182..2391495 | - | Na+/H+ antiporter NhaC family protein |
| FP | CDS | NC_007795.1 | 2420344..2420622 | - | hypothetical protein |
| FP | CDS | NC_007795.1 | 2420613..2421551 | - | multidrug effflux MFS transporter |
| FP | CDS | NC_007795.1 | 2423442..2423579 | - | hypothetical protein |
| FP | CDS | NC_007795.1 | 2423661..2423897 | - | hypothetical protein |
| FP | CDS | NC_007795.1 | 2425635..2425991 | + | hypothetical protein |
| FP | CDS | NC_007795.1 | 2457007..2457645 | + | hypothetical protein |
| FP | CDS | NC_007795.1 | 2457715..2457921 | + | hypothetical protein |
| FP | CDS | NC_007795.1 | 2558843..2559955 | - | membrane protein |
| FP | CDS | NC_007795.1 | 2591200..2591868 | - | hypothetical protein |
| FP | CDS | NC_007795.1 | 2591843..2593024 | - | hypothetical protein |
| FP | CDS | NC_007795.1 | 2671337..2671855 | - | hypothetical protein |
| FP | CDS | NC_007795.1 | 2671852..2672676 | - | hypothetical protein |
| FP | CDS | NC_007795.1 | 2721168..2721272 | - | SE2200 family small protein |
| FP | tRNA | NC_007795.1 | 1899027..1899101 | - | tRNA-Gly |
| FP | tRNA | NC_007795.1 | 1899115..1899189 | - | tRNA-Gly |
| FP | tRNA | NC_007795.1 | 2340343..2340435 | - | tRNA-Gln |

> Primeros pares con coordenadas no exactas

| tipo | seqid | ref | Δstart | Δend | solape |
|---|---|---|---|---|---|
| CDS | NC_007795.1 | 3670..3915 | +9 | +0 | 0.9634 |
| CDS | NC_007795.1 | 9755..10456 | +0 | +111 | 0.8635 |
| CDS | NC_007795.1 | 17365..18303 | +15 | +0 | 0.9840 |
| CDS | NC_007795.1 | 27515..28798 | -51 | +0 | 0.9618 |
| CDS | NC_007795.1 | 42681..43613 | +0 | +54 | 0.9453 |
| CDS | NC_007795.1 | 48313..49056 | +9 | +0 | 0.9879 |
| CDS | NC_007795.1 | 56284..57054 | -33 | +0 | 0.9590 |
| CDS | NC_007795.1 | 58855..61053 | -39 | +0 | 0.9826 |
| CDS | NC_007795.1 | 64265..65926 | +0 | +6 | 0.9964 |
| CDS | NC_007795.1 | 76521..77435 | +0 | +84 | 0.9159 |
| CDS | NC_007795.1 | 79885..80730 | -135 | +0 | 0.8624 |
| CDS | NC_007795.1 | 83505..84761 | +84 | +0 | 0.9332 |
| CDS | NC_007795.1 | 86396..88246 | +72 | +0 | 0.9611 |
| CDS | NC_007795.1 | 88221..88997 | -3 | +0 | 0.9962 |
| CDS | NC_007795.1 | 93570..94550 | +9 | +0 | 0.9908 |
| CDS | NC_007795.1 | 101370..102125 | +0 | -15 | 0.9802 |
| CDS | NC_007795.1 | 103075..104427 | -6 | +0 | 0.9956 |
| CDS | NC_007795.1 | 144157..145188 | -3 | +0 | 0.9971 |
| CDS | NC_007795.1 | 157205..157684 | +0 | +15 | 0.9697 |
| CDS | NC_007795.1 | 159987..161018 | +0 | +3 | 0.9971 |
| CDS | NC_007795.1 | 164125..164622 | +0 | +60 | 0.8925 |
| CDS | NC_007795.1 | 166663..168660 | +0 | +48 | 0.9765 |
| CDS | NC_007795.1 | 169309..170349 | -15 | +0 | 0.9858 |
| CDS | NC_007795.1 | 172714..173589 | -3 | +0 | 0.9966 |
| CDS | NC_007795.1 | 181951..183429 | +45 | +0 | 0.9696 |
| CDS | NC_007795.1 | 199583..200155 | +0 | -78 | 0.8639 |
| CDS | NC_007795.1 | 220761..222338 | +0 | -15 | 0.9905 |
| CDS | NC_007795.1 | 222779..223921 | -87 | +0 | 0.9293 |
| CDS | NC_007795.1 | 224235..225647 | +0 | +63 | 0.9573 |
| CDS | NC_007795.1 | 228482..229291 | -144 | +0 | 0.8491 |
| CDS | NC_007795.1 | 239969..240643 | -42 | +0 | 0.9414 |
| CDS | NC_007795.1 | 241836..243377 | -153 | +0 | 0.9097 |
| CDS | NC_007795.1 | 251527..253266 | -15 | +0 | 0.9915 |
| CDS | NC_007795.1 | 254122..254565 | +6 | +0 | 0.9865 |
| CDS | NC_007795.1 | 274790..275683 | +0 | +9 | 0.9900 |
| CDS | NC_007795.1 | 279728..280009 | +39 | +0 | 0.8617 |
| CDS | NC_007795.1 | 286569..287243 | -3 | +0 | 0.9956 |
| CDS | NC_007795.1 | 287243..287560 | +6 | +0 | 0.9811 |
| CDS | NC_007795.1 | 287570..289414 | -9 | +0 | 0.9951 |
| CDS | NC_007795.1 | 297408..298685 | +0 | +30 | 0.9771 |
| CDS | NC_007795.1 | 302466..303416 | -81 | +0 | 0.9215 |
| CDS | NC_007795.1 | 312595..313908 | +0 | +48 | 0.9648 |
| CDS | NC_007795.1 | 339971..341653 | -30 | +0 | 0.9825 |
| CDS | NC_007795.1 | 341746..342402 | +0 | -51 | 0.9224 |
| CDS | NC_007795.1 | 342464..342688 | +0 | -9 | 0.9600 |
| CDS | NC_007795.1 | 362010..362948 | +0 | +36 | 0.9631 |
| CDS | NC_007795.1 | 362987..363298 | +0 | -21 | 0.9327 |
| CDS | NC_007795.1 | 368020..368403 | +0 | +3 | 0.9922 |
| CDS | NC_007795.1 | 369357..369650 | +0 | -36 | 0.8776 |
| CDS | NC_007795.1 | 372329..373084 | -51 | +0 | 0.9368 |
| CDS | NC_007795.1 | 399449..400921 | -84 | +0 | 0.9461 |
| CDS | NC_007795.1 | 407345..408082 | -78 | +0 | 0.9044 |
| CDS | NC_007795.1 | 410343..410603 | +0 | -66 | 0.7471 |
| CDS | NC_007795.1 | 431034..432152 | +0 | -6 | 0.9946 |
| CDS | NC_007795.1 | 445284..446957 | -24 | +0 | 0.9859 |
| CDS | NC_007795.1 | 457344..458183 | -87 | +0 | 0.9061 |
| CDS | NC_007795.1 | 503169..503630 | -9 | +0 | 0.9809 |
| CDS | NC_007795.1 | 503727..504215 | -78 | +0 | 0.8624 |
| CDS | NC_007795.1 | 504205..505212 | -3 | +0 | 0.9970 |
| CDS | NC_007795.1 | 509574..510629 | -18 | +0 | 0.9832 |
| CDS | NC_007795.1 | 513063..513710 | -33 | +0 | 0.9515 |
| CDS | NC_007795.1 | 525875..529471 | -27 | +0 | 0.9925 |
| CDS | NC_007795.1 | 574840..575490 | +0 | +114 | 0.8510 |
| CDS | NC_007795.1 | 578134..578898 | -72 | +0 | 0.9140 |
| CDS | NC_007795.1 | 584361..584777 | +0 | +24 | 0.9456 |
| CDS | NC_007795.1 | 586908..587465 | -69 | +0 | 0.8900 |
| CDS | NC_007795.1 | 588392..588937 | -72 | +0 | 0.8835 |
| CDS | NC_007795.1 | 590589..591230 | +9 | +0 | 0.9860 |
| CDS | NC_007795.1 | 593183..593752 | -72 | +0 | 0.8879 |
| CDS | NC_007795.1 | 594375..594635 | +15 | +0 | 0.9425 |
| CDS | NC_007795.1 | 599248..599769 | +15 | +0 | 0.9713 |
| CDS | NC_007795.1 | 601586..601969 | -33 | +0 | 0.9209 |
| CDS | NC_007795.1 | 605090..605917 | -60 | +0 | 0.9324 |
| CDS | NC_007795.1 | 606065..607015 | -36 | +0 | 0.9635 |
| CDS | NC_007795.1 | 612623..613516 | +0 | +36 | 0.9613 |
| CDS | NC_007795.1 | 620306..620764 | +21 | +0 | 0.9542 |
| CDS | NC_007795.1 | 621143..623146 | -39 | +0 | 0.9809 |
| CDS | NC_007795.1 | 685710..686057 | +0 | -63 | 0.8190 |
| CDS | NC_007795.1 | 712364..714247 | +6 | +0 | 0.9968 |
| CDS | NC_007795.1 | 725114..725848 | +0 | +132 | 0.8478 |
| CDS | NC_007795.1 | 727082..729247 | +60 | +0 | 0.9723 |
| CDS | NC_007795.1 | 731044..732015 | -42 | +0 | 0.9586 |
| CDS | NC_007795.1 | 747288..748226 | -144 | +0 | 0.8670 |
| CDS | NC_007795.1 | 753021..753896 | -117 | +0 | 0.8822 |
| CDS | NC_007795.1 | 757511..758527 | +36 | +0 | 0.9646 |
| CDS | NC_007795.1 | 773243..774145 | +0 | -33 | 0.9635 |
| CDS | NC_007795.1 | 804682..805242 | -21 | +0 | 0.9639 |
| CDS | NC_007795.1 | 812195..812581 | -3 | +0 | 0.9923 |
| CDS | NC_007795.1 | 819904..821154 | +9 | +0 | 0.9928 |
| CDS | NC_007795.1 | 832811..833695 | -75 | +0 | 0.9219 |
| CDS | NC_007795.1 | 850790..851119 | +0 | +12 | 0.9649 |
| CDS | NC_007795.1 | 855577..855921 | -33 | +0 | 0.9127 |
| CDS | NC_007795.1 | 867744..868211 | -108 | +0 | 0.8125 |
| CDS | NC_007795.1 | 879806..880072 | -42 | +0 | 0.8641 |
| CDS | NC_007795.1 | 912569..912886 | +0 | +48 | 0.8689 |
| CDS | NC_007795.1 | 925707..926915 | +0 | -18 | 0.9851 |
| CDS | NC_007795.1 | 929948..930208 | +9 | +0 | 0.9655 |
| CDS | NC_007795.1 | 933094..935418 | +15 | +0 | 0.9935 |
| CDS | NC_007795.1 | 936932..938443 | -3 | +0 | 0.9980 |
| CDS | NC_007795.1 | 938930..939460 | +0 | +39 | 0.9316 |


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
/opt/conda/bin/python benchmarks/external/external_benchmark.py --genome data/reference_genomes/saureus/saureus_genomic.fna.gz --reference data/reference_genomes/saureus/saureus_genomic.gff.gz --tools bakta,prokka,dfast --run --threads 1 --json reports/external_benchmark/benchmark.json --md reports/external_benchmark/benchmark.md --bakta-db /data/bakta_db
```

Para comparar salidas producidas en otro host, usar `--query-gff TOOL=RUTA` sin `--run`; el origen queda como `precomputed`.
