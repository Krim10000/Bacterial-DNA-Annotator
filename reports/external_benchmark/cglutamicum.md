# Benchmark externo CPU: Bakta, Prokka y DFAST

**Generado:** 2026-11-25T21:03:00Z

**Esquema:** `external-benchmark-run/v1`

**Alcance:** Comparación local de CDS/tRNA/rRNA de Bakta, Prokka y DFAST contra una referencia GFF3/GBFF. CPU, sin descargas de datos.

**Comando reproducible:**

```bash
/opt/conda/bin/python /opt/benchmark/external/external_benchmark.py --genome /work/data/reference_genomes/cglutamicum/cglutamicum_genomic.fna.gz --reference /work/data/reference_genomes/cglutamicum/cglutamicum_genomic.gff.gz --tools all --run --threads 1 --bakta-db /data/bakta_db --extra-arg dfast=--dbroot --extra-arg dfast=/data/dfast_db --save-gffs /work/reports/external_benchmark/multi_v2/cglutamicum/evidence --json /work/reports/external_benchmark/multi_v2/cglutamicum.json --md /work/reports/external_benchmark/multi_v2/cglutamicum.md
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

- FASTA: `data/reference_genomes/cglutamicum/cglutamicum_genomic.fna.gz` (SHA256 `8aff3c6bb0d77233…`, 1 contig, 3309401 bp).
- Referencia: `data/reference_genomes/cglutamicum/cglutamicum_genomic.gff.gz` (formato `gff3`, SHA256 `b19bed9cb77d6291…`).
- Conteos de referencia: CDS **2948**, tRNA **59**, rRNA **18**.
- CDS pseudo excluidos en referencia: **51**.
- SHA256 FASTA: `8aff3c6bb0d772337e40f1790722e5c71eb92dd3be2d24faf1b959c6b3bca9ec`.
- SHA256 referencia: `b19bed9cb77d6291fe97c0a409ded73d0a80abf024ebc49f539794cd47cf43f9`.

## Validación

- Inicio: `2026-11-25T19:46:06Z`; fin: `2026-11-25T21:03:00Z`; total: **4613.85 s**.
- JSON declarado: `reports/external_benchmark/multi_v2/cglutamicum.json`.
- Markdown declarado: `reports/external_benchmark/multi_v2/cglutamicum.md`.

| check | estado | detalle |
|---|---|---|
| genome_path | ok | FASTA existe y es legible |
| genome_sha256 | ok | SHA256 8aff3c6bb0d772337e40f1790722e5c71eb92dd3be2d24faf1b959c6b3bca9ec |
| reference_path | ok | referencia gff3 con 3025 features |
| reference_sha256 | ok | SHA256 b19bed9cb77d6291fe97c0a409ded73d0a80abf024ebc49f539794cd47cf43f9 |
| reference_features | ok | CDS=2948 tRNA=59 rRNA=18 pseudo_excluidos=51 |
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
| Bakta | executed | ok | 0.9656 | 0.9916 | 1.0000 | 0.9481 | 0.9851 | 0.9663 | 2052.2 |
| Prokka | executed | ok | 0.9656 | 0.9916 | 1.0000 | 0.9481 | 0.9851 | 0.9663 | 1779.5 |
| DFAST | executed | ok | 0.9656 | 0.9916 | 1.0000 | 0.9481 | 0.9851 | 0.9663 | 740.3 |

La precisión/recall/F1 se calculan contra la referencia local. Si una herramienta no está instalada, su estado es `not_available` y no hay métricas que mostrar.

## Bakta

- Estado: **ok**.
- Fuente: `executed`.
- Bakta DB: `ok` (válida: **sí**, ruta `/data/bakta_db`, versión `0`).
- Versión: `bakta 1.12.2` (estado `ok`).
- Ejecución: **ok**, 2052.18 s, RSS máx. 2653.2 MiB (`2026-11-25T19:46:12Z` → `2026-11-25T20:20:25Z`).
- Salida: `/tmp/external_benchmark_bakta_1u85bdfq/out/benchmark.gff3` (3995064 bytes, SHA256 `8bdb5e84d919e587…`).
- Comando: `/usr/local/bin/bakta --output /tmp/external_benchmark_bakta_1u85bdfq/out --prefix benchmark --threads 1 --force --skip-plot --db /data/bakta_db /tmp/benchmark_genome_c12m6ik_.fna`.
- Anotación comparada: `reports/external_benchmark/multi_v2/cglutamicum/evidence/bakta/benchmark.gff3` (formato `gff3`, SHA256 `8bdb5e84d919e587…`).

### Features (Bakta)

| tipo | ref | query | TP | FP | FN | precisión | recall | F1 | exactos | borde distinto | solape medio |
|---|---|---|---|---|---|---|---|---|---|---|---|
| CDS | 2948 | 3065 | 2903 | 162 | 45 | 0.9471 | 0.9847 | 0.9656 | 2503 | 400 | 0.9909 |
| tRNA | 59 | 60 | 59 | 1 | 0 | 0.9833 | 1.0000 | 0.9916 | 59 | 0 | 1.0000 |
| rRNA | 18 | 18 | 18 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 18 | 0 | 1.0000 |

### Atributos (Bakta)

| tipo | pares | product en ambos | product literal | product tokens | gene literal | locus literal |
|---|---|---|---|---|---|---|
| CDS | 2903 | 2903 | 0.2949 | 0.4697 | 0.9295 | 0.0000 |
| tRNA | 59 | 59 | 0.0000 | 0.6455 | — | 0.0000 |
| rRNA | 18 | 18 | 1.0000 | 1.0000 | 1.0000 | 0.0000 |

### Hipotéticas (Bakta)

| tipo | ref CDS | ref con product | ref hyp | % ref | query CDS | query con product | query hyp | % query |
|---|---|---|---|---|---|---|---|---|
| CDS | 2948 | 2948 | 625 | 21.2 | 3065 | 3065 | 905 | 29.5 |
| tRNA | 0 | 0 | 0 | — | 0 | 0 | 0 | — |
| rRNA | 0 | 0 | 0 | — | 0 | 0 | 0 | — |

#### Diferencias

- Falsos negativos (referencia sin pareja): **45**.
- Falsos positivos (herramienta sin pareja): **163**.
- Pares solo con hebra distinta: **5**.
- Pares emparejados con algún borde distinto: **400**.

Las listas de diferencias están truncadas a 100 elementos por lista; los conteos son completos.

> Primeras diferencias no emparejadas

| tipo | feature | seqid | coordenadas | hebra | product/locus |
|---|---|---|---|---|---|
| FN | CDS | NC_003450.3 | 23612..23779 | - | hypothetical protein |
| FN | CDS | NC_003450.3 | 75745..76038 | + | CHY zinc finger protein |
| FN | CDS | NC_003450.3 | 221929..222210 | + | 2Fe-2S iron-sulfur cluster-binding protein |
| FN | CDS | NC_003450.3 | 377683..377835 | + | transposase |
| FN | CDS | NC_003450.3 | 378296..378514 | + | transposase |
| FN | CDS | NC_003450.3 | 674796..675077 | - | putative quinol monooxygenase |
| FN | CDS | NC_003450.3 | 805346..805420 | - | 50S ribosomal protein bL37 |
| FN | CDS | NC_003450.3 | 932305..932490 | + | hypothetical protein |
| FN | CDS | NC_003450.3 | 951556..951654 | + | hypothetical protein |
| FN | CDS | NC_003450.3 | 994650..994841 | + | hypothetical protein |
| FN | CDS | NC_003450.3 | 1015982..1016563 | + | hypothetical protein |
| FN | CDS | NC_003450.3 | 1103951..1104106 | + | hypothetical protein |
| FN | CDS | NC_003450.3 | 1122872..1123030 | + | hypothetical protein |
| FN | CDS | NC_003450.3 | 1280967..1281254 | + | hypothetical protein |
| FN | CDS | NC_003450.3 | 1367167..1367289 | + | hypothetical protein |
| FN | CDS | NC_003450.3 | 1367885..1368070 | - | DoxX family protein |
| FN | CDS | NC_003450.3 | 1371976..1372326 | - | Abi family protein |
| FN | CDS | NC_003450.3 | 1402272..1402706 | + | MepB family protein |
| FN | CDS | NC_003450.3 | 1555835..1556182 | - | hypothetical protein |
| FN | CDS | NC_003450.3 | 1564237..1564485 | + | Fic family protein |
| FN | CDS | NC_003450.3 | 1780504..1780830 | - | alpha/beta hydrolase |
| FN | CDS | NC_003450.3 | 1797182..1797340 | + | hypothetical protein |
| FN | CDS | NC_003450.3 | 1888228..1888599 | - | hypothetical protein |
| FN | CDS | NC_003450.3 | 2098409..2098945 | - | ribosome maturation factor RimP |
| FN | CDS | NC_003450.3 | 2357328..2357711 | + | hypothetical protein |
| FN | CDS | NC_003450.3 | 2385911..2386192 | - | hypothetical protein |
| FN | CDS | NC_003450.3 | 2386210..2386584 | + | hypothetical protein |
| FN | CDS | NC_003450.3 | 2393165..2393407 | + | hypothetical protein |
| FN | CDS | NC_003450.3 | 2609336..2609479 | - | hypothetical protein |
| FN | CDS | NC_003450.3 | 2612456..2612593 | - | hypothetical protein |
| FN | CDS | NC_003450.3 | 2612628..2612897 | - | hypothetical protein |
| FN | CDS | NC_003450.3 | 2627777..2628088 | + | LacI family DNA-binding transcriptional regulator |
| FN | CDS | NC_003450.3 | 2690151..2690441 | + | fluoride efflux transporter family protein |
| FN | CDS | NC_003450.3 | 2840756..2840929 | - | DUF4236 domain-containing protein |
| FN | CDS | NC_003450.3 | 2997479..2997688 | - | hypothetical protein |
| FN | CDS | NC_003450.3 | 3030360..3030539 | + | ABC transporter ATP-binding protein |
| FN | CDS | NC_003450.3 | 3144493..3144633 | - | hypothetical protein |
| FN | CDS | NC_003450.3 | 3160629..3160802 | - | helix-turn-helix domain-containing protein |
| FN | CDS | NC_003450.3 | 3160721..3160939 | - | recombinase family protein |
| FN | CDS | NC_003450.3 | 3161680..3162015 | - | cupin domain-containing protein |
| FN | CDS | NC_003450.3 | 3166265..3166438 | - | hypothetical protein |
| FN | CDS | NC_003450.3 | 3172908..3173090 | - | hypothetical protein |
| FN | CDS | NC_003450.3 | 3233152..3233205 | + | trp operon leader peptide |
| FN | CDS | NC_003450.3 | 3282782..3283105 | + | MGMT family protein |
| FN | CDS | NC_003450.3 | 3307689..3307988 | - | membrane protein insertion efficiency factor YidD |
| FP | CDS | NC_003450.3 | 21426..21527 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 50751..50972 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 75727..75981 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 85426..85623 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 105054..105254 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 106429..106590 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 116207..116434 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 155874..156170 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 186769..187305 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 201338..201580 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 201617..201730 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 222097..222210 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 239793..239948 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 273296..273484 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 332919..334565 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 334516..334956 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 345811..345975 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 378141..378230 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 378284..378538 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 387418..387528 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 452006..452212 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 527604..527759 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 532008..532526 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 671141..671770 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 671736..672656 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 674853..674975 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 748368..748637 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 847727..848029 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 855839..855949 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 895441..895596 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 953747..954148 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 954277..954756 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 954974..955357 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 955354..955500 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 959855..960256 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 960385..960864 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 982291..982674 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1014925..1015119 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 1079037..1079183 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 1084261..1084368 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 1088205..1088507 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1089862..1090008 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1095384..1095677 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1096185..1096331 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1099012..1099164 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1192129..1192257 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1227890..1228639 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 1228718..1229098 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 1341957..1342457 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1342458..1342727 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1342791..1343675 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1366270..1366500 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1371855..1371965 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1375930..1376088 | - | Transposase |
| FP | CDS | NC_003450.3 | 1400881..1401336 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 1402365..1402874 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1404867..1405301 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1408250..1408453 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 1409141..1409284 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 1555835..1555999 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1609839..1610192 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1639370..1639504 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1639656..1639820 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 1640153..1640536 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1652891..1653577 | - | ArsR family transcriptional regulator |
| FP | CDS | NC_003450.3 | 1663933..1664115 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1671270..1671677 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1673263..1674105 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1743843..1743971 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 1744025..1744522 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 1777332..1777649 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 1777698..1777823 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 1780498..1780671 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 1790831..1790956 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1793548..1793703 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 1798791..1799042 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 1799368..1799505 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1799624..1799743 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 1826641..1826937 | - | Transposase |
| FP | CDS | NC_003450.3 | 1829685..1829900 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1887932..1888129 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1888228..1889445 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1917561..1918070 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1928534..1928911 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 1930987..1931190 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1937202..1937414 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 1943095..1943313 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 1943345..1943656 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 1945592..1945882 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1967114..1967233 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1971401..1971550 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1987504..1987821 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1987884..1988306 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1988367..1988483 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1988527..1988664 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1994302..1994430 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1994434..1994568 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1995780..1996019 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 2013996..2014100 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 2025218..2025388 | + | hypothetical protein |

> Primeros pares con coordenadas no exactas

| tipo | seqid | ref | Δstart | Δend | solape |
|---|---|---|---|---|---|
| CDS | NC_003450.3 | 7830..8798 | +54 | +0 | 0.9443 |
| CDS | NC_003450.3 | 10104..11177 | +0 | -6 | 0.9944 |
| CDS | NC_003450.3 | 11520..11729 | +0 | +39 | 0.8434 |
| CDS | NC_003450.3 | 11891..14401 | -60 | +0 | 0.9767 |
| CDS | NC_003450.3 | 21564..22127 | +33 | +0 | 0.9415 |
| CDS | NC_003450.3 | 24295..24732 | -81 | +0 | 0.8439 |
| CDS | NC_003450.3 | 34366..34902 | -36 | +0 | 0.9372 |
| CDS | NC_003450.3 | 36244..37014 | +0 | +228 | 0.7718 |
| CDS | NC_003450.3 | 53094..54011 | +3 | +0 | 0.9967 |
| CDS | NC_003450.3 | 60699..62324 | -30 | +0 | 0.9819 |
| CDS | NC_003450.3 | 86399..87244 | -81 | +0 | 0.9126 |
| CDS | NC_003450.3 | 89615..90448 | -42 | +0 | 0.9521 |
| CDS | NC_003450.3 | 103503..104912 | -24 | +0 | 0.9833 |
| CDS | NC_003450.3 | 105170..105892 | +0 | -270 | 0.6266 |
| CDS | NC_003450.3 | 107444..110893 | -9 | +0 | 0.9974 |
| CDS | NC_003450.3 | 115578..115946 | -42 | +0 | 0.8978 |
| CDS | NC_003450.3 | 127189..127998 | +0 | +6 | 0.9926 |
| CDS | NC_003450.3 | 129486..130055 | +0 | +63 | 0.9005 |
| CDS | NC_003450.3 | 132968..134113 | +0 | +21 | 0.9820 |
| CDS | NC_003450.3 | 135515..136201 | +0 | -168 | 0.7555 |
| CDS | NC_003450.3 | 140329..141792 | +3 | +0 | 0.9980 |
| CDS | NC_003450.3 | 144725..145483 | +6 | +0 | 0.9921 |
| CDS | NC_003450.3 | 153223..155613 | +0 | +39 | 0.9840 |
| CDS | NC_003450.3 | 157614..158141 | -30 | +0 | 0.9462 |
| CDS | NC_003450.3 | 169038..169994 | -63 | +0 | 0.9382 |
| CDS | NC_003450.3 | 170044..170919 | -48 | +0 | 0.9481 |
| CDS | NC_003450.3 | 193262..194467 | -405 | +0 | 0.7486 |
| CDS | NC_003450.3 | 220151..221107 | +0 | +24 | 0.9755 |
| CDS | NC_003450.3 | 224354..225247 | +27 | +0 | 0.9698 |
| CDS | NC_003450.3 | 238218..239528 | -42 | +0 | 0.9690 |
| CDS | NC_003450.3 | 241880..242902 | +0 | +18 | 0.9827 |
| CDS | NC_003450.3 | 244466..244819 | +18 | +0 | 0.9492 |
| CDS | NC_003450.3 | 247319..248575 | -39 | +0 | 0.9699 |
| CDS | NC_003450.3 | 256599..257897 | +18 | +0 | 0.9861 |
| CDS | NC_003450.3 | 257900..258532 | +81 | +0 | 0.8720 |
| CDS | NC_003450.3 | 274366..275874 | -42 | +0 | 0.9729 |
| CDS | NC_003450.3 | 278385..278732 | +0 | +102 | 0.7733 |
| CDS | NC_003450.3 | 281401..282927 | +0 | +6 | 0.9961 |
| CDS | NC_003450.3 | 299723..301261 | +0 | -6 | 0.9961 |
| CDS | NC_003450.3 | 305876..306703 | -18 | +0 | 0.9787 |
| CDS | NC_003450.3 | 307724..307909 | +0 | -15 | 0.9194 |
| CDS | NC_003450.3 | 310132..311328 | +3 | +0 | 0.9975 |
| CDS | NC_003450.3 | 313622..314464 | +0 | +126 | 0.8700 |
| CDS | NC_003450.3 | 319333..321690 | +0 | +39 | 0.9837 |
| CDS | NC_003450.3 | 325919..326617 | -15 | +0 | 0.9790 |
| CDS | NC_003450.3 | 330573..330812 | +18 | +0 | 0.9250 |
| CDS | NC_003450.3 | 337539..338771 | +21 | +0 | 0.9830 |
| CDS | NC_003450.3 | 352693..353640 | +153 | +0 | 0.8386 |
| CDS | NC_003450.3 | 357234..359240 | +0 | +114 | 0.9463 |
| CDS | NC_003450.3 | 359759..360307 | +0 | +9 | 0.9839 |
| CDS | NC_003450.3 | 370685..371776 | -72 | +0 | 0.9381 |
| CDS | NC_003450.3 | 378665..379837 | +0 | -117 | 0.9003 |
| CDS | NC_003450.3 | 390784..392208 | +0 | -15 | 0.9895 |
| CDS | NC_003450.3 | 398278..399582 | -30 | +0 | 0.9775 |
| CDS | NC_003450.3 | 400039..400344 | +12 | +0 | 0.9608 |
| CDS | NC_003450.3 | 404505..405419 | +0 | +9 | 0.9903 |
| CDS | NC_003450.3 | 433062..433991 | +21 | +0 | 0.9774 |
| CDS | NC_003450.3 | 438516..439907 | +42 | +0 | 0.9698 |
| CDS | NC_003450.3 | 452637..454433 | -57 | +0 | 0.9693 |
| CDS | NC_003450.3 | 460020..461096 | -48 | +0 | 0.9573 |
| CDS | NC_003450.3 | 462557..463870 | +9 | +0 | 0.9932 |
| CDS | NC_003450.3 | 469370..470173 | +60 | +0 | 0.9254 |
| CDS | NC_003450.3 | 478089..478973 | +0 | -84 | 0.9051 |
| CDS | NC_003450.3 | 480621..480995 | +0 | +6 | 0.9843 |
| CDS | NC_003450.3 | 481391..482668 | +0 | +108 | 0.9221 |
| CDS | NC_003450.3 | 490580..491941 | +60 | +0 | 0.9559 |
| CDS | NC_003450.3 | 511536..512510 | -33 | +0 | 0.9673 |
| CDS | NC_003450.3 | 512913..516410 | +18 | +0 | 0.9949 |
| CDS | NC_003450.3 | 516494..520495 | -27 | +0 | 0.9933 |
| CDS | NC_003450.3 | 544754..547042 | +0 | -30 | 0.9869 |
| CDS | NC_003450.3 | 550696..551808 | +0 | +15 | 0.9867 |
| CDS | NC_003450.3 | 554449..554973 | +0 | -39 | 0.9257 |
| CDS | NC_003450.3 | 599759..599935 | -60 | +0 | 0.7468 |
| CDS | NC_003450.3 | 614800..615180 | +0 | +42 | 0.9007 |
| CDS | NC_003450.3 | 618090..618995 | +0 | +18 | 0.9805 |
| CDS | NC_003450.3 | 623800..624942 | -129 | +0 | 0.8986 |
| CDS | NC_003450.3 | 624985..625677 | -42 | +0 | 0.9429 |
| CDS | NC_003450.3 | 626067..626510 | +0 | +33 | 0.9308 |
| CDS | NC_003450.3 | 628563..630143 | -12 | +0 | 0.9925 |
| CDS | NC_003450.3 | 633076..633525 | +0 | -51 | 0.8867 |
| CDS | NC_003450.3 | 636086..638224 | +0 | +54 | 0.9754 |
| CDS | NC_003450.3 | 638314..639237 | +0 | +183 | 0.8347 |
| CDS | NC_003450.3 | 641037..642560 | -84 | +0 | 0.9478 |
| CDS | NC_003450.3 | 645173..647470 | +0 | +120 | 0.9504 |
| CDS | NC_003450.3 | 654676..655125 | -42 | +0 | 0.9146 |
| CDS | NC_003450.3 | 657202..657978 | +0 | +45 | 0.9453 |
| CDS | NC_003450.3 | 661166..662020 | -45 | +0 | 0.9500 |
| CDS | NC_003450.3 | 664123..664947 | +0 | +141 | 0.8540 |
| CDS | NC_003450.3 | 668189..670468 | -39 | +0 | 0.9832 |
| CDS | NC_003450.3 | 695074..696279 | +0 | +105 | 0.9199 |
| CDS | NC_003450.3 | 705211..708633 | -51 | +0 | 0.9853 |
| CDS | NC_003450.3 | 708839..709711 | +3 | +0 | 0.9966 |
| CDS | NC_003450.3 | 716283..716660 | +0 | -9 | 0.9762 |
| CDS | NC_003450.3 | 730436..731302 | +57 | +0 | 0.9343 |
| CDS | NC_003450.3 | 733060..734946 | +12 | +0 | 0.9936 |
| CDS | NC_003450.3 | 791476..793011 | +36 | +0 | 0.9766 |
| CDS | NC_003450.3 | 794810..795304 | -96 | +0 | 0.8376 |
| CDS | NC_003450.3 | 807943..809217 | +0 | -48 | 0.9624 |
| CDS | NC_003450.3 | 809528..810397 | +21 | +0 | 0.9759 |
| CDS | NC_003450.3 | 834927..835391 | +63 | +0 | 0.8645 |

- Aviso: 2 feature(s) GFF multi-segmento fusionados por ID; se comparan por su intervalo envolvente.

## Prokka

- Estado: **ok**.
- Fuente: `executed`.
- Versión: `prokka 1.15.6` (estado `ok`).
- Ejecución: **ok**, 1779.54 s, RSS máx. 2653.2 MiB (`2026-11-25T20:20:37Z` → `2026-11-25T20:50:16Z`).
- Salida: `/tmp/external_benchmark_prokka_rj52mipi/out/benchmark.gff` (4056346 bytes, SHA256 `134a6cacbf7455e5…`).
- Comando: `/usr/local/bin/prokka --outdir /tmp/external_benchmark_prokka_rj52mipi/out --prefix benchmark --cpus 1 --force --quiet /tmp/benchmark_genome_c12m6ik_.fna`.
- Anotación comparada: `reports/external_benchmark/multi_v2/cglutamicum/evidence/prokka/benchmark.gff` (formato `gff3`, SHA256 `134a6cacbf7455e5…`).

### Features (Prokka)

| tipo | ref | query | TP | FP | FN | precisión | recall | F1 | exactos | borde distinto | solape medio |
|---|---|---|---|---|---|---|---|---|---|---|---|
| CDS | 2948 | 3065 | 2903 | 162 | 45 | 0.9471 | 0.9847 | 0.9656 | 2503 | 400 | 0.9909 |
| tRNA | 59 | 60 | 59 | 1 | 0 | 0.9833 | 1.0000 | 0.9916 | 13 | 46 | 0.9890 |
| rRNA | 18 | 18 | 18 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 18 | 0 | 1.0000 |

### Atributos (Prokka)

| tipo | pares | product en ambos | product literal | product tokens | gene literal | locus literal |
|---|---|---|---|---|---|---|
| CDS | 2903 | 2903 | 0.2573 | 0.4205 | 0.7486 | 0.0000 |
| tRNA | 59 | 59 | 0.0000 | 0.6596 | — | 0.0000 |
| rRNA | 18 | 18 | 1.0000 | 1.0000 | — | 0.0000 |

### Hipotéticas (Prokka)

| tipo | ref CDS | ref con product | ref hyp | % ref | query CDS | query con product | query hyp | % query |
|---|---|---|---|---|---|---|---|---|
| CDS | 2948 | 2948 | 625 | 21.2 | 3065 | 3065 | 1408 | 45.9 |
| tRNA | 0 | 0 | 0 | — | 0 | 0 | 0 | — |
| rRNA | 0 | 0 | 0 | — | 0 | 0 | 0 | — |

#### Diferencias

- Falsos negativos (referencia sin pareja): **45**.
- Falsos positivos (herramienta sin pareja): **163**.
- Pares solo con hebra distinta: **5**.
- Pares emparejados con algún borde distinto: **446**.

Las listas de diferencias están truncadas a 100 elementos por lista; los conteos son completos.

> Primeras diferencias no emparejadas

| tipo | feature | seqid | coordenadas | hebra | product/locus |
|---|---|---|---|---|---|
| FN | CDS | NC_003450.3 | 23612..23779 | - | hypothetical protein |
| FN | CDS | NC_003450.3 | 75745..76038 | + | CHY zinc finger protein |
| FN | CDS | NC_003450.3 | 221929..222210 | + | 2Fe-2S iron-sulfur cluster-binding protein |
| FN | CDS | NC_003450.3 | 377683..377835 | + | transposase |
| FN | CDS | NC_003450.3 | 378296..378514 | + | transposase |
| FN | CDS | NC_003450.3 | 674796..675077 | - | putative quinol monooxygenase |
| FN | CDS | NC_003450.3 | 805346..805420 | - | 50S ribosomal protein bL37 |
| FN | CDS | NC_003450.3 | 932305..932490 | + | hypothetical protein |
| FN | CDS | NC_003450.3 | 951556..951654 | + | hypothetical protein |
| FN | CDS | NC_003450.3 | 994650..994841 | + | hypothetical protein |
| FN | CDS | NC_003450.3 | 1015982..1016563 | + | hypothetical protein |
| FN | CDS | NC_003450.3 | 1103951..1104106 | + | hypothetical protein |
| FN | CDS | NC_003450.3 | 1122872..1123030 | + | hypothetical protein |
| FN | CDS | NC_003450.3 | 1280967..1281254 | + | hypothetical protein |
| FN | CDS | NC_003450.3 | 1367167..1367289 | + | hypothetical protein |
| FN | CDS | NC_003450.3 | 1367885..1368070 | - | DoxX family protein |
| FN | CDS | NC_003450.3 | 1371976..1372326 | - | Abi family protein |
| FN | CDS | NC_003450.3 | 1402272..1402706 | + | MepB family protein |
| FN | CDS | NC_003450.3 | 1555835..1556182 | - | hypothetical protein |
| FN | CDS | NC_003450.3 | 1564237..1564485 | + | Fic family protein |
| FN | CDS | NC_003450.3 | 1780504..1780830 | - | alpha/beta hydrolase |
| FN | CDS | NC_003450.3 | 1797182..1797340 | + | hypothetical protein |
| FN | CDS | NC_003450.3 | 1888228..1888599 | - | hypothetical protein |
| FN | CDS | NC_003450.3 | 2098409..2098945 | - | ribosome maturation factor RimP |
| FN | CDS | NC_003450.3 | 2357328..2357711 | + | hypothetical protein |
| FN | CDS | NC_003450.3 | 2385911..2386192 | - | hypothetical protein |
| FN | CDS | NC_003450.3 | 2386210..2386584 | + | hypothetical protein |
| FN | CDS | NC_003450.3 | 2393165..2393407 | + | hypothetical protein |
| FN | CDS | NC_003450.3 | 2609336..2609479 | - | hypothetical protein |
| FN | CDS | NC_003450.3 | 2612456..2612593 | - | hypothetical protein |
| FN | CDS | NC_003450.3 | 2612628..2612897 | - | hypothetical protein |
| FN | CDS | NC_003450.3 | 2627777..2628088 | + | LacI family DNA-binding transcriptional regulator |
| FN | CDS | NC_003450.3 | 2690151..2690441 | + | fluoride efflux transporter family protein |
| FN | CDS | NC_003450.3 | 2840756..2840929 | - | DUF4236 domain-containing protein |
| FN | CDS | NC_003450.3 | 2997479..2997688 | - | hypothetical protein |
| FN | CDS | NC_003450.3 | 3030360..3030539 | + | ABC transporter ATP-binding protein |
| FN | CDS | NC_003450.3 | 3144493..3144633 | - | hypothetical protein |
| FN | CDS | NC_003450.3 | 3160629..3160802 | - | helix-turn-helix domain-containing protein |
| FN | CDS | NC_003450.3 | 3160721..3160939 | - | recombinase family protein |
| FN | CDS | NC_003450.3 | 3161680..3162015 | - | cupin domain-containing protein |
| FN | CDS | NC_003450.3 | 3166265..3166438 | - | hypothetical protein |
| FN | CDS | NC_003450.3 | 3172908..3173090 | - | hypothetical protein |
| FN | CDS | NC_003450.3 | 3233152..3233205 | + | trp operon leader peptide |
| FN | CDS | NC_003450.3 | 3282782..3283105 | + | MGMT family protein |
| FN | CDS | NC_003450.3 | 3307689..3307988 | - | membrane protein insertion efficiency factor YidD |
| FP | CDS | NC_003450.3 | 21426..21527 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 50751..50972 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 75727..75981 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 85426..85623 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 105054..105254 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 106429..106590 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 116207..116434 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 155874..156170 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 186769..187305 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 201338..201580 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 201617..201730 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 222097..222210 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 239793..239948 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 273296..273484 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 332919..334565 | + | Beta-hexosaminidase |
| FP | CDS | NC_003450.3 | 334516..334956 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 345811..345975 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 378141..378230 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 378284..378538 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 387418..387528 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 452006..452212 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 527604..527759 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 532008..532526 | + | Succinyl-CoA:acetate CoA-transferase |
| FP | CDS | NC_003450.3 | 671141..671770 | + | 2-methylcitrate dehydratase 2 |
| FP | CDS | NC_003450.3 | 671736..672656 | + | 2-methylcitrate dehydratase 2 |
| FP | CDS | NC_003450.3 | 674853..674975 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 748368..748637 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 847727..848029 | + | IS5 family transposase ISCef8 |
| FP | CDS | NC_003450.3 | 855839..855949 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 895441..895596 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 953747..954148 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 954277..954756 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 954974..955357 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 955354..955500 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 959855..960256 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 960385..960864 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 982291..982674 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1014925..1015119 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 1079037..1079183 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 1084261..1084368 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 1088205..1088507 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1089862..1090008 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1095384..1095677 | - | ISL3 family transposase IS1207 |
| FP | CDS | NC_003450.3 | 1096185..1096331 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1099012..1099164 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1192129..1192257 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1227890..1228639 | + | Acetyl esterase |
| FP | CDS | NC_003450.3 | 1228718..1229098 | + | Putative short-chain fatty acid transporter |
| FP | CDS | NC_003450.3 | 1341957..1342457 | - | Nitrate import ATP-binding protein NrtD |
| FP | CDS | NC_003450.3 | 1342458..1342727 | - | Bicarbonate transport ATP-binding protein CmpD |
| FP | CDS | NC_003450.3 | 1342791..1343675 | - | Bicarbonate transport ATP-binding protein CmpC |
| FP | CDS | NC_003450.3 | 1366270..1366500 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1371855..1371965 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1375930..1376088 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1400881..1401336 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 1402365..1402874 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1404867..1405301 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1408250..1408453 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 1409141..1409284 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 1555835..1555999 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1609839..1610192 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1639370..1639504 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1639656..1639820 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 1640153..1640536 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1652891..1653577 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1663933..1664115 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1671270..1671677 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1673263..1674105 | - | Monomeric sarcosine oxidase |
| FP | CDS | NC_003450.3 | 1743843..1743971 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 1744025..1744522 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 1777332..1777649 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 1777698..1777823 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 1780498..1780671 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 1790831..1790956 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1793548..1793703 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 1798791..1799042 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 1799368..1799505 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1799624..1799743 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 1826641..1826937 | - | IS3 family transposase IS1206 |
| FP | CDS | NC_003450.3 | 1829685..1829900 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1887932..1888129 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1888228..1889445 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1917561..1918070 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1928534..1928911 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 1930987..1931190 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1937202..1937414 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 1943095..1943313 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 1943345..1943656 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 1945592..1945882 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1967114..1967233 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1971401..1971550 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1987504..1987821 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1987884..1988306 | - | IS3 family transposase IS1206 |
| FP | CDS | NC_003450.3 | 1988367..1988483 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1988527..1988664 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1994302..1994430 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1994434..1994568 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1995780..1996019 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 2013996..2014100 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 2025218..2025388 | + | hypothetical protein |

> Primeros pares con coordenadas no exactas

| tipo | seqid | ref | Δstart | Δend | solape |
|---|---|---|---|---|---|
| CDS | NC_003450.3 | 7830..8798 | +54 | +0 | 0.9443 |
| CDS | NC_003450.3 | 10104..11177 | +0 | -6 | 0.9944 |
| CDS | NC_003450.3 | 11520..11729 | +0 | +39 | 0.8434 |
| CDS | NC_003450.3 | 11891..14401 | -60 | +0 | 0.9767 |
| CDS | NC_003450.3 | 21564..22127 | +33 | +0 | 0.9415 |
| CDS | NC_003450.3 | 24295..24732 | -81 | +0 | 0.8439 |
| CDS | NC_003450.3 | 34366..34902 | -36 | +0 | 0.9372 |
| CDS | NC_003450.3 | 36244..37014 | +0 | +228 | 0.7718 |
| CDS | NC_003450.3 | 53094..54011 | +3 | +0 | 0.9967 |
| CDS | NC_003450.3 | 60699..62324 | -30 | +0 | 0.9819 |
| CDS | NC_003450.3 | 86399..87244 | -81 | +0 | 0.9126 |
| CDS | NC_003450.3 | 89615..90448 | -42 | +0 | 0.9521 |
| CDS | NC_003450.3 | 103503..104912 | -24 | +0 | 0.9833 |
| CDS | NC_003450.3 | 105170..105892 | +0 | -270 | 0.6266 |
| CDS | NC_003450.3 | 107444..110893 | -9 | +0 | 0.9974 |
| CDS | NC_003450.3 | 115578..115946 | -42 | +0 | 0.8978 |
| CDS | NC_003450.3 | 127189..127998 | +0 | +6 | 0.9926 |
| CDS | NC_003450.3 | 129486..130055 | +0 | +63 | 0.9005 |
| CDS | NC_003450.3 | 132968..134113 | +0 | +21 | 0.9820 |
| CDS | NC_003450.3 | 135515..136201 | +0 | -168 | 0.7555 |
| CDS | NC_003450.3 | 140329..141792 | +3 | +0 | 0.9980 |
| CDS | NC_003450.3 | 144725..145483 | +6 | +0 | 0.9921 |
| CDS | NC_003450.3 | 153223..155613 | +0 | +39 | 0.9840 |
| CDS | NC_003450.3 | 157614..158141 | -30 | +0 | 0.9462 |
| CDS | NC_003450.3 | 169038..169994 | -63 | +0 | 0.9382 |
| CDS | NC_003450.3 | 170044..170919 | -48 | +0 | 0.9481 |
| CDS | NC_003450.3 | 193262..194467 | -405 | +0 | 0.7486 |
| CDS | NC_003450.3 | 220151..221107 | +0 | +24 | 0.9755 |
| CDS | NC_003450.3 | 224354..225247 | +27 | +0 | 0.9698 |
| CDS | NC_003450.3 | 238218..239528 | -42 | +0 | 0.9690 |
| CDS | NC_003450.3 | 241880..242902 | +0 | +18 | 0.9827 |
| CDS | NC_003450.3 | 244466..244819 | +18 | +0 | 0.9492 |
| CDS | NC_003450.3 | 247319..248575 | -39 | +0 | 0.9699 |
| CDS | NC_003450.3 | 256599..257897 | +18 | +0 | 0.9861 |
| CDS | NC_003450.3 | 257900..258532 | +81 | +0 | 0.8720 |
| CDS | NC_003450.3 | 274366..275874 | -42 | +0 | 0.9729 |
| CDS | NC_003450.3 | 278385..278732 | +0 | +102 | 0.7733 |
| CDS | NC_003450.3 | 281401..282927 | +0 | +6 | 0.9961 |
| CDS | NC_003450.3 | 299723..301261 | +0 | -6 | 0.9961 |
| CDS | NC_003450.3 | 305876..306703 | -18 | +0 | 0.9787 |
| CDS | NC_003450.3 | 307724..307909 | +0 | -15 | 0.9194 |
| CDS | NC_003450.3 | 310132..311328 | +3 | +0 | 0.9975 |
| CDS | NC_003450.3 | 313622..314464 | +0 | +126 | 0.8700 |
| CDS | NC_003450.3 | 319333..321690 | +0 | +39 | 0.9837 |
| CDS | NC_003450.3 | 325919..326617 | -15 | +0 | 0.9790 |
| CDS | NC_003450.3 | 330573..330812 | +18 | +0 | 0.9250 |
| CDS | NC_003450.3 | 337539..338771 | +21 | +0 | 0.9830 |
| CDS | NC_003450.3 | 352693..353640 | +153 | +0 | 0.8386 |
| CDS | NC_003450.3 | 357234..359240 | +0 | +114 | 0.9463 |
| CDS | NC_003450.3 | 359759..360307 | +0 | +9 | 0.9839 |
| CDS | NC_003450.3 | 370685..371776 | -72 | +0 | 0.9381 |
| CDS | NC_003450.3 | 378665..379837 | +0 | -117 | 0.9003 |
| CDS | NC_003450.3 | 390784..392208 | +0 | -15 | 0.9895 |
| CDS | NC_003450.3 | 398278..399582 | -30 | +0 | 0.9775 |
| CDS | NC_003450.3 | 400039..400344 | +12 | +0 | 0.9608 |
| CDS | NC_003450.3 | 404505..405419 | +0 | +9 | 0.9903 |
| CDS | NC_003450.3 | 433062..433991 | +21 | +0 | 0.9774 |
| CDS | NC_003450.3 | 438516..439907 | +42 | +0 | 0.9698 |
| CDS | NC_003450.3 | 452637..454433 | -57 | +0 | 0.9693 |
| CDS | NC_003450.3 | 460020..461096 | -48 | +0 | 0.9573 |
| CDS | NC_003450.3 | 462557..463870 | +9 | +0 | 0.9932 |
| CDS | NC_003450.3 | 469370..470173 | +60 | +0 | 0.9254 |
| CDS | NC_003450.3 | 478089..478973 | +0 | -84 | 0.9051 |
| CDS | NC_003450.3 | 480621..480995 | +0 | +6 | 0.9843 |
| CDS | NC_003450.3 | 481391..482668 | +0 | +108 | 0.9221 |
| CDS | NC_003450.3 | 490580..491941 | +60 | +0 | 0.9559 |
| CDS | NC_003450.3 | 511536..512510 | -33 | +0 | 0.9673 |
| CDS | NC_003450.3 | 512913..516410 | +18 | +0 | 0.9949 |
| CDS | NC_003450.3 | 516494..520495 | -27 | +0 | 0.9933 |
| CDS | NC_003450.3 | 544754..547042 | +0 | -30 | 0.9869 |
| CDS | NC_003450.3 | 550696..551808 | +0 | +15 | 0.9867 |
| CDS | NC_003450.3 | 554449..554973 | +0 | -39 | 0.9257 |
| CDS | NC_003450.3 | 599759..599935 | -60 | +0 | 0.7468 |
| CDS | NC_003450.3 | 614800..615180 | +0 | +42 | 0.9007 |
| CDS | NC_003450.3 | 618090..618995 | +0 | +18 | 0.9805 |
| CDS | NC_003450.3 | 623800..624942 | -129 | +0 | 0.8986 |
| CDS | NC_003450.3 | 624985..625677 | -42 | +0 | 0.9429 |
| CDS | NC_003450.3 | 626067..626510 | +0 | +33 | 0.9308 |
| CDS | NC_003450.3 | 628563..630143 | -12 | +0 | 0.9925 |
| CDS | NC_003450.3 | 633076..633525 | +0 | -51 | 0.8867 |
| CDS | NC_003450.3 | 636086..638224 | +0 | +54 | 0.9754 |
| CDS | NC_003450.3 | 638314..639237 | +0 | +183 | 0.8347 |
| CDS | NC_003450.3 | 641037..642560 | -84 | +0 | 0.9478 |
| CDS | NC_003450.3 | 645173..647470 | +0 | +120 | 0.9504 |
| CDS | NC_003450.3 | 654676..655125 | -42 | +0 | 0.9146 |
| CDS | NC_003450.3 | 657202..657978 | +0 | +45 | 0.9453 |
| CDS | NC_003450.3 | 661166..662020 | -45 | +0 | 0.9500 |
| CDS | NC_003450.3 | 664123..664947 | +0 | +141 | 0.8540 |
| CDS | NC_003450.3 | 668189..670468 | -39 | +0 | 0.9832 |
| CDS | NC_003450.3 | 695074..696279 | +0 | +105 | 0.9199 |
| CDS | NC_003450.3 | 705211..708633 | -51 | +0 | 0.9853 |
| CDS | NC_003450.3 | 708839..709711 | +3 | +0 | 0.9966 |
| CDS | NC_003450.3 | 716283..716660 | +0 | -9 | 0.9762 |
| CDS | NC_003450.3 | 730436..731302 | +57 | +0 | 0.9343 |
| CDS | NC_003450.3 | 733060..734946 | +12 | +0 | 0.9936 |
| CDS | NC_003450.3 | 791476..793011 | +36 | +0 | 0.9766 |
| CDS | NC_003450.3 | 794810..795304 | -96 | +0 | 0.8376 |
| CDS | NC_003450.3 | 807943..809217 | +0 | -48 | 0.9624 |
| CDS | NC_003450.3 | 809528..810397 | +21 | +0 | 0.9759 |
| CDS | NC_003450.3 | 834927..835391 | +63 | +0 | 0.8645 |

- Aviso: 2 feature(s) GFF multi-segmento fusionados por ID; se comparan por su intervalo envolvente.

## DFAST

- Estado: **ok**.
- Fuente: `executed`.
- Versión: `DFAST ver. 1.6.2` (estado `ok`).
- Ejecución: **ok**, 740.33 s, RSS máx. 2653.2 MiB (`2026-11-25T20:50:28Z` → `2026-11-25T21:02:48Z`).
- Salida: `/tmp/external_benchmark_dfast_yxng_5ev/out/genome.gff` (5815130 bytes, SHA256 `ed90432ba4cf8f71…`).
- Comando: `/usr/local/bin/dfast --genome /tmp/benchmark_genome_c12m6ik_.fna --out /tmp/external_benchmark_dfast_yxng_5ev/out --cpu 1 --force --dbroot /data/dfast_db`.
- Anotación comparada: `reports/external_benchmark/multi_v2/cglutamicum/evidence/dfast/genome.gff` (formato `gff3`, SHA256 `ed90432ba4cf8f71…`).

### Features (DFAST)

| tipo | ref | query | TP | FP | FN | precisión | recall | F1 | exactos | borde distinto | solape medio |
|---|---|---|---|---|---|---|---|---|---|---|---|
| CDS | 2948 | 3065 | 2903 | 162 | 45 | 0.9471 | 0.9847 | 0.9656 | 2503 | 400 | 0.9909 |
| tRNA | 59 | 60 | 59 | 1 | 0 | 0.9833 | 1.0000 | 0.9916 | 13 | 46 | 0.9890 |
| rRNA | 18 | 18 | 18 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 0 | 18 | 0.9861 |

### Atributos (DFAST)

| tipo | pares | product en ambos | product literal | product tokens | gene literal | locus literal |
|---|---|---|---|---|---|---|
| CDS | 2903 | 2903 | 0.8440 | 0.8974 | 0.9966 | 0.0000 |
| tRNA | 59 | 59 | 0.9831 | 0.9887 | — | 0.0000 |
| rRNA | 18 | 18 | 1.0000 | 1.0000 | — | 0.0000 |

### Hipotéticas (DFAST)

| tipo | ref CDS | ref con product | ref hyp | % ref | query CDS | query con product | query hyp | % query |
|---|---|---|---|---|---|---|---|---|
| CDS | 2948 | 2948 | 625 | 21.2 | 3065 | 3065 | 817 | 26.7 |
| tRNA | 0 | 0 | 0 | — | 0 | 0 | 0 | — |
| rRNA | 0 | 0 | 0 | — | 0 | 0 | 0 | — |

#### Diferencias

- Falsos negativos (referencia sin pareja): **45**.
- Falsos positivos (herramienta sin pareja): **163**.
- Pares solo con hebra distinta: **5**.
- Pares emparejados con algún borde distinto: **464**.

Las listas de diferencias están truncadas a 100 elementos por lista; los conteos son completos.

> Primeras diferencias no emparejadas

| tipo | feature | seqid | coordenadas | hebra | product/locus |
|---|---|---|---|---|---|
| FN | CDS | NC_003450.3 | 23612..23779 | - | hypothetical protein |
| FN | CDS | NC_003450.3 | 75745..76038 | + | CHY zinc finger protein |
| FN | CDS | NC_003450.3 | 221929..222210 | + | 2Fe-2S iron-sulfur cluster-binding protein |
| FN | CDS | NC_003450.3 | 377683..377835 | + | transposase |
| FN | CDS | NC_003450.3 | 378296..378514 | + | transposase |
| FN | CDS | NC_003450.3 | 674796..675077 | - | putative quinol monooxygenase |
| FN | CDS | NC_003450.3 | 805346..805420 | - | 50S ribosomal protein bL37 |
| FN | CDS | NC_003450.3 | 932305..932490 | + | hypothetical protein |
| FN | CDS | NC_003450.3 | 951556..951654 | + | hypothetical protein |
| FN | CDS | NC_003450.3 | 994650..994841 | + | hypothetical protein |
| FN | CDS | NC_003450.3 | 1015982..1016563 | + | hypothetical protein |
| FN | CDS | NC_003450.3 | 1103951..1104106 | + | hypothetical protein |
| FN | CDS | NC_003450.3 | 1122872..1123030 | + | hypothetical protein |
| FN | CDS | NC_003450.3 | 1280967..1281254 | + | hypothetical protein |
| FN | CDS | NC_003450.3 | 1367167..1367289 | + | hypothetical protein |
| FN | CDS | NC_003450.3 | 1367885..1368070 | - | DoxX family protein |
| FN | CDS | NC_003450.3 | 1371976..1372326 | - | Abi family protein |
| FN | CDS | NC_003450.3 | 1402272..1402706 | + | MepB family protein |
| FN | CDS | NC_003450.3 | 1555835..1556182 | - | hypothetical protein |
| FN | CDS | NC_003450.3 | 1564237..1564485 | + | Fic family protein |
| FN | CDS | NC_003450.3 | 1780504..1780830 | - | alpha/beta hydrolase |
| FN | CDS | NC_003450.3 | 1797182..1797340 | + | hypothetical protein |
| FN | CDS | NC_003450.3 | 1888228..1888599 | - | hypothetical protein |
| FN | CDS | NC_003450.3 | 2098409..2098945 | - | ribosome maturation factor RimP |
| FN | CDS | NC_003450.3 | 2357328..2357711 | + | hypothetical protein |
| FN | CDS | NC_003450.3 | 2385911..2386192 | - | hypothetical protein |
| FN | CDS | NC_003450.3 | 2386210..2386584 | + | hypothetical protein |
| FN | CDS | NC_003450.3 | 2393165..2393407 | + | hypothetical protein |
| FN | CDS | NC_003450.3 | 2609336..2609479 | - | hypothetical protein |
| FN | CDS | NC_003450.3 | 2612456..2612593 | - | hypothetical protein |
| FN | CDS | NC_003450.3 | 2612628..2612897 | - | hypothetical protein |
| FN | CDS | NC_003450.3 | 2627777..2628088 | + | LacI family DNA-binding transcriptional regulator |
| FN | CDS | NC_003450.3 | 2690151..2690441 | + | fluoride efflux transporter family protein |
| FN | CDS | NC_003450.3 | 2840756..2840929 | - | DUF4236 domain-containing protein |
| FN | CDS | NC_003450.3 | 2997479..2997688 | - | hypothetical protein |
| FN | CDS | NC_003450.3 | 3030360..3030539 | + | ABC transporter ATP-binding protein |
| FN | CDS | NC_003450.3 | 3144493..3144633 | - | hypothetical protein |
| FN | CDS | NC_003450.3 | 3160629..3160802 | - | helix-turn-helix domain-containing protein |
| FN | CDS | NC_003450.3 | 3160721..3160939 | - | recombinase family protein |
| FN | CDS | NC_003450.3 | 3161680..3162015 | - | cupin domain-containing protein |
| FN | CDS | NC_003450.3 | 3166265..3166438 | - | hypothetical protein |
| FN | CDS | NC_003450.3 | 3172908..3173090 | - | hypothetical protein |
| FN | CDS | NC_003450.3 | 3233152..3233205 | + | trp operon leader peptide |
| FN | CDS | NC_003450.3 | 3282782..3283105 | + | MGMT family protein |
| FN | CDS | NC_003450.3 | 3307689..3307988 | - | membrane protein insertion efficiency factor YidD |
| FP | CDS | NC_003450.3 | 21426..21527 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 50751..50972 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 75727..75981 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 85426..85623 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 105054..105254 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 106429..106590 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 116207..116434 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 155874..156170 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 186769..187305 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 201338..201580 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 201617..201730 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 222097..222210 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 239793..239948 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 273296..273484 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 332919..334565 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 334516..334956 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 345811..345975 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 378141..378230 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 378284..378538 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 387418..387528 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 452006..452212 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 527604..527759 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 532008..532526 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 671141..671770 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 671736..672656 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 674853..674975 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 748368..748637 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 847727..848029 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 855839..855949 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 895441..895596 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 953747..954148 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 954277..954756 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 954974..955357 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 955354..955500 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 959855..960256 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 960385..960864 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 982291..982674 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1014925..1015119 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 1079037..1079183 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 1084261..1084368 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 1088205..1088507 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1089862..1090008 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1095384..1095677 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1096185..1096331 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1099012..1099164 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1192129..1192257 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1227890..1228639 | + | alpha/beta hydrolase |
| FP | CDS | NC_003450.3 | 1228718..1229098 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 1341957..1342457 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1342458..1342727 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1342791..1343675 | - | ABC transporter substrate-binding protein |
| FP | CDS | NC_003450.3 | 1366270..1366500 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1371855..1371965 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1375930..1376088 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1400881..1401336 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 1402365..1402874 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1404867..1405301 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1408250..1408453 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 1409141..1409284 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 1555835..1555999 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1609839..1610192 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1639370..1639504 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1639656..1639820 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 1640153..1640536 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1652891..1653577 | - | transcriptional regulator |
| FP | CDS | NC_003450.3 | 1663933..1664115 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1671270..1671677 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1673263..1674105 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1743843..1743971 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 1744025..1744522 | + | CE1759 family FMN reductase |
| FP | CDS | NC_003450.3 | 1777332..1777649 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 1777698..1777823 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 1780498..1780671 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 1790831..1790956 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1793548..1793703 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 1798791..1799042 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 1799368..1799505 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1799624..1799743 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 1826641..1826937 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1829685..1829900 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1887932..1888129 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1888228..1889445 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1917561..1918070 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1928534..1928911 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 1930987..1931190 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1937202..1937414 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 1943095..1943313 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 1943345..1943656 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 1945592..1945882 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1967114..1967233 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1971401..1971550 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1987504..1987821 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1987884..1988306 | - | IS3 family transposase |
| FP | CDS | NC_003450.3 | 1988367..1988483 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1988527..1988664 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1994302..1994430 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1994434..1994568 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 1995780..1996019 | - | hypothetical protein |
| FP | CDS | NC_003450.3 | 2013996..2014100 | + | hypothetical protein |
| FP | CDS | NC_003450.3 | 2025218..2025388 | + | hypothetical protein |

> Primeros pares con coordenadas no exactas

| tipo | seqid | ref | Δstart | Δend | solape |
|---|---|---|---|---|---|
| CDS | NC_003450.3 | 7830..8798 | +54 | +0 | 0.9443 |
| CDS | NC_003450.3 | 10104..11177 | +0 | -6 | 0.9944 |
| CDS | NC_003450.3 | 11520..11729 | +0 | +39 | 0.8434 |
| CDS | NC_003450.3 | 11891..14401 | -60 | +0 | 0.9767 |
| CDS | NC_003450.3 | 21564..22127 | +33 | +0 | 0.9415 |
| CDS | NC_003450.3 | 24295..24732 | -81 | +0 | 0.8439 |
| CDS | NC_003450.3 | 34366..34902 | -36 | +0 | 0.9372 |
| CDS | NC_003450.3 | 36244..37014 | +0 | +228 | 0.7718 |
| CDS | NC_003450.3 | 53094..54011 | +3 | +0 | 0.9967 |
| CDS | NC_003450.3 | 60699..62324 | -30 | +0 | 0.9819 |
| CDS | NC_003450.3 | 86399..87244 | -81 | +0 | 0.9126 |
| CDS | NC_003450.3 | 89615..90448 | -42 | +0 | 0.9521 |
| CDS | NC_003450.3 | 103503..104912 | -24 | +0 | 0.9833 |
| CDS | NC_003450.3 | 105170..105892 | +0 | -270 | 0.6266 |
| CDS | NC_003450.3 | 107444..110893 | -9 | +0 | 0.9974 |
| CDS | NC_003450.3 | 115578..115946 | -42 | +0 | 0.8978 |
| CDS | NC_003450.3 | 127189..127998 | +0 | +6 | 0.9926 |
| CDS | NC_003450.3 | 129486..130055 | +0 | +63 | 0.9005 |
| CDS | NC_003450.3 | 132968..134113 | +0 | +21 | 0.9820 |
| CDS | NC_003450.3 | 135515..136201 | +0 | -168 | 0.7555 |
| CDS | NC_003450.3 | 140329..141792 | +3 | +0 | 0.9980 |
| CDS | NC_003450.3 | 144725..145483 | +6 | +0 | 0.9921 |
| CDS | NC_003450.3 | 153223..155613 | +0 | +39 | 0.9840 |
| CDS | NC_003450.3 | 157614..158141 | -30 | +0 | 0.9462 |
| CDS | NC_003450.3 | 169038..169994 | -63 | +0 | 0.9382 |
| CDS | NC_003450.3 | 170044..170919 | -48 | +0 | 0.9481 |
| CDS | NC_003450.3 | 193262..194467 | -405 | +0 | 0.7486 |
| CDS | NC_003450.3 | 220151..221107 | +0 | +24 | 0.9755 |
| CDS | NC_003450.3 | 224354..225247 | +27 | +0 | 0.9698 |
| CDS | NC_003450.3 | 238218..239528 | -42 | +0 | 0.9690 |
| CDS | NC_003450.3 | 241880..242902 | +0 | +18 | 0.9827 |
| CDS | NC_003450.3 | 244466..244819 | +18 | +0 | 0.9492 |
| CDS | NC_003450.3 | 247319..248575 | -39 | +0 | 0.9699 |
| CDS | NC_003450.3 | 256599..257897 | +18 | +0 | 0.9861 |
| CDS | NC_003450.3 | 257900..258532 | +81 | +0 | 0.8720 |
| CDS | NC_003450.3 | 274366..275874 | -42 | +0 | 0.9729 |
| CDS | NC_003450.3 | 278385..278732 | +0 | +102 | 0.7733 |
| CDS | NC_003450.3 | 281401..282927 | +0 | +6 | 0.9961 |
| CDS | NC_003450.3 | 299723..301261 | +0 | -6 | 0.9961 |
| CDS | NC_003450.3 | 305876..306703 | -18 | +0 | 0.9787 |
| CDS | NC_003450.3 | 307724..307909 | +0 | -15 | 0.9194 |
| CDS | NC_003450.3 | 310132..311328 | +3 | +0 | 0.9975 |
| CDS | NC_003450.3 | 313622..314464 | +0 | +126 | 0.8700 |
| CDS | NC_003450.3 | 319333..321690 | +0 | +39 | 0.9837 |
| CDS | NC_003450.3 | 325919..326617 | -15 | +0 | 0.9790 |
| CDS | NC_003450.3 | 330573..330812 | +18 | +0 | 0.9250 |
| CDS | NC_003450.3 | 337539..338771 | +21 | +0 | 0.9830 |
| CDS | NC_003450.3 | 352693..353640 | +153 | +0 | 0.8386 |
| CDS | NC_003450.3 | 357234..359240 | +0 | +114 | 0.9463 |
| CDS | NC_003450.3 | 359759..360307 | +0 | +9 | 0.9839 |
| CDS | NC_003450.3 | 370685..371776 | -72 | +0 | 0.9381 |
| CDS | NC_003450.3 | 378665..379837 | +0 | -117 | 0.9003 |
| CDS | NC_003450.3 | 390784..392208 | +0 | -15 | 0.9895 |
| CDS | NC_003450.3 | 398278..399582 | -30 | +0 | 0.9775 |
| CDS | NC_003450.3 | 400039..400344 | +12 | +0 | 0.9608 |
| CDS | NC_003450.3 | 404505..405419 | +0 | +9 | 0.9903 |
| CDS | NC_003450.3 | 433062..433991 | +21 | +0 | 0.9774 |
| CDS | NC_003450.3 | 438516..439907 | +42 | +0 | 0.9698 |
| CDS | NC_003450.3 | 452637..454433 | -57 | +0 | 0.9693 |
| CDS | NC_003450.3 | 460020..461096 | -48 | +0 | 0.9573 |
| CDS | NC_003450.3 | 462557..463870 | +9 | +0 | 0.9932 |
| CDS | NC_003450.3 | 469370..470173 | +60 | +0 | 0.9254 |
| CDS | NC_003450.3 | 478089..478973 | +0 | -84 | 0.9051 |
| CDS | NC_003450.3 | 480621..480995 | +0 | +6 | 0.9843 |
| CDS | NC_003450.3 | 481391..482668 | +0 | +108 | 0.9221 |
| CDS | NC_003450.3 | 490580..491941 | +60 | +0 | 0.9559 |
| CDS | NC_003450.3 | 511536..512510 | -33 | +0 | 0.9673 |
| CDS | NC_003450.3 | 512913..516410 | +18 | +0 | 0.9949 |
| CDS | NC_003450.3 | 516494..520495 | -27 | +0 | 0.9933 |
| CDS | NC_003450.3 | 544754..547042 | +0 | -30 | 0.9869 |
| CDS | NC_003450.3 | 550696..551808 | +0 | +15 | 0.9867 |
| CDS | NC_003450.3 | 554449..554973 | +0 | -39 | 0.9257 |
| CDS | NC_003450.3 | 599759..599935 | -60 | +0 | 0.7468 |
| CDS | NC_003450.3 | 614800..615180 | +0 | +42 | 0.9007 |
| CDS | NC_003450.3 | 618090..618995 | +0 | +18 | 0.9805 |
| CDS | NC_003450.3 | 623800..624942 | -129 | +0 | 0.8986 |
| CDS | NC_003450.3 | 624985..625677 | -42 | +0 | 0.9429 |
| CDS | NC_003450.3 | 626067..626510 | +0 | +33 | 0.9308 |
| CDS | NC_003450.3 | 628563..630143 | -12 | +0 | 0.9925 |
| CDS | NC_003450.3 | 633076..633525 | +0 | -51 | 0.8867 |
| CDS | NC_003450.3 | 636086..638224 | +0 | +54 | 0.9754 |
| CDS | NC_003450.3 | 638314..639237 | +0 | +183 | 0.8347 |
| CDS | NC_003450.3 | 641037..642560 | -84 | +0 | 0.9478 |
| CDS | NC_003450.3 | 645173..647470 | +0 | +120 | 0.9504 |
| CDS | NC_003450.3 | 654676..655125 | -42 | +0 | 0.9146 |
| CDS | NC_003450.3 | 657202..657978 | +0 | +45 | 0.9453 |
| CDS | NC_003450.3 | 661166..662020 | -45 | +0 | 0.9500 |
| CDS | NC_003450.3 | 664123..664947 | +0 | +141 | 0.8540 |
| CDS | NC_003450.3 | 668189..670468 | -39 | +0 | 0.9832 |
| CDS | NC_003450.3 | 695074..696279 | +0 | +105 | 0.9199 |
| CDS | NC_003450.3 | 705211..708633 | -51 | +0 | 0.9853 |
| CDS | NC_003450.3 | 708839..709711 | +3 | +0 | 0.9966 |
| CDS | NC_003450.3 | 716283..716660 | +0 | -9 | 0.9762 |
| CDS | NC_003450.3 | 730436..731302 | +57 | +0 | 0.9343 |
| CDS | NC_003450.3 | 733060..734946 | +12 | +0 | 0.9936 |
| CDS | NC_003450.3 | 791476..793011 | +36 | +0 | 0.9766 |
| CDS | NC_003450.3 | 794810..795304 | -96 | +0 | 0.8376 |
| CDS | NC_003450.3 | 807943..809217 | +0 | -48 | 0.9624 |
| CDS | NC_003450.3 | 809528..810397 | +21 | +0 | 0.9759 |
| CDS | NC_003450.3 | 834927..835391 | +63 | +0 | 0.8645 |

- Aviso: 2 feature(s) GFF multi-segmento fusionados por ID; se comparan por su intervalo envolvente.

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
/opt/conda/bin/python benchmarks/external/external_benchmark.py --genome data/reference_genomes/cglutamicum/cglutamicum_genomic.fna.gz --reference data/reference_genomes/cglutamicum/cglutamicum_genomic.gff.gz --tools bakta,prokka,dfast --run --threads 1 --json reports/external_benchmark/benchmark.json --md reports/external_benchmark/benchmark.md --bakta-db /data/bakta_db
```

Para comparar salidas producidas en otro host, usar `--query-gff TOOL=RUTA` sin `--run`; el origen queda como `precomputed`.
