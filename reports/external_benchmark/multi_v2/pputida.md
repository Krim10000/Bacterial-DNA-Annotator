# Benchmark externo CPU: Bakta, Prokka y DFAST

**Generado:** 2026-11-25T22:53:01Z

**Esquema:** `external-benchmark-run/v1`

**Alcance:** Comparación local de CDS/tRNA/rRNA de Bakta, Prokka y DFAST contra una referencia GFF3/GBFF. CPU, sin descargas de datos.

**Comando reproducible:**

```bash
/opt/conda/bin/python /opt/benchmark/external/external_benchmark.py --genome /work/data/reference_genomes/pputida/pputida_genomic.fna.gz --reference /work/data/reference_genomes/pputida/pputida_genomic.gff.gz --tools all --run --threads 4 --bakta-db /data/bakta_db --extra-arg dfast=--dbroot --extra-arg dfast=/data/dfast_db --save-gffs /work/reports/external_benchmark/multi_v2/pputida/evidence --json /work/reports/external_benchmark/multi_v2/pputida.json --md /work/reports/external_benchmark/multi_v2/pputida.md
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

- FASTA: `data/reference_genomes/pputida/pputida_genomic.fna.gz` (SHA256 `3905190150e291a5…`, 1 contig, 6156701 bp).
- Referencia: `data/reference_genomes/pputida/pputida_genomic.gff.gz` (formato `gff3`, SHA256 `8f8794b91cce6066…`).
- Conteos de referencia: CDS **5377**, tRNA **75**, rRNA **22**.
- CDS pseudo excluidos en referencia: **94**.
- SHA256 FASTA: `3905190150e291a5fa5c04d491d29bf8507c2d09f4d0471732c0102c6d3844d7`.
- SHA256 referencia: `8f8794b91cce6066124a70c2af3a2e181b6d753b7d5281a581c41957878bc7ba`.

## Validación

- Inicio: `2026-11-25T22:16:39Z`; fin: `2026-11-25T22:53:01Z`; total: **2181.98 s**.
- JSON declarado: `reports/external_benchmark/multi_v2/pputida.json`.
- Markdown declarado: `reports/external_benchmark/multi_v2/pputida.md`.

| check | estado | detalle |
|---|---|---|
| genome_path | ok | FASTA existe y es legible |
| genome_sha256 | ok | SHA256 3905190150e291a5fa5c04d491d29bf8507c2d09f4d0471732c0102c6d3844d7 |
| reference_path | ok | referencia gff3 con 5474 features |
| reference_sha256 | ok | SHA256 8f8794b91cce6066124a70c2af3a2e181b6d753b7d5281a581c41957878bc7ba |
| reference_features | ok | CDS=5377 tRNA=75 rRNA=22 pseudo_excluidos=94 |
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
| Bakta | executed | ok | 0.9806 | 1.0000 | 1.0000 | 0.9761 | 0.9858 | 0.9809 | 966.4 |
| Prokka | executed | ok | 0.9800 | 0.9677 | 1.0000 | 0.9750 | 0.9848 | 0.9799 | 762.6 |
| DFAST | executed | ok | 0.9806 | 0.9677 | 1.0000 | 0.9752 | 0.9858 | 0.9805 | 398.8 |

La precisión/recall/F1 se calculan contra la referencia local. Si una herramienta no está instalada, su estado es `not_available` y no hay métricas que mostrar.

## Bakta

- Estado: **ok**.
- Fuente: `executed`.
- Bakta DB: `ok` (válida: **sí**, ruta `/data/bakta_db`, versión `0`).
- Versión: `bakta 1.12.2` (estado `ok`).
- Ejecución: **ok**, 966.41 s, RSS máx. 2676.5 MiB (`2026-11-25T22:16:46Z` → `2026-11-25T22:32:52Z`).
- Salida: `/tmp/external_benchmark_bakta_n6yla9bx/out/benchmark.gff3` (7429429 bytes, SHA256 `f99061371e6bfa0f…`).
- Comando: `/usr/local/bin/bakta --output /tmp/external_benchmark_bakta_n6yla9bx/out --prefix benchmark --threads 4 --force --skip-plot --db /data/bakta_db /tmp/benchmark_genome_b9aq49k3.fna`.
- Anotación comparada: `reports/external_benchmark/multi_v2/pputida/evidence/bakta/benchmark.gff3` (formato `gff3`, SHA256 `f99061371e6bfa0f…`).

### Features (Bakta)

| tipo | ref | query | TP | FP | FN | precisión | recall | F1 | exactos | borde distinto | solape medio |
|---|---|---|---|---|---|---|---|---|---|---|---|
| CDS | 5377 | 5431 | 5299 | 132 | 78 | 0.9757 | 0.9855 | 0.9806 | 4977 | 322 | 0.9960 |
| tRNA | 75 | 75 | 75 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 75 | 0 | 1.0000 |
| rRNA | 22 | 22 | 22 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 22 | 0 | 1.0000 |

### Atributos (Bakta)

| tipo | pares | product en ambos | product literal | product tokens | gene literal | locus literal |
|---|---|---|---|---|---|---|
| CDS | 5299 | 5299 | 0.2200 | 0.4202 | 0.8792 | 0.0000 |
| tRNA | 75 | 75 | 0.0000 | 0.6389 | — | 0.0000 |
| rRNA | 22 | 22 | 1.0000 | 1.0000 | 1.0000 | 0.0000 |

### Hipotéticas (Bakta)

| tipo | ref CDS | ref con product | ref hyp | % ref | query CDS | query con product | query hyp | % query |
|---|---|---|---|---|---|---|---|---|
| CDS | 5377 | 5377 | 677 | 12.6 | 5431 | 5431 | 746 | 13.7 |
| tRNA | 0 | 0 | 0 | — | 0 | 0 | 0 | — |
| rRNA | 0 | 0 | 0 | — | 0 | 0 | 0 | — |

#### Diferencias

- Falsos negativos (referencia sin pareja): **78**.
- Falsos positivos (herramienta sin pareja): **132**.
- Pares solo con hebra distinta: **1**.
- Pares emparejados con algún borde distinto: **322**.

Las listas de diferencias están truncadas a 100 elementos por lista; los conteos son completos.

> Primeras diferencias no emparejadas

| tipo | feature | seqid | coordenadas | hebra | product/locus |
|---|---|---|---|---|---|
| FN | CDS | NC_021505.1 | 17890..18099 | + | hypothetical protein |
| FN | CDS | NC_021505.1 | 20903..21241 | + | zinc ribbon domain-containing protein |
| FN | CDS | NC_021505.1 | 29460..29879 | - | GtrA family protein |
| FN | CDS | NC_021505.1 | 333665..333862 | + | hypothetical protein |
| FN | CDS | NC_021505.1 | 344909..345034 | + | type II toxin-antitoxin system HigB family toxin |
| FN | CDS | NC_021505.1 | 347261..347839 | + | hypothetical protein |
| FN | CDS | NC_021505.1 | 483140..483211 | - | pyrroloquinoline quinone precursor peptide PqqA |
| FN | CDS | NC_021505.1 | 488761..489042 | - | hypothetical protein |
| FN | CDS | NC_021505.1 | 489540..489701 | + | hypothetical protein |
| FN | CDS | NC_021505.1 | 490974..491216 | + | DUF4160 domain-containing protein |
| FN | CDS | NC_021505.1 | 657757..658194 | + | antitoxin Xre/MbcA/ParS toxin-binding domain-containing protein |
| FN | CDS | NC_021505.1 | 659624..659848 | - | helix-turn-helix domain-containing protein |
| FN | CDS | NC_021505.1 | 702006..702368 | - | hypothetical protein |
| FN | CDS | NC_021505.1 | 702416..702622 | - | hypothetical protein |
| FN | CDS | NC_021505.1 | 769837..770604 | - | methyl-accepting chemotaxis protein |
| FN | CDS | NC_021505.1 | 776910..777107 | - | hypothetical protein |
| FN | CDS | NC_021505.1 | 822983..823846 | - | methyl-accepting chemotaxis protein |
| FN | CDS | NC_021505.1 | 1054460..1055224 | - | methyl-accepting chemotaxis protein |
| FN | CDS | NC_021505.1 | 1107034..1107246 | - | cysteine-rich CWC family protein |
| FN | CDS | NC_021505.1 | 1121505..1121963 | + | protein YgfX |
| FN | CDS | NC_021505.1 | 1308343..1308987 | + | hypothetical protein |
| FN | CDS | NC_021505.1 | 1403753..1404031 | + | hypothetical protein |
| FN | CDS | NC_021505.1 | 1408226..1408459 | + | hypothetical protein |
| FN | CDS | NC_021505.1 | 1408691..1409125 | + | FliH/SctL family protein |
| FN | CDS | NC_021505.1 | 1461733..1461894 | + | hypothetical protein |
| FN | CDS | NC_021505.1 | 1477667..1477921 | + | type II toxin-antitoxin system HicA family toxin |
| FN | CDS | NC_021505.1 | 1579480..1580889 | - | amino acid permease |
| FN | CDS | NC_021505.1 | 1580979..1581233 | - | hypothetical protein |
| FN | CDS | NC_021505.1 | 1738514..1738708 | - | hypothetical protein |
| FN | CDS | NC_021505.1 | 1887143..1887406 | + | ogr/Delta-like zinc finger family protein |
| FN | CDS | NC_021505.1 | 1941868..1942266 | + | type III secretion system chaperone |
| FN | CDS | NC_021505.1 | 2051692..2052342 | + | methyl-accepting chemotaxis protein |
| FN | CDS | NC_021505.1 | 2086549..2086770 | + | hypothetical protein |
| FN | CDS | NC_021505.1 | 2090570..2090794 | - | hypothetical protein |
| FN | CDS | NC_021505.1 | 2091737..2092078 | + | hypothetical protein |
| FN | CDS | NC_021505.1 | 2092409..2092888 | + | hypothetical protein |
| FN | CDS | NC_021505.1 | 2097312..2097554 | - | Pathogenicity locus |
| FN | CDS | NC_021505.1 | 2253759..2254148 | + | hypothetical protein |
| FN | CDS | NC_021505.1 | 2286392..2286691 | + | hypothetical protein |
| FN | CDS | NC_021505.1 | 2455375..2455701 | + | winged helix-turn-helix domain-containing protein |
| FN | CDS | NC_021505.1 | 2516256..2516543 | - | DUF4189 domain-containing protein |
| FN | CDS | NC_021505.1 | 2587063..2587818 | - | methyl-accepting chemotaxis protein |
| FN | CDS | NC_021505.1 | 2623632..2623856 | - | DUF7683 domain-containing protein |
| FN | CDS | NC_021505.1 | 2801116..2801475 | + | ADP-ribosylglycohydrolase family protein |
| FN | CDS | NC_021505.1 | 2802766..2803287 | + | hypothetical protein |
| FN | CDS | NC_021505.1 | 2970688..2971005 | - | GNAT family N-acetyltransferase |
| FN | CDS | NC_021505.1 | 2974652..2974864 | - | hypothetical protein |
| FN | CDS | NC_021505.1 | 2979306..2979488 | + | DUF2256 domain-containing protein |
| FN | CDS | NC_021505.1 | 3206194..3206688 | - | methyl-accepting chemotaxis protein |
| FN | CDS | NC_021505.1 | 3516873..3517016 | + | hypothetical protein |
| FN | CDS | NC_021505.1 | 3588358..3588429 | - | pyrroloquinoline quinone precursor peptide PqqA |
| FN | CDS | NC_021505.1 | 3645207..3645590 | - | APC family permease |
| FN | CDS | NC_021505.1 | 3653738..3653983 | + | hypothetical protein |
| FN | CDS | NC_021505.1 | 3829786..3830190 | + | hypothetical protein |
| FN | CDS | NC_021505.1 | 4096371..4096757 | - | hypothetical protein |
| FN | CDS | NC_021505.1 | 4142347..4142481 | + | hypothetical protein |
| FN | CDS | NC_021505.1 | 4154132..4154683 | - | DUF1090 family protein |
| FN | CDS | NC_021505.1 | 4238829..4239596 | - | methyl-accepting chemotaxis protein |
| FN | CDS | NC_021505.1 | 4406574..4406801 | - | hypothetical protein |
| FN | CDS | NC_021505.1 | 4498270..4499151 | - | methyl-accepting chemotaxis protein |
| FN | CDS | NC_021505.1 | 5110048..5110479 | - | methyl-accepting chemotaxis protein |
| FN | CDS | NC_021505.1 | 5114761..5115063 | + | DUF3565 domain-containing protein |
| FN | CDS | NC_021505.1 | 5233828..5234139 | + | hypothetical protein |
| FN | CDS | NC_021505.1 | 5235958..5236188 | - | hypothetical protein |
| FN | CDS | NC_021505.1 | 5395796..5396197 | - | hypothetical protein |
| FN | CDS | NC_021505.1 | 5402778..5403152 | + | (R)-mandelonitrile lyase |
| FN | CDS | NC_021505.1 | 5431418..5431687 | - | hypothetical protein |
| FN | CDS | NC_021505.1 | 5442530..5443009 | - | hypothetical protein |
| FN | CDS | NC_021505.1 | 5467841..5468236 | - | hypothetical protein |
| FN | CDS | NC_021505.1 | 5469103..5469498 | - | DUF6957 family protein |
| FN | CDS | NC_021505.1 | 5601999..5602883 | - | methyl-accepting chemotaxis protein |
| FN | CDS | NC_021505.1 | 5672247..5672666 | - | hypothetical protein |
| FN | CDS | NC_021505.1 | 5672675..5672989 | - | hypothetical protein |
| FN | CDS | NC_021505.1 | 5758049..5758807 | + | methyl-accepting chemotaxis protein |
| FN | CDS | NC_021505.1 | 5939374..5939829 | + | REP-associated tyrosine transposase |
| FN | CDS | NC_021505.1 | 6011288..6011566 | - | type II toxin-antitoxin system RelE/ParE family toxin |
| FN | CDS | NC_021505.1 | 6034453..6034803 | - | RHS repeat-associated core domain-containing protein |
| FN | CDS | NC_021505.1 | 6155699..6155944 | - | membrane protein insertion efficiency factor YidD |
| FP | CDS | NC_021505.1 | 81549..82484 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 125613..126467 | + | Endonuclease/exonuclease/phosphatase domain-containing protein |
| FP | CDS | NC_021505.1 | 180064..180327 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 180352..180819 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 186522..186668 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 213738..214403 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 347588..347839 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 456764..457207 | + | Response regulator |
| FP | CDS | NC_021505.1 | 658640..658804 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 705296..705565 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 715252..715350 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 741009..742970 | - | Methyl-accepting chemotaxis protein CtpL |
| FP | CDS | NC_021505.1 | 769837..771780 | - | Methyl-accepting chemotaxis protein |
| FP | CDS | NC_021505.1 | 822983..824902 | - | Methyl-accepting chemotaxis protein McpS |
| FP | CDS | NC_021505.1 | 950118..950207 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 981528..982448 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 982412..982789 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 1051215..1051448 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 1054460..1056334 | - | Methyl-accepting chemotaxis protein PctA |
| FP | CDS | NC_021505.1 | 1081662..1081835 | - | Lipoprotein |
| FP | CDS | NC_021505.1 | 1155617..1157155 | - | glutamate synthase (NADPH) |
| FP | CDS | NC_021505.1 | 1271388..1272137 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 1293060..1293191 | - | Alternative ribosome-rescue factor |
| FP | CDS | NC_021505.1 | 1296350..1296922 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 1308373..1308681 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 1308712..1308948 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 1315936..1316133 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 1436805..1440158 | + | Carrier domain-containing protein |
| FP | CDS | NC_021505.1 | 1489722..1489922 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 1580950..1581111 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 1657104..1657202 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 1683265..1683447 | + | AbrB/MazE/SpoVT family DNA-binding domain-containing protein |
| FP | CDS | NC_021505.1 | 1923228..1923554 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 1923559..1923867 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 1939322..1939822 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 1958442..1958573 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 1995395..1995490 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 2050261..2052342 | + | Methyl-accepting chemotaxis protein |
| FP | CDS | NC_021505.1 | 2053347..2053487 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 2087527..2087616 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 2088715..2088936 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 2097642..2098154 | + | putative metallopeptidase |
| FP | CDS | NC_021505.1 | 2211640..2212002 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 2234371..2234643 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 2470127..2470336 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 2587063..2589207 | - | Methyl-accepting chemotaxis protein McpU |
| FP | CDS | NC_021505.1 | 2591859..2592104 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 2624598..2624987 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 2638781..2638897 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 2676068..2676319 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 2680237..2680590 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 2765658..2765777 | - | Lipoprotein |
| FP | CDS | NC_021505.1 | 2803544..2803642 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 2815790..2816932 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 2818252..2819424 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 2903466..2904620 | + | Elongation factor GreAB |
| FP | CDS | NC_021505.1 | 2932339..2932446 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 2945929..2946105 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 2954493..2956055 | - | Acetate permease ActP |
| FP | CDS | NC_021505.1 | 2987476..2987784 | - | Transposase |
| FP | CDS | NC_021505.1 | 2995639..2995782 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 3162432..3162758 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 3173980..3174222 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 3206194..3207513 | - | Biofilm dispersion protein BdlA |
| FP | CDS | NC_021505.1 | 3257315..3257413 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 3332713..3333984 | + | Major facilitator superfamily (MFS) profile domain-containing protein |
| FP | CDS | NC_021505.1 | 3414088..3414246 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 3420611..3421051 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 3424128..3424478 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 3435041..3435205 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 3483868..3484614 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 3487051..3487371 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 3487311..3487604 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 3515175..3515303 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 3559390..3559857 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 3559891..3560814 | - | Amino acid permease/ SLC12A domain-containing protein |
| FP | CDS | NC_021505.1 | 3584428..3584769 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 3584868..3585329 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 3644767..3645210 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 3652974..3653726 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 3674140..3674805 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 3674802..3675497 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 3683260..3685575 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 3685928..3686188 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 3686400..3686882 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 3686711..3687328 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 3738942..3739184 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 3774303..3776261 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 3776251..3776877 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 3803439..3803816 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 3851781..3852041 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 3889254..3889907 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 3889962..3893675 | - | RING-type E3 ubiquitin transferase |
| FP | CDS | NC_021505.1 | 3916682..3916795 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 4017861..4019426 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 4142341..4142631 | - | DUF4258 domain-containing protein |
| FP | CDS | NC_021505.1 | 4170392..4170742 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 4170696..4170881 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 4238829..4240505 | - | Methyl-accepting chemotaxis protein |
| FP | CDS | NC_021505.1 | 4353945..4354175 | - | hypothetical protein |

> Primeros pares con coordenadas no exactas

| tipo | seqid | ref | Δstart | Δend | solape |
|---|---|---|---|---|---|
| CDS | NC_021505.1 | 14880..15563 | +216 | +0 | 0.6842 |
| CDS | NC_021505.1 | 15566..16435 | +84 | +0 | 0.9034 |
| CDS | NC_021505.1 | 30137..31105 | +150 | +0 | 0.8452 |
| CDS | NC_021505.1 | 59393..60178 | +0 | -15 | 0.9809 |
| CDS | NC_021505.1 | 72329..72643 | -3 | +0 | 0.9906 |
| CDS | NC_021505.1 | 76477..78042 | +0 | +27 | 0.9831 |
| CDS | NC_021505.1 | 93346..95433 | +36 | +0 | 0.9828 |
| CDS | NC_021505.1 | 95430..95711 | +15 | +0 | 0.9468 |
| CDS | NC_021505.1 | 117345..118301 | +48 | +0 | 0.9498 |
| CDS | NC_021505.1 | 119858..122605 | -51 | +0 | 0.9818 |
| CDS | NC_021505.1 | 175007..177421 | +0 | +72 | 0.9710 |
| CDS | NC_021505.1 | 182924..184282 | +0 | +18 | 0.9869 |
| CDS | NC_021505.1 | 226648..227286 | +0 | -48 | 0.9249 |
| CDS | NC_021505.1 | 280509..282404 | +30 | +0 | 0.9842 |
| CDS | NC_021505.1 | 315010..315492 | +0 | +57 | 0.8944 |
| CDS | NC_021505.1 | 324068..324898 | +0 | -45 | 0.9458 |
| CDS | NC_021505.1 | 372257..373105 | +0 | -3 | 0.9965 |
| CDS | NC_021505.1 | 399978..400850 | +0 | +753 | 0.5369 |
| CDS | NC_021505.1 | 415236..415997 | +0 | -3 | 0.9961 |
| CDS | NC_021505.1 | 421026..421979 | +6 | +0 | 0.9937 |
| CDS | NC_021505.1 | 448341..449453 | +0 | +39 | 0.9661 |
| CDS | NC_021505.1 | 462819..463583 | +0 | +63 | 0.9239 |
| CDS | NC_021505.1 | 490335..490889 | -33 | +0 | 0.9439 |
| CDS | NC_021505.1 | 513126..514445 | +0 | +6 | 0.9955 |
| CDS | NC_021505.1 | 519858..520856 | +0 | -30 | 0.9700 |
| CDS | NC_021505.1 | 523347..523979 | +3 | +0 | 0.9953 |
| CDS | NC_021505.1 | 555490..556449 | -174 | +0 | 0.8466 |
| CDS | NC_021505.1 | 559885..560418 | +54 | +0 | 0.8989 |
| CDS | NC_021505.1 | 640311..642206 | +0 | +102 | 0.9489 |
| CDS | NC_021505.1 | 644613..645839 | +21 | +0 | 0.9829 |
| CDS | NC_021505.1 | 663802..665064 | -36 | +0 | 0.9723 |
| CDS | NC_021505.1 | 666485..667252 | -36 | +0 | 0.9552 |
| CDS | NC_021505.1 | 705971..707173 | -18 | +0 | 0.9853 |
| CDS | NC_021505.1 | 728793..729857 | +0 | +24 | 0.9780 |
| CDS | NC_021505.1 | 781766..782686 | +0 | -3 | 0.9967 |
| CDS | NC_021505.1 | 837135..837842 | +0 | +51 | 0.9328 |
| CDS | NC_021505.1 | 863512..864294 | +0 | +6 | 0.9924 |
| CDS | NC_021505.1 | 876353..876772 | -6 | +0 | 0.9859 |
| CDS | NC_021505.1 | 889273..889818 | +0 | -51 | 0.9066 |
| CDS | NC_021505.1 | 935362..935607 | -3 | +0 | 0.9880 |
| CDS | NC_021505.1 | 954278..955198 | +0 | +78 | 0.9219 |
| CDS | NC_021505.1 | 958379..959590 | -3 | +0 | 0.9975 |
| CDS | NC_021505.1 | 967254..967709 | +51 | +0 | 0.8882 |
| CDS | NC_021505.1 | 968650..968943 | +45 | +0 | 0.8469 |
| CDS | NC_021505.1 | 968943..970682 | -3 | +0 | 0.9983 |
| CDS | NC_021505.1 | 985118..985573 | +0 | +261 | 0.6360 |
| CDS | NC_021505.1 | 1025956..1026909 | +0 | -6 | 0.9937 |
| CDS | NC_021505.1 | 1053100..1054221 | +0 | -9 | 0.9920 |
| CDS | NC_021505.1 | 1060660..1061469 | -66 | +0 | 0.9247 |
| CDS | NC_021505.1 | 1092383..1093702 | +0 | +96 | 0.9322 |
| CDS | NC_021505.1 | 1097167..1099458 | +0 | -3 | 0.9987 |
| CDS | NC_021505.1 | 1126525..1127925 | -33 | +0 | 0.9770 |
| CDS | NC_021505.1 | 1133415..1134137 | -18 | +0 | 0.9757 |
| CDS | NC_021505.1 | 1175885..1177414 | +0 | +42 | 0.9733 |
| CDS | NC_021505.1 | 1188638..1189285 | +0 | +27 | 0.9600 |
| CDS | NC_021505.1 | 1201632..1202585 | +6 | +0 | 0.9937 |
| CDS | NC_021505.1 | 1219842..1220318 | +0 | +3 | 0.9938 |
| CDS | NC_021505.1 | 1234281..1235888 | +0 | +24 | 0.9853 |
| CDS | NC_021505.1 | 1260499..1261572 | +0 | +3 | 0.9972 |
| CDS | NC_021505.1 | 1261703..1262794 | +33 | +0 | 0.9698 |
| CDS | NC_021505.1 | 1309045..1309839 | -51 | +0 | 0.9397 |
| CDS | NC_021505.1 | 1320969..1322027 | -6 | +0 | 0.9944 |
| CDS | NC_021505.1 | 1329479..1329700 | +42 | +0 | 0.8108 |
| CDS | NC_021505.1 | 1339634..1340440 | +180 | +0 | 0.7770 |
| CDS | NC_021505.1 | 1435648..1436592 | +6 | +0 | 0.9937 |
| CDS | NC_021505.1 | 1477918..1478250 | +9 | +0 | 0.9730 |
| CDS | NC_021505.1 | 1496842..1497969 | -39 | +0 | 0.9666 |
| CDS | NC_021505.1 | 1509108..1509410 | +0 | -96 | 0.6832 |
| CDS | NC_021505.1 | 1521625..1522038 | +0 | -48 | 0.8841 |
| CDS | NC_021505.1 | 1531447..1532508 | +81 | +0 | 0.9237 |
| CDS | NC_021505.1 | 1534707..1536047 | +0 | -6 | 0.9955 |
| CDS | NC_021505.1 | 1537748..1538341 | +0 | +39 | 0.9384 |
| CDS | NC_021505.1 | 1539283..1539912 | -36 | +0 | 0.9459 |
| CDS | NC_021505.1 | 1540038..1540763 | -129 | +0 | 0.8491 |
| CDS | NC_021505.1 | 1543505..1544923 | +15 | +0 | 0.9894 |
| CDS | NC_021505.1 | 1545105..1545989 | +0 | -27 | 0.9695 |
| CDS | NC_021505.1 | 1578049..1578375 | +0 | -9 | 0.9725 |
| CDS | NC_021505.1 | 1588884..1590122 | -81 | +0 | 0.9386 |
| CDS | NC_021505.1 | 1619147..1620673 | -3 | +0 | 0.9980 |
| CDS | NC_021505.1 | 1620675..1622504 | +93 | +0 | 0.9492 |
| CDS | NC_021505.1 | 1624121..1625200 | +36 | +0 | 0.9667 |
| CDS | NC_021505.1 | 1627787..1628728 | -6 | +0 | 0.9937 |
| CDS | NC_021505.1 | 1709527..1710312 | +0 | +69 | 0.9193 |
| CDS | NC_021505.1 | 1739201..1739944 | +0 | +24 | 0.9688 |
| CDS | NC_021505.1 | 1753247..1753438 | +0 | +9 | 0.9552 |
| CDS | NC_021505.1 | 1795414..1796706 | +0 | +21 | 0.9840 |
| CDS | NC_021505.1 | 1807714..1809396 | +0 | -18 | 0.9893 |
| CDS | NC_021505.1 | 1816361..1816966 | -3 | +0 | 0.9951 |
| CDS | NC_021505.1 | 1861660..1862895 | +12 | +0 | 0.9903 |
| CDS | NC_021505.1 | 1885612..1885860 | -66 | +0 | 0.7905 |
| CDS | NC_021505.1 | 1905917..1906681 | +0 | -30 | 0.9608 |
| CDS | NC_021505.1 | 1928607..1929131 | -246 | +0 | 0.6809 |
| CDS | NC_021505.1 | 1929238..1930797 | -102 | +0 | 0.9386 |
| CDS | NC_021505.1 | 1937608..1938270 | +15 | +0 | 0.9774 |
| CDS | NC_021505.1 | 1940427..1941146 | +54 | +0 | 0.9250 |
| CDS | NC_021505.1 | 1957845..1958501 | -9 | +0 | 0.9865 |
| CDS | NC_021505.1 | 1980999..1981283 | +0 | -3 | 0.9895 |
| CDS | NC_021505.1 | 2017226..2017654 | +0 | +54 | 0.8882 |
| CDS | NC_021505.1 | 2041217..2042218 | +0 | +15 | 0.9853 |
| CDS | NC_021505.1 | 2072492..2073916 | +0 | -15 | 0.9895 |

- Aviso: 1 feature(s) GFF multi-segmento fusionados por ID; se comparan por su intervalo envolvente.

## Prokka

- Estado: **ok**.
- Fuente: `executed`.
- Versión: `prokka 1.15.6` (estado `ok`).
- Ejecución: **ok**, 762.61 s, RSS máx. 3958.8 MiB (`2026-11-25T22:33:07Z` → `2026-11-25T22:45:50Z`).
- Salida: `/tmp/external_benchmark_prokka_g_li_hpw/out/benchmark.gff` (7531815 bytes, SHA256 `f6763eadc657b957…`).
- Comando: `/usr/local/bin/prokka --outdir /tmp/external_benchmark_prokka_g_li_hpw/out --prefix benchmark --cpus 4 --force --quiet /tmp/benchmark_genome_b9aq49k3.fna`.
- Anotación comparada: `reports/external_benchmark/multi_v2/pputida/evidence/prokka/benchmark.gff` (formato `gff3`, SHA256 `f6763eadc657b957…`).

### Features (Prokka)

| tipo | ref | query | TP | FP | FN | precisión | recall | F1 | exactos | borde distinto | solape medio |
|---|---|---|---|---|---|---|---|---|---|---|---|
| CDS | 5377 | 5427 | 5294 | 133 | 83 | 0.9755 | 0.9846 | 0.9800 | 4969 | 325 | 0.9959 |
| tRNA | 75 | 80 | 75 | 5 | 0 | 0.9375 | 1.0000 | 0.9677 | 75 | 0 | 1.0000 |
| rRNA | 22 | 22 | 22 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 22 | 0 | 1.0000 |

### Atributos (Prokka)

| tipo | pares | product en ambos | product literal | product tokens | gene literal | locus literal |
|---|---|---|---|---|---|---|
| CDS | 5294 | 5294 | 0.1653 | 0.3493 | 0.6628 | 0.0000 |
| tRNA | 75 | 75 | 0.0000 | 0.6611 | — | 0.0000 |
| rRNA | 22 | 22 | 1.0000 | 1.0000 | — | 0.0000 |

### Hipotéticas (Prokka)

| tipo | ref CDS | ref con product | ref hyp | % ref | query CDS | query con product | query hyp | % query |
|---|---|---|---|---|---|---|---|---|
| CDS | 5377 | 5377 | 677 | 12.6 | 5427 | 5427 | 1918 | 35.3 |
| tRNA | 0 | 0 | 0 | — | 0 | 0 | 0 | — |
| rRNA | 0 | 0 | 0 | — | 0 | 0 | 0 | — |

#### Diferencias

- Falsos negativos (referencia sin pareja): **83**.
- Falsos positivos (herramienta sin pareja): **138**.
- Pares solo con hebra distinta: **1**.
- Pares emparejados con algún borde distinto: **325**.

Las listas de diferencias están truncadas a 100 elementos por lista; los conteos son completos.

> Primeras diferencias no emparejadas

| tipo | feature | seqid | coordenadas | hebra | product/locus |
|---|---|---|---|---|---|
| FN | CDS | NC_021505.1 | 17890..18099 | + | hypothetical protein |
| FN | CDS | NC_021505.1 | 20903..21241 | + | zinc ribbon domain-containing protein |
| FN | CDS | NC_021505.1 | 29460..29879 | - | GtrA family protein |
| FN | CDS | NC_021505.1 | 333665..333862 | + | hypothetical protein |
| FN | CDS | NC_021505.1 | 344909..345034 | + | type II toxin-antitoxin system HigB family toxin |
| FN | CDS | NC_021505.1 | 347261..347839 | + | hypothetical protein |
| FN | CDS | NC_021505.1 | 483140..483211 | - | pyrroloquinoline quinone precursor peptide PqqA |
| FN | CDS | NC_021505.1 | 488761..489042 | - | hypothetical protein |
| FN | CDS | NC_021505.1 | 489540..489701 | + | hypothetical protein |
| FN | CDS | NC_021505.1 | 490974..491216 | + | DUF4160 domain-containing protein |
| FN | CDS | NC_021505.1 | 657757..658194 | + | antitoxin Xre/MbcA/ParS toxin-binding domain-containing protein |
| FN | CDS | NC_021505.1 | 659624..659848 | - | helix-turn-helix domain-containing protein |
| FN | CDS | NC_021505.1 | 702006..702368 | - | hypothetical protein |
| FN | CDS | NC_021505.1 | 702416..702622 | - | hypothetical protein |
| FN | CDS | NC_021505.1 | 769837..770604 | - | methyl-accepting chemotaxis protein |
| FN | CDS | NC_021505.1 | 776910..777107 | - | hypothetical protein |
| FN | CDS | NC_021505.1 | 822983..823846 | - | methyl-accepting chemotaxis protein |
| FN | CDS | NC_021505.1 | 1054460..1055224 | - | methyl-accepting chemotaxis protein |
| FN | CDS | NC_021505.1 | 1107034..1107246 | - | cysteine-rich CWC family protein |
| FN | CDS | NC_021505.1 | 1121505..1121963 | + | protein YgfX |
| FN | CDS | NC_021505.1 | 1308343..1308987 | + | hypothetical protein |
| FN | CDS | NC_021505.1 | 1403753..1404031 | + | hypothetical protein |
| FN | CDS | NC_021505.1 | 1408226..1408459 | + | hypothetical protein |
| FN | CDS | NC_021505.1 | 1408691..1409125 | + | FliH/SctL family protein |
| FN | CDS | NC_021505.1 | 1461733..1461894 | + | hypothetical protein |
| FN | CDS | NC_021505.1 | 1477667..1477921 | + | type II toxin-antitoxin system HicA family toxin |
| FN | CDS | NC_021505.1 | 1579480..1580889 | - | amino acid permease |
| FN | CDS | NC_021505.1 | 1580979..1581233 | - | hypothetical protein |
| FN | CDS | NC_021505.1 | 1738514..1738708 | - | hypothetical protein |
| FN | CDS | NC_021505.1 | 1887143..1887406 | + | ogr/Delta-like zinc finger family protein |
| FN | CDS | NC_021505.1 | 1941868..1942266 | + | type III secretion system chaperone |
| FN | CDS | NC_021505.1 | 1978693..1980846 | - | glycogen debranching protein GlgX |
| FN | CDS | NC_021505.1 | 2051692..2052342 | + | methyl-accepting chemotaxis protein |
| FN | CDS | NC_021505.1 | 2086549..2086770 | + | hypothetical protein |
| FN | CDS | NC_021505.1 | 2090570..2090794 | - | hypothetical protein |
| FN | CDS | NC_021505.1 | 2091737..2092078 | + | hypothetical protein |
| FN | CDS | NC_021505.1 | 2092409..2092888 | + | hypothetical protein |
| FN | CDS | NC_021505.1 | 2097312..2097554 | - | Pathogenicity locus |
| FN | CDS | NC_021505.1 | 2253759..2254148 | + | hypothetical protein |
| FN | CDS | NC_021505.1 | 2286392..2286691 | + | hypothetical protein |
| FN | CDS | NC_021505.1 | 2455375..2455701 | + | winged helix-turn-helix domain-containing protein |
| FN | CDS | NC_021505.1 | 2516256..2516543 | - | DUF4189 domain-containing protein |
| FN | CDS | NC_021505.1 | 2572868..2574136 | - | carbohydrate porin |
| FN | CDS | NC_021505.1 | 2587063..2587818 | - | methyl-accepting chemotaxis protein |
| FN | CDS | NC_021505.1 | 2623632..2623856 | - | DUF7683 domain-containing protein |
| FN | CDS | NC_021505.1 | 2674496..2675056 | + | hypothetical protein |
| FN | CDS | NC_021505.1 | 2801116..2801475 | + | ADP-ribosylglycohydrolase family protein |
| FN | CDS | NC_021505.1 | 2802766..2803287 | + | hypothetical protein |
| FN | CDS | NC_021505.1 | 2970688..2971005 | - | GNAT family N-acetyltransferase |
| FN | CDS | NC_021505.1 | 2974652..2974864 | - | hypothetical protein |
| FN | CDS | NC_021505.1 | 2979306..2979488 | + | DUF2256 domain-containing protein |
| FN | CDS | NC_021505.1 | 3206194..3206688 | - | methyl-accepting chemotaxis protein |
| FN | CDS | NC_021505.1 | 3516873..3517016 | + | hypothetical protein |
| FN | CDS | NC_021505.1 | 3588358..3588429 | - | pyrroloquinoline quinone precursor peptide PqqA |
| FN | CDS | NC_021505.1 | 3645207..3645590 | - | APC family permease |
| FN | CDS | NC_021505.1 | 3653738..3653983 | + | hypothetical protein |
| FN | CDS | NC_021505.1 | 3829786..3830190 | + | hypothetical protein |
| FN | CDS | NC_021505.1 | 4096371..4096757 | - | hypothetical protein |
| FN | CDS | NC_021505.1 | 4142347..4142481 | + | hypothetical protein |
| FN | CDS | NC_021505.1 | 4154132..4154683 | - | DUF1090 family protein |
| FN | CDS | NC_021505.1 | 4238829..4239596 | - | methyl-accepting chemotaxis protein |
| FN | CDS | NC_021505.1 | 4406574..4406801 | - | hypothetical protein |
| FN | CDS | NC_021505.1 | 4498270..4499151 | - | methyl-accepting chemotaxis protein |
| FN | CDS | NC_021505.1 | 4642854..4646333 | + | hybrid sensor histidine kinase/response regulator |
| FN | CDS | NC_021505.1 | 4979432..4980838 | + | NAD(P)/FAD-dependent oxidoreductase |
| FN | CDS | NC_021505.1 | 5110048..5110479 | - | methyl-accepting chemotaxis protein |
| FN | CDS | NC_021505.1 | 5114761..5115063 | + | DUF3565 domain-containing protein |
| FN | CDS | NC_021505.1 | 5233828..5234139 | + | hypothetical protein |
| FN | CDS | NC_021505.1 | 5235958..5236188 | - | hypothetical protein |
| FN | CDS | NC_021505.1 | 5395796..5396197 | - | hypothetical protein |
| FN | CDS | NC_021505.1 | 5402778..5403152 | + | (R)-mandelonitrile lyase |
| FN | CDS | NC_021505.1 | 5431418..5431687 | - | hypothetical protein |
| FN | CDS | NC_021505.1 | 5442530..5443009 | - | hypothetical protein |
| FN | CDS | NC_021505.1 | 5467841..5468236 | - | hypothetical protein |
| FN | CDS | NC_021505.1 | 5469103..5469498 | - | DUF6957 family protein |
| FN | CDS | NC_021505.1 | 5601999..5602883 | - | methyl-accepting chemotaxis protein |
| FN | CDS | NC_021505.1 | 5672247..5672666 | - | hypothetical protein |
| FN | CDS | NC_021505.1 | 5672675..5672989 | - | hypothetical protein |
| FN | CDS | NC_021505.1 | 5758049..5758807 | + | methyl-accepting chemotaxis protein |
| FN | CDS | NC_021505.1 | 5939374..5939829 | + | REP-associated tyrosine transposase |
| FN | CDS | NC_021505.1 | 6011288..6011566 | - | type II toxin-antitoxin system RelE/ParE family toxin |
| FN | CDS | NC_021505.1 | 6034453..6034803 | - | RHS repeat-associated core domain-containing protein |
| FN | CDS | NC_021505.1 | 6155699..6155944 | - | membrane protein insertion efficiency factor YidD |
| FP | CDS | NC_021505.1 | 81549..82484 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 125613..126467 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 180064..180327 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 180352..180819 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 186522..186668 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 213738..214403 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 347588..347839 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 456764..457207 | + | Response regulator protein PilR |
| FP | CDS | NC_021505.1 | 601742..602332 | + | Formate dehydrogenase-O major subunit |
| FP | CDS | NC_021505.1 | 658640..658804 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 705296..705565 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 715252..715350 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 741009..742970 | - | Methyl-accepting chemotaxis protein CtpL |
| FP | CDS | NC_021505.1 | 769837..771780 | - | Methyl-accepting chemotaxis protein McpH |
| FP | CDS | NC_021505.1 | 822983..824902 | - | Methyl-accepting chemotaxis protein McpS |
| FP | CDS | NC_021505.1 | 950118..950207 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 981528..982448 | + | Cell division protein FtsA |
| FP | CDS | NC_021505.1 | 982412..982789 | + | Cell division protein FtsA |
| FP | CDS | NC_021505.1 | 1051215..1051448 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 1054460..1056334 | - | Methyl-accepting chemotaxis protein McpG |
| FP | CDS | NC_021505.1 | 1081662..1081835 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 1155617..1157155 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 1271388..1272137 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 1293060..1293191 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 1296350..1296922 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 1308373..1308681 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 1308712..1308948 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 1315936..1316133 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 1436805..1440158 | + | Linear gramicidin synthase subunit B |
| FP | CDS | NC_021505.1 | 1489722..1489922 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 1580950..1581111 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 1657104..1657202 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 1683265..1683447 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 1923228..1923554 | + | Response regulator GacA |
| FP | CDS | NC_021505.1 | 1923559..1923867 | + | Response regulator GacA |
| FP | CDS | NC_021505.1 | 1939322..1939822 | + | Surface presentation of antigens protein SpaS |
| FP | CDS | NC_021505.1 | 1958442..1958573 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 1995395..1995490 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 2050261..2052342 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 2053347..2053487 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 2087527..2087616 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 2088715..2088936 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 2097642..2098154 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 2211640..2212002 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 2234371..2234643 | - | Transcriptional regulator SlyA |
| FP | CDS | NC_021505.1 | 2470127..2470336 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 2587063..2589207 | - | Methyl-accepting chemotaxis protein TlpQ |
| FP | CDS | NC_021505.1 | 2591859..2592104 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 2624598..2624987 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 2638781..2638897 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 2676068..2676319 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 2680237..2680590 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 2765658..2765777 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 2803544..2803642 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 2815790..2816932 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 2818252..2819424 | + | Double-stranded DNA deaminase toxin A |
| FP | CDS | NC_021505.1 | 2903466..2904620 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 2932339..2932446 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 2945929..2946105 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 2954493..2956055 | - | Cation/acetate symporter ActP |
| FP | CDS | NC_021505.1 | 2987476..2987784 | - | IS3 family transposase ISPsy24 |
| FP | CDS | NC_021505.1 | 2995639..2995782 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 3162432..3162758 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 3173980..3174222 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 3206194..3207513 | - | Biofilm dispersion protein BdlA |
| FP | CDS | NC_021505.1 | 3257315..3257413 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 3332713..3333984 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 3414088..3414246 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 3420611..3421051 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 3424128..3424478 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 3435041..3435205 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 3483868..3484614 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 3487051..3487371 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 3487311..3487604 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 3515175..3515303 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 3559390..3559857 | - | D-serine/D-alanine/glycine transporter |
| FP | CDS | NC_021505.1 | 3559891..3560814 | - | D-serine/D-alanine/glycine transporter |
| FP | CDS | NC_021505.1 | 3584428..3584769 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 3584868..3585329 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 3644767..3645210 | - | 2-succinylbenzoate--CoA ligase |
| FP | CDS | NC_021505.1 | 3652974..3653726 | - | 4-hydroxy-tetrahydrodipicolinate synthase |
| FP | CDS | NC_021505.1 | 3674140..3674805 | + | Hydrogen cyanide synthase subunit HcnB |
| FP | CDS | NC_021505.1 | 3674802..3675497 | + | Hydrogen cyanide synthase subunit HcnB |
| FP | CDS | NC_021505.1 | 3683260..3685575 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 3685928..3686188 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 3686400..3686882 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 3686711..3687328 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 3738942..3739184 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 3774303..3776261 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 3776251..3776877 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 3803439..3803816 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 3851781..3852041 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 3889254..3889907 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 3889962..3893675 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 3916682..3916795 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 4017861..4019426 | - | Methyl-accepting chemotaxis protein Aer |
| FP | CDS | NC_021505.1 | 4142341..4142631 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 4170392..4170742 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 4170696..4170881 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 4238829..4240505 | - | Methyl-accepting chemotaxis protein CtpH |

> Primeros pares con coordenadas no exactas

| tipo | seqid | ref | Δstart | Δend | solape |
|---|---|---|---|---|---|
| CDS | NC_021505.1 | 14880..15563 | +216 | +0 | 0.6842 |
| CDS | NC_021505.1 | 15566..16435 | +84 | +0 | 0.9034 |
| CDS | NC_021505.1 | 30137..31105 | +150 | +0 | 0.8452 |
| CDS | NC_021505.1 | 59393..60178 | +0 | -15 | 0.9809 |
| CDS | NC_021505.1 | 72329..72643 | -3 | +0 | 0.9906 |
| CDS | NC_021505.1 | 76477..78042 | +0 | +27 | 0.9831 |
| CDS | NC_021505.1 | 93346..95433 | +36 | +0 | 0.9828 |
| CDS | NC_021505.1 | 95430..95711 | +15 | +0 | 0.9468 |
| CDS | NC_021505.1 | 117345..118301 | +48 | +0 | 0.9498 |
| CDS | NC_021505.1 | 119858..122605 | -51 | +0 | 0.9818 |
| CDS | NC_021505.1 | 175007..177421 | +0 | +72 | 0.9710 |
| CDS | NC_021505.1 | 182924..184282 | +0 | +18 | 0.9869 |
| CDS | NC_021505.1 | 226648..227286 | +0 | -48 | 0.9249 |
| CDS | NC_021505.1 | 280509..282404 | +30 | +0 | 0.9842 |
| CDS | NC_021505.1 | 315010..315492 | +0 | +57 | 0.8944 |
| CDS | NC_021505.1 | 324068..324898 | +0 | -45 | 0.9458 |
| CDS | NC_021505.1 | 372257..373105 | +0 | -3 | 0.9965 |
| CDS | NC_021505.1 | 399978..400850 | +0 | +753 | 0.5369 |
| CDS | NC_021505.1 | 415236..415997 | +0 | -3 | 0.9961 |
| CDS | NC_021505.1 | 421026..421979 | +6 | +0 | 0.9937 |
| CDS | NC_021505.1 | 448341..449453 | +0 | +39 | 0.9661 |
| CDS | NC_021505.1 | 462819..463583 | +0 | +63 | 0.9239 |
| CDS | NC_021505.1 | 490335..490889 | -33 | +0 | 0.9439 |
| CDS | NC_021505.1 | 513126..514445 | +0 | +6 | 0.9955 |
| CDS | NC_021505.1 | 519858..520856 | +0 | -30 | 0.9700 |
| CDS | NC_021505.1 | 523347..523979 | +3 | +0 | 0.9953 |
| CDS | NC_021505.1 | 555490..556449 | -174 | +0 | 0.8466 |
| CDS | NC_021505.1 | 559885..560418 | +54 | +0 | 0.8989 |
| CDS | NC_021505.1 | 601742..604810 | +639 | +0 | 0.7918 |
| CDS | NC_021505.1 | 640311..642206 | +0 | +102 | 0.9489 |
| CDS | NC_021505.1 | 644613..645839 | +21 | +0 | 0.9829 |
| CDS | NC_021505.1 | 663802..665064 | -36 | +0 | 0.9723 |
| CDS | NC_021505.1 | 666485..667252 | -36 | +0 | 0.9552 |
| CDS | NC_021505.1 | 705971..707173 | -18 | +0 | 0.9853 |
| CDS | NC_021505.1 | 728793..729857 | +0 | +24 | 0.9780 |
| CDS | NC_021505.1 | 781766..782686 | +0 | -3 | 0.9967 |
| CDS | NC_021505.1 | 837135..837842 | +0 | +51 | 0.9328 |
| CDS | NC_021505.1 | 863512..864294 | +0 | +6 | 0.9924 |
| CDS | NC_021505.1 | 876353..876772 | -6 | +0 | 0.9859 |
| CDS | NC_021505.1 | 889273..889818 | +0 | -51 | 0.9066 |
| CDS | NC_021505.1 | 935362..935607 | -3 | +0 | 0.9880 |
| CDS | NC_021505.1 | 954278..955198 | +0 | +78 | 0.9219 |
| CDS | NC_021505.1 | 958379..959590 | -3 | +0 | 0.9975 |
| CDS | NC_021505.1 | 967254..967709 | +51 | +0 | 0.8882 |
| CDS | NC_021505.1 | 968650..968943 | +45 | +0 | 0.8469 |
| CDS | NC_021505.1 | 968943..970682 | -3 | +0 | 0.9983 |
| CDS | NC_021505.1 | 985118..985573 | +0 | +261 | 0.6360 |
| CDS | NC_021505.1 | 1025956..1026909 | +0 | -6 | 0.9937 |
| CDS | NC_021505.1 | 1053100..1054221 | +0 | -9 | 0.9920 |
| CDS | NC_021505.1 | 1060660..1061469 | -66 | +0 | 0.9247 |
| CDS | NC_021505.1 | 1092383..1093702 | +0 | +96 | 0.9322 |
| CDS | NC_021505.1 | 1097167..1099458 | +0 | -3 | 0.9987 |
| CDS | NC_021505.1 | 1126525..1127925 | -33 | +0 | 0.9770 |
| CDS | NC_021505.1 | 1133415..1134137 | -18 | +0 | 0.9757 |
| CDS | NC_021505.1 | 1175885..1177414 | +0 | +42 | 0.9733 |
| CDS | NC_021505.1 | 1188638..1189285 | +0 | +27 | 0.9600 |
| CDS | NC_021505.1 | 1201632..1202585 | +6 | +0 | 0.9937 |
| CDS | NC_021505.1 | 1219842..1220318 | +0 | +3 | 0.9938 |
| CDS | NC_021505.1 | 1234281..1235888 | +0 | +24 | 0.9853 |
| CDS | NC_021505.1 | 1260499..1261572 | +0 | +3 | 0.9972 |
| CDS | NC_021505.1 | 1261703..1262794 | +33 | +0 | 0.9698 |
| CDS | NC_021505.1 | 1309045..1309839 | -51 | +0 | 0.9397 |
| CDS | NC_021505.1 | 1320969..1322027 | -6 | +0 | 0.9944 |
| CDS | NC_021505.1 | 1329479..1329700 | +42 | +0 | 0.8108 |
| CDS | NC_021505.1 | 1339634..1340440 | +180 | +0 | 0.7770 |
| CDS | NC_021505.1 | 1435648..1436592 | +6 | +0 | 0.9937 |
| CDS | NC_021505.1 | 1477918..1478250 | +9 | +0 | 0.9730 |
| CDS | NC_021505.1 | 1496842..1497969 | -39 | +0 | 0.9666 |
| CDS | NC_021505.1 | 1509108..1509410 | +0 | -96 | 0.6832 |
| CDS | NC_021505.1 | 1521625..1522038 | +0 | -48 | 0.8841 |
| CDS | NC_021505.1 | 1531447..1532508 | +81 | +0 | 0.9237 |
| CDS | NC_021505.1 | 1534707..1536047 | +0 | -6 | 0.9955 |
| CDS | NC_021505.1 | 1537748..1538341 | +0 | +39 | 0.9384 |
| CDS | NC_021505.1 | 1539283..1539912 | -36 | +0 | 0.9459 |
| CDS | NC_021505.1 | 1540038..1540763 | -129 | +0 | 0.8491 |
| CDS | NC_021505.1 | 1543505..1544923 | +15 | +0 | 0.9894 |
| CDS | NC_021505.1 | 1545105..1545989 | +0 | -27 | 0.9695 |
| CDS | NC_021505.1 | 1578049..1578375 | +0 | -9 | 0.9725 |
| CDS | NC_021505.1 | 1588884..1590122 | -81 | +0 | 0.9386 |
| CDS | NC_021505.1 | 1619147..1620673 | -3 | +0 | 0.9980 |
| CDS | NC_021505.1 | 1620675..1622504 | +93 | +0 | 0.9492 |
| CDS | NC_021505.1 | 1624121..1625200 | +36 | +0 | 0.9667 |
| CDS | NC_021505.1 | 1627787..1628728 | -6 | +0 | 0.9937 |
| CDS | NC_021505.1 | 1709527..1710312 | +0 | +69 | 0.9193 |
| CDS | NC_021505.1 | 1721962..1724361 | +0 | -36 | 0.9850 |
| CDS | NC_021505.1 | 1739201..1739944 | +0 | +24 | 0.9688 |
| CDS | NC_021505.1 | 1753247..1753438 | +0 | +9 | 0.9552 |
| CDS | NC_021505.1 | 1795414..1796706 | +0 | +21 | 0.9840 |
| CDS | NC_021505.1 | 1807714..1809396 | +0 | -18 | 0.9893 |
| CDS | NC_021505.1 | 1816361..1816966 | -3 | +0 | 0.9951 |
| CDS | NC_021505.1 | 1861660..1862895 | +12 | +0 | 0.9903 |
| CDS | NC_021505.1 | 1885612..1885860 | -66 | +0 | 0.7905 |
| CDS | NC_021505.1 | 1905917..1906681 | +0 | -30 | 0.9608 |
| CDS | NC_021505.1 | 1928607..1929131 | -246 | +0 | 0.6809 |
| CDS | NC_021505.1 | 1929238..1930797 | -102 | +0 | 0.9386 |
| CDS | NC_021505.1 | 1937608..1938270 | +15 | +0 | 0.9774 |
| CDS | NC_021505.1 | 1940427..1941146 | +54 | +0 | 0.9250 |
| CDS | NC_021505.1 | 1957845..1958501 | -9 | +0 | 0.9865 |
| CDS | NC_021505.1 | 1980999..1981283 | +0 | -3 | 0.9895 |
| CDS | NC_021505.1 | 2017226..2017654 | +0 | +54 | 0.8882 |

- Aviso: 1 feature(s) GFF multi-segmento fusionados por ID; se comparan por su intervalo envolvente.

## DFAST

- Estado: **ok**.
- Fuente: `executed`.
- Versión: `DFAST ver. 1.6.2` (estado `ok`).
- Ejecución: **ok**, 398.80 s, RSS máx. 3958.8 MiB (`2026-11-25T22:46:05Z` → `2026-11-25T22:52:44Z`).
- Salida: `/tmp/external_benchmark_dfast_csjormus/out/genome.gff` (10798204 bytes, SHA256 `fae16713b9fb11cb…`).
- Comando: `/usr/local/bin/dfast --genome /tmp/benchmark_genome_b9aq49k3.fna --out /tmp/external_benchmark_dfast_csjormus/out --cpu 4 --force --dbroot /data/dfast_db`.
- Anotación comparada: `reports/external_benchmark/multi_v2/pputida/evidence/dfast/genome.gff` (formato `gff3`, SHA256 `fae16713b9fb11cb…`).

### Features (DFAST)

| tipo | ref | query | TP | FP | FN | precisión | recall | F1 | exactos | borde distinto | solape medio |
|---|---|---|---|---|---|---|---|---|---|---|---|
| CDS | 5377 | 5431 | 5299 | 132 | 78 | 0.9757 | 0.9855 | 0.9806 | 4977 | 322 | 0.9960 |
| tRNA | 75 | 80 | 75 | 5 | 0 | 0.9375 | 1.0000 | 0.9677 | 75 | 0 | 1.0000 |
| rRNA | 22 | 22 | 22 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 0 | 22 | 0.9830 |

### Atributos (DFAST)

| tipo | pares | product en ambos | product literal | product tokens | gene literal | locus literal |
|---|---|---|---|---|---|---|
| CDS | 5299 | 5299 | 0.7420 | 0.8116 | 0.9843 | 0.0000 |
| tRNA | 75 | 75 | 0.9867 | 0.9911 | — | 0.0000 |
| rRNA | 22 | 22 | 1.0000 | 1.0000 | — | 0.0000 |

### Hipotéticas (DFAST)

| tipo | ref CDS | ref con product | ref hyp | % ref | query CDS | query con product | query hyp | % query |
|---|---|---|---|---|---|---|---|---|
| CDS | 5377 | 5377 | 677 | 12.6 | 5431 | 5431 | 1006 | 18.5 |
| tRNA | 0 | 0 | 0 | — | 0 | 0 | 0 | — |
| rRNA | 0 | 0 | 0 | — | 0 | 0 | 0 | — |

#### Diferencias

- Falsos negativos (referencia sin pareja): **78**.
- Falsos positivos (herramienta sin pareja): **137**.
- Pares solo con hebra distinta: **1**.
- Pares emparejados con algún borde distinto: **344**.

Las listas de diferencias están truncadas a 100 elementos por lista; los conteos son completos.

> Primeras diferencias no emparejadas

| tipo | feature | seqid | coordenadas | hebra | product/locus |
|---|---|---|---|---|---|
| FN | CDS | NC_021505.1 | 17890..18099 | + | hypothetical protein |
| FN | CDS | NC_021505.1 | 20903..21241 | + | zinc ribbon domain-containing protein |
| FN | CDS | NC_021505.1 | 29460..29879 | - | GtrA family protein |
| FN | CDS | NC_021505.1 | 333665..333862 | + | hypothetical protein |
| FN | CDS | NC_021505.1 | 344909..345034 | + | type II toxin-antitoxin system HigB family toxin |
| FN | CDS | NC_021505.1 | 347261..347839 | + | hypothetical protein |
| FN | CDS | NC_021505.1 | 483140..483211 | - | pyrroloquinoline quinone precursor peptide PqqA |
| FN | CDS | NC_021505.1 | 488761..489042 | - | hypothetical protein |
| FN | CDS | NC_021505.1 | 489540..489701 | + | hypothetical protein |
| FN | CDS | NC_021505.1 | 490974..491216 | + | DUF4160 domain-containing protein |
| FN | CDS | NC_021505.1 | 657757..658194 | + | antitoxin Xre/MbcA/ParS toxin-binding domain-containing protein |
| FN | CDS | NC_021505.1 | 659624..659848 | - | helix-turn-helix domain-containing protein |
| FN | CDS | NC_021505.1 | 702006..702368 | - | hypothetical protein |
| FN | CDS | NC_021505.1 | 702416..702622 | - | hypothetical protein |
| FN | CDS | NC_021505.1 | 769837..770604 | - | methyl-accepting chemotaxis protein |
| FN | CDS | NC_021505.1 | 776910..777107 | - | hypothetical protein |
| FN | CDS | NC_021505.1 | 822983..823846 | - | methyl-accepting chemotaxis protein |
| FN | CDS | NC_021505.1 | 1054460..1055224 | - | methyl-accepting chemotaxis protein |
| FN | CDS | NC_021505.1 | 1107034..1107246 | - | cysteine-rich CWC family protein |
| FN | CDS | NC_021505.1 | 1121505..1121963 | + | protein YgfX |
| FN | CDS | NC_021505.1 | 1308343..1308987 | + | hypothetical protein |
| FN | CDS | NC_021505.1 | 1403753..1404031 | + | hypothetical protein |
| FN | CDS | NC_021505.1 | 1408226..1408459 | + | hypothetical protein |
| FN | CDS | NC_021505.1 | 1408691..1409125 | + | FliH/SctL family protein |
| FN | CDS | NC_021505.1 | 1461733..1461894 | + | hypothetical protein |
| FN | CDS | NC_021505.1 | 1477667..1477921 | + | type II toxin-antitoxin system HicA family toxin |
| FN | CDS | NC_021505.1 | 1580979..1581233 | - | hypothetical protein |
| FN | CDS | NC_021505.1 | 1738514..1738708 | - | hypothetical protein |
| FN | CDS | NC_021505.1 | 1887143..1887406 | + | ogr/Delta-like zinc finger family protein |
| FN | CDS | NC_021505.1 | 1941868..1942266 | + | type III secretion system chaperone |
| FN | CDS | NC_021505.1 | 2051692..2052342 | + | methyl-accepting chemotaxis protein |
| FN | CDS | NC_021505.1 | 2086549..2086770 | + | hypothetical protein |
| FN | CDS | NC_021505.1 | 2090570..2090794 | - | hypothetical protein |
| FN | CDS | NC_021505.1 | 2091737..2092078 | + | hypothetical protein |
| FN | CDS | NC_021505.1 | 2092409..2092888 | + | hypothetical protein |
| FN | CDS | NC_021505.1 | 2097312..2097554 | - | Pathogenicity locus |
| FN | CDS | NC_021505.1 | 2253759..2254148 | + | hypothetical protein |
| FN | CDS | NC_021505.1 | 2286392..2286691 | + | hypothetical protein |
| FN | CDS | NC_021505.1 | 2455375..2455701 | + | winged helix-turn-helix domain-containing protein |
| FN | CDS | NC_021505.1 | 2516256..2516543 | - | DUF4189 domain-containing protein |
| FN | CDS | NC_021505.1 | 2587063..2587818 | - | methyl-accepting chemotaxis protein |
| FN | CDS | NC_021505.1 | 2623632..2623856 | - | DUF7683 domain-containing protein |
| FN | CDS | NC_021505.1 | 2674496..2675056 | + | hypothetical protein |
| FN | CDS | NC_021505.1 | 2801116..2801475 | + | ADP-ribosylglycohydrolase family protein |
| FN | CDS | NC_021505.1 | 2802766..2803287 | + | hypothetical protein |
| FN | CDS | NC_021505.1 | 2970688..2971005 | - | GNAT family N-acetyltransferase |
| FN | CDS | NC_021505.1 | 2974652..2974864 | - | hypothetical protein |
| FN | CDS | NC_021505.1 | 2979306..2979488 | + | DUF2256 domain-containing protein |
| FN | CDS | NC_021505.1 | 3206194..3206688 | - | methyl-accepting chemotaxis protein |
| FN | CDS | NC_021505.1 | 3516873..3517016 | + | hypothetical protein |
| FN | CDS | NC_021505.1 | 3588358..3588429 | - | pyrroloquinoline quinone precursor peptide PqqA |
| FN | CDS | NC_021505.1 | 3645207..3645590 | - | APC family permease |
| FN | CDS | NC_021505.1 | 3653738..3653983 | + | hypothetical protein |
| FN | CDS | NC_021505.1 | 3829786..3830190 | + | hypothetical protein |
| FN | CDS | NC_021505.1 | 4096371..4096757 | - | hypothetical protein |
| FN | CDS | NC_021505.1 | 4142347..4142481 | + | hypothetical protein |
| FN | CDS | NC_021505.1 | 4154132..4154683 | - | DUF1090 family protein |
| FN | CDS | NC_021505.1 | 4238829..4239596 | - | methyl-accepting chemotaxis protein |
| FN | CDS | NC_021505.1 | 4406574..4406801 | - | hypothetical protein |
| FN | CDS | NC_021505.1 | 4498270..4499151 | - | methyl-accepting chemotaxis protein |
| FN | CDS | NC_021505.1 | 5110048..5110479 | - | methyl-accepting chemotaxis protein |
| FN | CDS | NC_021505.1 | 5114761..5115063 | + | DUF3565 domain-containing protein |
| FN | CDS | NC_021505.1 | 5233828..5234139 | + | hypothetical protein |
| FN | CDS | NC_021505.1 | 5235958..5236188 | - | hypothetical protein |
| FN | CDS | NC_021505.1 | 5395796..5396197 | - | hypothetical protein |
| FN | CDS | NC_021505.1 | 5402778..5403152 | + | (R)-mandelonitrile lyase |
| FN | CDS | NC_021505.1 | 5431418..5431687 | - | hypothetical protein |
| FN | CDS | NC_021505.1 | 5442530..5443009 | - | hypothetical protein |
| FN | CDS | NC_021505.1 | 5467841..5468236 | - | hypothetical protein |
| FN | CDS | NC_021505.1 | 5469103..5469498 | - | DUF6957 family protein |
| FN | CDS | NC_021505.1 | 5601999..5602883 | - | methyl-accepting chemotaxis protein |
| FN | CDS | NC_021505.1 | 5672247..5672666 | - | hypothetical protein |
| FN | CDS | NC_021505.1 | 5672675..5672989 | - | hypothetical protein |
| FN | CDS | NC_021505.1 | 5758049..5758807 | + | methyl-accepting chemotaxis protein |
| FN | CDS | NC_021505.1 | 5939374..5939829 | + | REP-associated tyrosine transposase |
| FN | CDS | NC_021505.1 | 6011288..6011566 | - | type II toxin-antitoxin system RelE/ParE family toxin |
| FN | CDS | NC_021505.1 | 6034453..6034803 | - | RHS repeat-associated core domain-containing protein |
| FN | CDS | NC_021505.1 | 6155699..6155944 | - | membrane protein insertion efficiency factor YidD |
| FP | CDS | NC_021505.1 | 81549..82484 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 125613..126467 | + | endonuclease/exonuclease/phosphatase family protein |
| FP | CDS | NC_021505.1 | 180064..180327 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 180352..180819 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 186522..186668 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 213738..214403 | + | RHS repeat-associated core domain-containing protein |
| FP | CDS | NC_021505.1 | 347588..347839 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 456764..457207 | + | response regulator |
| FP | CDS | NC_021505.1 | 658640..658804 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 705296..705565 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 715252..715350 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 741009..742970 | - | methyl-accepting chemotaxis protein CtpL |
| FP | CDS | NC_021505.1 | 769837..771780 | - | methyl-accepting chemotaxis protein |
| FP | CDS | NC_021505.1 | 822983..824902 | - | methyl-accepting chemotaxis protein |
| FP | CDS | NC_021505.1 | 950118..950207 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 981528..982448 | + | cell division protein FtsA |
| FP | CDS | NC_021505.1 | 982412..982789 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 1051215..1051448 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 1054460..1056334 | - | methyl-accepting chemotaxis protein PctA |
| FP | CDS | NC_021505.1 | 1081662..1081835 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 1155617..1157155 | - | FMN-binding glutamate synthase family protein |
| FP | CDS | NC_021505.1 | 1271388..1272137 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 1293060..1293191 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 1296350..1296922 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 1308373..1308681 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 1308712..1308948 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 1315936..1316133 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 1436805..1440158 | + | amino acid adenylation domain-containing protein |
| FP | CDS | NC_021505.1 | 1489722..1489922 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 1580950..1581111 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 1657104..1657202 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 1683265..1683447 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 1923228..1923554 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 1923559..1923867 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 1939322..1939822 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 1958442..1958573 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 1995395..1995490 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 2050261..2052342 | + | methyl-accepting chemotaxis protein |
| FP | CDS | NC_021505.1 | 2053347..2053487 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 2087527..2087616 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 2088715..2088936 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 2097642..2098154 | + | putative metallopeptidase |
| FP | CDS | NC_021505.1 | 2211640..2212002 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 2234371..2234643 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 2470127..2470336 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 2587063..2589207 | - | methyl-accepting chemotaxis protein |
| FP | CDS | NC_021505.1 | 2591859..2592104 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 2624598..2624987 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 2638781..2638897 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 2676068..2676319 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 2680237..2680590 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 2765658..2765777 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 2803544..2803642 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 2815790..2816932 | + | RHS repeat-associated core domain-containing protein |
| FP | CDS | NC_021505.1 | 2818252..2819424 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 2903466..2904620 | + | group II intron reverse transcriptase/maturase |
| FP | CDS | NC_021505.1 | 2932339..2932446 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 2945929..2946105 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 2954493..2956055 | - | cation/acetate symporter ActP |
| FP | CDS | NC_021505.1 | 2987476..2987784 | - | IS3 family transposase |
| FP | CDS | NC_021505.1 | 2995639..2995782 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 3162432..3162758 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 3173980..3174222 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 3206194..3207513 | - | methyl-accepting chemotaxis protein |
| FP | CDS | NC_021505.1 | 3257315..3257413 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 3332713..3333984 | + | MFS transporter |
| FP | CDS | NC_021505.1 | 3414088..3414246 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 3420611..3421051 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 3424128..3424478 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 3435041..3435205 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 3483868..3484614 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 3487051..3487371 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 3487311..3487604 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 3515175..3515303 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 3559390..3559857 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 3559891..3560814 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 3584428..3584769 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 3584868..3585329 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 3644767..3645210 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 3652974..3653726 | - | 4-hydroxy-tetrahydrodipicolinate synthase |
| FP | CDS | NC_021505.1 | 3674140..3674805 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 3674802..3675497 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 3683260..3685575 | - | type VI secretion system membrane subunit TssM |
| FP | CDS | NC_021505.1 | 3685928..3686188 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 3686400..3686882 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 3686711..3687328 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 3738942..3739184 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 3774303..3776261 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 3776251..3776877 | + | hypothetical protein |
| FP | CDS | NC_021505.1 | 3803439..3803816 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 3851781..3852041 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 3889254..3889907 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 3889962..3893675 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 3916682..3916795 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 4017861..4019426 | - | PAS domain-containing methyl-accepting chemotaxis protein |
| FP | CDS | NC_021505.1 | 4142341..4142631 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 4170392..4170742 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 4170696..4170881 | - | hypothetical protein |
| FP | CDS | NC_021505.1 | 4238829..4240505 | - | methyl-accepting chemotaxis protein |
| FP | CDS | NC_021505.1 | 4353945..4354175 | - | hypothetical protein |

> Primeros pares con coordenadas no exactas

| tipo | seqid | ref | Δstart | Δend | solape |
|---|---|---|---|---|---|
| CDS | NC_021505.1 | 14880..15563 | +216 | +0 | 0.6842 |
| CDS | NC_021505.1 | 15566..16435 | +84 | +0 | 0.9034 |
| CDS | NC_021505.1 | 30137..31105 | +150 | +0 | 0.8452 |
| CDS | NC_021505.1 | 59393..60178 | +0 | -15 | 0.9809 |
| CDS | NC_021505.1 | 72329..72643 | -3 | +0 | 0.9906 |
| CDS | NC_021505.1 | 76477..78042 | +0 | +27 | 0.9831 |
| CDS | NC_021505.1 | 93346..95433 | +36 | +0 | 0.9828 |
| CDS | NC_021505.1 | 95430..95711 | +15 | +0 | 0.9468 |
| CDS | NC_021505.1 | 117345..118301 | +48 | +0 | 0.9498 |
| CDS | NC_021505.1 | 119858..122605 | -51 | +0 | 0.9818 |
| CDS | NC_021505.1 | 175007..177421 | +0 | +72 | 0.9710 |
| CDS | NC_021505.1 | 182924..184282 | +0 | +18 | 0.9869 |
| CDS | NC_021505.1 | 226648..227286 | +0 | -48 | 0.9249 |
| CDS | NC_021505.1 | 280509..282404 | +30 | +0 | 0.9842 |
| CDS | NC_021505.1 | 315010..315492 | +0 | +57 | 0.8944 |
| CDS | NC_021505.1 | 324068..324898 | +0 | -45 | 0.9458 |
| CDS | NC_021505.1 | 372257..373105 | +0 | -3 | 0.9965 |
| CDS | NC_021505.1 | 399978..400850 | +0 | +753 | 0.5369 |
| CDS | NC_021505.1 | 415236..415997 | +0 | -3 | 0.9961 |
| CDS | NC_021505.1 | 421026..421979 | +6 | +0 | 0.9937 |
| CDS | NC_021505.1 | 448341..449453 | +0 | +39 | 0.9661 |
| CDS | NC_021505.1 | 462819..463583 | +0 | +63 | 0.9239 |
| CDS | NC_021505.1 | 490335..490889 | -33 | +0 | 0.9439 |
| CDS | NC_021505.1 | 513126..514445 | +0 | +6 | 0.9955 |
| CDS | NC_021505.1 | 519858..520856 | +0 | -30 | 0.9700 |
| CDS | NC_021505.1 | 523347..523979 | +3 | +0 | 0.9953 |
| CDS | NC_021505.1 | 555490..556449 | -174 | +0 | 0.8466 |
| CDS | NC_021505.1 | 559885..560418 | +54 | +0 | 0.8989 |
| CDS | NC_021505.1 | 640311..642206 | +0 | +102 | 0.9489 |
| CDS | NC_021505.1 | 644613..645839 | +21 | +0 | 0.9829 |
| CDS | NC_021505.1 | 663802..665064 | -36 | +0 | 0.9723 |
| CDS | NC_021505.1 | 666485..667252 | -36 | +0 | 0.9552 |
| CDS | NC_021505.1 | 705971..707173 | -18 | +0 | 0.9853 |
| CDS | NC_021505.1 | 728793..729857 | +0 | +24 | 0.9780 |
| CDS | NC_021505.1 | 781766..782686 | +0 | -3 | 0.9967 |
| CDS | NC_021505.1 | 837135..837842 | +0 | +51 | 0.9328 |
| CDS | NC_021505.1 | 863512..864294 | +0 | +6 | 0.9924 |
| CDS | NC_021505.1 | 876353..876772 | -6 | +0 | 0.9859 |
| CDS | NC_021505.1 | 889273..889818 | +0 | -51 | 0.9066 |
| CDS | NC_021505.1 | 935362..935607 | -3 | +0 | 0.9880 |
| CDS | NC_021505.1 | 954278..955198 | +0 | +78 | 0.9219 |
| CDS | NC_021505.1 | 958379..959590 | -3 | +0 | 0.9975 |
| CDS | NC_021505.1 | 967254..967709 | +51 | +0 | 0.8882 |
| CDS | NC_021505.1 | 968650..968943 | +45 | +0 | 0.8469 |
| CDS | NC_021505.1 | 968943..970682 | -3 | +0 | 0.9983 |
| CDS | NC_021505.1 | 985118..985573 | +0 | +261 | 0.6360 |
| CDS | NC_021505.1 | 1025956..1026909 | +0 | -6 | 0.9937 |
| CDS | NC_021505.1 | 1053100..1054221 | +0 | -9 | 0.9920 |
| CDS | NC_021505.1 | 1060660..1061469 | -66 | +0 | 0.9247 |
| CDS | NC_021505.1 | 1092383..1093702 | +0 | +96 | 0.9322 |
| CDS | NC_021505.1 | 1097167..1099458 | +0 | -3 | 0.9987 |
| CDS | NC_021505.1 | 1126525..1127925 | -33 | +0 | 0.9770 |
| CDS | NC_021505.1 | 1133415..1134137 | -18 | +0 | 0.9757 |
| CDS | NC_021505.1 | 1175885..1177414 | +0 | +42 | 0.9733 |
| CDS | NC_021505.1 | 1188638..1189285 | +0 | +27 | 0.9600 |
| CDS | NC_021505.1 | 1201632..1202585 | +6 | +0 | 0.9937 |
| CDS | NC_021505.1 | 1219842..1220318 | +0 | +3 | 0.9938 |
| CDS | NC_021505.1 | 1234281..1235888 | +0 | +24 | 0.9853 |
| CDS | NC_021505.1 | 1260499..1261572 | +0 | +3 | 0.9972 |
| CDS | NC_021505.1 | 1261703..1262794 | +33 | +0 | 0.9698 |
| CDS | NC_021505.1 | 1309045..1309839 | -51 | +0 | 0.9397 |
| CDS | NC_021505.1 | 1320969..1322027 | -6 | +0 | 0.9944 |
| CDS | NC_021505.1 | 1329479..1329700 | +42 | +0 | 0.8108 |
| CDS | NC_021505.1 | 1339634..1340440 | +180 | +0 | 0.7770 |
| CDS | NC_021505.1 | 1435648..1436592 | +6 | +0 | 0.9937 |
| CDS | NC_021505.1 | 1477918..1478250 | +9 | +0 | 0.9730 |
| CDS | NC_021505.1 | 1496842..1497969 | -39 | +0 | 0.9666 |
| CDS | NC_021505.1 | 1509108..1509410 | +0 | -96 | 0.6832 |
| CDS | NC_021505.1 | 1521625..1522038 | +0 | -48 | 0.8841 |
| CDS | NC_021505.1 | 1531447..1532508 | +81 | +0 | 0.9237 |
| CDS | NC_021505.1 | 1534707..1536047 | +0 | -6 | 0.9955 |
| CDS | NC_021505.1 | 1537748..1538341 | +0 | +39 | 0.9384 |
| CDS | NC_021505.1 | 1539283..1539912 | -36 | +0 | 0.9459 |
| CDS | NC_021505.1 | 1540038..1540763 | -129 | +0 | 0.8491 |
| CDS | NC_021505.1 | 1543505..1544923 | +15 | +0 | 0.9894 |
| CDS | NC_021505.1 | 1545105..1545989 | +0 | -27 | 0.9695 |
| CDS | NC_021505.1 | 1578049..1578375 | +0 | -9 | 0.9725 |
| CDS | NC_021505.1 | 1588884..1590122 | -81 | +0 | 0.9386 |
| CDS | NC_021505.1 | 1619147..1620673 | -3 | +0 | 0.9980 |
| CDS | NC_021505.1 | 1620675..1622504 | +93 | +0 | 0.9492 |
| CDS | NC_021505.1 | 1624121..1625200 | +36 | +0 | 0.9667 |
| CDS | NC_021505.1 | 1627787..1628728 | -6 | +0 | 0.9937 |
| CDS | NC_021505.1 | 1709527..1710312 | +0 | +69 | 0.9193 |
| CDS | NC_021505.1 | 1739201..1739944 | +0 | +24 | 0.9688 |
| CDS | NC_021505.1 | 1753247..1753438 | +0 | +9 | 0.9552 |
| CDS | NC_021505.1 | 1795414..1796706 | +0 | +21 | 0.9840 |
| CDS | NC_021505.1 | 1807714..1809396 | +0 | -18 | 0.9893 |
| CDS | NC_021505.1 | 1816361..1816966 | -3 | +0 | 0.9951 |
| CDS | NC_021505.1 | 1861660..1862895 | +12 | +0 | 0.9903 |
| CDS | NC_021505.1 | 1885612..1885860 | -66 | +0 | 0.7905 |
| CDS | NC_021505.1 | 1905917..1906681 | +0 | -30 | 0.9608 |
| CDS | NC_021505.1 | 1928607..1929131 | -246 | +0 | 0.6809 |
| CDS | NC_021505.1 | 1929238..1930797 | -102 | +0 | 0.9386 |
| CDS | NC_021505.1 | 1937608..1938270 | +15 | +0 | 0.9774 |
| CDS | NC_021505.1 | 1940427..1941146 | +54 | +0 | 0.9250 |
| CDS | NC_021505.1 | 1957845..1958501 | -9 | +0 | 0.9865 |
| CDS | NC_021505.1 | 1980999..1981283 | +0 | -3 | 0.9895 |
| CDS | NC_021505.1 | 2017226..2017654 | +0 | +54 | 0.8882 |
| CDS | NC_021505.1 | 2041217..2042218 | +0 | +15 | 0.9853 |
| CDS | NC_021505.1 | 2072492..2073916 | +0 | -15 | 0.9895 |

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
/opt/conda/bin/python benchmarks/external/external_benchmark.py --genome data/reference_genomes/pputida/pputida_genomic.fna.gz --reference data/reference_genomes/pputida/pputida_genomic.gff.gz --tools bakta,prokka,dfast --run --threads 1 --json reports/external_benchmark/benchmark.json --md reports/external_benchmark/benchmark.md --bakta-db /data/bakta_db
```

Para comparar salidas producidas en otro host, usar `--query-gff TOOL=RUTA` sin `--run`; el origen queda como `precomputed`.
