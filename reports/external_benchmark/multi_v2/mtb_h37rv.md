# Benchmark externo CPU: Bakta, Prokka y DFAST

**Generado:** 2026-11-25T22:16:38Z

**Esquema:** `external-benchmark-run/v1`

**Alcance:** Comparación local de CDS/tRNA/rRNA de Bakta, Prokka y DFAST contra una referencia GFF3/GBFF. CPU, sin descargas de datos.

**Comando reproducible:**

```bash
/opt/conda/bin/python /opt/benchmark/external/external_benchmark.py --genome /work/data/reference_genomes/mtb_h37rv/mtb_h37rv_genomic.fna.gz --reference /work/data/reference_genomes/mtb_h37rv/mtb_h37rv_genomic.gff.gz --tools all --run --threads 4 --bakta-db /data/bakta_db --extra-arg dfast=--dbroot --extra-arg dfast=/data/dfast_db --save-gffs /work/reports/external_benchmark/multi_v2/mtb_h37rv/evidence --json /work/reports/external_benchmark/multi_v2/mtb_h37rv.json --md /work/reports/external_benchmark/multi_v2/mtb_h37rv.md
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

- FASTA: `data/reference_genomes/mtb_h37rv/mtb_h37rv_genomic.fna.gz` (SHA256 `35911742cea55cd5…`, 1 contig, 4411532 bp).
- Referencia: `data/reference_genomes/mtb_h37rv/mtb_h37rv_genomic.gff.gz` (formato `gff3`, SHA256 `fbd1a13d212e1d04…`).
- Conteos de referencia: CDS **3906**, tRNA **45**, rRNA **3**.
- CDS pseudo excluidos en referencia: **0**.
- SHA256 FASTA: `35911742cea55cd556d54948f56e3efa757e40d44e18ebc8742cacc59b2a976e`.
- SHA256 referencia: `fbd1a13d212e1d042059de7e44dfe5d7dbd96bcc71a53f2e5040ea0e50a4a2ed`.

## Validación

- Inicio: `2026-11-25T21:30:55Z`; fin: `2026-11-25T22:16:38Z`; total: **2742.53 s**.
- JSON declarado: `reports/external_benchmark/multi_v2/mtb_h37rv.json`.
- Markdown declarado: `reports/external_benchmark/multi_v2/mtb_h37rv.md`.

| check | estado | detalle |
|---|---|---|
| genome_path | ok | FASTA existe y es legible |
| genome_sha256 | ok | SHA256 35911742cea55cd556d54948f56e3efa757e40d44e18ebc8742cacc59b2a976e |
| reference_path | ok | referencia gff3 con 3954 features |
| reference_sha256 | ok | SHA256 fbd1a13d212e1d042059de7e44dfe5d7dbd96bcc71a53f2e5040ea0e50a4a2ed |
| reference_features | ok | CDS=3906 tRNA=45 rRNA=3 pseudo_excluidos=0 |
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
| Bakta | executed | ok | 0.9479 | 1.0000 | 1.0000 | 0.9287 | 0.9691 | 0.9485 | 1200.0 |
| Prokka | executed | ok | 0.9465 | 0.9184 | 1.0000 | 0.9265 | 0.9666 | 0.9462 | 960.6 |
| DFAST | executed | ok | 0.9471 | 0.9677 | 1.0000 | 0.9278 | 0.9679 | 0.9474 | 521.2 |

La precisión/recall/F1 se calculan contra la referencia local. Si una herramienta no está instalada, su estado es `not_available` y no hay métricas que mostrar.

## Bakta

- Estado: **ok**.
- Fuente: `executed`.
- Bakta DB: `ok` (válida: **sí**, ruta `/data/bakta_db`, versión `0`).
- Versión: `bakta 1.12.2` (estado `ok`).
- Ejecución: **ok**, 1200.00 s, RSS máx. 2654.8 MiB (`2026-11-25T21:31:03Z` → `2026-11-25T21:51:03Z`).
- Salida: `/tmp/external_benchmark_bakta_kkpgkd7g/out/benchmark.gff3` (5361711 bytes, SHA256 `5dcbb66831cd9da8…`).
- Comando: `/usr/local/bin/bakta --output /tmp/external_benchmark_bakta_kkpgkd7g/out --prefix benchmark --threads 4 --force --skip-plot --db /data/bakta_db /tmp/benchmark_genome_azi2xbyp.fna`.
- Anotación comparada: `reports/external_benchmark/multi_v2/mtb_h37rv/evidence/bakta/benchmark.gff3` (formato `gff3`, SHA256 `5dcbb66831cd9da8…`).

### Features (Bakta)

| tipo | ref | query | TP | FP | FN | precisión | recall | F1 | exactos | borde distinto | solape medio |
|---|---|---|---|---|---|---|---|---|---|---|---|
| CDS | 3906 | 4078 | 3784 | 294 | 122 | 0.9279 | 0.9688 | 0.9479 | 2927 | 857 | 0.9819 |
| tRNA | 45 | 45 | 45 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 32 | 13 | 0.9890 |
| rRNA | 3 | 3 | 3 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 1 | 2 | 0.9931 |

### Atributos (Bakta)

| tipo | pares | product en ambos | product literal | product tokens | gene literal | locus literal |
|---|---|---|---|---|---|---|
| CDS | 3784 | 3784 | 0.2315 | 0.4416 | 0.8674 | 0.0000 |
| tRNA | 45 | 45 | 0.0000 | 0.6481 | 0.0000 | 0.0000 |
| rRNA | 3 | 3 | 1.0000 | 1.0000 | 1.0000 | 0.0000 |

### Hipotéticas (Bakta)

| tipo | ref CDS | ref con product | ref hyp | % ref | query CDS | query con product | query hyp | % query |
|---|---|---|---|---|---|---|---|---|
| CDS | 3906 | 3906 | 1052 | 26.9 | 4078 | 4078 | 934 | 22.9 |
| tRNA | 0 | 0 | 0 | — | 0 | 0 | 0 | — |
| rRNA | 0 | 0 | 0 | — | 0 | 0 | 0 | — |

#### Diferencias

- Falsos negativos (referencia sin pareja): **122**.
- Falsos positivos (herramienta sin pareja): **294**.
- Pares solo con hebra distinta: **13**.
- Pares emparejados con algún borde distinto: **872**.

Las listas de diferencias están truncadas a 100 elementos por lista; los conteos son completos.

> Primeras diferencias no emparejadas

| tipo | feature | seqid | coordenadas | hebra | product/locus |
|---|---|---|---|---|---|
| FN | CDS | NC_000962.3 | 9914..10828 | + | membrane protein |
| FN | CDS | NC_000962.3 | 148491..148859 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 309699..310073 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 399535..400050 | - | PE family protein PE6 |
| FN | CDS | NC_000962.3 | 476394..476642 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 547076..547357 | - | toxin MazF1 |
| FN | CDS | NC_000962.3 | 547344..547517 | - | antitoxin MazE1 |
| FN | CDS | NC_000962.3 | 703830..704057 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 759136..759309 | + | membrane protein |
| FN | CDS | NC_000962.3 | 789157..789411 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 798833..799444 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 799629..799850 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 835154..835681 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 937593..937835 | + | transmembrane protein |
| FN | CDS | NC_000962.3 | 992598..993455 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 1052360..1052638 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 1076196..1076678 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 1094670..1094864 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 1114293..1114724 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 1131128..1131421 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 1214513..1214947 | + | PE family protein PE9 |
| FN | CDS | NC_000962.3 | 1215225..1215407 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 1265034..1265438 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 1371529..1371702 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 1376230..1376826 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 1484982..1485278 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 1547072..1547530 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 1612256..1612393 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 1612342..1612950 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 1690407..1690805 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 1697356..1697859 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 1710733..1711002 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 1746919..1747146 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 1779930..1780241 | + | phage protein |
| FN | CDS | NC_000962.3 | 1788162..1789163 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 1895725..1896087 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 1981130..1981579 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 2047687..2048034 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 2159191..2159943 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 2181262..2181906 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 2182460..2183239 | - | transcriptional regulator |
| FN | CDS | NC_000962.3 | 2197508..2197960 | - | lipoprotein |
| FN | CDS | NC_000962.3 | 2201223..2201744 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 2231026..2231454 | + | chitinase |
| FN | CDS | NC_000962.3 | 2315174..2316709 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 2321057..2321467 | + | mRNA interferase MazF7 |
| FN | CDS | NC_000962.3 | 2334295..2334594 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 2338065..2338505 | - | transmembrane protein |
| FN | CDS | NC_000962.3 | 2343027..2343332 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 2497225..2497446 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 2505736..2506161 | - | ribonuclease VapC16 |
| FN | CDS | NC_000962.3 | 2506207..2506383 | - | antitoxin VapB16 |
| FN | CDS | NC_000962.3 | 2528980..2529174 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 2544289..2544621 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 2546488..2546805 | - | toxin MazF8 |
| FN | CDS | NC_000962.3 | 2546839..2547087 | - | antitoxin MazE8 |
| FN | CDS | NC_000962.3 | 2561328..2561705 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 2562599..2563114 | + | lipoprotein LppO |
| FN | CDS | NC_000962.3 | 2579228..2579419 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 2579504..2579935 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 2580028..2580210 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 2583045..2583332 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 2604740..2605078 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 2649946..2650248 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 2692224..2692439 | + | acid and phagosome regulated protein AprA |
| FN | CDS | NC_000962.3 | 2692551..2692715 | + | acid and phagosome regulated protein AprB |
| FN | CDS | NC_000962.3 | 2697728..2698057 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 2705762..2706736 | + | PE family protein PE24 |
| FN | CDS | NC_000962.3 | 2728437..2728847 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 2736709..2736987 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 2752262..2752660 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 2752848..2752994 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 2775272..2775565 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 2786575..2786898 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 2800846..2801145 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 2808083..2808304 | + | antitoxin VapB38 |
| FN | CDS | NC_000962.3 | 2854683..2854907 | - | antitoxin VapB39 |
| FN | CDS | NC_000962.3 | 2943376..2943603 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 2962470..2962712 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 2968533..2968850 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 2980963..2981190 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 2981187..2981576 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 2981482..2981754 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 3020200..3020457 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 3034635..3034883 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 3042542..3043018 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 3055515..3056348 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 3076367..3076720 | - | membrane protein |
| FN | CDS | NC_000962.3 | 3082352..3082756 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 3112465..3113094 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 3147959..3148258 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 3174747..3174995 | + | antitoxin VapB23 |
| FN | CDS | NC_000962.3 | 3183138..3183395 | + | antitoxin VapB43 |
| FN | CDS | NC_000962.3 | 3200266..3201015 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 3242970..3243281 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 3289705..3290235 | + | transposase |
| FN | CDS | NC_000962.3 | 3312953..3313201 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 3325934..3326104 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 3467210..3467662 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 3477060..3477500 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 9963..10160 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 10523..10828 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 11421..11528 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 11555..11692 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 13024..13176 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 31824..31967 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 99379..99573 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 148677..148859 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 180805..181029 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 223607..223909 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 226700..226819 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 257783..258856 | - | Acyl-CoA dehydrogenase |
| FP | CDS | NC_000962.3 | 274710..274904 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 339107..339367 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 399487..400098 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 418901..419104 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 466672..467406 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 467459..468001 | - | PPE family protein |
| FP | CDS | NC_000962.3 | 472914..473105 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 544736..545059 | - | HNH endonuclease |
| FP | CDS | NC_000962.3 | 545188..545418 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 591475..591576 | + | Ribosomal protein mS38 C-terminal domain-containing protein |
| FP | CDS | NC_000962.3 | 610935..611051 | - | Transmembrane protein |
| FP | CDS | NC_000962.3 | 612255..612605 | + | Methyltransferase/methylase |
| FP | CDS | NC_000962.3 | 612598..612903 | + | Methyltransferase |
| FP | CDS | NC_000962.3 | 664923..665045 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 665434..665601 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 688032..688859 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 688808..689062 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 697020..697193 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 702016..702759 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 704475..704672 | - | Uncharacterized protein Rv0609B |
| FP | CDS | NC_000962.3 | 711536..712231 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 712249..712719 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 759294..759524 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 767428..767541 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 789172..789420 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 799264..800019 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 832534..832848 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 835059..835205 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 850342..850512 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 874233..874943 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 874882..876390 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 876342..876503 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 889072..889398 | + | IS3 family IS987 transposase ORF A |
| FP | CDS | NC_000962.3 | 889449..890333 | + | IS3 family transposase |
| FP | CDS | NC_000962.3 | 908015..908113 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 931578..931991 | - | Transmembrane protein |
| FP | CDS | NC_000962.3 | 947312..947644 | + | Transposase |
| FP | CDS | NC_000962.3 | 986674..986823 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 993435..993752 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 1014397..1014555 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 1016820..1017050 | - | Glyoxalase |
| FP | CDS | NC_000962.3 | 1023865..1024122 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 1075111..1075251 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 1114173..1114382 | - | Methyltransferase |
| FP | CDS | NC_000962.3 | 1148080..1148172 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 1151777..1151923 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 1158918..1159409 | - | Transposase |
| FP | CDS | NC_000962.3 | 1159375..1160061 | - | Transposase |
| FP | CDS | NC_000962.3 | 1160095..1160418 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 1164164..1164436 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 1172110..1172313 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 1174995..1175621 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 1176011..1176265 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 1176372..1176485 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 1178604..1178750 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 1192235..1192375 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 1214049..1214360 | + | (2Z,6E)-farnesyl diphosphate synthase |
| FP | CDS | NC_000962.3 | 1214826..1215131 | + | PE family protein |
| FP | CDS | NC_000962.3 | 1215599..1216054 | + | putative cellulase Cel12b |
| FP | CDS | NC_000962.3 | 1231301..1231990 | + | Para-nitrobenzyl esterase |
| FP | CDS | NC_000962.3 | 1232311..1232826 | + | Para-nitrobenzyl esterase |
| FP | CDS | NC_000962.3 | 1241115..1241390 | - | PE family protein |
| FP | CDS | NC_000962.3 | 1242703..1242912 | + | Secreted protein |
| FP | CDS | NC_000962.3 | 1264314..1264556 | + | Beta-ketoadipyl-CoA thiolase |
| FP | CDS | NC_000962.3 | 1264606..1264947 | + | Enoyl-CoA hydratase |
| FP | CDS | NC_000962.3 | 1266990..1267439 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 1272036..1272167 | - | DUF551 domain-containing protein |
| FP | CDS | NC_000962.3 | 1278383..1278820 | + | Transposase IS4-like domain-containing protein |
| FP | CDS | NC_000962.3 | 1282890..1283036 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 1301248..1301355 | + | Transposase |
| FP | CDS | NC_000962.3 | 1351064..1351309 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 1376230..1376526 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 1451779..1452003 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 1526236..1526376 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 1530973..1531209 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 1541994..1542878 | - | IS3 family transposase |
| FP | CDS | NC_000962.3 | 1542929..1543255 | - | IS3 family IS987 transposase ORF A |
| FP | CDS | NC_000962.3 | 1612198..1612872 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 1616365..1616592 | - | Transcriptional regulator Rv0078-like C-terminal domain-containing protein |
| FP | CDS | NC_000962.3 | 1617565..1617729 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 1622768..1623007 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 1623124..1623237 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 1693996..1694436 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 1694545..1695108 | - | UDP-4-amino-4-deoxy-L-arabinose--oxoglutarate aminotransferase |
| FP | CDS | NC_000962.3 | 1701098..1701346 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 1712057..1712257 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 1746958..1747191 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 1753606..1754037 | + | hypothetical protein |

> Primeros pares con coordenadas no exactas

| tipo | seqid | ref | Δstart | Δend | solape |
|---|---|---|---|---|---|
| CDS | NC_000962.3 | 4434..4997 | +48 | +0 | 0.9149 |
| CDS | NC_000962.3 | 5240..7267 | -117 | +0 | 0.9455 |
| CDS | NC_000962.3 | 13133..13558 | +0 | -66 | 0.8451 |
| CDS | NC_000962.3 | 14089..14877 | +45 | +0 | 0.9430 |
| CDS | NC_000962.3 | 21637..23181 | +0 | -9 | 0.9942 |
| CDS | NC_000962.3 | 31514..31819 | +78 | +0 | 0.7451 |
| CDS | NC_000962.3 | 32057..33154 | +222 | +0 | 0.7978 |
| CDS | NC_000962.3 | 39877..41202 | +0 | -15 | 0.9887 |
| CDS | NC_000962.3 | 41304..41912 | -24 | +0 | 0.9621 |
| CDS | NC_000962.3 | 43562..46471 | +33 | +0 | 0.9887 |
| CDS | NC_000962.3 | 53663..55699 | -426 | +0 | 0.8270 |
| CDS | NC_000962.3 | 57410..57973 | -75 | +0 | 0.8826 |
| CDS | NC_000962.3 | 65012..65350 | +0 | +315 | 0.5183 |
| CDS | NC_000962.3 | 65552..66694 | +219 | +0 | 0.8084 |
| CDS | NC_000962.3 | 68418..68579 | +6 | +0 | 0.9630 |
| CDS | NC_000962.3 | 95414..96892 | +3 | +0 | 0.9980 |
| CDS | NC_000962.3 | 96927..97601 | +195 | +0 | 0.7111 |
| CDS | NC_000962.3 | 102815..103663 | +0 | -138 | 0.8375 |
| CDS | NC_000962.3 | 103710..104663 | +0 | +6 | 0.9938 |
| CDS | NC_000962.3 | 104805..105215 | +0 | -45 | 0.8905 |
| CDS | NC_000962.3 | 108156..109778 | +24 | +0 | 0.9852 |
| CDS | NC_000962.3 | 109783..110019 | -12 | +0 | 0.9518 |
| CDS | NC_000962.3 | 117714..119699 | +63 | +0 | 0.9683 |
| CDS | NC_000962.3 | 130895..131104 | +0 | +33 | 0.8642 |
| CDS | NC_000962.3 | 137319..137909 | -33 | +0 | 0.9471 |
| CDS | NC_000962.3 | 144049..145626 | -51 | +0 | 0.9687 |
| CDS | NC_000962.3 | 173238..174758 | +48 | +0 | 0.9684 |
| CDS | NC_000962.3 | 179319..180896 | +0 | -105 | 0.9335 |
| CDS | NC_000962.3 | 186495..186623 | +0 | -3 | 0.9767 |
| CDS | NC_000962.3 | 190607..191956 | -24 | +0 | 0.9825 |
| CDS | NC_000962.3 | 191984..193135 | +0 | -33 | 0.9714 |
| CDS | NC_000962.3 | 193626..194111 | +45 | +0 | 0.9074 |
| CDS | NC_000962.3 | 194144..194815 | +0 | +123 | 0.8453 |
| CDS | NC_000962.3 | 206814..207455 | -66 | +0 | 0.9068 |
| CDS | NC_000962.3 | 213028..214140 | +0 | -114 | 0.8976 |
| CDS | NC_000962.3 | 214088..214927 | -57 | +0 | 0.9365 |
| CDS | NC_000962.3 | 218390..218551 | +0 | -15 | 0.9074 |
| CDS | NC_000962.3 | 218705..219367 | +6 | +0 | 0.9910 |
| CDS | NC_000962.3 | 223564..224664 | +336 | +0 | 0.6948 |
| CDS | NC_000962.3 | 230899..231534 | -132 | +0 | 0.8281 |
| CDS | NC_000962.3 | 238392..241292 | +0 | -3 | 0.9990 |
| CDS | NC_000962.3 | 241976..243214 | +0 | +45 | 0.9650 |
| CDS | NC_000962.3 | 247384..248112 | +0 | +6 | 0.9918 |
| CDS | NC_000962.3 | 248115..248906 | +0 | -198 | 0.7500 |
| CDS | NC_000962.3 | 254637..255950 | +0 | +3 | 0.9977 |
| CDS | NC_000962.3 | 256064..257677 | +96 | +0 | 0.9405 |
| CDS | NC_000962.3 | 271574..272839 | +0 | +15 | 0.9883 |
| CDS | NC_000962.3 | 273055..274278 | -75 | +0 | 0.9423 |
| CDS | NC_000962.3 | 279605..281140 | +0 | -162 | 0.8945 |
| CDS | NC_000962.3 | 288428..289042 | +48 | +0 | 0.9220 |
| CDS | NC_000962.3 | 302866..305427 | +51 | +0 | 0.9801 |
| CDS | NC_000962.3 | 315783..316415 | +0 | +21 | 0.9679 |
| CDS | NC_000962.3 | 322764..323273 | +0 | -252 | 0.5059 |
| CDS | NC_000962.3 | 330933..331658 | +0 | +9 | 0.9878 |
| CDS | NC_000962.3 | 356862..358064 | +0 | +9 | 0.9926 |
| CDS | NC_000962.3 | 383602..384486 | +30 | +0 | 0.9661 |
| CDS | NC_000962.3 | 385173..385943 | +0 | -33 | 0.9572 |
| CDS | NC_000962.3 | 386305..387099 | +0 | -18 | 0.9774 |
| CDS | NC_000962.3 | 398254..398628 | -60 | +0 | 0.8621 |
| CDS | NC_000962.3 | 400192..401703 | +69 | +0 | 0.9544 |
| CDS | NC_000962.3 | 405950..408448 | +0 | -48 | 0.9808 |
| CDS | NC_000962.3 | 415050..415460 | -6 | +0 | 0.9856 |
| CDS | NC_000962.3 | 417304..418290 | +174 | +0 | 0.8237 |
| CDS | NC_000962.3 | 418949..419608 | +243 | +0 | 0.6318 |
| CDS | NC_000962.3 | 438302..438739 | +0 | -33 | 0.9247 |
| CDS | NC_000962.3 | 445314..446525 | +0 | +9 | 0.9926 |
| CDS | NC_000962.3 | 446531..447046 | +0 | +87 | 0.8557 |
| CDS | NC_000962.3 | 452294..453154 | +0 | +42 | 0.9535 |
| CDS | NC_000962.3 | 455637..455858 | -3 | +0 | 0.9867 |
| CDS | NC_000962.3 | 456268..456819 | +0 | +99 | 0.8479 |
| CDS | NC_000962.3 | 462135..463307 | -105 | +0 | 0.9178 |
| CDS | NC_000962.3 | 469591..470013 | +42 | +0 | 0.9007 |
| CDS | NC_000962.3 | 472781..474106 | +339 | +0 | 0.7443 |
| CDS | NC_000962.3 | 475350..475742 | -24 | +0 | 0.9424 |
| CDS | NC_000962.3 | 483977..485734 | +552 | +0 | 0.6860 |
| CDS | NC_000962.3 | 499713..500366 | +60 | +0 | 0.9083 |
| CDS | NC_000962.3 | 506561..506971 | +0 | -69 | 0.8321 |
| CDS | NC_000962.3 | 507758..508555 | +0 | +30 | 0.9638 |
| CDS | NC_000962.3 | 519073..519567 | +21 | +0 | 0.9576 |
| CDS | NC_000962.3 | 520324..521454 | -21 | +0 | 0.9818 |
| CDS | NC_000962.3 | 526143..527360 | +0 | -39 | 0.9680 |
| CDS | NC_000962.3 | 530751..532214 | +0 | -51 | 0.9652 |
| CDS | NC_000962.3 | 532396..532911 | -42 | +0 | 0.9247 |
| CDS | NC_000962.3 | 534445..535215 | +0 | +12 | 0.9847 |
| CDS | NC_000962.3 | 536504..537169 | +0 | +66 | 0.9098 |
| CDS | NC_000962.3 | 537229..538548 | +0 | -24 | 0.9818 |
| CDS | NC_000962.3 | 541488..541910 | +0 | -6 | 0.9858 |
| CDS | NC_000962.3 | 542142..542852 | -51 | +0 | 0.9331 |
| CDS | NC_000962.3 | 547586..549607 | +0 | +9 | 0.9956 |
| CDS | NC_000962.3 | 551749..551988 | -39 | +0 | 0.8602 |
| CDS | NC_000962.3 | 552026..552550 | -45 | +0 | 0.9211 |
| CDS | NC_000962.3 | 561854..562294 | +10 | -146 | 0.6463 |
| CDS | NC_000962.3 | 562723..563427 | +0 | +66 | 0.9144 |
| CDS | NC_000962.3 | 567921..568967 | +0 | -297 | 0.7163 |
| CDS | NC_000962.3 | 575348..576790 | +9 | +0 | 0.9938 |
| CDS | NC_000962.3 | 577664..578269 | -21 | +0 | 0.9665 |
| CDS | NC_000962.3 | 578426..579175 | -21 | +0 | 0.9728 |
| CDS | NC_000962.3 | 583701..584690 | +0 | -39 | 0.9606 |
| CDS | NC_000962.3 | 584695..585423 | +12 | +0 | 0.9835 |
| CDS | NC_000962.3 | 585424..586314 | +0 | -120 | 0.8653 |


## Prokka

- Estado: **ok**.
- Fuente: `executed`.
- Versión: `prokka 1.15.6` (estado `ok`).
- Ejecución: **ok**, 960.60 s, RSS máx. 2654.8 MiB (`2026-11-25T21:51:21Z` → `2026-11-25T22:07:21Z`).
- Salida: `/tmp/external_benchmark_prokka_y_o1bbxr/out/benchmark.gff` (5433507 bytes, SHA256 `11ee07681d4af59a…`).
- Comando: `/usr/local/bin/prokka --outdir /tmp/external_benchmark_prokka_y_o1bbxr/out --prefix benchmark --cpus 4 --force --quiet /tmp/benchmark_genome_azi2xbyp.fna`.
- Anotación comparada: `reports/external_benchmark/multi_v2/mtb_h37rv/evidence/prokka/benchmark.gff` (formato `gff3`, SHA256 `11ee07681d4af59a…`).

### Features (Prokka)

| tipo | ref | query | TP | FP | FN | precisión | recall | F1 | exactos | borde distinto | solape medio |
|---|---|---|---|---|---|---|---|---|---|---|---|
| CDS | 3906 | 4069 | 3774 | 295 | 132 | 0.9275 | 0.9662 | 0.9465 | 2920 | 854 | 0.9819 |
| tRNA | 45 | 53 | 45 | 8 | 0 | 0.8491 | 1.0000 | 0.9184 | 0 | 45 | 0.9781 |
| rRNA | 3 | 3 | 3 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 1 | 2 | 0.9931 |

### Atributos (Prokka)

| tipo | pares | product en ambos | product literal | product tokens | gene literal | locus literal |
|---|---|---|---|---|---|---|
| CDS | 3774 | 3774 | 0.3397 | 0.5401 | 0.7206 | 0.0000 |
| tRNA | 45 | 45 | 0.0000 | 0.6667 | — | 0.0000 |
| rRNA | 3 | 3 | 1.0000 | 1.0000 | — | 0.0000 |

### Hipotéticas (Prokka)

| tipo | ref CDS | ref con product | ref hyp | % ref | query CDS | query con product | query hyp | % query |
|---|---|---|---|---|---|---|---|---|
| CDS | 3906 | 3906 | 1052 | 26.9 | 4069 | 4069 | 1298 | 31.9 |
| tRNA | 0 | 0 | 0 | — | 0 | 0 | 0 | — |
| rRNA | 0 | 0 | 0 | — | 0 | 0 | 0 | — |

#### Diferencias

- Falsos negativos (referencia sin pareja): **132**.
- Falsos positivos (herramienta sin pareja): **303**.
- Pares solo con hebra distinta: **13**.
- Pares emparejados con algún borde distinto: **901**.

Las listas de diferencias están truncadas a 100 elementos por lista; los conteos son completos.

> Primeras diferencias no emparejadas

| tipo | feature | seqid | coordenadas | hebra | product/locus |
|---|---|---|---|---|---|
| FN | CDS | NC_000962.3 | 9914..10828 | + | membrane protein |
| FN | CDS | NC_000962.3 | 148491..148859 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 309699..310073 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 399535..400050 | - | PE family protein PE6 |
| FN | CDS | NC_000962.3 | 476394..476642 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 547076..547357 | - | toxin MazF1 |
| FN | CDS | NC_000962.3 | 547344..547517 | - | antitoxin MazE1 |
| FN | CDS | NC_000962.3 | 703830..704057 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 759136..759309 | + | membrane protein |
| FN | CDS | NC_000962.3 | 789157..789411 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 798833..799444 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 799629..799850 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 835154..835681 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 927837..930485 | - | PE-PGRS family protein PE_PGRS14 |
| FN | CDS | NC_000962.3 | 937593..937835 | + | transmembrane protein |
| FN | CDS | NC_000962.3 | 992598..993455 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 1052360..1052638 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 1076196..1076678 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 1094670..1094864 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 1114293..1114724 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 1131128..1131421 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 1214513..1214947 | + | PE family protein PE9 |
| FN | CDS | NC_000962.3 | 1215225..1215407 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 1265034..1265438 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 1371529..1371702 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 1376230..1376826 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 1484982..1485278 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 1547072..1547530 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 1581145..1581969 | - | methyltransferase |
| FN | CDS | NC_000962.3 | 1612256..1612393 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 1612342..1612950 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 1626959..1628080 | - | transaldolase |
| FN | CDS | NC_000962.3 | 1690407..1690805 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 1697356..1697859 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 1710733..1711002 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 1746919..1747146 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 1779930..1780241 | + | phage protein |
| FN | CDS | NC_000962.3 | 1788162..1789163 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 1862347..1865382 | - | PE-PGRS family protein PE_PGRS30 |
| FN | CDS | NC_000962.3 | 1895725..1896087 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 1981130..1981579 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 2047687..2048034 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 2159191..2159943 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 2181262..2181906 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 2182460..2183239 | - | transcriptional regulator |
| FN | CDS | NC_000962.3 | 2197508..2197960 | - | lipoprotein |
| FN | CDS | NC_000962.3 | 2201223..2201744 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 2231026..2231454 | + | chitinase |
| FN | CDS | NC_000962.3 | 2315174..2316709 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 2321057..2321467 | + | mRNA interferase MazF7 |
| FN | CDS | NC_000962.3 | 2334295..2334594 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 2338065..2338505 | - | transmembrane protein |
| FN | CDS | NC_000962.3 | 2343027..2343332 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 2497225..2497446 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 2505736..2506161 | - | ribonuclease VapC16 |
| FN | CDS | NC_000962.3 | 2506207..2506383 | - | antitoxin VapB16 |
| FN | CDS | NC_000962.3 | 2528980..2529174 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 2544289..2544621 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 2546488..2546805 | - | toxin MazF8 |
| FN | CDS | NC_000962.3 | 2546839..2547087 | - | antitoxin MazE8 |
| FN | CDS | NC_000962.3 | 2561328..2561705 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 2562599..2563114 | + | lipoprotein LppO |
| FN | CDS | NC_000962.3 | 2579228..2579419 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 2579504..2579935 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 2580028..2580210 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 2583045..2583332 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 2604740..2605078 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 2649946..2650248 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 2692224..2692439 | + | acid and phagosome regulated protein AprA |
| FN | CDS | NC_000962.3 | 2692551..2692715 | + | acid and phagosome regulated protein AprB |
| FN | CDS | NC_000962.3 | 2697728..2698057 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 2705762..2706736 | + | PE family protein PE24 |
| FN | CDS | NC_000962.3 | 2728437..2728847 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 2736709..2736987 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 2752262..2752660 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 2752848..2752994 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 2775272..2775565 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 2786575..2786898 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 2800846..2801145 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 2808083..2808304 | + | antitoxin VapB38 |
| FN | CDS | NC_000962.3 | 2854683..2854907 | - | antitoxin VapB39 |
| FN | CDS | NC_000962.3 | 2943376..2943603 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 2962470..2962712 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 2968533..2968850 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 2980963..2981190 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 2981187..2981576 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 2981482..2981754 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 3020200..3020457 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 3034635..3034883 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 3042542..3043018 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 3055515..3056348 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 3076367..3076720 | - | membrane protein |
| FN | CDS | NC_000962.3 | 3082352..3082756 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 3112465..3113094 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 3142309..3143628 | - | DNA-damage-inducible protein DinF |
| FN | CDS | NC_000962.3 | 3147959..3148258 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 3174747..3174995 | + | antitoxin VapB23 |
| FN | CDS | NC_000962.3 | 3183138..3183395 | + | antitoxin VapB43 |
| FN | CDS | NC_000962.3 | 3200266..3201015 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 3242970..3243281 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 9963..10160 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 10523..10828 | + | putative protein |
| FP | CDS | NC_000962.3 | 11421..11528 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 11555..11692 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 13024..13176 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 31824..31967 | - | putative proteinA |
| FP | CDS | NC_000962.3 | 99379..99573 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 148677..148859 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 180805..181029 | - | putative PE-PGRS family protein PE_PGRS24 |
| FP | CDS | NC_000962.3 | 223607..223909 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 226700..226819 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 257783..258856 | - | (2S)-methylsuccinyl-CoA dehydrogenase |
| FP | CDS | NC_000962.3 | 274710..274904 | - | Putative antitoxin VapB51 |
| FP | CDS | NC_000962.3 | 339107..339367 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 399487..400098 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 418901..419104 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 466672..467406 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 467459..468001 | - | PPE family protein PPE32 |
| FP | CDS | NC_000962.3 | 472914..473105 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 544736..545059 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 545188..545418 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 591475..591576 | + | putative proteinB |
| FP | CDS | NC_000962.3 | 610935..611051 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 612255..612605 | + | Sarcosine/dimethylglycine N-methyltransferase |
| FP | CDS | NC_000962.3 | 612598..612903 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 664923..665045 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 665434..665601 | - | putative proteinA |
| FP | CDS | NC_000962.3 | 688032..688859 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 688808..689062 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 697020..697193 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 702016..702759 | + | IS607 family transposase IS1536 |
| FP | CDS | NC_000962.3 | 704475..704672 | - | putative proteinB |
| FP | CDS | NC_000962.3 | 711536..712231 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 712249..712719 | + | Galactose-1-phosphate uridylyltransferase |
| FP | CDS | NC_000962.3 | 759294..759524 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 767428..767541 | + | Putative amidase AmiD |
| FP | CDS | NC_000962.3 | 789172..789420 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 799264..800019 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 832534..832848 | + | ISL3 family transposase IS1557 |
| FP | CDS | NC_000962.3 | 835059..835205 | - | putative 3-hydroxyisobutyrate dehydrogenase |
| FP | CDS | NC_000962.3 | 850342..850512 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 874233..874943 | + | Dipeptidyl aminopeptidase BI |
| FP | CDS | NC_000962.3 | 874882..876390 | + | Dipeptidyl aminopeptidase BI |
| FP | CDS | NC_000962.3 | 876342..876503 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 889072..889398 | + | IS3 family transposase IS987 |
| FP | CDS | NC_000962.3 | 889449..890333 | + | IS3 family transposase IS987 |
| FP | CDS | NC_000962.3 | 908015..908113 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 931578..931991 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 947312..947644 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 986674..986823 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 993435..993752 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 1014397..1014555 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 1016820..1017050 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 1023865..1024122 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 1075111..1075251 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 1114173..1114382 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 1148080..1148172 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 1151777..1151923 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 1158918..1159409 | - | ISNCY family transposase ISRjo3 |
| FP | CDS | NC_000962.3 | 1159375..1160061 | - | ISNCY family transposase ISRjo3 |
| FP | CDS | NC_000962.3 | 1160095..1160418 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 1164164..1164436 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 1172110..1172313 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 1174995..1175621 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 1176011..1176265 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 1176372..1176485 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 1178604..1178750 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 1192235..1192375 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 1214049..1214360 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 1214826..1215131 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 1215599..1216054 | + | putative cellulase Cel12b |
| FP | CDS | NC_000962.3 | 1231301..1231990 | + | Para-nitrobenzyl esterase |
| FP | CDS | NC_000962.3 | 1232311..1232826 | + | Para-nitrobenzyl esterase |
| FP | CDS | NC_000962.3 | 1241115..1241390 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 1242703..1242912 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 1264314..1264556 | + | 3-ketoacyl-CoA thiolase |
| FP | CDS | NC_000962.3 | 1264606..1264947 | + | Carnitinyl-CoA dehydratase |
| FP | CDS | NC_000962.3 | 1266990..1267439 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 1272036..1272167 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 1278383..1278820 | + | IS5 family transposase ISMt1 |
| FP | CDS | NC_000962.3 | 1282890..1283036 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 1301248..1301355 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 1351064..1351309 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 1376230..1376526 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 1451779..1452003 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 1526236..1526376 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 1530973..1531209 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 1541994..1542878 | - | IS3 family transposase IS987 |
| FP | CDS | NC_000962.3 | 1542929..1543255 | - | IS3 family transposase IS987 |
| FP | CDS | NC_000962.3 | 1612198..1612872 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 1616365..1616592 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 1617565..1617729 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 1622768..1623007 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 1623124..1623237 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 1693996..1694436 | - | dTDP-4-amino-4,6-dideoxygalactose transaminase |
| FP | CDS | NC_000962.3 | 1694545..1695108 | - | dTDP-4-amino-4,6-dideoxygalactose transaminase |
| FP | CDS | NC_000962.3 | 1701098..1701346 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 1712057..1712257 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 1746958..1747191 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 1753606..1754037 | + | Long-chain-fatty-acid--CoA ligase FadD13 |

> Primeros pares con coordenadas no exactas

| tipo | seqid | ref | Δstart | Δend | solape |
|---|---|---|---|---|---|
| CDS | NC_000962.3 | 1..1524 | +33 | +0 | 0.9783 |
| CDS | NC_000962.3 | 4434..4997 | +48 | +0 | 0.9149 |
| CDS | NC_000962.3 | 5240..7267 | -117 | +0 | 0.9455 |
| CDS | NC_000962.3 | 13133..13558 | +0 | -66 | 0.8451 |
| CDS | NC_000962.3 | 14089..14877 | +45 | +0 | 0.9430 |
| CDS | NC_000962.3 | 21637..23181 | +0 | -9 | 0.9942 |
| CDS | NC_000962.3 | 31514..31819 | +78 | +0 | 0.7451 |
| CDS | NC_000962.3 | 32057..33154 | +222 | +0 | 0.7978 |
| CDS | NC_000962.3 | 39877..41202 | +0 | -15 | 0.9887 |
| CDS | NC_000962.3 | 41304..41912 | -24 | +0 | 0.9621 |
| CDS | NC_000962.3 | 43562..46471 | +33 | +0 | 0.9887 |
| CDS | NC_000962.3 | 53663..55699 | -426 | +0 | 0.8270 |
| CDS | NC_000962.3 | 57410..57973 | -75 | +0 | 0.8826 |
| CDS | NC_000962.3 | 65012..65350 | +0 | +315 | 0.5183 |
| CDS | NC_000962.3 | 65552..66694 | +219 | +0 | 0.8084 |
| CDS | NC_000962.3 | 68418..68579 | +6 | +0 | 0.9630 |
| CDS | NC_000962.3 | 95414..96892 | +3 | +0 | 0.9980 |
| CDS | NC_000962.3 | 96927..97601 | +195 | +0 | 0.7111 |
| CDS | NC_000962.3 | 102815..103663 | +0 | -138 | 0.8375 |
| CDS | NC_000962.3 | 103710..104663 | +0 | +6 | 0.9938 |
| CDS | NC_000962.3 | 104805..105215 | +0 | -45 | 0.8905 |
| CDS | NC_000962.3 | 108156..109778 | +24 | +0 | 0.9852 |
| CDS | NC_000962.3 | 109783..110019 | -12 | +0 | 0.9518 |
| CDS | NC_000962.3 | 117714..119699 | +63 | +0 | 0.9683 |
| CDS | NC_000962.3 | 130895..131104 | +0 | +33 | 0.8642 |
| CDS | NC_000962.3 | 137319..137909 | -33 | +0 | 0.9471 |
| CDS | NC_000962.3 | 144049..145626 | -51 | +0 | 0.9687 |
| CDS | NC_000962.3 | 173238..174758 | +48 | +0 | 0.9684 |
| CDS | NC_000962.3 | 179319..180896 | +0 | -105 | 0.9335 |
| CDS | NC_000962.3 | 186495..186623 | +0 | -3 | 0.9767 |
| CDS | NC_000962.3 | 190607..191956 | -24 | +0 | 0.9825 |
| CDS | NC_000962.3 | 191984..193135 | +0 | -33 | 0.9714 |
| CDS | NC_000962.3 | 193626..194111 | +45 | +0 | 0.9074 |
| CDS | NC_000962.3 | 194144..194815 | +0 | +123 | 0.8453 |
| CDS | NC_000962.3 | 206814..207455 | -66 | +0 | 0.9068 |
| CDS | NC_000962.3 | 213028..214140 | +0 | -114 | 0.8976 |
| CDS | NC_000962.3 | 214088..214927 | -57 | +0 | 0.9365 |
| CDS | NC_000962.3 | 218390..218551 | +0 | -15 | 0.9074 |
| CDS | NC_000962.3 | 218705..219367 | +6 | +0 | 0.9910 |
| CDS | NC_000962.3 | 223564..224664 | +336 | +0 | 0.6948 |
| CDS | NC_000962.3 | 230899..231534 | -132 | +0 | 0.8281 |
| CDS | NC_000962.3 | 238392..241292 | +0 | -3 | 0.9990 |
| CDS | NC_000962.3 | 241976..243214 | +0 | +45 | 0.9650 |
| CDS | NC_000962.3 | 247384..248112 | +0 | +6 | 0.9918 |
| CDS | NC_000962.3 | 248115..248906 | +0 | -198 | 0.7500 |
| CDS | NC_000962.3 | 254637..255950 | +0 | +3 | 0.9977 |
| CDS | NC_000962.3 | 256064..257677 | +96 | +0 | 0.9405 |
| CDS | NC_000962.3 | 271574..272839 | +0 | +15 | 0.9883 |
| CDS | NC_000962.3 | 273055..274278 | -75 | +0 | 0.9423 |
| CDS | NC_000962.3 | 279605..281140 | +0 | -162 | 0.8945 |
| CDS | NC_000962.3 | 288428..289042 | +48 | +0 | 0.9220 |
| CDS | NC_000962.3 | 302866..305427 | +51 | +0 | 0.9801 |
| CDS | NC_000962.3 | 315783..316415 | +0 | +21 | 0.9679 |
| CDS | NC_000962.3 | 322764..323273 | +0 | -252 | 0.5059 |
| CDS | NC_000962.3 | 330933..331658 | +0 | +9 | 0.9878 |
| CDS | NC_000962.3 | 356862..358064 | +0 | +9 | 0.9926 |
| CDS | NC_000962.3 | 383602..384486 | +30 | +0 | 0.9661 |
| CDS | NC_000962.3 | 385173..385943 | +0 | -33 | 0.9572 |
| CDS | NC_000962.3 | 386305..387099 | +0 | -18 | 0.9774 |
| CDS | NC_000962.3 | 398254..398628 | -60 | +0 | 0.8621 |
| CDS | NC_000962.3 | 400192..401703 | +69 | +0 | 0.9544 |
| CDS | NC_000962.3 | 405950..408448 | +0 | -48 | 0.9808 |
| CDS | NC_000962.3 | 415050..415460 | -6 | +0 | 0.9856 |
| CDS | NC_000962.3 | 417304..418290 | +174 | +0 | 0.8237 |
| CDS | NC_000962.3 | 418949..419608 | +243 | +0 | 0.6318 |
| CDS | NC_000962.3 | 438302..438739 | +0 | -33 | 0.9247 |
| CDS | NC_000962.3 | 445314..446525 | +0 | +9 | 0.9926 |
| CDS | NC_000962.3 | 446531..447046 | +0 | +87 | 0.8557 |
| CDS | NC_000962.3 | 452294..453154 | +0 | +42 | 0.9535 |
| CDS | NC_000962.3 | 455637..455858 | -3 | +0 | 0.9867 |
| CDS | NC_000962.3 | 456268..456819 | +0 | +99 | 0.8479 |
| CDS | NC_000962.3 | 462135..463307 | -105 | +0 | 0.9178 |
| CDS | NC_000962.3 | 469591..470013 | +42 | +0 | 0.9007 |
| CDS | NC_000962.3 | 472781..474106 | +339 | +0 | 0.7443 |
| CDS | NC_000962.3 | 475350..475742 | -24 | +0 | 0.9424 |
| CDS | NC_000962.3 | 483977..485734 | +552 | +0 | 0.6860 |
| CDS | NC_000962.3 | 499713..500366 | +60 | +0 | 0.9083 |
| CDS | NC_000962.3 | 506561..506971 | +0 | -69 | 0.8321 |
| CDS | NC_000962.3 | 507758..508555 | +0 | +30 | 0.9638 |
| CDS | NC_000962.3 | 519073..519567 | +21 | +0 | 0.9576 |
| CDS | NC_000962.3 | 520324..521454 | -21 | +0 | 0.9818 |
| CDS | NC_000962.3 | 526143..527360 | +0 | -39 | 0.9680 |
| CDS | NC_000962.3 | 530751..532214 | +0 | -51 | 0.9652 |
| CDS | NC_000962.3 | 532396..532911 | -42 | +0 | 0.9247 |
| CDS | NC_000962.3 | 534445..535215 | +0 | +12 | 0.9847 |
| CDS | NC_000962.3 | 536504..537169 | +0 | +66 | 0.9098 |
| CDS | NC_000962.3 | 537229..538548 | +0 | -24 | 0.9818 |
| CDS | NC_000962.3 | 541488..541910 | +0 | -6 | 0.9858 |
| CDS | NC_000962.3 | 542142..542852 | -51 | +0 | 0.9331 |
| CDS | NC_000962.3 | 547586..549607 | +0 | +9 | 0.9956 |
| CDS | NC_000962.3 | 551749..551988 | -39 | +0 | 0.8602 |
| CDS | NC_000962.3 | 552026..552550 | -45 | +0 | 0.9211 |
| CDS | NC_000962.3 | 561854..562294 | +10 | -146 | 0.6463 |
| CDS | NC_000962.3 | 562723..563427 | +0 | +66 | 0.9144 |
| CDS | NC_000962.3 | 567921..568967 | +0 | -297 | 0.7163 |
| CDS | NC_000962.3 | 575348..576790 | +9 | +0 | 0.9938 |
| CDS | NC_000962.3 | 577664..578269 | -21 | +0 | 0.9665 |
| CDS | NC_000962.3 | 578426..579175 | -21 | +0 | 0.9728 |
| CDS | NC_000962.3 | 583701..584690 | +0 | -39 | 0.9606 |
| CDS | NC_000962.3 | 584695..585423 | +12 | +0 | 0.9835 |


## DFAST

- Estado: **ok**.
- Fuente: `executed`.
- Versión: `DFAST ver. 1.6.2` (estado `ok`).
- Ejecución: **ok**, 521.16 s, RSS máx. 2654.8 MiB (`2026-11-25T22:07:42Z` → `2026-11-25T22:16:23Z`).
- Salida: `/tmp/external_benchmark_dfast_1qstp86m/out/genome.gff` (7856459 bytes, SHA256 `404142262124aaf4…`).
- Comando: `/usr/local/bin/dfast --genome /tmp/benchmark_genome_azi2xbyp.fna --out /tmp/external_benchmark_dfast_1qstp86m/out --cpu 4 --force --dbroot /data/dfast_db`.
- Anotación comparada: `reports/external_benchmark/multi_v2/mtb_h37rv/evidence/dfast/genome.gff` (formato `gff3`, SHA256 `404142262124aaf4…`).

### Features (DFAST)

| tipo | ref | query | TP | FP | FN | precisión | recall | F1 | exactos | borde distinto | solape medio |
|---|---|---|---|---|---|---|---|---|---|---|---|
| CDS | 3906 | 4074 | 3779 | 295 | 127 | 0.9276 | 0.9675 | 0.9471 | 2923 | 856 | 0.9819 |
| tRNA | 45 | 48 | 45 | 3 | 0 | 0.9375 | 1.0000 | 0.9677 | 0 | 45 | 0.9781 |
| rRNA | 3 | 3 | 3 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 0 | 3 | 0.9714 |

### Atributos (DFAST)

| tipo | pares | product en ambos | product literal | product tokens | gene literal | locus literal |
|---|---|---|---|---|---|---|
| CDS | 3779 | 3779 | 0.2315 | 0.4444 | 0.8436 | 0.0000 |
| tRNA | 45 | 45 | 1.0000 | 1.0000 | — | 0.0000 |
| rRNA | 3 | 3 | 1.0000 | 1.0000 | — | 0.0000 |

### Hipotéticas (DFAST)

| tipo | ref CDS | ref con product | ref hyp | % ref | query CDS | query con product | query hyp | % query |
|---|---|---|---|---|---|---|---|---|
| CDS | 3906 | 3906 | 1052 | 26.9 | 4074 | 4074 | 874 | 21.4 |
| tRNA | 0 | 0 | 0 | — | 0 | 0 | 0 | — |
| rRNA | 0 | 0 | 0 | — | 0 | 0 | 0 | — |

#### Diferencias

- Falsos negativos (referencia sin pareja): **127**.
- Falsos positivos (herramienta sin pareja): **298**.
- Pares solo con hebra distinta: **13**.
- Pares emparejados con algún borde distinto: **904**.

Las listas de diferencias están truncadas a 100 elementos por lista; los conteos son completos.

> Primeras diferencias no emparejadas

| tipo | feature | seqid | coordenadas | hebra | product/locus |
|---|---|---|---|---|---|
| FN | CDS | NC_000962.3 | 9914..10828 | + | membrane protein |
| FN | CDS | NC_000962.3 | 148491..148859 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 309699..310073 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 333437..336310 | - | PE-PGRS family protein PE_PGRS3 |
| FN | CDS | NC_000962.3 | 399535..400050 | - | PE family protein PE6 |
| FN | CDS | NC_000962.3 | 476394..476642 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 547076..547357 | - | toxin MazF1 |
| FN | CDS | NC_000962.3 | 547344..547517 | - | antitoxin MazE1 |
| FN | CDS | NC_000962.3 | 703830..704057 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 759136..759309 | + | membrane protein |
| FN | CDS | NC_000962.3 | 789157..789411 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 798833..799444 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 799629..799850 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 835154..835681 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 927837..930485 | - | PE-PGRS family protein PE_PGRS14 |
| FN | CDS | NC_000962.3 | 937593..937835 | + | transmembrane protein |
| FN | CDS | NC_000962.3 | 992598..993455 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 1052360..1052638 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 1076196..1076678 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 1094670..1094864 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 1114293..1114724 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 1131128..1131421 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 1214513..1214947 | + | PE family protein PE9 |
| FN | CDS | NC_000962.3 | 1215225..1215407 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 1265034..1265438 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 1371529..1371702 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 1376230..1376826 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 1484982..1485278 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 1547072..1547530 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 1612256..1612393 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 1612342..1612950 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 1690407..1690805 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 1697356..1697859 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 1710733..1711002 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 1746919..1747146 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 1779930..1780241 | + | phage protein |
| FN | CDS | NC_000962.3 | 1788162..1789163 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 1862347..1865382 | - | PE-PGRS family protein PE_PGRS30 |
| FN | CDS | NC_000962.3 | 1895725..1896087 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 1981130..1981579 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 2047687..2048034 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 2087971..2089518 | - | PE-PGRS family protein PE_PGRS34 |
| FN | CDS | NC_000962.3 | 2159191..2159943 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 2181262..2181906 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 2182460..2183239 | - | transcriptional regulator |
| FN | CDS | NC_000962.3 | 2197508..2197960 | - | lipoprotein |
| FN | CDS | NC_000962.3 | 2201223..2201744 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 2231026..2231454 | + | chitinase |
| FN | CDS | NC_000962.3 | 2315174..2316709 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 2321057..2321467 | + | mRNA interferase MazF7 |
| FN | CDS | NC_000962.3 | 2334295..2334594 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 2338065..2338505 | - | transmembrane protein |
| FN | CDS | NC_000962.3 | 2343027..2343332 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 2497225..2497446 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 2505736..2506161 | - | ribonuclease VapC16 |
| FN | CDS | NC_000962.3 | 2506207..2506383 | - | antitoxin VapB16 |
| FN | CDS | NC_000962.3 | 2528980..2529174 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 2544289..2544621 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 2546488..2546805 | - | toxin MazF8 |
| FN | CDS | NC_000962.3 | 2546839..2547087 | - | antitoxin MazE8 |
| FN | CDS | NC_000962.3 | 2561328..2561705 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 2562599..2563114 | + | lipoprotein LppO |
| FN | CDS | NC_000962.3 | 2579228..2579419 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 2579504..2579935 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 2580028..2580210 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 2583045..2583332 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 2604740..2605078 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 2649946..2650248 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 2692224..2692439 | + | acid and phagosome regulated protein AprA |
| FN | CDS | NC_000962.3 | 2692551..2692715 | + | acid and phagosome regulated protein AprB |
| FN | CDS | NC_000962.3 | 2697728..2698057 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 2705762..2706736 | + | PE family protein PE24 |
| FN | CDS | NC_000962.3 | 2728437..2728847 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 2736709..2736987 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 2752262..2752660 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 2752848..2752994 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 2775272..2775565 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 2786575..2786898 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 2800846..2801145 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 2808083..2808304 | + | antitoxin VapB38 |
| FN | CDS | NC_000962.3 | 2854683..2854907 | - | antitoxin VapB39 |
| FN | CDS | NC_000962.3 | 2943376..2943603 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 2962470..2962712 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 2968533..2968850 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 2980963..2981190 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 2981187..2981576 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 2981482..2981754 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 3020200..3020457 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 3034635..3034883 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 3042542..3043018 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 3055515..3056348 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 3076367..3076720 | - | membrane protein |
| FN | CDS | NC_000962.3 | 3082352..3082756 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 3112465..3113094 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 3147959..3148258 | - | hypothetical protein |
| FN | CDS | NC_000962.3 | 3174747..3174995 | + | antitoxin VapB23 |
| FN | CDS | NC_000962.3 | 3183138..3183395 | + | antitoxin VapB43 |
| FN | CDS | NC_000962.3 | 3200266..3201015 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 3242970..3243281 | + | hypothetical protein |
| FN | CDS | NC_000962.3 | 3289705..3290235 | + | transposase |
| FP | CDS | NC_000962.3 | 9963..10160 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 10523..10828 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 11421..11528 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 11555..11692 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 13024..13176 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 31824..31967 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 99379..99573 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 148677..148859 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 180805..181029 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 223607..223909 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 226700..226819 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 257783..258856 | - | acyl-CoA dehydrogenase family protein |
| FP | CDS | NC_000962.3 | 274710..274904 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 339107..339367 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 399487..400098 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 418901..419104 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 466672..467406 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 467459..468001 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 472914..473105 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 544736..545059 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 545188..545418 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 591475..591576 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 610935..611051 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 612255..612605 | + | class I SAM-dependent methyltransferase |
| FP | CDS | NC_000962.3 | 612598..612903 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 664923..665045 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 665434..665601 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 688032..688859 | + | MCE family protein |
| FP | CDS | NC_000962.3 | 688808..689062 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 697020..697193 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 702016..702759 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 704475..704672 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 711536..712231 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 712249..712719 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 759294..759524 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 767428..767541 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 789172..789420 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 799264..800019 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 832534..832848 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 835059..835205 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 850342..850512 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 874233..874943 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 874882..876390 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 876342..876503 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 889072..889398 | + | transposase |
| FP | CDS | NC_000962.3 | 889449..890333 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 908015..908113 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 931578..931991 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 947312..947644 | + | transposase |
| FP | CDS | NC_000962.3 | 986674..986823 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 993435..993752 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 1014397..1014555 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 1016820..1017050 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 1023865..1024122 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 1075111..1075251 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 1114173..1114382 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 1148080..1148172 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 1151777..1151923 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 1158918..1159409 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 1159375..1160061 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 1160095..1160418 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 1164164..1164436 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 1172110..1172313 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 1174995..1175621 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 1176011..1176265 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 1176372..1176485 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 1178604..1178750 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 1192235..1192375 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 1214049..1214360 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 1214826..1215131 | + | PE family protein |
| FP | CDS | NC_000962.3 | 1215599..1216054 | + | glycosyl hydrolase family 12 |
| FP | CDS | NC_000962.3 | 1231301..1231990 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 1232311..1232826 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 1241115..1241390 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 1242703..1242912 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 1264314..1264556 | + | acetyl-CoA acetyltransferase |
| FP | CDS | NC_000962.3 | 1264606..1264947 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 1266990..1267439 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 1272036..1272167 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 1278383..1278820 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 1282890..1283036 | + | DUF1059 domain-containing protein |
| FP | CDS | NC_000962.3 | 1301248..1301355 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 1351064..1351309 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 1376230..1376526 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 1451779..1452003 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 1526236..1526376 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 1530973..1531209 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 1541994..1542878 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 1542929..1543255 | - | transposase |
| FP | CDS | NC_000962.3 | 1612198..1612872 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 1616365..1616592 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 1617565..1617729 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 1622768..1623007 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 1623124..1623237 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 1693996..1694436 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 1694545..1695108 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 1701098..1701346 | + | hypothetical protein |
| FP | CDS | NC_000962.3 | 1712057..1712257 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 1746958..1747191 | - | hypothetical protein |
| FP | CDS | NC_000962.3 | 1753606..1754037 | + | hypothetical protein |

> Primeros pares con coordenadas no exactas

| tipo | seqid | ref | Δstart | Δend | solape |
|---|---|---|---|---|---|
| CDS | NC_000962.3 | 4434..4997 | +48 | +0 | 0.9149 |
| CDS | NC_000962.3 | 5240..7267 | -117 | +0 | 0.9455 |
| CDS | NC_000962.3 | 13133..13558 | +0 | -66 | 0.8451 |
| CDS | NC_000962.3 | 14089..14877 | +45 | +0 | 0.9430 |
| CDS | NC_000962.3 | 21637..23181 | +0 | -9 | 0.9942 |
| CDS | NC_000962.3 | 31514..31819 | +78 | +0 | 0.7451 |
| CDS | NC_000962.3 | 32057..33154 | +222 | +0 | 0.7978 |
| CDS | NC_000962.3 | 39877..41202 | +0 | -15 | 0.9887 |
| CDS | NC_000962.3 | 41304..41912 | -24 | +0 | 0.9621 |
| CDS | NC_000962.3 | 43562..46471 | +33 | +0 | 0.9887 |
| CDS | NC_000962.3 | 53663..55699 | -426 | +0 | 0.8270 |
| CDS | NC_000962.3 | 57410..57973 | -75 | +0 | 0.8826 |
| CDS | NC_000962.3 | 65012..65350 | +0 | +315 | 0.5183 |
| CDS | NC_000962.3 | 65552..66694 | +219 | +0 | 0.8084 |
| CDS | NC_000962.3 | 68418..68579 | +6 | +0 | 0.9630 |
| CDS | NC_000962.3 | 95414..96892 | +3 | +0 | 0.9980 |
| CDS | NC_000962.3 | 96927..97601 | +195 | +0 | 0.7111 |
| CDS | NC_000962.3 | 102815..103663 | +0 | -138 | 0.8375 |
| CDS | NC_000962.3 | 103710..104663 | +0 | +6 | 0.9938 |
| CDS | NC_000962.3 | 104805..105215 | +0 | -45 | 0.8905 |
| CDS | NC_000962.3 | 108156..109778 | +24 | +0 | 0.9852 |
| CDS | NC_000962.3 | 109783..110019 | -12 | +0 | 0.9518 |
| CDS | NC_000962.3 | 117714..119699 | +63 | +0 | 0.9683 |
| CDS | NC_000962.3 | 130895..131104 | +0 | +33 | 0.8642 |
| CDS | NC_000962.3 | 137319..137909 | -33 | +0 | 0.9471 |
| CDS | NC_000962.3 | 144049..145626 | -51 | +0 | 0.9687 |
| CDS | NC_000962.3 | 173238..174758 | +48 | +0 | 0.9684 |
| CDS | NC_000962.3 | 179319..180896 | +0 | -105 | 0.9335 |
| CDS | NC_000962.3 | 186495..186623 | +0 | -3 | 0.9767 |
| CDS | NC_000962.3 | 190607..191956 | -24 | +0 | 0.9825 |
| CDS | NC_000962.3 | 191984..193135 | +0 | -33 | 0.9714 |
| CDS | NC_000962.3 | 193626..194111 | +45 | +0 | 0.9074 |
| CDS | NC_000962.3 | 194144..194815 | +0 | +123 | 0.8453 |
| CDS | NC_000962.3 | 206814..207455 | -66 | +0 | 0.9068 |
| CDS | NC_000962.3 | 213028..214140 | +0 | -114 | 0.8976 |
| CDS | NC_000962.3 | 214088..214927 | -57 | +0 | 0.9365 |
| CDS | NC_000962.3 | 218390..218551 | +0 | -15 | 0.9074 |
| CDS | NC_000962.3 | 218705..219367 | +6 | +0 | 0.9910 |
| CDS | NC_000962.3 | 223564..224664 | +336 | +0 | 0.6948 |
| CDS | NC_000962.3 | 230899..231534 | -132 | +0 | 0.8281 |
| CDS | NC_000962.3 | 238392..241292 | +0 | -3 | 0.9990 |
| CDS | NC_000962.3 | 241976..243214 | +0 | +45 | 0.9650 |
| CDS | NC_000962.3 | 247384..248112 | +0 | +6 | 0.9918 |
| CDS | NC_000962.3 | 248115..248906 | +0 | -198 | 0.7500 |
| CDS | NC_000962.3 | 254637..255950 | +0 | +3 | 0.9977 |
| CDS | NC_000962.3 | 256064..257677 | +96 | +0 | 0.9405 |
| CDS | NC_000962.3 | 271574..272839 | +0 | +15 | 0.9883 |
| CDS | NC_000962.3 | 273055..274278 | -75 | +0 | 0.9423 |
| CDS | NC_000962.3 | 279605..281140 | +0 | -162 | 0.8945 |
| CDS | NC_000962.3 | 288428..289042 | +48 | +0 | 0.9220 |
| CDS | NC_000962.3 | 302866..305427 | +51 | +0 | 0.9801 |
| CDS | NC_000962.3 | 315783..316415 | +0 | +21 | 0.9679 |
| CDS | NC_000962.3 | 322764..323273 | +0 | -252 | 0.5059 |
| CDS | NC_000962.3 | 330933..331658 | +0 | +9 | 0.9878 |
| CDS | NC_000962.3 | 356862..358064 | +0 | +9 | 0.9926 |
| CDS | NC_000962.3 | 383602..384486 | +30 | +0 | 0.9661 |
| CDS | NC_000962.3 | 385173..385943 | +0 | -33 | 0.9572 |
| CDS | NC_000962.3 | 386305..387099 | +0 | -18 | 0.9774 |
| CDS | NC_000962.3 | 398254..398628 | -60 | +0 | 0.8621 |
| CDS | NC_000962.3 | 400192..401703 | +69 | +0 | 0.9544 |
| CDS | NC_000962.3 | 405950..408448 | +0 | -48 | 0.9808 |
| CDS | NC_000962.3 | 415050..415460 | -6 | +0 | 0.9856 |
| CDS | NC_000962.3 | 417304..418290 | +174 | +0 | 0.8237 |
| CDS | NC_000962.3 | 418949..419608 | +243 | +0 | 0.6318 |
| CDS | NC_000962.3 | 438302..438739 | +0 | -33 | 0.9247 |
| CDS | NC_000962.3 | 445314..446525 | +0 | +9 | 0.9926 |
| CDS | NC_000962.3 | 446531..447046 | +0 | +87 | 0.8557 |
| CDS | NC_000962.3 | 452294..453154 | +0 | +42 | 0.9535 |
| CDS | NC_000962.3 | 455637..455858 | -3 | +0 | 0.9867 |
| CDS | NC_000962.3 | 456268..456819 | +0 | +99 | 0.8479 |
| CDS | NC_000962.3 | 462135..463307 | -105 | +0 | 0.9178 |
| CDS | NC_000962.3 | 469591..470013 | +42 | +0 | 0.9007 |
| CDS | NC_000962.3 | 472781..474106 | +339 | +0 | 0.7443 |
| CDS | NC_000962.3 | 475350..475742 | -24 | +0 | 0.9424 |
| CDS | NC_000962.3 | 483977..485734 | +552 | +0 | 0.6860 |
| CDS | NC_000962.3 | 499713..500366 | +60 | +0 | 0.9083 |
| CDS | NC_000962.3 | 506561..506971 | +0 | -69 | 0.8321 |
| CDS | NC_000962.3 | 507758..508555 | +0 | +30 | 0.9638 |
| CDS | NC_000962.3 | 519073..519567 | +21 | +0 | 0.9576 |
| CDS | NC_000962.3 | 520324..521454 | -21 | +0 | 0.9818 |
| CDS | NC_000962.3 | 526143..527360 | +0 | -39 | 0.9680 |
| CDS | NC_000962.3 | 530751..532214 | +0 | -51 | 0.9652 |
| CDS | NC_000962.3 | 532396..532911 | -42 | +0 | 0.9247 |
| CDS | NC_000962.3 | 534445..535215 | +0 | +12 | 0.9847 |
| CDS | NC_000962.3 | 536504..537169 | +0 | +66 | 0.9098 |
| CDS | NC_000962.3 | 537229..538548 | +0 | -24 | 0.9818 |
| CDS | NC_000962.3 | 541488..541910 | +0 | -6 | 0.9858 |
| CDS | NC_000962.3 | 542142..542852 | -51 | +0 | 0.9331 |
| CDS | NC_000962.3 | 547586..549607 | +0 | +9 | 0.9956 |
| CDS | NC_000962.3 | 551749..551988 | -39 | +0 | 0.8602 |
| CDS | NC_000962.3 | 552026..552550 | -45 | +0 | 0.9211 |
| CDS | NC_000962.3 | 561854..562294 | +10 | -146 | 0.6463 |
| CDS | NC_000962.3 | 562723..563427 | +0 | +66 | 0.9144 |
| CDS | NC_000962.3 | 567921..568967 | +0 | -297 | 0.7163 |
| CDS | NC_000962.3 | 575348..576790 | +9 | +0 | 0.9938 |
| CDS | NC_000962.3 | 577664..578269 | -21 | +0 | 0.9665 |
| CDS | NC_000962.3 | 578426..579175 | -21 | +0 | 0.9728 |
| CDS | NC_000962.3 | 583701..584690 | +0 | -39 | 0.9606 |
| CDS | NC_000962.3 | 584695..585423 | +12 | +0 | 0.9835 |
| CDS | NC_000962.3 | 585424..586314 | +0 | -120 | 0.8653 |


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
/opt/conda/bin/python benchmarks/external/external_benchmark.py --genome data/reference_genomes/mtb_h37rv/mtb_h37rv_genomic.fna.gz --reference data/reference_genomes/mtb_h37rv/mtb_h37rv_genomic.gff.gz --tools bakta,prokka,dfast --run --threads 1 --json reports/external_benchmark/benchmark.json --md reports/external_benchmark/benchmark.md --bakta-db /data/bakta_db
```

Para comparar salidas producidas en otro host, usar `--query-gff TOOL=RUTA` sin `--run`; el origen queda como `precomputed`.
