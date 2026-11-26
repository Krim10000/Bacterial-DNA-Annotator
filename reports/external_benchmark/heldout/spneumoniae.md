# Benchmark externo CPU: Bakta, Prokka y DFAST

**Generado:** 2026-11-26T02:43:57Z

**Esquema:** `external-benchmark-run/v1`

**Alcance:** Comparación local de CDS/tRNA/rRNA de Bakta, Prokka y DFAST contra una referencia GFF3/GBFF. CPU, sin descargas de datos.

**Comando reproducible:**

```bash
/opt/conda/bin/python /opt/benchmark/external/external_benchmark.py --genome /work/data/reference_genomes/spneumoniae/spneumoniae_genomic.fna.gz --reference /work/data/reference_genomes/spneumoniae/spneumoniae_genomic.gff.gz --tools all --run --threads 4 --timeout 5400 --bakta-db /data/bakta_db --extra-arg dfast=--dbroot --extra-arg dfast=/data/dfast_db --save-gffs /work/reports/external_benchmark/heldout/spneumoniae/evidence --json /work/reports/external_benchmark/heldout/spneumoniae.json --md /work/reports/external_benchmark/heldout/spneumoniae.md
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

- FASTA: `data/reference_genomes/spneumoniae/spneumoniae_genomic.fna.gz` (SHA256 `e5e580af13c7a8ff…`, 1 contig, 2110968 bp).
- Referencia: `data/reference_genomes/spneumoniae/spneumoniae_genomic.gff.gz` (formato `gff3`, SHA256 `5724f05787121edc…`).
- Conteos de referencia: CDS **1931**, tRNA **58**, rRNA **12**.
- CDS pseudo excluidos en referencia: **186**.
- SHA256 FASTA: `e5e580af13c7a8ff81e795fe656615ddc0495982e6afc6c7284a7df8c08a67c6`.
- SHA256 referencia: `5724f05787121edc09f41212753e002d930258c4a35b76f0cf1843204ba607ae`.

## Validación

- Inicio: `2026-11-26T02:22:33Z`; fin: `2026-11-26T02:43:57Z`; total: **1283.54 s**.
- JSON declarado: `reports/external_benchmark/heldout/spneumoniae.json`.
- Markdown declarado: `reports/external_benchmark/heldout/spneumoniae.md`.

| check | estado | detalle |
|---|---|---|
| genome_path | ok | FASTA existe y es legible |
| genome_sha256 | ok | SHA256 e5e580af13c7a8ff81e795fe656615ddc0495982e6afc6c7284a7df8c08a67c6 |
| reference_path | ok | referencia gff3 con 2001 features |
| reference_sha256 | ok | SHA256 5724f05787121edc09f41212753e002d930258c4a35b76f0cf1843204ba607ae |
| reference_features | ok | CDS=1931 tRNA=58 rRNA=12 pseudo_excluidos=186 |
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
| Bakta | executed | ok | 0.9180 | 1.0000 | 1.0000 | 0.8866 | 0.9575 | 0.9207 | 556.9 |
| Prokka | executed | ok | 0.9180 | 1.0000 | 1.0000 | 0.8866 | 0.9575 | 0.9207 | 420.1 |
| DFAST | executed | ok | 0.9180 | 1.0000 | 1.0000 | 0.8866 | 0.9575 | 0.9207 | 283.6 |

La precisión/recall/F1 se calculan contra la referencia local. Si una herramienta no está instalada, su estado es `not_available` y no hay métricas que mostrar.

## Bakta

- Estado: **ok**.
- Fuente: `executed`.
- Bakta DB: `ok` (válida: **sí**, ruta `/data/bakta_db`, versión `0`).
- Versión: `bakta 1.12.2` (estado `ok`).
- Ejecución: **ok**, 556.86 s, RSS máx. 2652.1 MiB (`2026-11-26T02:22:39Z` → `2026-11-26T02:31:55Z`).
- Salida: `/tmp/external_benchmark_bakta_6p7zn2r7/out/benchmark.gff3` (2598267 bytes, SHA256 `c03996e15821ac42…`).
- Comando: `/usr/local/bin/bakta --output /tmp/external_benchmark_bakta_6p7zn2r7/out --prefix benchmark --threads 4 --force --skip-plot --db /data/bakta_db /tmp/benchmark_genome_otq2w3ys.fna`.
- Anotación comparada: `reports/external_benchmark/heldout/spneumoniae/evidence/bakta/benchmark.gff3` (formato `gff3`, SHA256 `c03996e15821ac42…`).

### Features (Bakta)

| tipo | ref | query | TP | FP | FN | precisión | recall | F1 | exactos | borde distinto | solape medio |
|---|---|---|---|---|---|---|---|---|---|---|---|
| CDS | 1931 | 2091 | 1846 | 245 | 85 | 0.8828 | 0.9560 | 0.9180 | 1764 | 82 | 0.9967 |
| tRNA | 58 | 58 | 58 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 58 | 0 | 1.0000 |
| rRNA | 12 | 12 | 12 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 12 | 0 | 1.0000 |

### Atributos (Bakta)

| tipo | pares | product en ambos | product literal | product tokens | gene literal | locus literal |
|---|---|---|---|---|---|---|
| CDS | 1846 | 1846 | 0.2541 | 0.4592 | 0.9050 | 0.0000 |
| tRNA | 58 | 58 | 0.0000 | 0.6451 | — | 0.0000 |
| rRNA | 12 | 12 | 1.0000 | 1.0000 | 1.0000 | 0.0000 |

### Hipotéticas (Bakta)

| tipo | ref CDS | ref con product | ref hyp | % ref | query CDS | query con product | query hyp | % query |
|---|---|---|---|---|---|---|---|---|
| CDS | 1931 | 1931 | 268 | 13.9 | 2091 | 2091 | 339 | 16.2 |
| tRNA | 0 | 0 | 0 | — | 0 | 0 | 0 | — |
| rRNA | 0 | 0 | 0 | — | 0 | 0 | 0 | — |

#### Diferencias

- Falsos negativos (referencia sin pareja): **85**.
- Falsos positivos (herramienta sin pareja): **245**.
- Pares solo con hebra distinta: **1**.
- Pares emparejados con algún borde distinto: **82**.

Las listas de diferencias están truncadas a 100 elementos por lista; los conteos son completos.

> Primeras diferencias no emparejadas

| tipo | feature | seqid | coordenadas | hebra | product/locus |
|---|---|---|---|---|---|
| FN | CDS | NZ_LN831051.1 | 1588..2076 | + | phosphoglycerate kinase |
| FN | CDS | NZ_LN831051.1 | 5195..5314 | - | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 7280..7417 | + | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 11755..11892 | - | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 13858..13977 | + | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 26991..27116 | - | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 52213..52335 | + | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 76489..76674 | - | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 78278..78514 | - | CPBP family intramembrane glutamic endopeptidase |
| FN | CDS | NZ_LN831051.1 | 115640..115978 | + | lanthionine synthetase LanC family protein |
| FN | CDS | NZ_LN831051.1 | 238743..238835 | + | type I toxin-antitoxin system Fst family toxin |
| FN | CDS | NZ_LN831051.1 | 253201..253398 | + | helix-turn-helix domain-containing protein |
| FN | CDS | NZ_LN831051.1 | 277204..277344 | - | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 316038..316190 | - | AzlD domain-containing protein |
| FN | CDS | NZ_LN831051.1 | 325239..325631 | - | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 328845..329030 | - | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 337699..337857 | - | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 344487..344726 | - | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 347913..348095 | - | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 351359..351595 | - | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 361731..361898 | + | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 386482..386634 | + | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 405670..405786 | - | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 416534..416644 | - | phosphoribosylaminoimidazolesuccinocarboxamide synthase |
| FN | CDS | NZ_LN831051.1 | 421632..421757 | - | PncF family bacteriocin immunity protein |
| FN | CDS | NZ_LN831051.1 | 463328..463537 | - | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 556922..557119 | + | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 561107..561238 | + | teichoic acid D-Ala incorporation-associated protein DltX |
| FN | CDS | NZ_LN831051.1 | 641113..641319 | + | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 657535..657669 | + | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 675941..676120 | + | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 707936..708067 | + | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 740073..740225 | + | 50S ribosomal protein L33 |
| FN | CDS | NZ_LN831051.1 | 810117..810278 | + | helix-turn-helix domain-containing protein |
| FN | CDS | NZ_LN831051.1 | 824247..824438 | + | cell wall-binding repeat protein |
| FN | CDS | NZ_LN831051.1 | 824463..824525 | + | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 876879..876989 | - | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 917884..918072 | + | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 953070..953249 | + | type II toxin-antitoxin system HicA family toxin |
| FN | CDS | NZ_LN831051.1 | 970783..970926 | + | IS630 transposase-related protein |
| FN | CDS | NZ_LN831051.1 | 971597..976975 | + | accessory Sec-dependent serine-rich glycoprotein adhesin |
| FN | CDS | NZ_LN831051.1 | 1077873..1078064 | + | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 1098064..1098228 | + | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 1124411..1124587 | - | MazG-like protein |
| FN | CDS | NZ_LN831051.1 | 1134273..1134419 | + | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 1155403..1155549 | - | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 1232933..1233226 | + | IS630 transposase-related protein |
| FN | CDS | NZ_LN831051.1 | 1233868..1234110 | + | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 1273357..1273551 | + | IS66 family transposase |
| FN | CDS | NZ_LN831051.1 | 1286893..1287204 | + | CHY zinc finger protein |
| FN | CDS | NZ_LN831051.1 | 1295430..1295792 | + | G5 domain-containing protein |
| FN | CDS | NZ_LN831051.1 | 1322554..1322793 | + | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 1330351..1330611 | + | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 1416356..1416478 | - | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 1452346..1452471 | - | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 1476860..1477012 | - | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 1516001..1516225 | + | HAD family hydrolase |
| FN | CDS | NZ_LN831051.1 | 1519503..1519949 | - | tyrosine-type recombinase/integrase |
| FN | CDS | NZ_LN831051.1 | 1577534..1577713 | - | DUF3884 family protein |
| FN | CDS | NZ_LN831051.1 | 1619020..1619091 | - | SHP2/SHP3 family peptide pheromone |
| FN | CDS | NZ_LN831051.1 | 1635940..1636383 | - | restriction endonuclease subunit S |
| FN | CDS | NZ_LN831051.1 | 1673056..1673256 | - | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 1693058..1693207 | - | 50S ribosomal protein L33 |
| FN | CDS | NZ_LN831051.1 | 1707418..1707606 | - | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 1719152..1719424 | - | nucleotide-binding domain-containing protein |
| FN | CDS | NZ_LN831051.1 | 1732075..1732272 | + | AraC family transcriptional regulator |
| FN | CDS | NZ_LN831051.1 | 1751301..1751441 | + | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 1757310..1757558 | + | site-specific integrase |
| FN | CDS | NZ_LN831051.1 | 1773388..1773513 | - | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 1869636..1869782 | + | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 1874117..1874326 | + | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 1931896..1932018 | - | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 1981566..1981712 | - | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 1993132..1993314 | + | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 1997131..1997376 | - | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 1997775..1997984 | - | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 1998150..1998296 | + | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 2000396..2000533 | + | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 2001370..2002174 | - | IS5 family transposase |
| FN | CDS | NZ_LN831051.1 | 2028442..2028525 | - | putative autolysis signal peptide Pep27 |
| FN | CDS | NZ_LN831051.1 | 2068225..2068377 | + | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 2080905..2081024 | - | PncF family bacteriocin immunity protein |
| FN | CDS | NZ_LN831051.1 | 2081128..2081277 | - | bacteriocin-like peptide BlpO |
| FN | CDS | NZ_LN831051.1 | 2107831..2107968 | - | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 2109956..2110054 | + | glutamine synthetase |
| FP | CDS | NZ_LN831051.1 | 558..758 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 991..1263 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 54689..54811 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 55009..55557 | + | Transposase_20 domain protein |
| FP | CDS | NZ_LN831051.1 | 69576..70061 | - | Glycosyltransferase family 8 protein |
| FP | CDS | NZ_LN831051.1 | 106480..106674 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 106686..107156 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 114648..114977 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 117857..118105 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 118056..118799 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 136563..136874 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 137171..137527 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 137615..137923 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 140716..141498 | + | Transposase |
| FP | CDS | NZ_LN831051.1 | 141432..141614 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 141672..142022 | + | Transposase IS204/IS1001/IS1096/IS1165 DDE domain-containing protein |
| FP | CDS | NZ_LN831051.1 | 154662..155060 | + | Transposase Synechocystis PCC 6803 domain-containing protein |
| FP | CDS | NZ_LN831051.1 | 155188..155400 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 155499..156281 | + | Transposase |
| FP | CDS | NZ_LN831051.1 | 156215..156397 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 156455..156805 | + | Transposase IS204/IS1001/IS1096/IS1165 DDE domain-containing protein |
| FP | CDS | NZ_LN831051.1 | 161196..161546 | - | Group II intron reverse transcriptase/maturase |
| FP | CDS | NZ_LN831051.1 | 164179..164268 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 181656..181871 | - | Transposase IS200-like domain-containing protein |
| FP | CDS | NZ_LN831051.1 | 189038..189163 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 227465..228421 | - | Alkanal monooxygenase alpha chain |
| FP | CDS | NZ_LN831051.1 | 253615..253716 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 254029..254223 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 254287..254736 | + | Iron ABC transporter iron-binding protein |
| FP | CDS | NZ_LN831051.1 | 256088..256456 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 256453..257460 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 281258..281401 | - | Transposase |
| FP | CDS | NZ_LN831051.1 | 291487..291720 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 291860..292432 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 316038..316652 | - | AzlC family protein |
| FP | CDS | NZ_LN831051.1 | 325072..325209 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 325682..325915 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 327654..328514 | - | Glycosyltransferase |
| FP | CDS | NZ_LN831051.1 | 330982..331911 | + | Transposase |
| FP | CDS | NZ_LN831051.1 | 331902..332288 | + | Transposase IS204/IS1001/IS1096/IS1165 DDE domain-containing protein |
| FP | CDS | NZ_LN831051.1 | 335381..335488 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 349014..349829 | - | Leader peptidase |
| FP | CDS | NZ_LN831051.1 | 352124..352288 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 353330..353749 | - | Argininosuccinate synthase, truncation |
| FP | CDS | NZ_LN831051.1 | 422193..422486 | - | Mobile genetic element |
| FP | CDS | NZ_LN831051.1 | 430719..430976 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 431300..431533 | + | Adenylosuccinate lyase |
| FP | CDS | NZ_LN831051.1 | 431660..432139 | + | Degenerate transposase (Orf1) |
| FP | CDS | NZ_LN831051.1 | 433188..433772 | - | Repair protein |
| FP | CDS | NZ_LN831051.1 | 462098..462451 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 462487..462612 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 463721..463825 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 465018..465176 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 470467..470625 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 476619..476948 | - | Transposase |
| FP | CDS | NZ_LN831051.1 | 503511..504278 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 504363..504701 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 526186..526455 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 620828..622057 | + | Ribosomal protein L11 methyltransferase |
| FP | CDS | NZ_LN831051.1 | 639777..639923 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 657230..657610 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 657784..658209 | - | Transposase DDE domain-containing protein |
| FP | CDS | NZ_LN831051.1 | 664650..664952 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 665290..665679 | + | Transposase |
| FP | CDS | NZ_LN831051.1 | 705086..705268 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 733741..734091 | + | Transposase |
| FP | CDS | NZ_LN831051.1 | 759147..760040 | + | Gram-positive cocci surface proteins LPxTG domain-containing protein |
| FP | CDS | NZ_LN831051.1 | 762918..764600 | + | Bacteriocin-associated integral membrane family protein |
| FP | CDS | NZ_LN831051.1 | 764720..764947 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 765658..765876 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 768468..768854 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 795549..795638 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 799283..799471 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 799482..799586 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 818743..818871 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 818881..819459 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 819583..819813 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 819823..820032 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 829575..829904 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 830004..830330 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 840476..840913 | - | Transposase IS204/IS1001/IS1096/IS1165 DDE domain-containing protein |
| FP | CDS | NZ_LN831051.1 | 859001..859504 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 859619..859894 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 859939..860379 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 860388..860816 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 899345..899572 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 899526..899801 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 925005..925145 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 925182..925736 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 955587..955745 | + | 7,8-dihydro-8-oxoguanine-triphosphatase |
| FP | CDS | NZ_LN831051.1 | 959185..959628 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 959638..960204 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 970371..970514 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 972928..976980 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 980105..981178 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 981399..981734 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 991213..991395 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 991399..992115 | + | Glycosyl transferase, family 2 |
| FP | CDS | NZ_LN831051.1 | 1003775..1004734 | - | Serine/threonine transporter SstT |
| FP | CDS | NZ_LN831051.1 | 1024890..1025279 | + | Transposase |

> Primeros pares con coordenadas no exactas

| tipo | seqid | ref | Δstart | Δend | solape |
|---|---|---|---|---|---|
| CDS | NZ_LN831051.1 | 1991..2302 | -33 | +0 | 0.9043 |
| CDS | NZ_LN831051.1 | 9733..10338 | +60 | +0 | 0.9010 |
| CDS | NZ_LN831051.1 | 78527..79000 | +0 | -39 | 0.9177 |
| CDS | NZ_LN831051.1 | 105789..106127 | +0 | +54 | 0.8626 |
| CDS | NZ_LN831051.1 | 192006..193985 | +0 | -6 | 0.9970 |
| CDS | NZ_LN831051.1 | 201815..202789 | +0 | -30 | 0.9692 |
| CDS | NZ_LN831051.1 | 237611..238609 | +0 | -27 | 0.9730 |
| CDS | NZ_LN831051.1 | 307879..309549 | +0 | +12 | 0.9929 |
| CDS | NZ_LN831051.1 | 323503..323931 | +0 | +171 | 0.7150 |
| CDS | NZ_LN831051.1 | 363902..365041 | +0 | -93 | 0.9184 |
| CDS | NZ_LN831051.1 | 392700..393104 | +0 | -66 | 0.8370 |
| CDS | NZ_LN831051.1 | 427737..430370 | +0 | +36 | 0.9865 |
| CDS | NZ_LN831051.1 | 434582..435943 | +0 | -99 | 0.9273 |
| CDS | NZ_LN831051.1 | 435957..436472 | +0 | +42 | 0.9247 |
| CDS | NZ_LN831051.1 | 437853..438002 | -6 | +0 | 0.9615 |
| CDS | NZ_LN831051.1 | 440964..444182 | +0 | -48 | 0.9851 |
| CDS | NZ_LN831051.1 | 459977..460519 | +0 | +48 | 0.9188 |
| CDS | NZ_LN831051.1 | 461321..461479 | +19 | +4 | 0.8805 |
| CDS | NZ_LN831051.1 | 463027..463356 | +0 | -3 | 0.9909 |
| CDS | NZ_LN831051.1 | 465164..465991 | +0 | +30 | 0.9650 |
| CDS | NZ_LN831051.1 | 473252..474538 | +0 | +42 | 0.9684 |
| CDS | NZ_LN831051.1 | 497386..498144 | +0 | -15 | 0.9802 |
| CDS | NZ_LN831051.1 | 549958..551385 | +18 | +0 | 0.9874 |
| CDS | NZ_LN831051.1 | 556182..556790 | +0 | +51 | 0.9227 |
| CDS | NZ_LN831051.1 | 588251..589135 | -15 | +0 | 0.9833 |
| CDS | NZ_LN831051.1 | 605967..606836 | -12 | +0 | 0.9864 |
| CDS | NZ_LN831051.1 | 623548..624012 | +150 | +0 | 0.6774 |
| CDS | NZ_LN831051.1 | 624022..624228 | -9 | +0 | 0.9583 |
| CDS | NZ_LN831051.1 | 629101..629382 | +0 | +75 | 0.7899 |
| CDS | NZ_LN831051.1 | 633828..635081 | +0 | +69 | 0.9478 |
| CDS | NZ_LN831051.1 | 702128..703144 | +144 | +0 | 0.8584 |
| CDS | NZ_LN831051.1 | 708430..709254 | -6 | +0 | 0.9928 |
| CDS | NZ_LN831051.1 | 748284..748490 | +15 | +0 | 0.9275 |
| CDS | NZ_LN831051.1 | 774575..775825 | -894 | +0 | 0.5832 |
| CDS | NZ_LN831051.1 | 782878..783108 | +33 | +0 | 0.8571 |
| CDS | NZ_LN831051.1 | 784408..784596 | -6 | +0 | 0.9692 |
| CDS | NZ_LN831051.1 | 820935..821918 | -48 | +0 | 0.9535 |
| CDS | NZ_LN831051.1 | 893522..894328 | +36 | +0 | 0.9554 |
| CDS | NZ_LN831051.1 | 899954..900250 | +51 | +0 | 0.8283 |
| CDS | NZ_LN831051.1 | 907206..909365 | +33 | +0 | 0.9847 |
| CDS | NZ_LN831051.1 | 956951..959134 | -60 | +0 | 0.9733 |
| CDS | NZ_LN831051.1 | 981898..982218 | -42 | +0 | 0.8843 |
| CDS | NZ_LN831051.1 | 1011074..1011814 | -6 | +0 | 0.9920 |
| CDS | NZ_LN831051.1 | 1025347..1025595 | -39 | +0 | 0.8646 |
| CDS | NZ_LN831051.1 | 1045290..1045625 | -9 | +0 | 0.9739 |
| CDS | NZ_LN831051.1 | 1054009..1054371 | -6 | +0 | 0.9837 |
| CDS | NZ_LN831051.1 | 1074967..1075833 | -21 | +0 | 0.9764 |
| CDS | NZ_LN831051.1 | 1077585..1077863 | -54 | +0 | 0.8378 |
| CDS | NZ_LN831051.1 | 1081851..1083908 | +15 | +0 | 0.9927 |
| CDS | NZ_LN831051.1 | 1135543..1136508 | -18 | +0 | 0.9817 |
| CDS | NZ_LN831051.1 | 1136559..1136765 | +0 | +6 | 0.9718 |
| CDS | NZ_LN831051.1 | 1187212..1187700 | +57 | +0 | 0.8834 |
| CDS | NZ_LN831051.1 | 1243754..1245145 | -51 | +0 | 0.9647 |
| CDS | NZ_LN831051.1 | 1288573..1289379 | +87 | +0 | 0.8922 |
| CDS | NZ_LN831051.1 | 1315850..1316581 | +60 | +0 | 0.9180 |
| CDS | NZ_LN831051.1 | 1341610..1342455 | +9 | +0 | 0.9894 |
| CDS | NZ_LN831051.1 | 1350369..1350740 | -60 | +0 | 0.8611 |
| CDS | NZ_LN831051.1 | 1355496..1356755 | -162 | +0 | 0.8861 |
| CDS | NZ_LN831051.1 | 1383594..1384172 | +33 | +0 | 0.9430 |
| CDS | NZ_LN831051.1 | 1427900..1429216 | -27 | +0 | 0.9799 |
| CDS | NZ_LN831051.1 | 1456991..1457245 | -9 | +0 | 0.9659 |
| CDS | NZ_LN831051.1 | 1472042..1472473 | +39 | +0 | 0.9097 |
| CDS | NZ_LN831051.1 | 1492688..1493287 | +21 | +0 | 0.9650 |
| CDS | NZ_LN831051.1 | 1570351..1572084 | +0 | -9 | 0.9948 |
| CDS | NZ_LN831051.1 | 1618047..1618916 | +0 | -39 | 0.9552 |
| CDS | NZ_LN831051.1 | 1699481..1701646 | +0 | -57 | 0.9737 |
| CDS | NZ_LN831051.1 | 1701688..1702674 | +0 | -48 | 0.9514 |
| CDS | NZ_LN831051.1 | 1704769..1705299 | +0 | +27 | 0.9516 |
| CDS | NZ_LN831051.1 | 1743852..1744526 | +84 | +0 | 0.8756 |
| CDS | NZ_LN831051.1 | 1777895..1779037 | +0 | +60 | 0.9501 |
| CDS | NZ_LN831051.1 | 1779648..1781900 | +0 | +51 | 0.9779 |
| CDS | NZ_LN831051.1 | 1803732..1806203 | +0 | -36 | 0.9854 |
| CDS | NZ_LN831051.1 | 1845314..1845877 | +0 | +9 | 0.9843 |
| CDS | NZ_LN831051.1 | 1891134..1892229 | +0 | -118 | 0.8923 |
| CDS | NZ_LN831051.1 | 1929576..1930274 | +0 | -21 | 0.9700 |
| CDS | NZ_LN831051.1 | 1940539..1940838 | +0 | +48 | 0.8621 |
| CDS | NZ_LN831051.1 | 1945275..1946183 | -60 | +0 | 0.9381 |
| CDS | NZ_LN831051.1 | 1974083..1974529 | -177 | +0 | 0.7163 |
| CDS | NZ_LN831051.1 | 1993814..1994092 | +0 | +3 | 0.9894 |
| CDS | NZ_LN831051.1 | 2032743..2033108 | +0 | +12 | 0.9683 |
| CDS | NZ_LN831051.1 | 2043911..2044177 | +0 | -69 | 0.7416 |
| CDS | NZ_LN831051.1 | 2073903..2075039 | +0 | -12 | 0.9894 |

- Aviso: 2 feature(s) GFF multi-segmento fusionados por ID; se comparan por su intervalo envolvente.

## Prokka

- Estado: **ok**.
- Fuente: `executed`.
- Versión: `prokka 1.15.6` (estado `ok`).
- Ejecución: **ok**, 420.14 s, RSS máx. 2989.6 MiB (`2026-11-26T02:32:01Z` → `2026-11-26T02:39:01Z`).
- Salida: `/tmp/external_benchmark_prokka_2ou8quep/out/benchmark.gff` (2642554 bytes, SHA256 `adfc51a42afd9848…`).
- Comando: `/usr/local/bin/prokka --outdir /tmp/external_benchmark_prokka_2ou8quep/out --prefix benchmark --cpus 4 --force --quiet /tmp/benchmark_genome_otq2w3ys.fna`.
- Anotación comparada: `reports/external_benchmark/heldout/spneumoniae/evidence/prokka/benchmark.gff` (formato `gff3`, SHA256 `adfc51a42afd9848…`).

### Features (Prokka)

| tipo | ref | query | TP | FP | FN | precisión | recall | F1 | exactos | borde distinto | solape medio |
|---|---|---|---|---|---|---|---|---|---|---|---|
| CDS | 1931 | 2091 | 1846 | 245 | 85 | 0.8828 | 0.9560 | 0.9180 | 1764 | 82 | 0.9967 |
| tRNA | 58 | 58 | 58 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 35 | 23 | 0.9892 |
| rRNA | 12 | 12 | 12 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 12 | 0 | 1.0000 |

### Atributos (Prokka)

| tipo | pares | product en ambos | product literal | product tokens | gene literal | locus literal |
|---|---|---|---|---|---|---|
| CDS | 1846 | 1846 | 0.2010 | 0.4034 | 0.7103 | 0.0000 |
| tRNA | 58 | 58 | 0.0000 | 0.6595 | — | 0.0000 |
| rRNA | 12 | 12 | 1.0000 | 1.0000 | — | 0.0000 |

### Hipotéticas (Prokka)

| tipo | ref CDS | ref con product | ref hyp | % ref | query CDS | query con product | query hyp | % query |
|---|---|---|---|---|---|---|---|---|
| CDS | 1931 | 1931 | 268 | 13.9 | 2091 | 2091 | 786 | 37.6 |
| tRNA | 0 | 0 | 0 | — | 0 | 0 | 0 | — |
| rRNA | 0 | 0 | 0 | — | 0 | 0 | 0 | — |

#### Diferencias

- Falsos negativos (referencia sin pareja): **85**.
- Falsos positivos (herramienta sin pareja): **245**.
- Pares solo con hebra distinta: **1**.
- Pares emparejados con algún borde distinto: **105**.

Las listas de diferencias están truncadas a 100 elementos por lista; los conteos son completos.

> Primeras diferencias no emparejadas

| tipo | feature | seqid | coordenadas | hebra | product/locus |
|---|---|---|---|---|---|
| FN | CDS | NZ_LN831051.1 | 1588..2076 | + | phosphoglycerate kinase |
| FN | CDS | NZ_LN831051.1 | 5195..5314 | - | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 7280..7417 | + | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 11755..11892 | - | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 13858..13977 | + | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 26991..27116 | - | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 52213..52335 | + | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 76489..76674 | - | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 78278..78514 | - | CPBP family intramembrane glutamic endopeptidase |
| FN | CDS | NZ_LN831051.1 | 115640..115978 | + | lanthionine synthetase LanC family protein |
| FN | CDS | NZ_LN831051.1 | 238743..238835 | + | type I toxin-antitoxin system Fst family toxin |
| FN | CDS | NZ_LN831051.1 | 253201..253398 | + | helix-turn-helix domain-containing protein |
| FN | CDS | NZ_LN831051.1 | 277204..277344 | - | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 316038..316190 | - | AzlD domain-containing protein |
| FN | CDS | NZ_LN831051.1 | 325239..325631 | - | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 328845..329030 | - | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 337699..337857 | - | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 344487..344726 | - | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 347913..348095 | - | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 351359..351595 | - | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 361731..361898 | + | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 386482..386634 | + | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 405670..405786 | - | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 416534..416644 | - | phosphoribosylaminoimidazolesuccinocarboxamide synthase |
| FN | CDS | NZ_LN831051.1 | 421632..421757 | - | PncF family bacteriocin immunity protein |
| FN | CDS | NZ_LN831051.1 | 463328..463537 | - | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 556922..557119 | + | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 561107..561238 | + | teichoic acid D-Ala incorporation-associated protein DltX |
| FN | CDS | NZ_LN831051.1 | 641113..641319 | + | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 657535..657669 | + | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 675941..676120 | + | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 707936..708067 | + | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 740073..740225 | + | 50S ribosomal protein L33 |
| FN | CDS | NZ_LN831051.1 | 810117..810278 | + | helix-turn-helix domain-containing protein |
| FN | CDS | NZ_LN831051.1 | 824247..824438 | + | cell wall-binding repeat protein |
| FN | CDS | NZ_LN831051.1 | 824463..824525 | + | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 876879..876989 | - | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 917884..918072 | + | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 953070..953249 | + | type II toxin-antitoxin system HicA family toxin |
| FN | CDS | NZ_LN831051.1 | 970783..970926 | + | IS630 transposase-related protein |
| FN | CDS | NZ_LN831051.1 | 971597..976975 | + | accessory Sec-dependent serine-rich glycoprotein adhesin |
| FN | CDS | NZ_LN831051.1 | 1077873..1078064 | + | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 1098064..1098228 | + | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 1124411..1124587 | - | MazG-like protein |
| FN | CDS | NZ_LN831051.1 | 1134273..1134419 | + | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 1155403..1155549 | - | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 1232933..1233226 | + | IS630 transposase-related protein |
| FN | CDS | NZ_LN831051.1 | 1233868..1234110 | + | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 1273357..1273551 | + | IS66 family transposase |
| FN | CDS | NZ_LN831051.1 | 1286893..1287204 | + | CHY zinc finger protein |
| FN | CDS | NZ_LN831051.1 | 1295430..1295792 | + | G5 domain-containing protein |
| FN | CDS | NZ_LN831051.1 | 1322554..1322793 | + | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 1330351..1330611 | + | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 1416356..1416478 | - | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 1452346..1452471 | - | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 1476860..1477012 | - | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 1516001..1516225 | + | HAD family hydrolase |
| FN | CDS | NZ_LN831051.1 | 1519503..1519949 | - | tyrosine-type recombinase/integrase |
| FN | CDS | NZ_LN831051.1 | 1577534..1577713 | - | DUF3884 family protein |
| FN | CDS | NZ_LN831051.1 | 1619020..1619091 | - | SHP2/SHP3 family peptide pheromone |
| FN | CDS | NZ_LN831051.1 | 1635940..1636383 | - | restriction endonuclease subunit S |
| FN | CDS | NZ_LN831051.1 | 1673056..1673256 | - | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 1693058..1693207 | - | 50S ribosomal protein L33 |
| FN | CDS | NZ_LN831051.1 | 1707418..1707606 | - | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 1719152..1719424 | - | nucleotide-binding domain-containing protein |
| FN | CDS | NZ_LN831051.1 | 1732075..1732272 | + | AraC family transcriptional regulator |
| FN | CDS | NZ_LN831051.1 | 1751301..1751441 | + | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 1757310..1757558 | + | site-specific integrase |
| FN | CDS | NZ_LN831051.1 | 1773388..1773513 | - | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 1869636..1869782 | + | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 1874117..1874326 | + | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 1931896..1932018 | - | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 1981566..1981712 | - | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 1993132..1993314 | + | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 1997131..1997376 | - | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 1997775..1997984 | - | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 1998150..1998296 | + | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 2000396..2000533 | + | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 2001370..2002174 | - | IS5 family transposase |
| FN | CDS | NZ_LN831051.1 | 2028442..2028525 | - | putative autolysis signal peptide Pep27 |
| FN | CDS | NZ_LN831051.1 | 2068225..2068377 | + | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 2080905..2081024 | - | PncF family bacteriocin immunity protein |
| FN | CDS | NZ_LN831051.1 | 2081128..2081277 | - | bacteriocin-like peptide BlpO |
| FN | CDS | NZ_LN831051.1 | 2107831..2107968 | - | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 2109956..2110054 | + | glutamine synthetase |
| FP | CDS | NZ_LN831051.1 | 558..758 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 991..1263 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 54689..54811 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 55009..55557 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 69576..70061 | - | Glycosyltransferase GlyE |
| FP | CDS | NZ_LN831051.1 | 106480..106674 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 106686..107156 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 114648..114977 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 117857..118105 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 118056..118799 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 136563..136874 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 137171..137527 | - | dTDP-4-dehydrorhamnose reductase |
| FP | CDS | NZ_LN831051.1 | 137615..137923 | - | dTDP-4-dehydrorhamnose reductase |
| FP | CDS | NZ_LN831051.1 | 140716..141498 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 141432..141614 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 141672..142022 | + | ISL3 family transposase ISSpn14 |
| FP | CDS | NZ_LN831051.1 | 154662..155060 | + | IS630 family transposase IS630-Spn1 |
| FP | CDS | NZ_LN831051.1 | 155188..155400 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 155499..156281 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 156215..156397 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 156455..156805 | + | ISL3 family transposase ISSpn14 |
| FP | CDS | NZ_LN831051.1 | 161196..161546 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 164179..164268 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 181656..181871 | - | IS200/IS605 family transposase ISStin10 |
| FP | CDS | NZ_LN831051.1 | 189038..189163 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 227465..228421 | - | N-acetyl-S-alkylcysteine monooxygenase |
| FP | CDS | NZ_LN831051.1 | 253615..253716 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 254029..254223 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 254287..254736 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 256088..256456 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 256453..257460 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 281258..281401 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 291487..291720 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 291860..292432 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 316038..316652 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 325072..325209 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 325682..325915 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 327654..328514 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 330982..331911 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 331902..332288 | + | ISL3 family transposase ISSth1 |
| FP | CDS | NZ_LN831051.1 | 335381..335488 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 349014..349829 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 352124..352288 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 353330..353749 | - | Argininosuccinate synthase |
| FP | CDS | NZ_LN831051.1 | 422193..422486 | - | IS5 family transposase ISSpn7 |
| FP | CDS | NZ_LN831051.1 | 430719..430976 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 431300..431533 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 431660..432139 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 433188..433772 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 462098..462451 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 462487..462612 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 463721..463825 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 465018..465176 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 470467..470625 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 476619..476948 | - | IS630 family transposase IS630-Spn1 |
| FP | CDS | NZ_LN831051.1 | 503511..504278 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 504363..504701 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 526186..526455 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 620828..622057 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 639777..639923 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 657230..657610 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 657784..658209 | - | IS1380 family transposase ISSpn5 |
| FP | CDS | NZ_LN831051.1 | 664650..664952 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 665290..665679 | + | IS5 family transposase IS1381 |
| FP | CDS | NZ_LN831051.1 | 705086..705268 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 733741..734091 | + | IS630 family transposase IS630-Spn1 |
| FP | CDS | NZ_LN831051.1 | 759147..760040 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 762918..764600 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 764720..764947 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 765658..765876 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 768468..768854 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 795549..795638 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 799283..799471 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 799482..799586 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 818743..818871 | + | putative transcriptional regulatory protein YeeN |
| FP | CDS | NZ_LN831051.1 | 818881..819459 | + | putative transcriptional regulatory protein YeeN |
| FP | CDS | NZ_LN831051.1 | 819583..819813 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 819823..820032 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 829575..829904 | + | Epoxyqueuosine reductase QueH |
| FP | CDS | NZ_LN831051.1 | 830004..830330 | + | Epoxyqueuosine reductase QueH |
| FP | CDS | NZ_LN831051.1 | 840476..840913 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 859001..859504 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 859619..859894 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 859939..860379 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 860388..860816 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 899345..899572 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 899526..899801 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 925005..925145 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 925182..925736 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 955587..955745 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 959185..959628 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 959638..960204 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 970371..970514 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 972928..976980 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 980105..981178 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 981399..981734 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 991213..991395 | + | Glycosyltransferase GlyG |
| FP | CDS | NZ_LN831051.1 | 991399..992115 | + | Glycosyltransferase GlyG |
| FP | CDS | NZ_LN831051.1 | 1003775..1004734 | - | Serine/threonine transporter SstT |
| FP | CDS | NZ_LN831051.1 | 1024890..1025279 | + | IS5 family transposase IS1381 |

> Primeros pares con coordenadas no exactas

| tipo | seqid | ref | Δstart | Δend | solape |
|---|---|---|---|---|---|
| CDS | NZ_LN831051.1 | 1991..2302 | -33 | +0 | 0.9043 |
| CDS | NZ_LN831051.1 | 9733..10338 | +60 | +0 | 0.9010 |
| CDS | NZ_LN831051.1 | 78527..79000 | +0 | -39 | 0.9177 |
| CDS | NZ_LN831051.1 | 105789..106127 | +0 | +54 | 0.8626 |
| CDS | NZ_LN831051.1 | 192006..193985 | +0 | -6 | 0.9970 |
| CDS | NZ_LN831051.1 | 201815..202789 | +0 | -30 | 0.9692 |
| CDS | NZ_LN831051.1 | 237611..238609 | +0 | -27 | 0.9730 |
| CDS | NZ_LN831051.1 | 307879..309549 | +0 | +12 | 0.9929 |
| CDS | NZ_LN831051.1 | 323503..323931 | +0 | +171 | 0.7150 |
| CDS | NZ_LN831051.1 | 363902..365041 | +0 | -93 | 0.9184 |
| CDS | NZ_LN831051.1 | 392700..393104 | +0 | -66 | 0.8370 |
| CDS | NZ_LN831051.1 | 427737..430370 | +0 | +36 | 0.9865 |
| CDS | NZ_LN831051.1 | 434582..435943 | +0 | -99 | 0.9273 |
| CDS | NZ_LN831051.1 | 435957..436472 | +0 | +42 | 0.9247 |
| CDS | NZ_LN831051.1 | 437853..438002 | -6 | +0 | 0.9615 |
| CDS | NZ_LN831051.1 | 440964..444182 | +0 | -48 | 0.9851 |
| CDS | NZ_LN831051.1 | 459977..460519 | +0 | +48 | 0.9188 |
| CDS | NZ_LN831051.1 | 461321..461479 | +19 | +4 | 0.8805 |
| CDS | NZ_LN831051.1 | 463027..463356 | +0 | -3 | 0.9909 |
| CDS | NZ_LN831051.1 | 465164..465991 | +0 | +30 | 0.9650 |
| CDS | NZ_LN831051.1 | 473252..474538 | +0 | +42 | 0.9684 |
| CDS | NZ_LN831051.1 | 497386..498144 | +0 | -15 | 0.9802 |
| CDS | NZ_LN831051.1 | 549958..551385 | +18 | +0 | 0.9874 |
| CDS | NZ_LN831051.1 | 556182..556790 | +0 | +51 | 0.9227 |
| CDS | NZ_LN831051.1 | 588251..589135 | -15 | +0 | 0.9833 |
| CDS | NZ_LN831051.1 | 605967..606836 | -12 | +0 | 0.9864 |
| CDS | NZ_LN831051.1 | 623548..624012 | +150 | +0 | 0.6774 |
| CDS | NZ_LN831051.1 | 624022..624228 | -9 | +0 | 0.9583 |
| CDS | NZ_LN831051.1 | 629101..629382 | +0 | +75 | 0.7899 |
| CDS | NZ_LN831051.1 | 633828..635081 | +0 | +69 | 0.9478 |
| CDS | NZ_LN831051.1 | 702128..703144 | +144 | +0 | 0.8584 |
| CDS | NZ_LN831051.1 | 708430..709254 | -6 | +0 | 0.9928 |
| CDS | NZ_LN831051.1 | 748284..748490 | +15 | +0 | 0.9275 |
| CDS | NZ_LN831051.1 | 774575..775825 | -894 | +0 | 0.5832 |
| CDS | NZ_LN831051.1 | 782878..783108 | +33 | +0 | 0.8571 |
| CDS | NZ_LN831051.1 | 784408..784596 | -6 | +0 | 0.9692 |
| CDS | NZ_LN831051.1 | 820935..821918 | -48 | +0 | 0.9535 |
| CDS | NZ_LN831051.1 | 893522..894328 | +36 | +0 | 0.9554 |
| CDS | NZ_LN831051.1 | 899954..900250 | +51 | +0 | 0.8283 |
| CDS | NZ_LN831051.1 | 907206..909365 | +33 | +0 | 0.9847 |
| CDS | NZ_LN831051.1 | 956951..959134 | -60 | +0 | 0.9733 |
| CDS | NZ_LN831051.1 | 981898..982218 | -42 | +0 | 0.8843 |
| CDS | NZ_LN831051.1 | 1011074..1011814 | -6 | +0 | 0.9920 |
| CDS | NZ_LN831051.1 | 1025347..1025595 | -39 | +0 | 0.8646 |
| CDS | NZ_LN831051.1 | 1045290..1045625 | -9 | +0 | 0.9739 |
| CDS | NZ_LN831051.1 | 1054009..1054371 | -6 | +0 | 0.9837 |
| CDS | NZ_LN831051.1 | 1074967..1075833 | -21 | +0 | 0.9764 |
| CDS | NZ_LN831051.1 | 1077585..1077863 | -54 | +0 | 0.8378 |
| CDS | NZ_LN831051.1 | 1081851..1083908 | +15 | +0 | 0.9927 |
| CDS | NZ_LN831051.1 | 1135543..1136508 | -18 | +0 | 0.9817 |
| CDS | NZ_LN831051.1 | 1136559..1136765 | +0 | +6 | 0.9718 |
| CDS | NZ_LN831051.1 | 1187212..1187700 | +57 | +0 | 0.8834 |
| CDS | NZ_LN831051.1 | 1243754..1245145 | -51 | +0 | 0.9647 |
| CDS | NZ_LN831051.1 | 1288573..1289379 | +87 | +0 | 0.8922 |
| CDS | NZ_LN831051.1 | 1315850..1316581 | +60 | +0 | 0.9180 |
| CDS | NZ_LN831051.1 | 1341610..1342455 | +9 | +0 | 0.9894 |
| CDS | NZ_LN831051.1 | 1350369..1350740 | -60 | +0 | 0.8611 |
| CDS | NZ_LN831051.1 | 1355496..1356755 | -162 | +0 | 0.8861 |
| CDS | NZ_LN831051.1 | 1383594..1384172 | +33 | +0 | 0.9430 |
| CDS | NZ_LN831051.1 | 1427900..1429216 | -27 | +0 | 0.9799 |
| CDS | NZ_LN831051.1 | 1456991..1457245 | -9 | +0 | 0.9659 |
| CDS | NZ_LN831051.1 | 1472042..1472473 | +39 | +0 | 0.9097 |
| CDS | NZ_LN831051.1 | 1492688..1493287 | +21 | +0 | 0.9650 |
| CDS | NZ_LN831051.1 | 1570351..1572084 | +0 | -9 | 0.9948 |
| CDS | NZ_LN831051.1 | 1618047..1618916 | +0 | -39 | 0.9552 |
| CDS | NZ_LN831051.1 | 1699481..1701646 | +0 | -57 | 0.9737 |
| CDS | NZ_LN831051.1 | 1701688..1702674 | +0 | -48 | 0.9514 |
| CDS | NZ_LN831051.1 | 1704769..1705299 | +0 | +27 | 0.9516 |
| CDS | NZ_LN831051.1 | 1743852..1744526 | +84 | +0 | 0.8756 |
| CDS | NZ_LN831051.1 | 1777895..1779037 | +0 | +60 | 0.9501 |
| CDS | NZ_LN831051.1 | 1779648..1781900 | +0 | +51 | 0.9779 |
| CDS | NZ_LN831051.1 | 1803732..1806203 | +0 | -36 | 0.9854 |
| CDS | NZ_LN831051.1 | 1845314..1845877 | +0 | +9 | 0.9843 |
| CDS | NZ_LN831051.1 | 1891134..1892229 | +0 | -118 | 0.8923 |
| CDS | NZ_LN831051.1 | 1929576..1930274 | +0 | -21 | 0.9700 |
| CDS | NZ_LN831051.1 | 1940539..1940838 | +0 | +48 | 0.8621 |
| CDS | NZ_LN831051.1 | 1945275..1946183 | -60 | +0 | 0.9381 |
| CDS | NZ_LN831051.1 | 1974083..1974529 | -177 | +0 | 0.7163 |
| CDS | NZ_LN831051.1 | 1993814..1994092 | +0 | +3 | 0.9894 |
| CDS | NZ_LN831051.1 | 2032743..2033108 | +0 | +12 | 0.9683 |
| CDS | NZ_LN831051.1 | 2043911..2044177 | +0 | -69 | 0.7416 |
| CDS | NZ_LN831051.1 | 2073903..2075039 | +0 | -12 | 0.9894 |
| tRNA | NZ_LN831051.1 | 480353..480425 | -1 | +1 | 0.9733 |
| tRNA | NZ_LN831051.1 | 681531..681603 | -1 | +1 | 0.9733 |
| tRNA | NZ_LN831051.1 | 684904..684974 | -1 | +1 | 0.9726 |
| tRNA | NZ_LN831051.1 | 685009..685082 | -1 | +1 | 0.9737 |
| tRNA | NZ_LN831051.1 | 685430..685510 | -1 | +1 | 0.9759 |
| tRNA | NZ_LN831051.1 | 685515..685585 | -1 | +1 | 0.9726 |
| tRNA | NZ_LN831051.1 | 685600..685672 | -1 | +1 | 0.9733 |
| tRNA | NZ_LN831051.1 | 685686..685757 | -1 | +1 | 0.9730 |
| tRNA | NZ_LN831051.1 | 685767..685850 | -1 | +1 | 0.9767 |
| tRNA | NZ_LN831051.1 | 743619..743691 | -1 | +1 | 0.9733 |
| tRNA | NZ_LN831051.1 | 747164..747245 | -1 | +1 | 0.9762 |
| tRNA | NZ_LN831051.1 | 747622..747695 | -1 | +1 | 0.9737 |
| tRNA | NZ_LN831051.1 | 787279..787349 | -2 | +1 | 0.9595 |
| tRNA | NZ_LN831051.1 | 835544..835616 | -1 | +1 | 0.9733 |
| tRNA | NZ_LN831051.1 | 839089..839170 | -1 | +1 | 0.9762 |
| tRNA | NZ_LN831051.1 | 839547..839620 | -1 | +1 | 0.9737 |
| tRNA | NZ_LN831051.1 | 839626..839699 | -1 | +1 | 0.9737 |
| tRNA | NZ_LN831051.1 | 840175..840248 | -1 | +1 | 0.9737 |

- Aviso: 2 feature(s) GFF multi-segmento fusionados por ID; se comparan por su intervalo envolvente.

## DFAST

- Estado: **ok**.
- Fuente: `executed`.
- Versión: `DFAST ver. 1.6.2` (estado `ok`).
- Ejecución: **ok**, 283.60 s, RSS máx. 2989.6 MiB (`2026-11-26T02:39:08Z` → `2026-11-26T02:43:51Z`).
- Salida: `/tmp/external_benchmark_dfast_azkc9rt9/out/genome.gff` (3795676 bytes, SHA256 `9eabdd92070599cb…`).
- Comando: `/usr/local/bin/dfast --genome /tmp/benchmark_genome_otq2w3ys.fna --out /tmp/external_benchmark_dfast_azkc9rt9/out --cpu 4 --force --dbroot /data/dfast_db`.
- Anotación comparada: `reports/external_benchmark/heldout/spneumoniae/evidence/dfast/genome.gff` (formato `gff3`, SHA256 `9eabdd92070599cb…`).

### Features (DFAST)

| tipo | ref | query | TP | FP | FN | precisión | recall | F1 | exactos | borde distinto | solape medio |
|---|---|---|---|---|---|---|---|---|---|---|---|
| CDS | 1931 | 2091 | 1846 | 245 | 85 | 0.8828 | 0.9560 | 0.9180 | 1764 | 82 | 0.9967 |
| tRNA | 58 | 58 | 58 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 35 | 23 | 0.9892 |
| rRNA | 12 | 12 | 12 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 0 | 12 | 0.9762 |

### Atributos (DFAST)

| tipo | pares | product en ambos | product literal | product tokens | gene literal | locus literal |
|---|---|---|---|---|---|---|
| CDS | 1846 | 1846 | 0.8115 | 0.8638 | 0.9862 | 0.0000 |
| tRNA | 58 | 58 | 0.9828 | 0.9885 | — | 0.0000 |
| rRNA | 12 | 12 | 1.0000 | 1.0000 | — | 0.0000 |

### Hipotéticas (DFAST)

| tipo | ref CDS | ref con product | ref hyp | % ref | query CDS | query con product | query hyp | % query |
|---|---|---|---|---|---|---|---|---|
| CDS | 1931 | 1931 | 268 | 13.9 | 2091 | 2091 | 502 | 24.0 |
| tRNA | 0 | 0 | 0 | — | 0 | 0 | 0 | — |
| rRNA | 0 | 0 | 0 | — | 0 | 0 | 0 | — |

#### Diferencias

- Falsos negativos (referencia sin pareja): **85**.
- Falsos positivos (herramienta sin pareja): **245**.
- Pares solo con hebra distinta: **1**.
- Pares emparejados con algún borde distinto: **117**.

Las listas de diferencias están truncadas a 100 elementos por lista; los conteos son completos.

> Primeras diferencias no emparejadas

| tipo | feature | seqid | coordenadas | hebra | product/locus |
|---|---|---|---|---|---|
| FN | CDS | NZ_LN831051.1 | 1588..2076 | + | phosphoglycerate kinase |
| FN | CDS | NZ_LN831051.1 | 5195..5314 | - | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 7280..7417 | + | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 11755..11892 | - | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 13858..13977 | + | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 26991..27116 | - | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 52213..52335 | + | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 76489..76674 | - | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 78278..78514 | - | CPBP family intramembrane glutamic endopeptidase |
| FN | CDS | NZ_LN831051.1 | 115640..115978 | + | lanthionine synthetase LanC family protein |
| FN | CDS | NZ_LN831051.1 | 238743..238835 | + | type I toxin-antitoxin system Fst family toxin |
| FN | CDS | NZ_LN831051.1 | 253201..253398 | + | helix-turn-helix domain-containing protein |
| FN | CDS | NZ_LN831051.1 | 277204..277344 | - | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 316038..316190 | - | AzlD domain-containing protein |
| FN | CDS | NZ_LN831051.1 | 325239..325631 | - | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 328845..329030 | - | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 337699..337857 | - | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 344487..344726 | - | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 347913..348095 | - | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 351359..351595 | - | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 361731..361898 | + | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 386482..386634 | + | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 405670..405786 | - | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 416534..416644 | - | phosphoribosylaminoimidazolesuccinocarboxamide synthase |
| FN | CDS | NZ_LN831051.1 | 421632..421757 | - | PncF family bacteriocin immunity protein |
| FN | CDS | NZ_LN831051.1 | 463328..463537 | - | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 556922..557119 | + | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 561107..561238 | + | teichoic acid D-Ala incorporation-associated protein DltX |
| FN | CDS | NZ_LN831051.1 | 641113..641319 | + | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 657535..657669 | + | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 675941..676120 | + | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 707936..708067 | + | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 740073..740225 | + | 50S ribosomal protein L33 |
| FN | CDS | NZ_LN831051.1 | 810117..810278 | + | helix-turn-helix domain-containing protein |
| FN | CDS | NZ_LN831051.1 | 824247..824438 | + | cell wall-binding repeat protein |
| FN | CDS | NZ_LN831051.1 | 824463..824525 | + | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 876879..876989 | - | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 917884..918072 | + | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 953070..953249 | + | type II toxin-antitoxin system HicA family toxin |
| FN | CDS | NZ_LN831051.1 | 970783..970926 | + | IS630 transposase-related protein |
| FN | CDS | NZ_LN831051.1 | 971597..976975 | + | accessory Sec-dependent serine-rich glycoprotein adhesin |
| FN | CDS | NZ_LN831051.1 | 1077873..1078064 | + | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 1098064..1098228 | + | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 1124411..1124587 | - | MazG-like protein |
| FN | CDS | NZ_LN831051.1 | 1134273..1134419 | + | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 1155403..1155549 | - | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 1232933..1233226 | + | IS630 transposase-related protein |
| FN | CDS | NZ_LN831051.1 | 1233868..1234110 | + | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 1273357..1273551 | + | IS66 family transposase |
| FN | CDS | NZ_LN831051.1 | 1286893..1287204 | + | CHY zinc finger protein |
| FN | CDS | NZ_LN831051.1 | 1295430..1295792 | + | G5 domain-containing protein |
| FN | CDS | NZ_LN831051.1 | 1322554..1322793 | + | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 1330351..1330611 | + | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 1416356..1416478 | - | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 1452346..1452471 | - | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 1476860..1477012 | - | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 1516001..1516225 | + | HAD family hydrolase |
| FN | CDS | NZ_LN831051.1 | 1519503..1519949 | - | tyrosine-type recombinase/integrase |
| FN | CDS | NZ_LN831051.1 | 1577534..1577713 | - | DUF3884 family protein |
| FN | CDS | NZ_LN831051.1 | 1619020..1619091 | - | SHP2/SHP3 family peptide pheromone |
| FN | CDS | NZ_LN831051.1 | 1635940..1636383 | - | restriction endonuclease subunit S |
| FN | CDS | NZ_LN831051.1 | 1673056..1673256 | - | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 1693058..1693207 | - | 50S ribosomal protein L33 |
| FN | CDS | NZ_LN831051.1 | 1707418..1707606 | - | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 1719152..1719424 | - | nucleotide-binding domain-containing protein |
| FN | CDS | NZ_LN831051.1 | 1732075..1732272 | + | AraC family transcriptional regulator |
| FN | CDS | NZ_LN831051.1 | 1751301..1751441 | + | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 1757310..1757558 | + | site-specific integrase |
| FN | CDS | NZ_LN831051.1 | 1773388..1773513 | - | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 1869636..1869782 | + | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 1874117..1874326 | + | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 1931896..1932018 | - | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 1981566..1981712 | - | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 1993132..1993314 | + | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 1997131..1997376 | - | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 1997775..1997984 | - | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 1998150..1998296 | + | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 2000396..2000533 | + | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 2001370..2002174 | - | IS5 family transposase |
| FN | CDS | NZ_LN831051.1 | 2028442..2028525 | - | putative autolysis signal peptide Pep27 |
| FN | CDS | NZ_LN831051.1 | 2068225..2068377 | + | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 2080905..2081024 | - | PncF family bacteriocin immunity protein |
| FN | CDS | NZ_LN831051.1 | 2081128..2081277 | - | bacteriocin-like peptide BlpO |
| FN | CDS | NZ_LN831051.1 | 2107831..2107968 | - | hypothetical protein |
| FN | CDS | NZ_LN831051.1 | 2109956..2110054 | + | glutamine synthetase |
| FP | CDS | NZ_LN831051.1 | 558..758 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 991..1263 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 54689..54811 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 55009..55557 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 69576..70061 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 106480..106674 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 106686..107156 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 114648..114977 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 117857..118105 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 118056..118799 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 136563..136874 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 137171..137527 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 137615..137923 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 140716..141498 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 141432..141614 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 141672..142022 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 154662..155060 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 155188..155400 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 155499..156281 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 156215..156397 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 156455..156805 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 161196..161546 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 164179..164268 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 181656..181871 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 189038..189163 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 227465..228421 | - | LLM class flavin-dependent oxidoreductase |
| FP | CDS | NZ_LN831051.1 | 253615..253716 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 254029..254223 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 254287..254736 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 256088..256456 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 256453..257460 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 281258..281401 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 291487..291720 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 291860..292432 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 316038..316652 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 325072..325209 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 325682..325915 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 327654..328514 | - | glycosyltransferase family 2 protein |
| FP | CDS | NZ_LN831051.1 | 330982..331911 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 331902..332288 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 335381..335488 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 349014..349829 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 352124..352288 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 353330..353749 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 422193..422486 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 430719..430976 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 431300..431533 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 431660..432139 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 433188..433772 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 462098..462451 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 462487..462612 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 463721..463825 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 465018..465176 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 470467..470625 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 476619..476948 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 503511..504278 | + | YhgE/Pip family protein |
| FP | CDS | NZ_LN831051.1 | 504363..504701 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 526186..526455 | + | transposase |
| FP | CDS | NZ_LN831051.1 | 620828..622057 | + | transporter |
| FP | CDS | NZ_LN831051.1 | 639777..639923 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 657230..657610 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 657784..658209 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 664650..664952 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 665290..665679 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 705086..705268 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 733741..734091 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 759147..760040 | + | DUF1542 domain-containing protein |
| FP | CDS | NZ_LN831051.1 | 762918..764600 | + | bacteriocin-associated integral membrane family protein |
| FP | CDS | NZ_LN831051.1 | 764720..764947 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 765658..765876 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 768468..768854 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 795549..795638 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 799283..799471 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 799482..799586 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 818743..818871 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 818881..819459 | + | YebC/PmpR family DNA-binding transcriptional regulator |
| FP | CDS | NZ_LN831051.1 | 819583..819813 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 819823..820032 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 829575..829904 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 830004..830330 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 840476..840913 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 859001..859504 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 859619..859894 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 859939..860379 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 860388..860816 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 899345..899572 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 899526..899801 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 925005..925145 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 925182..925736 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 955587..955745 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 959185..959628 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 959638..960204 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 970371..970514 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 972928..976980 | - | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 980105..981178 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 981399..981734 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 991213..991395 | + | hypothetical protein |
| FP | CDS | NZ_LN831051.1 | 991399..992115 | + | glycosyltransferase family 2 protein |
| FP | CDS | NZ_LN831051.1 | 1003775..1004734 | - | serine/threonine transporter SstT |
| FP | CDS | NZ_LN831051.1 | 1024890..1025279 | + | hypothetical protein |

> Primeros pares con coordenadas no exactas

| tipo | seqid | ref | Δstart | Δend | solape |
|---|---|---|---|---|---|
| CDS | NZ_LN831051.1 | 1991..2302 | -33 | +0 | 0.9043 |
| CDS | NZ_LN831051.1 | 9733..10338 | +60 | +0 | 0.9010 |
| CDS | NZ_LN831051.1 | 78527..79000 | +0 | -39 | 0.9177 |
| CDS | NZ_LN831051.1 | 105789..106127 | +0 | +54 | 0.8626 |
| CDS | NZ_LN831051.1 | 192006..193985 | +0 | -6 | 0.9970 |
| CDS | NZ_LN831051.1 | 201815..202789 | +0 | -30 | 0.9692 |
| CDS | NZ_LN831051.1 | 237611..238609 | +0 | -27 | 0.9730 |
| CDS | NZ_LN831051.1 | 307879..309549 | +0 | +12 | 0.9929 |
| CDS | NZ_LN831051.1 | 323503..323931 | +0 | +171 | 0.7150 |
| CDS | NZ_LN831051.1 | 363902..365041 | +0 | -93 | 0.9184 |
| CDS | NZ_LN831051.1 | 392700..393104 | +0 | -66 | 0.8370 |
| CDS | NZ_LN831051.1 | 427737..430370 | +0 | +36 | 0.9865 |
| CDS | NZ_LN831051.1 | 434582..435943 | +0 | -99 | 0.9273 |
| CDS | NZ_LN831051.1 | 435957..436472 | +0 | +42 | 0.9247 |
| CDS | NZ_LN831051.1 | 437853..438002 | -6 | +0 | 0.9615 |
| CDS | NZ_LN831051.1 | 440964..444182 | +0 | -48 | 0.9851 |
| CDS | NZ_LN831051.1 | 459977..460519 | +0 | +48 | 0.9188 |
| CDS | NZ_LN831051.1 | 461321..461479 | +19 | +4 | 0.8805 |
| CDS | NZ_LN831051.1 | 463027..463356 | +0 | -3 | 0.9909 |
| CDS | NZ_LN831051.1 | 465164..465991 | +0 | +30 | 0.9650 |
| CDS | NZ_LN831051.1 | 473252..474538 | +0 | +42 | 0.9684 |
| CDS | NZ_LN831051.1 | 497386..498144 | +0 | -15 | 0.9802 |
| CDS | NZ_LN831051.1 | 549958..551385 | +18 | +0 | 0.9874 |
| CDS | NZ_LN831051.1 | 556182..556790 | +0 | +51 | 0.9227 |
| CDS | NZ_LN831051.1 | 588251..589135 | -15 | +0 | 0.9833 |
| CDS | NZ_LN831051.1 | 605967..606836 | -12 | +0 | 0.9864 |
| CDS | NZ_LN831051.1 | 623548..624012 | +150 | +0 | 0.6774 |
| CDS | NZ_LN831051.1 | 624022..624228 | -9 | +0 | 0.9583 |
| CDS | NZ_LN831051.1 | 629101..629382 | +0 | +75 | 0.7899 |
| CDS | NZ_LN831051.1 | 633828..635081 | +0 | +69 | 0.9478 |
| CDS | NZ_LN831051.1 | 702128..703144 | +144 | +0 | 0.8584 |
| CDS | NZ_LN831051.1 | 708430..709254 | -6 | +0 | 0.9928 |
| CDS | NZ_LN831051.1 | 748284..748490 | +15 | +0 | 0.9275 |
| CDS | NZ_LN831051.1 | 774575..775825 | -894 | +0 | 0.5832 |
| CDS | NZ_LN831051.1 | 782878..783108 | +33 | +0 | 0.8571 |
| CDS | NZ_LN831051.1 | 784408..784596 | -6 | +0 | 0.9692 |
| CDS | NZ_LN831051.1 | 820935..821918 | -48 | +0 | 0.9535 |
| CDS | NZ_LN831051.1 | 893522..894328 | +36 | +0 | 0.9554 |
| CDS | NZ_LN831051.1 | 899954..900250 | +51 | +0 | 0.8283 |
| CDS | NZ_LN831051.1 | 907206..909365 | +33 | +0 | 0.9847 |
| CDS | NZ_LN831051.1 | 956951..959134 | -60 | +0 | 0.9733 |
| CDS | NZ_LN831051.1 | 981898..982218 | -42 | +0 | 0.8843 |
| CDS | NZ_LN831051.1 | 1011074..1011814 | -6 | +0 | 0.9920 |
| CDS | NZ_LN831051.1 | 1025347..1025595 | -39 | +0 | 0.8646 |
| CDS | NZ_LN831051.1 | 1045290..1045625 | -9 | +0 | 0.9739 |
| CDS | NZ_LN831051.1 | 1054009..1054371 | -6 | +0 | 0.9837 |
| CDS | NZ_LN831051.1 | 1074967..1075833 | -21 | +0 | 0.9764 |
| CDS | NZ_LN831051.1 | 1077585..1077863 | -54 | +0 | 0.8378 |
| CDS | NZ_LN831051.1 | 1081851..1083908 | +15 | +0 | 0.9927 |
| CDS | NZ_LN831051.1 | 1135543..1136508 | -18 | +0 | 0.9817 |
| CDS | NZ_LN831051.1 | 1136559..1136765 | +0 | +6 | 0.9718 |
| CDS | NZ_LN831051.1 | 1187212..1187700 | +57 | +0 | 0.8834 |
| CDS | NZ_LN831051.1 | 1243754..1245145 | -51 | +0 | 0.9647 |
| CDS | NZ_LN831051.1 | 1288573..1289379 | +87 | +0 | 0.8922 |
| CDS | NZ_LN831051.1 | 1315850..1316581 | +60 | +0 | 0.9180 |
| CDS | NZ_LN831051.1 | 1341610..1342455 | +9 | +0 | 0.9894 |
| CDS | NZ_LN831051.1 | 1350369..1350740 | -60 | +0 | 0.8611 |
| CDS | NZ_LN831051.1 | 1355496..1356755 | -162 | +0 | 0.8861 |
| CDS | NZ_LN831051.1 | 1383594..1384172 | +33 | +0 | 0.9430 |
| CDS | NZ_LN831051.1 | 1427900..1429216 | -27 | +0 | 0.9799 |
| CDS | NZ_LN831051.1 | 1456991..1457245 | -9 | +0 | 0.9659 |
| CDS | NZ_LN831051.1 | 1472042..1472473 | +39 | +0 | 0.9097 |
| CDS | NZ_LN831051.1 | 1492688..1493287 | +21 | +0 | 0.9650 |
| CDS | NZ_LN831051.1 | 1570351..1572084 | +0 | -9 | 0.9948 |
| CDS | NZ_LN831051.1 | 1618047..1618916 | +0 | -39 | 0.9552 |
| CDS | NZ_LN831051.1 | 1699481..1701646 | +0 | -57 | 0.9737 |
| CDS | NZ_LN831051.1 | 1701688..1702674 | +0 | -48 | 0.9514 |
| CDS | NZ_LN831051.1 | 1704769..1705299 | +0 | +27 | 0.9516 |
| CDS | NZ_LN831051.1 | 1743852..1744526 | +84 | +0 | 0.8756 |
| CDS | NZ_LN831051.1 | 1777895..1779037 | +0 | +60 | 0.9501 |
| CDS | NZ_LN831051.1 | 1779648..1781900 | +0 | +51 | 0.9779 |
| CDS | NZ_LN831051.1 | 1803732..1806203 | +0 | -36 | 0.9854 |
| CDS | NZ_LN831051.1 | 1845314..1845877 | +0 | +9 | 0.9843 |
| CDS | NZ_LN831051.1 | 1891134..1892229 | +0 | -118 | 0.8923 |
| CDS | NZ_LN831051.1 | 1929576..1930274 | +0 | -21 | 0.9700 |
| CDS | NZ_LN831051.1 | 1940539..1940838 | +0 | +48 | 0.8621 |
| CDS | NZ_LN831051.1 | 1945275..1946183 | -60 | +0 | 0.9381 |
| CDS | NZ_LN831051.1 | 1974083..1974529 | -177 | +0 | 0.7163 |
| CDS | NZ_LN831051.1 | 1993814..1994092 | +0 | +3 | 0.9894 |
| CDS | NZ_LN831051.1 | 2032743..2033108 | +0 | +12 | 0.9683 |
| CDS | NZ_LN831051.1 | 2043911..2044177 | +0 | -69 | 0.7416 |
| CDS | NZ_LN831051.1 | 2073903..2075039 | +0 | -12 | 0.9894 |
| rRNA | NZ_LN831051.1 | 477135..477250 | +4 | -4 | 0.9310 |
| rRNA | NZ_LN831051.1 | 477327..480228 | +2 | -1 | 0.9990 |
| rRNA | NZ_LN831051.1 | 480477..482022 | +0 | -2 | 0.9987 |
| rRNA | NZ_LN831051.1 | 679934..681479 | +2 | +0 | 0.9987 |
| rRNA | NZ_LN831051.1 | 681728..684629 | +1 | -2 | 0.9990 |
| rRNA | NZ_LN831051.1 | 684706..684821 | +4 | -4 | 0.9310 |
| rRNA | NZ_LN831051.1 | 742022..743567 | +2 | +0 | 0.9987 |
| rRNA | NZ_LN831051.1 | 743816..746717 | +1 | -2 | 0.9990 |
| rRNA | NZ_LN831051.1 | 746794..746909 | +4 | -4 | 0.9310 |
| rRNA | NZ_LN831051.1 | 833947..835492 | +2 | +0 | 0.9987 |
| rRNA | NZ_LN831051.1 | 835741..838642 | +1 | -2 | 0.9990 |
| rRNA | NZ_LN831051.1 | 838719..838834 | +4 | -4 | 0.9310 |
| tRNA | NZ_LN831051.1 | 480353..480425 | -1 | +1 | 0.9733 |
| tRNA | NZ_LN831051.1 | 681531..681603 | -1 | +1 | 0.9733 |
| tRNA | NZ_LN831051.1 | 684904..684974 | -1 | +1 | 0.9726 |
| tRNA | NZ_LN831051.1 | 685009..685082 | -1 | +1 | 0.9737 |
| tRNA | NZ_LN831051.1 | 685430..685510 | -1 | +1 | 0.9759 |
| tRNA | NZ_LN831051.1 | 685515..685585 | -1 | +1 | 0.9726 |

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
/opt/conda/bin/python benchmarks/external/external_benchmark.py --genome data/reference_genomes/spneumoniae/spneumoniae_genomic.fna.gz --reference data/reference_genomes/spneumoniae/spneumoniae_genomic.gff.gz --tools bakta,prokka,dfast --run --threads 1 --json reports/external_benchmark/benchmark.json --md reports/external_benchmark/benchmark.md --bakta-db /data/bakta_db
```

Para comparar salidas producidas en otro host, usar `--query-gff TOOL=RUTA` sin `--run`; el origen queda como `precomputed`.
