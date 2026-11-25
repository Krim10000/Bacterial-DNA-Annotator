# Inventario de catálogos ComplexPortal

Generado: `2026-11-25T00:16:28Z`

Alcance: genomas con `data/regulation/*/genes.tsv`; especie/taxid desde `data/taxonomy/full/*.json` (fallback `data/taxonomy/*.json` y GBFF) y catálogos locales en `data/complexes/`. Este informe no descarga nada ni consulta la red.

## Resumen

| métrica | valor |
|---|---:|
| genomas | 401 |
| panel GCF | 400 |
| referencias no-GCF | 1 |
| con especie | 401 |
| con taxid | 401 |
| panel con taxid | 400 |
| con catálogo local | 1 |
| panel con catálogo local | 0 |
| catálogo gold listo | 1 |
| TSV crudo sin mapear | 0 |
| sin catálogo local | 400 |
| sin taxid | 0 |
| sin taxonomía | 0 |

`not_available` significa que el dato no está en las fuentes locales; no significa que el taxid o el catálogo no existan en ComplexPortal.

## Estados

| status | n | significado |
|---|---:|---|
| `ready` | 1 | catálogo gold local + taxid |
| `raw_only` | 0 | TSV crudo local, sin mapear a genes |
| `catalog_missing` | 400 | taxid disponible, sin catálogo local |
| `taxid_missing` | 0 | sin taxid; no se puede resolver catálogo |
| `taxonomy_missing` | 0 | sin especie, taxid ni JSON de taxonomía |

## Catálogos locales

| genome | species | taxid | local_catalog | status |
|---|---|---|---|---|
| `ecoli_k12` | Escherichia coli str. K-12 substr. MG1655 | 83333 | `data/complexes/gold_k12.tsv` | `ready` |

## Inventario completo

| genome | species | taxid | local_catalog | status | reason |
|---|---|---|---|---|---|
| `GCF_000006685.1_ASM668v1` | Chlamydia muridarum | 243161 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 243161; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000006985.1_ASM698v1` | Chlorobaculum tepidum | 194439 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 194439; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000007205.1_ASM720v1` | Chlamydia pneumoniae | 182082 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 182082; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000007605.1_ASM760v1` | Chlamydia caviae | 227941 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 227941; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000007985.2_ASM798v2` | Geobacter sulfurreducens | 243231 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 243231; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000008625.1_ASM862v1` | Aquifex aeolicus | 224324 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 224324; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000008725.1_ASM872v1` | Chlamydia trachomatis | 272561 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 272561; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000010305.1_ASM1030v1` | Gemmatimonas aurantiaca | 379066 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 379066; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000010985.1_ASM1098v1` | Deferribacter desulfuricans | 639282 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 639282; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000011905.1_ASM1190v1` | Dehalococcoides mccartyi | 243164 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 243164; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000012485.1_ASM1248v1` | Pelodictyon luteolum | 319225 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 319225; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000014025.1_ASM1402v1` | Aliarcobacter butzleri | 367737 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 367737; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000014125.1_ASM1412v1` | Clostridium novyi | 386415 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 386415; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000015045.1_ASM1504v1` | Pelobacter propionicus | 338966 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 338966; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000015125.1_ASM1512v1` | Chlorobium phaeobacteroides | 290317 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 290317; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000017685.1_ASM1768v1` | Leptospira biflexa | 456481 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 456481; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000017805.1_ASM1780v1` | Roseiflexus castenholzii | 383372 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 383372; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000017865.1_ASM1786v1` | Pseudothermotoga lettingae | 416591 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 416591; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000018605.1_ASM1860v1` | Petrotoga mobilis | 403833 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 403833; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000018865.1_ASM1886v1` | Chloroflexus aurantiacus | 324602 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 324602; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000018945.1_ASM1894v1` | Thermotoga neapolitana | 309803 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 309803; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000019965.1_ASM1996v1` | Opitutus terrae | 452637 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 452637; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000020145.1_ASM2014v1` | Elusimicrobium minutum | 445932 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 445932; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000020465.1_ASM2046v1` | Chlorobium limicola | 290315 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 290315; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000020505.1_ASM2050v1` | Chlorobaculum parvum | 517417 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 517417; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000020525.1_ASM2052v1` | Chloroherpeton thalassium | 517418 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 517418; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000020625.1_ASM2062v1` | Prosthecochloris aestuarii | 290512 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 290512; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000020645.1_ASM2064v1` | Pelodictyon phaeoclathratiforme | 324925 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 324925; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000020945.1_ASM2094v1` | Coprothermobacter proteolyticus | 309798 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 309798; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000020965.1_ASM2096v1` | Dictyoglomus thermophilum | 309799 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 309799; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000020985.1_ASM2098v1` | Thermodesulfovibrio yellowstonii | 289376 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 289376; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000021545.1_ASM2154v1` | Sulfurihydrogenibium azorense | 204536 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 204536; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000021565.1_ASM2156v1` | Persephonella marina | 123214 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 123214; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000021645.1_ASM2164v1` | Dictyoglomus turgidum | 515635 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 515635; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000021685.1_ASM2168v1` | Thermomicrobium roseum | 309801 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 309801; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000021905.1_ASM2190v1` | Desulfatibacillum aliphaticivorans | 218208 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 218208; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000021945.1_ASM2194v1` | Chloroflexus aggregans | 326427 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 326427; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000022145.1_ASM2214v1` | Anaeromyxobacter dehalogenans | 455488 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 455488; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000022565.1_ASM2256v1` | Acidobacterium capsulatum | 240015 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 240015; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000023325.1_ASM2332v1` | Kosmotoga olearia | 521045 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 521045; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000024405.1_ASM2440v1` | Sebaldella termitidis | 526218 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 526218; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000024425.1_ASM2442v1` | Meiothermus ruber | 504728 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 504728; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000024845.1_ASM2484v1` | Rhodothermus marinus | 518766 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 518766; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000024905.1_ASM2490v1` | Thermanaerovibrio acidaminovorans | 525903 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 525903; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000024985.1_ASM2498v1` | Sphaerobacter thermophilus | 479434 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 479434; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000025005.1_ASM2500v1` | Thermobaculum terrenum | 525904 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 525904; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000025605.1_ASM2560v1` | Thermocrinis albus | 638303 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 638303; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000025725.1_ASM2572v1` | Denitrovibrio acetiphilus | 522772 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 522772; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000025885.1_ASM2588v1` | Aminobacterium colombiense | 572547 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 572547; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000025905.1_ASM2590v1` | Coraliomargarita akajimensis | 583355 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 583355; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000025945.1_ASM2594v1` | Desulfotalea psychrophila | 177439 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 177439; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000091545.1_ASM9154v1` | Thermus thermophilus | 300852 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 300852; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000092425.1_ASM9242v1` | Truepera radiovictrix | 649638 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 649638; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000092785.1_ASM9278v1` | Waddlia chondrophila | 716544 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 716544; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000143725.1_ASM14372v1` | Brachyspira pilosicoli | 759914 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 759914; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000146065.2_ASM14606v1` | Cloacimonas acidaminivorans (SeqCode) | 459349 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 459349; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000146505.1_ASM14650v1` | Fibrobacter succinogenes | 59374 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 59374; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000147715.2_ASM14771v3` | Mesotoga prima | 660470 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 660470; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000164905.1_ASM16490v1` | Hydrogenobacter thermophilus | 608538 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 608538; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000172155.1_ASM17215v1` | Verrucomicrobium spinosum | 240016 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 240016; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000177635.2_ASM17763v2` | Desulfurispirillum indicum | 653733 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 653733; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000179915.2_ASM17991v2` | Terriglobus saanensis | 401053 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 401053; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000183405.1_ASM18340v1` | Calditerrivibrio nitroreducens | 768670 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 768670; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000185805.1_ASM18580v1` | Thermovibrio ammonificans | 648996 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 648996; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000186245.1_ASM18624v1` | Nitratifractor salsuginis | 749222 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 749222; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000191045.1_ASM19104v1` | Desulfurobacterium thermolithotrophum | 868864 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 868864; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000195335.1_ASM19533v1` | Marinithermus hydrothermalis | 869210 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 869210; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000195975.1_ASM19597v1` | Parasynechococcus marenigrum | 84588 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 84588; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000196115.1_ASM19611v1` | Rhodopirellula baltica | 243090 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 243090; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000196175.1_ASM19617v1` | Bdellovibrio bacteriovorus | 264462 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 264462; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000199675.1_ASM19967v1` | Anaerolinea thermophila | 926569 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 926569; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000204255.1_ASM20425v1` | Chlamydia psittaci | 331636 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 331636; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000210915.2_ASM21091v2` | Halobacteriovorax marinus | 862908 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 862908; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000212395.1_ASM21239v1` | Thermodesulfobium narugense | 747365 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 747365; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000212675.2_ASM21267v2` | Desulfocurvibacter africanus | 690850 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 690850; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000214155.1_ASM21415v1` | Mycolicibacter sinensis | 875328 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 875328; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000215975.1_ASM21597v1` | Thermodesulfobacterium geofontis | 795359 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 795359; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000218625.1_ASM21862v1` | Flexistipes sinusarabici | 717231 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 717231; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000219725.1_ASM21972v1` | Gracilinema caldarium | 744872 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 744872; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000225345.1_ASM22534v1` | Roseburia hominis | 585394 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 585394; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000226295.1_ASM22629v1` | Chloracidobacterium thermophilum | 981222 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 981222; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000226315.1_ASM22631v1` | Micavibrio aeruginosavorus | 856793 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 856793; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000235405.2_ASM23540v3` | Fervidobacterium pennivorans | 771875 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 771875; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000237205.1_ASM23720v1` | Simkania negevensis | 331113 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 331113; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000239475.1_ASM23947v1` | Borreliella garinii | 1081646 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 1081646; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000252445.1_ASM25244v1` | Deinococcus gobiensis | 745776 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 745776; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000253035.1_ASM25303v1` | Parachlamydia acanthamoebae | 765952 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 765952; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000255135.1_ASM25513v1` | Marinitoga piezophila | 443254 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 443254; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000258405.1_ASM25840v1` | Ignavibacterium album | 945713 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 945713; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000259275.1_ASM25927v1` | Helicobacter cetorum | 1163745 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 1163745; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000265295.1_ASM26529v1` | Sulfurospirillum barnesii | 760154 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 760154; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000265425.1_ASM26542v1` | Terriglobus roseus | 926566 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 926566; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000279145.1_ASM27914v1` | Melioribacter roseus | 1191523 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 1191523; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000284115.1_ASM28411v1` | Phycisphaera mikurensis | 1142394 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 1142394; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000284315.1_ASM28431v1` | Leptospirillum ferrooxidans | 1162668 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 1162668; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000284335.1_ASM28433v1` | Caldisericum exile | 511051 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 511051; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000309885.1_ASM30988v1` | Thermus oshimai | 751945 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 751945; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000317105.1_ASM31710v1` | Oxynema acuminatum | 56110 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 56110; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000317575.1_ASM31757v1` | Stanieria cyanosphaera | 111780 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 111780; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000348725.1_ASM34872v1` | Pseudobdellovibrio exovorus | 1184267 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 1184267; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000392435.1_ASM39243v1` | Saccharimonas aalborgensis (SeqCode) | 1332188 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 1332188; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000439435.1_ASM43943v1` | Spiroplasma taiwanense | 1276220 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 1276220; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000455605.1_GPTSA100-9` | Flavobacterium indicum | 1094466 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 1094466; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000471025.2_ASM47102v2` | Chlamydia gallinacea | 1143323 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 1143323; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000479335.1_ASM47933v1` | Streptococcus ilei | 1156431 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 1156431; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000484535.1_ASM48453v1` | Gloeobacter kilaueensis | 1183438 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 1183438; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000487995.2_ASM48799v2` | Mucispirillum schaedleri | 1379858 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 1379858; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000504085.1_ASM50408v1` | Pseudothermotoga elfii | 1123383 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 1123383; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000508245.1_ASM50824v1` | Mycoplasma ovis | 1415773 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 1415773; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000512735.1_ASM51273v1` | Thermocrinis ruber | 75906 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 75906; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000513475.1_Bbla` | Candidatus Babela massiliensis | 673862 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 673862; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000522985.1_ASM52298v1` | Gemmatirosa kalamazoonensis | 861299 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 861299; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000582535.1_ASM58253v1` | Hepatoplasma crinochetorum (SeqCode) | 1427984 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 1427984; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000695095.2_ASM69509v2` | Gemmatimonas phototrophica | 1379270 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 1379270; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000724625.1_ASM72462v1` | Fimbriimonas ginsengisoli | 661478 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 661478; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000767275.3_ASM76727v4` | Fervidobacterium islandicum | 2423 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 2423; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000800395.1_ASM80039v1` | Deinococcus radiopugnans | 57497 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 57497; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000803625.1_ASM80362v1` | Nanosynbacter lyticus | 2093824 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 2093824; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000828675.1_ASM82867v1` | Pseudothermotoga profunda | 1408160 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 1408160; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000829235.1_ASM82923v1` | cyanobacterium endosymbiont of Epithemia turgida | 1228987 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 1228987; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_000941075.1_ASM94107v1` | Mesomycoplasma dispar | 86660 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 86660; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_001017655.1_ASM101765v1` | Kiritimatiella glycovorans | 1307763 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 1307763; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_001027025.1_ASM102702v1` | Kosmotoga pacifica | 1330330 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 1330330; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_001027285.1_ASM102728v1` | Archangium gephyra | 48 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 48; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_001027545.1_ASM102754v1` | Endomicrobium proavitum | 1408281 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 1408281; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_001273775.1_ASM127377v1` | Nitrospira moscoviensis | 42253 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 42253; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_001296145.1_ASM129614v1` | Fusobacterium animalis | 76859 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 76859; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_001458695.1_NiCh1` | Candidatus Nitrospira inopinata | 1715989 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 1715989; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_001547735.1_ASM154773v1` | Thermosulfidibacter takaii | 1298851 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 1298851; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_001610955.1_ASM161095v1` | Geobacillus thermoleovorans | 33941 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 33941; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_001611795.1_ASM161179v1` | Streptomyces qaidamensis | 1783515 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 1783515; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_001618865.1_ASM161886v1` | Luteitalea pratensis | 1855912 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 1855912; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_001676785.2_ASM167678v2` | Brachyspira hyodysenteriae | 1266923 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 1266923; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_001701045.1_ASM170104v1` | Cloacibacillus porcorum | 1197717 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 1197717; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_001712875.1_ASM171287v1` | Candidatus Phytoplasma zeae | 202462 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 202462; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_001747405.1_ASM174740v1` | Chlorobaculum limnaeum | 274537 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 274537; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_001747425.1_ASM174742v1` | Actinoalloteichus hymeniacidonis | 340345 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 340345; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_001886435.1_ASM188643v1` | Alteromonas mediterranea | 314275 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 314275; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_001886815.1_ASM188681v1` | Caldithrix abyssi | 880073 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 880073; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_001889025.1_ASM188902v1` | Allosphingorhabdus lutea | 1913578 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 1913578; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_001953175.1_ASM195317v1` | Dehalogenimonas formicexedens | 1839801 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 1839801; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_002139875.1_ASM213987v1` | Campylobacter vicugnae | 1660076 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 1660076; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_002157165.1_ASM215716v1` | Acidovorax carolinensis | 553814 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 553814; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_002269385.1_ASM226938v1` | Maribacter cobaltidurans | 1178778 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 1178778; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_002302565.1_ASM230256v1` | Capnocytophaga canimorsus | 28188 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 28188; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_002305895.1_ASM230589v1` | Corallococcus macrosporus | 1189310 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 1189310; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_002310495.1_ASM231049v1` | Nibricoccus aquaticus | 2576891 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 2576891; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_002355835.1_ASM235583v1` | Candidatus Endomicrobiellum trichonymphae | 1408204 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 1408204; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_002804025.1_ASM280402v1` | Mesoplasma tabanidae | 219745 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 219745; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_002813575.1_ASM281357v1` | Nostoc flagelliforme | 2038116 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 2038116; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_002872415.1_ASM287241v1` | Bacteriovorax stolpii | 960 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 960; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_003019315.1_ASM301931v1` | Fusobacterium mortiferum | 469616 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 469616; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_003019655.1_ASM301965v1` | Fusobacterium varium | 469618 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 469618; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_003019675.1_ASM301967v1` | Fusobacterium ulcerans | 861 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 861; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_003057965.1_ASM305796v1` | Thermodesulfobium acidiphilum | 1794699 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 1794699; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_003173015.1_ASM317301v1` | Deinococcus irradiatisoli | 2202254 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 2202254; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_003258315.1_ASM325831v1` | Bradymonas sediminis | 1548548 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 1548548; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_003325455.1_ASM332545v1` | Acidipropionibacterium virtanenii | 2057246 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 2057246; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_003330725.1_ASM333072v1` | Acidisarcina polymorpha | 2211140 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 2211140; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_003351905.1_ASM335190v1` | Cardinium endosymbiont of Sogatella furcifera | 650378 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 650378; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_003353065.1_ASM335306v1` | Cyclonatronum proteinivorum | 1457365 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 1457365; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_003355515.1_ASM335551v1` | Aliarcobacter trophiarum | 1032241 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 1032241; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_003366055.1_ASM336605v1` | Candidatus Chromulinivorax destructor | 2066483 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 2066483; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_003491405.1_ASM349140v1` | Salinibacter ruber | 146919 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 146919; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_003544815.1_ASM354481v1` | Arcobacter suis | 663365 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 663365; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_003574315.2_ASM357431v2` | Capsulimonas corticalis | 2219043 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 2219043; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_003721155.3_ASM372115v4` | Janibacter melonis | 262209 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 262209; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_003815875.1_ASM381587v1` | Chryseobacterium lactis | 1241981 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 1241981; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_003966735.1_ASM396673v1` | Desulfovibrio ferrophilus | 241368 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 241368; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_003966975.1_ASM396697v1` | Pelolinea submarina | 913107 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 913107; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_004087915.1_ASM408791v1` | Geovibrio thiophilus | 139438 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 139438; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_004102945.1_ASM410294v1` | Candidatus Velamenicoccus archaeovorus | 1930593 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 1930593; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_004103715.1_ASM410371v1` | Acidilutibacter cellobiosedens | 2507161 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 2507161; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_004214795.1_ASM421479v1` | Malaciobacter pacificus | 1080223 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 1080223; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_004376055.1_ASM437605v2` | Candidatus Magnetomonas plexicatena | 2552947 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 2552947; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_004803815.1_ASM480381v1` | Campylobacter showae | 204 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 204; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_005144435.2_ASM514443v2` | Agrobacterium leguminum | 2792015 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 2792015; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_005697565.1_ASM569756v1` | Candidatus Nanosynbacter featherlites | 2572088 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 2572088; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_006352065.1_ASM635206v1` | Georgenia yuyongxinii | 2589797 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 2589797; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_007475525.1_ASM747552v1` | Methylacidiphilum kamchatkense | 1202785 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 1202785; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_007747655.1_ASM774765v1` | Anatilimnocola aggregata | 2528021 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 2528021; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_007750655.1_ASM775065v1` | Rohdeia mirabilis | 2528008 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 2528008; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_007904085.1_ASM790408v1` | Euhalothece natronophila | 522448 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 522448; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_007990365.1_ASM799036v1` | Leptotrichia trevisanii | 109328 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 109328; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_007990505.1_ASM799050v1` | Pseudoleptotrichia goodfellowii | 157692 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 157692; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_007993755.1_ASM799375v1` | Allomicrovenator marinus | 2600177 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 2600177; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_007993775.1_ASM799377v1` | Persicimonas caeni | 2292766 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 2292766; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_008086545.1_ASM808654v1` | Spiroplasma chinense | 216932 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 216932; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_008087625.1_ASM808762v1` | Aquisphaera giovannonii | 406548 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 406548; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_008087665.1_ASM808766v1` | Mariniblastus fucicola | 980251 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 980251; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_009002475.1_ASM900247v1` | Candidatus Uabimicrobium amorphum | 2596890 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 2596890; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_009258225.1_ASM925822v1` | Sulfurimonas lithotrophica | 2590022 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 2590022; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_009688965.1_ASM968896v1` | Desulfosarcina widdelii | 947919 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 947919; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_009689585.1_ASM968958v1` | Aquipluma nitroreducens | 2010828 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 2010828; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_009731575.1_ASM973157v1` | Akkermansia muciniphila | 239935 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 239935; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_009734425.1_ASM973442v1` | Corynebacterium occultum | 2675219 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 2675219; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_009905215.1_ASM990521v1` | Tichowtungia aerotolerans | 2697043 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 2697043; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_010726765.1_ASM1072676v1` | Mycobacterium gallinarum | 39689 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 39689; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_011751885.2_ASM1175188v2` | Ruicaihuangia jiafui | 2708338 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 2708338; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_013004105.1_ASM1300410v1` | Gemmatimonas groenlandica | 2732249 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 2732249; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_013128195.2_ASM1312819v2` | Frigoriglobus tundricola | 2774151 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 2774151; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_013177655.1_ASM1317765v1` | Campylobacter canadensis | 449520 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 449520; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_013340845.1_ASM1334084v1` | Hydrogenovibrio marinus | 28885 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 28885; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_014274685.1_ASM1427468v1` | Edaphobacter albus | 2763071 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 2763071; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_014337215.1_ASM1433721v1` | Wansuia hejianensis | 2763667 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 2763667; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_014337255.1_ASM1433725v1` | Fusobacterium hominis | 2764326 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 2764326; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_014529625.1_ASM1452962v2` | Candidatus Pelagisphaera phototrophica | 2684113 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 2684113; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_014672715.2_ASM1467271v2` | Dehalogenimonas etheniformans | 1536648 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 1536648; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_014701235.1_ASM1470123v1` | Dissulfurispira thermophila | 2715679 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 2715679; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_015277775.1_ASM1527777v1` | Paludibaculum fermentans | 1473598 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 1473598; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_015356815.2_ASM1535681v2` | Candidatus Rhabdochlamydia porcellionis | 225148 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 225148; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_015500095.1_ASM1550009v1` | Phototrophicus methaneseepsis | 2710758 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 2710758; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_015775515.1_ASM1577551v1` | Atribacter laminatus | 2847778 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 2847778; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_016027375.1_ASM1602737v1` | Clostridium perfringens | 1502 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 1502; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_016446285.1_ASM1644628v1` | Desulfobulbus oligotrophicus | 1909699 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 1909699; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_016724785.1_ASM1672478v1` | Fusobacterium canifelinum | 285729 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 285729; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_017498545.1_ASM1749854v1` | Sulfidibacter corallicola | 2818388 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 2818388; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_017893965.1_ASM1789396v1` | Treponema parvum | 138851 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 138851; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_018069705.1_ASM1806970v1` | Aminithiophilus ramosus | 3029084 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 3029084; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_018127985.1_ASM1812798v1` | Prevotella multiformis | 282402 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 282402; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_018137965.1_ASM1813796v1` | Luteolibacter ambystomatis | 2824561 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 2824561; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_018243235.1_ASM1824323v1` | Ralstonia nicotianae | 3037696 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 3037696; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_018292205.1_ASM1829220v1` | Bacteroides caccae | 997873 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 997873; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_018304765.1_ASM1830476v1` | Candidatus Chloracidobacterium aggregatum | 2821540 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 2821540; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_018343815.1_ASM1834381v1` | Chlamydia crocodili | 2766982 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 2766982; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_018398915.1_ASM1839891v1` | Aminirod propionatiphilus | 3415223 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 3415223; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_018771645.1_ASM1877164v1` | Thermus antranikianii | 88190 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 88190; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_018866205.1_ASM1886620v1` | Candidatus Minimicrobia vallesae | 2841264 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 2841264; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_019048065.1_ASM1904806v1` | Citrobacter farmeri | 67824 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 67824; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_019048105.1_ASM1904810v1` | Providencia rettgeri | 587 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 587; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_019056495.1_ASM1905649v1` | Chlamydia buteonis | 2494525 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 2494525; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_019063885.1_ASM1906388v1` | Geomonas diazotrophica | 2843197 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 2843197; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_019552005.1_ASM1955200v1` | Fusobacterium nucleatum | 851 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 851; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_020091585.1_ASM2009158v1` | Kovacikia minuta | 2653194 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 2653194; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_020549965.1_ASM2054996v1` | Chlamydia pecorum | 85991 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 85991; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_020616335.1_ASM2061633v1` | Corynebacterium uberis | 2883169 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 2883169; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_021057185.1_ASM2105718v1` | Candidatus Absconditicoccus praedator | 2735562 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 2735562; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_021117425.1_ASM2111742v1` | Brevibacterium zhoupengii | 2898795 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 2898795; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_021379005.1_ASM2137900v1` | Leptodesmis sichuanensis | 2835870 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 2835870; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_022026155.1_ASM2202615v1` | Dissulfurimicrobium hydrothermale | 1750598 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 1750598; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_022376835.1_ASM2237683v1` | Brucepastera parasyntrophica | 2880008 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 2880008; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_022813675.1_ASM2281367v1` | Deferrivibrio metallireducens | 2782161 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 2782161; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_022814725.1_ASM2281472v1` | Bdellovibrio reynosensis | 2835041 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 2835041; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_022846135.1_ASM2284613v1` | Pyramidobacter piscolens | 638849 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 638849; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_022869165.1_ASM2286916v1` | Sphaerochaeta associata | 1129264 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 1129264; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_022870745.1_ASM2287074v1` | Silvibacterium acidisoli | 2883246 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 2883246; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_023035795.1_ASM2303579v1` | Borrelia nietonii | 1293576 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 1293576; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_023035875.1_ASM2303587v1` | Borrelia puertoricensis | 2756107 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 2756107; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_023169805.1_ASM2316980v1` | Arboribacterium acericorticis | 2562437 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 2562437; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_023921165.1_ASM2392116v1` | Aquirufa antheringensis | 2516559 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 2516559; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_023935825.1_ASM2393582v1` | Sphingomonas aerolata | 185951 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 185951; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_024205945.1_ASM2420594v1` | Macrococcus equipercicus | 69967 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 69967; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_024338845.1_ASM2433884v1` | Mycoplasma bradburyae | 2963128 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 2963128; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_024584745.1_ASM2458474v1` | Dolichospermum heterosporum | 747523 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 747523; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_025148635.1_ASM2514863v1` | Butyricimonas virosa | 544645 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 544645; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_025273655.1_ASM2527365v1` | Candidatus Mycosynbacter amalyticus | 2665156 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 2665156; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_025558325.1_ASM2555832v1` | Staphylococcus schleiferi | 1295 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 1295; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_025725725.1_ASM2572572v1` | Brachybacterium huguangmaarense | 1652028 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 1652028; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_025997855.1_ASM2599785v1` | Deinococcus aetherius | 200252 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 200252; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_025997995.1_ASM2599799v1` | Sediminibacterium planctonicum | 575378 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 575378; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_026184435.1_ASM2618443v1` | Zhurongbacter thermophilus | 2681877 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 2681877; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_026240635.1_ASM2624063v1` | Trichothermofontia sichuanensis | 2016064 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 2016064; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_026965555.1_ASM2696555v1` | Nannocystis punicea | 2995304 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 2995304; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_027594505.1_ASM2759450v1` | Tepidiforma flava | 3004094 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 3004094; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_027625475.1_ASM2762547v1` | Sphingomonas abietis | 3012344 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 3012344; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_027923555.1_ASM2792355v1` | Vulcanimicrobium alpinum | 3016050 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 3016050; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_028622095.1_ASM2862209v1` | Deinococcus aquaticus | 328692 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 328692; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_028728065.1_ASM2872806v1` | Lentisphaera profundi | 1658616 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 1658616; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_028885455.1_ASM2888545v1` | Pseudoalteromonas rhizosphaerae | 2518973 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 2518973; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_029023725.1_ASM2902372v1` | Actomonas lomoniglobus | 2909285 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 2909285; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_029677705.1_ASM2967770v1` | Mycoplasma feriruminatoris | 1179777 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 1179777; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_029919255.1_ASM2991925v1` | cyanobacterium endosymbiont of Epithemia clementina EcSB | 3034674 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 3034674; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_030169785.1_ASM3016978v1` | Bdellovibrio bacteriopascens | 3044166 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 3044166; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_030253535.1_ASM3025353v1` | Candidatus Southlakia epibionticum | 3043284 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 3043284; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_030295385.1_ASM3029538v1` | Geothrix oryzae | 2927975 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 2927975; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_030389965.1_ASM3038996v1` | Candidatus Chlorohelix allophototropha | 3003348 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 3003348; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_030438115.2_ASM3043811v2` | Borreliella japonica | 34095 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 34095; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_030439285.2_ASM3043928v2` | [Borrelia] turdi | 57863 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 57863; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_030444335.1_ASM3044433v1` | Polaribacter huanghezhanensis | 1354726 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 1354726; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_030503875.2_ASM3050387v2` | Marivirga arenosa | 3059076 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 3059076; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_030520105.1_ASM3052010v1` | Fervidibacter sacchari | 1448929 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 1448929; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_030585795.1_ASM3058579v1` | Leuconostoc carnosum | 1252 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 1252; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_030623445.1_ASM3062344v1` | Pseudogemmatithrix spongiicola | 3062599 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 3062599; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_030866825.1_ASM3086682v1` | Rathayibacter soli | 3144168 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 3144168; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_032050975.1_ASM3205097v1` | Candidatus Nitrospira allomarina | 3020900 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 3020900; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_032051675.1_ASM3205167v1` | Candidatus Nitrospira neomarina | 3020899 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 3020899; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_032102965.1_ASM3210296v1` | Demequina capsici | 3075620 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 3075620; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_032595795.2_ASM3259579v2` | Borreliella americana | 478807 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 478807; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_032698475.1_ASM3269847v1` | Xanthomonas dyei | 743699 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 743699; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_033547035.1_ASM3354703v1` | Archangium violaceum | 83451 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 83451; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_033878955.1_ASM3387895v1` | Coraliomargarita algicola | 3092156 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 3092156; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_033969665.2_ASM3396966v2` | Borreliella sinica | 87162 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 87162; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_033971325.1_ASM3397132v1` | Pampinifervens diazotrophicum | 1632018 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 1632018; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_034259205.1_ASM3425920v1` | Peredibacter starrii | 28202 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 28202; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_034424215.1_ASM3442421v1` | Paraniabella yanshanensis | 577386 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 577386; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_037076345.1_ASM3707634v1` | Haloferula helveola | 490095 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 490095; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_037076455.1_ASM3707645v1` | Desulfoferula mesophila | 3058419 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 3058419; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_037081835.1_ASM3708183v1` | Spirobacillus cienkowskii | 495820 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 495820; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_037157775.1_ASM3715777v1` | Pendulispora brunnea | 2905690 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 2905690; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_037477945.1_MHOCP` | Neomoorella humiferrea | 676965 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 676965; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_038095605.1_ASM3809560v1` | Gaopeijia maritima | 3119007 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 3119007; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_038396615.1_ASM3839661v1` | Thermatribacter velox | 3039681 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 3039681; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_038400315.1_ASM3840031v1` | Kueselia aquadivae | 3139350 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 3139350; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_038443385.1_ASM3844338v1` | Stieleria varia | 2528005 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 2528005; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_038725785.1_ASM3872578v1` | Termitidicoccus mucosus | 1184151 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 1184151; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_039555295.1_ASM3955529v1` | Helicobacter mastomyrinus | 287948 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 287948; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_039645985.1_ASM3964598v1` | Sulfonitrofixus jiaomeiensis | 3131938 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 3131938; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_039994875.1_ASM3999487v1` | Occallatibacter bavaricus | 3036704 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 3036704; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_040256395.1_GCC_DsonKR87_v1` | Deinococcus sonorensis | 694439 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 694439; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_040450955.1_ASM4045095v1` | Thermodesulfovibrio autotrophicus | 3118333 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 3118333; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_040451425.1_ASM4045142v1` | Thermodesulfovibrio obliviosus | 3118332 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 3118332; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_040564635.1_ASM4056463v1` | Planktothricoides raciborskii | 2792601 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 2792601; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_040833985.2_ASM4083398v1` | Leptospira mgodei | 3235172 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 3235172; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_041026625.1_ASM4102662v1` | Metamycoplasma spumans | 92406 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 92406; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_041154365.1_ASM4115436v1` | Fidelibacter multiformis | 3377529 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 3377529; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_041200005.1_ASM4120000v1` | Leptotrichia rugosa | 3239302 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 3239302; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_041505935.1_Fervidobacterium_thailandense_FC2004` | Fervidobacterium thailandense | 1008305 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 1008305; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_042033325.1_ASM4203332v1` | Sneathia sanguinegens | 40543 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 40543; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_042847215.1_ASM4284721v1` | Gilvimarinus japonicus | 1796469 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 1796469; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_044336245.1_ASM4433624v1` | Schlesneria sphaerica | 3373610 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 3373610; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_044789115.1_ASM4478911v1` | Lotipaludicoccus stagni | 3034021 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 3034021; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_045784255.1_ASM4578425v1` | Deinococcus altitudinis | 468914 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 468914; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_045784305.1_DRML_24019` | Deinococcus radiomollis | 468916 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 468916; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_045864955.1_ASM4586495v1` | Atrimonas thermophila | 3064161 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 3064161; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_046267995.1_ASM4626799v1` | Helcococcus bovis | 3153252 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 3153252; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_046774955.1_ASM4677495v1` | Deinococcus yunweiensis | 367282 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 367282; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_046845665.1_ASM4684566v1` | Rubrivirga aquatilis | 3391575 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 3391575; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_046845755.1_ASM4684575v1` | Rubrivirga halophila | 3391614 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 3391614; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_047159735.1_ASM4715973v1` | Minisyncoccus archaeiphilus | 3238481 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 3238481; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_048565275.1_ASM4856527v1` | Desulfosarcina variabilis | 859321 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 859321; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_049669535.1_sGff_cult_auto_v1` | Spiroplasma endosymbiont of Glossina fuscipes fuscipes | 2004463 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 2004463; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_049676565.1_Dip510` | Elusimicrobium posterum | 3116653 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 3116653; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_050039385.1_ASM5003938v1` | Myxococcus faecalis | 3115646 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 3115646; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_050267035.1_ASM5026703v1` | Candidatus Caldatribacterium saccharofermentans | 1454753 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 1454753; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_050389785.1_ASM5038978v1` | Parelusimicrobium proximum | 3228953 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 3228953; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_050500235.1_ASM5050023v1` | Microcystis protocystis | 3422496 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 3422496; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_050563455.1_ASM5056345v1` | Crocosphaera watsonii | 165597 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 165597; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_050780155.1_QI0013` | Bilophila wadsworthia | 35833 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 35833; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_051175385.1_ASM5117538v1` | Tundrisphaera lichenicola | 2029860 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 2029860; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_052161605.1_ASM5216160v1` | Brakhagea baltica | 3454458 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 3454458; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_054054625.1_ASM5405462v1` | Pigmentibacter ruber | 2683196 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 2683196; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_054917915.1_ASM5491791v1` | Sphingomonas citri | 2862499 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 2862499; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_055378285.1_ASM5537828v1` | Desulfurispirillum alkaliphilum | 393030 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 393030; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_055378975.1_ASM5537897v1` | Desulfitobacterium elongatum | 3375637 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 3375637; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_055379415.1_ASM5537941v1` | Leucogemmata humicola | 3417965 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 3417965; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_055379825.1_ASM5537982v1` | Longimonas haloalkaliphila | 1807135 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 1807135; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_055379935.1_ASM5537993v1` | Natronotalea proteinilytica | 1807134 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 1807134; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_055380895.1_ASM5538089v1` | Thermovibrio ruber | 171869 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 171869; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_055383505.1_ASM5538350v1` | Piscibacter tructae | 3455867 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 3455867; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_055384375.1_ASM5538437v1` | Thiofractor thiocamini | 1095805 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 1095805; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_055384745.1_ASM5538474v1` | Desulfurobacterium pacificum | 240166 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 240166; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_055385515.1_ASM5538551v1` | Hydrogenivirga okinawensis | 345592 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 345592; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_055387075.1_ASM5538707v1` | Marinitoga okinawensis | 389480 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 389480; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_055387105.1_ASM5538710v1` | Caldifarcimen microaerophilum | 3388968 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 3388968; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_055387525.1_ASM5538752v1` | Caldilinea tarbellica | 859243 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 859243; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_055388355.1_ASM5538835v1` | Prosthecobacter fluviatilis | 445931 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 445931; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_055390195.1_ASM5539019v1` | Hydrogenobacter subterraneus | 92642 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 92642; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_055390205.1_ASM5539020v1` | Mesotoga infera | 1236046 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 1236046; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_055390565.1_ASM5539056v1` | Aminobacterium thunnarium | 1482737 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 1482737; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_055390845.1_ASM5539084v1` | Thermodesulfovibrio aggregans | 86166 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 86166; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_055390895.1_ASM5539089v1` | Lactivibrio alcoholicus | 760261 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 760261; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_055391325.1_ASM5539132v1` | Defluviitoga tunisiensis | 1006576 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 1006576; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_055396605.1_ASM5539660v1` | Clostridium budayi | 83338 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 83338; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_055694045.1_ASM5569404v1` | Fusobacterium abscessus | 3475676 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 3475676; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_055763695.1_ASM5576369v1` | Nocardiopsis moroccensis | 3461151 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 3461151; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_056388415.1_ASM5638841v1` | Mycoplasmopsis alligatoris | 747682 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 747682; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_058362665.1_ASM5836266v1` | Corallococcus interemptor | 2316720 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 2316720; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_059632245.1_ASM5963224v1` | Bythopirellula superficialaquae | 3454465 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 3454465; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_059696725.1_ASM5969672v1` | Thermogemmatispora argillosa | 2045280 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 2045280; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_059705355.1_ASM5970535v1` | Sulfurimonas lacustris | 3686813 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 3686813; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_059726255.1_ASM5972625v1` | Koromonas humicola | 3421943 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 3421943; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_059726345.1_ASM5972634v1` | Humiphilus diazotrophicus | 3421940 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 3421940; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_060026435.1_ASM6002643v1` | Sorangium bulgaricum | 2093423 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 2093423; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_060027055.1_ASM6002705v1` | Racemicystis crocea | 1707966 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 1707966; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_060030665.1_ASM6003066v1` | Vitiosangium cumulatum | 1867796 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 1867796; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_060030705.1_ASM6003070v1` | Simulacricoccus ruber | 2303410 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 2303410; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_060030825.1_ASM6003082v1` | Sorangium ambruticinum | 2093422 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 2093422; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_060031025.1_ASM6003102v1` | Chondromyces robustus | 83448 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 83448; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_060138195.1_ASM6013819v1` | Nannocystis konarekensis | 2066023 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 2066023; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_900169565.1_NSJP_Ch1` | Nitrospira japonica | 1325564 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 1325564; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_900184705.1_CAMBI-1` | Candidatus Brevifilum fermentans | 1986204 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 1986204; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_900217235.1_DPRO_PRJEB22548_v1` | Pseudodesulfovibrio profundus | 57320 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 57320; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_900239975.1_Chlamydia_sp._S15-834K` | Chlamydia poikilotherma | 1967783 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 1967783; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_900475595.1_44343_G01` | Streptococcus uberis | 1349 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 1349; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_900660485.1_50648_E01-3` | Mesomycoplasma neurolyticum | 2120 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 2120; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_900660665.1_57612_D01-3` | Mycoplasmopsis maculosa | 114885 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 114885; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_901538265.1_ASM90153826v1` | Gemmata massiliana | 1210884 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 1210884; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_902806995.1_15-2067_O50` | Chlamydiifrater volucris | 2681470 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 2681470; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_907165195.1_Cryomorphaceae_CECT30217` | Parvicella tangerina | 2829795 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 2829795; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_918320885.1_Harlequin_endosymbiont1` | unidentified bacterial endosymbiont | 2355 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 2355; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_946900835.1_NSPWAT` | Nitrospina watsonii | 1323948 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 1323948; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_947241125.1_NsDNF` | Nitrospira tepida | 2973512 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 2973512; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_949774925.1_MFUM-new` | Candidatus Methylacidiphilum fumarolicum | 591154 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 591154; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_963264975.1_NZ_220527_assembly_043` | Fusobacterium periodonticum | 860 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 860; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_964030865.1_idDioLine1.Spiroplasma_sp_1.1` | Spiroplasma endosymbiont of Dioctria linearis | 3066290 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 3066290; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_964030885.1_icBarPell3.Wolbachia_sp_1.1` | Wolbachia endosymbiont (group A) of Barypeithes pellucidus | 3139322 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 3139322; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_964642095.1_Ff3` | Fretibacterium fastidiosum | 651822 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 651822; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_965152565.1_NS_bockiana` | Candidatus Nitrospira bockiana | 463288 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 463288; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_976986435.1_NS_calida2` | Nitrospira calida | 876851 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 876851; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `GCF_978936825.1_DTS69` | Dissulfuribacter thermophilus | 1156395 | `not_available` | `catalog_missing` | sin catálogo ComplexPortal local para taxid 1156395; mapeo UniProt→locus: functional.tsv (UniProt→gene_id) |
| `ecoli_k12` | Escherichia coli str. K-12 substr. MG1655 | 83333 | `data/complexes/gold_k12.tsv` | `ready` | catálogo gold local para taxid 83333; taxid leído del catálogo ComplexPortal local; mapeo UniProt→locus: genes.tsv acc (UniProt→locus_tag/gene_id) |

## Regenerar

```bash
python scripts/complexes/inventory_complexportal_400.py
```

## Lectura

- Un `catalog_missing` solo dice que no hay catálogo local. No prueba que ComplexPortal tenga o deje de tener una entrada para ese taxid; eso se comprueba contra el listado remoto (ver `docs/COMPLEXPORTAL_400.md`).
- Los `model_<genoma>.complexes.tsv` y `transfer_<genoma>.tsv` no son catálogos ComplexPortal y no se cuentan como tales.
- `raw_only` todavía no sirve para anotar: falta el mapeo UniProt→gen y, en su caso, el filtro de complejos incompletos.
