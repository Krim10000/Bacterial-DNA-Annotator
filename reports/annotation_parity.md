# Paridad de anotación base

**Generado:** `2026-11-24T13:20:57Z`  
**Genoma:** `GCF_000006685.1_ASM668v1`  
**Estado:** **evaluado**

## Alcance

Comparación de anotación base contra un GFF3 de referencia local. Mide CDS/tRNA/rRNA/ncRNA, coordenadas, productos/hipotéticas, ORFs pequeños (<150 aa) y diferencias explicables. No mide identidad de secuencia ni función.

```bash
python tools/annotation_parity.py --genome GCF_000006685.1_ASM668v1
```

## Entradas

| rol | ruta | sha256 | CDS | productos |
|---|---|---|---|---|
| nuestra anotación | `data/annotations_ab_initio/GCF_000006685.1_ASM668v1/annotation.gff3` | f08245f2c17e | 902 | TSV vecino: functional.tsv |
| referencia | `data/annotations/GCF_000006685.1_ASM668v1.gff.gz` | 50ece152af6e | 891 | atributo product del GFF |

Parámetros: umbral de solapamiento recíproco **0,500**, alerta de borde **30 nt**, ORF pequeño **<150 aa**, pseudo-CDS excluidos: **sí**.

## Contigs

| lado | n | contigs |
|---|---|---|
| referencia | 2 | NC_002182.1, NC_002620.2 |
| nuestra | 2 | NC_002182.1, NC_002620.2 |
| compartidos | 2 | NC_002182.1, NC_002620.2 |
| solo referencia | 0 | — |
| solo nuestra | 0 | — |

## Conteos por clase

| clase | nuestra | referencia | delta | en alcance (nuestra/ref) | estado |
|---|---|---|---|---|---|
| CDS | 902 | 891 | 11 | 902 / 891 | evaluado |
| tRNA | 37 | 37 | 0 | 37 / 37 | evaluado |
| rRNA | 6 | 6 | 0 | 6 / 6 | evaluado |
| ncRNA | 0 | 3 | -3 | 0 / 3 | ausente en nuestra anotación |

## Paridad de coordenadas

| clase | emparejados | recall | precision | F1 | exactos | % exactos | alertas de borde | largo discordante |
|---|---|---|---|---|---|---|---|---|
| CDS | 884 | 0,992 | 0,980 | 0,986 | 803 | 90,8 % | 37 | 4 |
| tRNA | 37 | 1,000 | 1,000 | 1,000 | 37 | 100,0 % | 0 | 0 |
| rRNA | 6 | 1,000 | 1,000 | 1,000 | 6 | 100,0 % | 0 | 0 |
| ncRNA | — | — | — | — | — | — | — | — |

Micro recall **0,993**, micro precision **0,981**, micro F1 **0,987** (solo clases evaluadas); tasa de coordenadas exactas **91,3 %** sobre 927 emparejamientos.

## Productos e hipotéticas

| lado | fuente | CDS | con product | sin product | hipotéticas | % hypothetical |
|---|---|---|---|---|---|---|
| nuestra | TSV vecino: functional.tsv | 902 | 638 | 264 | 58 | 9,1 % |
| referencia | atributo product del GFF | 891 | 891 | 0 | 154 | 17,3 % |

Pares CDS con producto en ambos lados: **630**; coincidencia literal **20,9 %**; Jaccard de tokens medio **0,460**, mediana **0,333**.
Desacuerdos hipotética/no-hipotética: **37**; posibles nombres distintos: **268**.

## ORFs pequeños (<150 aa)

| lado | CDS | <150 aa | % del total |
|---|---|---|---|
| nuestra | 902 | 156 | 17,3 % |
| referencia | 891 | 153 | 17,2 % |

Emparejados pequeños: **144**; ausentes en nuestra anotación: **7**; extra en nuestra anotación: **9**.

## ncRNA

| lado | total | tipos |
|---|---|---|
| nuestra | 0 | — |
| referencia | 3 | RNase_P_RNA×1, SRP_RNA×1, tmRNA×1 |
| estado | ausente en nuestra anotación | — |

## Diferencias accionables

**Total:** **98**. En este Markdown se muestran 40; el JSON incluye los 98.

| categoría | n |
|---|---|
| boundary_shift | 37 |
| hypothetical_disagreement | 37 |
| extra_cds | 12 |
| missing_cds | 7 |
| length_discordance | 4 |
| missing_ncrna | 1 |

| severidad | categoría | clase | referencia | nuestra | detalle |
|---|---|---|---|---|---|
| high | boundary_shift | CDS | NC_002620.2:833828-837736(+) | NC_002620.2:833354-837736(+) | Bordes desplazados: 5' 474 nt, 3' 0 nt respecto de la referencia (umbral 30 nt). |
| high | boundary_shift | CDS | NC_002620.2:92932-94260(-) | NC_002620.2:92932-94506(-) | Bordes desplazados: 5' 246 nt, 3' 0 nt respecto de la referencia (umbral 30 nt). |
| high | boundary_shift | CDS | NC_002620.2:378645-380018(-) | NC_002620.2:378645-379799(-) | Bordes desplazados: 5' 219 nt, 3' 0 nt respecto de la referencia (umbral 30 nt). |
| high | boundary_shift | CDS | NC_002620.2:622331-623551(-) | NC_002620.2:622331-623368(-) | Bordes desplazados: 5' 183 nt, 3' 0 nt respecto de la referencia (umbral 30 nt). |
| high | boundary_shift | CDS | NC_002620.2:887950-889057(-) | NC_002620.2:887950-888612(-) | Bordes desplazados: 5' 445 nt, 3' 0 nt respecto de la referencia (umbral 30 nt). |
| high | boundary_shift | CDS | NC_002620.2:173453-174511(+) | NC_002620.2:173588-174511(+) | Bordes desplazados: 5' 135 nt, 3' 0 nt respecto de la referencia (umbral 30 nt). |
| high | boundary_shift | CDS | NC_002620.2:177525-178328(-) | NC_002620.2:177525-178193(-) | Bordes desplazados: 5' 135 nt, 3' 0 nt respecto de la referencia (umbral 30 nt). |
| high | boundary_shift | CDS | NC_002620.2:658528-659187(-) | NC_002620.2:658528-658956(-) | Bordes desplazados: 5' 231 nt, 3' 0 nt respecto de la referencia (umbral 30 nt). |
| high | boundary_shift | CDS | NC_002620.2:690363-690899(-) | NC_002620.2:690363-690779(-) | Bordes desplazados: 5' 120 nt, 3' 0 nt respecto de la referencia (umbral 30 nt). |
| high | boundary_shift | CDS | NC_002620.2:413993-414436(-) | NC_002620.2:413993-414277(-) | Bordes desplazados: 5' 159 nt, 3' 0 nt respecto de la referencia (umbral 30 nt). |
| high | missing_ncrna | ncRNA | — | — | La referencia anota 3 ncRNA y nuestra anotación no tiene ninguno. |
| medium | hypothetical_disagreement | CDS | NC_002620.2:883492-886515(+) | NC_002620.2:883492-886515(+) | La referencia y nuestra anotación discrepan en si el producto es hipotético: ref=type III secretion system actin-recruiting effector Tarp; ours=Uncharacterized protein TC_0741. |
| medium | boundary_shift | CDS | NC_002620.2:610516-613314(+) | NC_002620.2:610558-613314(+) | Bordes desplazados: 5' 42 nt, 3' 0 nt respecto de la referencia (umbral 30 nt). |
| medium | hypothetical_disagreement | CDS | NC_002620.2:760882-762996(-) | NC_002620.2:760882-762996(-) | La referencia y nuestra anotación discrepan en si el producto es hipotético: ref=thioredoxin domain-containing protein; ours=Uncharacterized protein YyaL. |
| medium | boundary_shift | CDS | NC_002620.2:846547-848412(+) | NC_002620.2:846499-848412(+) | Bordes desplazados: 5' 48 nt, 3' 0 nt respecto de la referencia (umbral 30 nt). |
| medium | hypothetical_disagreement | CDS | NC_002620.2:321624-323441(+) | NC_002620.2:321624-323441(+) | La referencia y nuestra anotación discrepan en si el producto es hipotético: ref=type III secretion system early effector TepP; ours=Uncharacterized protein TC_0268. |
| medium | hypothetical_disagreement | CDS | NC_002620.2:288801-290606(+) | NC_002620.2:288801-290606(+) | La referencia y nuestra anotación discrepan en si el producto es hipotético: ref=protease-like activity factor CPAF; ours=Uncharacterized protein TC_0248. |
| medium | boundary_shift | CDS | NC_002620.2:203306-205102(-) | NC_002620.2:203306-205045(-) | Bordes desplazados: 5' 57 nt, 3' 0 nt respecto de la referencia (umbral 30 nt). |
| medium | hypothetical_disagreement | CDS | NC_002620.2:418784-420523(-) | NC_002620.2:418784-420523(-) | La referencia y nuestra anotación discrepan en si el producto es hipotético: ref=menaquinone biosynthesis decarboxylase; ours=Uncharacterized protein TC_0359. |
| medium | hypothetical_disagreement | CDS | NC_002620.2:671920-673596(+) | NC_002620.2:671920-673596(+) | La referencia y nuestra anotación discrepan en si el producto es hipotético: ref=inclusion membrane protein IncM; ours=Uncharacterized protein TC_0561. |
| medium | boundary_shift | CDS | NC_002620.2:791260-792894(-) | NC_002620.2:791260-792840(-) | Bordes desplazados: 5' 54 nt, 3' 0 nt respecto de la referencia (umbral 30 nt). |
| medium | boundary_shift | CDS | NC_002620.2:421064-422647(-) | NC_002620.2:421064-422677(-) | Bordes desplazados: 5' 30 nt, 3' 0 nt respecto de la referencia (umbral 30 nt). |
| medium | hypothetical_disagreement | CDS | NC_002620.2:1071329-1072918(-) | NC_002620.2:1071329-1072918(-) | La referencia y nuestra anotación discrepan en si el producto es hipotético: ref=FAD-dependent thymidylate synthase; ours=UPF0159 protein TC_0921. |
| medium | boundary_shift | CDS | NC_002620.2:739987-741519(+) | NC_002620.2:740059-741519(+) | Bordes desplazados: 5' 72 nt, 3' 0 nt respecto de la referencia (umbral 30 nt). |
| medium | boundary_shift | CDS | NC_002620.2:659189-660700(+) | NC_002620.2:659249-660700(+) | Bordes desplazados: 5' 60 nt, 3' 0 nt respecto de la referencia (umbral 30 nt). |
| medium | boundary_shift | CDS | NC_002620.2:566521-567933(-) | NC_002620.2:566521-567993(-) | Bordes desplazados: 5' 60 nt, 3' 0 nt respecto de la referencia (umbral 30 nt). |
| medium | hypothetical_disagreement | CDS | NC_002620.2:607956-609428(+) | NC_002620.2:607956-609428(+) | La referencia y nuestra anotación discrepan en si el producto es hipotético: ref=sodium-dependent transporter; ours=Uncharacterized sodium-dependent transporter MJ1319. |
| medium | boundary_shift | CDS | NC_002620.2:927596-929014(+) | NC_002620.2:927665-929014(+) | Bordes desplazados: 5' 69 nt, 3' 0 nt respecto de la referencia (umbral 30 nt). |
| medium | boundary_shift | CDS | NC_002620.2:778341-779669(+) | NC_002620.2:778254-779669(+) | Bordes desplazados: 5' 87 nt, 3' 0 nt respecto de la referencia (umbral 30 nt). |
| medium | hypothetical_disagreement | CDS | NC_002620.2:593061-594461(+) | NC_002620.2:593061-594461(+) | La referencia y nuestra anotación discrepan en si el producto es hipotético: ref=amino acid permease; ours=Uncharacterized transporter lpg1691. |
| medium | hypothetical_disagreement | CDS | NC_002620.2:127321-128688(+) | NC_002620.2:127321-128688(+) | La referencia y nuestra anotación discrepan en si el producto es hipotético: ref=alanine/glycine:cation symporter family protein; ours=Uncharacterized transporter HI_0883. |
| medium | boundary_shift | CDS | NC_002620.2:1062106-1063398(-) | NC_002620.2:1062106-1063440(-) | Bordes desplazados: 5' 42 nt, 3' 0 nt respecto de la referencia (umbral 30 nt). |
| medium | hypothetical_disagreement | CDS | NC_002620.2:1004471-1005781(+) | NC_002620.2:1004471-1005781(+) | La referencia y nuestra anotación discrepan en si el producto es hipotético: ref=type III secretion system translocon subunit SctE; ours=Uncharacterized protein TC_0868. |
| medium | hypothetical_disagreement | CDS | NC_002620.2:1037164-1038396(-) | NC_002620.2:1037164-1038396(-) | La referencia y nuestra anotación discrepan en si el producto es hipotético: ref=DUF1343 domain-containing protein; ours=Peptidoglycan beta-N-acetylmuramidase NamZ. |
| medium | boundary_shift | CDS | NC_002620.2:504322-505440(-) | NC_002620.2:504322-505506(-) | Bordes desplazados: 5' 66 nt, 3' 0 nt respecto de la referencia (umbral 30 nt). |
| medium | boundary_shift | CDS | NC_002620.2:848610-849719(+) | NC_002620.2:848544-849719(+) | Bordes desplazados: 5' 66 nt, 3' 0 nt respecto de la referencia (umbral 30 nt). |
| medium | hypothetical_disagreement | CDS | NC_002620.2:800547-801686(-) | NC_002620.2:800547-801686(-) | La referencia y nuestra anotación discrepan en si el producto es hipotético: ref=CT392 family protein; ours=Uncharacterized protein TC_0671. |
| medium | length_discordance | CDS | NC_002620.2:887950-889057(-) | NC_002620.2:887950-888612(-) | Largo discordante: referencia 1110 nt vs nuestra 663 nt (ratio 0.5973). |
| medium | hypothetical_disagreement | CDS | NC_002620.2:327372-328463(-) | NC_002620.2:327372-328463(-) | La referencia y nuestra anotación discrepan en si el producto es hipotético: ref=IncV family inclusion membrane protein; ours=Uncharacterized protein TC_0273. |
| medium | hypothetical_disagreement | CDS | NC_002620.2:998746-999816(-) | NC_002620.2:998746-999816(-) | La referencia y nuestra anotación discrepan en si el producto es hipotético: ref=M24 family metallopeptidase; ours=Uncharacterized peptidase YqhT. |

## Revisión (productos y pseudogenes)

**Total:** **274**. Son diferencias de nomenclatura o limitaciones de la referencia; no se cuentan como errores de anotación. En este Markdown se muestran 20; el JSON guarda 200.

| categoría | referencia | nuestra | detalle |
|---|---|---|---|
| product_mismatch | NC_002620.2:506330-516097(+) | NC_002620.2:506330-516097(+) | Productos con poco solapamiento de tokens (Jaccard 0.17): ref=LifA/Efa1-related large cytotoxin; ours=Cytotoxin-L. |
| product_mismatch | NC_002620.2:526720-536376(+) | NC_002620.2:526720-536376(+) | Productos con poco solapamiento de tokens (Jaccard 0.17): ref=LifA/Efa1-related large cytotoxin; ours=Cytotoxin-L. |
| product_mismatch | NC_002620.2:729938-735301(+) | NC_002620.2:729938-735301(+) | Productos con poco solapamiento de tokens (Jaccard 0.00): ref=excinuclease ABC subunit UvrA; ours=UvrABC system protein A. |
| product_mismatch | NC_002620.2:828166-833175(+) | NC_002620.2:828166-833175(+) | Productos con poco solapamiento de tokens (Jaccard 0.33): ref=polymorphic outer membrane protein middle domain-containing protein; ours=Probable outer membrane protein PmpB. |
| product_mismatch | NC_002620.2:229856-234418(+) | NC_002620.2:229856-234418(+) | Productos con poco solapamiento de tokens (Jaccard 0.33): ref=polymorphic outer membrane protein middle domain-containing protein; ours=Probable outer membrane protein PmpD. |
| product_mismatch | NC_002620.2:833828-837736(+) | NC_002620.2:833354-837736(+) | Productos con poco solapamiento de tokens (Jaccard 0.33): ref=polymorphic outer membrane protein middle domain-containing protein; ours=Probable outer membrane protein PmpC. |
| product_mismatch | NC_002620.2:977680-981279(+) | NC_002620.2:977680-981279(+) | Productos con poco solapamiento de tokens (Jaccard 0.11): ref=DEAD/DEAH box helicase; ours=Putative ATP-dependent helicase MPN_020. |
| product_mismatch | NC_002620.2:94618-98109(+) | NC_002620.2:94618-98109(+) | Productos con poco solapamiento de tokens (Jaccard 0.00): ref=DEAD/DEAH box helicase; ours=TATA-binding protein-associated factor MOT1. |
| product_mismatch | NC_002620.2:60401-63643(+) | NC_002620.2:60401-63643(+) | Productos con poco solapamiento de tokens (Jaccard 0.11): ref=penicillin-binding transpeptidase domain-containing protein; ours=Peptidoglycan D,D-transpeptidase MrdA. |
| product_mismatch | NC_002620.2:7797-10877(-) | NC_002620.2:7797-10877(-) | Productos con poco solapamiento de tokens (Jaccard 0.14): ref=exodeoxyribonuclease V subunit beta; ours=RecBCD enzyme subunit RecB. |
| product_mismatch | NC_002620.2:308277-311354(-) | NC_002620.2:308277-311354(-) | Productos con poco solapamiento de tokens (Jaccard 0.33): ref=polymorphic outer membrane protein middle domain-containing protein; ours=Probable outer membrane protein PmpF. |
| product_mismatch | NC_002620.2:311529-314489(+) | NC_002620.2:311529-314489(+) | Productos con poco solapamiento de tokens (Jaccard 0.33): ref=polymorphic outer membrane protein middle domain-containing protein; ours=Probable outer membrane protein PmpG. |
| product_mismatch | NC_002620.2:314520-317462(+) | NC_002620.2:314520-317462(+) | Productos con poco solapamiento de tokens (Jaccard 0.33): ref=polymorphic outer membrane protein middle domain-containing protein; ours=Probable outer membrane protein PmpH. |
| product_mismatch | NC_002620.2:305344-308274(-) | NC_002620.2:305344-308274(-) | Productos con poco solapamiento de tokens (Jaccard 0.33): ref=polymorphic outer membrane protein middle domain-containing protein; ours=Probable outer membrane protein PmpE. |
| product_mismatch | NC_002620.2:825098-828028(+) | NC_002620.2:825098-828028(+) | Productos con poco solapamiento de tokens (Jaccard 0.33): ref=polymorphic outer membrane protein middle domain-containing protein; ours=Probable outer membrane protein PmpA. |
| product_mismatch | NC_002620.2:247799-250726(+) | NC_002620.2:247799-250726(+) | Productos con poco solapamiento de tokens (Jaccard 0.00): ref=insulinase family protein; ours=Presequence protease 2, chloroplastic/mitochondrial. |
| product_mismatch | NC_002620.2:221882-224701(-) | NC_002620.2:221882-224701(-) | Productos con poco solapamiento de tokens (Jaccard 0.00): ref=insulinase family protein; ours=Insulin-degrading enzyme. |
| product_mismatch | NC_002620.2:50282-53038(+) | NC_002620.2:50282-53038(+) | Productos con poco solapamiento de tokens (Jaccard 0.12): ref=secretin N-terminal domain-containing protein; ours=Secretin GspD 2. |
| product_mismatch | NC_002620.2:382015-384732(+) | NC_002620.2:382015-384732(+) | Productos con poco solapamiento de tokens (Jaccard 0.25): ref=2-oxoglutarate dehydrogenase E1 component; ours=Multifunctional 2-oxoglutarate metabolism enzyme. |
| product_mismatch | NC_002620.2:318745-321348(-) | NC_002620.2:318745-321348(-) | Productos con poco solapamiento de tokens (Jaccard 0.33): ref=polymorphic outer membrane protein middle domain-containing protein; ours=Probable outer membrane protein PmpI. |

## Límites

- La comparación es de anotación (coordenadas, clases y etiquetas), no de función, homología ni identidad de secuencia.
- El emparejamiento es 1:1 voraz por solapamiento recíproco; no es un óptimo bipartito ni un alineamiento global.
- Las coordenadas GFF3 son 1-based inclusive. No se hace liftover ni se alinean secuencias.
- El largo en aminoácidos es una estimación desde el ADN (span/3 menos un codón de parada), no una traducción.
- Las diferencias de nombre de producto entre bases de datos no son errores de anotación; se separan como 'revisión'.
- No se descargan referencias ni se ejecutan Bakta/Prokka/DFAST; solo se leen GFF locales.
- Un GFF de referencia puede omitir tRNA/rRNA/ncRNA o incluir pseudogenes; 'reference_absent' no indica un error nuestro.
- La ausencia de referencia no es un resultado: el estado es 'not_evaluable' y no se calculan precision/recall.
- Nuestros productos vienen de un TSV funcional (p. ej. Swiss-Prot), no del GFF; la comparación de nombres con la referencia es orientativa.

## Reproducir

```bash
python tools/annotation_parity.py --genome GCF_000006685.1_ASM668v1
.venv/bin/python -m pytest tests/test_annotation_parity.py -q
```

Este informe no declara superioridad entre anotadores. Sin una verdad independiente y ejecuciones pareadas, solo describe diferencias.

