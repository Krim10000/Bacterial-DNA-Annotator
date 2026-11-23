# Catálogos de complejos disponibles

## Lo que hay hoy

| catálogo | especie | estado local | uso |
|---|---|---|---|
| ComplexPortal `83333` | *E. coli* K-12 | `data/complexes/raw/complexportal_83333.tsv` | gold K-12 |
| STRING `511145` | *E. coli* K-12 | `data/complexes/raw/string_511145_full.txt.gz` | canales de interacción |
| gold K-12 | *E. coli* K-12 | `data/complexes/gold_k12.tsv` | 273 complejos |
| transferencia K-12 | 400 genomas | `data/complexes/transfer_<genoma>.tsv` | ortología, nivel B |

No hay, hoy, un catálogo curado de complejos para cada una de las 400 especies.

## ¿Se pueden añadir otras especies?

Sí, pero no es descargar un TSV y ya:

- **ComplexPortal** tiene entradas para muchas especies. Para cada una hace falta:
  - `complextab/<taxid>.tsv`;
  - mapeo UniProt → locus_tag del genoma;
  - decidir si se usa como anotación o como validación (circularidad).
- **EcoCyc/MetaCyc** son de *E. coli* y tienen licencia BioCyc; no se
  redistribuyen.
- **KEGG MODULE** son módulos funcionales, no complejos físicos; no sirven para
  “quién se une con quién”.
- **CORUM** es de mamíferos.
- **PDB/IntAct/BioGRID** aportan interacciones/estructuras, no un catálogo
  completo por especie.

## Estrategia actual

1. K-12 como catálogo curado.
2. Transferencia por ortología K-12 → 400.
3. Modelo genómico K-12 como prior.
4. Evidencia de masa, perfil filogenético y estructura.
5. Todo lo transferido se marca `documentado_inferido` (`ECO:0005544`), no
   experimental.

## Qué haría para ampliar

1. Listar taxids de las 400 con entrada en ComplexPortal.
2. Descargar `complextab/<taxid>.tsv`.
3. Mapear UniProt a locus_tag usando `functional.tsv`/GenBank.
4. Guardar `gold_<taxid>.tsv` y usarlo como **anotación**, no como validación.
5. Medir cobertura; no todas las especies tendrán catálogo.

Mientras no exista ese mapeo, la transferencia desde K-12 es la mejor
aproximación disponible, y siempre es una sospecha.
