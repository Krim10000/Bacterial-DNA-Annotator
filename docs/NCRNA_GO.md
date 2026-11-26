# GO para ARN (rRNA, tRNA, ncRNA)

## Independiente de proteína

Esta capa no usa accesión UniProt ni ortología. Usa:

```text
GFF de anotación (Infernal/cmsearch, tRNAscan-SE)
results/<genoma>/<genoma>.ncrna.tsv
Rfam (rfam_models/bac_ncRNA.cm, bac_rRNA.cm)
GBFF para K-12
```

## Términos asignados

### rRNA

```text
GO:0005840  ribosome                            CC
GO:0003735  structural constituent of ribosome  MF
GO:0006412  translation                         BP
```

### tRNA

```text
GO:0006412  translation            BP
GO:0006399  tRNA metabolic process BP
```

### ncRNA con Rfam conocido

```text
RF00011 RNase P RNA → GO:0001682 tRNA 5'-leader removal
RF00169 SRP RNA     → GO:0006614 SRP-dependent cotranslational protein targeting
RF00023 tmRNA       → GO:0070929 trans-translation
```

ncRNA sin Rfam conocido **no recibe GO**: `no_data`.

## Salida

```text
data/labels/ncrna_go/<genoma>/ncrna_go.tsv
data/labels/ncrna_go/<genoma>/ncrna_summary.json
```

Columnas:

```text
genome, feature_id, contig, start, end, strand, type, family,
rfam, go_id, go_name, aspect, evidence, source, confidence
```

## Cobertura

```text
402 genomas procesados
26.596 filas GO
23.546 features
```

K-12: 108 features, 130 filas.
B. subtilis: 117 features, 147 filas.

## Fuente

Cada fila declara:

```text
source = GFF:<archivo> | GBFF:<archivo>
evidence = curated_feature_type | Rfam_curated
confidence = high | medium
```

No se presenta como experimento.

## Límites

- Depende de que la anotación GFF tenga rRNA/tRNA/ncRNA tipados.
- No descubre ncRNA nuevos por sí sola.
- Rfam desconocido queda fuera.
- No es prueba de función.

## GUI

Pista:

```text
GO ARN
```

colores:

```text
rRNA  violeta
tRNA  naranja
tmRNA verde
```

Tooltip: tipo, familia, GO y archivo fuente.

## Módulo

```text
modules/ncrna_go/module.py
order = 29
```

Para una especie nueva:

```bash
python scripts/annotation/00_annotate_single.py ...
```

y luego ejecutar el módulo `ncrna_go`.
