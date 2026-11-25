# Bacillus subtilis como genoma de referencia

`subtilis` se añade como **genoma de referencia curado**, separado de los 400
genomas libres.

## Archivos

```text
data/raw/subtilis_genomic.fna.gz
data/annotations/subtilis.gff.gz
data/regulation/subtilis/genes.tsv
data/operon_finder/subtilis_hairpins.tsv
data/dbtbs/subtilis_tss_tts.tsv
data/taxonomy/full/subtilis.json
```

## Pipeline ejecutado

```text
bp: 4.215.606
GC: 43,5 %
cobertura anotada: 88,5 %
CDS: 4.242
tRNA: 86
rRNA: 20
ncRNA: 6
funcional: 4.134 hits (97,45 %)
operones: 1.334
promotores: 1.334
complejos: 841 candidatos
TSS: 3.390
TTS: 3.009
oriC: 3.237 pb (consenso ORCA + dnaA)
```

Salidas:

```text
results/subtilis/
```

## Límites

- No forma parte de las métricas de los 400.
- Taxonomía curada, no ANI.
- Complejos: candidatos; no estructura cuaternaria.
- TSS/TTS: precalculados BSGatlas; no medidos en este ensamblado.
- No hay catálogo ComplexPortal independiente.
