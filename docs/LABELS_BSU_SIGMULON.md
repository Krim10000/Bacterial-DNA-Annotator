# Símbulón derivado de B. subtilis

No existía una tabla σ→TU explícita como `NetworkSigmaTU` de RegulonDB.
Se derivó de:

```text
data/sigma/gold/subtilis/tss_sigma.tsv
data/dbtbs/gt_operons.tsv
```

Método:

```text
TSS gold + strand
→ operón BSGatlas que lo contiene o empieza ≤1.500 pb aguas abajo
→ asignación σ→TU
```

Resultado:

```text
TSS:        2.226
Mapeados:   2.167
Sin mapear:    59
Operones:   1.560
```

Por sigma:

| Sigma | TSS |
|---|---:|
| sigma70 | 1.900 |
| sigma38 | 177 |
| sigma28 | 68 |
| sigma32 | 11 |
| sigma54 | 11 |

Estado:

```text
derived_bsu_sigmulon
no independiente
```

Usa la misma anotación BSGatlas/DBTBS que el gold. No es una validación
experimental independiente; sirve para perfiles y evaluación descriptiva.

Salida:

```text
reports/labels_pilot/bsu_derived_sigmulon.tsv
reports/labels_pilot/bsu_derived_sigmulon.json
```

Script:

```bash
python scripts/labels/11_build_bsu_sigmulon.py
```
