# Modelo K-12 de complejos en los 400

`tools/run_complex_model_400.py` aplica el modelo genómico K-12 a los 400
genomas, sin STRING ni sintenia de K-12.

## Resultado

| métrica | valor |
|---|---:|
| genomas | 400 |
| clusters candidatos | **281.090** |
| confianza alta | 106.547 |
| confianza media | 9.126 |
| confianza baja | 165.417 |
| fallos | 0 |

Cada fila:

```text
status = complex_candidate
flag   = yellow
source = modelo K-12 (build_candidates.py)
```

`expected`/`missing` van vacíos: el modelo no conoce la composición esperada.
No es una predicción de estructura cuaternaria ni cierra brechas.

Salidas:

```text
data/complexes/model_<genoma>.complexes.tsv
reports/complex_model_400_summary.tsv/json
```

El navegador suma estos clusters a las firmas/transferencias y los muestra
como amarillo. El pipeline nuevo usa el mismo modelo mediante
`webapp/complex_model.py`.

## Límite

- K-12 es un prior, no una calibración para otra especie.
- Muchos clusters son `nuevo`/baja confianza.
- Es una hipótesis de asociación, no interacción física.
- No se usa para cerrar brechas.
