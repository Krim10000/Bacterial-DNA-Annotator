# Benchmark por capas

Herramienta: `tools/benchmark_layers.py`
Informe: `reports/layer_benchmark.md` · JSON: `reports/layer_benchmark.json`

Ejecuta:

```bash
python tools/benchmark_layers.py
```

## Resultados

| capa | estado | resultado |
|---|---|---|
| A/B/C | parcial | 400 genomas, 1.337.745 proteínas; A 566.957, B 123.200, C 647.588 |
| Membrana | parcial | 400 genomas, 319.715 filas; 54.951 high, 52.485 medium, 212.279 low |
| Complejos | interno, no ciego | gold K-12 273 complejos; J≥0,8: 86 pares, recall gold 31,5 %, pureza genes 100 % |
| oriC | parcial | 271 validated, 77 consensus, 42 anomalous, 10 lineage_specific; DoriC 306/400 |

## Qué NO mide

- A/B/C no tiene gold local de función/dominios: solo conteos y cobertura.
- Membrana valida formato/coordenadas, no verdad biológica.
- Complejos es un benchmark **interno**: el modelo se ajustó con el mismo gold
  K-12; el solapamiento está favorecido.
- oriC `validated` usa DoriC como ancla; no es exactitud independiente.

## Regla

Un benchmark interno no es una métrica de producción. Para afirmar precisión
hace falta un holdout externo y, para función/estructura, verdad experimental.
