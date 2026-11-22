# 05_trace_compare.py — comparación de trazas entre genomas (v1.0)

Una tarea: leer los reportes `{genome}_network_trace.json` que ya escribió
`04_network_trace.py` (v3.3 o superior) y separar las brechas por frecuencia
entre genomas. NO re-analiza redes ni Rhea: solo consume los JSON de reporte.

## Cubetas por frecuencia

- 🎯 **Sistemáticas** — presentes en `>= --min-freq` de los genomas (por
  defecto, en TODOS). Aparecen igual en cualquier genoma → artefacto
  estructural probable (clases genéricas de Rhea, cofactores tipo
  NTP/acil-CoA) o biología central que el anotador no cubre. No discriminan
  entre genomas: revisar una vez, no por genoma.
- 🤝 **Compartidas** — en más de un genoma pero no en todos.
- 🔬 **Específicas** — solo en un genoma → candidatas reales de enzima
  faltante (o de anotación perdida en ese genoma).

Además agrega la **tría** de 04 v3.3 (campo `triage`, por voto mayoritario
entre los reportes que mencionan la brecha); con reportes de 04 < v3.3 la
tría sale como `unknown`.

## Uso

```
# 1) generar/actualizar los reportes (batch, el RDF se carga una vez):
python scripts/metabolism/04_network_trace.py

# 2) comparar:
python scripts/metabolism/05_trace_compare.py
# opciones:
python scripts/metabolism/05_trace_compare.py \
    --reports-dir data/metabolism \
    --out data/metabolism/trace_compare.json \
    --min-freq 0.9 --top 15
```

Con un solo reporte no hay comparación (todo queda como específica y se
avisa). Si el directorio mezcla redes sintéticas de prueba con redes reales,
la comparación se contamina: comparar solo genomas del mismo pipeline.

## Salida

- **Consola**: una línea por genoma (reacciones matched, brechas con su tría,
  sinks, fuentes) y el top de cada cubeta de frecuencia.
- **`trace_compare.json`** (por defecto en `--reports-dir`): `per_genome`
  (counts por genoma), `gap_frequency` (chebi, nombre, tría mayoritaria, freq,
  genomes), `systematic_gaps`, `shared_gaps`, `specific_gaps` y el umbral
  usado (`systematic_threshold`).
