# Genomas de referencia

Panel **separado** de los 400 libres. Se usa para benchmark y comparación
sistema vs anotación curada; no entrena el modelo ni cambia las métricas de los
400.

## Manifiesto

`data/reference_genomes.json`

Campos:

```text
id, species, strain, level, resources, status, genome_path,
annotation_path, note
```

Estado:

- `available`: hay FASTA/anotación y resultados del sistema;
- `planned`: pendiente de descarga.

## Panel

Nivel 1: 12 organismos modelo.
Nivel 2: 30 organismos descritos.

Actualmente disponibles:

```text
subtilis
ecoli_k12 (anotación curada; resultados del sistema pueden faltar)
```

El resto está como `planned`.

## Comparación

`webapp/reference_genomes.py` compara:

```text
sistema (summary.json) vs referencia GFF (CDS/tRNA/rRNA)
```

No evalúa función ni exactitud de coordenadas. El gold hereda los criterios de
la fuente.

## Añadir un genoma

1. Descargar FASTA y anotación.
2. Crear/ajustar `data/regulation/<id>/`.
3. Ejecutar el pipeline.
4. Cambiar `status` a `available` en el manifiesto.
5. Añadir `genome_path` y `annotation_path`.
6. No usar el mismo gold para calibrar y evaluar.
