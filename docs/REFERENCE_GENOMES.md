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

Estado de la descarga (2026-11-25):

```text
42 genomas en el manifiesto
39 descargados (FASTA + GFF)
2 disponibles y anotados en la GUI: subtilis, ecoli_k12
2 anotados por el sistema: mtb_h37rv, cglutamicum
1 planificado: rsphaeroides (el ensamblado no publica FASTA/GFF)
```

Los 38 restantes tienen FASTA/GFF descargados y están pendientes de ejecutar el
pipeline. No se declaran como `available` hasta que exista `results/<id>/`.

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
