# Cobertura GO no directa

## Objetivo

Cubrir genes sin GO directo usando evidencia no directa, sin inventar.

## Capas

```text
1. Homología / ortología      ISS
2. Dominio InterPro/Pfam      IEA
3. Ruta completa              IGC
4. Estructura/Foldseek        ISS baja
5. Contexto/operón            IGC baja
```

Cada fila debe declarar:

```text
source
evidence
with_from
confidence
```

## Estado actual

### Donantes

```text
277.510 proteínas con GO directo
```

### Consultas

```text
1.068.777 proteínas sin GO directo
```

### Base DIAMOND

```text
donors.faa       121 MB
unannotated.faa  389 MB
donors.dmnd      125 MB
```

### Ejecución

El `diamond blastp` completo de 1,07 M consultas contra 277.510 donantes
no terminó en un tiempo razonable y se detuvo antes de agotar memoria.

```text
no se reporta cobertura no directa todavía
```

## Estrategia correcta

```text
1. Partir las consultas por genoma/lote.
2. Ejecutar DIAMOND por lote.
3. Transferir solo con:
   - identidad ≥ 30 %
   - cobertura ≥ 70 %
   - donante con GO directo
4. Marcar ISS y `with_from`.
5. Medir cobertura por lote.
```

No usar el cluster global de 1,35 M: tampoco cupo en 23 GB.

## Scripts

```text
scripts/labels/19_build_orthogroups.py
scripts/labels/20_transfer_homology_go.py
```

## Lo que falta

- ejecución por lotes;
- criterio de ortología 1:1 para ISO;
- dominio InterPro/Pfam;
- ruta completa IGC;
- validación.

## Límites

- La homología no es ortología 1:1.
- No transferir IEA como si fuera experimental.
- Un hit de dominio no demuestra proceso.
