# Especie nueva: qué corre de verdad

## GUI y CLI

Desde el sprint actual, la GUI usa `webapp/pipeline_adapter.py`, que ejecuta
la misma secuencia de capas que el CLI:

```text
annotate
 → ncRNA
 → taxonomía local
 → functional
 → dominios C
 → contexto C
 → B estructural
 → complejos (firma/modelo/transferencia)
 → regulación
 → membrana
 → oriC
 → metabolismo opcional
 → informe
```

Cada capa devuelve `ok/status/not_available/error`; si falta una dependencia no
se inventan filas.

## Qué significa “especie nueva”

- **Genoma de los 400:** puede tener taxonomía, transferencia K-12 y datos
  precalculados.
- **FASTA realmente nuevo:** no tiene taxonomía local ni transferencia K-12.
  - `taxonomy = not_available`;
  - no se identifica especie por ANI;
  - se puede calcular un **vecino por k-mers** contra los 400;
  - ese vecino es una aproximación, no ANI ni identificación taxonómica.
  - se puede marcar su rama en el árbol como `nuevo/externo`.

## Metabolismo

Solo se ejecuta si se marca la opción y existen Swiss-Prot/Rhea:

```text
proteínas → EC/Rhea → red → brechas
```

No depende de que la especie esté en los 400, pero sí de las bases de datos.

## Complejos

- firmas por vecindad;
- modelo K-12;
- transferencia solo si el genoma coincide con los 400;
- **no** predice estructura cuaternaria ni interfaz;
- el modelo es un prior, no una calibración para la especie nueva.

## Límites

- especie nueva ≠ especie identificada;
- vecino k-mer ≠ ANI;
- metabolismo reconstruido ≠ metabolismo validado;
- complejo candidato ≠ interacción física;
- no hay cierre de brechas.
