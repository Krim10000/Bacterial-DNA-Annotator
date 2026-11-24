# Genes solapados y marcos de lectura

## Estado actual

`run_cds()` usa `pyrodigal.GeneFinder(meta=True)` y conserva **todas** las
predicciones de ambas hebras.

No hay:
- resolución automática de solapamientos;
- elección de un gen ganador;
- fase/marco calculado por defecto.

Por eso dos CDS pueden compartir coordenadas y ambos aparecer en el GFF,
proteínas, función y metabolismo.

## Navegador

- sección **HEBRA +** con todos sus elementos;
- sección **HEBRA −** con todos sus elementos;
- sección **SIN HEBRA** para complejos sin hebra única, TSS/TTS sin dato y sospechas;
- flechas en sentido de lectura;
- complejos: la hebra se deriva de sus subunidades; si hay subunidades en ambas hebras, quedan sin hebra;
- horquillas: usan su `strand` cuando existe; las `intergénica` quedan sin hebra;
- tabla “Solapamientos en pantalla”;
- `bp_solape`, hebras, `mismo_marco` y `locus_compartido`;
- no se elimina ni decide cuál es correcto.

## Límites

- el marco solo se muestra si viene en `frame`/`phase`; si no, `unknown`;
- un solapamiento no implica error: puede ser un gen alternativo, un
  pseudogén, un ORF embebido o un artefacto;
- la decisión requiere evidencia adicional.

## Futuro

- carril de marcos por CDS;
- fase GFF real;
- prioridad configurable;
- alerta de locus duplicado.
