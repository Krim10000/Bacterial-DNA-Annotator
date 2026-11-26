# Piloto de dominios Pfam

## Ejecución

```text
10.000 proteínas no_data
hmmscan contra Pfam
pfam2go / interpro2go
```

## Resultado

```text
proteínas con dominio:        262 / 10.000  (2,6 %)
proteínas con dominio significativo: 212
proteínas con GO por dominio: 133 / 10.000  (1,3 %)
dominios significativos:      136 Pfam
términos GO:                  147
```

## Lectura

El rendimiento fue bajo en esta muestra.

Causa probable:

```text
las primeras 10.000 secuencias de unannotated.faa
no son una muestra estratificada por genoma
```

No justifica lanzar 476.000 sin estratificar.

## Siguiente

1. Muestra estratificada: 20–50 proteínas por genoma.
2. Comparar `hmmscan` normal vs `--cut_ga`.
3. Medir por filo.
4. Si el rendimiento sube, escalar por lotes.
5. Si sigue bajo, priorizar ortología/estructura.

## Regla

```text
dominio Pfam ≠ función confirmada
IEA / confianza media-baja
```
