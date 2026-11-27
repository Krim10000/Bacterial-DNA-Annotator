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

## Piloto estratificado

Se construyó una muestra de 20 `no_data` por genoma:

```text
query.faa: 2,6 MB
hmmscan --cut_ga contra Pfam
```

El `hmmscan` no terminó y fue detenido antes de producir el TSV.

```text
no se reporta cobertura por filo
no se escala
```

Causa probable: 8.800 secuencias contra Pfam completo en un proceso que
agotó el tiempo/memoria disponible.

Siguiente intento:

```text
5.000 secuencias
--cpu 2
Pfam-A
medir tiempo y memoria antes de escalar
```

## Intento reducido

```text
10 no_data por genoma (~4.400 secuencias)
hmmscan --cut_ga --cpu 2
```

Resultado: el `domtblout` quedó incompleto y no se generó
`gene_domain_evidence.tsv`.

Conclusión: Pfam completo no es viable en este entorno con este tamaño de
muestra sin un proceso por lotes más pequeño o más memoria.

Siguiente:

```text
500 secuencias por lote
medir RSS/tiempo
o InterProScan externo
```

## Estado: PILOTO COMPLETADO, ESCALADO DIFERIDO

El intento reducido terminó:

```text
muestra: 4.010 no_data (401 genomas, 10/genoma)
hmmscan --cut_ga --cpu 2
tiempo: 1.625 s (~27 min)
RSS: 1,68 GB
dominios: 1.286 filas
genes con GO por dominio: 655 / 4.010 = 16,3 %
Pfam: 278
términos GO: 211
```

El resultado es válido para el piloto, pero el escalado a 476.531 requiere
mucho más tiempo/lotes: a este ritmo son ~50 h de CPU.

No se usa como cobertura global todavía.

### Retomar en otro entorno

Requisitos recomendados:

```text
RAM: 32–64 GB
CPU: 8–16
Disco: >300 GB
```

Comando base:

```bash
# muestra estratificada
python scripts/labels/29_domains_stratified.py \
  --per-genome 10 \
  --threads 8

# o por lotes de 500
```

Antes de escalar:

```text
1. medir RSS de hmmscan
2. confirmar que el domtblout termina
3. validar 100 hits contra Pfam/InterPro
4. recién entonces correr 476.531
```

No se debe consumir ningún resultado de dominio hasta que exista
`gene_domain_evidence.tsv` completo.

### Aviso de mezcla

La corrida vieja de 20/genoma terminó después y reescribió el resumen con
`sampled=8020`, pero reutilizó el `domtblout` del piloto de 10/genoma.
Los hits son los mismos (655 genes), así que **no aporta cobertura nueva**.

Resultado válido: 4.010 muestreados, 655 genes.
El resumen de 8.020 está marcado `mixed_stale_do_not_use`.
