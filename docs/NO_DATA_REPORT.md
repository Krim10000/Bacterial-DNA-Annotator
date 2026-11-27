# Reporte de no_data

## Total

```text
no_data: 859.533
```

## Por contexto

```text
con contexto de operón:          205.801  (23,9 %)
  en medio del operón:            83.793
  en borde:                      122.008
sin contexto:                    653.732
  sin archivo de operones:       152.719
  operón sin vecino anotado:     476.531
  singleton/sin contexto:         24.482
```

## Qué significa cada bloque

### 205.801 con contexto

Tienen vecinos anotados. Se pueden usar como:

```text
IGC bajo
priorización
feature del clasificador
```

No son GO directo.

### 476.531 en operón sin vecino anotado

El contexto local no ayuda. Rutas:

```text
clusters/ortología
contexto conservado entre genomas
dominios
estructura
no_data
```

### 152.719 sin archivo de operones

Falta generar operones para esos genomas. Es un problema de pipeline, no de biología.

### 24.482 singleton

Sin contexto. Dependen de dominio/estructura/ortología.

## Cobertura actual

```text
genes con evidencia: 665.386
no_data:              859.533
total:              1.524.919 filas de genes
```

## Clasificador

```text
full:    557.763
partial: 107.623
no_data: 859.532
```

## Siguiente

```text
1. Transferencia de contexto desde clusters para los 476.531.
2. Generar operones faltantes.
3. Dominios por lotes en entorno con más RAM.
4. Estructura/Foldseek para el residuo.
5. Reporte por especie y por filo.
```

## Contexto conservado por cluster

```text
con contexto por cluster:  64.330
sin cluster:              177.319
cluster sin donante:      617.883
```

Archivo:

```text
data/labels/context/gene_cluster_context.tsv
```
