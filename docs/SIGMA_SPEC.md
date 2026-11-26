# Especificación canónica de sigma

**Estado:** normativo. Este documento manda sobre los demás `SIGMA_*.md`.
Fecha: 2026-11-26.

## 1. Objetivo

Dado un genoma bacteriano:

1. detectar qué factores sigma codifica;
2. separar los que tienen ortólogo/familia conocida de los linaje-específicos;
3. proponer, con evidencia, qué promotores/regiones podrían depender de cada uno;
4. asociar hipótesis funcionales a rutas/procesos;
5. abstenerse cuando no hay evidencia.

No es objetivo:

- asignar un sigma a cada TSS;
- forzar el número de sigma de la literatura;
- confirmar función sin experimento;
- transferir un modelo de una especie a otra sin auditoría.

## 2. Definiciones

```text
sigma candidato
  proteína detectada por dominio/ortología, completa o parcial.

familia
  sigma70, sigma54, ECF, esporulación, etc.

orthogroup
  grupo de ortología; unidad para transferir función.

native_name
  nombre local (SigA, SigB, RpoD, RpoS...). No es una equivalencia.

regulón
  genes cuya expresión depende de un sigma, con evidencia.

promotor dependiente de sigma
  región con evidencia de reconocimiento por un sigma.

TSS
  inicio de transcripción. No es sinónimo de promotor.

hipótesis contextual
  asociación inferida por vecindad/regulón/ruta. No es confirmación.
```

## 3. Regla de oro

Nunca:

```text
SigB = sigma38
sigmaE = sigmaE
```

Siempre:

```text
SigB de B. subtilis → orthogroup rpoS-like
sigmaE de E. coli   → orthogroup rpoE/ECF
sigmaE de B. subtilis → orthogroup sigE/esporulación
```

## 4. Arquitectura por capas

### Capa 0 — Inventario

Entrada: proteoma.

```text
HMMER/Pfam + arquitectura de dominios + ortología
```

Dominios sigma válidos:

```text
PF00140  Sigma70_r1_2
PF04542  Sigma70_r2
PF04539  Sigma70_r3
PF04545  Sigma70_r4
PF00309  Sigma54
```

**PF04997 no es sigma**: es `RNA_pol_Rpb1_1`. No usarlo.

Salida:

```text
sigma_id, locus, protein, length, domains, family, orthogroup,
copy_number, status, evidence, confidence
```

Estados:

```text
complete
partial
pseudogene
lineage_specific
uncharacterized
```

### Capa 1 — Identidad y ortología

- DIAMOND/HMMER;
- OrthoFinder/árbol;
- orthogroup estable;
- native_name local;
- literatura por orthogroup.

Nunca se transfiere función por nombre.

### Capa 2 — Evidencia de promotor

- PWM por familia/orthogroup;
- arquitectura: −35/−10, extended −10, −24/−12, ECF;
- fondo nulo específico de especie;
- ventanas de operón si no hay TSS;
- `model_source` y `window_source`.

Salida:

```text
score, z, pvalue, architecture, model_source, status=candidate
```

No `confirmed`.

### Capa 3 — Contexto

Features permitidas:

- regulón curado;
- operón/vecindad;
- rutas curadas;
- co-expresión;
- conservación;
- sigma/anti-sigma.

Toda feature debe declarar:

```text
source, independence, weight, ablation
```

Prohibido:

- usar la misma anotación como feature y gold;
- convertir un gen vecino en prueba;
- llamar función a una asociación de ruta.

### Capa 4 — Estructura

La estructura sirve para:

- familia/dominio;
- integridad;
- ECF;
- anti-sigma/interfaz.

No sirve para decidir el promotor.

### Capa 5 — Calibración

Solo con gold de la especie o de un clado validado.

Salida:

```text
P(sigma | evidencia)
```

con calibración, cobertura y abstención.

## 5. Salida por sigma

```text
sigma_id
locus_tag
native_name
family
orthogroup
domains
status
literature_status: known | family_only | undescribed
hypothesized_role
pathways
regulon_evidence
promoter_evidence
context_evidence
structure_evidence
calibration_source
confidence: A|B|C|D|E
rationale
```

## 6. Niveles de confianza

```text
A  ortholog + literature + functional evidence
B  ortholog + domains + context
C  known family, unknown function
D  lineage-specific, hypothesis only
E  partial/pseudogene
```

Ningún nivel A equivale a confirmación experimental de un promotor.

## 7. Dossier por especie

Ejemplo:

```text
Sigma detectados: 10
  completos: 8
  parciales: 2
  ortólogo conocido: 7
  linaje-específicos: 2
  sin literatura: 6
  con hipótesis contextual: 4
```

Tabla:

| ID | Locus | Nombre | Familia | Orthogroup | Literatura | Contexto | Hipótesis | Confianza |
|---|---|---|---|---|---|---|---|---|

El dossier debe explicar qué falta y por qué.

## 8. Escala a 400

### Tier A

Inventario por dominios/ortología para todas.

### Tier B

Evidencia de promotor por familia para todas, con abstención.

### Tier C

Calibración solo para especies con gold.

Nunca reportar accuracy global en especies sin gold.

## 9. Validación

- GroupKFold por locus/operón;
- leave-one-species-out;
- leave-one-clade-out;
- negativos reales;
- macro-F1, AUPRC, cobertura, calibración;
- no top-1 aislado.

## 10. Estado actual

Implementado:

```text
inventario por anotación
políticas argmax/two-stage/hybrid
calibrador k-mer B. subtilis
inventario masivo 442 genomas
```

Experimental:

```text
modelo jerárquico K-12
contexto sigma
k-mer K-12
```

No implementado:

```text
HMMER/ortología para 400
PWM por familia
dossier por especie
LOSO/LOCO
```

## 11. Documentos fuente

Canónicos:

```text
SIGMA_SPEC.md
SIGMA_STRATEGIES.md
SIGMA_INVENTORY.md
SIGMA_KMER_CALIBRATOR.md
SIGMA_GOLD.md
SIGMA_BENCHMARK.md
SIGMA_400_PLAN.md
SIGMA_META_CLASSIFICATION.md
LABELS_MODULE.md
literature/README.md
```

Históricos o parciales: no usar como fuente de verdad sin contrastar.

## 12. Control de cambios

Toda afirmación nueva debe indicar:

```text
fuente
dato
método
límite
```

Si no hay fuente, se marca `hipótesis`.
