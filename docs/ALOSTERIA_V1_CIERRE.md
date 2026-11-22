# Cierre de la línea de alosteria — versión 1.0

> Decisión: **la línea de energía libre y magnitud del efector se cierra para v1.0.**
> Motivo: requiere recursos de cómputo (GPU) que no están disponibles en este entorno.
>
> Lo que sí entra en v1.0 es el sistema que **encuentra y anota represores** en genomas
> bacterianos desconocidos, con preparación estructural local y puntuación de contacto.
>
> Fecha: 2026-11-22.

---

## 0 · Objetivo de la línea (palabras textuales del proyecto)

> *"Poder llegar a decir: este represor tiene esta secuencia de afinidad con esta fuerza,
> y este otro factor se une con este efecto para que se una con más fuerza al ADN o lo
> suelte."*

O sea: un motor que, dado un represor bacteriano, entregue **cuatro cosas**:

| # | Qué debe decir | En una frase |
|---|---|---|
| 1 | a qué secuencia de ADN se une | "reconoce este operador" |
| 2 | con qué fuerza | "con esta afinidad" |
| 3 | qué molécula lo regula | "responde a este metabolito" |
| 4 | **en qué dirección y cuánto** | "lo aprieta / lo suelta con este efecto" |

### Criterios de éxito que se fijaron

| Criterio | Medida | Estado |
|---|---|---|
| Recuperar reguladores conocidos en un genoma conocido | recall contra RegulonDB | **cumplido** (92 %) |
| Predecir sitios mejor que el azar | percentil frente al genoma | **cumplido** (98,7 %) |
| Predecir la **dirección** del efector | signo correcto en ≥8/10 casos | **no cumplido** (1/6, 2/6) |
| Predecir la **magnitud** | ΔΔG con error ≤2 kcal/mol | **no abordado** (sin GPU ni banco) |

### Lo que esta línea NO era

- No era diseño de novo de represores (RFdiffusion/ProteinMPNN).
- No era predecir plegamiento (se parte de estructura existente o AFDB).
- No era sustituir la experimentación: era priorizar qué medir.

### Por qué se cierra

Los criterios 1 y 2 se cumplieron y son lo que entra en v1.0. Los criterios 3 y 4
dependen de cómputo (GPU) o de banco, y ninguno de los dos está disponible. La línea se
cierra **con los dos primeros cumplidos y los dos segundos con criterio de reapertura
definido**, no por falta de ideas.

---

## 1 · Qué entra en v1.0

| Pieza | Qué hace | Medido |
|---|---|---|
| **Pipeline genoma → represores** (`18_genome_pipeline.py`) | genes → dominios → familia → homólogo DIAMOND → efector heredado → sitios PWM | **recall 92 %** en *E. coli* (63 % de novo) · ~22 s/genoma de 2 Mb · 130 MB RAM |
| **Campaña** (`22_campaign.py`) | corre el pipeline sobre lotes y agrega | 85 represores en 5 genomas · 22,2 s/genoma |
| **Benchmark curado** | 43 represores con efector, 36 decisivos, PMID y evidencia | `data/allostery/` |
| **Evaluador GBA** (`00_c0_gba.py`) | ¿el efector está en la química del regulón? | **80 % vs 5,3 %** a distancia 0 (**15×**) |
| **Detector de bisagra** (`06`-`08`) | modo DBD-dominado y antiparalelo | sensibilidad 1/1 · **especificidad 5/8** |
| **Preparación local** (`25_prepare_local.py`) | protonar pH 7 + minimizar | 10,9 s prep + 1,9 s/paso |
| **Puntuación de contacto** (`26_contact_score.py`) | pose real vs decoy en el bolsillo | **separa** (PurR 1,50 vs 0,00; LacI 1,37 vs 0,00) |
| **GUI** (`app_alloreg.py`) | genoma → represores · preparar estructura | 2 páginas |
| **Documentación** | métodos, límites y decisiones | `PIPELINE_GENOMA.md`, `PROTEINAS_REGULADORAS.md`, `PLAN_ALOSTERIA.md`, `MODAL.md` |

## 2 · Qué queda fuera de v1.0, y por qué

| Pieza | Estado | Motivo |
|---|---|---|
| MM/GBSA del complejo | **fuera** | requiere MD con solvente explícito: 10 ns = meses en local |
| ABFE / ΔΔG del efector | **fuera** | 100-150 GPU-h por sistema |
| Magnitud del efecto del efector | **fuera** | depende de lo anterior |
| Energía MM local | **parcial** | funciona, pero en vacío: Coulomb sin apantallar, no ordena afinidades |

Resumen de los cuatro métodos ya probados para "¿inductor o correpresor?" — todos
negativos o inválidos con los datos y el cómputo disponibles:

| Método | Resultado |
|---|---|
| Secuencia (familia, ESM-2) | 0,42 / 0,59 balanced acc; p=0,24 |
| Docking ciego (13 casos) | top-1 0,31 vs 0,25 azar |
| Co-folding Boltz-2 (2 diseños) | 1/6 y 2/6 signos correctos |
| Gating de cavidad | control positivo inválido |

## 3 · Qué haría falta para reabrir la línea

| Requisito | Detalle | Coste estimado |
|---|---|---|
| **Presupuesto GPU** | MM/GBSA de 3 controles en NeuroSnap | **$1-9** |
| | ABFE de 1 caso, solo si lo anterior ordena | **$7-82** |
| **o** presupuesto Modal | MD 10 ns × 3 + MM/GBSA local | **$51-135** |
| **o** banco | EMSA/ITC con y sin ligando en 12-20 casos | — |

El criterio de reapertura es explícito: **primero MM/GBSA de LacI, PurR y TrpR**; si no
ordena las direcciones conocidas, no se gasta en ABFE y la línea queda cerrada con
evidencia, no por falta de intentos.

## 4 · Cómo se retoma

Todo está en el repo y es reproducible:

```bash
# pipeline
python scripts/allostery/18_genome_pipeline.py genoma.fna --out salida --find-repressors

# campaña
python scripts/allostery/22_campaign.py --n 5

# preparación local y contacto
python scripts/allostery/25_prepare_local.py receptor.pdb --max-iter 100
python scripts/allostery/26_contact_score.py receptor_preparado.pdb ligando.sdf --decoy

# GUI
streamlit run scripts/allostery/app_alloreg.py
```

El entorno de energía MM (Python 3.12 + OpenFF, **3,8 GB**) está en
`.mamba_root/envs/.mamba_mm312`. Si se borra para liberar disco, se recrea con lo
documentado en `PIPELINE_GENOMA.md` §5.8 (~15 min, $0).

## 5 · Lo que se aprendió

Lo biológico y lo metodológico está en [`PROTEINAS_REGULADORAS.md`](PROTEINAS_REGULADORAS.md):
siete hechos establecidos, cinco métodos con sus números, ocho errores que costaron
tiempo y cuatro preguntas abiertas.

La frase que resume v1.0:

> El sistema puede decir **qué regulador** actúa en **qué región**, con **qué afinidad
> relativa** y **qué efector es probable por familia u homología**. No puede decir
> **cuánto** cambia el efector la unión. Esa segunda mitad es un problema de datos
> experimentales y de cómputo, no de código, y así queda registrado.
