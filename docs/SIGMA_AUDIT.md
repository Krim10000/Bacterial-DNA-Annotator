# Auditoría: ¿el sistema actual detecta los diferentes tipos de sigma?

> **Fecha:** 2026-11-25 · **Ámbito:** `dna_poc/scripts/regulation/`,
> `data/regulation/`, `data/regulondb/gt/`.
>
> **Respuesta corta:** el sistema **previo** (`12_promoters_sigma.py`) detecta
> **genes** de factor sigma en un inventario, pero **no asigna promotores de
> sigma alternativos**. La capa de promotores solo busca el consenso σ70 y,
> cuando encuentra algo, etiqueta `sigma70` por defecto (`assumed_default`). No
> usa las PWMs de σ24/σ28/σ32/σ38/σ54 que ya están en el repositorio.
>
> **Actualización:** el módulo nuevo [`sigma`](SIGMA_MODULE.md) corrige esta
> capa con escaneo por arquitectura, competición entre sigmas y evidencia
> trazable. Esta auditoría queda como diagnóstico del sistema previo y como
> base de calibración.

---

## 1. Qué existía en el sistema previo

### 1.1 Inventario de genes sigma

`scripts/regulation/12_promoters_sigma.py` construye
`regulators_inventory.tsv` con dos señales:

1. **Nombre/producto** (`SIGMA_PATTERNS`): `rpoD`, `rpoH`, `rpoS`, `rpoN`,
   `fliA`/`rpoF`, `rpoE`, `sigA/B/E/F/H/K`, y "sigma factor" genérico.
2. **Dominio Pfam** (`hmmsearch` contra el panel de regulación), con los modelos
   `PF00140`, `PF04542`, `PF04539`, `PF04545`, `PF04997`.

Salida: `data/regulation/<genoma>/regulators_inventory.tsv`.

#### Cobertura agregada (401 genomas: 400 + K-12)

| Familia | Filas |
|---|---:|
| `sigma_unresolved` | 3.765 |
| `sigma70` | 727 |
| `sigmaE` | 650 |
| `sigma54` | 331 |
| `sigmaB` | 124 |
| `sigmaH` | 101 |
| `sigma28` | 91 |
| `sigmaK` | 74 |
| `sigmaF` | 71 |
| `sigma32` | 34 |
| `sigma38` | 25 |
| **Total** | **5.993** |

**Lectura:** 3.765/5.993 = **62,8 %** queda en `sigma_unresolved`. El inventario
sirve para saber que hay genes parecidos a sigma, no para asignar familia con
confianza. Además, `sigmaE`, `sigmaB`, `sigmaH`, `sigmaK` y `sigmaF` son etiquetas
de nombre, sin modelo de promotor asociado.

#### K-12

`data/regulation/ecoli_k12/regulators_inventory.tsv` tiene 8 filas:

| Gen | Familia | Estado | Comentario |
|---|---|---|---|
| `fliA` | `sigma28` | validated | Dominio + nombre |
| `rpoE` | `sigmaE` | validated | Dominio + nombre |
| `rpoS` | `sigma38` | validated | Dominio + nombre |
| `rpoD` | `sigma70` | validated | Dominio + nombre |
| `rpoN` | `sigma54` | consensus | Sin hit HMM; solo nombre |
| `rpoH` | `sigma32` | validated | Dominio + nombre |
| `fecI` | `sigma_unresolved` | validated | Dominio + nombre, pero sin familia |
| `rseD` | `sigmaE` | novel | **Falso positivo**: es un péptido líder de `rpoE`, no un sigma |

El caso `rseD` muestra el límite de las reglas por texto: "rpoE leader peptide"
contiene `rpoE` y se clasifica como `sigmaE`. El estado `novel`/`low` lo marca,
pero sigue apareciendo en el inventario.

### 1.2 Asignación de sigma a promotores

`12_promoters_sigma.py` escanea una ventana de 150 pb aguas arriba del primer gen
del operón y busca:

- caja −35 `TTGACA` (≤2–3 mismatches);
- caja −10 `TATAAT` (≤2–3 mismatches);
- espaciador 14–21 pb;
- −10 extendida `TG`.

Si encuentra un par compatible, escribe:

```text
sigma = sigma70
sigma_status = assumed_default
```

No compara contra PWMs de σ alternativos. No calcula un score por sigma. No
resuelve competencia entre sigma.

#### Cobertura agregada de promotores (401 genomas)

| Campo | Filas |
|---|---:|
| Total de promotores/operones | 786.805 |
| `sigma=sigma70` | 305.897 |
| `sigma=""` | 480.908 |
| `sigma_status=assumed_default` | 305.897 |
| `sigma_status=unresolved` | 417.456 |
| `sigma_status=""` | 63.452 |

**No hay una sola fila con `sigma24`, `sigma28`, `sigma32`, `sigma38` o
`sigma54`.** El 100 % de las asignaciones es σ70 por defecto.

#### K-12: lo que hay vs lo que debería haber

`data/regulation/ecoli_k12/promoters.tsv`:

- 2.544 operones;
- 1.677 con par −35/−10;
- 867 `unresolved`;
- **todas** las asignaciones con promotor son `sigma70` / `assumed_default`.

Pero `data/regulondb/gt/gt_promoters.tsv` (verdad curada de K-12) tiene:

| Sigma | Promotores curados |
|---|---:|
| `sigma70` | 2.001 |
| `sigma24` | 523 |
| `sigma32` | 326 |
| `sigma38` | 242 |
| `sigma28` | 147 |
| `sigma54` | 101 |
| `sigma19` | 1 |
| sin etiqueta | 716 |
| **Total** | **4.057** |

Es decir: el sistema previo **no recuperaba** 1.340 promotores de sigma
alternativo que están curados en K-12. No es que falle el 100 %; es que ni
siquiera los busca.

### 1.3 PWMs de sigma disponibles pero no usadas

`data/regulondb/gt/promoter_pwms.tsv` contiene 14 matrices:

| Sigma | Cajas | Comentario |
|---|---|---|
| `sigma24` | −10 (7 pb), −35 (7 pb) | ECF |
| `sigma28` | −10 (8 pb), −35 (6 pb) | Flagelar |
| `sigma32` | −10 (6 pb), −35 (7 pb) | Heat shock |
| `sigma38` | −10 (11 pb), −35 (6 pb) | Estacionaria |
| `sigma54` | −10 (5 pb), −35 (6 pb) | En realidad −24/−12; arquitectura distinta |
| `sigma70` | −10 (6 pb), −35 (6 pb) | Housekeeping |

`scripts/regulation/find_operons.py` tiene `load_promoter_pwm(sigma="sigma70")`,
así que **puede cargar** las matrices alternativas. Pero:

- `rule_operons.py` llama `F.load_promoter_pwm()` sin argumento: siempre σ70;
- `--pwm` permite cambiar de organismo, no de sigma;
- `promoter_scan` asume arquitectura −35/−10 con espaciador 14–21 pb; para
  σ54 (−24/−12, espaciador 10–11) esa geometría no aplica.

**Conclusión:** el dato existía; el consumidor previo no lo usaba.

### 1.4 La capa TSS no es un predictor de sigma

`webapp/tss_layer.py` puede leer etiquetas de sigma de TSS curados (RegulonDB,
DBTBS, BSGatlas) y normalizarlas. Eso es **ingesta de anotación**, no predicción:
si el TSS no trae sigma, no se infiere. No sustituye un escáner de promotores.

---

## 2. Tabla de estado por tipo de sigma

| Sigma | ¿Inventario de genes? | ¿Asignación a promotores? | ¿PWM disponible? | Estado real |
|---|---|---|---|---|
| σ70 | Sí (727) | Sí, por defecto | Sí, usada | Funciona, pero sin competencia |
| σ54 | Sí (331, mostly name) | **No** | Sí, no usada | Falta escáner −24/−12 + bEBP/IHF |
| σ32 | Sí (34) | **No** | Sí, no usada | Falta escáner |
| σ38 | Sí (25) | **No** | Sí, no usada | Falta escáner |
| σ28 | Sí (91) | **No** | Sí, no usada | Falta escáner |
| σ24/ECF | **No como familia** | **No** | Sí, no usada | `fecI` cae en `sigma_unresolved`; falta ECF |
| σ19 | **No** | **No** | No | Solo 1 promotor curado; falta |
| σE | Sí (650) | **No** | No | Etiqueta de nombre; sin modelo |
| σB/H/K/F | Sí (name) | **No** | No | Etiquetas de nombre; sin modelo |
| Anti-sigma | Excluido | No aplica | No | Correcto: no es sigma |

---

## 3. Por qué importa

El flujo propuesto depende de saber **qué sigma usa cada promotor**:

- un promotor “débil” para σ70 puede ser fuerte para σ54 o σ32;
- un promotor con spacer 19–20 puede ser MerR, pero también puede ser σ54;
- buscar activadores en [−100,−35] no sirve para σ54, que usa −24/−12 y bEBP;
- la ruta de rescate “sigma alternativa” no se puede ejecutar si el scanner no
  distingue sigma.

Hoy el sistema no puede responder “¿este promotor es σ70, σ54 o σ32?”. Solo
responde “hay un par −35/−10 compatible, asumo σ70”.

---

## 4. Qué implementar

### 4.1 Separar inventario de asignación

- Inventario: genes sigma, con evidencia y familia.
- Asignación: promotor → sigma, con score y margen.

No usar el inventario como si fuera asignación.

### 4.2 Escáner por arquitectura

| Sigma | Arquitectura | Regla |
|---|---|---|
| σ70 | −35/−10 | Espaciador 15–19; −10 extendida; UP |
| σ54 | −24/−12 | Espaciador 10–11; bEBP; IHF |
| σ32 | −35/−10 | Matriz propia; espaciador 13–15 |
| σ38 | −35/−10 | Matriz propia; −10 extendida |
| σ28 | −35/−10 | Matriz propia; espaciador 15–19 |
| σ24/ECF | −35/−10 | Matriz propia; ECF |

No reutilizar el scanner σ70 para σ54.

### 4.3 Score y competencia

Para cada promotor:

```text
score_sigma70
score_sigma54
score_sigma32
score_sigma38
score_sigma28
score_sigma24
```

Asignar el máximo solo si:

- supera un umbral calibrado;
- el margen sobre el segundo es suficiente;
- la arquitectura es compatible.

Si no, `sigma_status=unresolved`, no `sigma70` por defecto.

### 4.4 Calibración

Usar `data/regulondb/gt/gt_promoters.tsv` (4.057 promotores curados) para:

- calibrar umbrales por sigma;
- medir recall/precisión por sigma;
- construir la matriz de confusión σ70/σ54/σ32/σ38/σ28/σ24;
- reportar por clase, no solo global.

### 4.5 Salida

Añadir a `promoters.tsv`:

```text
sigma
sigma_score
sigma_margin
sigma_status      # assigned | ambiguous | unresolved | not_applicable
sigma_evidence
```

---

## 5. Criterios de aceptación

- [ ] Ningún promotor se etiqueta `sigma70` por defecto.
- [ ] σ54 se escanea con −24/−12 y bEBP/IHF, no con −35/−10.
- [ ] σ24/ECF y σ19 están en el inventario o se declaran `no_model`.
- [ ] Cada asignación tiene score, margen y evidencia.
- [ ] K-12 recupera al menos una fracción medida de los 1.340 promotores
      alternativos curados; el número real se reporta, no se estima.
- [ ] La matriz de confusión por sigma está en el reporte.
- [ ] `rseD` no se cuenta como sigmaE.

---

## 6. Implementación en este repositorio

La auditoría se convirtió en un módulo ejecutable:

| Pieza | Ruta |
|---|---|
| Núcleo | `scripts/sigma/sigma_annotator.py` |
| CLI | `scripts/sigma/10_annotate_sigma.py` |
| Módulo GUI | `modules/sigma/module.py` |
| Documentación | [`SIGMA_MODULE.md`](SIGMA_MODULE.md) |
| Agrupación entre especies | [`SIGMA_GROUPS.md`](SIGMA_GROUPS.md) |
| Tests | `tests/test_sigma_annotator.py` |

Uso:

```bash
python scripts/sigma/10_annotate_sigma.py --genome ecoli_k12
python scripts/sigma/10_annotate_sigma.py --all
python -m modules --selftest sigma
```

El módulo aparece en **🧩 Módulos** y muestra una tabla de evidencia con la fila
completa al seleccionarla.

## 7. Respuesta en una frase

**El sistema previo detectaba genes sigma de forma parcial (62,8 % queda
`sigma_unresolved`) y no asigna promotores de sigma alternativos: solo busca
σ70 y lo pone por defecto. Las PWMs de σ24/σ28/σ32/σ38/σ54 ya existen, pero no
se usaban; el módulo nuevo implementa el escáner por arquitectura.**
