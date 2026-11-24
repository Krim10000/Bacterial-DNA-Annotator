# Criterio de cierre e interpretación de brechas (frente CPU)

Este documento define **cuándo una brecha metabólica puede considerarse
asociada, candidata o no cerrable**, y cuándo —nunca— puede considerarse
cerrada. Acompaña a `tools/audit_gap_closure.py`, que audita de solo lectura
`data/metabolism/gap_interpretation.tsv` y `data/complexes/transfer_*` /
`model_*` sin GPU y sin red.

## Alcance

- Frente **CPU/offline**: no usa GPU, AlphaFold-Multimer, Boltz, Foldseek ni
  APIs externas.
- **No modifica** la red metabólica, `gap_interpretation.tsv` ni las tablas de
  complejos. Solo lee y escribe `reports/gap_closure_audit.{md,json}`.
- **No declara cierres.** El estado `confirmed` no existe en esta taxonomía.
- El resultado es una **interpretación con evidencia y límites**, no una
  validación experimental.

## Regla central: asociada ≠ cerrada

Una brecha puede estar **asociada** a una enzima homóloga, a un dominio, a un
complejo transferido o a una reacción de transporte y seguir **abierta**. La
asociación solo cambia la explicación y la confianza.

Un cierre real requeriría, según el caso:

- actividad enzimática medida para el sustrato/producto correcto;
- complementación de un mutante o reconstitución in vitro;
- ensayo de transporte/importación;
- estructura + bioquímica que demuestre la interfaz y la estequiometría;
- evidencia genética equivalente en el organismo.

Ninguna de estas cosas se obtiene de las tablas actuales. Por eso:

```text
homología fuerte      -> evidencia de candidato, no cierre
dominio HMM           -> hipótesis de familia, no cierre
complejo transferido  -> sospecha de ensamblaje, no cierre
transportador         -> hipótesis de importación, no cierre
modelo estructural    -> hipótesis geométrica, no cierre
```

La regla también aplica a GPU: un modelo de AlphaFold-Multimer o un ipTM alto
no demuestra que el complejo exista ni que cierre la brecha.

## Estados

Cada brecha recibe **un** estado, por precedencia. Los estados son mutuamente
excluyentes; las señales de evidencia, complejo, transporte y dirección se
reportan aparte como cobertura.

| estado | significado | ¿cierra? |
|---|---|---|
| `gap_evidence_supported` | hay evidencia directa de un cierre candidato: mapeo Rhea validado, homólogo fuerte a un productor o dominio confirmado/reanotado | **no** |
| `gap_candidate` | hay evidencia parcial: dominio HMM, homólogo débil, revisión de soporte, multifunción misma EC o asociación débil a complejo | **no** |
| `complex_suspected` | hay asociación a un complejo transferido fuerte (`frac_mapped >= 0.8`) no declarado incompleto | **no** |
| `complex_incomplete` | no hay asociación fuerte no incompleta y sí al menos una asociación a un complejo incompleto (faltan subunidades o estequiometría) | **no** |
| `not_evaluated` | no hay evidencia, complejo, transporte ni bloqueador explícito: no hay hipótesis de cierre que evaluar | **no** |
| `not_closable_cpu` | hay una capa evaluada que impide cerrar con CPU/offline; incluye no-brechas y brechas reales que requieren experimento, ΔG, base externa o GPU | **no** |

Precedencia:

```text
evidencia fuerte
  > evidencia débil
  > complejo fuerte no incompleto
  > complejo incompleto
  > asociación débil a complejo
  > bloqueador explícito
  > not_evaluated
```

`not_closable_cpu` agrupa dos cosas que el reporte separa:

- **no es una brecha biosintética** (catabolismo, falso positivo por
  dirección/moneda, punto de entrada, reacción espontánea);
- **brecha real sin cierre CPU/offline** (dirección no resuelta, transporte sin
  homólogo, sin homólogo, química novel, dependencia del hospedador).

Agruparlas en un solo estado es deliberado: en ambos casos **no hay cierre
biosintético que emitir** con esta evidencia. El motivo primario y el grupo
quedan explícitos en el JSON y el Markdown.

La evidencia directa se evalúa antes que los bloqueadores. Si una fila marcada
como catabólica o `false_positive` tiene homología fuerte o mapeo validado, se
conserva como evidencia: la etiqueta de rol/auditoría es revisable y no debe
descartar una señal positiva. En la instantánea actual no hay filas fuertes en
esos grupos; solo unas pocas débiles.

## Evidencia directa

Se calcula solo con columnas de `gap_interpretation.tsv`:

| nivel | criterio |
|---|---|
| fuerte | `mapping_validated == 1`; o `homolog_status == homolog_strong` con `homolog_acc`; o `review_status` en {`reannotated_hmm`, `domain_confirmed`} con `hmm_evidence` |
| débil | `hmm_evidence`; o `homolog_status` en {`homolog_weak`, `homolog_low`} con `homolog_acc`; o `status == candidate`; o `review_status` en el conjunto de soporte parcial; o `origin == promiscuous` (misma EC); o `mapping_rhea` presente |
| ninguna | no cumple ninguno de los anteriores |

Conjunto de revisión parcial (`WEAK_REVIEW_STATUSES` en el código):

```text
evidence_high
mapping_context_supported
mapping_rhea_and_ec_supported
mapping_rhea_supported_ec_unknown
weak_rhea_supported
high_priority_homology
high_priority_low
evidence_high_no_rhea
mapping_context_no_producer
```

No se usan como evidencia `keep_missing`, `no_domain_data` ni
`one_sided_not_transport`: son ausencia de evidencia o dudas, no cierres.

**Límite:** estas etiquetas siguen siendo inferencia. `homologo_fuerte` no
prueba que el producto sea el correcto en este genoma, y un dominio HMM puede
tener funciones amplias. Por eso el estado se llama `supported`, nunca
`confirmed`.

## Complejos

La asociación usa el mismo puente que el resto del proyecto:

```text
EC del complejo -> rhea2ec -> reacción Rhea -> productos (ChEBI)
```

Una brecha se asocia a un complejo si:

1. el genoma del complejo es el mismo que el de la brecha;
2. el ChEBI de la brecha está **exactamente** entre los productos de alguna EC
   del complejo;
3. el complejo tiene EC/Rhea utilizable.

No se asocia por nombre de complejo, texto libre, vecindad génica ni similitud
de ligandos.

### Fuerza y completitud

| fuente | fuerte | incompleto |
|---|---|---|
| `transfer_<genoma>.tsv` | `frac_mapped >= 0.8` | `stoich_status == incomplete`, o `observed_n < expected_n` |
| `model_<genoma>.complexes.tsv` | no aplica (`complex_candidate`) | `missing` no vacío |

`frac_mapped` mide ortología desde K-12, no interacción física. Una
estequiometría `unknown` **no** es `complete`. El modelo K-12 es un prior
genómico, no una calibración para otra especie.

### Estados de complejo

- `complex_suspected`: hay al menos una asociación fuerte no incompleta.
- `complex_incomplete`: no hay asociación fuerte no incompleta y sí una
  asociación incompleta.
- `gap_candidate`: solo hay asociaciones débiles (`frac_mapped < 0.8`) no
  incompletas; es una hipótesis débil, no un cierre.
- Si una brecha tiene a la vez un complejo fuerte y otro incompleto, prevalece
  `complex_suspected`, pero el reporte conserva la señal de incompletitud.

## Transporte y dirección

Se reportan como **señales**, no como cierres:

| campo | valores | lectura |
|---|---|---|
| `transport_status` | `transport_reaction_only` | hay reacción de transporte en Rhea, pero sin homólogo de transportador; importación no confirmada |
| `transport_status` | `t1_bacterial` | hay un transportador homologado; la importación sigue siendo una hipótesis |
| `transport_status` | `no_transport` | se evaluó y no se encontró transporte; no es una prueba de ausencia |
| `direction_status` | `thermo_forward`, `thermo_reversible` | dirección resuelta (señal positiva) |
| `direction_status` | `direction_unspecified`, `thermo_unknown` | dirección no resuelta; Rhea master `UN` deja los lados arbitrarios |
| `direction_status` | vacío | no evaluado |

`thermo_status` está vacío en la instantánea actual. No se calcula ΔG ni se
ejecuta Equilibrator en esta auditoría. Por eso una dirección no resuelta es un
bloqueador: sin ΔG o experimento no se puede afirmar que la reacción **produce**
el metabolito de la brecha.

## Bloqueadores (`not_closable_cpu`)

| código | significado |
|---|---|
| `catabolic` | catabolismo/degradación; no requiere productor biosintético |
| `false_positive_direction` | falso positivo por dirección/moneda; la reacción consume el metabolito |
| `entry_point` | punto de entrada/importado; no es una brecha biosintética |
| `spontaneous` | posible reacción espontánea/no enzimática; falta base curada |
| `host_dependent` | endosimbionte/endoparásito; el hospedador puede proveer el metabolito |
| `no_homolog` | sin homólogo detectado; hace falta base mayor o experimento |
| `novel_no_reference` | química sin referencia en los 400; requiere caracterización experimental |
| `transport_reaction_only` | reacción de transporte sin homólogo; importación no confirmada |
| `transporter_homolog` | transportador presente; importación no medida |
| `direction_unresolved` | dirección no resuelta; requiere ΔG/experimento |

Una fila puede tener varios bloqueadores; el reporte guarda el primario y el
conteo por todos los motivos.

## Resultado de la auditoría (instantánea 2026-11-24)

Entrada: 8.655 brechas en 309 genomas. Puente local disponible: 75.208
reacciones Rhea y 6.161 EC con productos.

### Estados

| estado | brechas | % |
|---|---:|---:|
| `gap_evidence_supported` | 1.022 | 11,8 % |
| `gap_candidate` | 1.147 | 13,3 % |
| `complex_suspected` | 14 | 0,2 % |
| `complex_incomplete` | 0 | 0,0 % |
| `not_evaluated` | 7 | 0,1 % |
| `not_closable_cpu` | 6.465 | 74,7 % |
| **total** | **8.655** | 100 % |

### Cobertura de señales

| señal | brechas |
|---|---:|
| evidencia directa (fuerte o débil) | 2.161 |
| evidencia fuerte | 1.022 |
| evidencia débil | 1.139 |
| asociación a complejo | 22 |
| señal de transporte | 907 |
| dirección resuelta | 439 |
| dirección no resuelta | 4.379 |
| hipótesis de multifunción (misma EC) | 522 |
| **sin señal positiva de ningún tipo** | **5.830** |
| sin ninguna capa de interpretación | 0 |

### Complejos

| insumo | archivos | filas | con EC | con productos | fuertes | incompletos |
|---|---:|---:|---:|---:|---:|---:|
| `transfer_*` | 400 | 19.085 | 12.996 | 11.361 | 8.408 | 1.785 |
| `model_*` | 400 | 281.090 | 0 | 0 | — | 0 |

Asociaciones encontradas: **22 brechas**, **24 pares brecha–complejo**.

- 14 brechas con complejo fuerte (`frac_mapped >= 0.8`);
- 2 brechas con complejo incompleto;
- 2 brechas con complejo fuerte **y** además incompleto;
- 2 brechas con solapamiento de IDs Rhea (señal de apoyo, no estado distinto);
- 0 asociaciones desde `model_*`, porque esas 281.090 filas no traen EC.

Las 14 asociaciones fuertes se concentran en:

- CPX-1092 (succinil-CoA sintetasa) y succinil-CoA: 6 brechas;
- CPX-3446 (TrpAB, triptófano sintasa) e indol: 6 brechas;
- CPX-4443 (transportador ABC de Na/Li, EC 7.6.2.2) y daunorrubicina: 2
  brechas. La EC es de transporte y la asociación es solo por intersección de
  producto; debe revisarse antes de interpretarla como producción.
- Los transportadores ABC de poliaminas CPX-4383/CPX-4384 (EC 7.6.2.11) dan
  asociaciones débiles (`frac_mapped` 0,5–0,75) y quedan como `gap_candidate`,
  no como sospecha fuerte. Asociar una EC de transporte no convierte una
  importación en una biosíntesis.

### Bloqueadores

| motivo primario | brechas |
|---|---:|
| `direction_unresolved` | 2.574 |
| `catabolic` | 2.390 |
| `false_positive_direction` | 603 |
| `transport_reaction_only` | 600 |
| `entry_point` | 246 |
| `novel_no_reference` | 43 |
| `spontaneous` | 9 |
| **total `not_closable_cpu`** | **6.465** |

De las 6.465: **3.248** no son brechas biosintéticas (catabolismo, falso
positivo, entrada, espontánea) y **3.217** son brechas reales sin cierre
CPU/offline (dirección no resuelta, transporte sin homólogo, novel).

### Hallazgos

1. **No hay ninguna brecha cerrada ni confirmada.** El reporte emite
   `confirmed = 0` y `closed = 0`; las 22 asociaciones quedan abiertas.
2. **La mayoría de las brechas no tienen señal positiva** (5.830/8.655). No es
   lo mismo que `not_evaluated`: muchas tienen un resultado negativo
   (`no_transport`) o un bloqueador.
3. **`model_*` no puede asociar brechas** en esta instantánea: 281.090
   candidatos, 0 con EC. Es un límite del modelo, no evidencia de que los
   complejos no existan.
4. **La transferencia K-12 solo asocia 22 brechas** (0,25 %). La estructura
   cuaternaria no parece la explicación mayoritaria de las brechas con estos
   datos.
5. **7 filas quedan `not_evaluated`**: todas son `propanoyl-CoA` con un hit de
   transportador unilateral (`review_status = one_sided_not_transport`), sin
   evidencia de cierre ni bloqueador concluyente.
6. **La dirección no resuelta es el bloqueador individual más frecuente**
   (2.574 como motivo primario; 3.173 entre las `not_closable_cpu`). Sin
   ΔG/experimento, la mitad “productora” de Rhea master no se puede afirmar.

## Cómo ejecutar

```bash
.venv/bin/python tools/audit_gap_closure.py \
    --data-root data \
    --out-md reports/gap_closure_audit.md \
    --out-json reports/gap_closure_audit.json

.venv/bin/python -m pytest tests/test_audit_gap_closure.py -q
```

El JSON incluye definiciones, conteos por estado, señales, estadísticas de
complejos, lista de asociaciones, detalle de `not_evaluated` y conteos por
genoma. El Markdown es el resumen legible.

## Límites

- **Puente limitado a EC/Rhea/ChEBI.** Una subunidad sin EC no puede asociarse
  aunque sea esencial; una EC compartida puede producir metabolitos moneda y
  generar una asociación espuria. Las asociaciones se reportan como tales.
- **`frac_mapped` no es interacción física.** Es ortología recíproca desde
  K-12, con identidad y cobertura incompletas en los TSV actuales.
- **Estequiometría `unknown` no es completa.** El modelo K-12 no conoce la
  composición esperada y sus 281.090 filas no tienen EC.
- **Evidencia directa no es función demostrada.** `mapping_rhea_validada`,
  `homologo_fuerte` y `dominio_confirmado` son candidatos.
- **No se calcula termodinámica.** `direction_unresolved` queda como bloqueador;
  Equilibrator/ΔG sería otra ruta CPU, no ejecutada aquí.
- **No hay GPU ni estructura.** AlphaFold-Multimer, Boltz, co-folding o
  Foldseek-Multimer podrían priorizar hipótesis, pero no cerrar brechas.
- **La instantánea puede cambiar.** Los conteos corresponden a la fecha
  indicada; rerun del script los actualiza.
- **La asociación no modifica la red.** El script no escribe en `data/` ni en
  el pipeline/GUI; cualquier integración posterior es un cambio separado.

## Relación con otros documentos

- [`data/metabolism/GAP_LOGIC.md`](../data/metabolism/GAP_LOGIC.md): estados
  `validated`/`consensus`/`novel`/`anomalous`/`lineage_specific` y soporte 400.
  No son estados de cierre; son interpretación de la brecha.
- [`data/metabolism/GAP_AUDIT.md`](../data/metabolism/GAP_AUDIT.md): auditoría
  de falsos positivos, degradación, transporte y dirección.
- [`COMPLEX_GAPS.md`](COMPLEX_GAPS.md): por qué un complejo no cierra una
  brecha.
- [`COMPLEX_TRANSFER.md`](COMPLEX_TRANSFER.md): transferencia K-12 y umbral
  `frac_mapped >= 0.8`.
- [`COMPLEX_CPU.md`](COMPLEX_CPU.md): estados de completitud de complejos.
- [`CPU_GPU_LIMITS.md`](CPU_GPU_LIMITS.md): qué puede y no puede CPU.
