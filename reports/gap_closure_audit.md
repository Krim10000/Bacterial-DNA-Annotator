# Auditoría de cierre/interpretación de brechas

Generado por `tools/audit_gap_closure.py`. Frente CPU, sin GPU, sin red y sin modificar redes ni tablas de entrada.

## Regla

- **Asociada no es cerrada.** Una brecha puede tener evidencia, un
  complejo o un transportador asociado y seguir abierta.
- **No se emite `confirmed`.** Un cierre real exige experimento
  (actividad enzimática, complementación, transporte medido) o una
  validación equivalente.
- Cierres declarados en esta auditoría: **0**.
- Brechas asociadas a algo pero no cerradas: **22**.

## Instantánea

- Brechas auditadas: **8655** en **309** genomas.
- Generado (UTC): `2026-11-24T13:21:04.112881+00:00`.
- Entrada: `data/metabolism/gap_interpretation.tsv`.
- Puente EC→Rhea→productos: disponible; 75208 reacciones Rhea y 6161 EC con productos.

## Estados

| estado | brechas | % | criterio |
|---|---:|---:|---|
| `gap_evidence_supported` | 1022 | 11.8% | hay evidencia directa (mapeo Rhea validado, homólogo fuerte o dominio confirmado) de un cierre candidato; NO es un cierre |
| `gap_candidate` | 1147 | 13.3% | hay evidencia parcial (dominio HMM, homólogo débil, revisión, multifunción same-EC o asociación débil a complejo); NO es un cierre |
| `complex_suspected` | 14 | 0.2% | asociación a un complejo transferido fuerte (frac_mapped>=0.8); asociación, no cierre |
| `complex_incomplete` | 0 | 0.0% | asociación a un complejo declarado incompleto (faltan subunidades o estequiometría); no puede cerrar la brecha |
| `not_evaluated` | 7 | 0.1% | sin evidencia, sin complejo, sin transporte y sin bloqueador explícito: no hay hipótesis de cierre que evaluar |
| `not_closable_cpu` | 6465 | 74.7% | hay una capa evaluada que impide cerrar con CPU/offline: no es una brecha biosintética o requiere experimento/ΔG/base externa/GPU |
| **total** | **8655** | 100.0% | |

Precedencia: evidencia fuerte > evidencia débil > complejo fuerte > complejo incompleto > asociación débil > bloqueador > sin evaluar.

### Evidencia directa

| nivel | brechas | base |
|---|---:|---|
| fuerte | 1022 | dominio_confirmado=43, mapping_rhea_validada=257, homologo_fuerte=980 |
| débil | 1139 | dominio_hmm=535, multifuncion_misma_ec=522, homologo_debil=85, revision_soporte=74 |
| ninguna | 6494 | — |

`fuerte` significa candidato con mapeo validado, homólogo fuerte o dominio confirmado; **no** significa función demostrada.

### Filas `not_evaluated` (7)

| genoma | brecha | ChEBI | review | motivo |
|---|---|---|---|---|
| `GCF_000455605.1_GPTSA100-9` | propanoyl-CoA | `CHEBI:57392` | `one_sided_not_transport` | sin evidencia, complejo, transporte ni bloqueador explícito |
| `GCF_000695095.2_ASM69509v2` | propanoyl-CoA | `CHEBI:57392` | `one_sided_not_transport` | sin evidencia, complejo, transporte ni bloqueador explícito |
| `GCF_003353065.1_ASM335306v1` | propanoyl-CoA | `CHEBI:57392` | `one_sided_not_transport` | sin evidencia, complejo, transporte ni bloqueador explícito |
| `GCF_007750655.1_ASM775065v1` | propanoyl-CoA | `CHEBI:57392` | `one_sided_not_transport` | sin evidencia, complejo, transporte ni bloqueador explícito |
| `GCF_030623445.1_ASM3062344v1` | propanoyl-CoA | `CHEBI:57392` | `one_sided_not_transport` | sin evidencia, complejo, transporte ni bloqueador explícito |
| `GCF_055387105.1_ASM5538710v1` | propanoyl-CoA | `CHEBI:57392` | `one_sided_not_transport` | sin evidencia, complejo, transporte ni bloqueador explícito |
| `GCF_907165195.1_Cryomorphaceae_CECT30217` | propanoyl-CoA | `CHEBI:57392` | `one_sided_not_transport` | sin evidencia, complejo, transporte ni bloqueador explícito |

## Cobertura de señales

| señal | brechas |
|---|---:|
| evidencia directa (fuerte o débil) | 2161 |
| evidencia fuerte | 1022 |
| evidencia débil | 1139 |
| asociación a complejo (transfer/model) | 22 |
| señal de transporte (reacción o transportador) | 907 |
| transporte evaluado sin homólogo | 3555 |
| dirección resuelta (forward/reversible) | 439 |
| dirección no resuelta (UN/unknown) | 4379 |
| hipótesis de multifunción (misma EC) | 522 |
| sin señal positiva de ningún tipo | 5830 |
| sin ninguna capa de interpretación | 0 |

`sin señal positiva` no es lo mismo que `not_evaluated`: la fila puede tener un resultado negativo (p. ej. `no_transport`) o un bloqueador.

## Complejos

| insumo | archivos | filas | con EC | con productos | fuertes | incompletos |
|---|---:|---:|---:|---:|---:|---:|
| `transfer_*` | 400 | 19085 | 12996 | 11361 | 8408 | 1785 |
| `model_*` | 400 | 281090 | 0 | 0 | — | 0 (con `missing`) |

| asociación | brechas |
|---|---:|
| con cualquier complejo | 22 |
| con complejo fuerte (frac_mapped≥0.8) | 14 |
| con complejo incompleto | 2 |
| fuerte **y** además incompleto | 2 |
| con solapamiento Rhea | 2 |
| asociaciones totales (una brecha puede tener varias) | 24 |

### Brechas asociadas a complejos

| genoma | brecha | ChEBI | complejo | fuente | frac | estequiometría | estado |
|---|---|---|---|---|---:|---|---|
| `GCF_000179915.2_ASM17991v2` | succinyl-CoA | `CHEBI:57292` | `CPX-1092` | transfer | 1.00 | unknown | `complex_suspected` |
| `GCF_000265425.1_ASM26542v1` | succinyl-CoA | `CHEBI:57292` | `CPX-1092` | transfer | 1.00 | unknown | `complex_suspected` |
| `GCF_001027285.1_ASM102728v1` | indole | `CHEBI:16881` | `CPX-3446` | transfer | 1.00 | unknown | `complex_suspected` |
| `GCF_001676785.2_ASM167678v2` | spermidine | `CHEBI:57834` | `CPX-4384` | transfer | 0.50 | unknown | `gap_candidate` |
| `GCF_002269385.1_ASM226938v1` | indole | `CHEBI:16881` | `CPX-3446` | transfer | 1.00 | unknown | `complex_suspected` |
| `GCF_003721155.3_ASM372115v4` | spermidine | `CHEBI:57834` | `CPX-4384` | transfer | 0.50 | unknown | `gap_candidate` |
| `GCF_006352065.1_ASM635206v1` | spermidine | `CHEBI:57834` | `CPX-4383` | transfer | 0.75 | unknown | `gap_candidate` |
| `GCF_009688965.1_ASM968896v1` | daunorubicin | `CHEBI:64677` | `CPX-4443` | transfer | 1.00 | unknown | `complex_suspected` |
| `GCF_014274685.1_ASM1427468v1` | succinyl-CoA | `CHEBI:57292` | `CPX-1092` | transfer | 1.00 | unknown | `complex_suspected` |
| `GCF_016446285.1_ASM1644628v1` | daunorubicin | `CHEBI:64677` | `CPX-4443` | transfer | 1.00 | unknown | `complex_suspected` |
| `GCF_022870745.1_ASM2287074v1` | succinyl-CoA | `CHEBI:57292` | `CPX-1092` | transfer | 1.00 | unknown | `complex_suspected` |
| `GCF_027594505.1_ASM2759450v1` | indole | `CHEBI:16881` | `CPX-3446` | transfer | 1.00 | unknown | `complex_suspected` |
| `GCF_030866825.1_ASM3086682v1` | spermidine | `CHEBI:57834` | `CPX-4384` | transfer | 0.50 | unknown | `gap_candidate` |
| `GCF_032102965.1_ASM3210296v1` | spermidine | `CHEBI:57834` | `CPX-4384` | transfer | 0.50 | unknown | `gap_candidate` |
| `GCF_033547035.1_ASM3354703v1` | indole | `CHEBI:16881` | `CPX-3446` | transfer | 1.00 | unknown | `complex_suspected` |
| `GCF_037081835.1_ASM3708183v1` | succinyl-CoA | `CHEBI:57292` | `CPX-1092` | transfer | 1.00 | unknown | `complex_suspected` |
| `GCF_037081835.1_ASM3708183v1` | succinyl-CoA | `CHEBI:57292` | `CPX-3921` | transfer | 0.67 | incompleto | `complex_suspected` |
| `GCF_045784305.1_DRML_24019` | spermidine | `CHEBI:57834` | `CPX-4383` | transfer | 0.50 | unknown | `gap_candidate` |
| `GCF_054054625.1_ASM5405462v1` | succinyl-CoA | `CHEBI:57292` | `CPX-1092` | transfer | 1.00 | unknown | `complex_suspected` |
| `GCF_054054625.1_ASM5405462v1` | succinyl-CoA | `CHEBI:57292` | `CPX-3921` | transfer | 0.67 | incompleto | `complex_suspected` |
| `GCF_055396605.1_ASM5539660v1` | spermidine | `CHEBI:57834` | `CPX-4383` | transfer | 0.50 | unknown | `gap_candidate` |
| `GCF_060026435.1_ASM6002643v1` | indole | `CHEBI:16881` | `CPX-3446` | transfer | 1.00 | unknown | `complex_suspected` |
| `GCF_060030825.1_ASM6003082v1` | indole | `CHEBI:16881` | `CPX-3446` | transfer | 1.00 | unknown | `complex_suspected` |
| `GCF_900217235.1_DPRO_PRJEB22548_v1` | spermidine | `CHEBI:57834` | `CPX-4383` | transfer | 0.75 | unknown | `gap_candidate` |

**Hallazgo:** los `model_*` no aportan asociaciones en esta instantánea porque tienen 0 filas con EC. No se fuerza una asociación por vecindad: sin EC/Rhea no hay producto que comparar con el ChEBI de la brecha.

## Bloqueadores de `not_closable_cpu`

| motivo primario | brechas | significado |
|---|---:|---|
| `direction_unresolved` | 2574 | dirección no resuelta (Rhea master UN): sin ΔG/experimento no se puede afirmar que la reacción produce el metabolito |
| `catabolic` | 2390 | catabolismo/degradación: no requiere un productor biosintético; no hay cierre de brecha que buscar |
| `false_positive_direction` | 603 | falso positivo por dirección/moneda: la reacción consume el metabolito en sentido invertido; no es una brecha biosintética |
| `transport_reaction_only` | 600 | reacción de transporte en Rhea pero sin homólogo de transportador: la importación no está confirmada |
| `entry_point` | 246 | punto de entrada/importado: el metabolito puede venir del medio; no es una brecha biosintética |
| `novel_no_reference` | 43 | química sin referencia en los 400: requiere caracterización experimental o una base externa |
| `spontaneous` | 9 | posible reacción espontánea/no enzimática; falta una base curada de reacciones no enzimáticas |

Conteo por todos los motivos (una brecha puede tener varios):

| motivo | brechas |
|---|---:|
| `direction_unresolved` | 3173 |
| `catabolic` | 2390 |
| `false_positive_direction` | 603 |
| `transport_reaction_only` | 602 |
| `entry_point` | 248 |
| `novel_no_reference` | 45 |
| `transporter_homolog` | 40 |
| `host_dependent` | 16 |
| `spontaneous` | 13 |

`catabolic`, `false_positive_direction`, `entry_point` y `spontaneous` no son brechas biosintéticas: no hay cierre que buscar. El resto son brechas reales sin cierre CPU/offline posible.

Resumen: **3248** filas no son brechas biosintéticas y **3217** son brechas reales sin cierre CPU/offline.

## Límites

- El puente es **EC del complejo → Rhea → productos (ChEBI)**, con
  intersección exacta de ChEBI. No usa texto libre, estructura ni
  similitud de ligandos.
- La asociación puede incluir metabolitos moneda si una EC del complejo
  los produce; se reporta como asociación, nunca como cierre.
- `frac_mapped` es ortología K-12, no interacción física. `unknown` en
  estequiometría no es `complete`.
- `gap_evidence_supported` agrupa mapeo Rhea validado, homólogo fuerte y
  dominio confirmado; ninguno demuestra actividad en este genoma.
- No se calcula ΔG ni se ejecuta Equilibrator: `direction_unresolved`
  queda como bloqueador de esta auditoría.
- No hay GPU, AlphaFold-Multimer, co-folding ni dinámica; una predicción
  de estructura tampoco cerraría una brecha.
- `not_closable_cpu` no significa imposible: significa que con las capas
  CPU/offline auditadas no se puede emitir un cierre.
- La instantánea es de solo lectura; `gap_interpretation.tsv` y
  `data/complexes/` no se modifican.

## Reproducir

```bash
.venv/bin/python tools/audit_gap_closure.py \
    --data-root data \
    --out-md reports/gap_closure_audit.md \
    --out-json reports/gap_closure_audit.json
.venv/bin/python -m pytest tests/test_audit_gap_closure.py -q
```
