> **NOTA IMPORTANTE:** la evaluación de subunidades queda congelada hasta terminar las 2.570 estructuras. Ver `IMPORTANT_NOTE.md`.

# Pares para AlphaFold-Multimer por vecindad génica

## Lista cruda

`MULTIMER_PAIRS.tsv`: 1.823 pares query–vecino (±3 genes).

Incluye mucho ruido: genes cercanos no siempre forman complejo.

`MULTIMER_PAIRS_COMPLEX.tsv`: 60 pares con descripción compatible con complejo/subunidad.

## Lista curada

`MULTIMER_TOP_PAIRS.tsv`: 8 pares con contexto y estructura coherentes.

| query | partner | lectura |
|---|---|---|
| P44104 | P44105 | sistema tipo III Res/Mod |
| P44104 | P71366 | Mod de tipo III |
| P77795 | P0AFR9 YdcV | permeasa ABC |
| P77795 | P77156 YdcU | permeasa ABC |
| P77795 | P76108 YdcS | proteína periplásmica ABC |
| Q58207 | Q8RFV0 LolD | ATPasa LolCDE |
| P37494 | P59653 ComA | ATPasa de transporte |
| P55680 | P55679 Y4WA | proteasa de zinc |

FASTA:

`data/metabolism/targets/multimer_fastas/*.fasta`

## Descarte por vecindad

- P75128 PheT: sus vecinos son TyrS, no PheS.
- P38423 UvrC: no hay UvrA/UvrB cerca.
- Q57860 PolC: no hay subunidades de Pol III cerca.
- P41883 Rpn11: no hay proteasoma cerca.
- Q5UQ92 RecD: no hay RecB/RecC cerca.

La estructura sola no basta; sin pareja génica no se justifica el complejo.

## AlphaFold-Multimer en Kaggle

Sí, con ColabFold y T4.

```python
!pip install -q "colabfold[alphafold]"
!colabfold_batch pair.fasta out \
  --model-type alphafold2_multimer_v3 \
  --num-models 1 --num-recycle 3
```

- Internet activo.
- MSA vía servidor MMseqs2.
- 2 cadenas pequeñas.
- T4 16 GB: evita complejos >900 aa totales.
- P37494+ComA (926 aa) y P55680+Y4WA (1.005 aa) pueden dar OOM.
- Revisar ipTM/PAE, no solo pTM.

## Regla

- Vecindad no es prueba de complejo.
- Coexpresión/operón conservado sí ayuda.
- AlphaFold-Multimer da una hipótesis de interfaz.
- No confirma interacción in vivo.
