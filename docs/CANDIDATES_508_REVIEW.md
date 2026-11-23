# Revisión priorizada de 508 cierres candidatos

## Fuentes

- `domain_closure_candidates_top1000.tsv`: dominio específico + dirección.
- `closure_final_top1000.tsv`: cruce dominio × homología.
- `homology_closure_top1000.tsv`: DIAMOND independiente.
- `closure_review.tsv`: revisión manual disponible.

## Tiers

| tier | filas | significado |
|---|---:|---|
| A_validated_same | 50 | dominio y homología en la **misma proteína** |
| A_validated_genome | 101 | dominio + homología en el genoma, distinta proteína |
| B_consensus | 314 | dominio específico + dirección; sin homología |
| D_anomalous | 41 | dominio/homología en conflicto |
| D_reject | 2 | dominio amplio; revisión manual rechazó |

Total: **508**.

## Interpretación

- `validated` aquí significa **dos métodos independientes coinciden**, no experimento.
- `same_protein` es la señal más fuerte: 50 casos.
- `consensus` es dominio + dirección; sigue siendo hipótesis.
- `anomalous` no debe anotarse sin revisión.
- Ninguno es `confirmed`.

## Ejemplos top

| tier | brecha | dominio | homología | pident | lectura |
|---|---|---|---:|---:|---|
| A_same | L-cysteinylglycine | G_glu_transpept | P15557 | 32,8% | γ-glutamiltranspeptidasa probable |
| A_genome | D-threo-isocitrate | Aconitase | O13966 | 36,8% | aconitasa probable |
| A_genome | UDP-GlcA | UDPG_MGDP_dh_N | Q9RHD7 | 56,0% | deshidrogenasa probable |
| B_consensus | D-threo-isocitrate | Aconitase | — | — | dominio + dirección |
| D_anomalous | D-serine | PALP | A2XWA9 | 38,7% | conflicto; revisar |

## Archivo

`data/metabolism/CANDIDATES_508_REVIEW.tsv`

Columnas:

```text
tier genome gap chebi role support_400 domain domain_score domain_evalue
direction rhea_master query homology_acc homology_pident homology_qcov
homology_protein same_protein domain_verdict verdict final_status
manual_verdict manual_reason
```

## Regla

- A_same: prioridad máxima de revisión.
- A_genome: revisar contexto génico.
- B_consensus: revisar dominio/dirección.
- D: no usar sin validación.
- Nunca `confirmed`.
