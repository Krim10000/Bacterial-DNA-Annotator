# Revisión de los 50 candidatos `A_validated_same`

Criterio: dominio específico y homología DIAMOND en **la misma proteína**, dirección Rhea direccional.

Archivo: `data/metabolism/CANDIDATES_50_SAME.tsv`

## Grupos

| candidatos | brecha | dominio | homólogo | lectura |
|---:|---|---|---|---|
| 10 | L-cysteinylglycine | G_glu_transpept | P15557 | **conflicto**: producto anotado oxamate amidohidrolasa/acilasa; no GGT clásica |
| 10 | (2R,3R)-2,3-dihidroxi-3-metilpentanoato | KARI_N | P06168 | **KARI consistente**; mejor candidato |
| 7 | adenylyl-molybdopterin | MoCF_biosynth | Q03555 | biosíntesis de MoCo |
| 4 | adenylyl-molybdopterin | MoCF_biosynth | Q39054 | biosíntesis de MoCo |
| 2 | adenylyl-molybdopterin | MoCF_biosynth | Q54KM0 | biosíntesis de MoCo |
| 4 | D-threo-isocitrate | Aconitase | O13966 | aconitasa probable |
| 4 | molybdopterin | MoaE | A2X0R4 | MoaE probable |
| 2 | quinolinate | QRPTase_C | A0A1S4D475 | quinolato fosforribosiltransferasa |
| 1 | quinolinate | QRPTase_C | A0A1S4CL59 | idem |
| 1 | quinolinate | QRPTase_C | B2RFS9 | idem |
| 1 | 3-(4-hidroxifenil)piruvato | 2-Hacid_dh_C | J7SHB8 | oxidorreductasa probable |
| 1 | 3-(4-hidroxifenil)piruvato | Ldh_1_C | P40925 | deshidrogenasa probable |
| 1 | molybdopterin | MoaE | A4FUY7 | MoaE probable |
| 1 | 3',3'-c-di-GMP | GGDEF | A0A0H3LVU7 | diguanilato ciclasa probable |

## Top 10

Hay dos grupos de 10:

### L-cysteinylglycine / G_glu_transpept / P15557

- domain score 489–532;
- e-value 1e-150 a 1e-163;
- identidad 31–35%;
- cobertura ~1.0;
- misma proteína;
- **pero** el producto anotado es oxamate amidohidrolasa/acilasa, no GGT clásica.

Es un **conflicto de anotación**, no un cierre limpio.

### KARI / KARI_N / P06168

- identidad 55–60%;
- cobertura ~1.0;
- producto anotado: KARI;
- dominio y homólogo consistentes.

Es el **mejor candidato** del conjunto.

## Interpretación

- `A_validated_same` es la evidencia más fuerte de los 508.
- Sigue siendo **candidato**, no `confirmed`.
- La recurrencia en varios genomas refuerza la hipótesis.
- No se ha validado experimentalmente ni por estructura/docking.

## Regla

1. Revisar los 10 de L-cysteinylglycine.
2. Revisar los 10 KARI.
3. Revisar los 13 MoCF/MoaE.
4. Luego Aconitase, quinolinate y los demás.
5. No anotar como confirmado sin validación.
