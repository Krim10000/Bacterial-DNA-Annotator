# Candidatos a subunidades de complejos

## Límite

ESMFold predice **monómeros**. No modela interfaz, ensamblaje ni estequiometría.

Un modelo monomérico no descarta que la proteína sea subunidad.

## Búsqueda

Sobre los 304 hits no-self de los lotes 01 y 02:

- 58 tienen descripción compatible con complejo/subunidad/transporte;
- 42 tienen pLDDT ≥0,70 y pTM ≥0,50.

## Candidatos fuertes

| query | hit | lectura |
|---|---|---|
| P75128 | PheT | fenilalanil-tRNA sintetasa, subunidad beta |
| P38423 | UvrC | complejo UvrABC de reparación |
| Q57860 | PolC | DNA polimerasa III |
| Q58207 | Rpo5 | RNA polimerasa |
| P41883 | Rpn11 | proteasoma 26S |
| Q5UQ92 | RecD | complejo RecBCD |
| P44104 | EcoP15I Res | restricción-modificación |
| Q45062 | GlyS | glicil-tRNA sintetasa beta (baja confianza) |

## Ribosomales

| query | hit | pLDDT |
|---|---|---:|
| O27952 | uS11c | 0,73 |
| O05398 | uL18 | 0,79 |
| P76463 | uL6 | 0,76 |
| C0H410 | bL31 | 0,54 |
| P92540 | bS16 | 0,41 |
| O83479 | uS5 | 0,54 |

Son subunidades ribosomales probables, pero varias tienen baja confianza.

## Transporte

| query | hit | lectura |
|---|---|---|
| P77795 | CysA | ATPasa de importe sulfato |
| P37494 | CysA | idem |
| P39911 | ModC | ATPasa de importe molibdato |
| O33683 | MnhG1 | subunidad del antiportador Na+/H+ |
| O30083 | EfeU | permeasa de hierro |
| O34649 | RibX | permeasa de riboflavina |

## PPE

- P9WI44, P9WI08, P9WI02.
- Familia PPE de micobacterias.
- Posible relación con sistemas ESX.
- pLDDT/pTM bajos: no concluir.

## Cómo validar

1. Identificar posible pareja por operón/vecindad.
2. Buscar anotación de complejo.
3. AlphaFold-Multimer/ColabFold con la pareja.
4. Medir ipTM/PAE.
5. Comparar con PDB biological assembly.
6. No confundir pTM monomérico con confianza de complejo.
