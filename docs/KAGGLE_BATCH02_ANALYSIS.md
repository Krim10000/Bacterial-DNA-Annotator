# Análisis b02_c01 y b02_c02

## Integridad

| métrica | valor |
|---|---:|
| Proteínas | **50** |
| PDB | **50** |
| Errores | **0** |
| PDB con NaN | **0** |
| PDB sin ATOM | **0** |

## Tiempos

- Total: **4.401,5 s** = 1 h 13,4 min.
- Media: **88,0 s/proteína**.
- Mediana: 80,0 s.
- Máximo: 161,2 s.

## Longitudes

- mínima: 386 aa
- media: 466 aa
- mediana: 456 aa
- máxima: 600 aa

## Calidad

| pLDDT | proteínas |
|---|---:|
| ≥0,80 | 31 |
| 0,70–0,80 | 13 |
| 0,60–0,70 | 3 |
| 0,50–0,60 | 3 |
| <0,50 | 0 |

Media pLDDT: **0,793**.

| pTM | proteínas |
|---|---:|
| ≥0,80 | 33 |
| 0,70–0,80 | 5 |
| 0,60–0,70 | 4 |
| 0,50–0,60 | 2 |
| <0,50 | 6 |

Media pTM: **0,805**.

## Cobertura

- 50 accesiones únicas.
- **184 clusters**.
- **400 proteínas**.
- **120 genomas**.

## Hallazgos estructurales

- 50/50 tienen hit no-self.
- 42/50 son de alta confianza (pLDDT ≥0,70 y pTM ≥0,50).

Categorías de los 42:

| categoría | cantidad |
|---|---:|
| otros | 23 |
| uncharacterized/DUF | 8 |
| transporte | 5 |
| RNA metiltransferasa | 4 |
| aminotransferasa | 1 |
| peptidasa | 1 |

## Casos destacados

| query | target | identidad | lectura |
|---|---|---:|---|
| Q8FCN5 | P59195 YhhS | **93%** | transportador MFS |
| Q8CRU6 | Q5HN37 | 64% | RNA metiltransferasa |
| P9WJW7 | P9WJW6 | 56% | transportador MFS |
| P74311 | A6TP80 Acr3 | 45% | resistencia a arsenicales |
| O94562 | A0A098DDI1 | 45% | aminotransferasa |
| P38423 | Q2S280 UvrC | 42% | reparación de ADN |
| O34708 | Q6LX42 AgcS | 40% | symporter sodio/alanina |
| P32703 | A8FHE3 | 36% | antiportador |

## Casos de baja confianza

- P9WI44: pLDDT 0,56, pTM 0,41.
- P9WI08: pLDDT 0,60, pTM 0,49.
- P9WI02: pLDDT 0,58, pTM 0,41.

Son PPE family; no usar para conclusiones fuertes.

## Interpretación

- El lote b02 tiene mejor pLDDT que el lote 01.
- Predominan transportadores, RNA metiltransferasas y proteínas PPE.
- Q8FCN5 es el caso más claro: MFS transporter con 93% de identidad.
- Siguen siendo hipótesis, no funciones confirmadas.
