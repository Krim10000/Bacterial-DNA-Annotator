# Análisis acumulado: 587 estructuras

## Lotes

| lote | proteínas |
|---|---:|
| Lote 01 | 257 |
| b02_c01–c08 | 200 |
| `esmfold_300_short` parcial | **130** |
| **Total único** | **587** |

No hay duplicados entre ZIP.

## Integridad

| métrica | valor |
|---|---:|
| PDB | 587 |
| `done` | 587 |
| Errores | **0** |
| PDB con NaN | **0** |
| PDB sin ATOM | **0** |

## Tiempo

- Total: **20.380 s** = **5 h 39,7 min**.
- Media: **34,7 s/proteína**.
- Mediana: **20,9 s**.
- Máximo: 164,5 s.

## Longitudes

- mínima: 37 aa
- media: 279 aa
- mediana: 265 aa
- máxima: 600 aa

## Calidad

| pLDDT | proteínas | % |
|---|---:|---:|
| ≥0,80 | 299 | 50,9% |
| 0,70–0,80 | 168 | 28,6% |
| 0,60–0,70 | 74 | 12,6% |
| 0,50–0,60 | 30 | 5,1% |
| <0,50 | **16** | **2,7%** |

Media pLDDT: **0,767**.

| pTM | proteínas | % |
|---|---:|---:|
| ≥0,80 | 346 | 58,9% |
| 0,70–0,80 | 69 | 11,7% |
| 0,60–0,70 | 70 | 11,9% |
| 0,50–0,60 | 34 | 5,8% |
| <0,50 | **68** | **11,6%** |

Media pTM: **0,775**.

## Cobertura

- 587 accesiones únicas;
- **1.939 clusters**;
- **3.845 proteínas**;
- **330 genomas**.

Porcentajes:

- 587/2.570 = **22,8%** de las cortas;
- 1.939/25.058 = **7,7%** de clusters;
- 3.845/54.031 = **7,1%** de proteínas uncharacterized.

## Casos de baja confianza

No usar para conclusiones fuertes:

```text
P56728  pLDDT 0,31  pTM 0,19
Q9S4W9  pLDDT 0,32  pTM 0,20
O67076  pLDDT 0,34  pTM 0,32
O67138  pLDDT 0,37  pTM 0,21
Q58674  pLDDT 0,41  pTM 0,31
Q00247  pLDDT 0,41  pTM 0,30
P92540  pLDDT 0,41  pTM 0,25
O83451  pLDDT 0,42  pTM 0,36
Q9Z7M9  pLDDT 0,42  pTM 0,26
O83729  pLDDT 0,44  pTM 0,28
```

## Interpretación

- El 79,5% tiene pLDDT ≥0,70.
- El 11,6% tiene pTM <0,50.
- Los 587 PDB son válidos, pero no todos sirven para conclusiones finas.
- Faltan 1.983 cortas y 169 largas.
- Los hallazgos estructurales previos (transportadores, oxidorreductasas, reguladores) siguen siendo hipótesis.

## Siguiente

- Terminar `esmfold_300_short`;
- luego 800;
- luego 1.013;
- después 169 largas.

## Hallazgos estructurales

Foldseek local sobre las 587:

- 578/587 tienen hit no-self.
- 458/587 son de alta confianza (pLDDT ≥0,70 y pTM ≥0,50).
- 9 no tienen hit no-self.

Categorías de los 458:

| categoría | cantidad |
|---|---:|
| otros | 198 |
| uncharacterized/DUF | 113 |
| transporte | **46** |
| oxidorreductasa | **33** |
| RNA metiltransferasa | **25** |
| regulador | **24** |
| peptidasa | 8 |
| descarboxilasa | 5 |
| aminotransferasa | 3 |
| ribosomal | 3 |

## Casos destacados

| query | hit | identidad | lectura |
|---|---|---:|---|
| Q8FCN5 | P59195 YhhS | 93% | transportador MFS |
| P0AA74 | P0AA73 YhbE | 93% | transportador de membrana |
| A8AII1 | Q3Z3M2 YcaD | 93% | transportador MFS |
| P38423 | Q2S280 UvrC | 42% | reparación UvrABC |
| P77795 | CysA | — | ATPasa ABC |
| P39911 | ModC | — | ATPasa ABC |
| P0ADD8 | P0ADD7 YjjQ | 69% | posible regulador |
| Q9ZCD6 | Q7NQN3 UbiD | 61% | descarboxilasa |
| Q53196 | E1V7V7 | 59% | aminotransferasa |
| Q4L749 | Q49YD7 | 69% | peptidasa |

## Lectura

- Predominan transporte, oxidorreductasas, RNA metiltransferasas y reguladores.
- Muchos hits siguen siendo “uncharacterized”.
- Identidad alta apoya familia, no función exacta.
- No hay cierres confirmados.
