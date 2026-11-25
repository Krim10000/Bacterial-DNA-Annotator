# Análisis del lote 01 de Kaggle

## Integridad

| métrica | valor |
|---|---:|
| Proteínas | **257** |
| PDB | **257** |
| `done` | **257** |
| Errores | **0** |
| PDB con `nan` | **0** |
| PDB sin ATOM | **0** |
| PDB pequeños | **0** |

## Tiempos

| métrica | valor |
|---|---:|
| Tiempo total | **7.423,7 s** = 2 h 3,7 min |
| Media | **28,9 s/proteína** |
| Mediana | 16,1 s |
| Mínimo | 1,81 s |
| Máximo | 156,4 s |

## Longitudes

- mínima: 39 aa
- media: 257 aa
- mediana: 239 aa
- máxima: 600 aa

## pLDDT

| rango | proteínas | % |
|---|---:|---:|
| ≥0,80 | 125 | 48,6% |
| 0,70–0,80 | 78 | 30,4% |
| 0,60–0,70 | 30 | 11,7% |
| 0,50–0,60 | 16 | 6,2% |
| <0,50 | **8** | **3,1%** |

Media: **0,761**  
Mediana: **0,79**

## pTM

| rango | proteínas | % |
|---|---:|---:|
| ≥0,80 | 150 | 58,4% |
| 0,70–0,80 | 34 | 13,2% |
| 0,60–0,70 | 30 | 11,7% |
| 0,50–0,60 | 15 | 5,8% |
| <0,50 | **28** | **10,9%** |

Media: **0,775**  
Mediana: **0,842**

## Casos de baja confianza

No usar para conclusiones estructurales fuertes:

```text
P56728  pLDDT 0.31  pTM 0.19
Q9S4W9  pLDDT 0.32  pTM 0.20
O67076  pLDDT 0.34  pTM 0.32
O67138  pLDDT 0.37  pTM 0.21
P92540  pLDDT 0.41  pTM 0.25
O67369  pLDDT 0.47  pTM 0.34
Q52968  pLDDT 0.47  pTM 0.21
Q5XD68  pLDDT 0.49  pTM 0.25
```

## Lectura

- 79% tienen pLDDT ≥0,70.
- 90% tienen pLDDT ≥0,60.
- 11% tienen pTM <0,50.
- La mayoría sirve para hipótesis, no para estructura fina.
- El lote 01 está terminado y es consistente.

## Proyección

- 257 en 2,06 h.
- 10 lotes ≈ **20,6 h**.
- Con Kaggle 30 h/semana, factible en **1–1,5 semanas**.

## Cobertura

- 257 accesiones únicas.
- 792 clusters.
- 1.659 proteínas.
- 262 genomas.

## Hallazgo estructural

- 254/257 tienen al menos un hit no-self.
- 201/257 tienen pLDDT ≥0,70 y pTM ≥0,50.
- Los mejores hits no-self:

| query | target | identidad | lectura |
|---|---|---:|---|
| P0AA74 | P0AA73 YhbE | 93,4% | transportador de membrana |
| A8AII1 | Q3Z3M2 YcaD | 92,9% | transportador MFS |
| P0ADA6 | P0ADA5 YajG | 91,6% | lipoproteína |
| P0AAL7 | P0AAL8 YdhY | 89,9% | ferredoxina-like |
| P0ADK2 | P0ADK3 YiaF | 89,4% | proteína no caracterizada |
| Q4L749 | Q49YD7 | 69,3% | peptidasa |
| P0ADD8 | P0ADD7 YjjQ | 69,2% | posible regulador |
| P9WGR3 | P9WGR2 | 67,6% | oxidorreductasa |
| Q9ZCD6 | Q7NQN3 UbiD | 60,9% | descarboxilasa |
| Q53196 | E1V7V7 | 58,8% | aminotransferasa |


Categorías de los 201 de alta confianza:

- 78 otros;
- 50 uncharacterized/DUF;
- 24 transporte;
- 21 oxidorreductasa;
- 11 regulador;
- 10 RNA metiltransferasa;
- 4 peptidasa;
- 2 descarboxilasa;
- 1 aminotransferasa.

## Interpretación

- Los hits no-self evitan el auto-match.
- Identidad alta (≥90%) apoya familia, no función exacta.
- Muchos hits siguen siendo “uncharacterized”.
- Son hipótesis estructurales, no cierres.
