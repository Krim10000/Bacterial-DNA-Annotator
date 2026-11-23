# Reducir lo desconocido sin inventar función

## Estado de implementación (2026-11-23)

Implementado sin GPU:

- filtros baratos y deduplicación/clustering;
- F1 Foldseek local y F2 EnzyMM/M-CSA;
- gating por pLDDT/pTM y early-exit si F1/F2 no da señal
  (`webapp/structural_filters.py`);
- priorización y propagación a clusters;
- documentación de límites y categorías.

Pendiente:

- completar estructuras (2.570 ≤600 aa y 169 >600 aa);
- separar las colas baratas/caras en comandos independientes;
- asignación automática de categorías amplias;
- GPU de primer filtro en Kaggle para casos que lo justifiquen.

## Regla

No se “cierra” una proteína por una sola señal. Se reduce la incertidumbre por capas.

```text
secuencia
  → dominio/familia
  → topología/membrana
  → contexto génico
  → estructura
  → filtros catalíticos
  → hipótesis
  → revisión
```

## 1. Terminar estructura

- 2.570 proteínas ≤600 aa en curso.
- 169 >600 aa aparte.
- Control: pLDDT ≥0,70 y pTM ≥0,50 para interpretar.
- pLDDT/pTM bajos → no concluir.

## 2. Filtros baratos antes de estructura

Sin GPU:

- Pfam/HMMER;
- InterPro/DUF;
- TMHMM/Phobius;
- péptido señal;
- coiled-coil;
- bajo complejidad;
- taxonomía;
- vecindad génica/operón.

## 3. Categorías amplias

| categoría | evidencia mínima |
|---|---|
| enzima probable | dominio catalítico + F2/F3/F4 |
| transportador probable | ≥2 TM + familia de transporte |
| regulador probable | dominio de unión a ADN + contexto |
| proteína de membrana | TM/singal péptido |
| desordenada | pLDDT bajo + baja complejidad |
| desconocida | sin señal suficiente |

Son **categorías**, no funciones concretas.

## 4. Enzimas

- F2 EnzyMM/M-CSA: plantilla/EC.
- F3 P2Rank: bolsillo.
- F4 AlphaFill: cofactor.
- F5: solo con mapa propio.
- Contexto: operón/ruta.

## 5. Transportadores

- TMHMM/Phobius: número de hélices.
- Pfam: MFS, APC, ABC_tran, NBC, etc.
- Foldseek: familia estructural.
- Contexto: operón, sustrato cercano.
- Límite: ESMFold no modela membrana; pLDDT puede ser bajo.
- No se puede asignar sustrato solo por estructura.

## 6. Reguladores

- Pfam: HTH, LysR, TetR, AraC, CRP, LuxR, CheY, sigma, etc.
- Foldseek: homólogo estructural.
- Contexto: genes vecinos, promotor, operón.
- Límite: un dominio HTH no confirma regulón ni gen diana.

## 7. Propagación

- Se calcula en el representante.
- Se propaga al cluster.
- Se marca `inferred_from_cluster`.
- Nunca `confirmed`.

## 8. Estados

```text
unknown
domain_only
membrane_only
enzyme_hypothesis
transporter_hypothesis
regulator_hypothesis
supported_hypothesis
conflict
```

## 9. Límite honesto

- No todas las proteínas se podrán anotar.
- Algunas quedarán `unknown`.
- ESMFold no confirma función.
- AlphaFill no confirma sustrato.
- P2Rank no confirma sitio activo.
- F5 ausente no es negativo.

## 10. Orden de trabajo

```text
1. Completar 2.570 estructuras.
2. QC pLDDT/pTM.
3. Dominios/Pfam y TM.
4. F1 local.
5. F2.
6. F3/F4 si hay señal.
7. Categoría amplia.
8. Revisión.
9. Propagar.
10. 169 largas.
```
