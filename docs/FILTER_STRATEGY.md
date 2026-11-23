# Estrategia de filtros: agrupar y abaratar antes

## Principio

No ejecutar filtros caros sobre proteínas duplicadas ni sobre modelos mal plegados.

```text
deduplicar
  → agrupar
  → filtrar barato
  → evaluar estructura
  → filtros estructurales
  → revisión
```

## 1. Agrupar antes de calcular

Orden:

1. deduplicar accesión UniProt;
2. clusterizar secuencia 90/70/50/30%;
3. agrupar por dominio/familia;
4. estratificar por taxonomía/genoma;
5. conservar un representante por grupo.

Resultados ya obtenidos:

| conjunto | secuencias | clusters 50% |
|---|---:|---:|
| sin hit Swiss-Prot | 647.588 | 446.927 |
| uncharacterized con UniProt | 54.031 | 25.058 |
| únicos UniProt de esos clusters | — | 3.931 |
| con AFDB | — | 1.192 |

Regla:

> Se filtra el representante, no cada copia genómica.

## 2. Filtros baratos primero

| orden | filtro | coste | cuándo |
|---:|---|---|---|
| 0 | anotación existente | barato | siempre |
| 1 | dominio/Pfam/familia | barato | antes de GPU |
| 2 | contexto génico/operón | barato | candidatos |
| 3 | F1 Foldseek local | barato-medio | estructura disponible |
| 4 | F2 EnzyMM/M-CSA | medio | tras F1 |
| 5 | F3 P2Rank | medio | tras F1 |
| 6 | F4 AlphaFill | medio | tras F1 |
| 7 | F5 sitio activo propio | caro/escaso | solo si hay mapa |

## 3. Control de calidad estructural

Antes de interpretar:

| pLDDT | pTM | uso |
|---:|---:|---|
| ≥0,70 | ≥0,50 | usar para hipótesis |
| 0,60–0,70 | 0,40–0,50 | revisar |
| <0,60 | <0,40 | no concluir |

- pLDDT bajo en un dominio: no interpretar ese dominio.
- pTM bajo: no concluir topología global.
- Recortar colas de baja confianza antes de P2Rank/AlphaFill.

## 4. Regla de evidencia

```text
F2 + F3 + F4 + contexto + pLDDT alto → hipótesis fuerte
F2 solo → hipótesis débil
pLDDT bajo → no concluir
F5 solo → no cierra
```

Ningún filtro cierra una brecha. Ninguno es `confirmed`.

## 5. Estado real del código

Implementado:

- clustering y manifiestos;
- F1 local;
- F2/F3/F4/F5 en `run_filters`;
- propagación a clusters;
- priorización;
- **corte automático por pLDDT/pTM** (`quality_gate`);
- **early-exit** si F1/F2 no da señal (`gating=no_structural_signal`);
- **gating de F3/F4** tras F1/F2 y calidad estructural;
- registro de `plddt`, `ptm`, `quality_gate` y `gating` en la tabla de filtros.

Pendiente de ingeniería, no de ciencia:

- separar físicamente las colas baratas y caras en etapas de CLI independientes;
- exponer el gating en un comando dedicado, sin recalcular F1/F2.

La regla ya no es solo documental: el pipeline no gasta P2Rank/AlphaFill con
pLDDT/pTM bajos ni cuando F1/F2 no tienen señal.

## 6. Uso práctico

```text
1. elegir representantes
2. pLDDT/pTM
3. F1 local
4. F2
5. F3/F4 solo si F2 o F1 da señal
6. F5 solo con mapa propio
7. revisión
```

## 7. Límites

- ESMFold no confirma función.
- pLDDT/pTM no son evidencia funcional.
- F3 bolsillo no es sitio activo.
- F4 cofactor no define reacción.
- F5 ausente no es negativo.
