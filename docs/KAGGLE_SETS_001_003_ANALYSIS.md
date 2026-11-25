# Análisis sets 001–003

## Integridad

| métrica | valor |
|---|---:|
| Proteínas nuevas | **755** |
| PDB | 755 |
| Errores | **0** |
| PDB con NaN | **0** |
| PDB sin ATOM | **0** |

## Tiempos

| set | proteínas | tiempo |
|---|---:|---:|
| set_001 | 290 | 12,0 min |
| set_002 | 282 | 21,2 min |
| set_003 | 183 | 22,6 min |
| **Total** | **755** | **55,8 min** |

## Calidad

- pLDDT medio: **0,718**.
- pTM medio: **0,710**.
- pLDDT <0,50: **81**.
- pTM <0,50: **228**.

Son proteínas cortas: la confianza global baja, pero muchas tienen pLDDT/pTM aceptables.

## Cobertura acumulada

- 1.342 proteínas.
- **3.825 clusters**.
- **8.116 proteínas**.
- **358 genomas**.

Progreso:

- 1.342/2.570 = **52,2%** de las cortas.
- Restan 1.228 cortas.
- 169 largas siguen aparte.

## Hallazgos estructurales nuevos

- 709/755 tienen hit no-self.
- 400/755 son de alta confianza (pLDDT ≥0,70 y pTM ≥0,50).

Categorías de los 400:

| categoría | cantidad |
|---|---:|
| otros | 224 |
| **uncharacterized/DUF** | **120** |
| oxidorreductasa | 17 |
| transporte | 14 |
| regulador | 13 |
| descarboxilasa | 4 |
| ribosomal | 3 |
| RNA metiltransferasa | 2 |
| peptidasa | 2 |
| aminotransferasa | 1 |

## Hallazgo principal

Los sets 001–003 están dominados por **proteínas cortas conservadas no caracterizadas**:

- alta identidad con homólogos `uncharacterized` de E. coli/Shigella;
- familias pequeñas, muchas veces de 60–150 aa;
- pLDDT/pTM bajos por tamaño;
- no son basura: son familias conservadas.

Ejemplos:

```text
P64623  YhdV  100%  pLDDT 0,45
P0ADU9  YqiB   96%  pLDDT 0,66
P64543  YpeC   95%  pLDDT 0,44
P0AD11  YecJ   94%  pLDDT 0,68
P0AAN7  YaiA   94%  pLDDT 0,73
P65293  YgdI   92%  pLDDT 0,65
```

No son funciones; son **familias conservadas desconocidas**.

## Casos con posible función

| query | hit | identidad | lectura |
|---|---|---:|---|
| P04394 | NuoE | 41% | NADH-quinona oxidorreductasa |
| P26946 | NatA | 39% | ATPasa ABC |
| P44580 | Lrp | 36% | regulador tipo Lrp |
| C0H3S9 | ArsR | 31% | represor ArsR |
| Q58906 | 4-carboximuconolactona descarboxilasa | 34% | posible descarboxilasa |
| P72581 | pseudouridina sintasa E | 35% | posible modificadora de rRNA |

Identidad baja: solo hipótesis de familia.

## Interpretación

- Este bloque refuerza el valor de los `uncharacterized/DUF`.
- Hay muchas familias cortas conservadas, no anotadas.
- La estructura ayuda a agruparlas, no a nombrarlas.
- No usar los pLDDT bajos para conclusiones finas.
- La señal más sólida es la conservación de secuencia, no la estructura.
