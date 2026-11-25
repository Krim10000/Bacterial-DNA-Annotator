# 647.588 sin hit Swiss-Prot: plan CPU

## Aclaración

Sin hit Swiss-Prot no significa sin homología.

Swiss-Prot es una base curada pequeña. Un ORF puede tener:

- homología en TrEMBL/UniRef;
- dominio Pfam;
- familia estructural;
- contexto génico;
- señales de membrana;
- taxonomía.

## Estado

| métrica | valor |
|---|---:|
| Sin hit Swiss-Prot | 647.588 |
| Clusters 50% | 446.927 |
| Singletons | 360.741 |
| Clusters >1 | 86.186 |
| Clusters >5 | 9.521 |

No se pueden anotar profundamente las 647k con CPU en tiempo razonable. Hay que reducir y priorizar.

## Opciones CPU, de barata a cara

### 1. Reclustering y representantes

- 30% o 50%;
- dominio/familia;
- taxonomía;
- elegir representantes.

Objetivo: reducir antes de HMMER.

### 2. Búsqueda contra bases mayores

- MMseqs2/DIAMOND contra **UniRef90/TrEMBL**;
- contra **PDB seq**;
- contra **eggNOG**;
- no solo Swiss-Prot.

Salida: homología, familia, posible EC.

### 3. Dominios

- Pfam/InterPro/pyhmmer;
- TIGRFAM/NCBIfam;
- CDD.

Coste alto si se corre sobre 647k. Hacerlo sobre representantes y sin hit previo.

### 4. Topología

- TMHMM/Phobius;
- péptido señal;
- lipoproteína;
- coiled-coil;
- baja complejidad.

Salida: membrana, secretada, soluble.

### 5. Contexto génico

- operón;
- vecindad conservada;
- sintenía;
- clusters biosintéticos;
- sistemas de defensa/móviles.

Salida: hipótesis por vecindad, no función confirmada.

### 6. Taxonomía y conservación

- distribución en 400 genomas;
- filo;
- core/variable;
- ganancia/pérdida.

### 7. Embeddings

- ESM-2/ProtT5 CPU solo para representantes;
- clustering y outliers;
- no para 647k completas.

### 8. Estructura

- solo después de priorizar;
- ESMFold/AFDB/AlphaFold para top;
- no para 647k.

## Estrategia realista

```text
647k
 → 446k clusters
 → representantes por dominio/taxonomía
 → MMseqs2 vs UniRef/PDB/eggNOG
 → Pfam/TMHMM en los que sigan sin señal
 → contexto génico
 → catálogo de hipótesis
 → top para estructura
```

## Anotación “completa”

No significa función confirmada. Significa:

```text
familia
dominio
membrana/señal
homología
contexto
estructura si aplica
confianza
limitación
```

Algunas quedarán `unknown`.

## Coste

- CPU local: $0.
- Tiempo: días/semanas según profundidad.
- GPU solo para el top.

## Regla

- No correr HMMER profundo sobre 647k sin reducir.
- No llamar “completa” a una hipótesis.
- No plegar 647k.


## Caminos

Ver `docs/THREE_PROTEIN_PATHS.md`: A función nombrada, B uncharacterized/DUF, C sin hit Swiss-Prot.
