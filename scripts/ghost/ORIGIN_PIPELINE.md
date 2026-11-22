# Pipeline de origen de replicación

## Principio

Cada anotación debe declarar **qué se intentó, con qué métodos, con qué evidencia
y por qué se acepta o se marca**. El pipeline no es un script para 400 genomas:
es un módulo que recibe secuencia, marcadores y contexto, y devuelve un estado
explicable.

## Árbol de decisión

1. **Contexto**: filo, tamaño, estilo de vida (endosimbionte/parásito), replicones.
2. **Marcadores**: `dnaA`/`dnaN` (HMM + homología), RPA (Fusobacteriaceae),
   cajas DnaA.
3. **Métodos**: ORCA multi-configuración (ambos, dnaA, dnaN, ninguno), GC skew
   acumulado, Z-curve, densidad de cajas dnaA, DoriC como verdad externa.
4. **Consenso**:
   - `validated`: ≥2 métodos independientes coinciden (ORCA + DoriC,
     DoriC + GC skew, GC skew + cajas).
   - `consensus`: ORCA + proximidad a dnaA, sin verdad externa.
   - `lineage_specific`: desviación esperada del linaje (Fusobacteriaceae sin
     `dnaA`, con RPA/iterones).
   - `anomalous`: los métodos discrepan; se declara el motivo.
   - `unresolved`: ningún método dio candidato; se listan los intentos.
5. **Salida**: posición, estado, `validated_by`, `reason`, `attempted`,
   confianza; atributos en el GFF.

## Mensajes al usuario (ejemplos)

- «Origen validado por **ORCA + DoriC** (Δ=0 pb)».
- «Origen **anómalo**: ORCA discrepa de DoriC (Δ=1,3 Mb); se intentó GC skew;
  GC skew apoya DoriC → **validado por DoriC + GC skew**».
- «**Anómalo esperado**: sin `dnaA` (Fusobacteriaceae); se espera RPA/iterones;
  se intentó GC skew y búsqueda de cajas».

## Bibliografía

- Brezellec & Pasek (2026, bioRxiv): *A plasmid-derived machinery drives
  chromosomal replication initiation in Fusobacterium nucleatum and its close
  relatives*. **Fusobacteriaceae perdió `dnaA`**; usa RPA e iterones.
- Karpathy et al. (2007, PLoS One): genoma de *F. nucleatum* subsp.
  *polymorphum*; origen tentativo junto a un ORF con baja identidad a DnaA,
  `recF`/`gyrB` aguas abajo.
- Gil et al. (2003, PNAS): *Blochmannia floridanus*; oriC por **GC skew** porque
  no hay un cluster diagnóstico de cajas DnaA.
- *Mycoplasma conjunctivae* (2009): replicación simplificada; oriC por
  comparación de genomas.
- ORCA (2024): disparidades de nucleótidos, cajas DnaA y posición de genes.
- DoriC (2007): Z-curve, cajas DnaA, genes adyacentes y filogenia.

## Otras familias/linajes con iniciación independiente de DnaA

- **Fusobacteriaceae** (Fusobacterium, Leptotrichiaceae): pérdida de `dnaA`;
  RPA/iterones derivados de plásmido ([Brezellec & Pasek 2026](https://www.biorxiv.org/content/10.64898/2026.10.08.757545)).
- **Endosimbiontes**: *Blochmannia*, *Buchnera*, *Wigglesworthia* y otros genomas
  reducidos; sin cajas DnaA diagnósticas, oriC por **GC skew**
  ([Gil et al. 2003](https://www.pnas.org/doi/10.1073/pnas.1533499100)).
- **Cianobacterias**: la dependencia de DnaA varía; *Synechococcus elongatus*
  puede replicar sin DnaA ([Ohbayashi et al.](https://kaken.nii.ac.jp/report/KAKENHI-PROJECT-21K15142/21K15142seika/)).
- **Mollicutes** (Mycoplasma/Ureaplasma): complejo de replicación simplificado
  ([MMBR 1998](https://journals.asm.org/doi/full/10.1128/mmbr.62.4.1094-1156.1998)).
- **Plástidos** (no bacterias) también carecen de dnaA.

La tabla `DNAA_LESS` es un **prior de interpretación**, no una rama
hardcodeada: el pipeline intenta los mismos métodos en todos los genomas y solo
cambia cómo se explica el resultado.

## Iterones: alcance y límite

Se implementó `iteron_scan` (k=12, ±10 kb, ≥4 copias). En Fusobacterium
(27% GC) los repetidos dominantes son AT-ricos de baja complejidad, no iterones
claros; por eso el estado queda `lineage_specific` con "iterones no detectados"
y el siguiente paso es un **HMM de RPA/Rep** o una base de iniciadores.

## Integración

- `webapp/origin_pipeline.py`: árbol y estados.
- `oriC_status.tsv`: estado, motivo, métodos, intentos.
- `oriC_all.gff3`: atributos `status`, `validated_by`, `reason`.
- GUI: pestaña **🧭 Origen** con posición, estado y explicación.

## Próximo

- HMM de DnaA/RPA/Rep para los 17 sin marcadores; validar iterones con HMM.
- Detección de iterones (RPA) en Fusobacteriaceae.
- Revisión automática de los 41 discordantes con GC skew.
