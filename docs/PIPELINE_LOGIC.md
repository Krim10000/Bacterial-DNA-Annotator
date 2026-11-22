# Lógica del pipeline: estados, evidencia e intentos

## Principios

1. **No hay ADN basura.** Después de miles de millones de años de evolución, una
   región sin función conocida es **desconocida**, no descartable. Todo desierto
   y toda región intergénica se exploran y se etiquetan.
2. **Evidencia primero.** Cada resultado declara métodos, evidencia, confianza y
   qué se intentó. Una etiqueta sin fuente no es válida.
3. **Generalidad.** Los métodos no dependen de un caso particular. Los linajes
   solo cambian la **interpretación** del resultado, no qué se intenta.
4. **Fallbacks.** Si el método principal falla, se prueban alternativas y se
   registra el intento; nunca se devuelve "vacío" sin explicación.
5. **Trazabilidad.** `label` dice qué es; `status` dice cómo se decidió.
6. **Reproducibilidad.** Versiones, parámetros y hashes quedan en el manifiesto.

## Vocabulario de estados

| estado | significado |
|---|---|
| `validated` | ≥2 métodos independientes coinciden. |
| `consensus` | Los métodos coinciden, sin verdad externa. |
| `lineage_specific` | Desviación esperada del linaje (p. ej. sin dnaA). |
| `anomalous` | Los métodos discrepan; se declara el motivo. |
| `unresolved` | Ningún método dio candidato; se listan los intentos. |
| `novel` | Sin referencia conocida; candidato exploratorio. |

Cada fila puede llevar: `status`, `status_validated_by`, `status_reason`,
`status_attempted`, `status_confidence`, `status_evidence` (módulo
`ghost/status.py`).

## Qué debe declarar cada paso

| paso | método principal | validación | estado esperado |
|---|---|---|---|
| CDS | pyrodigal | modelo + ORF completo + no solape | `validated` / `anomalous` (partial) |
| tRNA | tRNAscan-SE | isotipo + anticodón | `validated` / `lineage_specific` (sin aa) |
| rRNA | cmsearch (Rfam) | Rfam + longitud completa | `validated` / `anomalous` (partial/dup) |
| Genes fantasma (T1) | 6 marcos + scoring | codón/RBS/ESM/conservación | `novel` / `consensus` |
| Tier 2 | hypothetical/DUF | dominio + conservación | `consensus` / `novel` |
| Tier 3 | modelos alternativos | evidencia vs canónico | `anomalous` / `consensus` |
| Desiertos | ≥3 kb sin señales | ORFs conservados | `novel` / `consensus` |
| Red metabólica | DIAMOND→Rhea | reacciones + soporte | `validated` / `consensus` |
| Brechas | traza + demanda | soporte en referencia | `validated` / `anomalous` |
| Cierres | dominio/homología | dirección + reacción | `validated` / `anomalous` |
| Origen | ORCA + GC skew + DoriC | ≥2 métodos | ya implementado |
| Defensa/CRISPR | paneles/HMM | dominio + contexto | `validated` / `novel` |

## Aplicación por paso

### CDS
- Principal: pyrodigal.
- `validated`: CDS completo con inicio/stop y sin solape problemático.
- `anomalous`: parcial, solapante o sin start/stop; se intenta el motor de ORFs.
- `attempted`: pyrodigal; ORFs 6 marcos.

### tRNA
- tRNAscan-SE da `product=tRNA-Ala` (isotipo).
- `validated`: isotipo + anticodón.
- `lineage_specific`: genoma reducido sin un tRNA (p. ej. Glu/Asn por GatCAB).
- `anomalous`: `tRNA-Undet`/`tRNA-Sup`; revisar.
- Se reporta cobertura 20/20 y los ausentes.

### rRNA
- cmsearch (Rfam) como principal; barrnap como prefiltro rápido.
- `validated`: 16S/23S/5S con Rfam y longitud completa.
- `anomalous`: parcial o duplicado (hoy `rna_annotations.json` tiene duplicados).
- `attempted`: cmsearch; barrnap.

### Genes fantasma y tiers
- La evidencia ya existe (codón, RBS, composición, ESM, conservación).
- Añadir `status`: `novel` si no hay hit/dominio; `consensus` si hay evidencia
  múltiple; `anomalous` si las señales discrepan.
- `attempted`: 6 marcos; DIAMOND; Pfam; ESM.

### Desiertos
- `novel`: desierto sin ORFs conservados.
- `consensus`: desierto con ORFs conservados (candidato a gen/operón).
- `validated`: si además hay dominio/homología.
- `attempted`: ORFs, conservación, dominios.

### Red metabólica y brechas
- `validated`: reacción con enzima y demanda, soporte en referencia.
- `anomalous`: ión/transporte/ROS o sin demanda (se descarta con motivo).
- `novel`: brecha sin soporte en los 400 (p. ej. pigmento raro).
- `attempted`: DIAMOND, Rhea, soporte, filo.

### Cierres de brecha
- `validated`: dominio + reacción + dirección verificada.
- `anomalous`: dominio correcto pero dirección no verificada.
- `attempted`: dominio, homología, dirección.

### Origen de replicación
- Ya implementado: `validated`, `consensus`, `lineage_specific`, `anomalous`,
  `unresolved`.
- Ejemplo: Fusobacteriaceae sin dnaA → `lineage_specific` con RPA/iterones.

## Prioridad para adoptar el esquema

| prioridad | paso | por qué |
|---|---|---|
| P1 | Brechas y cierres | hoy hay soporte/confianza, falta motivo/intentos |
| P1 | Desiertos | falta distinguir novel/consensus/validated |
| P1 | tRNA/rRNA | falta 20/20, duplicados y undetermined |
| P2 | Genes fantasma/Tier 2/3 | la evidencia existe; falta estado |
| P2 | Defensa/CRISPR | paneles y contexto |
| P3 | Estructura 3D | opcional |

## Cómo extender un paso

1. Definir el método principal y al menos un fallback.
2. Definir la evidencia que se guarda.
3. Mapear el resultado a un estado del vocabulario.
4. Declarar `validated_by`, `reason`, `attempted`.
5. Añadir una prueba y documentar el árbol.

## Implementación

- `ghost/status.py`: `STATUSES`, `Decision`, `decide`, `apply_status`.
- `ghost/labels.py`: etiqueta + fuente + evidencia + confianza + interpretación.
- `webapp/origin_pipeline.py`: primer paso completo con el vocabulario.

## Referencias

- Fusobacteriaceae sin dnaA y RPA/iterones: Brezellec & Pasek 2026.
- Endosimbiontes y GC skew: Gil et al. 2003.
- Cianobacterias DnaA-independientes: Ohbayashi et al.
- Mollicutes: MMBR 1998.
