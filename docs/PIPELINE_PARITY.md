# Paridad: los 400 genomas vs un FASTA nuevo

Detalle de oriC: [`ORIC.md`](ORIC.md).

**Fecha:** 2026-11-23. **Pregunta:** ¿todo lo que tienen los 400 se puede
reproducir en un genoma nuevo?

## Inventario real

| capacidad | evidencia en los 400 | pipeline de FASTA nuevo | estado |
|---|---|---|---|
| CDS, tRNA, rRNA, GFF3, proteínas | 400/400 | `annotate` + flags | ✅ |
| ORFs no anotados y desiertos | 400/400 | `ghost` + desiertos | ✅ |
| Anotación funcional Swiss-Prot | 400 `functional.tsv` | ahora automática | ✅ nuevo |
| Caminos A/B/C | 400/400 `protein_paths.tsv` | `protein_paths.tsv` | ✅ clasificación; B/C estructural parcial |
| Complejos por vecindad/ortología | 397/400 firmas; 400/400 transferencia; 281.090 clusters del modelo | firma + modelo K-12 + transferencia | 🟡 candidatos, sin predicción cuaternaria |
| Red metabólica, reacciones y brechas | 400 redes | automática si hay DB | ✅ condicionada a datos |
| oriC ORCA | 400 filas de estado | `webapp/oric.py` | ✅ nuevo |
| DoriC | 306 genomas con verdad | lookup con `--accession` | 🟡 opcional |
| Operones/promotores/regulación | 402 directorios | no automático | ❌ |
| Móviles/defensa | 10 pilotos | módulo de resumen | 🟡 fuera de los 400 |
| Membrana/transporte | 400 TSV `full/` + 10 GFF piloto | `webapp/membrane_layer.py` | ✅ automático si hmmscan + panel |
| Taxonomía/ANI | 400 JSON `full/` | no automático | 🟡 datos 400; no pipeline nuevo |
| TSS/TTS/horquillas | precalculado | no automático | ❌ |
| Estructura/filtros | subconjunto | flags opcionales | 🟡 |
| Targets/allosteria | precalculado | no automático | 🟡 |

## Tasa de éxito de oriC en los 400

- 400 genomas con estado: **271 validated**, **77 consensus**,
  **10 lineage_specific**, **42 anomalous**.
- 306 genomas tienen verdad DoriC; de los 67 de baja confianza con DoriC,
  **63 quedaron a ≤50 kb** (mediana 0).
- De los 95 de baja confianza, **94/95 resueltos** al combinar ORCA, dnaA,
  GC skew y cajas; queda **1 en revisión**.
- `validated` no es una medida independiente de precisión: muchos casos usan
  DoriC como ancla. La evaluación honesta separa DoriC, consenso y skew.

## Validación en un genoma de los 400

`GCF_000006985.1_ASM698v1` (2.154.946 pb) procesado como FASTA nuevo:

- 1.301 proteínas con hit Swiss-Prot; caminos A/B/C = 1.100/201/708;
- ORCA + DoriC + dnaA → **validated**, posición **32.396**;
- la tabla precalculada de los 400 también tiene **32.396**, `validated`;
- cobertura de anotación 87,4 %.

Esto demuestra paridad para el núcleo anotación + función + oriC. No demuestra
paridad para regulación/operones, que sigue pendiente.

## Lo que ya cierra la brecha principal

El pipeline nuevo ahora produce:

```text
*.gff3
*.proteins.faa
*.functional.tsv
*.functional_annotation.json
*.protein_paths.tsv      # A / B / C
*.oric.json
*.oric.gff3
*.summary.json
*.report.md
```

Y ejecuta automáticamente, si existen las dependencias:

1. anotación CDS/tRNA/rRNA;
2. ghost/desiertos;
3. DIAMOND/Swiss-Prot;
4. clasificación A/B/C;
5. metabolismo si no se usa `--no-metabolic`;
6. ORCA + GC skew + cajas dnaA;
7. DoriC si se pasa `--accession`.

## Lo que sigue faltando para paridad total

### 1. Operones y regulación

Los 400 tienen `data/regulation/<genoma>/` con tablas de genes, operones,
promotores y evidencia. El pipeline nuevo no las construye. Implementarlo exige:

- tabla de genes y sintenía;
- modelos PWM y reglas;
- datos curados (RegulonDB/DBTBS) y sus licencias;
- benchmark contra RegulonDB.

No es un flag: es una etapa nueva.

### 2. Móviles, membrana y taxonomía

Membrana y taxonomía ya tienen salidas `full` para 400 genomas (TSV/JSON), pero
solo 10 GFF piloto de membrana y **no se ejecutan automáticamente** en un FASTA
nuevo. Móviles sigue siendo piloto. No son paridad de pipeline nuevo todavía.

### 3. TSS/TTS/horquillas

Requieren datos precalculados y modelos de ARN. No se ejecutan en un FASTA
arbitrario sin un pipeline de TSS/TTS.

### 4. Estructura y alosteria

Dependen de PDB/AFDB/ESMFold y, para dirección de efectores, GPU y datos
experimentales. No son parte del núcleo CPU.

## Conclusión honesta

No, todavía **no hay paridad total**. Sí hay paridad en el núcleo que define el
producto: anotación, proteínas, función, metabolismo, brechas y oriC.

La brecha más grande es **regulación/operones**. La segunda es **estructura
masiva**. Membrana ya se ejecuta en un FASTA nuevo si hmmscan y el panel Pfam
están instalados; taxonomía sigue precalculada. Móviles/defensa sigue siendo
piloto de 10 genomas.

## Orden para cerrar la paridad

1. Estabilizar functional + A/B/C + ORCA en 3–5 genomas públicos.
2. Implementar el módulo `oric` en la GUI con entrada FASTA.
3. Llevar regulación/operones a un comando reproducible.
4. Integrar membrana/taxonomía full 400 en el pipeline nuevo y ampliar
   móviles/defensa a los 400.
5. Estructura selectiva y GPU.
