# Hoja de ruta hacia la versión 1.0

## Qué significa "1.0"

Un usuario instala el proyecto, descarga los recursos, sube un FASTA y obtiene:

1. anotación estructural con **% de cobertura**;
2. genes nuevos, desiertos y su evidencia etiquetada;
3. red metabólica y **brechas** priorizadas;
4. **origen de replicación** con calidad declarada;
5. **informe descargable** y trazable.

Y el equipo puede **reproducir los 400 genomas** con los mismos resultados.

## Estado actual

| componente | estado | evidencia |
|---|---|---|
| Anotación CDS (pyrodigal) | ✅ | F1 0,98 en benchmark |
| tRNA / rRNA | ✅ | F1 0,99 / 0,80 |
| Genes fantasma (T1) | ✅ | recall 0,99 · FPR 0,000 |
| Tier 2 (hypothetical) | ✅ | 53.650 candidatas · 1.164 clusters conservados |
| Tier 3 (modelos alternativos) | ✅ | 281.869 candidatos |
| Red metabólica dinámica | ✅ | DIAMOND→Rhea; reutiliza red de referencia |
| Brechas refinadas + soporte | ✅ | 8.655 brechas · soporte en 400 genomas |
| Catálogo KEGG | ✅ | 194 rutas · 5.361 reacciones |
| GUI (flujo, cobertura, informe) | ✅ | capturas v6–v8 |
| **ORCA / DoriC** | ⚠️ | **95/384 con confianza <0,5** |
| Móviles / defensa / CRISPR | ⚠️ | precalculado para 400, **no** en GFF ni GUI |
| Multidioma | ❌ | solo español |
| Estructura 3D | ⚠️ | Foldseek validado; ESMFold no integrado |
| Referencia NCBI para GCF externos | ⚠️ | solo GCF del proyecto |
| KEGG en redes subidas | ⚠️ | hoy solo genomas del proyecto |
| Reproducibilidad | ⚠️ | faltan versiones/hashes/manifiesto |
| Packaging | ⚠️ | setup/Docker/CI sí; falta `pyproject` y versión |
| Tests | ⚠️ | 41 unitarios; sin tests de GUI |
| Licencias (KEGG/Pfam) | ⚠️ | revisar términos de distribución |

## Vía de regulación L2/L3 (nueva)

- L2 ya tiene `genes.tsv`, `operons.tsv` y `promoters.tsv` para K-12 y un genoma
  de ejemplo.
- **Gate L2: PASA — 18/20 operones exactos** (mismos genes y orden), según
  `data/regulation/GATE.md`.
- L3 **desbloqueado**; en curso por `regulon-builder`.
- Detalle y checklist en `docs/REGULATION_PLAN.md` y `PLAN_REGULON_BUILDER.md`.
- Nota: este bloque reemplaza la afirmación antigua de “no pasa el gate” y
  “L3 bloqueado”.

## Bloqueantes para 1.0 (por prioridad)

### P0 — ORCA circular (credibilidad científica)

- 384 predicciones: **214 alta, 74 media, 95 baja, 1 sin dato**.
- **364 genomas son circulares**; **93 de las 95 bajas son circulares**.
- 67 de las bajas tienen DoriC disponible para anclar.
- Causa: ORCA asume linearidad; el origen en el borde del contig y el gate de
  marcadores penalizan la predicción.
- **Hallazgo nuevo:** inyectar `dnaA`+`dnaN` juntos puede arruinar la predicción
  (caso 0,02). Con **solo `dnaA`** sube a 0,42 y el error baja de 1,3 Mb a 19 kb.

**Plan:**
1. Correr ORCA en 4 configuraciones (ambos, dnaA, dnaN, ninguno) y elegir por
   confianza + consistencia con dnaA.
2. Rotar/duplicar flancos del contig circular y mapear de vuelta.
3. Consenso ORCA + DoriC + GC skew + densidad de cajas dnaA.
4. Relajar el gate de marcadores (dnaA, dnaN, gyrB, recA).
5. Marcar `quality` por genoma y "no resoluble" en ensamblajes fragmentados.
6. Reporte de calidad y delta a DoriC.

**Aceptación:** ≥95% en alta/media y delta mediano a DoriC <5 kb.

**Estado 2026-11-20:** 42 anómalos triados (23 skew-ORCA, 16 sin
convergencia, 3 solo GC skew); el fix circular completo sigue pendiente.

### P1 — Cobertura de genomas externos

- Descargar la anotación de referencia NCBI para un GCF dado (no solo los 400).
- Comparar cobertura y usar la referencia como plantilla de red.

### P1 — KEGG en redes subidas

- Llevar el catálogo KEGG a la red dinámica de un FASTA subido.
- Mostrar rutas con reacciones presentes/faltantes y soporte.

### P1 — Estados en todos los pasos

- Adoptar `ghost/status.py` (validated / consensus / lineage_specific /
  anomalous / unresolved / novel) en brechas, cierres, desiertos y ARN.
- Cada salida declara `status_validated_by`, `status_reason`, `status_attempted`.
- Detalle y prioridades en `PIPELINE_LOGIC.md`.

### P1 — Higiene de ARN

- Corregir los duplicados de `rna_annotations.json` (400/400 genomas).
- Reportar cobertura 20/20 de tRNA y marcar `undetermined`.

### P1 — Completitud de la anotación (no hay ADN basura)

Premisa del proyecto: después de miles de millones de años de evolución no hay
ADN sin función; si no conocemos su función, es **desconocida**, no basura.

- Integrar **móviles, defensa, CRISPR y virulencia** (hoy en JSON para 400) al
  GFF, al JSON y a la GUI.
- **Estado 2026-11-20:** JSON/GUI/TSV hechos; falta tipado Cas/IS/AMR y GFF por
  contig. Delegado a `mobile-builder`.
- Para genomas subidos: paneles HMM de defensa/CRISPR + detección de móviles.
- Explorar sistemáticamente todo desierto y toda región intergénica; etiquetar
  lo desconocido como `candidate_novel`, nunca descartarlo.

### P1 — Multidioma

- Español/inglés con selector en la GUI y mensajes de CLI.
- Textos, informes y documentación traducidos.
- Estructura i18n (gettext o diccionarios) desde el inicio para no reescribir.

### P1 — Reproducibilidad

- Fijar versiones (pyrodigal, DIAMOND, Pfam, Rhea, KEGG).
- Hashes de las bases y manifiesto por corrida.
- Semillas y parámetros en un `pipeline.yaml`.

### P2 — Packaging y distribución

- `pyproject.toml`, entry points (`dna-annotator`), `--version`.
- Docker verificado y CI que corra tests + smoke de GUI.
- Release notes y licencia.
- **Estado 2026-11-20:** `pyproject.toml`, `setup.py`, entrypoint
  `dna-annotate` y CI hechos; faltan wheel autónomo, Docker, i18n y licencia.

### P2 — Estructura 3D (opcional, no bloquea 1.0)

- Lo actual (dominios, homología, ESM-2) alcanza para priorizar.
- Añadir plegado **solo para el top N** con ESMFold (API/GPU) y Foldseek para
  comparar; no plegar todo el genoma.

### P2 — ARN rápido

- tRNAscan-SE y barrnap como primera pasada; cmsearch solo para confirmar o
  regiones sin hit.
- Caché de las 400 anotaciones como referencia rápida; no sustituye la
  anotación de un genoma nuevo.

### P2 — Validación y documentación

- Tests de la GUI (Playwright) en CI.
- Tutorial de 10 minutos y dataset de ejemplo.
- Página de límites y términos de uso.

## Lo que ya está listo

- Motor de anotación y descubrimiento validado.
- Catálogo de genes nuevos y desiertos con etiquetas de procedencia.
- Red metabólica dinámica y brechas con soporte en 400 genomas.
- Catálogo KEGG y rutas curadas con vista acotada.
- Flujo de usuario con cobertura, informe y campana.
- 41 pruebas unitarias y setup de recursos para GitHub.

## Recomendación

1. **Cerrar ORCA circular** (P0): es lo único que impide decir "toda bacteria
   tiene su origen con calidad declarada".
2. KEGG en redes subidas y referencia NCBI (P1).
3. Reproducibilidad y packaging (P1/P2).
4. Tests de GUI y documentación (P2).
