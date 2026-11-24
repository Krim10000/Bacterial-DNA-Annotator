# Plan de trabajo — 12 semanas

**Fecha base:** 2026-11-24  
**Objetivo:** convertir el anotador en el mejor **CPU-first, multi-genoma,
explicable y honesto** para bacterias, con benchmark externo y release
reproducible.

**Restricción:** no hay GPU. El release 1.x debe cerrarse **sin GPU**. Las tareas
GPU pasan a un backlog separado y no bloquean la release.

## 1. Posición

No competir en gene calling con Bakta/Prokka. Ganar en:

- integración local CPU;
- metabolismo y brechas;
- regulación/operones;
- complejos y estructura;
- evidencia/confianza explícitas;
- GUI ES/EN/ZH.

## 2. Estado de partida

- 400 genomas con A/B/C (`protein_paths.tsv`).
- 281.090 clusters de complejos del modelo K-12.
- 397/400 firmas; 400/400 transferencia.
- Membrana 400 y automática en FASTA nuevo.
- Regulación por reglas en FASTA nuevo.
- 196 tests.
- Benchmark externo: `not_evaluable` (Bakta/Prokka/DFAST no instalados).
- Complejos: holdout independiente `not_evaluable` fuera de K-12.
- Estructura: visor 3D; Foldseek-Multimer sin DB de ensamblajes.
- Móviles/defensa: piloto.
- Taxonomía: 400 precomputada; no automática.

## 3. Principios

1. Ninguna etiqueta se convierte en `confirmed` sin experimento.
2. `not_evaluable` no es resultado negativo.
3. Un benchmark interno no es precisión de producción.
4. No entrenar y validar con el mismo gold.
5. No redistribuir datos restringidos.
6. Cada tarea deja test, doc y salida reproducible.
7. Sin GPU no se pliega masivamente ni se hace coevolución. El núcleo CPU debe
   ser útil por sí solo.
8. La GPU, si aparece, solo añade estructura/coevolución; nunca “cierra” función.

## 4. Líneas de trabajo

| ID | línea | prioridad | dueño |
|---|---|---|---|
| A | Benchmark externo y calibración | P0 | benchmark-layers |
| B | Paridad de anotación base | P0 | nuevo `annotation-base` |
| C | Función remota y bases curadas | P0 | functional-remote |
| D | Regulación/operones/TSS | P0 | regulation-auto |
| E | Complejos: holdout y priorización | P1 | cluster-calibration |
| F | Estructura CPU (AFDB/PDB/Foldseek local) | P1 | b-structural + root |
| G | Brechas: criterio de cierre | P1 | root + nuevo `gap-closure` |
| H | Móviles/defensa/AMR/BGC | P1 | nuevo `mobile-defense` |
| I | Taxonomía automática | P1 | root |
| J | Empaquetado y adopción | P1 | packaging |
| K | GUI: comparar/explicar | P1 | gui-compare + root |
| L | Datos/licencias/reproducibilidad | P0 | license-audit |

## 5. Calendario

### Semana 1 — Preparar la verdad

- [ ] Instalar Bakta, Prokka y DFAST en entorno controlado (Docker/Conda).
- [ ] Congelar genomas held-out (no K-12).
- [ ] Congelar versiones de DB (Swiss-Prot, Rhea, Pfam, ComplexPortal).
- [ ] Ejecutar `tools/benchmark_external.py --run`.
- [ ] Publicar `reports/external_benchmark.*`.
- [ ] CI con benchmark estructural, no solo tests.

**Entregable:** primer benchmark externo real.

### Semanas 2–3 — Paridad base

- [ ] Comparar CDS/tRNA/rRNA/% hypothetical/tiempo/RAM contra Bakta/Prokka/DFAST.
- [ ] Corregir diferencias de GFF3/atributos.
- [ ] Añadir ncRNA/small ORF/plásmidos al benchmark.
- [ ] Documentar qué no igualamos y por qué.

**Métrica:** CDS F1 ≥0,98; tRNA/rRNA F1 ≥0,95.

### Semanas 4–5 — Función remota

- [ ] Integrar UniRef/TrEMBL/InterPro/eggNOG locales.
- [ ] Medir cobertura de B/C y EC.
- [ ] Añadir `remote_status` al pipeline/GUI.
- [ ] No cambiar A/B/C automáticamente: solo evidencia adicional.

**Métrica:** cobertura B/C y EC medida; precisión por etiqueta.

### Semanas 6–7 — Regulación

- [ ] Benchmark ciego de operones contra RegulonDB.
- [ ] TSS/TTS y sigma factors.
- [ ] Promotores: definir emparejamiento TSS→ventana.
- [ ] Separar `rule_based` de `experimental`.

**Métrica:** F1 de operones publicada; TSS `not_evaluable` hasta tener gold.

### Semanas 8–9 — Complejos y estructura

- [ ] Conseguir ComplexPortal de otras especies.
- [ ] Holdout real por especie/complejo.
- [ ] Revisar P1 (207) y P2 (1.350).
- [ ] Foldseek local con DB de monómeros/ensamblajes disponibles.
- [ ] AFDB/PDB para B/C priorizados, sin plegar.
- [ ] Visor 3D y B-factor/RMSX.
- [ ] No prometer pLDDT/ipTM sin GPU.

**Métrica:** precisión de alta confianza ≥0,80 en holdout; cobertura estructural
CPU medida; no cerrar brechas.

### Semanas 10–11 — Cobertura y producto

- [ ] Taxonomía automática (normalizador local ya existe; falta ANI/red).
- [ ] Móviles/defensa/AMR: 400 TSV tipados; 400 FASTA (90 en `raw_excluded`).
- [ ] BGC: decidir integrar antiSMASH o no.
- [ ] Brechas: dirección, transporte, termodinámica.
- [ ] SBML borrador (opcional).
- [ ] Bioconda construida y Galaxy probado con Planemo.
- [ ] GUI: comparar genomas y explicar por capa.

### Semana 12 — Release y credibilidad

- [ ] Release con DOI.
- [ ] Benchmarks publicados.
- [ ] `bio.tools`, Bioconda, Galaxy.
- [ ] Post/JOSS/paper.
- [ ] Revisión de licencias y datos.

## 6. Tareas por prioridad

### P0

1. **Benchmark externo**
   - salida: reporte y CI.
   - criterio: 5–6 genomas, mismas entradas.

2. **Paridad base**
   - salida: tabla de diferencias y fixes.
   - criterio: F1 objetivo.

3. **Función remota**
   - salida: DBs locales y cobertura.
   - criterio: `remote_status` medido.

4. **Regulación**
   - salida: benchmark operones/TSS.
   - criterio: F1 ciego.

5. **Licencias/reproducibilidad**
   - salida: inventario y versiones.
   - criterio: release reproducible sin datos restringidos.

### P1

6. Holdout de complejos.
7. Estructura CPU: AFDB/PDB, Foldseek local, visor; GPU fuera de 1.x.
8. Brechas con criterio de cierre.
9. Móviles/defensa/AMR.
10. Taxonomía.
11. Bioconda/Galaxy.
12. GUI comparar/explicar.

### P2 — CPU

13. FBA/SBML.
14. BGC.
15. Rendimiento/offline extremo.
16. Flipbook/RMSX solo si hay trayectoria y R (extra opcional).

### Backlog GPU (no bloquea 1.x)

- ESMFold local.
- AlphaFold-Multimer/ColabFold.
- Boltz-2.
- Foldseek-Multimer con DB de ensamblajes.
- Coevolución con MSA.
- MD/RMSX a escala.

Se activan solo si hay GPU y datos; no se anuncian como parte de 1.x.

## 7. KPIs

| capa | KPI | objetivo |
|---|---|---|
| CDS | F1 vs Bakta/Prokka | ≥0,98 |
| tRNA/rRNA | F1 | ≥0,95 |
| Función | cobertura EC/UniProt | medida |
| A/B/C | calibración | por etiqueta |
| Operones | F1 ciego | publicar |
| oriC | mediana a DoriC | ≤5 kb |
| Membrana | acuerdo TMHMM | ≥0,90 |
| Complejos | precisión alta confianza | ≥0,80 holdout |
| Estructura CPU | cobertura AFDB/PDB | medida |
| Brechas | cierres con evidencia | 0 sin prueba |
| Tiempo 5 Mb | CPU | <5 min |
| Reproducibilidad | release + DOI | cada release |
| GPU | no requerida | 1.x CPU-only |

## 7b. Extra: dinámica (RMSX/Flipbook)

No es núcleo de anotación. Se añade como extra opcional para proteínas
priorizadas con trayectoria MD.

- D1 visor PDB multimodelo + B-factor;
- D2 adaptador `dynamics_layer.py`;
- D3 validación contra RMSF;
- D4 docs/licencias.

Sin R/MDAnalysis/trayectoria: `not_configured`. No cierra brechas.

## 8. Riesgos

| riesgo | mitigación |
|---|---|
| Bakta/Prokka no instalables | Docker/Conda; reportar `not_evaluable` |
| Datos restringidos | no redistribuir; versionar fuentes |
| Circularidad | holdout por especie/complejo |
| Sin GPU | no plegar; usar AFDB/PDB; GPU en backlog |
| Exceso de candidatos | P1/P2 y filtros de evidencia |
| Sobreafirmar | etiquetas y `closure=no` |
| Benchmarks caros | lotes y CI semanal |

## 9. Definición de “mejor”

No es tener más columnas. Es:

```text
precisión medida + cobertura + calibración + explicabilidad
+ reproducibilidad + adopción
```

## 9b. Definición de terminado CPU-only

Una release 1.x CPU-only está terminada cuando:

- benchmark externo ejecutado o `not_evaluable` documentado;
- A/B/C, membrana, regulación, complejos y oriC reproducibles en CPU;
- B/C tienen estado/dominio/contexto, no estructura masiva;
- complejos están priorizados y no se presentan como validados;
- packaging y GUI funcionan sin GPU;
- ningún resultado GPU es requisito;
- las limitaciones están en README/docs.

## 10. Regla final

```text
sin benchmark no hay “mejor”
sin holdout no hay precisión
sin experimento no hay “confirmed”
sin evidencia no hay cierre
```
