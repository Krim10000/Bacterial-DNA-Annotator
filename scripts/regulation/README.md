# Capa regulatoria (L2 + L3)

Convierte la anotación de un genoma bacteriano en **unidades de transcripción**
(L2) y en un **regulón** con reglas ejecutables en pseudo-código (L3). Es la
capa que faltaba entre "bolsa de CDS" y "programa": el `if` de la bacteria vive
en el promotor y en el operón, no en el gen.

Todo es **determinista** (orden estable, sin aleatoriedad), **sin red** y escribe
solo bajo `data/regulation/`.

## Submódulos fusionados

Esta carpeta absorbió tres módulos de la reorganización de 2026:

- [`README.operon_finder.md`](README.operon_finder.md) — buscador de operones
  determinista y terminadores.
- [`README.tss_tts.md`](README.tss_tts.md) — fronteras TSS/TTS.
- [`README.chains.md`](README.chains.md) — cadenas alternativas de reacción.

## Orden de ejecución

```bash
P=.venv/bin/python
$P scripts/regulation/10_build_gene_table.py --gbff data/test_ecoli.gbff   # K-12
$P scripts/regulation/11_predict_operons.py --build-synteny               # 1 vez
$P scripts/regulation/11_predict_operons.py --genome ecoli_k12
$P scripts/regulation/12_promoters_sigma.py --genome ecoli_k12
$P scripts/regulation/13_build_regulon.py  --genome ecoli_k12
$P scripts/regulation/13_build_regulon.py  --validate --genome ecoli_k12
$P scripts/regulation/14_validate.py                                     # gates + VALIDATION.md

# lote completo (400 genomas del set ab initio)
$P scripts/regulation/10_build_gene_table.py --all --with-seq
$P scripts/regulation/11_predict_operons.py --all
$P scripts/regulation/12_promoters_sigma.py --all
$P scripts/regulation/13_build_regulon.py  --all
$P scripts/regulation/14_validate.py
```

Cada script acepta `--selftest` (unittest implícito con `assert`, sin pytest).

## 10 — tabla génica (`<G>/genes.tsv`)

Una fila por gen. Es la interfaz estable hacia L2/L3.

| columna | contenido |
|---|---|
| `contig,start,end,strand` | coordenadas 1-based, hebra `1`/`-1` |
| `gene_id` | ID de pyrodigal (set de 400) o `protein_id`/`locus_tag` (GenBank) |
| `locus_tag` | b-number en K-12 (es el ID que usa RegulonDB) |
| `gene_name,product,ec` | nombre, producto, EC si existe |
| `acc,fam` | accesión y familia (stem) de UniProt/Swiss-Prot |
| `partial,rbs_motif,rbs_source,start_type,conf` | señales de traducción |
| `rna,gene_type` | `CDS`/`tRNA`/`rRNA`/`ncRNA` |
| `status,evidence,attempted,confidence` | decisión y procedencia |

`rbs_motif` viene de pyrodigal cuando existe; si no, se calcula
(`rbs_source=computed`, heurística Shine-Dalgarno simple).

## 11 — operones (`operons.tsv` + `operon_genes.tsv`)

Método por **combinación de evidencia**, nunca una señal sola:

1. **Distancia** con umbral calibrado por genoma: histograma de huecos
   (bins de 5 pb), suavizado gaussiano y **primer mínimo local** tras el modo
   pequeño → `t_lo`. El bulto secundario (operones laxos) da `t_mid`; `t_hi` es
   el techo de la zona ambigua. Nada de 50 pb fijos.
2. **Acoplamiento traduccional** (solapamiento o hueco ≤ 3 pb).
3. **RBS** del gen de abajo (ausente/débil ⇒ compatible con cotranscripción).
4. **Sintenía**: mismo par de ortólogos adyacente en otros genomas
   (`synteny_pairs.tsv`, construido con `--build-synteny`). Clave = accesión
   UniProt; reportado como proxy de ortología.
5. **Coherencia funcional** (prefijo génico o palabra de ruta; EC a 3 niveles).

Decisión en tres bandas:

| hueco | regla |
|---|---|
| `≤ t_lo` | se une |
| `(t_lo, t_mid]` | se une si `score ≥ 0.35` (basta una evidencia secundaria) |
| `(t_mid, t_hi]` | se une si `score ≥ 0.50` (sintenía + coherencia) |
| `> t_hi` | no se une |

`score = 0.45·sintenía(≥3) + 0.35·coherencia + 0.25·acoplamiento + 0.20·sin RBS`.

Restricciones duras: un gen en un solo operón, coordenadas monótonas, un ARN
intercalado o un cambio de hebra corta la unidad. **No** se usa la ausencia de
terminador (no hay predictor instalado).

> Los ARN **anidados** dentro de un CDS (p. ej. `xylZ` dentro de `xylA`) no
> cortan: son elementos internos del transcrito. Sin esta regla, operones reales
> se fragmentan.

## 12 — promotores y σ (`promoters.tsv`)

Escaneo posicional en 120 pb río arriba: `-35` TTGACA (≤2 cambios), `-10`
TATAAT (≤3), espaciador 15–19 pb, `-10` extendida (TG) con el espaciador
medido desde el TG. `score` = 6−mm35 + 6−mm10 + bonus espaciador + extendida
(máx. 15); umbral 9. **Es evidencia posicional, no un PWM**: `consensus` como
máximo. `sigma=sigma70` es el valor por defecto declarado (`sigma_status`).

`regulators_inventory.tsv` lista factores sigma por producto/nombre de gen.
La confirmación por Pfam **no** se pudo hacer: el panel local
(`data/db/pfam/panel_focus.hmm`, 18 modelos) no contiene PF00140/PF04542/
PF04539/PF04545/PF04997. Se declara `pfam_checked=false`, no se finge.

## 13 — regulón (`regulon.json` + `expression_rules.tsv`)

- Aristas `documented`: tabla curada TF→genes con cita (única fuente que permite
  `logic_status=documented`).
- Aristas `motif`: escaneo del motivo del TF en 300 pb río arriba, restringido a
  la **zona de arquitectura** del promotor (`ZONES`). Signo por posición
  (operador sobre el promotor ⇒ represión; río arriba del -35 ⇒ activación).
- Lógica por operón **asumida**:
  `ON IF (algún activador) AND NOT (algún represor)`.
- Motivos incluidos: CRP, FNR, Fur, LexA, LacI, TrpR, ArcA, IHF (cada uno con
  cita). Excluidos y declarados: H-NS, Lrp, NarL, OxyR, SoxS, Fis (sin consenso
  fiable) y RpoS (es factor sigma, va en 12).
- `effector` (cAMP, triptófano, Fe²⁺…) se toma de `data/regulondb/gt/gt_tfs.tsv`
  cuando existe: es la entrada del `if`.

## Límites conocidos

- **La distancia no basta**: el techo medido de pair-level F1 es ~0.79 con este
  conjunto de evidencia (el gate pedía 0.80). Ver `VALIDATION.md`.
- **La sintenía depende del set**: K-12 está representado por 4 genomas
  Escherichia en los 400, así que cubre ~10% de los enlaces.
- **Sin expresión**: no hay RNA-seq; nada de esto mide expresión, solo
  potencial regulatorio.
- **La lógica es asumida**, no medida. Un `if/else` "compilado" de aquí hereda
  esa suposición.
- **Motivos de TF**: 8 entradas curadas; P/R real medido en `VALIDATION.md`.
  No reconstruyen el regulón completo.

## Sprint 1 (2026-11-19): qué cambió

1. **ncRNA es miembro del operón, no frontera.** RegulonDB incluye `fimR2`,
   `xylZ`, `rbsZ` y `malH` dentro de sus operones; solo tRNA/rRNA cortan y un
   ncRNA que solapa el intervalo se absorbe como miembro.
2. **Coherencia ponderada por rareza**: `0.30 + 0.25·especificidad`. Un prefijo
   raro (`ara`, 9 genes) es evidencia fuerte; uno común, no.
3. **Veto suave a reguladores**: un TF adyacente suele tener promotor propio;
   en la banda lejana exige score ≥0.75.
4. **Techo de la zona ambigua** 168 → 300 pb (el `min()` lo recortaba).
5. **Promotores**: ventana 150 pb, espaciador 14–21 (distribución real de
   RegulonDB), -35 hasta 3 cambios, y **prior de distancia** (la -10 a 30–60 pb
   del inicio). Cobertura 24% → 65.9%.
6. **Sigma por HMM**: `data/db/pfam/panel_regulation.hmm` (PF00140/04542/04539/
   04545/04997) + `hmmsearch --cut_ga`. `validated` solo si dominio **y** nombre;
   el dominio sin nombre entra como `novel`. `Sigma70_r4` es un HTH compartido
   con reguladores de respuesta: solo no basta.
7. **`15_gate.py`** genera `GATE.md` con los 20 operones curados.
