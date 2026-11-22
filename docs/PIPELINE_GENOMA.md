# Pipeline de genoma crudo → represores (prototipo)

> Entrada: un FASTA de genoma sin anotar.
> Salida: "este gen se convierte en esta proteína, es un represor de esta familia, se une
> a estas zonas y responde a esto (o a aquello)" — con el estado de cada afirmación.
>
> **Es un prototipo.** Todo lo que dice es una **hipótesis con `status`**, no un hecho.
> Las cuatro limitaciones grandes están al final, y hay que leerlas antes de usar la salida.

---

## 1 · Flujo

```
FASTA crudo
   │
   ├─ 1. pyrodigal ──────────────► genes + proteínas
   │
   ├─ 2. hmmsearch (299 HMM) ────► candidatos con dominio de TF + familia
   │
   ├─ 3. referencia (RegulonDB) ─► homólogo más cercano por k-mers
   │                                └─ efector heredado si similitud ≥ 0,35
   │
   ├─ 4. PWM de la familia ──────► sitios en la región −250 pb
   │                                └─ percentil frente a 1.000-2.000 ventanas del genoma
   │
   └─ 5. informe ────────────────► TSV + Markdown + JSON
```

## 2 · Uso

**CLI**

```bash
python scripts/allostery/18_genome_pipeline.py genoma.fna \
    --out data/allostery/genome_run --max-tfs 40 --find-repressors
```

**GUI**

```bash
streamlit run scripts/allostery/app_alloreg.py
```

Sube el FASTA, marca **"Buscar solo posibles represores"** (activado por defecto), ajusta
el máximo de candidatos y pulsa *Ejecutar*. La GUI muestra tiempo total y RAM pico, y
despliega los tiempos por etapa. Para probar, un fragmento de 100-200 kb tarda segundos.

## 3 · Salidas

| Archivo | Contenido |
|---|---|
| `proteins.faa` | proteínas llamadas por pyrodigal |
| `tf.domtbl` | hits de HMM de dominio de TF |
| `represores.tsv` | una fila por candidato: gen, familia, homólogo, identidad, modo, efector, `status`, sitios |
| `summary.json` | conteos y familias |
| `informe.md` | tabla legible |

Ejemplo real (fragmento de 200 kb de *E. coli*):

| gen | familia | homólogo | sim. | modo | efector | status |
|---|---|---|---:|---|---|---|
| slice_200kb_61 | AraC | araC | 1,00 | inductor | D-fucose | inherited_homolog |
| slice_200kb_107 | GntR | pdhR | 1,00 | inductor | pyruvate | inherited_homolog |
| slice_200kb_76 | LacI/GalR | cra | 1,00 | correpresor | β-D-fructofuranosa-1-P | inherited_homolog |
| slice_200kb_73 | LysR | leuO | 1,00 | correpresor | 2-isopropilmalato | inherited_homolog |
| slice_200kb_18 | LysR | nhaR | 0,99 | correpresor | Na⁺ | inherited_homolog |
| slice_200kb_20 | HTH_Tnp_IS1 | insA1 | 0,99 | desconocido | — | unknown |

## 4 · Estados (`status`)

| status | Significado | Qué se puede afirmar |
|---|---|---|
| `inherited_homolog` | el homólogo de **RegulonDB** tiene efector anotado y cumple **identidad ≥ 40 %, cobertura ≥ 70 %, e ≤ 1e-5** | "probablemente responde a X, por homología" |
| `family_prior` (con `diamond_sprot`) | no hay efector heredable, pero sí homólogo nombrado en Swiss-Prot | "es la TF X de *Bacillus*, se desconoce su efector" |
| `family_prior` | no hay homólogo con efector; se usa el efector más común de la familia | "la familia suele responder a X" |
| `unknown` | ni homólogo con efector ni prior de familia | nada |

Ningún estado es `confirmed`. Para eso haría falta un ensayo.

## 5 · Métodos, en una línea cada uno

- **Genes**: **reutiliza** `webapp.annotator.annotate` del repo (no duplica la llamada de
  genes ni el GFF). Si el anotador cambia, este pipeline lo hereda.
- **TF**: hmmsearch contra un panel de 299 HMM de dominios de unión a ADN, umbral
  E ≤ 1e-5. Se toma el mejor dominio por proteína.
- **Familia**: mapa Pfam → familia canónica (`LacI/GalR`, `LysR`, `GntR`, `TetR`, `MarR`,
  `AraC`, `DeoR`, `ArsR/SmtB`, `MerR`, `Fur`, `IclR`, `RpiR`, `AsnC`, `PadR`, `Rrf2`…).
- **Homólogo**, en dos pasos:
  1. **DIAMOND** contra los 243 TFs de RegulonDB (DB de 65 KB, `ref_tfs.dmnd`). Si
     identidad ≥ 40 %, cobertura ≥ 70 % y e ≤ 1e-5 → `homolog_metodo=diamond` y el
     efector se puede heredar.
  2. Si no hay hit, **DIAMOND contra 8.300 TFs de Swiss-Prot** filtrados
     (`19_build_sprot_tf_db.py`, DB de 3,9 MB) → `homolog_metodo=diamond_sprot`, con
     nombre y organismo del hit, **sin efector** (Swiss-Prot no lo anota).
  3. Solo si el binario falta se usa Jaccard de 4-mers, y por debajo de 0,25 se reporta
     `sin_homologo` en lugar de un falso homólogo.
- **Sitios**: PWM de 14 pb (núcleo de máxima información) de los TFs **de la misma
  familia**; se escanea la región −250 pb de cada gen. Percentil frente a ventanas
  barajadas del propio genoma.
- **Efector**: herencia desde el homólogo (capa 1 del módulo, 955 pares validados).

## 5.0 · Caché

El pipeline guarda el resultado por **hash del FASTA + parámetros** en
`data/allostery/work/cache/<hash>/`. Repetir el mismo genoma baja de ~18 s a **0,5 s**.

## 5.1 · Tiempo y recursos medidos

| Genoma | Genes | Candidatos TF | Tiempo total | RAM pico | Dominio |
|---|---:|---:|---:|---:|---|
| 200 kb | 170 | 8 | **17,9 s** | 128 MB | 1 contig |
| 1 Mb | 912 | 64 | **37,5 s** | 129 MB | 1 contig |
| 2,05 Mb (*F. nucleatum*) | 1.911 | 48 | **43,4 s** | 134 MB | 1 cromosoma |
| 4,6 Mb (estimado) | ~4.200 | ~300 | **~2-3 min** | ~150 MB | 1 cromosoma |

Desglose en el caso de 1 Mb (con DIAMOND): genes 13,6 s · dominios 5,2 s ·
referencia + DB 3,3 s · asignación de homólogo/efector/sitios 15,5 s.

Recursos: **1-2 núcleos** (hmmsearch con `--cpu 2`), **sin GPU**, **sin red** en tiempo de
ejecución (la referencia y los PWM son locales). Disco de salida: ~0,3 MB por Mb de genoma.

Escalado: lineal con el número de genes; el cuello es la asignación (un candidato × hasta
8 PWM de su familia × 600 ventanas de fondo).

## 5.2 · Efectores multi-organismo: cobertura real de nuestros 400

Se aplicaron las dos vías pedidas sobre los 400 organismos del estudio.

| Vía | Fuente | Resultado |
|---|---|---|
| **1. Scrapeo RegPrecise** | `genome.jsp` de los organismos presentes allí | **19/400** organismos · **325 pares** TF–efector · **102 TFs** distintos |
| **2. Minado de texto UniProt** | comentarios de función de Swiss-Prot | **12 %** de las especies tienen ≥1 entrada revisada · **1 par** tras el filtro de precisión |

Datos duros del corpus: de 25 especies muestreadas, solo **3** tienen alguna entrada
revisada en Swiss-Prot. Nuestros 400 organismos son deliberadamente no-modelo.

El mecanismo funciona: en *F. nucleatum*, un candidato con homólogo **FruR** hereda
**fructosa-6-fosfato** (`status=regprecise_tf_name`). Pero la cobertura es baja por
escasez de datos, no por el código.

**Conclusión:** el cuello de la capa de efectores no es el pipeline, es que la anotación de
efectores está concentrada en organismos modelo. Para este corpus, `unknown` es la
respuesta honesta en la mayoría de los casos.

## 5.3 · Validación en genoma conocido y campaña

**Validación en *E. coli* K-12** (4,64 Mb, 2 min 1 s), contra los 243 TFs de RegulonDB
(comparando por gen, no por nombre):

| Estrategia | Recall | Candidatos |
|---|---:|---:|
| Solo dominio HMM (de novo) | **63 %** (153/243) | 261 |
| HMM ∪ homólogo DIAMOND | **92 %** (223/243) | 337 |

El panel HMM cubre ~63 % de las familias DBD; el resto se recupera por homología. En un
genoma sin parientes cercanos, el techo *de novo* es el 63 %.

**Campaña en 5 genomas desconocidos** (2,0-2,1 Mb cada uno, `22_campaign.py`):

| genoma | genes | candidatos TF | represores | tiempo |
|---|---:|---:|---:|---:|
| GCF_019552005.1 | 1.911 | 48 | 20 | 24,0 s |
| GCF_000020985.1 | 2.013 | 42 | 19 | 22,4 s |
| GCF_000091545.1 | 2.227 | 53 | 14 | 22,9 s |
| GCF_014337255.1 | 1.950 | 67 | 17 | 19,6 s |
| GCF_041200005.1 | 2.197 | 52 | 15 | 22,0 s |

**85 represores en 5 genomas, 22,2 s por genoma.** Familias dominantes: PadR 19, GntR 15,
HxlR 12, TetR 11, Rrf2 7, LacI/GalR 7, ArsR/SmtB 6, MarR 6.

Cobertura de efector en estos genomas: **0** — son organismos no modelo, como se documenta
en §5.2.

## 5.4 · Niveles de anotación: barato primero, caro bajo demanda

| Nivel | Qué añade | Coste medido (E. coli K-12, 4,6 Mb) | Quién lo paga |
|---|---|---|---|
| `--tier rapida` | dominios, familia, homólogo, efector | **55 s** | CPU local |
| `--tier estandar` | + sitios PWM en −250 pb | **66-84 s** | CPU local |
| *profunda* (bajo demanda) | conservación entre genomas (E2), co-folding (E4) | minutos-horas + GPU | Modal |
| *energía* (bajo demanda) | MM/GBSA → ABFE del efector | GPU o créditos | NeuroSnap / Modal |

El pipeline **solo corre los dos primeros niveles por defecto** y declara en
`summary.json` qué queda pendiente (`profunda_bajo_demanda`). Así una campaña de 400
genomas cuesta ~2,5 h de CPU y cualquier profundización se decide caso por caso.

Con `--find-repressors` no se recorta antes de filtrar (antes el tope sesgaba hacia
activadores de familias grandes): en *E. coli* K-12 el resultado pasa de 6 a **105
represores** con los dos niveles.

## 5.5 · Insumos del piloto de energía libre

`23_pilot_inputs.py` dejó preparados en `data/allostery/pilot/`:

| Sistema | Estructura holo | Ligando | Estado |
|---|---|---|---|
| LacI / IPTG | 1LBH | IPT | **listo**: química del CCD + pose cristalográfica |
| PurR / hipoxantina | 1PNR | HPA | **listo**: química del CCD + pose cristalográfica |
| TrpR / Trp | 1TRO | TRP | receptor listo; ligando **pendiente de alinear** (el cristal no tiene el carboxilato) |

Ver `data/allostery/pilot/README.md` para qué falta antes de enviar (protonación a pH 7,
unidad funcional, minimización) y los comandos de OpenFE / NeuroSnap.

## 5.6 · Preparación local de estructuras (OpenMM, sin GPU)

`25_prepare_local.py` protona a pH 7 y minimiza con OpenMM **en el equipo local**.
Integrado como segunda página de la GUI (`app_alloreg.py` → *Preparar estructura*).

```bash
python scripts/allostery/25_prepare_local.py receptor.pdb --out preparado.pdb \
    --ph 7.0 --max-iter 100
```

| Medición en este equipo (8 hilos, sin CUDA) | Valor |
|---|---|
| Cargar PDB + PDBFixer (PurR, 5.832 átomos) | **10,9 s** |
| Minimizar 20 pasos | 37 s (**1,9 s/paso**) |
| Minimizar 500 pasos | 455 s (7,6 min) |
| **Dinámica MD** (implícito, 5.832 átomos) | **0,25 ns/día** |

La consecuencia: localmente se puede **preparar**, no simular. Un ns de dinámica tarda
4 días; 10 ns, 40. Para dinámica hace falta GPU (Modal A10G: 20-50 ns/día) o créditos
de NeuroSnap.

**Ligandos**: la energía de interacción con un ligando requiere parámetros GAFF/OpenFF
(`openmmforcefields` + `openff-toolkit`), que no están instalados. El flag `--interaction`
existe y devuelve ese motivo en vez de fallar en silencio. La preparación del **receptor**
no depende de eso.

Verificado sobre los dos casos del piloto: LacI/IPTG (17.964 átomos con H) y
PurR/hipoxantina (5.832), sin contactos ligando–receptor < 2 Å tras minimizar.

## 5.7 · Puntuación de contacto (local, $0) en lugar de energía MM

**Bloqueo encontrado:** `openff-toolkit` y `openmmforcefields` no están disponibles en el
índice de paquetes de este entorno (Python 3.10; solo aparece `openmmforcefields 0.15.1`,
que exige ≥3.11). Sin parámetros de ligando no se puede calcular una energía MM del
complejo.

**Alternativa implementada** (`26_contact_score.py`): puntuación geométrica de contacto,
sin campo de fuerzas.

| Métrica | Qué mide |
|---|---|
| `contactos` | átomos pesados del ligando con algún átomo del receptor < 4,5 Å |
| `fraccion_enterrada` | fracción de átomos con ≥3 contactos |
| `dist_min_A` | choques (< 2,0 Å = mal) |
| `puentes_H_plausibles` | pares N/O a 2,5-3,5 Å |
| `complementariedad` | (contactos + 0,5·enterrados) / átomos |

**Validación en los dos pilotos** (pose cristalográfica vs la misma molécula desplazada
25 Å):

| Sistema | complementariedad pose | complementariedad decoy | separa |
|---|---:|---:|---|
| PurR/hipoxantina | **1,50** (10/10 contactos, 3 H) | 0,00 | sí |
| LacI/IPTG | **1,37** (15/15, 11 enterrados, 6 H) | 0,00 | sí |

La métrica separa una pose real de una desplazada, que es lo mínimo exigible. **No dice
cuánto une**: para eso hace falta energía (GAFF/OpenFF + MM/GBSA, o ABFE en NeuroSnap).

Integrado en la GUI (*Preparar estructura* → subir SDF) y como CLI:
```bash
python scripts/allostery/26_contact_score.py receptor_preparado.pdb ligand.sdf --decoy
```

## 5.8 · Entorno de energía MM local (Python 3.11 + OpenFF)

El sistema tiene Python 3.11 (`/usr/bin/python3.11`), pero **pip no puede instalar
`openff-toolkit`**: el índice solo ofrece 0.18.0 (yanked) y 0.19.0 (exige ≥3.12).
`openmmforcefields` también está yanked.

**Solución aplicada:** micromamba + conda-forge con **Python 3.12**, que sí tiene
disponible `openff-toolkit` 0.19.0 (con 3.11 el índice solo ofrece versiones yanked).
HOME y caché van dentro del workspace porque el sandbox no deja escribir en `~/.conda`:

```bash
export HOME=$PWD/.fakehome XDG_CACHE_HOME=$PWD/.cache_home \
       MAMBA_ROOT_PREFIX=$PWD/.mamba_root CONDA_PKGS_DIRS=$PWD/.mamba_pkgs
.tools/bin/micromamba -r $MAMBA_ROOT_PREFIX create -y \
    -p .mamba_root/envs/.mamba_mm312 -c conda-forge \
    python=3.12 openmm openmmforcefields openff-toolkit rdkit numpy ambertools
```

| Resultado | Valor |
|---|---|
| Python | **3.12.16** |
| OpenMM | 8.6.1 |
| openmmforcefields | 0.16.0 |
| openff-toolkit | **0.19.0** |
| AmberTools (AM1-BCC) | instalado |
| Disco | **3,8 GB** |
| Tiempo de instalación | ~15 min |
| Coste | **$0** |

`antechamber` debe estar en `PATH` para que OpenFF asigne cargas AM1-BCC; por eso el
script se lanza con `PATH=.mamba_root/envs/.mamba_mm312/bin:$PATH`.

Los entornos anteriores de 3.11 (venv y conda) se **borraron** tras verificar que el de
3.12 reproduce los mismos números (PurR: −94,5 kJ/mol de LJ, +318,8 de Coulomb).

### Energía de interacción MM (`27_interaction_mm.py`)

```
.mamba_root/envs/.mamba_mm312/bin/python scripts/allostery/27_interaction_mm.py \
    receptor_preparado.pdb ligand.sdf --steps 20
```

| Sistema | E interacción | LJ | Coulomb | t construcción |
|---|---:|---:|---:|---:|
| PurR/hipoxantina | **+224 kJ/mol** | −95 | +319 | 12-15 s |
| LacI/IPTG | **+146 kJ/mol** | −61 | +207 | 147 s |

**Aviso importante:** está calculada **en vacío**. El término de Coulomb no está
apantallado, así que el valor absoluto y el orden no son interpretables como afinidad.
Es un paso intermedio: falta solvente implícito (MM/GBSA), que requiere tratar la
energía de Born (no descomponible por pares) y por eso no se hace con la suma directa.

## 6 · Limitaciones (leer antes de usar)

1. **Homología**: con DIAMOND son identidad y cobertura alineadas, con umbrales
   estándar. El fallback por k-mers solo se usa si el binario falta, y se marca en
   `homolog_metodo`. Lo que sigue pendiente es **calibrar** esos umbrales con genomas
   curados: 40 %/70 % es la convención habitual, no una validación propia.
2. **Afinidad relativa, no kcal/mol.** Los percentiles indican "qué tan parecido es el
   sitio al motivo", no energía de unión. Ver `PLAN_ALOSTERIA.md` §8.9: no hay fuente
   pública de ΔΔG por TF sin rehacer el pipeline experimental de Reg-Seq.
3. **El efector es una hipótesis por homología, sin magnitud.** No dice si el ligando
   aprieta o suelta con qué fuerza. El detector de bisagra dice si el TF *puede* responder
   (especificidad 5/8), no cuánto.
4. **La llamada de genes en contigs cortos es frágil.** pyrodigal en modo meta funciona
   en fragmentos, pero los genes truncados en los bordes se pierden.

Limitaciones menores: no se modelan ARN reguladores ni sigma factors; los sitios se
buscan solo en −250 pb (no en regiones distales); los PWM cubren 82 TFs de *E. coli*, y
un TF de otra familia cae en el fallback (todos los PWM) — señalado en la salida.

## 7 · Qué haría falta para que deje de ser prototipo

| Pieza | Mejora | Coste |
|---|---|---|
| Homología | MMseqs2 o BLAST contra UniProt/RefSeq | horas |
| Sitios | PWM por familia aprendido de RegPrecise/CollecTF, no de *E. coli* | días |
| Afinidad | energía medida (Reg-Seq bien procesado) o experimental | semanas |
| Efector | magnitud: requiere ΔΔG; hoy solo dirección probable | fuera de alcance |
| Validación | correr sobre genomas con anotación curada y medir F1 por columna | días |

## 8 · Archivos

```
scripts/allostery/18_genome_pipeline.py   pipeline
scripts/allostery/app_alloreg.py          GUI Streamlit
data/allostery/genome_run/                   salida del ejemplo (200 kb)
data/allostery/work/tf_domains.hmm           panel de 299 HMM de TF
data/regulondb/gt/tf_pwms.tsv                82 PWM de RegulonDB
```

## 9 · Honestidad de la salida

El pipeline **no** dice "este gen es un represor que responde a xilosa". Dice:
"esta ORF tiene un dominio HTH de la familia XylR, su homólogo más cercano es XylR de
*E. coli* con similitud 0,98, y hereda su efector (xilosa) con estado
`inherited_homolog`; en su región −250 pb hay dos sitios con percentil 0,99 frente al
genoma".

La diferencia entre esas dos frases es toda la diferencia entre un resultado y una
promesa.
