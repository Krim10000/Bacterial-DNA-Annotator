# Homología remota B/C — auditoría CPU

## Alcance

"Remota" aquí significa **homología distante** (UniRef, TrEMBL, InterPro,
eggNOG, Pfam), no un servidor remoto. La auditoría es local y CPU-only:

- consulta solo bases ya presentes en `data_root`;
- **no usa red** ni descarga bases;
- audita **solo filas B/C**;
- añade columnas `remote_*` y escribe
  `<genoma>.remote_homology.tsv`;
- **no reclasifica A/B/C** y **no cierra función**.

La fila A no se toca. Un hit remoto es evidencia de homología/familia, no una
función confirmada. Un `no_hit` no prueba ausencia de homología.

## API

```python
from pathlib import Path
from webapp.remote_homology import available, run

info = available(Path("data"))
res = run(rows, proteins, outdir, "GCF_000005845.2", Path("data"),
          threads=4)
```

`rows` son las filas de `functional.annotate` (`protein_id`, `path`, etc.) y
`proteins` la lista con `id`/`aa`. `run` modifica en sitio solo las filas B/C.

`available(data_root)` devuelve un diccionario con las cinco fuentes:

```python
{
  "uniref":   {"available": True, "path": "...", "kind": "diamond",
               "engine": "diamond", "role": "sequence",
               "searchable": True, "note": ""},
  "trembl":   {...},
  "interpro": {...},
  "eggnog":   {...},
  "pfam":     {"available": True, "path": ".../Pfam-A.hmm",
               "kind": "hmmscan", "engine": "hmmscan",
               "role": "profile", "searchable": True, "note": ""},
}
```

| campo | significado |
|---|---|
| `available` | hay un archivo local de esa fuente |
| `path` | archivo detectado (o `""`) |
| `kind` | `diamond`, `hmmscan` o `data` |
| `engine` | motor previsto (`diamond`/`hmmscan`) o `""` |
| `role` | `sequence` (UniRef/TrEMBL/eggNOG) o `profile` (Pfam/InterPro) |
| `searchable` | hay motor, índice y, para HMM, `hmmpress` |
| `note` | motivo cuando no es consultable |

`available` no toca la red; solo lista archivos y comprueba `PATH`.

## Detección local

Se mira en `data_root/db`, `data_root` y `data_root/databases`, con un descenso
limitado a carpetas cuyo nombre contiene `uniref`, `trembl`, `interpro`,
`eggnog`, `pfam`, `db` o `database`. No se recorre todo `data_root`.

| fuente | se detecta como | motor | requisito para buscar |
|---|---|---|---|
| UniRef | `uniref*.dmnd` (también en subcarpeta) | DIAMOND | `.dmnd` + `diamond` en `PATH` |
| TrEMBL | `*trembl*.dmnd` | DIAMOND | `.dmnd` + `diamond` |
| InterPro | `interpro*.hmm` o `interpro*.dmnd` | hmmscan / DIAMOND | HMM presionado o `.dmnd` + motor |
| eggNOG | `eggnog*.dmnd`; `eggnog.db`/`.tsv` se detectan | DIAMOND | `.dmnd` + `diamond` |
| Pfam | `Pfam-A.hmm` o `panel*.hmm` en carpeta `pfam` | hmmscan | HMM presionado + `hmmscan` |

Notas:

- Un FASTA de UniRef/TrEMBL sin `.dmnd` se detecta pero no se consulta.
- Un `eggnog.db` (SQLite/anotaciones) se detecta pero no se consulta sin un
  `.dmnd` de proteínas eggNOG.
- Un HMM sin `.h3f/.h3i/.h3m/.h3p` se detecta pero no se consulta.
- Si falta el binario, la fuente queda como no consultable; no se lanza nada.

## Flujo de `run`

```text
filas
  → solo B/C
  → available(data_root)
  → ¿hay DB local?
       no → remote_status=not_available (no falla)
       sí → búsqueda de secuencia (una base, por prioridad)
              UniRef > TrEMBL > eggNOG
          → filas sin hit de secuencia → búsqueda de perfil
              Pfam > InterPro
  → escribir <genoma>.remote_homology.tsv
```

Decisiones conservadoras:

- se usa **una** base de secuencia y **una** de perfil por ejecución, no todas;
- el perfil solo se busca en las filas que siguen sin hit de secuencia;
- DIAMOND: `--evalue 1e-5`, un solo target por query;
- hmmscan: `--cut_ga` si el HMM declara umbrales GA; si no,
  `-E 1e-5 --domE 1e-5`;
- cualquier fallo del motor se registra en `warnings` y no tumba la anotación.

Comandos equivalentes (los construye el módulo):

```bash
diamond blastp -d data/db/uniref90.dmnd -q query.faa -o hits.m8 \
  --outfmt 6 qseqid sseqid pident length evalue bitscore stitle \
  --evalue 1e-5 --max-target-seqs 1 --threads 4

hmmscan --domtblout hits.domtbl --noali --cpu 4 --cut_ga \
  data/db/pfam/Pfam-A.hmm query.faa
```

## Columnas del TSV

`<genoma>.remote_homology.tsv`, separado por tabuladores, solo con filas B/C:

| columna | contenido |
|---|---|
| `genome` | genoma auditado |
| `protein_id` | id de la proteína |
| `path` | `B` o `C` (sin cambios) |
| `remote_status` | estado de la auditoría |
| `remote_source` | `uniref`, `trembl`, `interpro`, `eggnog` o `pfam`; vacío sin hit |
| `remote_hit` | accesión del hit local; vacío sin hit |
| `remote_evalue` | e-value del hit; vacío si no hay hit |
| `remote_note` | procedencia, interpretación y límites |

Las filas A no aparecen en el TSV.

## Estados

| `remote_status` | cuándo | interpretación |
|---|---|---|
| `hit` | búsqueda local con hit | homología/familia; **no** función |
| `no_hit` | se buscó y no hubo hit | no prueba ausencia de homología |
| `not_available` | no hay ninguna DB local | auditoría no ejecutada |
| `not_searched` | hay DB pero sin motor/índice | no se buscó; ver `note` |
| `no_sequence` | la fila no tiene `aa` | no se buscó |
| `error` | el motor falló | revisar `warnings`/`note` |

Reglas:

- `hit` no cambia `path`, `label`, `evidence_strength` ni `next_step`;
- `hit` no cierra función, ni reacción, ni complejo;
- `no_hit` no se convierte en "sin función";
- las etiquetas A/B/C siguen siendo las de `functional.py`.

## Límites

- **Cobertura**: solo las cinco fuentes detectadas. No hay búsqueda contra PDB,
  AlphaFold ni bases no locales.
- **Sin red**: no hay APIs, descargas ni actualización automática de bases.
- **Una base por tipo**: si hay UniRef y TrEMBL, solo se consulta la primera
  según prioridad (no hay parámetro para forzar la fuente). Para auditar una
  base concreta, ejecuta el módulo con un `data_root` que solo contenga esa
  base.
- **Pfam de panel**: si no está `Pfam-A.hmm`, se usa el panel local disponible
  (por ejemplo `panel_focus.hmm`). Un panel cura familias, no todo Pfam.
- **InterPro**: si solo existe como HMM/`.dmnd` local; no se consulta
  InterProScan ni la web.
- **eggNOG**: `eggnog.db` sin `.dmnd` no permite buscar secuencias nuevas.
- **Sensibilidad**: DIAMOND con `1e-5` y un target no es una búsqueda
  exhaustiva; homólogos muy divergentes pueden quedar fuera.
- **Coste**: DIAMOND contra UniRef/TrEMBL y hmmscan contra Pfam-A pueden tardar
  horas. No ejecutar sobre 647k C sin reducir a representantes; ver
  [`NO_HIT_CPU_PLAN.md`](NO_HIT_CPU_PLAN.md).
- **No sustituye**: dominio completo, TMHMM/SignalP, contexto génico, sintenía,
  estructura, taxonomía ni experimento.
- **No integrado**: este módulo es independiente; no se llama desde
  `webapp/functional.py` ni desde el pipeline/app. El llamador decide cuándo
  ejecutarlo.

## Pruebas

```bash
.venv/bin/python -m pytest tests/test_remote_homology.py -q
```

Cubren: DB ausente, fila A intacta, contrato de columnas, detección local,
búsqueda DIAMOND simulada, hmmscan simulado, eggNOG sin índice y ausencia de
red (socket bloqueado durante una búsqueda local).

## Lectura recomendada

- [`THREE_PROTEIN_PATHS.md`](THREE_PROTEIN_PATHS.md): qué significa A/B/C.
- [`CPU_GPU_LIMITS.md`](CPU_GPU_LIMITS.md): qué puede y no puede CPU.
- [`NO_HIT_CPU_PLAN.md`](NO_HIT_CPU_PLAN.md): priorización antes de buscar.
