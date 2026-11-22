# ESMFold para cierre de brechas — local y Neurosnap

## Regla de transparencia

Cuando el usuario marca una alternativa de estructura, la GUI debe decir
**explícitamente que usa ESMFold**:

- `ESMFold local`
- `ESMFold vía API Neurosnap`

ESMFold predice estructura 3D desde secuencia. **No confirma función,
sustrato, dirección metabólica ni cierre de brecha.** Una estructura de
alta confianza es apoyo computacional, no experimento.

## Evaluación de la máquina local

La GUI evalúa la máquina donde corre Streamlit/pipeline, no el navegador:

- GPU NVIDIA y VRAM (`nvidia-smi` o PyTorch CUDA).
- RAM total.
- disco libre.
- PyTorch/CUDA, `esm`, `openfold`.

Niveles:

| nivel | criterio | uso |
|---|---|---|
| `recommended` | GPU ≥16 GB VRAM, RAM ≥32 GB, disco ≥20 GB, CUDA | local |
| `possible_slow` | GPU ≥10 GB, RAM ≥16 GB, disco ≥15 GB | local puntual, no lotes |
| `not_recommended` | sin GPU o recursos insuficientes | API Neurosnap |

Si no cumple, la GUI debe mostrar el motivo y ofrecer Neurosnap.

- Código: `webapp/esmfold_tools.py` (`detect_capability`).
- UI: `webapp/esmfold_ui.py`.
- Script local: `scripts/structure/04_local_esmfold.py`.

## Ruta local

Requisitos:

```bash
pip install 'fair-esm[esmfold]'
pip install 'openfold @ git+https://github.com/aqlaboratory/openfold.git@4b41059694619831a7db195b7e0988fc4ff3a307'
python scripts/structure/04_local_esmfold.py \
  --fasta data/metabolism/targets/structural_candidates_5.faa \
  --out data/metabolism/targets/structures/local
```

- No envía secuencias a la nube.
- Requiere pesos ESMFold (~GB) y GPU.
- OpenFold puede fallar por CUDA/compilación; no ocultar ese error.

## Ruta Neurosnap API

### Credencial

`ESM_API_KEY` está en `.env` (ignorado por git). No imprimir ni commitear.

```bash
cd <raíz-del-repositorio>
set -a; . ./.env; set +a
```

### Verificar

```bash
curl -s https://neurosnap.ai/api/jobs \
  -H "X-API-KEY: $ESM_API_KEY"
```

### Enviar ESMFold

Endpoint:

```text
POST https://neurosnap.ai/api/job/submit/ESMFold?note=...
```

Multipart:

```python
import json, os, requests

payload = json.dumps({"aa": {"protein": "MAQTTR..."}, "dna": {}, "rna": {}})
files = {
    "Input Sequences": (None, payload),
    "Number Recycles": (None, "6"),
}
r = requests.post(
    "https://neurosnap.ai/api/job/submit/ESMFold",
    params={"note": "dna_poc"},
    headers={"X-API-KEY": os.environ["ESM_API_KEY"]},
    files=files,
)
job_id = r.json()   # string JSON, no objeto
```

Reglas del nombre: máximo 50 caracteres, solo alfanumérico `_` `-`.
Sanear IDs con `.` o `|`.

### Poll y descarga

```text
GET /api/job/status/{job_id}   -> "pending" | "running" | "completed"
GET /api/job/data/{job_id}     -> {"out": [[filename, size], ...]}
GET /api/job/file/{job_id}/out/{filename}
```

El endpoint de archivo responde **307**; seguir redirect.

### Script del repo

```bash
set -a; . ./.env; set +a
.venv/bin/python scripts/structure/05_neurosnap_fold.py \
  --fasta data/metabolism/targets/structural_candidates_5.faa \
  --out data/metabolism/targets/structures/neuros_nap \
  --limit 1 --note dna_poc
```

Salidas: `*.pdb`, `*_metrics.json`, `*_job.json`.

## Coste

- **$1.00/crédito**; académico **$0.90**.
- ESMFold: `resource_cost = 0.0033325` créditos/segundo.
- ESMFold de 300 aa: **~0,04 créditos** (~$0.04).
- Free: **1 crédito**.
- **$10 = 10 créditos**.
- PayPal no documentado; la compra es por créditos.

## GUI

En `🧪 Brechas → 5 · Cierre candidato` y en el resultado subido
`Brechas → 5 · Cierre por estructura (ESMFold)`:

1. Radio: `No usar ESMFold` / `ESMFold local (usa ESMFold)` /
   `ESMFold vía API Neurosnap (usa ESMFold)`.
2. Al marcar una ruta, aviso explícito de que usa ESMFold.
3. Tabla de capacidad local y motivos.
4. Si local no es viable, recomendar Neurosnap.
5. Si hay ORF y `ESM_API_KEY`, botón de envío con aviso de coste.
6. Mostrar job, pLDDT, pTM y PDB.

## Pipeline

```bash
python -m webapp.pipeline genoma.fna \
  --outdir results/ \
  --esmfold "neuros nap" \
  --esmfold-gaps 3
```

- `--esmfold off` (default), `local`, `"neuros nap"`.
- `--esmfold-gaps N` limita el gasto (default 3).
- `--esmfold-out DIR` guarda PDB/métricas.
- Salida: `GENOMA.esmfold_closure.tsv` y `GENOMA.summary.json`.
- Si local no es viable, el pipeline avisa y no finge haber plegado.
- Si pLDDT ≥80 se etiqueta `structure_predicted`; **no** `structural_supported`.
  `structural_supported` exige mapear residuos catalíticos/sitio activo sobre la
  estructura y comprobar conservación. **Nunca** `confirmed` sin experimento.

## Evidencia ya obtenida

- Job Neurosnap: `6b021688c7790c6f67ecdfa7`
- pLDDT medio: **97,47**
- pTM: **0,97**
- Foldseek local: top hit **O06745**, bits **945**
- Interpretación: estructura plausible del ORF; no cierra la brecha.

## Seguridad y límites

- No imprimir `ESM_API_KEY`.
- Confirmar autorización antes de enviar secuencias a un servicio externo.
- ESMFold es monomer-only; para complejos usar AlphaFold/Boltz.
- pLDDT alto no prueba actividad, sustrato ni dirección.
- AFDB solo cubre homólogos UniProt; el ORF exacto requiere plegado.

## Archivos

- `webapp/esmfold_tools.py`
- `webapp/esmfold_ui.py`
- `scripts/structure/04_local_esmfold.py`
- `scripts/structure/05_neurosnap_fold.py`
- `data/metabolism/targets/STRUCTURE_PIPELINE.md`

## Evidencia multi-fuente (pipeline)

Después de tener PDB:

```bash
.venv/bin/python scripts/structure/07_foldseek_pdb.py
.venv/bin/python scripts/structure/08_alphafill.py
.venv/bin/python scripts/structure/09_consensus.py
```

- **07**: Foldseek API contra PDB100; fallback AFDB local.
- **08**: AlphaFill por UniProt esperado (ligandos/cofactores).
- **09**: consenso con sitios activos locales.

Salidas:
- `FOLDSEEK_PDB_HITS.tsv`
- `ALPHAFILL_HITS.tsv`
- `STRUCTURAL_EVIDENCE.tsv` / `.md`

Pipeline:

```bash
python -m webapp.pipeline genoma.fna --outdir results/ \
  --esmfold "neuros nap" --esmfold-gaps 3 --structure-evidence
```

Salida extra: `GENOMA.structural_evidence.tsv`.

GUI: en el cierre ESMFold aparece **“Cruzar Foldseek + AlphaFill”**. Es un botón explícito; no se ejecuta solo.

Reglas:
- `multi_source_support` = ≥2 fuentes coinciden; **no** es cierre.
- `closure=no` siempre.
- AlphaFill puede tener clashes; se filtra ≤20.
- Foldseek top puede no ser el homólogo esperado; se conserva el resultado real.
