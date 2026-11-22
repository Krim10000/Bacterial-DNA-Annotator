# Benchmark de modelos pretrained — capacidad del equipo y resultados

## 1. Capacidad real del equipo

| recurso | valor |
|---|---|
| CPU | Intel i7-8550U, 4 núcleos / 8 hilos, 1,8 GHz |
| RAM | 23 GB |
| GPU | **ninguna** (`torch.cuda.is_available() == False`) |
| Disco libre | 75 GB |
| Red | HuggingFace accesible (HTTP 200) |
| Stack | Python 3.10, torch 2.14 CPU, transformers 4.46.3, sentencepiece, sklearn |

Consecuencia: se pueden evaluar modelos de **27M–138M** en CPU con subsets;
**Evo2 1B/7B, NT-v2 250M/500M y fine-tuning completo** requieren GPU
(Modal A10G) para ser prácticos.

## 2. Modelos descargados y evaluados

| modelo | params | dim | contexto | tokenizer | licencia | estado |
|---|---:|---:|---:|---|---|---|
| **ProkBERT-mini-long** | 26,4M | 384 | 4 kb | 6-mer shift 2 | CC-BY-NC-4.0 | ✅ descargado/evaluado |
| **NT-v2 50M multi-species** | 55,9M | 512 | 12 kb | 6-mer solapado | CC-BY-NC-SA-4.0 | ✅ descargado/evaluado |
| **v2 epoch9 (propio)** | 12,2M | 384 | 32 kb efectivo | BPE 4.096 | propio | ✅ referencia |
| Evo2 1B-8K | 1B | 1.920 | 8 kb | carácter | Apache 2.0 | ⏸️ CPU no práctico |
| BacLM 350M | 350M | 960 | 2 kb | carácter | Apache 2.0 | ⏸️ CPU lento, contexto 2 kb |
| DNABERT-2 | 117M | 768 | 512 tokens | BPE | MIT-like | ⏸️ contexto corto |
| Mistral-DNA 138M bacteria | 138M | 768 | 512 tokens | BPE | — | ⏸️ contexto corto |
| Caduceus / HyenaDNA | 1–47M | — | 131 kb / 1M | carácter | Apache | ⏸️ preentrenados en humano |
| Bacformer | 27M | 480/960 | 6.000 proteínas | proteína | Apache | ❌ no es ADN |

Fuentes: [ProkBERT-mini-long](https://huggingface.co/neuralbioinfo/prokbert-mini-long),
[NT-v2 50M](https://huggingface.co/InstaDeepAI/nucleotide-transformer-v2-50m-multi-species),
[Evo2](https://pypi.org/project/evo2/), [BacLM](https://huggingface.co/macwiatrak/baclm-350m-masked),
[BacBench](https://github.com/macwiatrak/Bacbench).

## 3. Protocolo

- **600 ventanas de 256 pb**, 20 genomas, 20 filos, estratificadas:
  180 CDS, 120 tRNA, 60 rRNA, 150 intergénico, 60 upstream, 30 transposón.
- Embeddings **congelados** (mean pooling con attention mask).
- Clasificación **leave-one-genome-out** (LOGO) con `StandardScaler + SelectKBest(50) + LogisticRegression`.
- Ajuste por **GC/GC² y genoma** (efectos fijos) antes de las firmas.
- Baselines: **GC**, **3-mer**, **conteo de TATAAT**.
- Promotores: 400 ventanas de 256 pb (200 TATAAT + 200 dinucleótido-barajado)
  de 4 genomas no vistos, misma partición CV.

## 4. Resultados de componentes (600 ventanas)

| modelo | balanced acc | AUC CDS | AUC tRNA | AUC rRNA | AUC intergénico |
|---|---:|---:|---:|---:|---:|
| **ProkBERT-mini-long** | **0,825** | **0,981** | **0,992** | **0,967** | **0,909** |
| NT-v2 50M | 0,510 | 0,818 | 0,840 | 0,663 | 0,626 |
| v2 epoch9 | 0,376 | 0,669 | 0,708 | 0,701 | 0,528 |
| 3-mer baseline | 0,502 | 0,719 | 0,732 | 0,842 | 0,638 |
| GC baseline | 0,361 | 0,503 | 0,652 | 0,699 | 0,629 |

Firmas seleccionadas (dims con `q<0.05`, `|d|≥0.25`, consistencia ≥0.8):

| modelo | CDS | tRNA | rRNA | intergénico | upstream | transposón |
|---|---:|---:|---:|---:|---:|---:|
| **ProkBERT** | **188** | **223** | **149** | **99** | 9 | 51 |
| NT-v2 50M | 102 | 81 | 22 | 1 | 1 | 0 |
| v2 epoch9 | 49 | 193 | 7 | 0 | 0 | 0 |

Lectura:

- **ProkBERT es el mejor por un margen grande** en todos los componentes,
  incluido intergénico y upstream. No es solo composición: el baseline 3-mer
  se queda en 0,502 de balanced accuracy y ProkBERT llega a 0,825.
- **NT-v2 50M** mejora a v2 en CDS/tRNA/intergénico, pero está cerca del
  baseline 3-mer (0,510 vs 0,502). No es un salto cualitativo.
- **v2 epoch9** queda por debajo del baseline 3-mer; confirma que el modelo
  propio está subentrenado y que su embedding no separa componentes.

## 5. Promotores (TATAAT vs dinucleótido-barajado)

| modelo | AUROC | accuracy |
|---|---:|---:|
| **ProkBERT-mini-long** | **0,988** | 0,953 |
| NT-v2 50M | 0,908 | 0,833 |
| v2 epoch9 | 0,767 | 0,703 |
| Conteo de TATAAT (baseline) | 0,850 | — |
| GC (baseline) | 0,455 | — |

- ProkBERT **supera al conteo de TATAAT** (0,988 vs 0,850): no solo detecta el
  motivo, aprende contexto.
- NT-v2 también supera el baseline (0,908).
- v2 está **por debajo** del conteo de TATAAT: su AUROC de promotor no es
  mejor que contar el motivo.
- Caveat: los negativos son dinucleótido-barajados; parte de la señal es
  “secuencia natural vs barajada”. Para un test definitivo, usar negativos con
  motivo mutado y flancos intactos.

## 6. Coste de extracción en CPU

| modelo | s/ventana (256 pb) | 600 ventanas |
|---|---:|---:|
| v2 | 0,066 | 40 s |
| NT-v2 50M | 0,234 | 141 s |
| ProkBERT | 0,249 | 150 s |

Con 8 hilos y carga alta del host. Para 10.000 ventanas: 10–40 min por modelo;
para los 400 genomas completos, GPU.

## 7. Qué significa para el modelo V3

La propuesta V3 (conv-transformer 20–30M desde cero) **ya no es la primera
opción**. Los datos dicen:

1. **No entrenar desde cero para componentes locales.** ProkBERT, 27M y
   bacterial, supera al modelo propio por 0,45 de balanced accuracy y aprende
   intergénico/upstream. Adoptarlo o fine-tunearlo es más barato y mejor.
2. **ProkBERT tiene contexto 4 kb.** Sirve para CDS/tRNA/rRNA/promotor local,
   no para operones largos (>4 kb). Para eso, NT-v2 12 kb o un modelo
   long-context.
3. **NT-v2 50M no justifica el salto.** Es mejor que v2, pero no supera
   claramente al baseline 3-mer. NT-v2 250M/500M o Evo2 1B podrían cambiar el
   resultado; hay que probarlos en GPU.
4. **V3 revisado**:
   - **Opción A (recomendada):** fine-tunear ProkBERT (LoRA o full, 27M) para
     componentes/promotor; usar NT-v2 12 kb para contexto de operón; stack de
     dos modelos.
   - **Opción B:** destilar ProkBERT + NT-v2 en un encoder de 12–24 kb con
     tokenizer v3c 256/4. El conv-transformer propuesto sirve como estudiante,
     no como entrenamiento desde cero.
   - **Opción C (fallback):** si el fine-tuning de ProkBERT no generaliza a
     genomas no vistos, entrenar el conv-transformer con v3c y ProkBERT como
     teacher.
5. **El tokenizer v3c 256/4 sigue siendo útil** para el modelo propio o para
   destilación; no se usa si se adopta ProkBERT/NT-v2 tal cual, porque ellos
   traen su tokenizer.

## 8. Advertencias

- **Posible fuga de pretraining:** ProkBERT se entrenó con 206.650 millones de
  nucleótidos de microbioma; los 20 genomas pueden estar en su corpus. El
  0,825 es zero-shot optimista. Hay que validar con hold-out por género/especie
  que no esté en el pretraining.
- **Subset de 600 ventanas:** suficiente para una señal comparativa, no para
  conclusiones finas de estabilidad.
- **Promotores:** negativos barajados; el baseline TATAAT ya da 0,85.
- **Sin GPU:** Evo2 1B, NT-v2 250M y fine-tuning completo no se probaron.
- **BacBench 2026** advierte que la posición génica y el parentesco confunden
  benchmarks bacterianos; añadir baseline de posición/distancia antes de
  reclamar superioridad biológica ([Zenodo 23484461](https://zenodo.org/records/23484461),
  [Zenodo 23605241](https://zenodo.org/records/23605241)).

## 9. Reproducir

```bash
cd dna_poc
source .venv/bin/activate
export HF_HOME=$PWD/.cache/huggingface

# Subconjunto
python scripts/interpretability/05_pretrained_benchmark.py \
    --model prokbert \
    --windows data/embeddings/_pretrained_subset600/windows.tsv \
    --out data/embeddings/pretrained/prokbert_subset600

python scripts/interpretability/05_pretrained_benchmark.py \
    --model ntv2 \
    --windows data/embeddings/_pretrained_subset600/windows.tsv \
    --out data/embeddings/pretrained/ntv2_50m_subset600

# Firmas + LOGO
python scripts/interpretability/03_component_signatures.py \
    --run data/embeddings/pretrained/prokbert_subset600 --permutations 0

# Promotores
python scripts/interpretability/05_pretrained_benchmark.py \
    --model prokbert \
    --windows data/embeddings/_promoter_windows/windows.tsv \
    --out data/embeddings/pretrained/prokbert_promoter
python scripts/interpretability/06_promoter_probe.py \
    --run data/embeddings/pretrained/prokbert_promoter
```

## 10. Siguiente paso

1. **Fine-tunear ProkBERT** (27M) con nuestros 400 genomas y etiquetas de
   componentes/promotor, hold-out por género; evaluar con los mismos gates.
2. **Probar NT-v2 250M y Evo2 1B en Modal A10G** para contexto 12 kb/8 kb.
3. **Test de promotor con negativos de motivo mutado** (no barajado) para
   separar “natural vs barajado” de gramática de promotor.
4. **Decidir**: stack ProkBERT + NT-v2, destilación, o V3 desde cero solo si
   los dos primeros fallan.
