# Estado de la GPU y procesos que quedaron pendientes

🌐 **Español** · [English](GPU_STATUS.en.md)

**Resumen:** la anotación, el navegador, las brechas metabólicas, los operones y
la interfaz funcionan **sin GPU**. Lo que falta por GPU es opcional o quedó
cerrado en la nube: el entrenamiento completo del modelo de ADN, su evaluación y
ESMFold local.

Máquina de referencia de esta versión: `torch 2.14.1+cpu`,
`torch.cuda.is_available() == False`, sin `nvidia-smi` ni `/dev/nvidia*`.
La GUI detecta esto y lo muestra en **⏱️ ¿Qué necesito y cuánto tarda?**.

## 1 · Qué necesita GPU y qué no

| Etapa | Estado en v1.0 | ¿GPU? | Cómo se ejecuta | Impacto |
|---|---|---|---|---|
| Anotación (CDS, tRNA, rRNA, genes fantasma) | ✅ funciona | No | GUI o `dna-annotate` | núcleo |
| Brechas, cierres, transporte, operones | ✅ funciona | No | GUI / `run_pipeline.sh` | núcleo |
| Navegador de genoma | ✅ funciona | No | GUI | núcleo |
| Filtros estructurales (Foldseek, P2Rank, AlphaFill) | ✅ funciona | No | GUI / `scripts/structure/` | incluido |
| ESMFold por Neurosnap (API) | ☁️ opcional | No (remoto) | `--esmfold "neuros nap"` | ~0,04 créditos/proteína |
| ESMFold local | ⛔ bloqueado | ≥16 GB VRAM | `scripts/structure/04_local_esmfold.py` | no incluido |
| ESMFold en Modal | ⛔ nunca ejecutado | A10G | `scripts/structure/modal_esmfold.py` (deprecado) | no incluido |
| Entrenar modelo de ADN v2 | ⚠️ parcial (epoch 6 y 9 de 25) | A10G | `modal_train_v2.py` | no lo usa el anotador |
| Evaluación completa del modelo | ⛔ bloqueado | GPU | `scripts/13_comprehensive_eval.py --n 100` | no incluido |
| Boltz-2 a escala | ☁️ parcial (24 predicciones) | A10G | `scripts/allostery/04_boltz_modal.py` | evidencia recogida |
| MM/GBSA, ABFE, dinámica molecular | ⛔ fuera de v1.0 | GPU | `docs/ALOSTERIA_V1_CIERRE.md` | cerrado |

## 2 · Evidencia de lo que quedó bloqueado

- **Modelo de ADN.** `modal_train_v2.py:16` pide `gpu="A10G"`. Solo existen los
  checkpoints de las épocas 6 y 9 (`dna_model_v2_epoch6.pt`,
  `dna_model_v2_epoch9.pt`, 48,9 MB, ~13,8 M parámetros). El plan de 25 épocas
  no terminó: `checkpoints_descargados/checkpoints_v2/` está vacío.
  `logs/09_train_model.log` mide 0 bytes: nunca se intentó en local.
- **Evaluación.** `reports/README.md:29-35` documenta que la corrida local se
  canceló tras 45 minutos; `reports/eval_comprehensive.json` usa solo 8 ventanas
  y `device: "cpu"`, así que **no es una evaluación**. `ROADMAP.md` afirmaba
  97,3 M de parámetros y 25 épocas; los pesos reales son 11–14 M y 9–10 épocas.
  Esta versión corrige la afirmación y no presenta el modelo como evaluado.
- **ESMFold local.** `scripts/structure/04_local_esmfold.py:38-40` aborta sin
  CUDA; `webapp/esmfold_tools.py:78-88` informa «sin GPU NVIDIA visible ·
  PyTorch sin CUDA · falta openfold». La alternativa es la API de Neurosnap.
- **Boltz-2.** Se hicieron 24 predicciones y el resultado fue negativo
  (`docs/NEGATIVE_RESULTS.md:11-24`); no se lanzó la campaña de ~900 co-folds
  (20 h de A10G, 25–40 USD).
- **Logs.** No hay ni una ocurrencia de OOM/CUDA en los 913 logs: el trabajo
  pesado se externalizó a Modal o Neurosnap, no falló en local.

## 3 · Qué significa para quien usa la app

- Sin GPU obtienes **anotación y priorización completas**. Lo que no obtienes
  son predicciones de energía libre ni la evaluación del modelo de ADN.
- La app **no intenta** cargar CUDA ni descargar pesos: si no hay GPU, las
  opciones que la necesitan se muestran como no disponibles, con el motivo.
- ESMFold por Neurosnap envía la secuencia a un servicio externo y consume
  crédito. La app lo advierte antes de usarlo; no se activa solo.

## 4 · Cómo continuar en una máquina con GPU

```bash
# 1. Verificar que la GPU existe y que torch la ve
nvidia-smi
python -c "import torch; print(torch.cuda.is_available())"

# 2. Entrenamiento (Modal; revisar docs/MODAL_GPU_POLICY.md antes)
modal run modal_train_v2.py

# 3. Evaluación con la GPU local
python scripts/13_comprehensive_eval.py --n 100 --window 2048

# 4. ESMFold local (requiere pesos y OpenFold)
python scripts/structure/04_local_esmfold.py --help
```

Coste orientativo documentado: Boltz-2 a escala 25–40 USD; ABFE 100–150 GPU-h
por sistema; MD 10 ns × 3 réplicas 51–135 USD. Nada de esto es necesario para
usar el anotador.

## 5 · Decisión para v1.0

El modelo de ADN se publica como **experimental y no evaluado**. Los
checkpoints no se incluyen en el repositorio: se documenta cómo regenerarlos.
La anotación y la GUI son el producto estable de esta versión.
