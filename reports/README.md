# Reportes de evaluación

## Estado

| script | estado | motivo |
|---|---|---|
| `08_eval/11_eval_clustering.py` | import arreglado, no ejecutado | requiere GPU para ser útil |
| `08_eval/12_eval_probe.py` | import arreglado, no ejecutado | ruta de checkpoint v2/epoch25 inexistente |
| `13_comprehensive_eval.py` | **ejecutado en modo humo** | CPU sin GPU: n=8 no es una métrica |

## Lo que sí hay

`eval_comprehensive.json` — checkpoint `dna_model_v2_epoch6.pt`, tokenizer v2
(4096), CPU, `n=8`, ventana promotores 256.

**Estos números no son una evaluación.** Con 8 ventanas por clase el AUROC tiene
un error estándar enorme; sirven solo para demostrar que el arnés carga el
checkpoint, tokeniza, extrae embeddings y entrena la sonda logística. No los
cites como rendimiento del modelo.

| test | AUROC (n=8) | lectura |
|---|---|---|
| Promotores (TATAAT) | 0,667 | sin valor estadístico |
| RBS (AGGAGG) | 0,333 | sin valor estadístico |
| Start vs Stop | 0,333 | sin valor estadístico |

## El bloqueo real

`torch 2.14.1+cpu`, sin GPU. Una corrida con `--n 40 --window 512` no terminó en
~45 min y se canceló. El modelo es de 13,8 M de parámetros, pero la atención
sobre ventanas de 256–512 en CPU multiplicada por cientos de ventanas es
inviable en esta máquina.

No es un bug del script: es hardware. El entrenamiento se hace en Modal; la
evaluación debería hacerse allí con el mismo checkpoint.

## Reproducir

```bash
# humo local (segundos)
python scripts/13_comprehensive_eval.py \
    --checkpoint dna_model_v2_epoch6.pt --n 8 --window 256

# evaluación real (necesita GPU o Modal)
python scripts/13_comprehensive_eval.py \
    --checkpoint dna_model_v2_epoch6.pt --n 100 --window 2048
```

Los scripts `11` y `12` quedaron con el `sys.path` corregido tras la
reorganización de carpetas, pero siguen apuntando a artefactos que no están en
disco (`checkpoints_descargados/checkpoints_v2/..._epoch25.pt`). Para
ejecutarlos hay que pasarles el checkpoint existente.
