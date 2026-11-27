# P2Rank de baja prioridad

## Proceso

```text
scripts/sigma/114_p2rank_lowcpu.py
nice 19
1 hilo
arranca solo si CPU < 90 %
reanudable
```

## Cola

```text
1.972 estructuras AFDB/PDB sin P2Rank
data/sigma/structures/p2rank_lowcpu/queue.tsv
```

No incluye Kaggler ya procesadas ni los 169 sigma ya cubiertos.

## Salida

```text
data/sigma/structures/p2rank_lowcpu/<acc>/
data/sigma/structures/p2rank_lowcpu/status.json
data/sigma/structures/p2rank_lowcpu/daemon.log
```

## Reanudar

```bash
nice -n 19 .venv/bin/python scripts/sigma/114_p2rank_lowcpu.py --daemon
```

## Supervisión

Tarea programada cada 45 minutos:

```text
p2rank-supervisor
revisa status.json y daemon.log
si está atascado, relanza
```

No lanza P2Rank manualmente ni satura CPU.

## Regla

```text
pocket ≠ función
pocket ≠ promotor
```

Es un feature estructural, no una anotación.
