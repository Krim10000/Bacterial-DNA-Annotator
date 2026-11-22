# Política de uso de Modal y GPU

## Regla crítica

Modal permite **máximo 10 GPU activas simultáneas por cuenta**.
Un agente que use Modal **debe anunciarlo** y **cerrar todo lo que abrió**.

## Antes de usar

1. Anunciar: `USANDO MODAL GPU: <app> · <n> GPUs · <motivo>`.
2. Verificar que no haya apps/contenedores previos:
   ```bash
   modal app list
   modal container list
   ```
3. Ejecutar **sin `--detach`** y con `max_containers=1`.
4. Usar lotes pequeños (5 secuencias), no 10–20 a la vez.

## Después de usar

1. Confirmar que la app efímera terminó.
2. Si quedó algo:
   ```bash
   modal app list
   modal app stop <app-id>
   modal container list
   modal container stop <container-id>
   ```
3. Reportar: `MODAL CERRADO: <app> · GPUs liberadas: <n>`.
4. No dejar `modal deploy`, `--detach`, `keep_warm` ni contenedores idle.

## Reglas para agentes

- **Un agente = un lote GPU a la vez.**
- No compartir app persistente entre agentes sin coordinación.
- No usar GPU para Foldseek: Foldseek corre en CPU local.
- Si se alcanza el límite de 10 GPUs, **no reintentar en bucle**; esperar o reducir concurrencia.
- Documentar app, GPU, tiempo y coste.

## Estado actual

- `modal_esmfold.py` tiene `max_containers=1`, `timeout=1800`.
- No se ha ejecutado: falta `modal token new`.
- GPUs activas de este agente: **0**.
