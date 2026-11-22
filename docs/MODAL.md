# Uso de Modal en este proyecto

> **Límite duro: Modal admite 10 GPU activas simultáneas por workspace.** Un agente que
> deje contenedores abiertos bloquea a todos los demás. Por eso: **avisar antes, cerrar
> después**.

## Obligaciones de cualquier agente (humano o IA) que use Modal

1. **Anunciar** que va a usar Modal y para qué (en el mensaje o en el commit).
2. **No escalar**: usar 1 contenedor salvo que el trabajo lo exija. `.map()` sobre 35
   jobs abre hasta el límite de concurrencia si no se acota.
3. **Cerrar** al terminar: `bash scripts/allostery/modal_guard.sh stop-all`.
4. **Verificar**: `bash scripts/allostery/modal_guard.sh status` debe salir con
   código 0 y `0 activas`.

## Comandos

```bash
# ¿qué está activo?
bash scripts/allostery/modal_guard.sh status

# cerrar todo lo que no esté parado
bash scripts/allostery/modal_guard.sh stop-all

# inspección directa
modal app list
modal container list
modal app stop -y <app_id>          # una app
modal container stop <container_id> # un contenedor suelto
```

`status` devuelve código 1 si hay algo activo — sirve como guardia en CI o en un
`pre-commit`.

## Convenciones del repo

| Regla | Motivo |
|---|---|
| `modal run --detach` solo con **un** orquestador (`orchestrate`) | en detached, solo sobrevive la última función disparada |
| Preferir `.map()` acotado o `--serial` en trabajos cortos | evita abrir 30 contenedores sin querer |
| GPU por defecto: `A10G` | suficiente para ANM, Boltz-2 y docking; más GPU no acelera trabajo CPU |
| Los volúmenes (`boltz-weights`, `alloreg-results`, `alloreg-anm-results`) no consumen GPU | se pueden dejar montados sin coste |
| Cerrar también los contenedores sueltos (`modal container list`) | una app parada puede dejar un contenedor vivo unos minutos |

## Inventario de scripts de este módulo

| Script | GPU | Qué abre | Coste típico |
|---|---|---|---|
| `04_boltz_modal.py` | A10G | 1 orquestador + hasta 3 `run_boltz` | ~0,5 USD por 24 predicciones |
| `06_anm_modal.py` | A10G | 35 jobs `.map()` (ANM es CPU) | centavos |
| `07_laci_oligomer.py` | A10G | 1 job | centavos |
| `08_hinge_specificity.py` | A10G | 9 jobs `.map()` | centavos |

## Checklist (copiar en cada sesión)

- [ ] avisé que uso Modal
- [ ] limité la concurrencia (≤3 contenedores)
- [ ] `modal_guard.sh stop-all`
- [ ] `modal_guard.sh status` → `0 activas`, exit 0

## Historial de uso

| fecha | qué | apps | cerrado |
|---|---|---|---|
| 2026-11-21 | Boltz-2 (smoke, six, repaired) + ANM (35 TFs) + LacI oligómero + especificidad de bisagra | 11 | sí, 11/11 stopped, 0 contenedores |
