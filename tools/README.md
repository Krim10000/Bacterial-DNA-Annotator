# Herramientas de desarrollo

No forman parte del pipeline científico; sirven para mantener la interfaz, las
traducciones y el paquete publicable.

| herramienta | para qué |
|---|---|
| `extract_ui_strings.py` | extrae las cadenas visibles de la GUI a `webapp/locales/_source.json`; se ejecuta después de añadir textos |
| `build_demo_data.py` | regenera `webapp/demo/` desde `webapp/samples/demo.fna` (genes y operones ilustrativos) |
| `gui_smoke.py` | prueba con Playwright que la app carga en ES/EN/ZH, que la demo dibuja genes y que el clic abre la ficha |
| `make_release_repo.sh` | construye un repositorio limpio (sin datos, modelos ni historial) listo para publicar |
| `audit_script_paths.py` | comprueba que todas las rutas `scripts/...` citadas existan |
| `reorg_dry_run.py` | simula la reorganización de carpetas: colisiones y referencias a reescribir |
| `reorg_apply.py` | aplica el manifiesto con `--confirm` (mueve, renombra y reescribe); sin el flag solo simula |

## Flujo al añadir un texto nuevo

```bash
# 1. escribe el texto en español en webapp/
python tools/extract_ui_strings.py
# 2. completa las claves nuevas en webapp/locales/en.json y zh.json
python -m pytest tests/test_i18n.py -q
```

## Flujo antes de publicar

```bash
python -m pytest tests -q
cd scripts/ghost && PYTHONPATH=. python -m unittest discover -s tests
python tools/build_demo_data.py
python tools/make_release_repo.sh dist/DNA_TOKEN_BAC_v1.0 --force
```

El script de release verifica que no haya archivos >50 MB y deja un commit en la
rama `v1.0`.
