## Qué cambia

<!-- Una frase clara. Si toca ciencia, indica la tabla o cifra afectada. -->

## Por qué

<!-- El problema real que resuelve. -->

## Cómo se probó

- [ ] `python -m pytest tests -q`
- [ ] `python tools/audit_script_paths.py`
- [ ] `python -m compileall -q scripts webapp modules tools`
- [ ] `cd scripts/ghost && PYTHONPATH=. python -m unittest discover -s tests`
- [ ] Si toca datos: comparé salida antes/después

## Checklist

- [ ] No incluí datos restringidos, tokens ni rutas absolutas.
- [ ] Actualicé la documentación y las traducciones si añadí textos.
- [ ] Si es un módulo nuevo: `python -m modules --selftest <id>` pasa.
