# Cómo contribuir

Gracias por querer mejorar el anotador. Este proyecto se usa para investigación
y docencia, así que la prioridad es **no romper la trazabilidad**: cada número
debe poder rastrearse hasta su fuente.

## Antes de empezar

1. Lee `docs/ARCHITECTURE.md` y `scripts/README.md`.
2. Comprueba el entorno: `python scripts/00_doctor.py`.
3. Ejecuta las pruebas: `python -m pytest tests -q` y
   `cd scripts/ghost && PYTHONPATH=. python -m unittest discover -s tests`.
4. Si vas a tocar datos científicos, revisa `DATA_NOTICE.md` y no redistribuyas
   RegulonDB, KEGG, DoriC, STRING, DBTBS ni BATTER.

## Reglas de código

- **No cambies la profundidad** de `scripts/<etapa>/<archivo>.py`: 97 scripts
  calculan la raíz con `parents[2]`. Si renombras una carpeta, actualiza
  `run_pipeline.sh` y los `sys.path` en el mismo cambio.
- La lógica no importa Streamlit. Solo `app.py`, `*_ui.py`, `browser.py` y
  `requirements_ui.py` pueden hacerlo; así el resto se prueba sin interfaz.
- Las funciones nuevas llevan docstring en español y un ejemplo de uso cuando
  no es obvio.
- Los textos visibles van en español y se traducen: añade la cadena, ejecuta
  `python tools/extract_ui_strings.py` y completa `locales/en.json` y
  `locales/zh.json`. `tests/test_i18n.py` lo verifica.
- No uses rutas absolutas (`/home/...`). Usa `Path(__file__).resolve().parents[2]`
  o argumentos CLI.
- No apliques correcciones automáticas: marca `suggested` y deja la decisión al
  usuario.
- No declares `confirmed` sin evidencia EC, contexto o experimento.

## Cómo proponer un cambio

1. Abre un *issue* con: genoma, etapa, comando exacto, resultado esperado y
   resultado obtenido. Si es un error visual, adjunta una captura.
2. Haz un *fork* y una rama (`feat/…`, `fix/…`, `docs/…`).
3. Añade una prueba que falle antes y pase después.
4. Envía un *pull request* describiendo el cambio y su impacto científico.
5. Si el cambio modifica un resultado, indica qué tabla o cifra cambia y por
   qué.

## Estilo

- Python 3.10+, sin dependencias nuevas si se pueden evitar.
- Mensajes de interfaz claros, sin jerga innecesaria; piensa en un estudiante
  de secundaria.
- Nada de emojis en nombres de archivo ni en columnas de datos.
- Los commits explican el *porqué*, no solo el *qué*.

## Conducta

Al participar aceptas el `CODE_OF_CONDUCT.md`. No se toleran ataques
personales, acoso ni descalificaciones. La crítica científica es bienvenida;
la falta de respeto, no.

## Licencia de las contribuciones

Al enviar un *pull request* aceptas que tu contribución se publique bajo la
licencia MIT del proyecto (ver `LICENSE`), sin derechos adicionales.
