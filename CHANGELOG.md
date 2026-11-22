# Historial de cambios

Formato basado en [Keep a Changelog](https://keepachangelog.com/es/1.1.0/).
Versionado [SemVer](https://semver.org/lang/es/).

## [1.0.0] — 2026-11-22

Primera versión pensada para publicarse y para que la use alguien que no
conoce el proyecto.

### Añadido

- **Arquitectura de módulos autodescubiertos** (`modules/`): cada módulo declara
  sus opciones, estado, tiempo, autoprueba y ejecución. La vista **🧩 Módulos**
  se construye sola; añadir un módulo no requiere editar `app.py`. Incluye CLI
  (`python -m modules --list`, `--describe`, `--selftest`) y pruebas por módulo.
  Regulación, móviles, membrana, taxonomía, targets y validación ejecutan un
  **resumen verificable** de datos ya calculados (no fingen recalcular).
- **Interfaz en español, inglés y chino simplificado**, con selector en la
  barra lateral y enlaces `?lang=es|en|zh`. Traducciones completas (996 cadenas), con pruebas que fallan si falta una clave o se rompe un
  placeholder.
- **Demo incluida**: 150 kb de *Chlorobaculum tepidum* anotados desde cero
  (`webapp/demo/`), para que el navegador y las fichas funcionen en un clon
  recién descargado sin datos restringidos.
- **Navegador de genoma clicable**: al pulsar un gen, operón, promotor,
  horquilla o sospecha metabólica se abre una ficha con producto, coordenadas,
  hebra, EC, confianza, cajas −10/−35, ΔG o genes vecinos. Incluye respaldo por
  URL si el navegador no procesa el valor del componente.
- **Instalador para Ubuntu** (`install.sh`) con lanzador en el menú de
  aplicaciones, icono propio, opciones `--with-data`, `--with-tools`, `--full`,
  `--dry-run` y `--uninstall`.
- **Panel de requisitos, tiempos y glosario** para estudiantes: qué es DIAMOND,
  HMMER, Rhea o una brecha, y cuánto tarda cada paso.
- **Créditos de autoría**: Camilo Rojas Bravo (USACH) y declaración explícita
  de la asistencia de DeepSeek Harness con DeepSeek (noviembre de 2026).
- **Licencias y datos**: `LICENSE` (MIT), `THIRD_PARTY_LICENSES.md`,
  `DATA_NOTICE.md`, `NOTICE` y textos GPL-3.0/Apache-2.0.
- **Documentación de GPU**: `docs/GPU_STATUS.md` con las etapas bloqueadas,
  la evidencia y cómo continuar en una máquina con GPU.
- **Arquitectura y orden de scripts**: `docs/ARCHITECTURE.md` y
  `scripts/README.md` con índice por etapa, tiempos, GPU y estado.
- **Pruebas y CI**: pruebas de i18n, SVG, imports, instalador y anotación real;
  CI en Python 3.10/3.11/3.12 y arranque de la GUI.
- **Empaquetado 1.0.0**: `pyproject.toml` con datos de paquete (locales y
  componente), `MANIFEST.in` y `CITATION.cff`.

### Reorganización de scripts

- `scripts/` pasó de 27 a 17 carpetas por dominio: `data_prep`, `annotation`,
  `metabolism`, `regulation`, `structure`, `evidence`, `validation`, `ghost`,
  `allostery`, `taxonomy`, `membrane`, `mobile`, `targets`,
  `interpretability`, `oric`, `lib` y `tests`.
- `12_operon_finder`, `13_tss_tts` y `10_chains` se fusionaron en
  `regulation`; se actualizaron los 11 `sys.path` que los consumían.
- 23 movimientos/fusiones, 127 elementos movidos y 548 referencias reescritas.
- Se corrigieron rutas que el plan no detectaba: `webapp/annotator.py` al
  paquete `ghost`, rutas partidas `ROOT / "scripts" / "<carpeta>"` y dos
  `ANNOTATOR` que apuntaban a archivos inexistentes.
- El ML custom (`07_ml`, `08_eval` y scripts de la raíz) quedó fuera de la
  reorganización de v1.0.

### Corregido

- Las flechas de los desplegables mostraban la ligadura `keyboard_arrow_right`
  o `B8` cuando la fuente de iconos no cargaba; ahora son CSS puro.
- `_save_ui` ya no borra el idioma al guardar la vista o el genoma.
- El instalador y el lanzador detectan el puerto libre y no cierran la ventana
  ante un error.
- Se retiró `cobra` de los requisitos: no se usaba y añadía una dependencia
  LGPL innecesaria.
- Se corrigió la afirmación del `ROADMAP` sobre el modelo (97,3 M de parámetros
  y 25 épocas): los pesos reales son 11–14 M y 9–10 épocas, y no está evaluado.

### Cambiado

- README reescrito para tres públicos: estudiante, investigador y quien va a
  publicar el repositorio.
- Los scripts de legado (`scripts_orginales/`, `_review/`) quedan fuera de la
  versión publicada; se conservan en el archivo local.
- Los datos de terceros no se redistribuyen; se descargan con scripts.

### Limitaciones conocidas

- El modelo de ADN es experimental y no está evaluado (sin GPU).
- ESMFold local no está incluido; se ofrece la API de Neurosnap.
- Las predicciones de operones tienen F1 ≈ 0,78–0,81.
- La cobertura de TTS es baja (16 % de las fronteras).
- La app necesita datos locales para los modos que exploran los 400 genomas.

