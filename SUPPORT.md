# Soporte y colaboración

Este es un proyecto de investigación y docencia. La ayuda más valiosa no es
dinero: son **datos, revisión, cómputo y tiempo de revisión**.

## Cómo pedir ayuda

1. Lee [`docs/INSTALL.md`](docs/INSTALL.md) y [`docs/GPU_STATUS.md`](docs/GPU_STATUS.md).
2. Prueba la demo: `python -m webapp.pipeline webapp/samples/demo.fna --outdir results/`.
3. Abre un *issue* con:
   - genoma o secuencia pública (o la demo),
   - comando exacto,
   - resultado esperado y obtenido,
   - log sin datos sensibles.

No incluyas datos restringidos, tokens ni información personal.

## Cómo colaborar

| área | qué necesitamos | dificultad |
|---|---|---|
| Benchmark | correr Bakta, Prokka y DFAST en 5–6 genomas y comparar CDS/tRNA/rRNA/tiempo/RAM | media |
| Chino | revisión nativa de términos científicos y del README | baja |
| macOS / WSL2 | probar `install.sh` y escribir guía | media |
| `protein_query` | proteína → homólogos/función con DIAMOND | media |
| Operones | benchmark replicable contra RegulonDB | media |
| Desconocidas | filtros sin GPU y categorías amplias | alta |
| GPU | créditos de Kaggle/Modal o una corrida controlada | variable |
| Datos | genomas públicos, validación experimental o PWM curados | alta |

Abre un *issue* con la etiqueta `collaboration` o escribe a
`Camilo.Rojas@usach.cl` indicando qué área te interesa.

## Crédito y autoría

- Los contribuyentes aparecen en `AUTHORS.md` y en la release.
- Un aporte sustancial puede optar a coautoría en el paper de JOSS; se decide
  por contribución real, no por jerarquía.
- Los datos y herramientas de terceros se citan en
  `THIRD_PARTY_LICENSES.md` y en la página de créditos.
- No se ofrece autoría honoraria ni se oculta la asistencia de IA.

## Donaciones

**Todavía no hay un enlace de donaciones.** Si el proyecto consigue adopción y
se abre un canal, será uno solo y transparente:

- **GitHub Sponsors** (preferido) o **Open Collective** con host fiscal;
- el dinero se destinaría a créditos de GPU, archivo/DOI y hosting;
- no habría funciones pagas ni resultados reservados;
- se publicaría un informe de uso.

Antes de abrir donaciones hay que resolver el régimen tributario en Chile y la
política de conflictos de interés. Mientras tanto, **créditos de cómputo o
tiempo de revisión valen más que una donación pequeña**.

## Citar

Si usas el proyecto, cita el DOI:

```text
Rojas Bravo, C. (2026). Bacterial DNA Annotator (v1.1.3).
https://doi.org/10.5281/zenodo.24469732
```

Una estrella en GitHub también ayuda a que otras personas lo encuentren.
