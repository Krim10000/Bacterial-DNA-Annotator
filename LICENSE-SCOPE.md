# Alcance de la licencia / License scope

Esta licencia cubre **solo el código original** de este repositorio, es decir,
`webapp/**`, `scripts/**`, `tools/**`, `tests/**` y la documentación escrita
por el autor, salvo que un archivo indique otra cosa.

No se aplica a:

- las dependencias de terceros ni a sus binarios;
- las bases de datos y modelos descargados;
- los textos de licencia de terceros incluidos en `LICENSES/`.

Cada uno conserva su propia licencia. La lista completa, con las obligaciones
que implica, está en `THIRD_PARTY_LICENSES.md`; el régimen de los datos está en
`DATA_NOTICE.md`.

This license covers only the original source in `webapp/`, `scripts/`,
`tools/` and `tests/`. Third-party dependencies, data and models keep their own
licenses (see `THIRD_PARTY_LICENSES.md` and `DATA_NOTICE.md`) and are not
relicensed here.

## Aviso importante / Important notice

Varias dependencias de este proyecto son **GPL-3.0**: pyrodigal, orcapy,
DIAMOND, tRNAscan-SE y barrnap. El código original puede usarse bajo MIT, pero
**cualquier distribución que empaquete esas dependencias** (una imagen Docker,
un servicio desplegado o un instalador con todo incluido) forma una obra
combinada y debe cumplir la GPL-3.0. El instalador de este repositorio las
instala por separado, no las redistribuye.
