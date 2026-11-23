# DOI con Zenodo

El repositorio ya incluye `.zenodo.json` con título, autoría, licencia,
keywords y relación con GitHub. Falta conectar Zenodo y publicar una release
nueva para que Zenodo la archive y emita el DOI.

## Opción A — integración GitHub (recomendada)

1. Entra a https://zenodo.org/ e inicia sesión **con GitHub**.
2. Menú superior → **Settings** → **GitHub**.
3. Busca `Krim10000/Bacterial-DNA-Annotator` y activa el interruptor.
4. Vuelve a GitHub y crea una **release nueva** (por ejemplo `v1.1.3`).
   - Zenodo solo archiva releases creadas **después** de activar el repo.
   - No sirve volver a editar `v1.1.2`; el webhook ya pasó.
5. En Zenodo → **Upload** verás el depósito. Revisa autores, licencia y
   descripción; pulsa **Publish**.
6. Copia el DOI (`10.5281/zenodo...`) y pásalo para integrarlo en
   `CITATION.cff`, README y paper.

Al publicar, Zenodo genera:

- un DOI de **concepto** (siempre apunta a la última versión),
- un DOI de **versión** (específico de `v1.1.3`).

Usa el DOI de concepto en el README y el de versión en `CITATION.cff` si quieres
citar exactamente esta release.

## Opción B — subida manual

1. Descarga el ZIP de la release:
   `https://github.com/Krim10000/Bacterial-DNA-Annotator/releases`
2. Zenodo → **New upload** → **Upload files**.
3. Tipo: **Software**; licencia: **MIT**; acceso: **Open**.
4. Pega la descripción de `.zenodo.json` (o usa el formulario).
5. Añade el enlace a GitHub en **Related identifiers**.
6. Publica y copia el DOI.

## Opción C — API con token

Si prefieres que lo haga yo, necesito un token de Zenodo con permiso de
depósito. No lo pegues en el chat: créalo en
https://zenodo.org/account/settings/applications/tokens/new/ y guárdalo en un
archivo temporal, por ejemplo:

```bash
read -rsp "Zenodo token: " ZENODO_TOKEN; echo
printf '%s' "$ZENODO_TOKEN" > /tmp/zenodo_token
chmod 600 /tmp/zenodo_token
unset ZENODO_TOKEN
```

Después dime la ruta. Publicar un DOI es **irreversible**; antes de pulsar
Publish revisaré contigo título, autoría, licencia y versión.

## Después del DOI

- Añadir `doi:` a `CITATION.cff`.
- Añadir badge `DOI` al README.
- Añadir el DOI al `paper.md` de JOSS.
- Registrar la herramienta en bio.tools.
- Archivar en Software Heritage.

## Nota sobre versiones

`v1.0.0`, `v1.1.0`, `v1.1.1` y `v1.1.2` son releases históricas. Para el DOI
público usa la última versión (`v1.1.3` o superior) y evita citar un asset
antiguo con la afiliación previa.
