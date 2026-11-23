# Docker: Windows, macOS y Linux

La vía más sencilla para usar la app en Windows o macOS sin instalar Python ni
las herramientas del pipeline es Docker. La imagen funciona igual en los tres
sistemas.

## Requisitos

- Docker Desktop en Windows o macOS (backend WSL2 en Windows), o Docker Engine
  en Linux.
- 2 GB de RAM libre para la imagen; los datos se montan desde tu carpeta.
- Puerto 8501 libre.

## Arrancar

Desde la raíz del repositorio:

```bash
docker compose up --build
```

Abre:

```text
http://localhost:8501
```

Para detenerlo:

```bash
docker compose down
```

## Datos y resultados

`docker-compose.yml` monta:

```text
./data    → /app/data
./out_fl  → /app/out_fl
```

La imagen **no incluye** los 400 genomas, RegulonDB, KEGG ni Swiss-Prot. Si
tienes esos datos en `data/`, la app los lee. Si no, funcionan la demo, la
anotación de FASTA y los módulos cuyo estado esté disponible.

Para descargar Swiss-Prot y Rhea dentro del contenedor:

```bash
docker compose run --rm annotator python -m webapp.setup_data --data-root /app/data
```

## Herramientas incluidas

| herramienta | incluida | uso |
|---|---|---|
| pyrodigal | ✅ pip | CDS |
| DIAMOND | ✅ apt | homología |
| HMMER | ✅ apt | dominios |
| Infernal | ✅ apt | Rfam |
| tRNAscan-SE | ❌ | tRNA se desactiva |
| barrnap | ❌ | rRNA se desactiva |
| ViennaRNA | ✅ pip | horquillas |

La app detecta lo que falta y desactiva la casilla correspondiente; no falla.

## Limitaciones

- La imagen pesa ~650 MB.
- El pipeline completo con tRNA/rRNA requiere una imagen propia o WSL2/Ubuntu.
- En Apple Silicon, Docker usa emulación amd64 si la imagen no tiene variante
  arm64; para uso intensivo conviene construirla en la arquitectura nativa.
- Los datos restringidos no se redistribuyen; debes montarlos tú.

## Actualizar

```bash
git pull origin v1.0
docker compose build --pull
docker compose up
```

## Alternativa sin Docker

- **Linux/Ubuntu:** `bash install.sh`
- **macOS:** Homebrew + `bash install.sh` (script pendiente)
- **Windows nativo:** WSL2 con Ubuntu y `bash install.sh`
