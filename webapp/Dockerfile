FROM python:3.10-slim-bookworm

# Herramientas externas disponibles en Debian. tRNAscan-SE y barrnap no están
# empaquetados aquí; el núcleo (CDS, ghost, GUI) funciona sin ellos y la app
# desactiva las casillas que los necesitan.
RUN apt-get update && apt-get install -y --no-install-recommends \
        infernal hmmer diamond-aligner curl \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY requirements-annotator.txt /app/requirements.txt
RUN pip install --no-cache-dir -r /app/requirements.txt

COPY .streamlit /app/.streamlit
COPY webapp /app/webapp
COPY modules /app/modules
COPY scripts/ghost /app/scripts/ghost

ENV DNA_DATA_ROOT=/app/data \
    DNA_OUT_ROOT=/app/out_fl
VOLUME ["/app/data", "/app/out_fl"]

EXPOSE 8501
HEALTHCHECK CMD curl -f http://localhost:8501/_stcore/health || exit 1
CMD ["streamlit", "run", "webapp/app.py", "--server.address=0.0.0.0", "--server.port=8501"]
