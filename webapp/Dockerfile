FROM python:3.10-slim

# Herramientas externas opcionales (CDS y genes fantasma no las necesitan).
RUN apt-get update && apt-get install -y --no-install-recommends \
        trnascan-se infernal hmmer diamond-aligner curl \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY requirements-annotator.txt /app/requirements.txt
RUN pip install --no-cache-dir -r /app/requirements.txt

COPY .streamlit /app/.streamlit
COPY webapp /app/webapp
COPY scripts/ghost /app/scripts/ghost

EXPOSE 8501
HEALTHCHECK CMD curl -f http://localhost:8501/_stcore/health || exit 1
CMD ["streamlit", "run", "webapp/app.py", "--server.address=0.0.0.0", "--server.port=8501"]
