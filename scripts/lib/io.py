"""Utilidades compartidas: config, FASTA, manifest, logging."""
import gzip, json, logging, hashlib
from pathlib import Path
import yaml

def read_config(path="config.yaml"):
    with open(path) as f:
        return yaml.safe_load(f)

def read_fasta(path):
    """Generador: yield (name, seq). Streaming, no carga todo en RAM."""
    opener = gzip.open if str(path).endswith(".gz") else open
    name, buf = None, []
    with opener(path, "rt") as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(buf).upper()
                name, buf = line[1:].strip().split()[0], []
            else:
                buf.append(line.strip())
    if name is not None:
        yield name, "".join(buf).upper()

def genome_id(fasta_path):
    """ID estable y único por archivo (nombre + hash corto)."""
    p = Path(fasta_path)
    h = hashlib.md5(p.name.encode()).hexdigest()[:8]
    return f"{p.stem.replace('.fna','')}_{h}"

def manifest_append(path, record):
    """Añade una línea JSONL de forma atómica."""
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "a") as f:
        f.write(json.dumps(record) + "\n")

def manifest_load(path):
    """Carga manifest como dict {id: record}."""
    out = {}
    if not Path(path).exists():
        return out
    with open(path) as f:
        for line in f:
            r = json.loads(line)
            out[r.get("seq_id") or r.get("genome_id")] = r
    return out

def setup_logging(name, log_dir):
    Path(log_dir).mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        filename=Path(log_dir) / f"{name}.log",
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        force=True,
    )
    return logging.getLogger(name)
