"""Llamada de rRNA con modelos Rfam/Infernal y lectura de resultados.

Motivación (ver ``docs/REFERENCE_PARITY.md``): el ``barrnap`` 0.9 instalado en
esta máquina trae un ``bac.hmm`` sin el modelo 23S y usa por defecto un umbral
``-E 1e-6`` que descarta muchas copias de 5S. La referencia RefSeq/NCBI, en
cambio, anota 16S/23S/5S con Infernal y los modelos Rfam
(RF00177/RF02541/RF00001). Este módulo usa los mismos modelos locales
(``data/rfam_models/bac_rRNA.cm``) para que la llamada del sistema no dependa
del ``barrnap`` del sistema ni de sus umbrales.

Reglas:

- ``cmsearch`` es la fuente preferida; ``barrnap`` queda como respaldo si
  ``cmsearch`` no está, falla o no devuelve ninguna feature.
- Las coordenadas de ``cmsearch --tblout`` y de la salida GFF de ``barrnap``
  son 1-based e inclusivas (formato GFF3). Las coordenadas invertidas de hebra
  negativa se normalizan a ``[menor, mayor]`` conservando la hebra.
- Los duplicados se detectan por contención recíproca (>= 80 % del feature más
  corto) dentro del mismo contig, hebra y tipo (16S/23S/5S). No se fusionan
  tipos distintos ni hebras distintas. Así no se "inventan" ni se colapsan
  features legítimos.
- No se usan GFF de referencia ni coordenadas externas como entrada.
"""
from __future__ import annotations

import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Any, Iterable, Sequence

# Modelos Rfam incluidos en data/rfam_models/bac_rRNA.cm.
MODEL_PRODUCTS: dict[str, tuple[str, str]] = {
    "SSU_rRNA_bacteria": ("16S", "RF00177"),
    "LSU_rRNA_bacteria": ("23S", "RF02541"),
    "5S_rRNA": ("5S", "RF00001"),
}
DEFAULT_RRNA_CM = Path("data") / "rfam_models" / "bac_rRNA.cm"
DEFAULT_CONTAINMENT = 0.8
_RRNA_CLASS_RE = re.compile(r"(16S|23S|5S|18S|28S|5\.8S)", re.IGNORECASE)
_RFAM_CLASS = {
    "RF00177": "16S",
    "RF02541": "23S",
    "RF00001": "5S",
}


# ---------------------------------------------------------------------------
# Modelo de feature
# ---------------------------------------------------------------------------
def _feature(
    *,
    contig: str,
    start: int,
    end: int,
    strand: str,
    product: str,
    source: str,
    score: Any = None,
    evalue: Any = None,
    rfam: str = "",
    model: str = "",
) -> dict[str, Any]:
    start, end = (start, end) if start <= end else (end, start)
    strand = strand if strand in ("+", "-") else "."
    return {
        "contig": contig,
        "type": "rRNA",
        "start": int(start),
        "end": int(end),
        "strand": strand,
        "phase": 0,
        "product": product,
        "source": source,
        "score": score,
        "evalue": evalue,
        "rfam": rfam,
        "model": model,
    }


def rrna_class(feature: dict[str, Any]) -> str:
    """Devuelve 16S/23S/5S (u otra) a partir del modelo, Rfam o producto."""
    model = str(feature.get("model") or "")
    if model in MODEL_PRODUCTS:
        return MODEL_PRODUCTS[model][0]
    rfam = str(feature.get("rfam") or "")
    if rfam in _RFAM_CLASS:
        return _RFAM_CLASS[rfam]
    match = _RRNA_CLASS_RE.search(str(feature.get("product") or ""))
    return match.group(1) if match else "other"


def _product_from_barrnap(text: str) -> str:
    match = re.search(r"(?:^|;)product=([^;]+)", text or "")
    if match:
        return match.group(1).strip()
    match = re.search(r"(?:^|;)Name=([^;]+)", text or "")
    if match:
        value = match.group(1).strip()
        return value.replace("_", " ") if "RNA" in value else value
    return "rRNA"


# ---------------------------------------------------------------------------
# Parsers
# ---------------------------------------------------------------------------
def parse_cmsearch_tblout(text: str) -> list[dict[str, Any]]:
    """Parsea un ``--tblout`` de ``cmsearch`` (1-based, inclusivo).

    Columnas usadas: contig, modelo, acceso Rfam, ``seq from``, ``seq to``,
    hebra, score y E-value. Las líneas ``#`` y las filas incompletas se
    ignoran. Un modelo desconocido se conserva como ``rRNA`` (no se descarta
    silenciosamente un hit real), pero no se le asigna una clase inventada.
    """
    features: list[dict[str, Any]] = []
    for line in (text or "").splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        cols = line.split()
        if len(cols) < 16:
            continue
        contig = cols[0]
        model = cols[2]
        rfam = cols[3]
        try:
            start = int(float(cols[7]))
            end = int(float(cols[8]))
        except (TypeError, ValueError):
            continue
        if start < 1 or end < 1:
            continue
        strand = cols[9] if cols[9] in ("+", "-") else "."
        product, acc = MODEL_PRODUCTS.get(
            model, ("rRNA", rfam or "")
        )
        if model in MODEL_PRODUCTS:
            rfam = acc
        try:
            score = float(cols[13])
        except (TypeError, ValueError):
            score = None
        try:
            evalue = float(cols[15])
        except (TypeError, ValueError):
            evalue = None
        features.append(
            _feature(
                contig=contig,
                start=start,
                end=end,
                strand=strand,
                product=f"{product} ribosomal RNA" if product != "rRNA" else "rRNA",
                source="cmsearch",
                score=score,
                evalue=evalue,
                rfam=rfam,
                model=model,
            )
        )
    return features


def parse_barrnap_gff(gff_text: str) -> list[dict[str, Any]]:
    """Parsea la salida GFF3 de ``barrnap`` (1-based, inclusivo)."""
    features: list[dict[str, Any]] = []
    for line in (gff_text or "").splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        cols = line.rstrip("\n").split("\t")
        if len(cols) < 9 or cols[2] != "rRNA":
            continue
        try:
            start = int(cols[3])
            end = int(cols[4])
        except (TypeError, ValueError):
            continue
        if start < 1 or end < 1:
            continue
        strand = cols[6] if cols[6] in ("+", "-") else "."
        attrs = cols[8]
        product = _product_from_barrnap(attrs)
        name = re.search(r"(?:^|;)Name=([^;]+)", attrs or "")
        model = name.group(1).strip() if name else ""
        features.append(
            _feature(
                contig=cols[0],
                start=start,
                end=end,
                strand=strand,
                product=product,
                source="barrnap",
                score=None,
                evalue=None,
                rfam=MODEL_PRODUCTS.get(model, ("", ""))[1],
                model=model,
            )
        )
    return features


def read_text(path: Path | str) -> str:
    path = Path(path)
    if path.suffix.lower() == ".gz":
        import gzip

        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
            return fh.read()
    return path.read_text(encoding="utf-8", errors="replace")


# ---------------------------------------------------------------------------
# Deduplicación
# ---------------------------------------------------------------------------
def _overlap(left: dict[str, Any], right: dict[str, Any]) -> int:
    return max(
        0,
        min(left["end"], right["end"]) - max(left["start"], right["start"]) + 1,
    )


def _is_duplicate(
    left: dict[str, Any], right: dict[str, Any], containment: float
) -> bool:
    if (
        left["contig"] != right["contig"]
        or left["strand"] != right["strand"]
        or rrna_class(left) != rrna_class(right)
    ):
        return False
    overlap = _overlap(left, right)
    if overlap <= 0:
        return False
    shortest = min(
        left["end"] - left["start"] + 1, right["end"] - right["start"] + 1
    )
    return shortest > 0 and overlap / shortest >= containment


def _preference_key(feature: dict[str, Any]) -> tuple:
    source_rank = 0 if feature.get("source") == "cmsearch" else 1
    score = feature.get("score")
    score_rank = -(float(score) if score is not None else float("-inf"))
    length = feature["end"] - feature["start"] + 1
    return (source_rank, score_rank, -length, feature["contig"], feature["start"])


def dedupe_rrna(
    features: Iterable[dict[str, Any]],
    containment: float = DEFAULT_CONTAINMENT,
) -> list[dict[str, Any]]:
    """Elimina copias solapadas del mismo contig, hebra y tipo.

    Prefiere ``cmsearch`` sobre ``barrnap`` y, dentro de la misma fuente, el
    hit con mayor score y longitud. Los tipos distintos (16S/23S/5S) y las
    hebras distintas nunca se colapsan.
    """
    if not 0.0 < containment <= 1.0:
        raise ValueError("containment debe estar entre 0 y 1")
    kept: list[dict[str, Any]] = []
    for candidate in sorted(features, key=_preference_key):
        if any(_is_duplicate(candidate, other, containment) for other in kept):
            continue
        kept.append(candidate)
    kept.sort(key=lambda f: (f["contig"], f["start"], f["end"], f["strand"]))
    return kept


def assign_ids(features: Sequence[dict[str, Any]]) -> list[dict[str, Any]]:
    """Asigna ``rrna_1..n`` en orden genómico y devuelve los mismos dicts."""
    ordered = sorted(
        features, key=lambda f: (f["contig"], f["start"], f["end"], f["strand"])
    )
    for index, feature in enumerate(ordered, 1):
        feature["id"] = f"rrna_{index}"
    return ordered


# ---------------------------------------------------------------------------
# Selección de modelos del panel CM
# ---------------------------------------------------------------------------
def split_cm_records(text: str) -> list[str]:
    """Divide un archivo CM (INFERNAL1/a) en registros terminados en ``//``."""
    records: list[str] = []
    current: list[str] = []
    for line in (text or "").splitlines(keepends=True):
        current.append(line)
        if line.strip() == "//":
            records.append("".join(current))
            current = []
    if current and any(line.strip() for line in current):
        records.append("".join(current))
    return records


def extract_cm_models(text: str, models: Iterable[str]) -> str:
    """Devuelve un CM con solo los modelos pedidos (por ``NAME``).

    El panel ``bac_rRNA.cm`` contiene registros CM y también bloques HMMER3
    embebidos; aquí solo se conservan los registros ``INFERNAL1/a``. Si no se
    encuentra ningún modelo, devuelve cadena vacía para que el llamador pueda
    degradar sin inventar nada.
    """
    wanted = {str(model) for model in models}
    selected: list[str] = []
    for record in split_cm_records(text):
        match = re.search(r"^NAME\s+(\S+)", record, re.MULTILINE)
        if match and match.group(1) in wanted:
            # Se conservan tanto el registro CM (INFERNAL1/a) como su bloque
            # HMMER3 asociado: cmsearch usa el p7 HMM como filtro y falla si
            # se separa del CM.
            selected.append(record)
    return "".join(selected)


def _cm_models_text(cm_path: Path | str) -> str:
    return Path(cm_path).read_text(encoding="utf-8", errors="replace")


def run_cmsearch_models(
    fasta_path: Path | str,
    cm_path: Path | str,
    models: Iterable[str],
    *,
    cpu: int = 1,
    timeout: int = 1800,
    hmmonly: bool = False,
) -> list[dict[str, Any]]:
    """Ejecuta cmsearch sobre un subconjunto de modelos del panel."""
    text = extract_cm_models(_cm_models_text(cm_path), models)
    if not text.strip():
        return []
    with tempfile.TemporaryDirectory(prefix="rrna_models_") as td:
        subset = Path(td) / "subset.cm"
        subset.write_text(text, encoding="utf-8")
        return run_cmsearch(
            fasta_path,
            subset,
            cpu=cpu,
            timeout=timeout,
            hmmonly=hmmonly,
        )


# ---------------------------------------------------------------------------
# Ejecución de herramientas
# ---------------------------------------------------------------------------
def _run(cmd: Sequence[str], timeout: int) -> subprocess.CompletedProcess:
    return subprocess.run(
        list(cmd),
        check=True,
        capture_output=True,
        text=True,
        timeout=timeout,
        stdin=subprocess.DEVNULL,
    )


def run_cmsearch(
    fasta_path: Path | str,
    cm_path: Path | str,
    *,
    cpu: int = 1,
    timeout: int = 1800,
    hmmonly: bool = False,
) -> list[dict[str, Any]]:
    """Ejecuta ``cmsearch --cut_ga`` y devuelve las features rRNA.

    ``hmmonly=True`` usa solo la parte HMM del modelo (mucho más rápido) para
    localizar copias; los bordes quedan aproximados. El modo por defecto usa
    el CM completo, que reproduce los bordes de Infernal/Rfam.
    """
    fasta_path = Path(fasta_path)
    with tempfile.TemporaryDirectory(prefix="rrna_cm_") as td:
        tbl = Path(td) / "hits.tbl"
        cmd = [
            "cmsearch",
            "--cpu",
            str(max(1, int(cpu))),
            "--cut_ga",
            "--noali",
            "--tblout",
            str(tbl),
            str(cm_path),
            str(fasta_path),
        ]
        if hmmonly:
            cmd.insert(1, "--hmmonly")
        _run(cmd, timeout)
        if not tbl.is_file():
            return []
        return parse_cmsearch_tblout(tbl.read_text(encoding="utf-8"))


def run_barrnap(
    fasta_path: Path | str,
    *,
    kingdom: str = "bac",
    evalue: float = 10.0,
    cpu: int = 1,
    timeout: int = 1800,
) -> list[dict[str, Any]]:
    """Ejecuta ``barrnap`` con un E-value permisivo y parsea su GFF.

    El E-value por defecto de barrnap (1e-6) descarta copias de 5S que sí
    pasan el corte GA de Rfam. ``10`` recupera esas copias sin añadir
    coordenadas: solo acepta hits que el propio modelo reporta.
    """
    cmd = [
        "barrnap",
        "--kingdom",
        kingdom,
        "--threads",
        str(max(1, int(cpu))),
        "--evalue",
        str(float(evalue)),
        str(fasta_path),
    ]
    result = _run(cmd, timeout)
    return parse_barrnap_gff(result.stdout)


def call_rrna_fasta(
    fasta_path: Path | str,
    *,
    cm_path: Path | str | None = None,
    kingdom: str = "bac",
    cpu: int = 1,
    timeout: int = 1800,
    barrnap_evalue: float = 10.0,
    hmmonly: bool = False,
    strategy: str = "exact",
) -> dict[str, Any]:
    """Llama rRNA sobre un FASTA plano con Infernal y respaldo barrnap.

    Estrategias:

    - ``exact``: CM completo para 16S/23S/5S (bordes de Infernal/Rfam).
    - ``balanced``: CM completo para 16S y 5S, HMM para 23S. Detecta todas las
      copias con bordes exactos en 16S/5S y aproximados en 23S.
    - ``recall``: CM completo solo para 5S y HMM para 16S/23S. Es la opción
      rápida para corregir el recall de copias; los bordes de 16S/23S quedan
      aproximados.
    - ``fast``: HMM para los tres modelos (solo localización).

    Devuelve ``{"features", "tool", "reason", "strategy"}``. Nunca inventa
    features: si no hay herramienta, modelo o hits, la lista queda vacía con
    el motivo.
    """
    fasta_path = Path(fasta_path)
    features: list[dict[str, Any]] = []
    tool = "none"
    reasons: list[str] = []
    strategy = strategy if strategy in ("exact", "balanced", "recall", "fast") else "exact"

    def _safe_cmsearch(models: Sequence[str], *, hmm: bool = False) -> list[dict[str, Any]]:
        try:
            if set(models) == set(MODEL_PRODUCTS):
                return run_cmsearch(
                    fasta_path,
                    cm_path,
                    cpu=cpu,
                    timeout=timeout,
                    hmmonly=hmm,
                )
            return run_cmsearch_models(
                fasta_path,
                cm_path,
                models,
                cpu=cpu,
                timeout=timeout,
                hmmonly=hmm,
            )
        except (
            subprocess.CalledProcessError,
            subprocess.TimeoutExpired,
            OSError,
        ) as exc:
            reasons.append(f"cmsearch falló: {type(exc).__name__}")
            return []

    if cm_path is not None and Path(cm_path).is_file() and shutil.which("cmsearch"):
        if strategy == "exact":
            features = _safe_cmsearch(tuple(MODEL_PRODUCTS), hmm=hmmonly)
        elif strategy == "balanced":
            features = _safe_cmsearch(("SSU_rRNA_bacteria", "5S_rRNA"))
            features += _safe_cmsearch(("LSU_rRNA_bacteria",), hmm=True)
        elif strategy == "recall":
            features = _safe_cmsearch(("5S_rRNA",))
            features += _safe_cmsearch(
                ("SSU_rRNA_bacteria", "LSU_rRNA_bacteria"), hmm=True
            )
        else:  # fast
            features = _safe_cmsearch(tuple(MODEL_PRODUCTS), hmm=True)
        if features:
            tool = "cmsearch"
        else:
            reasons.append("cmsearch no devolvió hits")
    elif cm_path is not None:
        reasons.append(f"modelos/cmsearch no disponibles: {cm_path}")
    else:
        reasons.append("sin modelos Rfam")

    if not features and shutil.which("barrnap"):
        try:
            features = run_barrnap(
                fasta_path,
                kingdom=kingdom,
                evalue=barrnap_evalue,
                cpu=cpu,
                timeout=timeout,
            )
            if features:
                tool = "barrnap"
            else:
                reasons.append("barrnap no devolvió hits")
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired, OSError) as exc:
            reasons.append(f"barrnap falló: {type(exc).__name__}")

    features = assign_ids(dedupe_rrna(features))
    return {
        "features": features,
        "tool": tool,
        "reason": "; ".join(reasons),
        "strategy": strategy,
    }


def format_gff_line(feature: dict[str, Any]) -> str:
    """Serializa un hit rRNA como fila GFF3 (1-based, inclusivo)."""
    attrs = f"ID={feature.get('id', 'rrna')};product={feature.get('product', 'rRNA')}"
    if feature.get("rfam"):
        attrs += f";Dbxref=Rfam:{feature['rfam']}"
    if feature.get("source"):
        attrs += f";source={feature['source']}"
    return "\t".join(
        [
            str(feature["contig"]),
            str(feature.get("source") or "rrna"),
            "rRNA",
            str(int(feature["start"])),
            str(int(feature["end"])),
            ".",
            str(feature.get("strand") or "."),
            ".",
            attrs,
        ]
    ) + "\n"


def default_cm_path(project_root: Path | str) -> Path:
    return Path(project_root) / DEFAULT_RRNA_CM


def tool_available(name: str) -> bool:
    return shutil.which(name) is not None


def cpu_count(default: int = 2) -> int:
    try:
        return max(1, min(8, os.cpu_count() or default))
    except (TypeError, ValueError):
        return default
