"""Ficha reutilizable de un gen para la interfaz Streamlit.

La tarjeta **solo muestra lo que el registro declara** (identidad, coordenadas,
A/B/C, confianza, estado, EC y complejo) y deriva dos descargas a partir de la
secuencia disponible: ADN del gen y su traducción con el código genético
estándar. No infiere función, EC, confianza ni estructura.

Convenciones
------------
- ``start``/``end`` son coordenadas 1-based inclusivas, como en GFF3/TSV.
- ``sequences`` es ``{contig: secuencia}`` (el mismo diccionario que devuelve
  ``annotator.annotate``).
- ``gene["pdb"]`` puede ser una ruta absoluta o relativa a ``data_root``/raíz
  del proyecto. Si no existe, la ficha lo dice y no dibuja nada.
- Las etiquetas y explicaciones son bilingües español/inglés. En otros
  idiomas se usa español: no se inventan traducciones nuevas.

La ficha no reemplaza una revisión manual. Una traducción automática no
confirma que el marco de lectura, el inicio ni el producto sean correctos.
"""
from __future__ import annotations

import csv
import gzip
import re
from collections.abc import Mapping
from pathlib import Path
from typing import Any

import streamlit as st

try:  # importación plana cuando la app corre desde webapp/
    from i18n import get_language
except ImportError:  # pragma: no cover - paquete webapp.*
    from webapp.i18n import get_language

try:  # importación plana cuando la app corre desde webapp/
    import pae_utils
    import pae_view
except Exception:  # pragma: no cover
    try:
        from webapp import pae_utils
    except Exception:
        pae_utils = None
        pae_view = None

try:  # importación plana cuando la app corre desde webapp/
    import structure_viewer
except Exception:  # pragma: no cover - depende del entorno
    try:
        from webapp import structure_viewer
    except Exception:
        structure_viewer = None


# --------------------------------------------------------------------------
# Texto bilingüe
# --------------------------------------------------------------------------
def _language() -> str:
    """Idioma activo; ante cualquier fallo usa español (el de la app)."""
    try:
        return get_language()
    except Exception:  # pragma: no cover - sin Streamlit/i18n
        return "es"


def _t(es: str, en: str) -> str:
    """Etiqueta bilingüe: inglés si el idioma activo es ``en``, si no español."""
    return en if _language() == "en" else es


_ABC_EXPLANATION = {
    "A": (
        "Homología Swiss-Prot a una función nombrada. Es evidencia de "
        "homología, no una función confirmada experimentalmente.",
        "Swiss-Prot homology to a named function. This is homology evidence, "
        "not an experimentally confirmed function.",
    ),
    "B": (
        "Homología a una proteína no caracterizada, DUF o UPF. Indica familia, "
        "no una función resuelta.",
        "Homology to an uncharacterized, DUF or UPF protein. It points to a "
        "family, not to a resolved function.",
    ),
    "C": (
        "Sin hit Swiss-Prot. No prueba ausencia de función: la base puede no "
        "tener el homólogo o la proteína puede ser nueva.",
        "No Swiss-Prot hit. It does not prove absence of function: the "
        "database may lack the homolog, or the protein may be novel.",
    ),
    "NOT_EVALUATED": (
        "No evaluado: la búsqueda Swiss-Prot no se ejecutó o falló. No "
        "equivale al camino C.",
        "Not evaluated: the Swiss-Prot search did not run or failed. It is "
        "not equivalent to path C.",
    ),
}

_NOT_EVALUATED = {
    "not_evaluated", "no_evaluado", "no evaluado", "not evaluated",
    "none", "sin evaluar",
}

_CONFIDENCE_EXPLANATION = {
    "high": (
        "Confianza alta según el registro; no equivale a validación "
        "experimental.",
        "High confidence per the record; it is not equivalent to experimental "
        "validation.",
    ),
    "medium": (
        "Confianza media según el registro; conviene revisar la evidencia.",
        "Medium confidence per the record; the evidence should be reviewed.",
    ),
    "low": (
        "Confianza baja según el registro; trátala como hipótesis.",
        "Low confidence per the record; treat it as a hypothesis.",
    ),
    "strong": (
        "Evidencia fuerte según el registro; no confirma la función.",
        "Strong evidence per the record; it does not confirm the function.",
    ),
    "moderate": (
        "Evidencia moderada según el registro; no confirma la función.",
        "Moderate evidence per the record; it does not confirm the function.",
    ),
    "weak": (
        "Evidencia débil según el registro; no confirma la función.",
        "Weak evidence per the record; it does not confirm the function.",
    ),
    "none": (
        "Sin evidencia declarada; no prueba ausencia de función.",
        "No evidence declared; it does not prove absence of function.",
    ),
}

_COMPLEX_EXPLANATION = {
    "complex_suspected": (
        "Sospecha por ortología recíproca/vecindad; no demuestra interacción "
        "física ni estequiometría.",
        "Suspected from reciprocal orthology/neighbourhood; it does not prove "
        "physical interaction or stoichiometry.",
    ),
    "complex_candidate": (
        "Candidato por modelo o vecindad; no demuestra interacción física ni "
        "estequiometría.",
        "Candidate from a model or neighbourhood; it does not prove physical "
        "interaction or stoichiometry.",
    ),
    "complex_ambiguous": (
        "Varias asignaciones posibles; la evidencia no es concluyente.",
        "Several assignments are possible; the evidence is inconclusive.",
    ),
    "not_evaluated": (
        "Complejo no evaluado para este gen.",
        "Complex not evaluated for this gene.",
    ),
}


# --------------------------------------------------------------------------
# Utilidades de datos (puras y testeables)
# --------------------------------------------------------------------------
_COMPLEMENT = str.maketrans({
    "A": "T", "C": "G", "G": "C", "T": "A", "U": "A", "N": "N",
    "R": "Y", "Y": "R", "K": "M", "M": "K", "S": "S", "W": "W",
    "B": "V", "V": "B", "D": "H", "H": "D",
})

_MINUS_STRANDS = {"-1", "-", "−", "minus", "reverse", "r", "rev"}
_PLUS_STRANDS = {"1", "+", "+1", "plus", "forward", "f", "fwd"}

#: Código genético estándar (tabla 1 de NCBI). ``*`` = parada.
_STANDARD_CODE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def extract_gene_dna(seq: str, start: int, end: int,
                     strand: str | int = "+") -> str:
    """Extrae el ADN del gen en coordenadas 1-based inclusivas.

    En hebra negativa devuelve el reverso complementario. Si la secuencia está
    vacía, las coordenadas son inválidas o caen fuera del contig, devuelve
    ``""``: nunca inventa ni rellena bases. Acepta ``+``/``-``, ``1``/``-1`` y
    sinónimos habituales (``plus``/``minus``); cualquier otro valor se trata
    como hebra positiva, que es el valor por defecto.
    """
    if not seq:
        return ""
    clean = "".join(str(seq).split())
    if not clean:
        return ""
    try:
        start_i = int(start)
        end_i = int(end)
    except (TypeError, ValueError):
        try:
            start_i = int(float(start))
            end_i = int(float(end))
        except (TypeError, ValueError):
            return ""
    if start_i < 1 or end_i < start_i or start_i > len(clean) or end_i > len(clean):
        return ""
    fragment = clean[start_i - 1:end_i].upper()
    if _is_minus(strand):
        fragment = fragment.translate(_COMPLEMENT)[::-1]
    return fragment


def translate_dna(seq: str, to_stop: bool = False) -> str:
    """Traduce ADN con el código genético estándar.

    - Acepta minúsculas y ARN (``U`` se trata como ``T``).
    - Un codón con bases ambiguas o desconocidas se traduce a ``X``.
    - Un resto final incompleto (1–2 bases) se ignora, igual que en un CDS
      truncado; la ficha avisa si la longitud no es múltiplo de tres.
    - ``to_stop=True`` se detiene antes del primer codón de parada; por
      defecto la parada se representa como ``*``.
    """
    if not seq:
        return ""
    clean = "".join(str(seq).split()).upper().replace("U", "T")
    protein: list[str] = []
    for i in range(0, len(clean) - 2, 3):
        amino_acid = _STANDARD_CODE.get(clean[i:i + 3], "X")
        if to_stop and amino_acid == "*":
            break
        protein.append(amino_acid)
    return "".join(protein)


def _is_minus(strand: Any) -> bool:
    return str(strand).strip().lower() in _MINUS_STRANDS


def _strand_symbol(strand: Any) -> str:
    """Símbolo de hebra sin adivinar: ``+``, ``−`` o ``—`` si no se declara."""
    value = str(strand).strip().lower() if strand is not None else ""
    if value in _PLUS_STRANDS:
        return "+"
    if value in _MINUS_STRANDS:
        return "−"
    return "—"


def _text(value: Any) -> str:
    """Texto limpio para mostrar; convierte None/NaN/booleanos en ``""``."""
    if value is None or isinstance(value, bool):
        return ""
    text = str(value).strip()
    if text.lower() in {"nan", "none", "null", "<na>", "nat"}:
        return ""
    return text


def _first(gene: Mapping[str, Any], *keys: str) -> str:
    for key in keys:
        value = _text(gene.get(key))
        if value:
            return value
    return ""


def _int_or_none(value: Any) -> int | None:
    if value is None or isinstance(value, bool):
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        try:
            return int(float(value))
        except (TypeError, ValueError):
            return None


def _slug(value: Any) -> str:
    text = re.sub(r"[^0-9A-Za-z_.-]+", "_", str(value or "")).strip("_")
    return text


def _key(prefix: str, gene_id: str, suffix: str) -> str:
    base = "_".join(part for part in (_slug(prefix), _slug(gene_id)) if part)
    return f"{base or 'gene'}_{suffix}"


def _fasta(header: str, sequence: str, width: int = 60) -> str:
    lines = [f">{header}"]
    lines.extend(sequence[i:i + width] for i in range(0, len(sequence), width))
    return "\n".join(lines) + "\n"


def _abc_path(gene: Mapping[str, Any]) -> str:
    """Devuelve A/B/C (o ``not_evaluated``) solo si el registro lo declara.

    Revisa las claves que usan las capas funcionales y, si no hay ninguna,
    deduce el camino de etiquetas como ``A_named_function``. No asigna un
    camino por defecto: un gen sin clasificación devuelve ``""``.
    """
    for key in ("path", "abc", "abc_path", "evidence_path",
                "functional_path", "classification"):
        value = _text(gene.get(key))
        if not value:
            continue
        upper = value.upper()
        if upper in {"A", "B", "C"}:
            return upper
        match = re.match(r"^([ABC])(?:_|$)", value, re.IGNORECASE)
        if match:
            return match.group(1).upper()
        if value.strip().lower() in _NOT_EVALUATED:
            return "not_evaluated"
        return value  # declarado pero no reconocido: se muestra tal cual
    label = _first(gene, "label", "evidence_label", "abc_label")
    match = re.match(r"^([ABC])(?:_|$)", label, re.IGNORECASE)
    return match.group(1).upper() if match else ""


def _abc_explanation(path: str) -> str:
    key = path.strip().lower().replace(" ", "_")
    aliases = {
        "a": "A", "b": "B", "c": "C",
        "not_evaluated": "NOT_EVALUATED", "no_evaluado": "NOT_EVALUATED",
    }
    known = _ABC_EXPLANATION.get(aliases.get(key, ""))
    if known:
        return _t(*known)
    return _t(
        f"Valor declarado en el registro ({path}); la ficha no lo reinterpreta.",
        f"Value declared by the record ({path}); this card does not "
        "reinterpret it.",
    )


def _complex_info(gene: Mapping[str, Any]) -> dict[str, str]:
    """Reúne los campos de complejo presentes, sin completar los ausentes."""
    info: dict[str, str] = {}
    nested = gene.get("complex")
    if isinstance(nested, Mapping):
        aliases = {"name": "complex_name", "id": "complex_id",
                   "status": "complex_status", "role": "complex_role"}
        for key, value in nested.items():
            clean = _text(value)
            if clean:
                info[aliases.get(str(key), str(key))] = clean
    elif _text(nested):
        info["complex_name"] = _text(nested)
    for key in ("complex_id", "complex_name", "complex_status",
                "complex_role", "complex_flag"):
        value = _text(gene.get(key))
        if value:
            info[key] = value
    return info


def _complex_explanation(status: str) -> str:
    key = status.strip().lower().replace(" ", "_")
    if key in _COMPLEX_EXPLANATION:
        return _t(*_COMPLEX_EXPLANATION[key])
    return _t(
        "El complejo se muestra tal como lo declara el registro; no se infiere.",
        "The complex is shown as declared by the record; it is not inferred.",
    )


_REPO_ROOT = Path(__file__).resolve().parents[1]
_STRUCTURE_CATALOG = "structure_catalog.tsv"
_STRUCTURE_SUFFIXES = (
    ".pdb", ".pdb.gz", ".cif", ".cif.gz", ".mmcif", ".mmcif.gz",
)
_UNIPROT_PREFIXES = {"sp", "tr", "gb", "ref", "emb", "dbj", "pdb"}


def _resolve_pdb(value: Any, data_root: Path | str | None) -> Path | None:
    """Resuelve un PDB/mmCIF local; ``None`` si no existe el archivo.

    Acepta ``.pdb``, ``.pdb.gz``, ``.cif``, ``.cif.gz`` (y ``.mmcif``) sin
    descomprimir aquí: el visor y la descarga se encargan del gzip.
    """
    raw = _text(value)
    if not raw:
        return None
    candidate = Path(raw).expanduser()
    if candidate.is_absolute():
        return candidate if candidate.is_file() else None
    roots: list[Path] = []
    if data_root:
        roots.append(Path(data_root))
    roots.append(_REPO_ROOT)
    for root in roots:
        resolved = root / candidate
        if resolved.is_file():
            return resolved
    return None


def _structure_format(path: str | Path) -> str:
    """Devuelve ``"cif"`` o ``"pdb"`` a partir del nombre del archivo."""
    name = Path(path).name.lower()
    if name.endswith((".cif", ".cif.gz", ".mmcif", ".mmcif.gz")):
        return "cif"
    return "pdb"


def _is_gzip(path: str | Path) -> bool:
    return Path(path).name.lower().endswith(".gz")


def _read_structure_bytes(path: str | Path) -> bytes:
    """Lee coordenadas descomprimiendo gzip de forma transparente."""
    p = Path(path)
    if _is_gzip(p):
        with gzip.open(p, "rb") as fh:
            return fh.read()
    return p.read_bytes()


def _accession(gene: Mapping[str, Any]) -> str:
    """Accesión UniProt del gen, normalizando ``sp|ACC|NOMBRE``."""
    for key in ("uniprot", "uniprot_acc", "uniprot_id", "acc", "accession"):
        raw = _text(gene.get(key))
        if not raw:
            continue
        parts = [part.strip() for part in raw.split("|") if part.strip()]
        if not parts:
            continue
        if parts[0].lower() in _UNIPROT_PREFIXES and len(parts) > 1:
            return parts[1].upper()
        return parts[0].upper()
    return ""


def _catalog_files(data_root: Path | str | None) -> list[Path]:
    """Catálogos unificados candidatos, en orden y sin repetir."""
    roots = [Path(data_root)] if data_root else [_REPO_ROOT]
    candidates: list[Path] = []
    for root in roots:
        candidates.extend((
            root / "metabolism" / "targets" / _STRUCTURE_CATALOG,
            root / _STRUCTURE_CATALOG,
            root / "data" / "metabolism" / "targets" / _STRUCTURE_CATALOG,
        ))
    found: list[Path] = []
    seen: set[Path] = set()
    for candidate in candidates:
        try:
            key = candidate.resolve()
        except OSError:
            key = candidate
        if key in seen:
            continue
        seen.add(key)
        if candidate.is_file():
            found.append(candidate)
    return found


def _resolve_catalog_path(raw: Any, catalog: Path,
                          data_root: Path | str | None) -> Path | None:
    """Resuelve la ruta de una fila del catálogo sin asumir una sola raíz."""
    text = _text(raw)
    if not text:
        return None
    candidate = Path(text).expanduser()
    if candidate.is_absolute():
        return candidate if candidate.is_file() else None
    roots: list[Path] = [catalog.parent]
    roots.extend(catalog.parents[1:4])
    if data_root:
        base = Path(data_root)
        roots.extend((base, base.parent))
    roots.append(_REPO_ROOT)
    seen: set[Path] = set()
    for root in roots:
        resolved = root / candidate
        try:
            key = resolved.resolve()
        except OSError:
            key = resolved
        if key in seen:
            continue
        seen.add(key)
        if resolved.is_file():
            return resolved
    return None


def _catalog_structure(accession: str,
                       data_root: Path | str | None) -> Path | None:
    """Busca la accesión en ``structure_catalog.tsv`` y resuelve su archivo.

    Prefiere filas con estado ``ok`` y, dentro de ellas, la fuente AFDB. Una
    fila con archivo existente se usa aunque su estado no sea ``ok``; una ruta
    inexistente se ignora para que la ficha pueda caer a AFDB/ESMFold.
    """
    wanted = (accession or "").strip().upper()
    if not wanted:
        return None
    for catalog in _catalog_files(data_root):
        best_path: Path | None = None
        best_rank: tuple[int, int] | None = None
        try:
            with catalog.open(
                encoding="utf-8-sig", errors="replace", newline=""
            ) as fh:
                for row in csv.DictReader(fh, delimiter="\t"):
                    if _text(row.get("accession")).upper() != wanted:
                        continue
                    path = _resolve_catalog_path(
                        row.get("path") or row.get("pdb"), catalog, data_root)
                    if path is None:
                        continue
                    status = _text(row.get("status")).lower()
                    source = _text(row.get("source")).lower()
                    rank = (
                        0 if status in ("", "ok") else 1,
                        0 if source == "afdb" else 1 if source == "esmfold" else 2,
                    )
                    if best_rank is None or rank < best_rank:
                        best_rank, best_path = rank, path
        except OSError:
            continue
        if best_path is not None:
            return best_path
    return None


def _fallback_structure(accession: str,
                        data_root: Path | str | None) -> Path | None:
    """Caché local AFDB/ESMFold por accesión, sin descargar nada."""
    wanted = (accession or "").strip().upper()
    if not wanted or not data_root:
        return None
    targets = Path(data_root) / "metabolism" / "targets"
    for folder in ("afdb25k", "esmfold_kaggle"):
        for suffix in _STRUCTURE_SUFFIXES:
            candidate = targets / folder / f"{wanted}{suffix}"
            if candidate.is_file():
                return candidate
    return None


def _sequence_for_gene(gene: Mapping[str, Any],
                       sequences: Mapping[str, str] | None
                       ) -> tuple[str, str]:
    """Devuelve ``(adn, motivo)``; ``adn`` vacío si no se puede extraer."""
    if not sequences:
        return "", _t(
            "las secuencias del genoma no están disponibles en esta vista",
            "the genome sequences are not available in this view",
        )
    contig = _text(gene.get("contig"))
    if not contig:
        return "", _t(
            "el gen no declara un contig",
            "the gene does not declare a contig",
        )
    if contig not in sequences:
        return "", _t(
            f"el contig `{contig}` no está entre las secuencias cargadas",
            f"contig `{contig}` is not among the loaded sequences",
        )
    seq = sequences.get(contig)
    seq_text = "".join(str(seq).split()) if seq is not None else ""
    if not seq_text:
        return "", _t(
            f"la secuencia del contig `{contig}` está vacía",
            f"the sequence for contig `{contig}` is empty",
        )
    start = _int_or_none(gene.get("start"))
    end = _int_or_none(gene.get("end"))
    if start is None or end is None:
        return "", _t(
            "faltan coordenadas de inicio/fin válidas",
            "valid start/end coordinates are missing",
        )
    if start < 1 or end < start:
        return "", _t(
            f"las coordenadas {start}–{end} no son válidas",
            f"coordinates {start}–{end} are not valid",
        )
    if end > len(seq_text):
        return "", _t(
            f"las coordenadas {start}–{end} exceden el contig "
            f"({len(seq_text):,} pb)",
            f"coordinates {start}–{end} exceed the contig "
            f"({len(seq_text):,} bp)",
        )
    dna = extract_gene_dna(seq_text, start, end, gene.get("strand", "+"))
    if not dna:
        return "", _t(
            "no se pudo extraer ADN con las coordenadas y la hebra declaradas",
            "DNA could not be extracted with the declared coordinates and strand",
        )
    return dna, ""


# --------------------------------------------------------------------------
# Ficha
# --------------------------------------------------------------------------
def render_gene_card(gene: dict, genome: str, data_root: Path,
                     sequences: dict | None = None,
                     key_prefix: str = "") -> None:
    """Dibuja una ficha de gen reutilizable en Streamlit.

    Parameters
    ----------
    gene:
        Registro del gen. Se leen ``gene_name``/``locus_tag``/``gene_id``,
        ``product``, ``gene_type``/``rna``/``type``, ``contig``, ``start``,
        ``end``, ``strand``, ``path``/``label``, ``confidence``, ``status``,
        ``ec``, ``complex*`` y ``pdb``.
    genome:
        Nombre del genoma, solo para mostrar.
    data_root:
        Raíz de datos; se usa para resolver un ``pdb`` relativo.
    sequences:
        ``{contig: secuencia}``. Si falta, la ficha informa que el ADN y la
        traducción no están disponibles.
    key_prefix:
        Prefijo para claves de widgets; evita colisiones entre fichas.
    """
    if gene is None:
        empty = True
    else:
        try:
            empty = len(gene) == 0
        except (TypeError, ValueError):
            empty = False
    if empty:
        st.info(_t(
            "No hay datos del gen para mostrar.",
            "No gene data to display.",
        ))
        return

    name = _first(gene, "gene_name", "name", "locus_tag", "gene_id", "id")
    locus = _first(gene, "locus_tag", "gene_id")
    gene_id = _first(gene, "gene_id", "locus_tag", "gene_name", "id", "name")
    product = _text(gene.get("product"))
    gene_type = _first(gene, "gene_type", "rna", "type")
    contig = _text(gene.get("contig"))
    start = _int_or_none(gene.get("start"))
    end = _int_or_none(gene.get("end"))

    unit = _t("pb", "bp")
    if start is not None and end is not None:
        span = f"{start:,}–{end:,}"
        coords = f"{contig}:{span}" if contig else span
        length = end - start + 1 if end >= start else 0
        length_text = f"{length:,} {unit}"
    else:
        coords = contig or "—"
        length_text = "—"
    strand_symbol = _strand_symbol(gene.get("strand"))

    st.markdown(f"### 🧬 {name or _t('Gen sin nombre', 'Unnamed gene')}")
    st.caption(
        f"{_t('Genoma', 'Genome')}: `{genome or '—'}` · "
        f"{_t('Coordenadas', 'Coordinates')}: `{coords}` · "
        f"{_t('Longitud', 'Length')}: {length_text}"
    )

    identity = [
        f"- **{_t('Nombre / locus', 'Name / locus')}:** "
        f"`{name or '—'}`" + (f" · `{locus}`" if locus and locus != name else ""),
        f"- **{_t('Producto', 'Product')}:** "
        f"{product or _t('sin anotar', 'not annotated')}",
        f"- **{_t('Tipo', 'Type')}:** {gene_type or '—'}",
        f"- **{_t('Coordenadas', 'Coordinates')}:** `{coords}`",
        f"- **{_t('Hebra', 'Strand')}:** {strand_symbol}",
    ]
    st.markdown("\n".join(identity))

    # ------------------------------------------------ A/B/C + confianza
    st.markdown(f"**{_t('A/B/C y confianza', 'A/B/C and confidence')}**")
    abc = _abc_path(gene)
    label = _first(gene, "label", "evidence_label", "abc_label")
    if abc:
        shown = abc
        if label and label != abc:
            shown = f"{abc} · {label}"
        st.markdown(f"- **A/B/C:** `{shown}`")
        st.caption(_abc_explanation(abc))
    else:
        st.caption(_t(
            "A/B/C no declarado en este registro; no se infiere.",
            "A/B/C is not declared in this record; it is not inferred.",
        ))

    confidence = _text(gene.get("confidence"))
    status = _text(gene.get("status"))
    strength = _text(gene.get("evidence_strength"))
    source = _text(gene.get("source"))
    meta: list[str] = [
        f"{_t('Confianza', 'Confidence')}: `{confidence or '—'}`",
        f"{_t('Estado', 'Status')}: `{status or '—'}`",
    ]
    if strength:
        meta.append(f"{_t('Fuerza de evidencia', 'Evidence strength')}: `{strength}`")
    if source:
        meta.append(f"{_t('Fuente', 'Source')}: `{source}`")
    st.markdown(" · ".join(meta))
    if not confidence or not status:
        st.caption(_t(
            "Confianza/estado no declarados en el registro; no se inventan.",
            "Confidence/status are not declared in the record; they are not "
            "invented.",
        ))
    if confidence:
        conf_key = confidence.strip().lower()
        conf_note = _CONFIDENCE_EXPLANATION.get(
            {"alta": "high", "media": "medium", "baja": "low"}.get(
                conf_key, conf_key))
        if conf_note:
            st.caption(_t(*conf_note))
    st.caption(_t(
        "La confianza y el estado se muestran tal como los declara el registro; "
        "la ficha no los recalcula.",
        "Confidence and status are shown as declared by the record; this card "
        "does not recalculate them.",
    ))

    # ------------------------------------------------------ EC y complejo
    ec = _text(gene.get("ec"))
    complex_info = _complex_info(gene)
    if ec:
        st.markdown(f"- **EC:** `{ec}`")
    if complex_info:
        complex_name = (complex_info.get("complex_name")
                        or complex_info.get("complex_id")
                        or _t("complejo declarado", "declared complex"))
        st.markdown(f"- **{_t('Complejo', 'Complex')}:** {complex_name}")
        details: list[str] = []
        cid = complex_info.get("complex_id", "")
        if cid and cid != complex_name:
            details.append(f"ID: `{cid}`")
        if complex_info.get("complex_status"):
            details.append(
                f"{_t('estado', 'status')}: "
                f"`{complex_info['complex_status']}`"
            )
        if complex_info.get("complex_role"):
            details.append(
                f"{_t('rol', 'role')}: `{complex_info['complex_role']}`"
            )
        if details:
            st.caption(" · ".join(details))
        st.caption(_complex_explanation(
            complex_info.get("complex_status", "")))

    # ---------------------------------------------------------- descargas
    st.markdown(f"**{_t('Descargas', 'Downloads')}**")
    dna, reason = _sequence_for_gene(gene, sequences)
    if dna:
        stem = _slug(gene_id or name or "gene") or "gene"
        fasta_header = " ".join(
            part for part in (gene_id or name, product) if part
        ) or "gene"
        dna_fasta = _fasta(fasta_header, dna)
        st.download_button(
            _t("⬇️ ADN del gen (FASTA)", "⬇️ Gene DNA (FASTA)"),
            data=dna_fasta.encode("utf-8"),
            file_name=f"{stem}.fna",
            mime="text/x-fasta",
            key=_key(key_prefix, gene_id or name, "dna_dl"),
        )
        with st.expander(_t("Ver secuencia de ADN", "View DNA sequence")):
            st.code(dna, language=None)

        protein = translate_dna(dna)
        if protein:
            protein_fasta = _fasta(f"{fasta_header} [translation]", protein)
            st.download_button(
                _t("⬇️ Traducción a proteína (FASTA)",
                   "⬇️ Protein translation (FASTA)"),
                data=protein_fasta.encode("utf-8"),
                file_name=f"{stem}.faa",
                mime="text/x-fasta",
                key=_key(key_prefix, gene_id or name, "protein_dl"),
            )
            with st.expander(_t("Ver traducción", "View translation")):
                st.code(protein, language=None)
            st.caption(_t(
                "Traducción con el código genético estándar (tabla 1). Es "
                "automática y no está validada experimentalmente; no confirma "
                "el marco de lectura, el inicio ni el producto.",
                "Translation with the standard genetic code (table 1). It is "
                "automatic and not experimentally validated; it does not "
                "confirm the reading frame, start or product.",
            ))
            if len(dna) % 3:
                st.caption(_t(
                    f"Longitud de ADN no múltiplo de 3 ({len(dna):,} pb): el "
                    "último codón incompleto se omitió en la traducción.",
                    f"DNA length is not a multiple of 3 ({len(dna):,} bp): the "
                    "incomplete last codon was omitted from the translation.",
                ))
            if "*" in protein[:-1]:
                st.caption(_t(
                    "La traducción contiene un codón de parada interno; "
                    "revisa el marco de lectura y las coordenadas antes de "
                    "usarla.",
                    "The translation contains an internal stop codon; check "
                    "the reading frame and coordinates before using it.",
                ))
    else:
        st.info(_t(
            f"Secuencia de ADN no disponible: {reason}. No se genera "
            "traducción ni se inventa la secuencia.",
            f"DNA sequence not available: {reason}. No translation is "
            "generated and no sequence is invented.",
        ))

    # --------------------------------------------------------------- PDB
    # Orden: (1) ``pdb`` declarado por el registro, (2) catálogo unificado por
    # accesión, (3) caché local AFDB/ESMFold. Nunca se descarga ni se inventa.
    pdb_raw = _text(gene.get("pdb"))
    accession = _accession(gene)
    pdb_path = _resolve_pdb(pdb_raw, data_root) if pdb_raw else None
    if pdb_path is None:
        pdb_path = _catalog_structure(accession, data_root)
    if pdb_path is None:
        pdb_path = _fallback_structure(accession, data_root)
    if pdb_path is None:
        if pdb_raw:
            st.info(_t(
                f"El registro declara `pdb={pdb_raw}`, pero el archivo no está "
                "disponible; no se inventa una estructura.",
                f"The record declares `pdb={pdb_raw}`, but the file is not "
                "available; no structure is invented.",
            ))
        else:
            st.info(_t(
                "No hay PDB local para este gen. El proyecto no descarga una "
                "estructura para cada proteína: usa AFDB/PDB si ya está en la "
                "caché local o ESMFold explícitamente. No se genera ni inventa "
                "una estructura.",
                "No local PDB for this gene. The project does not download a "
                "structure for every protein: it uses AFDB/PDB from the local "
                "cache or explicit ESMFold. No structure is generated or "
                "invented.",
            ))
        return

    try:
        pdb_bytes = _read_structure_bytes(pdb_path)
    except (OSError, EOFError) as exc:
        st.info(_t(
            f"No se pudo leer el PDB `{pdb_path.name}`: {exc}. No se inventa "
            "una estructura.",
            f"The PDB `{pdb_path.name}` could not be read: {exc}. No "
            "structure is invented.",
        ))
        return

    fmt = _structure_format(pdb_path)
    download_name = pdb_path.name
    if _is_gzip(pdb_path):
        download_name = download_name[:-3]
    st.download_button(
        _t("⬇️ Descargar PDB", "⬇️ Download PDB"),
        data=pdb_bytes,
        file_name=download_name,
        mime="chemical/x-cif" if fmt == "cif" else "chemical/x-pdb",
        key=_key(key_prefix, gene_id or name, "pdb_dl"),
    )
    if st.button(
        _t("🧊 Ver estructura 3D", "🧊 View 3D structure"),
        key=_key(key_prefix, gene_id or name, "pdb_view"),
        help=_t(
            "Dibuja las coordenadas del PDB con py3Dmol. No predice ni valida "
            "la función.",
            "Draws the PDB coordinates with py3Dmol. It does not predict or "
            "validate function.",
        ),
    ):
        viewer = (getattr(structure_viewer, "render", None)
                  if structure_viewer is not None else None)
        if viewer is None:
            st.info(_t(
                "Visor 3D no disponible en este entorno; el PDB se puede "
                "descargar igualmente.",
                "3D viewer not available in this environment; the PDB can "
                "still be downloaded.",
            ))
        else:
            viewer(
                str(pdb_path),
                key=_key(key_prefix, gene_id or name, "viewer"),
                caption=_t(
                    "Vista de coordenadas. Una estructura predicha no confirma "
                    "función, sustrato ni interfaz.",
                    "View of coordinates. A predicted structure does not "
                    "confirm function, substrate or interface.",
                ),
            )

    # PAE local de AFDB, si se descargó junto al PDB.
    if pae_utils is not None:
        acc = _text(gene.get("uniprot") or gene.get("acc"))
        pae_path = None
        if gene.get("pae"):
            pae_path = Path(_text(gene.get("pae")))
        elif acc and data_root:
            pae_path = (Path(data_root) / "metabolism" / "targets"
                        / "afdb25k" / f"{acc.split('|')[-1]}.pae.json")
        info = pae_utils.summarize(pae_path) if pae_path else None
        if info and info.get("ok"):
            st.markdown("**PAE (Predicted Aligned Error)**")
            st.caption(
                _t(
                    f"PAE medio {info['mean']:.1f} Å · máximo {info['max']:.1f} Å. "
                    "Es error relativo del modelo, no distancia ni evidencia "
                    "experimental.",
                    f"Mean PAE {info['mean']:.1f} Å · maximum {info['max']:.1f} Å. "
                    "It is relative model error, not a distance or experimental "
                    "evidence.",
                )
            )
            if pae_view is not None and pae_path:
                try:
                    _pae = pae_utils.load_pae(pae_path)
                    _matrix = _pae.get("matrix") or []
                except Exception:  # noqa: BLE001
                    _matrix = []
                if _matrix:
                    with st.expander(_t("Ver matriz PAE", "View PAE matrix")):
                        st.markdown(
                            pae_view.heatmap_svg(_matrix),
                            unsafe_allow_html=True,
                        )
                        st.caption(_t(
                            "Verde = menor error relativo; rojo = mayor. "
                            "No es una distancia ni un mapa de contactos.",
                            "Green = lower relative error; red = higher. "
                            "It is not a distance or a contact map.",
                        ))
            with st.expander(_t("¿Cómo interpretar el PAE?", "How to read PAE")):
                st.markdown(_t(
                    "PAE bajo entre dos regiones indica que su posición "
                    "relativa es más fiable; PAE alto indica que pueden estar "
                    "bien plegadas pero mal orientadas entre sí. No demuestra "
                    "interacción.",
                    "Low PAE between two regions means their relative position "
                    "is more reliable; high PAE means they may be well folded "
                    "but misoriented relative to each other. It does not prove "
                    "interaction.",
                ))
        elif pae_path and Path(pae_path).exists() is False:
            st.caption(_t(
                "PAE no descargado para esta proteína; el PDB se puede usar "
                "igual, pero sin confianza de posición relativa.",
                "PAE not downloaded for this protein; the PDB can still be "
                "used, but without relative-position confidence.",
            ))
