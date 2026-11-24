"""Comparación de anotaciones GFF3 de dos genomas.

Alcance
-------
- Lee features ``CDS`` y usa la metadata de ``gene``/``mRNA`` (``Parent``,
  ``locus_tag``, ``Name``) solo para agrupar CDS y heredar nombres.
- Empareja genes por ``locus_tag``, después por ``ID`` y, para los que quedan,
  por coordenadas exactas en el mismo contig y hebra.
- No lee ADN, no traduce, no alinea y no hace peticiones de red.

Límite central: compartir ID, ``locus_tag`` o coordenadas no prueba ortología
ni función. Sin secuencia no se puede calcular identidad; cada fila de la
tabla lo declara y la interfaz lo repite en un caption visible.
"""
from __future__ import annotations

import gzip
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from urllib.parse import unquote

try:  # importación normal dentro de la app (webapp/ en sys.path)
    from i18n import tr  # type: ignore
except ImportError:  # importación como paquete: webapp.compare_ui
    try:
        from .i18n import tr  # type: ignore
    except ImportError:  # pragma: no cover - permite usar el módulo aislado
        def tr(text: Any) -> Any:  # type: ignore
            return text


__all__ = [
    "CAPTION_COMPARE",
    "COLUMNS",
    "CompareError",
    "LIMITS",
    "compare_gff",
    "parse_cds",
    "render_compare",
    "to_tsv",
]

#: Texto que debe acompañar siempre a la tabla en la interfaz.
CAPTION_COMPARE = (
    "Comparar coordenadas **no es ortología ni función**: dos genes en la "
    "misma posición pueden ser genes distintos, y dos ortólogos pueden estar "
    "en coordenadas distintas. Sin secuencia no se calcula identidad; esta "
    "tabla compara anotaciones, no biología."
)

#: Advertencias que la interfaz muestra junto al resultado.
LIMITS = (
    "La comparación usa solo la anotación GFF3 (locus_tag, ID, Parent y "
    "coordenadas); no lee ni alinea ADN o proteína y no consulta la red.",
    "No se calcula identidad de secuencia ni similitud. Un gen compartido "
    "por anotación no es necesariamente un ortólogo.",
    "Dos genes con las mismas coordenadas pueden ser distintos; dos ortólogos "
    "pueden tener coordenadas distintas. La coincidencia posicional no prueba "
    "ortología ni función.",
    "Un gen exclusivo puede ser un gen no anotado, un fragmento, un cambio de "
    "coordenadas o una diferencia real entre los genomas.",
    "Si los contigs tienen nombres distintos o los genomas no están alineados, "
    "las coordenadas no son comparables entre sí.",
)

#: Columnas estables de la tabla de diferencias (también las del TSV).
COLUMNS = [
    "estado",
    "coincidencia",
    "gen_a",
    "gen_b",
    "locus_tag_a",
    "locus_tag_b",
    "id_a",
    "id_b",
    "contig_a",
    "contig_b",
    "inicio_a",
    "fin_a",
    "hebra_a",
    "inicio_b",
    "fin_b",
    "hebra_b",
    "bp_a",
    "bp_b",
    "delta_bp",
    "delta_inicio",
    "delta_fin",
    "misma_hebra",
    "producto_a",
    "producto_b",
    "nota",
]

#: Etiqueta de coincidencia -> aclaración que se escribe en ``nota``.
_MATCH_NOTES = {
    "locus_tag": (
        "Mismo locus_tag en la anotación; no se verificó la secuencia ni la "
        "función."
    ),
    "ID": (
        "Mismo ID en la anotación; no se verificó la secuencia ni la función."
    ),
    "coordenadas": (
        "Mismas coordenadas en el mismo contig y hebra; es solo una "
        "coincidencia posicional, no implica ortología ni función."
    ),
}

_STATE_KEY = "compare_ui_result"
_SIG_KEY = "compare_ui_signature"
_UPLOAD_TYPES = ["gff3", "gff", "txt", "gz"]


class CompareError(ValueError):
    """Error legible de lectura o de validación de un GFF3."""


@dataclass
class _Feature:
    seqid: str
    source: str
    ftype: str
    start: int
    end: int
    strand: str
    attrs: dict[str, str]


class _Group:
    """Acumula los CDS que pertenecen al mismo gen dentro de un genoma."""

    def __init__(self, feat: _Feature, meta: dict[str, str]) -> None:
        self.seqid = feat.seqid
        self.source = feat.source
        self.start = feat.start
        self.end = feat.end
        self._strands: set[str] = {feat.strand}
        self.locus_tag = meta["locus_tag"]
        self.gene = meta["gene"]
        self.product = meta["product"]
        self.gene_id = meta["gene_id"]
        self.id = meta["id"]
        self.parent = meta["parent"]
        self.cds_ids: list[str] = []
        self.parts: list[tuple[int, int]] = []
        self.add(feat, meta)

    def add(self, feat: _Feature, meta: dict[str, str]) -> None:
        self.start = min(self.start, feat.start)
        self.end = max(self.end, feat.end)
        self._strands.add(feat.strand)
        cds_id = feat.attrs.get("ID", "")
        if cds_id and cds_id not in self.cds_ids:
            self.cds_ids.append(cds_id)
        self.parts.append((feat.start, feat.end))
        for name in ("locus_tag", "gene", "product", "gene_id", "id", "parent"):
            if not getattr(self, name):
                setattr(self, name, meta[name])

    @property
    def strand(self) -> str:
        return next(iter(self._strands)) if len(self._strands) == 1 else "."

    @property
    def cds_bp(self) -> int:
        return sum(end - start + 1 for start, end in self.parts)

    def to_dict(self) -> dict[str, Any]:
        parts = sorted(self.parts)
        overlap = any(b[0] <= a[1] for a, b in zip(parts, parts[1:]))
        return {
            "seqid": self.seqid,
            "start": self.start,
            "end": self.end,
            "strand": self.strand,
            "id": self.gene_id or self.id,
            "gene_id": self.gene_id,
            "parent": self.parent,
            "locus_tag": self.locus_tag,
            "gene": self.gene,
            "product": self.product,
            "source": self.source,
            "cds_count": len(self.parts),
            "cds_bp": self.cds_bp,
            "cds_ids": list(self.cds_ids),
            "cds_overlap": overlap,
        }


def _read_text(source: Any) -> str:
    """Lee una ruta, un ``Path`` o un archivo subido (texto o bytes)."""
    if isinstance(source, (str, Path)):
        path = Path(source)
        if not path.exists():
            raise CompareError(f"No existe el archivo GFF3: {path}")
        if not path.is_file():
            raise CompareError(f"La ruta no es un archivo: {path}")
        try:
            if path.suffix.lower() == ".gz":
                with gzip.open(path, "rt", encoding="utf-8-sig",
                               errors="replace") as fh:
                    return fh.read()
            return path.read_text(encoding="utf-8-sig", errors="replace")
        except OSError as exc:
            raise CompareError(f"No se pudo leer {path}: {exc}") from exc

    if source is None:
        raise CompareError("Falta la ruta o el archivo GFF3.")
    reader = getattr(source, "read", None)
    if reader is None:
        raise CompareError(
            "Se esperaba una ruta o un archivo subido con método read()."
        )
    seek = getattr(source, "seek", None)
    if callable(seek):
        try:
            seek(0)
        except (OSError, ValueError):
            pass
    try:
        data = reader()
    except Exception as exc:  # noqa: BLE001 - se informa al usuario
        raise CompareError(f"No se pudo leer el archivo GFF3: {exc}") from exc
    if isinstance(data, bytes):
        if data[:2] == b"\x1f\x8b":
            try:
                data = gzip.decompress(data)
            except OSError as exc:
                raise CompareError(f"El archivo gzip está dañado: {exc}") from exc
        return data.decode("utf-8-sig", errors="replace")
    if not isinstance(data, str):
        raise CompareError("El archivo GFF3 no se pudo interpretar como texto.")
    return data


def _decode(value: str) -> str:
    return unquote(value.strip())


def _parse_attributes(raw: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for chunk in raw.split(";"):
        chunk = chunk.strip()
        if not chunk or "=" not in chunk:
            continue
        key, value = chunk.split("=", 1)
        key = _decode(key)
        if key and key not in attrs:
            attrs[key] = _decode(value)
    return attrs


def _attr(attrs: dict[str, str], *names: str) -> str:
    for name in names:
        value = attrs.get(name, "")
        if value:
            return value
    lowered = {key.casefold(): value for key, value in attrs.items()}
    for name in names:
        value = lowered.get(name.casefold(), "")
        if value:
            return value
    return ""


def _first_parent(value: str) -> str:
    return value.split(",")[0].strip() if value else ""


def _iter_features(text: str) -> tuple[list[_Feature], list[str]]:
    features: list[_Feature] = []
    warnings: list[str] = []
    for lineno, line in enumerate(text.splitlines(), start=1):
        if line.strip().startswith("##FASTA"):
            break
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 9:
            warnings.append(
                f"Línea {lineno}: se esperaban 9 columnas y hay "
                f"{len(parts)}; se omite."
            )
            continue
        seqid, source, ftype, start, end, _score, strand, _phase = parts[:8]
        attrs_raw = "\t".join(parts[8:])
        try:
            start_i, end_i = int(start), int(end)
        except ValueError:
            warnings.append(
                f"Línea {lineno}: coordenadas no numéricas; se omite."
            )
            continue
        if start_i < 1 or end_i < start_i:
            warnings.append(
                f"Línea {lineno}: coordenadas inválidas "
                f"({start_i}-{end_i}); se omite."
            )
            continue
        features.append(
            _Feature(
                seqid=seqid.strip(),
                source=source.strip(),
                ftype=ftype.strip(),
                start=start_i,
                end=end_i,
                strand=(strand.strip() or "."),
                attrs=_parse_attributes(attrs_raw),
            )
        )
    return features, warnings


def _resolve_metadata(
    cds: _Feature, by_id: dict[str, _Feature]
) -> dict[str, str]:
    """Hereda locus_tag/Name/product de ``gene``/``mRNA`` vía ``Parent``."""
    attrs = cds.attrs
    locus_tag = _attr(attrs, "locus_tag", "locusTag")
    gene = _attr(attrs, "gene", "Name")
    product = _attr(attrs, "product", "description")
    parent = _first_parent(_attr(attrs, "Parent"))
    gene_id = ""

    ancestors: list[_Feature] = []
    current = parent
    seen: set[str] = set()
    while current and current not in seen:
        seen.add(current)
        ancestor = by_id.get(current)
        if ancestor is None:
            break
        ancestors.append(ancestor)
        locus_tag = locus_tag or _attr(
            ancestor.attrs, "locus_tag", "locusTag"
        )
        gene = gene or _attr(ancestor.attrs, "gene", "Name")
        product = product or _attr(
            ancestor.attrs, "product", "description"
        )
        current = _first_parent(_attr(ancestor.attrs, "Parent"))

    for ancestor in ancestors:
        if ancestor.ftype.casefold() in {"gene", "pseudogene"}:
            gene_id = _attr(ancestor.attrs, "ID")
            break
    if not gene_id and ancestors:
        gene_id = _attr(ancestors[-1].attrs, "ID")

    cds_id = _attr(attrs, "ID")
    return {
        "locus_tag": locus_tag,
        "gene": gene,
        "product": product,
        "gene_id": gene_id,
        "id": gene_id or cds_id,
        "parent": parent,
    }


def _group_key(
    feat: _Feature, meta: dict[str, str]
) -> tuple[tuple[Any, ...], str]:
    if meta["locus_tag"]:
        return ("locus_tag", meta["locus_tag"]), "locus_tag"
    if meta["parent"]:
        return ("parent", meta["parent"]), "parent"
    if meta["id"]:
        return ("id", meta["id"]), "id"
    return ("coords", feat.seqid, feat.start, feat.end, feat.strand), "coords"


def parse_cds(source: Any) -> tuple[list[dict[str, Any]], list[str]]:
    """Lee un GFF3 y devuelve ``(genes, advertencias)``.

    ``genes`` contiene un registro por gen (CDS agrupados por ``locus_tag``,
    ``Parent``, ``ID`` o coordenadas, en ese orden de prioridad). No se
    incluyen secuencias: solo anotación y coordenadas.
    """
    text = _read_text(source)
    features, warnings = _iter_features(text)

    by_id: dict[str, _Feature] = {}
    for feat in features:
        fid = _attr(feat.attrs, "ID")
        if fid and fid not in by_id:
            by_id[fid] = feat

    groups: dict[tuple[Any, ...], _Group] = {}
    seen_rows: set[tuple[Any, ...]] = set()
    without_key = 0
    for feat in features:
        if feat.ftype.casefold() != "cds":
            continue
        row_key = (
            feat.seqid,
            feat.start,
            feat.end,
            feat.strand,
            _attr(feat.attrs, "ID"),
            _attr(feat.attrs, "Parent"),
        )
        if row_key in seen_rows:
            warnings.append(
                f"CDS duplicado en {feat.seqid}:{feat.start}-{feat.end}; "
                "se cuenta una vez."
            )
            continue
        seen_rows.add(row_key)

        meta = _resolve_metadata(feat, by_id)
        key, kind = _group_key(feat, meta)
        if kind == "coords":
            without_key += 1
        group = groups.get(key)
        if group is None:
            groups[key] = _Group(feat, meta)
        else:
            group.add(feat, meta)

    genes: list[dict[str, Any]] = []
    for group in groups.values():
        if len(group._strands) > 1:
            warnings.append(
                f"CDS con hebras distintas en {group.seqid}:"
                f"{group.start}-{group.end}; se usa hebra '.'."
            )
        gene = group.to_dict()
        if gene["cds_overlap"]:
            warnings.append(
                f"CDS solapados dentro de {gene['seqid']}:"
                f"{gene['start']}-{gene['end']}; revisar la anotación."
            )
        genes.append(gene)
    genes.sort(key=lambda g: (g["seqid"], g["start"], g["end"], g["strand"]))

    if without_key:
        warnings.append(
            f"{without_key} CDS sin ID, locus_tag ni Parent: se agruparon "
            "solo por coordenadas exactas dentro del genoma."
        )
    return genes, list(dict.fromkeys(warnings))


def _identity_match(
    genes_a: list[dict[str, Any]],
    genes_b: list[dict[str, Any]],
    field_name: str,
    method: str,
    used_a: set[int],
    used_b: set[int],
    *,
    skip_conflicting_locus_tag: bool = False,
) -> list[tuple[str, dict[str, Any], dict[str, Any]]]:
    """Empareja uno a uno por un campo exacto (locus_tag o ID).

    ``skip_conflicting_locus_tag`` evita usar el ID como identidad cuando ambas
    anotaciones declaran ``locus_tag`` distintos: en ese caso el ID es
    contradictorio y no se debe inventar una pareja.
    """
    index_b: dict[str, list[int]] = {}
    for index, gene in enumerate(genes_b):
        value = gene.get(field_name)
        if value:
            index_b.setdefault(value, []).append(index)

    pairs: list[tuple[str, dict[str, Any], dict[str, Any]]] = []
    for index_a, gene_a in enumerate(genes_a):
        if index_a in used_a:
            continue
        value = gene_a.get(field_name)
        if not value:
            continue
        for match_b in index_b.get(value, []):
            if match_b in used_b:
                continue
            gene_b = genes_b[match_b]
            if skip_conflicting_locus_tag:
                tag_a = gene_a.get("locus_tag") or ""
                tag_b = gene_b.get("locus_tag") or ""
                if tag_a and tag_b and tag_a != tag_b:
                    continue
            used_a.add(index_a)
            used_b.add(match_b)
            pairs.append((method, gene_a, gene_b))
            break
    return pairs


def _coordinate_matches(
    genes_a: list[dict[str, Any]],
    genes_b: list[dict[str, Any]],
    used_a: set[int],
    used_b: set[int],
) -> list[tuple[str, dict[str, Any], dict[str, Any]]]:
    """Empareja por coordenadas exactas y deja constancia del método."""
    index_b: dict[tuple[Any, ...], list[int]] = {}
    for index, gene in enumerate(genes_b):
        if index in used_b:
            continue
        key = (
            gene["seqid"],
            gene["start"],
            gene["end"],
            gene["strand"],
        )
        index_b.setdefault(key, []).append(index)

    pairs: list[tuple[str, dict[str, Any], dict[str, Any]]] = []
    for index_a, gene_a in enumerate(genes_a):
        if index_a in used_a:
            continue
        key = (
            gene_a["seqid"],
            gene_a["start"],
            gene_a["end"],
            gene_a["strand"],
        )
        for match_b in index_b.get(key, []):
            if match_b in used_b:
                continue
            used_a.add(index_a)
            used_b.add(match_b)
            pairs.append(("coordenadas", gene_a, genes_b[match_b]))
            break
    return pairs


def _match_genes(
    genes_a: list[dict[str, Any]], genes_b: list[dict[str, Any]]
) -> tuple[
    list[tuple[str, dict[str, Any], dict[str, Any]]],
    list[dict[str, Any]],
    list[dict[str, Any]],
]:
    used_a: set[int] = set()
    used_b: set[int] = set()
    shared: list[tuple[str, dict[str, Any], dict[str, Any]]] = []

    shared.extend(
        _identity_match(genes_a, genes_b, "locus_tag", "locus_tag",
                        used_a, used_b)
    )
    shared.extend(
        _identity_match(
            genes_a, genes_b, "id", "ID", used_a, used_b,
            skip_conflicting_locus_tag=True,
        )
    )
    shared.extend(
        _coordinate_matches(genes_a, genes_b, used_a, used_b)
    )
    shared.sort(key=lambda p: (
        p[1]["seqid"], p[1]["start"], p[1]["end"], p[1]["strand"]
    ))

    exclusive_a = [
        gene for index, gene in enumerate(genes_a) if index not in used_a
    ]
    exclusive_b = [
        gene for index, gene in enumerate(genes_b) if index not in used_b
    ]
    return shared, exclusive_a, exclusive_b


def _gene_name(gene: dict[str, Any]) -> str:
    return (
        gene.get("locus_tag")
        or gene.get("id")
        or gene.get("gene")
        or f"{gene['seqid']}:{gene['start']}-{gene['end']}"
    )


def _blank_row() -> dict[str, Any]:
    return {column: None for column in COLUMNS}


def _row_for(
    estado: str,
    match: str | None,
    gene_a: dict[str, Any] | None,
    gene_b: dict[str, Any] | None,
    label_a: str,
    label_b: str,
) -> dict[str, Any]:
    row = _blank_row()
    row["estado"] = estado
    row["coincidencia"] = match
    if gene_a is not None:
        row["gen_a"] = _gene_name(gene_a)
        row["locus_tag_a"] = gene_a["locus_tag"] or None
        row["id_a"] = gene_a["id"] or None
        row["contig_a"] = gene_a["seqid"]
        row["inicio_a"] = gene_a["start"]
        row["fin_a"] = gene_a["end"]
        row["hebra_a"] = gene_a["strand"]
        row["bp_a"] = gene_a["cds_bp"]
        row["producto_a"] = gene_a["product"] or None
    if gene_b is not None:
        row["gen_b"] = _gene_name(gene_b)
        row["locus_tag_b"] = gene_b["locus_tag"] or None
        row["id_b"] = gene_b["id"] or None
        row["contig_b"] = gene_b["seqid"]
        row["inicio_b"] = gene_b["start"]
        row["fin_b"] = gene_b["end"]
        row["hebra_b"] = gene_b["strand"]
        row["bp_b"] = gene_b["cds_bp"]
        row["producto_b"] = gene_b["product"] or None
    if match and gene_a is not None and gene_b is not None:
        row["delta_bp"] = gene_b["cds_bp"] - gene_a["cds_bp"]
        row["delta_inicio"] = gene_b["start"] - gene_a["start"]
        row["delta_fin"] = gene_b["end"] - gene_a["end"]
        if {gene_a["strand"], gene_b["strand"]} <= {"+", "-"}:
            row["misma_hebra"] = gene_a["strand"] == gene_b["strand"]
        row["nota"] = _MATCH_NOTES[match]
    elif gene_a is not None:
        row["nota"] = (
            f"Sin pareja en {label_b} por locus_tag, ID ni coordenadas "
            "exactas."
        )
    else:
        row["nota"] = (
            f"Sin pareja en {label_a} por locus_tag, ID ni coordenadas "
            "exactas."
        )
    return row


def _rows(
    shared: list[tuple[str, dict[str, Any], dict[str, Any]]],
    exclusive_a: list[dict[str, Any]],
    exclusive_b: list[dict[str, Any]],
    label_a: str,
    label_b: str,
) -> list[dict[str, Any]]:
    rows = [
        _row_for("compartido", method, gene_a, gene_b, label_a, label_b)
        for method, gene_a, gene_b in shared
    ]
    rows.extend(
        _row_for(f"solo {label_a}", None, gene, None, label_a, label_b)
        for gene in exclusive_a
    )
    rows.extend(
        _row_for(f"solo {label_b}", None, None, gene, label_a, label_b)
        for gene in exclusive_b
    )
    return rows


def _table(rows: list[dict[str, Any]]) -> Any:
    import pandas as pd

    # ``convert_dtypes`` evita que una columna numérica con huecos se vuelva
    # float (2.0 en el TSV) y conserva los vacíos como NA.
    return pd.DataFrame(rows, columns=COLUMNS).convert_dtypes()


def _source_signature(source: Any) -> str:
    if source is None:
        return ""
    if isinstance(source, (str, Path)):
        return str(source)
    name = getattr(source, "name", "") or ""
    size = getattr(source, "size", "") or ""
    return f"upload:{name}:{size}"


def compare_gff(
    path_a: Any,
    path_b: Any,
    label_a: str = "A",
    label_b: str = "B",
) -> dict[str, Any]:
    """Compara dos GFF3 y devuelve conteos, parejas y tabla de diferencias.

    Parameters
    ----------
    path_a, path_b:
        Ruta (``str``/``Path``) o archivo con ``read()`` que devuelva texto o
        bytes. Se aceptan rutas/archivos ``.gz``.
    label_a, label_b:
        Etiquetas visibles para los genomas. Deben ser distintas.

    Returns
    -------
    dict
        ``labels``, ``counts`` (incluye ``genes_por_genoma``), ``genes_a``,
        ``genes_b``, ``shared``, ``exclusive_a``, ``exclusive_b``, ``table``
        (``pandas.DataFrame`` con ``COLUMNS``), ``warnings`` y ``limits``.

    Notes
    -----
    El emparejamiento es exacto y en cascada: ``locus_tag``, luego ``ID`` y
    luego coordenadas (``seqid`` + inicio + fin + hebra). No se calcula
    identidad de secuencia ni se consulta la red.
    """
    label_a = str(label_a).strip() or "A"
    label_b = str(label_b).strip() or "B"
    if label_a == label_b:
        raise CompareError(
            "Las etiquetas de los genomas deben ser distintas."
        )

    genes_a, warnings_a = parse_cds(path_a)
    genes_b, warnings_b = parse_cds(path_b)
    shared, exclusive_a, exclusive_b = _match_genes(genes_a, genes_b)

    by_method = {"locus_tag": 0, "ID": 0, "coordenadas": 0}
    for method, _gene_a, _gene_b in shared:
        by_method[method] += 1

    counts = {
        "genes_a": len(genes_a),
        "genes_b": len(genes_b),
        "cds_a": sum(gene["cds_count"] for gene in genes_a),
        "cds_b": sum(gene["cds_count"] for gene in genes_b),
        "compartidos": len(shared),
        "exclusivos_a": len(exclusive_a),
        "exclusivos_b": len(exclusive_b),
        "compartidos_locus_tag": by_method["locus_tag"],
        "compartidos_id": by_method["ID"],
        "compartidos_coordenadas": by_method["coordenadas"],
        "genes_por_genoma": {
            label_a: len(genes_a),
            label_b: len(genes_b),
        },
    }

    warnings = [f"{label_a}: {w}" for w in warnings_a]
    warnings.extend(f"{label_b}: {w}" for w in warnings_b)

    if by_method["coordenadas"]:
        warnings.append(
            f"{by_method['coordenadas']} genes se emparejaron solo por "
            "coordenadas exactas. Es una coincidencia posicional: no "
            "demuestra ortología ni función."
        )
    if by_method["ID"]:
        warnings.append(
            f"{by_method['ID']} genes se emparejaron por ID de anotación. "
            "Un ID puede ser un contador local (por ejemplo, cds_1) y no "
            "implica el mismo gen; no se verificó la secuencia."
        )
    if genes_a and genes_b and not shared:
        warnings.append(
            "No hay genes compartidos con las claves disponibles "
            "(locus_tag, ID o coordenadas exactas)."
        )
    if not genes_a:
        warnings.append(f"{label_a}: el GFF3 no contiene CDS.")
    if not genes_b:
        warnings.append(f"{label_b}: el GFF3 no contiene CDS.")

    seqids_a = {gene["seqid"] for gene in genes_a}
    seqids_b = {gene["seqid"] for gene in genes_b}
    if genes_a and genes_b and not (seqids_a & seqids_b):
        warnings.append(
            "Los genomas no comparten nombres de contig; la coincidencia por "
            "coordenadas solo se evalúa dentro del mismo contig, así que no "
            "puede haberla entre contigs distintos."
        )

    rows = _rows(shared, exclusive_a, exclusive_b, label_a, label_b)
    return {
        "labels": {"a": label_a, "b": label_b},
        "counts": counts,
        "genes_a": genes_a,
        "genes_b": genes_b,
        "shared": [
            {"metodo": method, "a": gene_a, "b": gene_b}
            for method, gene_a, gene_b in shared
        ],
        "exclusive_a": exclusive_a,
        "exclusive_b": exclusive_b,
        "table": _table(rows),
        "warnings": list(dict.fromkeys(warnings)),
        "limits": list(LIMITS),
    }


def to_tsv(result: dict[str, Any]) -> str:
    """Serializa la tabla de ``compare_gff`` como TSV con encabezado."""
    table = result.get("table")
    if table is None:
        raise CompareError("El resultado no contiene tabla para exportar.")
    return table.to_csv(sep="\t", index=False, lineterminator="\n")


def _render_result(st: Any, result: dict[str, Any]) -> None:
    counts = result["counts"]
    labels = result["labels"]
    label_a, label_b = labels["a"], labels["b"]

    m1, m2, m3 = st.columns(3)
    m1.metric(
        tr("Genes {genoma}").format(genoma=label_a), counts["genes_a"]
    )
    m2.metric(tr("Genes compartidos"), counts["compartidos"])
    m3.metric(
        tr("Genes {genoma}").format(genoma=label_b), counts["genes_b"]
    )
    m4, m5 = st.columns(2)
    m4.metric(
        tr("Exclusivos {genoma}").format(genoma=label_a),
        counts["exclusivos_a"],
    )
    m5.metric(
        tr("Exclusivos {genoma}").format(genoma=label_b),
        counts["exclusivos_b"],
    )
    st.caption(
        tr(
            "Coincidencias por método: locus_tag={locus_tag}, ID={id}, "
            "coordenadas={coordenadas}."
        ).format(
            locus_tag=counts["compartidos_locus_tag"],
            id=counts["compartidos_id"],
            coordenadas=counts["compartidos_coordenadas"],
        )
    )

    st.dataframe(result["table"], hide_index=True, width="stretch")
    st.caption(tr(CAPTION_COMPARE))

    safe = re.sub(
        r"[^0-9A-Za-z._-]+", "_", f"{label_a}_vs_{label_b}"
    ).strip("_") or "comparacion"
    st.download_button(
        tr("⬇️ Descargar tabla de diferencias (TSV)"),
        data=to_tsv(result),
        file_name=f"{safe}_diferencias.tsv",
        mime="text/tab-separated-values",
        key="compare_download",
    )

    if result["warnings"]:
        with st.expander(
            tr("Advertencias ({n})").format(n=len(result["warnings"]))
        ):
            for warning in result["warnings"]:
                st.warning(warning)

    with st.expander(tr("Qué no prueba esta comparación")):
        for limit in result["limits"]:
            st.write("• " + tr(limit))


def render_compare() -> None:
    """Renderiza en Streamlit dos GFF3, la tabla y la descarga TSV."""
    import streamlit as st

    st.subheader("🧬 " + tr("Comparar dos genomas"))
    st.caption(tr(CAPTION_COMPARE))
    st.caption(tr(
        "Sube un GFF3 por genoma o escribe su ruta local. Si hay archivo "
        "subido, se usa ese; la ruta se usa solo cuando no hay archivo. "
        "Formatos: .gff3, .gff, .txt y .gz."
    ))

    col_a, col_b = st.columns(2)
    with col_a:
        label_a = (
            st.text_input(tr("Etiqueta del genoma A"), value="A",
                          key="compare_label_a") or "A"
        ).strip() or "A"
        upload_a = st.file_uploader(
            tr("GFF3 del genoma A"), type=_UPLOAD_TYPES, key="compare_file_a"
        )
        path_a = (
            st.text_input(
                tr("… o ruta local al GFF3 del genoma A"),
                value="",
                key="compare_path_a",
                placeholder="/ruta/genoma_a.gff3",
            ) or ""
        ).strip()
        source_a = upload_a if upload_a is not None else (path_a or None)
    with col_b:
        label_b = (
            st.text_input(tr("Etiqueta del genoma B"), value="B",
                          key="compare_label_b") or "B"
        ).strip() or "B"
        upload_b = st.file_uploader(
            tr("GFF3 del genoma B"), type=_UPLOAD_TYPES, key="compare_file_b"
        )
        path_b = (
            st.text_input(
                tr("… o ruta local al GFF3 del genoma B"),
                value="",
                key="compare_path_b",
                placeholder="/ruta/genoma_b.gff3",
            ) or ""
        ).strip()
        source_b = upload_b if upload_b is not None else (path_b or None)

    signature = (
        label_a,
        label_b,
        _source_signature(source_a),
        _source_signature(source_b),
    )
    if st.session_state.get(_SIG_KEY) != signature:
        st.session_state[_SIG_KEY] = signature
        st.session_state.pop(_STATE_KEY, None)

    if st.button(tr("Comparar GFF3"), type="primary", key="compare_run"):
        if source_a is None or source_b is None:
            st.session_state.pop(_STATE_KEY, None)
            st.warning(tr(
                "Selecciona un GFF3 para A y otro para B, como archivo o "
                "como ruta."
            ))
        elif label_a == label_b:
            st.session_state.pop(_STATE_KEY, None)
            st.error(tr("Las etiquetas de los genomas deben ser distintas."))
        else:
            try:
                with st.spinner(tr("Comparando anotaciones…")):
                    result = compare_gff(
                        source_a, source_b, label_a=label_a, label_b=label_b
                    )
            except Exception as exc:  # noqa: BLE001 - se muestra al usuario
                st.session_state.pop(_STATE_KEY, None)
                st.error(
                    tr("No se pudo comparar: {error}").format(error=exc)
                )
            else:
                st.session_state[_STATE_KEY] = result

    result = st.session_state.get(_STATE_KEY)
    if result is None:
        st.info(tr(
            "Cuando compares, verás cuántos genes tiene cada genoma, cuáles "
            "comparten locus_tag/ID/coordenadas y cuáles son exclusivos. "
            "La tabla se puede descargar como TSV."
        ))
    else:
        _render_result(st, result)
