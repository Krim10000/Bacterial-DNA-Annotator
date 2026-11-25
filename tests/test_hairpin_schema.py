"""Tests de ``migrate_hairpin_schema`` (strand antiguo -> location + hebra).

Se usan TSV sintéticos: ni red, ni datos reales, ni dependencias externas.
Cubren hebras ``+``/``-``, herencia del gen izquierdo/derecho, ausencia de gen
flanqueante, idempotencia, ``--dry-run``, ``excluded/`` y fin de línea CRLF.
"""
from __future__ import annotations

import csv
import importlib.util
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "regulation" / "migrate_hairpin_schema.py"

HP_FIELDS = [
    "genome", "contig", "start", "end", "strand", "class", "stem", "gc",
    "loop", "polyU", "dG", "gap", "left_gene", "right_gene",
]
GENE_FIELDS = [
    "genome", "contig", "start", "end", "strand", "gene_id", "locus_tag",
    "gene_name",
]


def _load_module():
    spec = importlib.util.spec_from_file_location(
        "migrate_hairpin_schema_test", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module       # dataclasses necesita el módulo vivo
    spec.loader.exec_module(module)
    return module


migrate = _load_module()


def _write_tsv(path: Path, fields: list[str], rows: list[dict],
               newline: str = "\n") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t",
                                lineterminator=newline)
        writer.writeheader()
        writer.writerows(rows)


def _read_tsv(path: Path) -> tuple[list[str], list[dict]]:
    with path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def _gene(genome: str, start: int, strand: str, gene_id: str,
          locus_tag: str, gene_name: str) -> dict:
    return {
        "genome": genome, "contig": "chr", "start": str(start),
        "end": str(start + 99), "strand": strand, "gene_id": gene_id,
        "locus_tag": locus_tag, "gene_name": gene_name,
    }


def _genes_g1() -> list[dict]:
    return [
        _gene("g1", 100, "1", "idL", "gL", "alpha"),
        _gene("g1", 400, "-1", "idR", "gR", "beta"),
        _gene("g1", 700, "-1", "idR2", "gR2", "gamma"),
    ]


def _hairpin(start: int, left: str, right: str, cls: str) -> dict:
    return {
        "genome": "g1", "contig": "chr", "start": str(start),
        "end": str(start + 20), "strand": "intergénica", "class": cls,
        "stem": "12", "gc": "0.67", "loop": "4", "polyU": "9", "dG": "25.0",
        "gap": str(start - 1), "left_gene": left, "right_gene": right,
    }


def _hairpins_g1() -> list[dict]:
    return [
        _hairpin(250, "gL", "gR", "A_fuerte"),       # izquierdo gana
        _hairpin(350, "", "gR", "B_media"),          # solo derecho
        _hairpin(550, "no-existe", "gR2", "B_media"),  # izquierdo no está
        _hairpin(950, "", "", "C_debil"),            # sin gen flanqueante
        _hairpin(1050, "gL", "gR", "A_fuerte"),      # izquierdo +, derecho -
        _hairpin(1150, "idL", "", "A_fuerte"),       # identificado por gene_id
        _hairpin(1250, "beta", "", "B_media"),       # identificado por gene_name
    ]


def _tree(tmp_path: Path, rows: list[dict] | None = None) -> tuple[Path, Path]:
    base = tmp_path / "operon_finder"
    genes_root = tmp_path / "regulation"
    _write_tsv(genes_root / "g1" / "genes.tsv", GENE_FIELDS, _genes_g1())
    _write_tsv(base / "g1_hairpins.tsv", HP_FIELDS,
               _hairpins_g1() if rows is None else rows)
    return base, genes_root


def _run(base: Path, genes_root: Path, *extra: str) -> int:
    return migrate.main(["--dir", str(base), "--genes-root", str(genes_root),
                         *extra])


def test_migra_hebras_y_conserva_columnas(tmp_path: Path):
    base, genes_root = _tree(tmp_path)
    original = _hairpins_g1()

    assert _run(base, genes_root) == 0

    fields, rows = _read_tsv(base / "g1_hairpins.tsv")
    assert fields == HP_FIELDS + ["location", "strand_source"]
    assert len(rows) == len(original)

    expected = {
        "250": ("+", "gen izquierdo"),
        "350": ("-", "gen derecho"),
        "550": ("-", "gen derecho"),
        "950": ("", "sin gen flanqueante"),
        "1050": ("+", "gen izquierdo"),
        "1150": ("+", "gen izquierdo"),
        "1250": ("-", "gen izquierdo"),
    }
    by_start = {row["start"]: row for row in rows}
    assert set(by_start) == set(expected)
    for start, (strand, source) in expected.items():
        row = by_start[start]
        assert row["strand"] == strand
        assert row["strand_source"] == source
        assert row["location"] == "intergénica"

    # El resto de columnas y valores no se toca.
    old_by_start = {row["start"]: row for row in original}
    for start, new_row in by_start.items():
        for field in HP_FIELDS:
            if field != "strand":
                assert new_row[field] == old_by_start[start][field]


def test_idempotente(tmp_path: Path, capsys):
    base, genes_root = _tree(tmp_path)
    path = base / "g1_hairpins.tsv"

    assert _run(base, genes_root) == 0
    first = path.read_bytes()
    mtime = path.stat().st_mtime_ns
    assert b"location" in first and b"strand_source" in first

    assert _run(base, genes_root) == 0
    second = path.read_bytes()
    out = capsys.readouterr().out
    assert second == first
    assert path.stat().st_mtime_ns == mtime
    assert "0 reescritos" in out
    assert "0 cambiadas" in out
    assert not list(base.glob("*.tmp")) and not list(base.glob(".*.tmp"))


def test_dry_run_no_escribe(tmp_path: Path, capsys):
    base, genes_root = _tree(tmp_path)
    path = base / "g1_hairpins.tsv"
    before = path.read_bytes()
    before_files = sorted(p.name for p in base.iterdir())

    assert _run(base, genes_root, "--dry-run") == 0
    out = capsys.readouterr().out

    assert path.read_bytes() == before
    assert sorted(p.name for p in base.iterdir()) == before_files
    assert "dry-run" in out
    assert "7 total, 7 a cambiar" in out
    assert "hebra: + 3, - 3, sin hebra 1" in out
    assert ("fuente: gen izquierdo 4, gen derecho 2, "
            "sin gen flanqueante 1") in out


def test_incluye_excluded_y_no_inventa_sin_genes_tsv(tmp_path: Path, capsys):
    base, genes_root = _tree(tmp_path)
    excluded = base / "excluded"

    # g2 sí tiene tabla de genes: la hebra sale de su gen izquierdo.
    _write_tsv(genes_root / "g2" / "genes.tsv", GENE_FIELDS,
               [_gene("g2", 100, "-1", "idL2", "gL2", "delta")])
    g2_row = dict(_hairpin(250, "gL2", "", "B_media"), genome="g2")
    _write_tsv(excluded / "g2_hairpins.tsv", HP_FIELDS, [g2_row])

    # g3 no tiene genes.tsv: no se inventa hebra.
    g3_row = dict(_hairpin(250, "cualquiera", "otro", "C_debil"), genome="g3")
    _write_tsv(excluded / "g3_hairpins.tsv", HP_FIELDS, [g3_row])

    assert _run(base, genes_root) == 0
    out = capsys.readouterr().out

    _, rows2 = _read_tsv(excluded / "g2_hairpins.tsv")
    assert rows2[0]["strand"] == "-"
    assert rows2[0]["strand_source"] == "gen izquierdo"

    _, rows3 = _read_tsv(excluded / "g3_hairpins.tsv")
    assert rows3[0]["strand"] == ""
    assert rows3[0]["strand_source"] == "sin gen flanqueante"
    assert rows3[0]["location"] == "intergénica"

    assert "excluded" in out
    assert "principal: 1 archivos" in out
    assert "excluded: 2 archivos" in out
    assert "genomas sin genes.tsv (1): g3" in out


def test_gene_name_ambiguo_no_resuelve(tmp_path: Path):
    genes_root = tmp_path / "regulation"
    base = tmp_path / "operon_finder"
    _write_tsv(genes_root / "g4" / "genes.tsv", GENE_FIELDS, [
        _gene("g4", 100, "1", "idP", "gP", "dup"),
        _gene("g4", 400, "-1", "idM", "gM", "dup"),
    ])
    row = dict(_hairpin(250, "dup", "", "B_media"), genome="g4")
    _write_tsv(base / "g4_hairpins.tsv", HP_FIELDS, [row])

    assert _run(base, genes_root) == 0
    _, rows = _read_tsv(base / "g4_hairpins.tsv")
    assert rows[0]["strand"] == ""
    assert rows[0]["strand_source"] == "sin gen flanqueante"


def test_conserva_fin_de_linea_crlf(tmp_path: Path):
    base, genes_root = _tree(tmp_path)
    path = base / "g1_hairpins.tsv"
    _write_tsv(path, HP_FIELDS, _hairpins_g1(), newline="\r\n")
    assert b"\r\n" in path.read_bytes()

    assert _run(base, genes_root) == 0

    raw = path.read_bytes()
    assert b"\r\n" in raw
    assert b"\n" not in raw.replace(b"\r\n", b"")


def test_normalize_strand():
    assert migrate.normalize_strand("1") == "+"
    assert migrate.normalize_strand("-1") == "-"
    assert migrate.normalize_strand("+") == "+"
    assert migrate.normalize_strand("−") == "-"
    assert migrate.normalize_strand("intergénica") == ""
    assert migrate.normalize_strand("") == ""
