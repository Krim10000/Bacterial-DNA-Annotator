"""Integración de ``tss_layer`` en ``pipeline_adapter.run_in_memory``.

Las pruebas usan tablas **locales** de *B. subtilis* y *E. coli* K-12, y un
genoma nuevo sintético. No hay red, descargas ni predicción de posiciones:
solo se comprueba que el pipeline carga las tablas, enlaza las marcas a las
features y deja ``not_available`` cuando no hay datos para el genoma.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from webapp import tss_layer as T
from webapp.pipeline_adapter import run_in_memory

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
SUBTILIS_FASTA = DATA / "dbtbs" / "genome.fna"
SUBTILIS_TABLE = DATA / "dbtbs" / "subtilis_tss_tts.tsv"
K12_FASTA = DATA / "terminators" / "ecoli_k12.fna"

#: Todo apagado: el objetivo es la capa TSS/TTS, no Swiss-Prot/DIAMOND/Rfam.
#: Con ``cds=False`` no se ejecuta pyrodigal; las pruebas de enlace lo encienden.
ALL_LAYERS_OFF = {
    "cds": False,
    "trna": False,
    "rrna": False,
    "ghost": False,
    "deserts": False,
    "esm": False,
    "ncrna": False,
    "taxonomy": False,
    "functional": False,
    "domains": False,
    "context": False,
    "structural_b": False,
    "structural_b_filters": False,
    "complexes": False,
    "regulation": False,
    "membrane": False,
    "oric": False,
    "metabolic": False,
}


def _options(**overrides) -> dict:
    options = dict(ALL_LAYERS_OFF)
    options.update(overrides)
    return options


def _run(tmp_path: Path, filename: str, fasta: bytes,
         options: dict | None = None, data_root: Path | None = None):
    outdir = tmp_path / "out"
    root = data_root if data_root is not None else DATA
    res = run_in_memory(filename, fasta, outdir, root,
                        options=options or _options(), threads=1)
    return res, outdir


def _fasta_prefix(path: Path, bp: int) -> bytes:
    """Primer contig de un FASTA real, recortado a ``bp`` bases."""
    lines = path.read_text(encoding="utf-8").splitlines()
    header = next(line for line in lines if line.startswith(">"))
    chunks: list[str] = []
    length = 0
    for line in lines[lines.index(header) + 1:]:
        if line.startswith(">"):
            break
        chunk = line.strip()
        chunks.append(chunk)
        length += len(chunk)
        if length >= bp:
            break
    return f"{header}\n{''.join(chunks)[:bp]}\n".encode("utf-8")


# --------------------------------------------------------------------------
# Genomas reales: la tabla precalculada se carga tal cual
# --------------------------------------------------------------------------
def test_real_subtilis_loads_table_and_writes_output(tmp_path: Path):
    if not SUBTILIS_FASTA.is_file() or not SUBTILIS_TABLE.is_file():
        pytest.skip("faltan data/dbtbs/genome.fna o subtilis_tss_tts.tsv")

    res, outdir = _run(tmp_path, "subtilis.fna", SUBTILIS_FASTA.read_bytes())

    tss_tts = res["tss_tts"]
    assert tss_tts["status"] == T.STATUS_AVAILABLE
    assert tss_tts["ok"] is True
    assert tss_tts["available"] is True
    # La tabla real de BSGatlas de este repositorio.
    assert tss_tts["n_tss"] == 3390
    assert tss_tts["n_tts"] == 3009
    assert tss_tts["n_total"] == 6399
    assert len(tss_tts["tss"]) == tss_tts["n_tss"]
    assert len(tss_tts["tts"]) == tss_tts["n_tts"]
    assert tss_tts["n_duplicates"] == 0
    assert tss_tts["sources"][0]["source"] == T.SOURCE_BSGATLAS
    assert tss_tts["sources"][0]["path"].endswith("subtilis_tss_tts.tsv")
    assert tss_tts["reason"] == ""
    assert tss_tts["limitation"]
    assert tss_tts["limits"]["predicts_positions"] is False
    assert tss_tts["limits"]["downloads"] is False
    assert tss_tts["attach"]["n_features"] == 0
    assert res["layer_status"]["tss_tts"] == T.STATUS_AVAILABLE

    path = Path(res["outputs"]["tss_tts"])
    assert path.is_file()
    assert path == outdir / "subtilis.tss_tts.tsv"
    lines = path.read_text(encoding="utf-8").splitlines()
    assert lines[0].split("\t")[:4] == ["type", "contig", "pos", "end"]
    assert len(lines) - 1 == tss_tts["n_total"]
    # Solo coordenadas cargadas: la primera fila es el TSS real de dnaA.
    first = lines[1].split("\t")
    assert first[0] == "TSS"
    assert first[1] == "basu168"
    assert first[2] == "157"
    assert first[5] == T.SOURCE_BSGATLAS


def test_real_subtilis_subset_gets_nearest_marks_attached(tmp_path: Path):
    pytest.importorskip("pyrodigal")
    if not SUBTILIS_FASTA.is_file() or not SUBTILIS_TABLE.is_file():
        pytest.skip("faltan data/dbtbs/genome.fna o subtilis_tss_tts.tsv")

    fasta = _fasta_prefix(SUBTILIS_FASTA, 120_000)
    res, _ = _run(tmp_path, "subtilis_subset.fna",
                  fasta, options=_options(cds=True))

    assert res["genome"] == "subtilis_subset"
    features = res["features"]
    assert features, "pyrodigal debería encontrar CDS en el subconjunto real"

    attach = res["tss_tts"]["attach"]
    assert attach["ok"] is True
    assert attach["error"] == ""
    assert attach["n_features"] == len(features)
    assert attach["n_with_tss_upstream"] > 0
    assert attach["n_with_tts_downstream"] > 0
    assert attach["max_distance"] == T.MAX_ATTACH_DISTANCE

    linked = [f for f in features if isinstance(f.get("tss_upstream"), dict)]
    assert linked, "alguna feature debe tener un TSS aguas arriba"
    for feature in linked:
        mark = feature["tss_upstream"]
        start, end = sorted((feature["start"], feature["end"]))
        assert mark["contig"] == feature["contig"]
        assert mark["strand"] == feature["strand"]
        assert feature["tss_upstream_distance"] <= T.MAX_ATTACH_DISTANCE
        if feature["strand"] == "+":
            assert mark["pos"] < start
        else:
            assert mark["pos"] > end
    for feature in features:
        mark = feature.get("tts_downstream")
        if not isinstance(mark, dict):
            continue
        start, end = sorted((feature["start"], feature["end"]))
        assert feature["tts_downstream_distance"] <= T.MAX_ATTACH_DISTANCE
        if feature["strand"] == "+":
            assert mark["pos"] > end
        else:
            assert mark["pos"] < start


def test_real_k12_loads_regulondb_tables(tmp_path: Path):
    if not K12_FASTA.is_file():
        pytest.skip("falta data/terminators/ecoli_k12.fna")
    tss = DATA / "regulondb" / "tss_k12_union.tsv"
    tts = DATA / "regulondb" / "TerminatorSet.tsv"
    if not tss.is_file() or not tts.is_file():
        pytest.skip("faltan las tablas RegulonDB locales")

    res, _ = _run(tmp_path, "ecoli_k12.fna", K12_FASTA.read_bytes())

    tss_tts = res["tss_tts"]
    assert tss_tts["status"] == T.STATUS_AVAILABLE
    assert tss_tts["n_tss"] == 14989
    assert tss_tts["n_tts"] == 370
    assert {source["source"] for source in tss_tts["sources"]} == {
        T.SOURCE_REGULONDB}
    assert res["layer_status"]["tss_tts"] == T.STATUS_AVAILABLE
    assert "tss_tts" in res["outputs"]
    assert tss_tts["data_root"].endswith("data")
    # Los TSS de la unión K-12 no traen contig; la capa lo deja vacío en el
    # TSV en lugar de inventar un cromosoma para poder enlazarlos.
    lines = Path(res["outputs"]["tss_tts"]).read_text(
        encoding="utf-8").splitlines()
    tss_rows = [line.split("\t") for line in lines[1:]
                if line.startswith("TSS\t")]
    assert tss_rows
    assert all(row[1] == "" for row in tss_rows)
    assert tss_tts["attach"]["n_features"] == 0
    assert tss_tts["limitation"]


# --------------------------------------------------------------------------
# Genoma nuevo: ausencia explícita, sin heredar otra especie
# --------------------------------------------------------------------------
def _synthetic_orf_fasta() -> bytes:
    orf = "ATG" + "GCT" * 200 + "TAA"
    spacer = "A" * 60
    return (">contig1 sintético\n"
            + (spacer + orf + spacer) * 2 + "\n").encode("utf-8")


def test_new_synthetic_genome_is_not_available_and_gets_no_marks(
        tmp_path: Path):
    pytest.importorskip("pyrodigal")
    res, outdir = _run(
        tmp_path, "genoma_nuevo_sintetico.fna", _synthetic_orf_fasta(),
        options=_options(cds=True), data_root=DATA,
    )

    tss_tts = res["tss_tts"]
    assert tss_tts["status"] == T.STATUS_NOT_AVAILABLE
    assert tss_tts["available"] is False
    assert tss_tts["ok"] is False
    assert tss_tts["n_tss"] == 0
    assert tss_tts["n_tts"] == 0
    assert tss_tts["n_total"] == 0
    assert tss_tts["tss"] == []
    assert tss_tts["tts"] == []
    assert tss_tts["sources"] == []
    assert "no se predicen" in tss_tts["reason"]
    assert tss_tts["limits"]["cross_genome_transfer"] is False
    assert res["layer_status"]["tss_tts"] == T.STATUS_NOT_AVAILABLE
    assert "tss_tts" not in res["outputs"]
    assert not (outdir / "genoma_nuevo_sintetico.tss_tts.tsv").exists()

    # Las features existen, pero no se les asignó ninguna coordenada nueva.
    assert res["features"]
    for feature in res["features"]:
        assert feature["tss_upstream"] is None
        assert feature["tss_upstream_distance"] is None
        assert feature["tts_downstream"] is None
        assert feature["tts_downstream_distance"] is None
        assert feature["n_tss_inside"] == 0
        assert feature["n_tts_inside"] == 0
    assert tss_tts["attach"]["n_with_tss_upstream"] == 0
    assert tss_tts["attach"]["n_with_tts_downstream"] == 0


# --------------------------------------------------------------------------
# Status parcial y fallo de la capa base
# --------------------------------------------------------------------------
def test_partial_local_table_does_not_invent_the_missing_type(tmp_path: Path):
    data_root = tmp_path / "data"
    (data_root / "dbtbs").mkdir(parents=True)
    (data_root / "dbtbs" / "demo_genoma_tss_tts.tsv").write_text(
        "type\tcontig\tpos\tstrand\tclass\tsource\n"
        "TSS\tcontig1\t100\t+\tA\tdbtbs\n",
        encoding="utf-8",
    )

    res, outdir = _run(
        tmp_path, "demo_genoma.fna", ">contig1\n" + "ACGT" * 100 + "\n",
        data_root=data_root,
    )

    tss_tts = res["tss_tts"]
    assert tss_tts["status"] == T.STATUS_PARTIAL
    assert tss_tts["ok"] is True
    assert tss_tts["n_tss"] == 1
    assert tss_tts["n_tts"] == 0
    assert "no hay TTS" in tss_tts["reason"]
    assert res["layer_status"]["tss_tts"] == T.STATUS_PARTIAL
    path = outdir / "demo_genoma.tss_tts.tsv"
    assert path.is_file()
    assert len(path.read_text(encoding="utf-8").splitlines()) == 2


def test_failed_annotation_keeps_tss_tts_shape(monkeypatch, tmp_path: Path):
    def broken_annotate(*_args, **_kwargs):
        raise RuntimeError("pyrodigal falló")

    monkeypatch.setattr("webapp.pipeline_adapter.annotate", broken_annotate)
    res, _ = _run(tmp_path, "roto.fna", b">contig1\nACGT\n")

    assert res["status"] == "error"
    assert res["tss_tts"]["status"] == "error"
    assert res["tss_tts"]["ok"] is False
    assert res["tss_tts"]["n_tss"] == 0
    assert res["tss_tts"]["n_tts"] == 0
    assert res["layer_status"]["tss_tts"] == "error"
    assert "tss_tts" not in res["outputs"]
