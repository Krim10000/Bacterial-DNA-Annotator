"""Pruebas del lote del modelo genómico K-12 para los 400 genomas.

No se ejecuta el lote completo: se prueban descubrimiento, conversión de
coordenadas, orquestación (con un runner falso) y una integración mínima real
con `build_candidates.py --apply-only`.
"""
from __future__ import annotations

import csv
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))

import run_complex_model_400 as M  # noqa: E402

GENES_COLS = ["genome", "contig", "start", "end", "strand", "gene_id",
              "locus_tag", "gene_name", "ec", "product"]
CLUSTER_COLS = ["cluster_id", "genome", "n_genes", "genes", "gene_names",
                "mean_score", "min_score", "confidence", "evidence"]


def _write_genes(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=GENES_COLS,
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def _write_clusters(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=CLUSTER_COLS,
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def _gene(genome: str, tag: str, name: str, contig: str, start: int,
          end: int) -> dict:
    return {"genome": genome, "contig": contig, "start": start, "end": end,
            "strand": "1", "gene_id": tag, "locus_tag": "", "gene_name": name,
            "ec": "", "product": ""}


def _cluster(genome: str, cid: str, tags: list[str], names: list[str],
             confidence: str = "alta") -> dict:
    return {"cluster_id": cid, "genome": genome, "n_genes": len(tags),
            "genes": ",".join(tags), "gene_names": ",".join(names),
            "mean_score": "0.9900", "min_score": "0.9500",
            "confidence": confidence, "evidence": "operón;adyacente"}


def _make_data(tmp_path: Path, genomes: list[str]) -> Path:
    data = tmp_path / "data"
    for genome in genomes:
        _write_genes(data / "regulation" / genome / "genes.tsv", [
            _gene(genome, "m1", "atpA", "c1", 100, 300),
            _gene(genome, "m2", "atpB", "c1", 350, 600),
        ])
        _write_clusters(data / "regulation" / genome / "operons.tsv", [])
    # los operones se escriben con el formato mínimo que lee complexes_lib
    for genome in genomes:
        op = data / "regulation" / genome / "operons.tsv"
        with open(op, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, delimiter="\t",
                               fieldnames=["operon_id", "genes"])
            w.writeheader()
            w.writerow({"operon_id": f"{genome}:op1", "genes": "m1,m2"})
    model = data / "complexes" / M.MODEL_NAME
    model.parent.mkdir(parents=True, exist_ok=True)
    model.write_text("{}", encoding="utf-8")
    return data


def test_discover_genomes_excluye_k12_y_prefiere_operon_finder(tmp_path):
    data = tmp_path / "data"
    for genome in ("g1", "g2", "ecoli_k12"):
        _write_genes(data / "regulation" / genome / "genes.tsv",
                     [_gene(genome, "m1", "x", "c", 1, 100)])
        (data / "regulation" / genome / "operons.tsv").write_text(
            "operon_id\tgenes\nop\tm1\n", encoding="utf-8")
    # g2 además tiene el fichero histórico de operon_finder
    (data / "operon_finder").mkdir(parents=True)
    (data / "operon_finder" / "g2_operons.tsv").write_text(
        "operon_id\tgenes\nop\tm1\n", encoding="utf-8")
    # g3 tiene genes pero no operones: no debe entrar
    _write_genes(data / "regulation" / "g3" / "genes.tsv",
                 [_gene("g3", "m1", "x", "c", 1, 100)])

    jobs = M.discover_genomes(data)
    assert [j.genome for j in jobs] == ["g1", "g2"]
    assert jobs[0].operon_source == "regulation"
    assert jobs[1].operon_source == "operon_finder"
    assert jobs[1].operons == data / "operon_finder" / "g2_operons.tsv"


def test_convert_clusters_coordenadas_y_filtro(tmp_path):
    genes = tmp_path / "genes.tsv"
    _write_genes(genes, [
        _gene("g", "m1", "atpA", "c1", 100, 200),
        _gene("g", "m2", "atpB", "c1", 300, 400),
        _gene("g", "m3", "atpC", "c2", 500, 600),
        # sin coordenadas válidas: no debe contar
        {"genome": "g", "contig": "c1", "start": 0, "end": 50,
         "strand": "1", "gene_id": "m4", "locus_tag": "", "ec": "",
         "gene_name": "atpD", "product": ""},
    ])
    coords = M.load_gene_coords(genes)
    clusters = [
        _cluster("g", "g:cx1", ["m1", "m2", "m4"],
                 ["atpA", "atpB", "atpD"]),
        _cluster("g", "g:cx2", ["m1", "m4"], ["atpA", "atpD"]),
        _cluster("g", "g:cx3", ["m1", "m3"], ["atpA", "atpC"]),
        _cluster("g", "g:cx4", ["m2", "m3", "m1"],
                 ["atpB", "atpC", "atpA"]),
    ]
    rows, omitted = M.convert_clusters(clusters, coords)
    assert omitted == 2
    assert [r["complex_id"] for r in rows] == ["g:cx1", "g:cx4"]
    r = rows[0]
    assert r["contig"] == "c1"
    assert (r["start"], r["end"]) == (100, 400)
    assert r["subunits"] == "atpA:m1;atpB:m2"
    assert r["expected"] == r["missing"] == ""
    assert r["status"] == "complex_candidate"
    assert r["flag"] == "yellow"
    assert r["source"] == "modelo K-12 (build_candidates.py)"
    assert r["confidence"] == "alta"
    assert "modelo K-12" in r["evidence"]
    assert "operón" in r["evidence"]
    assert "no cierra brechas" in r["limitation"].lower()
    # cx4 conserva el contig con más genes localizados
    assert rows[1]["contig"] == "c1"
    assert rows[1]["subunits"] == "atpB:m2;atpA:m1"


def _fake_runner(calls: list[list[str]]):
    def run(cmd: list[str]):
        calls.append(list(cmd))
        genome = cmd[cmd.index("--genome") + 1]
        tmp = Path(cmd[cmd.index("--out-dir") + 1])
        _write_clusters(tmp / f"clusters_{genome}.tsv", [
            _cluster(genome, f"{genome}:cx1", ["m1", "m2"],
                     ["atpA", "atpB"]),
            _cluster(genome, f"{genome}:cx2", ["m1"], ["atpA"]),
        ])
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")
    return run


def test_main_limit_force_y_resumen(tmp_path, monkeypatch):
    data = _make_data(tmp_path, ["g1", "g2"])
    out = tmp_path / "out"
    calls: list[list[str]] = []
    monkeypatch.setattr(M, "_run_build", _fake_runner(calls))

    argv = ["--data-root", str(data), "--out-dir", str(out),
            "--limit", "1", "--threads", "1"]
    assert M.main(argv) == 0
    assert len(calls) == 1
    assert calls[0][calls[0].index("--string") + 1] == "/dev/null"
    assert calls[0][calls[0].index("--synteny") + 1] == "/dev/null"
    assert calls[0][calls[0].index("--model") + 1].endswith(
        M.MODEL_NAME)
    assert "--apply-only" in calls[0]
    assert calls[0][calls[0].index("--genome") + 1] == "g1"

    path = out / "model_g1.complexes.tsv"
    assert path.is_file()
    with open(path, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert list(rows[0]) == M.COLS
    assert len(rows) == 1  # cx2 tiene un solo gen y se omite

    summary_path = tmp_path / "reports" / "complex_model_400_summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    assert summary["genomas"] == 1
    assert summary["clusters"] == 2
    assert summary["altas"] == 1
    assert summary["fallos"] == 0
    assert summary["omitidos"] == 1
    tsv = tmp_path / "reports" / "complex_model_400_summary.tsv"
    with open(tsv, newline="", encoding="utf-8") as fh:
        assert list(csv.DictReader(fh, delimiter="\t"))[0]["genome"] == "g1"

    # sin --force no se vuelve a llamar y la salida existente se cuenta
    assert M.main(argv) == 0
    assert len(calls) == 1
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    assert summary["saltados"] == 1
    assert summary["altas"] == 1

    # con --force sí se rehace
    assert M.main(argv + ["--force"]) == 0
    assert len(calls) == 2


def test_main_registra_fallo(tmp_path, monkeypatch):
    data = _make_data(tmp_path, ["g1"])
    out = tmp_path / "out"

    def boom(cmd):
        raise subprocess.TimeoutExpired(cmd, 1)

    monkeypatch.setattr(M, "_run_build", boom)
    argv = ["--data-root", str(data), "--out-dir", str(out),
            "--limit", "1", "--threads", "1"]
    assert M.main(argv) == 1
    summary = json.loads(
        (tmp_path / "reports" / "complex_model_400_summary.json").read_text(
            encoding="utf-8"))
    assert summary["fallos"] == 1
    assert summary["resultados"][0]["status"] == "fallo"
    assert "TimeoutExpired" in summary["resultados"][0]["error"]
    assert not (out / "model_g1.complexes.tsv").exists()


def test_integracion_real_con_build_candidates(tmp_path):
    """Sin monkeypatch: el subproceso real corre sin STRING ni sintenia."""
    data = tmp_path / "data"
    genome = "g1"
    _write_genes(data / "regulation" / genome / "genes.tsv", [
        _gene(genome, "m1", "atpA", "c1", 100, 300),
        _gene(genome, "m2", "atpB", "c1", 350, 600),
        _gene(genome, "m3", "atpC", "c1", 650, 900),
    ])
    op = data / "regulation" / genome / "operons.tsv"
    with open(op, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=["operon_id", "genes"])
        w.writeheader()
        w.writerow({"operon_id": "op1", "genes": "m1,m2,m3"})
    model = data / "complexes" / M.MODEL_NAME
    model.parent.mkdir(parents=True, exist_ok=True)
    model.write_text(json.dumps({
        "version": 1,
        "features": ["same_operon", "adjacent", "log_gap", "log_dist",
                     "synteny", "same_ec"],
        "mean": [0, 0, 0, 0, 0, 0],
        "scale": [1, 1, 1, 1, 1, 1],
        "coef": [3, 2, 0, 0, 0, 0],
        "intercept": 1.0,
    }), encoding="utf-8")

    out = tmp_path / "out"
    assert M.main(["--data-root", str(data), "--out-dir", str(out),
                   "--limit", "1", "--threads", "1"]) == 0
    path = out / "model_g1.complexes.tsv"
    with open(path, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert rows, "el modelo sintético debería producir al menos un cluster"
    assert rows[0]["status"] == "complex_candidate"
    assert rows[0]["flag"] == "yellow"
    assert rows[0]["contig"] == "c1"
    assert int(rows[0]["start"]) == 100
    assert int(rows[0]["end"]) == 900
    assert "modelo K-12" in rows[0]["evidence"]
