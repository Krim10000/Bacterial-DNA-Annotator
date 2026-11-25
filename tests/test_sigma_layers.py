"""Pruebas del evaluador por capas y del cargador de gold."""
from __future__ import annotations

import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "sigma_eval_layers", ROOT / "scripts" / "sigma" / "50_eval_layers.py")
MOD = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MOD)


def test_selftest():
    assert "OK" in MOD.selftest()


def test_map_sigma():
    assert MOD.map_sigma("sigma70", "k12") == ("sigma70", "ok")
    assert MOD.map_sigma("A", "subtilis") == ("sigma70", "ok")
    assert MOD.map_sigma("E", "subtilis") == ("", "no_pwm")
    assert MOD.map_sigma("?", "k12") == ("", "unlabeled")


def test_gold_groups_multi_sigma(tmp_path: Path):
    path = tmp_path / "tss_sigma.tsv"
    path.write_text(
        "species\tgenome_id\tcontig\tpos\tstrand\tsigma\tgene\tevidence\t"
        "confidence\n"
        "E. coli\tecoli_k12\tNC_1\t100\t+\tsigma70\tthrL\ttest\tS\n"
        "E. coli\tecoli_k12\tNC_1\t100\t+\tsigma32\tthrL\ttest\tS\n",
        encoding="utf-8")
    rows = MOD.load_gold(path, "gold")
    assert len(rows) == 1
    rec = rows[0]
    assert rec["tss"] == 100
    assert rec["strand"] == "+"
    assert rec["sigmas"] == ["sigma32", "sigma70"]
    assert rec["n_rows"] == 2


def test_gold_accepts_pos_and_letters(tmp_path: Path):
    path = tmp_path / "tss_sigma.tsv"
    path.write_text(
        "contig\tpos\tstrand\tsigma\nbasu168\t157\t+\tA\n",
        encoding="utf-8")
    rows = MOD.load_gold(path, "gold")
    assert rows[0]["sigmas"] == ["sigma70"]
    assert rows[0]["strand"] == "+"
