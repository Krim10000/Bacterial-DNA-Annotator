#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E2E de 09_make_gene_map.py: los 4 modos de anotación (prokka GFF+FAA,
bakta TSV, eggNOG+GFF) + filtro por red + integración con 07.load_gene_map
y build_payload del portal."""
import csv
import importlib
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
M05 = ROOT / "scripts" / "metabolism"
ANN = ROOT / "scripts" / "testdata" / "annotation"
NET = ROOT / "scripts" / "testdata" / "net"
OUT = ANN / "out"

sys.path.insert(0, str(M05))
m07 = importlib.import_module("07_network_view")

shutil.rmtree(OUT, ignore_errors=True)
OUT.mkdir(parents=True, exist_ok=True)


def run(mode_args, name):
    cmd = [sys.executable, "-I", str(M05 / "09_make_gene_map.py"),
           "--genome", "GCF_TEST", "--rhea-dir", str(ANN),
           "--out-dir", str(OUT)] + mode_args
    r = subprocess.run(cmd, capture_output=True, text=True)
    assert r.returncode == 0, (name, r.stdout[-800:], r.stderr[-800:])
    p = OUT / "GCF_TEST_gene_map.tsv"
    with open(p, encoding="utf-8") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    return rows, r.stdout + r.stderr   # el log va a stderr


# --- modo prokka (GFF + FAA): EC del product, sequence, strand, hipotéticas
rows, out = run(["--gff", str(ANN / "GCF_TEST_prokka.gff"),
                 "--faa", str(ANN / "GCF_TEST_prokka.faa")], "prokka")
assert "sequence" in rows[0], "con FAA debe incluir columna sequence"
r1 = [r for r in rows if r["gene_id"] == "GTEST_00001"]
assert r1 and r1[0]["rhea_id"] == "RHEA:10001" and r1[0]["strand"] == "+" \
    and r1[0]["gene_length_bp"] == "901" and r1[0]["sequence"], r1
assert any(r["rhea_id"] == "RHEA:20002" and r["gene_id"] == "GLK1" for r in rows), \
    "GLK1 debe mapear a las variantes UN/BI de su EC"
hyp = [r for r in rows if r["gene_id"] == "PROK2"]
assert not hyp, "la hipotética sin EC no debe producir filas"
assert all(r["hypothetical"] == "false" for r in rows), rows
print(f"prokka: {len(rows)} filas, genes {sorted({r['gene_id'] for r in rows})}")

# --- modo bakta (EC en columna dedicada, sin FAA -> sin sequence)
rows_b, _ = run(["--bakta", str(ANN / "GCF_TEST_bakta.tsv")], "bakta")
assert "sequence" not in rows_b[0], "sin FAA no debe haber columna sequence"
assert any(r["gene_id"] == "WPROT1" and r["rhea_id"] == "RHEA:10001"
           for r in rows_b), rows_b
assert not any(r["gene_id"] == "WPROT2" for r in rows_b), \
    "hipotética bakta sin EC: sin filas"
print(f"bakta: {len(rows_b)} filas")

# --- modo eggNOG + GFF: eggNOG ENRIQUECE el gene del GFF (sin duplicar filas)
rows_e, _ = run(["--eggnog", str(ANN / "GCF_TEST_eggnog.emapper.annotations"),
                 "--gff", str(ANN / "GCF_TEST_prokka.gff"),
                 "--faa", str(ANN / "GCF_TEST_prokka.faa")], "eggnog")
glk = [r for r in rows_e if r["gene_id"] == "GLK1"]
assert glk and glk[0]["contig"] == "ctg1" and glk[0]["start"] == "4000" \
    and glk[0]["end"] == "5300" and glk[0]["strand"] == "+" \
    and glk[0]["sequence"] == "MEIKVMEIKV", glk
assert sorted(r["rhea_id"] for r in glk) == ["RHEA:20001", "RHEA:20002"], \
    "GLK1 debe aparecer una vez por variante Rhea de su EC"
assert len(rows_e) == len(rows), \
    "eggnog+gff no debe añadir genes: solo enriquece los del GFF"
print(f"eggnog: {len(rows_e)} filas, sin duplicados (GLK1 con coords del GFF)")

# --- filtro por la red: solo rhea presentes en GCF_TEST_network.json
import json
net = json.loads((NET / "GCF_TEST_network.json").read_text(encoding="utf-8"))
net_ids = {n["id"] for n in net.get("nodes", []) if n.get("type") == "reaction"}
rows_f, _ = run(["--gff", str(ANN / "GCF_TEST_prokka.gff"),
                 "--faa", str(ANN / "GCF_TEST_prokka.faa"),
                 "--net-dir", str(NET)], "filtro-red")
assert rows_f and all(f"RHEA:{r['rhea_id'].split(':')[1]}" in net_ids
                      or r["rhea_id"] in net_ids for r in rows_f), \
    (rows_f[:3], list(net_ids)[:5])
extra = [r for r in rows_f if r["rhea_id"] not in net_ids]
assert not extra, extra
print(f"filtro-red: {len(rows_f)}/{len(rows)} filas dentro de la red")

# --- integración con el portal: load_gene_map + build_payload
import importlib as _il
m04 = _il.import_module("04_network_trace")
gm = m07.load_gene_map(OUT / "GCF_TEST_gene_map.tsv")
assert gm and all(isinstance(k, int) for k in gm), gm
assert gm[10001][0]["gene_id"] == "GTEST_00001"
payload = m07.build_payload("GCF_TEST", NET, m04.RHEA_DIR, gm, None,
                            set(), False)   # rhea_dir real: cache TSV de Rhea
with_genes = [n for n in payload["nodes"] if n.get("genes")]
assert with_genes, "el payload debe incluir genes con el gene_map generado"
assert any(g["gene_id"] == "GTEST_00001" for n in with_genes
           for g in n["genes"]), payload["stats"]
print(f"portal: {len(with_genes)} reacciones con procedencia génica "
      f"({payload['stats']['hypothetical_rxns']} hipotéticas)")

# --- auto-descubrimiento con .gz (sin --gff/--faa): layout de data/annotations
import gzip as _gz
AUTO = ANN / "auto"
shutil.rmtree(AUTO, ignore_errors=True)
AUTO.mkdir()
for src in ("GCF_TEST_prokka.gff", "GCF_TEST_prokka.faa"):
    (AUTO / f"{src}.gz").write_bytes(
        _gz.compress((ANN / src).read_bytes()))
rows_a, out_a = run(["--annotations-dir", str(AUTO)], "auto-gz")
assert "auto-descubierta" in out_a and "GCF_TEST_prokka.gff.gz" in out_a, out_a
assert rows_a == rows, \
    "el modo auto con .gz debe dar las mismas filas que --gff/--faa"
print(f"auto-gz: {len(rows_a)} filas idénticas al modo explícito "
      f"(gff+faa descubiertos en {AUTO.name}/)")

# --- error claro si no hay anotación para el genoma
r_err = subprocess.run(
    [sys.executable, "-I", str(M05 / "09_make_gene_map.py"),
     "--genome", "GENOME_INEXISTENTE", "--rhea-dir", str(ANN),
     "--annotations-dir", str(AUTO), "--out-dir", str(OUT)],
    capture_output=True, text=True)
assert r_err.returncode != 0 and "sin anotación para" in r_err.stderr, \
    (r_err.returncode, r_err.stderr[-300:])
print("auto-gz: genoma inexistente -> error con mensaje útil")

# --- 08: lectores FASTA gz-aware y resolución {genome}.faa[.gz]
m08 = _il.import_module("08_fold_candidates")
p_faa = m08.find_seq_file(AUTO, "GCF_TEST", ".faa")
assert p_faa and p_faa.name == "GCF_TEST_prokka.faa.gz", p_faa
prot = m08.read_fasta_dict(p_faa)
assert set(prot) == {"GTEST_00001", "GLK1"}, prot
assert m08.read_fasta_first(p_faa).startswith("MEIKV")
assert m08.find_seq_file(AUTO, "OTRO_GENOMA", ".faa") is None
print("08: find_seq_file + read_fasta_dict leen .faa.gz sin descomprimir")

print("TEST 09_make_gene_map: OK — prokka/bakta/eggnog + auto-descubrimiento "
      ".gz + filtro por red + integración 07 + FASTA gz en 08")
