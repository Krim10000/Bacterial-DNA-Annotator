#!/usr/bin/env python3
"""Corre la anotacion de membrana en los 10 genomas piloto y escribe summary.json."""
from __future__ import annotations

import argparse
import csv
import json
import platform
import re
import shutil
import subprocess
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
TYPES = ["tm", "signal_peptide", "lipoprotein", "transporter", "secretion"]


def tool_version(cmd: str) -> str:
    try:
        out = subprocess.run([cmd, "-h"], capture_output=True, text=True, timeout=20).stdout
        for line in out.splitlines():
            if "HMMER" in line:
                return line.lstrip("# ").strip()
        for line in out.splitlines():
            if line.startswith("#"):
                return line.lstrip("# ").strip()
        return out.splitlines()[0].strip() if out else "desconocida"
    except Exception:
        return "no disponible"


def read_pilot(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def genome_stats(tsv: Path) -> dict:
    counts = Counter()
    conf = Counter()
    sources = Counter()
    substrates = Counter()
    families = Counter()
    sec_sub = Counter()
    substrate_rows = 0
    rows = 0
    with open(tsv, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            rows += 1
            counts[r["type"]] += 1
            conf[r["confidence"]] += 1
            ev = r["evidence"]
            if "hmmscan:" in ev:
                sources["hmmscan_pfam"] += 1
            if "kd_tm=" in ev:
                sources["kd_tm"] += 1
            if "sp_heur=" in ev:
                sources["sp_heur"] += 1
            if "lipobox=" in ev:
                sources["lipobox"] += 1
            if "swissprot=" in ev:
                sources["swissprot"] += 1
            if "swissprot_lipoprotein_text=" in ev:
                sources["swissprot_lipoprotein_text"] += 1
            if r["substrate"].strip():
                substrate_rows += 1
                for term in r["substrate"].split("|"):
                    substrates[term] += 1
            families[r["family"].split(" (")[0]] += 1
            if r["type"] == "secretion":
                m = re.search(r"\(([^;)]+)", r["family"])
                if m:
                    sec_sub[m.group(1)] += 1
    return dict(rows=rows, counts={t: counts.get(t, 0) for t in TYPES},
                confidence={c: conf.get(c, 0) for c in ("high", "medium", "low")},
                substrate_rows=substrate_rows,
                substrates=dict(substrates.most_common(12)),
                top_families=dict(families.most_common(12)),
                secretion_subgroups=dict(sec_sub.most_common()),
                evidence_sources=dict(sources))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pilot", default=str(ROOT / "data/membrane/pilot_genomes.tsv"))
    ap.add_argument("--outdir", default=str(ROOT / "data/membrane"))
    ap.add_argument("--annot-root", default=str(ROOT / "data/annotations_ab_initio"))
    ap.add_argument("--panel-dir", default=str(ROOT / "data/membrane/pfam_panel"))
    ap.add_argument("--cpu", type=int, default=4)
    ap.add_argument("--skip-existing", action="store_true")
    ap.add_argument("--summary", default=str(ROOT / "data/membrane/summary.json"))
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    pilot = read_pilot(Path(args.pilot))
    if len(pilot) != 10:
        sys.exit(f"[run] se esperaban 10 genomas piloto, hay {len(pilot)}")

    genomes = []
    failed = []
    for i, meta in enumerate(pilot, 1):
        g = meta["genome"]
        tsv = outdir / f"{g}.tsv"
        gff = outdir / "gff" / f"{g}.gff3"
        print(f"[run] {i}/10 {g} ({meta['organism']})", flush=True)
        if args.skip_existing and tsv.exists() and gff.exists():
            print("  ya existe; se reutiliza", flush=True)
        else:
            cmd = [sys.executable, str(HERE / "03_annotate_membrane.py"),
                   "--genome", g, "--annot-root", args.annot_root,
                   "--panel-dir", args.panel_dir, "--outdir", str(outdir),
                   "--cpu", str(args.cpu)]
            res = subprocess.run(cmd, capture_output=True, text=True)
            if res.returncode != 0:
                print(res.stderr[-3000:], file=sys.stderr)
                failed.append(g)
                continue
            print("  " + res.stdout.strip(), flush=True)
        st = genome_stats(tsv)
        genomes.append(dict(
            genome=g, phylum=meta.get("phylum", ""), organism=meta.get("organism", ""),
            rationale=meta.get("rationale", ""),
            contigs=int(meta.get("contigs", 0)), bp=int(meta.get("bp", 0)),
            cds=int(meta.get("cds", 0)), proteins=int(meta.get("proteins", 0)),
            tsv=str(tsv.relative_to(ROOT)), gff=str(gff.relative_to(ROOT)),
            **st,
        ))

    by_type = Counter()
    by_conf = Counter()
    for g in genomes:
        by_type.update(g["counts"])
        by_conf.update(g["confidence"])
    n_transporters = by_type.get("transporter", 0)
    substrate_rows = sum(g["substrate_rows"] for g in genomes)
    all_substrates = Counter()
    all_sec_sub = Counter()
    for g in genomes:
        all_substrates.update(g.get("substrates", {}))
        all_sec_sub.update(g.get("secretion_subgroups", {}))

    # validacion
    val_path = outdir / "validation.json"
    val_cmd = [sys.executable, str(HERE / "05_validate_outputs.py"),
               "--pilot", args.pilot, "--outdir", str(outdir), "--json", str(val_path)]
    val = subprocess.run(val_cmd, capture_output=True, text=True)
    validation = {"ok": False, "errors_total": None, "path": str(val_path.relative_to(ROOT))}
    if val_path.exists():
        v = json.load(open(val_path))
        validation = dict(ok=v.get("ok", False), errors_total=v.get("totals", {}).get("errors"),
                          genomes_ok=v.get("totals", {}).get("genomes_ok"),
                          path=str(val_path.relative_to(ROOT)))

    panel_models_path = Path(args.panel_dir) / "panel_models.tsv"
    n_panel = sum(1 for _ in csv.DictReader(open(panel_models_path), delimiter="\t"))

    pfam_first_date = ""
    try:
        with open(ROOT / ".cache/pfam/Pfam-A.hmm", errors="replace") as fh:
            for line in fh:
                if line.startswith("DATE "):
                    pfam_first_date = line[5:].strip()
                    break
    except OSError:
        pass

    summary = {
        "schema_version": 1,
        "module": "membrane",
        "agent": "membrane-builder",
        "generated_at": datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z"),
        "scope": "10 genomas piloto; capa de membrana/transporte/secrecion",
        "inputs": {
            "annotations": str((ROOT / "data/annotations_ab_initio").relative_to(ROOT)),
            "gff": "<genome>/annotation.gff3",
            "proteins": "<genome>/proteins.faa",
            "swissprot": "<genome>/functional.tsv (DIAMOND, solo lectura)",
            "pfam": ".cache/pfam/Pfam-A.hmm (solo lectura)",
            "pfam_first_model_date": pfam_first_date,
        },
        "method": {
            "pfam_panel": "%d modelos Pfam curados" % n_panel,
            "hmmscan": "hmmscan --cut_ga --noali --domtblout; 1 hit principal por proteina",
            "tm": "Kyte-Doolittle, ventana 19 aa, media >= 2.0, segmento >= 18 aa; >=2 segmentos, o SP + 1 segmento; NO es TMHMM",
            "signal_peptide": "h-region 9-13 aa en primeros 21 aa, media KD >= 2.2, residuos -1/-3 pequenos; NO es SignalP",
            "lipobox": "[LVI][ASTM][AG]C en primeros 40 aa + ventana de 7 aa en los 12 aa previos con media KD >= 1.8; NO es LipoP",
            "substrate": "solo texto Swiss-Prot explicito en transporte/SBP con contexto de transporte (pident >= 50, E <= 1e-5); max 2 terminos; si no, vacio",
            "classification_priority": ["secretion",
                                        "lipoprotein (lipobox + Pfam lipoprotein)",
                                        "transporter",
                                        "lipoprotein (lipobox solo)",
                                        "tm (>=2 segmentos KD)",
                                        "tm (SP + 1 TM fuera de la h-region)",
                                        "signal_peptide (SP sin TM)",
                                        "tm (Pfam membrane)"],
            "one_row_per_protein": True,
        },
        "tools": {
            "hmmer_hmmscan": tool_version("hmmscan"),
            "hmmer_hmmpress": tool_version("hmmpress"),
            "signalp": "no instalado" if shutil.which("signalp") is None else shutil.which("signalp"),
            "tmhmm": "no instalado" if shutil.which("tmhmm") is None else shutil.which("tmhmm"),
            "phobius": "no instalado" if shutil.which("phobius") is None else shutil.which("phobius"),
            "psortb": "no instalado" if shutil.which("psortb") is None else shutil.which("psortb"),
            "macsyfinder": "no instalado" if shutil.which("macsyfinder") is None else shutil.which("macsyfinder"),
            "python": platform.python_version(),
        },
        "pilot": [{k: m.get(k, "") for k in ("genome", "phylum", "organism", "rationale")} for m in pilot],
        "genomes": genomes,
        "totals": {
            "genomes": len(genomes),
            "proteins": sum(g["proteins"] for g in genomes),
            "rows": sum(g["rows"] for g in genomes),
            "by_type": {t: by_type.get(t, 0) for t in TYPES},
            "by_confidence": {c: by_conf.get(c, 0) for c in ("high", "medium", "low")},
            "substrate_rows": substrate_rows,
            "substrates_top": dict(all_substrates.most_common(15)),
            "secretion_subgroups": dict(all_sec_sub.most_common()),
            "transporter_with_substrate": None,
            "substrate_fraction_transporters": round(substrate_rows / n_transporters, 3) if n_transporters else 0.0,
        },
        "validation": validation,
        "limitations": [
            "SignalP, TMHMM, Phobius, PSORTb y MacSyFinder no estan instalados: TM/SP/lipobox son heuristicas, no predicciones de referencia.",
            "El panel Pfam es curado (255 modelos) y no cubre todo el espacio de transporte; familias fuera del panel no se anotan.",
            "No hay TCDB completo ni asignacion de sustrato por familia: el sustrato solo se toma de texto Swiss-Prot explicito.",
            "hmmscan usa --cut_ga de Pfam; proteinas divergentes o truncadas pueden quedar sin hit.",
            "tm puede tener falsos positivos en proteinas solubles con parches hidrofobicos; se marca low salvo >=4 segmentos.",
            "Una fila por proteina: una proteina con dos roles (p. ej. lipoproteina de un sistema de secrecion) recibe un solo type y el resto queda en family/evidence.",
            "No se distingue membrana interna/externa ni topologia; no hay prediccion de localizacion experimental.",
            "T1SS esta subrepresentado: el panel no incluye HlyD/TolC y solo usa ABC_export como aproximacion.",
            "secretion incluye componentes de sistemas y sustratos (T9SS, sortasa); el subgrupo en family/evidence distingue ambos.",
            "Un SP seguido de una TM fuera de la h-region se reporta como tm (ancla de membrana), no como signal_peptide.",
            "No se valido contra una referencia curada de membrana (gold standard); la validacion es estructural (coordenadas, GFF, campos).",
        ],
        "reproduce": [
            "python3 scripts/membrane/01_select_genomes.py",
            "python3 scripts/membrane/02_build_pfam_panel.py",
            "python3 scripts/membrane/04_run_pilot.py --cpu 4",
            "python3 scripts/membrane/05_validate_outputs.py",
        ],
    }
    # corregir el campo derivado de sustrato en transportadores
    n_tr_sub = 0
    for g in genomes:
        tsv = ROOT / g["tsv"]
        with open(tsv, newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r["type"] == "transporter" and r["substrate"].strip():
                    n_tr_sub += 1
    summary["totals"]["transporter_with_substrate"] = n_tr_sub
    summary["totals"]["substrate_fraction_transporters"] = round(n_tr_sub / n_transporters, 3) if n_transporters else 0.0

    Path(args.summary).parent.mkdir(parents=True, exist_ok=True)
    with open(args.summary, "w") as fh:
        json.dump(summary, fh, indent=2, ensure_ascii=False)
        fh.write("\n")
    print(f"[run] summary -> {args.summary}")
    print(f"[run] filas: {summary['totals']['rows']}; por tipo: {summary['totals']['by_type']}")
    print(f"[run] validacion: {validation}")
    if failed:
        print(f"[run] genomas fallidos: {failed}", file=sys.stderr)
        return 1
    return 0 if validation.get("ok") else 1


if __name__ == "__main__":
    raise SystemExit(main())
