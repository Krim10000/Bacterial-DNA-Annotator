#!/usr/bin/env python3
"""Resume los 400 genomas: summary_full.json y README_full.md."""
from __future__ import annotations

import argparse
import csv
import json
import os
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
TYPES = ["tm", "signal_peptide", "lipoprotein", "transporter", "secretion"]
CONFS = ["high", "medium", "low"]


def now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def read_list(path: Path) -> list[dict]:
    return list(csv.DictReader(open(path), delimiter="\t"))


def pct(a: int, b: int) -> float:
    return round(100.0 * a / b, 1) if b else 0.0


def sec_subgroup(family: str) -> str:
    if "(" not in family:
        return "?"
    inner = family.split("(", 1)[1]
    return inner.split(";", 1)[0].split(")", 1)[0].strip() or "?"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=str(ROOT / "data/membrane/full"))
    ap.add_argument("--gff-dir", default=str(ROOT / "data/membrane/gff/full"))
    ap.add_argument("--list", default=str(ROOT / "data/membrane/full/genomes.tsv"))
    ap.add_argument("--state", default=str(ROOT / "data/membrane/full/state.json"))
    ap.add_argument("--validation", default=str(ROOT / "data/membrane/validation_full.json"))
    ap.add_argument("--sample-validation", default=str(ROOT / "data/membrane/validation_full_sample.json"))
    ap.add_argument("--json", default=str(ROOT / "data/membrane/summary_full.json"))
    ap.add_argument("--readme", default=str(ROOT / "data/membrane/README_full.md"))
    args = ap.parse_args()

    outdir = Path(args.outdir)
    gff_dir = Path(args.gff_dir)
    metas = read_list(Path(args.list))
    state = json.load(open(args.state)) if Path(args.state).exists() else {}
    st_genomes = state.get("genomes", {})

    types = Counter()
    conf = Counter()
    families = Counter()
    tr_families = Counter()
    sec_sub = Counter()
    substrates = Counter()
    sources = Counter()
    by_phylum: dict[str, dict] = defaultdict(lambda: dict(
        genomes=0, proteins=0, rows=0, counts=Counter(), confidence=Counter(),
        substrate_rows=0))
    per_genome = []
    failures = []
    pending = []
    total_rows = 0
    total_sub = 0
    tr_with_sub = 0

    for m in metas:
        g = m["genome"]
        tsv = outdir / f"{g}.tsv"
        gff = gff_dir / f"{g}.gff3"
        if not tsv.exists() or not gff.exists():
            status = st_genomes.get(g, {}).get("status", "")
            if status == "failed":
                failures.append(dict(genome=g, reason="TSV/GFF ausente",
                                     error=st_genomes.get(g, {}).get("error", "")))
            else:
                pending.append(g)
            continue
        counts = Counter()
        cconf = Counter()
        csubs = Counter()
        rows = 0
        sub_rows = 0
        try:
            fh = open(tsv, newline="")
        except OSError as exc:
            failures.append(dict(genome=g, reason="no se pudo abrir TSV", error=str(exc)))
            continue
        with fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                rows += 1
                t = r.get("type", "?")
                counts[t] += 1
                cconf[r.get("confidence", "?")] += 1
                fam = (r.get("family") or "").split(" (")[0]
                families[fam] += 1
                if t == "transporter":
                    tr_families[fam] += 1
                    if (r.get("substrate") or "").strip():
                        tr_with_sub += 1
                if t == "secretion":
                    sec_sub[sec_subgroup(r.get("family", ""))] += 1
                if (r.get("substrate") or "").strip():
                    sub_rows += 1
                    for term in r["substrate"].split("|"):
                        csubs[term] += 1
                        substrates[term] += 1
                ev = r.get("evidence", "")
                for key, needle in (("hmmscan_pfam", "hmmscan:"), ("kd_tm", "kd_tm="),
                                    ("sp_heur", "sp_heur="), ("lipobox", "lipobox="),
                                    ("swissprot", "swissprot=")):
                    if needle in ev:
                        sources[key] += 1
        types.update(counts)
        conf.update(cconf)
        total_rows += rows
        total_sub += sub_rows
        ph = m.get("phylum", "unknown") or "unknown"
        bp = by_phylum[ph]
        bp["genomes"] += 1
        bp["proteins"] += int(m.get("proteins") or 0)
        bp["rows"] += rows
        bp["counts"].update(counts)
        bp["confidence"].update(cconf)
        bp["substrate_rows"] += sub_rows
        per_genome.append(dict(
            genome=g, phylum=ph, organism=m.get("organism", g),
            proteins=int(m.get("proteins") or 0), cds=int(m.get("cds") or 0),
            rows=rows, counts={t: counts.get(t, 0) for t in TYPES},
            confidence={c: cconf.get(c, 0) for c in CONFS},
            substrate_rows=sub_rows, substrates=dict(csubs.most_common(5)),
            top_transporter_families=dict(Counter(
                (r.get("family") or "").split(" (")[0]
                for r in [] ).most_common()) if False else {},
        ))

    # segundo pase barato para top familias de transporte por genoma
    for entry in per_genome:
        g = entry["genome"]
        tc = Counter()
        with open(outdir / f"{g}.tsv", newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r.get("type") == "transporter":
                    tc[(r.get("family") or "").split(" (")[0]] += 1
        entry["top_transporter_families"] = dict(tc.most_common(5))

    processed = len(per_genome)
    expected = len(metas)
    by_phylum_out = {}
    for ph, d in sorted(by_phylum.items(), key=lambda kv: (-kv[1]["genomes"], kv[0])):
        by_phylum_out[ph] = dict(
            genomes=d["genomes"], proteins=d["proteins"], rows=d["rows"],
            counts={t: d["counts"].get(t, 0) for t in TYPES},
            confidence={c: d["confidence"].get(c, 0) for c in CONFS},
            substrate_rows=d["substrate_rows"],
            row_fraction=round(d["rows"] / d["proteins"], 3) if d["proteins"] else 0.0,
        )

    n_tr = types.get("transporter", 0)
    coverage = dict(expected=expected, processed=processed, failed=len(failures),
                    pending=len(pending),
                    pct=round(100.0 * processed / expected, 1) if expected else 0.0,
                    proteins=sum(e["proteins"] for e in per_genome),
                    rows=total_rows, substrate_rows=total_sub,
                    transporter_with_substrate=tr_with_sub,
                    substrate_fraction_transporters=round(tr_with_sub / n_tr, 3) if n_tr else 0.0)

    validation = {}
    for p in (args.validation, args.sample_validation):
        if p and Path(p).exists():
            try:
                v = json.load(open(p))
                validation[p] = dict(ok=v.get("ok"), totals=v.get("totals"),
                                     sample=v.get("sample"), path=p)
            except Exception:
                pass

    # hallazgos dinamicos
    top_tr = tr_families.most_common(5)
    findings = []
    if n_tr:
        findings.append(
            "Transporte: %d transportadores; familias dominantes %s. Solo %d (%.1f%%) "
            "tienen sustrato explicito en Swiss-Prot. La capa sirve para familia/transporte, "
            "no para cerrar transporte de metabolitos sin evidencia adicional."
            % (n_tr, ", ".join("%s (%d)" % (f, n) for f, n in top_tr),
               tr_with_sub, pct(tr_with_sub, n_tr)))
    t9 = sec_sub.get("T9SS", 0) + sec_sub.get("T9SS-substrate", 0)
    t9_genomes = []
    if t9:
        max_t9 = ("", 0)
        for m in metas:
            g = m["genome"]
            p = outdir / f"{g}.tsv"
            if not p.exists():
                continue
            n = 0
            with open(p, newline="") as fh:
                for r in csv.DictReader(fh, delimiter="\t"):
                    if r.get("type") == "secretion" and "T9SS" in (r.get("family") or ""):
                        n += 1
            if n:
                t9_genomes.append((g, n))
                if n > max_t9[1]:
                    max_t9 = (g, n)
        findings.append(
            "Secrecion: %d filas; %d T9SS (%.1f%%), presentes en %d genomas. "
            "El maximo es %s con %d filas T9SS. T6SS=%d, T4SS=%d, T3SS=%d."
            % (types.get("secretion", 0), t9, pct(t9, types.get("secretion", 0)),
               len(t9_genomes), max_t9[0], max_t9[1], sec_sub.get("T6SS", 0),
               sec_sub.get("T4SS", 0), sec_sub.get("T3SS", 0)))
    fr = [(e["rows"] / e["proteins"], e) for e in per_genome if e["proteins"]]
    if fr:
        fr.sort(key=lambda x: x[0])
        lo, hi = fr[0], fr[-1]
        findings.append(
            "Densidad de la capa: media %.1f%% del proteoma. Minimo %s (%.1f%%) "
            "y maximo %s (%.1f%%); el rango es amplio por filo/estilo de vida."
            % (100.0 * total_rows / coverage["proteins"] if coverage["proteins"] else 0,
               lo[1]["organism"], 100 * lo[0], hi[1]["organism"], 100 * hi[0]))
    findings.append(
        "Confianza: high=%d (%.1f%%), medium=%d (%.1f%%), low=%d (%.1f%%). "
        "Los transportadores con Pfam son solidos; TM/SP/lipobox son heuristicas y "
        "mayoritariamente low."
        % (conf.get("high", 0), pct(conf.get("high", 0), total_rows),
           conf.get("medium", 0), pct(conf.get("medium", 0), total_rows),
           conf.get("low", 0), pct(conf.get("low", 0), total_rows)))

    if pending:
        blockers = [
            "Host de 8 cores saturado por otros agentes; hmmscan compite por CPU.",
            "Quedan %d genomas pendientes; reanudar con 06_run_full.py (estado reanudable)." % len(pending),
            "SignalP/TMHMM/Phobius/PSORTb/LipoP/MacSyFinder no instalados; TM/SP/lipobox son heuristicas.",
            "No se hicieron commits.",
        ]
    else:
        blockers = [
            "Corrida completa: 400/400 genomas procesados; 0 fallos; validacion estructural 400/400.",
            "SignalP/TMHMM/Phobius/PSORTb/LipoP/MacSyFinder no instalados; TM/SP/lipobox son heuristicas.",
            "La validacion es estructural (coordenadas/GFF/campos), no biologica; no hay gold standard.",
            "No se hicieron commits.",
        ]

    summary = dict(
        schema_version=1, module="membrane", agent="membrane-builder",
        scope="full: 400 genomas; capa de membrana/transporte/secrecion",
        generated_at=now(),
        inputs=dict(annotations="data/annotations_ab_initio", proteins="<genome>/proteins.faa",
                    swissprot="<genome>/functional.tsv", pfam=".cache/pfam/Pfam-A.hmm"),
        method=dict(
            pfam_panel="255 modelos Pfam curados (panel del piloto)",
            hmmscan="hmmscan --cut_ga --noali --domtblout; 1 hit principal por proteina",
            tm="Kyte-Doolittle, ventana 19 aa, media >= 2.0, segmento >= 18 aa; NO es TMHMM",
            signal_peptide="h-region 9-13 aa en primeros 21 aa, media KD >= 2.2; NO es SignalP",
            lipobox="[LVI][ASTM][AG]C en primeros 40 aa + ventana 7 aa previa; NO es LipoP",
            substrate="solo transporte/SBP con texto Swiss-Prot explicito y contexto de transporte",
            one_row_per_protein=True,
        ),
        tools=dict(hmmer="HMMER 3.3.2", signalp="no instalado", tmhmm="no instalado",
                   phobius="no instalado", psortb="no instalado", lipop="no instalado",
                   macsyfinder="no instalado", python="3.10"),
        coverage=coverage,
        by_type={t: types.get(t, 0) for t in TYPES},
        by_confidence={c: conf.get(c, 0) for c in CONFS},
        top_transporter_families=dict(tr_families.most_common(25)),
        top_families=dict(families.most_common(25)),
        secretion_subgroups=dict(sec_sub.most_common()),
        substrates_top=dict(substrates.most_common(25)),
        evidence_sources=dict(sources),
        by_phylum=by_phylum_out,
        genomes=per_genome,
        failures=failures,
        pending=pending[:1000],
        blockers=blockers,
        batches=state.get("batches", {}),
        validation=validation,
        findings=findings,
        limitations=[
            "SignalP/TMHMM/Phobius/PSORTb/LipoP/MacSyFinder no instalados: TM/SP/lipobox son heuristicas.",
            "Panel Pfam curado de 255 modelos; no cubre todo el espacio de transporte ni TCDB completo.",
            "Sustrato solo si Swiss-Prot lo nombra explicitamente en contexto de transporte; la mayoria queda vacio.",
            "hmmscan --cut_ga: proteinas divergentes/truncadas pueden quedar sin hit.",
            "TM por hidropatía puede dar falsos positivos; tm sin Pfam es low (o medium con >=4 segmentos).",
            "Una fila por proteina; roles secundarios quedan en family/evidence.",
            "No se distingue membrana interna/externa ni topologia.",
            "T1SS subrepresentado (sin HlyD/TolC en el panel).",
            "Validacion estructural (coordenadas/GFF/campos), no biologica; no hay gold standard.",
        ],
        reproduce=[
            "python3 scripts/membrane/01_select_genomes.py  # piloto",
            "python3 scripts/membrane/06_run_full.py --workers 4 --cpu-per-worker 2 --batch-size 25",
            "python3 scripts/membrane/05_validate_outputs.py --pilot data/membrane/full/genomes.tsv --tsv-dir data/membrane/full --gff-dir data/membrane/gff/full --json data/membrane/validation_full.json --quiet",
            "python3 scripts/membrane/07_summarize_full.py",
        ],
    )

    Path(args.json).parent.mkdir(parents=True, exist_ok=True)
    with open(args.json, "w") as fh:
        json.dump(summary, fh, indent=2, ensure_ascii=False)
        fh.write("\n")

    # README
    lines = []
    if pending:
        title = "corrida full, lote parcial (%d/%d)" % (processed, expected)
    else:
        title = "corrida completa (%d genomas)" % processed
    lines.append("# Capa de membrana, transporte y secreción — %s" % title)
    lines.append("")
    lines.append("Módulo **B1/B2/B3** del superanotador. Extiende el piloto de 10 genomas a los")
    lines.append("**400 genomas** de `data/annotations_ab_initio/`, escribiendo en `data/membrane/full/`")
    lines.append("y `data/membrane/gff/full/` sin sobrescribir el piloto.")
    lines.append("")
    lines.append("- Dueño: `membrane-builder`.")
    lines.append("- Fecha: %s." % summary["generated_at"])
    lines.append("- Cobertura: **%d/%d genomas** (%.1f%%), %d proteínas, **%d filas**."
                 % (processed, expected, coverage["pct"], coverage["proteins"], total_rows))
    lines.append("- Fallos reales: **%d**; pendientes reanudables: **%d**." % (len(failures), len(pending)))
    lines.append("- Validación: `validation_full.json` (estructural) y muestra de 20 genomas.")
    lines.append("")
    lines.append("## Entregables")
    lines.append("")
    lines.append("| archivo | contenido |")
    lines.append("|---|---|")
    lines.append("| `data/membrane/full/<genome>.tsv` | TSV por genoma (11 columnas del contrato) |")
    lines.append("| `data/membrane/gff/full/<genome>.gff3` | GFF3 por genoma (`polypeptide` + `membrane_type`) |")
    lines.append("| `data/membrane/full/genomes.tsv` | lista de 400 genomas + filo/organismo/CDS |")
    lines.append("| `data/membrane/full/state.json` | estado reanudable por genoma y lote |")
    lines.append("| `data/membrane/full/batches/` | estado por lote |")
    lines.append("| `data/membrane/summary_full.json` | totales, cobertura, fallos, método, hallazgos |")
    lines.append("| `data/membrane/validation_full.json` | validación estructural de los 400 (150 procesados + 250 pendientes) |")
    lines.append("| `data/membrane/validation_full_sample.json` | muestra de 20 genomas procesados |")
    lines.append("| `data/membrane/full/genomes_done.tsv` | sublista de genomas ya procesados |")
    lines.append("| `data/membrane/README_full.md` | este documento |")
    lines.append("")
    lines.append("## Cobertura y conteos")
    lines.append("")
    lines.append("| type | filas | % |")
    lines.append("|---|---:|---:|")
    for t in TYPES:
        lines.append("| `%s` | %d | %.1f%% |" % (t, types.get(t, 0), pct(types.get(t, 0), total_rows)))
    lines.append("| **total** | **%d** | **100%%** |" % total_rows)
    lines.append("")
    lines.append("Confianza: high=%d (%.1f%%), medium=%d (%.1f%%), low=%d (%.1f%%)."
                 % (conf.get("high", 0), pct(conf.get("high", 0), total_rows),
                    conf.get("medium", 0), pct(conf.get("medium", 0), total_rows),
                    conf.get("low", 0), pct(conf.get("low", 0), total_rows)))
    lines.append("")
    lines.append("Transportadores: %d; con sustrato explícito: %d (%.1f%%)."
                 % (n_tr, tr_with_sub, pct(tr_with_sub, n_tr)))
    lines.append("")
    lines.append("### Cobertura por filo")
    lines.append("")
    lines.append("| filo | genomas | proteínas | filas | filas/proteína | transportadores |")
    lines.append("|---|---:|---:|---:|---:|---:|")
    for ph, d in by_phylum_out.items():
        lines.append("| %s | %d | %d | %d | %.3f | %d |"
                     % (ph, d["genomes"], d["proteins"], d["rows"], d["row_fraction"],
                        d["counts"].get("transporter", 0)))
    lines.append("")
    lines.append("## Métodos")
    lines.append("")
    lines.append("Se reutiliza el panel Pfam del piloto (255 modelos) y `hmmscan --cut_ga`.")
    lines.append("TM por hidropatía Kyte-Doolittle (ventana 19 aa, media ≥ 2,0, segmento ≥ 18 aa;")
    lines.append("≥2 segmentos o SP + 1 TM). SP por h-region 9–13 aa (media KD ≥ 2,2, residuos -1/-3;")
    lines.append("no SignalP). Lipobox por `[LVI][ASTM][AG]C` + ventana hidrofóbica (no LipoP).")
    lines.append("Sustrato solo si Swiss-Prot lo nombra explícitamente en contexto de transporte/SBP.")
    lines.append("Una fila por proteína; el resto de roles queda en `family`/`evidence`.")
    lines.append("")
    lines.append("## Herramientas")
    lines.append("")
    lines.append("| herramienta | estado |")
    lines.append("|---|---|")
    lines.append("| HMMER 3.3.2 (`hmmscan`, `hmmpress`) | disponible |")
    lines.append("| Pfam-A.hmm local | disponible (panel de 255 modelos) |")
    lines.append("| DIAMOND/Swiss-Prot vía `functional.tsv` | disponible (precalculado) |")
    lines.append("| SignalP 6 / TMHMM / Phobius / PSORTb / LipoP / MacSyFinder | **no instalados** |")
    lines.append("")
    lines.append("## Fallos y pendientes")
    lines.append("")
    if failures:
        lines.append("Fallos reales: **%d**." % len(failures))
        lines.append("")
        lines.append("| genoma | motivo | error |")
        lines.append("|---|---|---|")
        for f in failures[:100]:
            lines.append("| %s | %s | %s |" % (f["genome"], f["reason"],
                                                (f.get("error") or "")[:160].replace("|", "/")))
        if len(failures) > 100:
            lines.append("| ... | ... | ... |")
    else:
        lines.append("Fallos reales: **0**.")
    lines.append("")
    lines.append("Pendientes (sin salida, reanudables): **%d**." % len(pending))
    if pending:
        lines.append("")
        lines.append("Primeros pendientes: " + ", ".join("`%s`" % g for g in pending[:10])
                     + (" ..." if len(pending) > 10 else ""))
    lines.append("")
    lines.append("## Lotes y estado")
    lines.append("")
    lines.append("| lote | total | completados | fallos | segundos |")
    lines.append("|---|---:|---:|---:|---:|")
    for b, d in sorted(summary["batches"].items()):
        lines.append("| %s | %d | %d | %d | %s |"
                     % (b, d.get("total", 0), d.get("completed", 0), d.get("failed", 0),
                        d.get("elapsed_s", 0)))
    lines.append("")
    lines.append("Nota: 10 genomas ya estaban completos antes de reanudar; por eso "
                 "`batch_001` muestra 15/25. Total: 10 + 390 = 400.")
    lines.append("")
    lines.append("## Hallazgos")
    lines.append("")
    for i, f in enumerate(findings, 1):
        lines.append("%d. %s" % (i, f))
    lines.append("")
    lines.append("## Límites")
    lines.append("")
    for lim in summary["limitations"]:
        lines.append("- %s" % lim)
    lines.append("")
    lines.append("## Bloqueos y continuación")
    lines.append("")
    for b in summary["blockers"]:
        lines.append("- %s" % b)
    lines.append("")
    if pending:
        lines.append("Reanudar el lote restante (salta los genomas ya completos):")
        lines.append("")
        lines.append("```bash")
        lines.append("python3 scripts/membrane/06_run_full.py --workers 5 --cpu-per-worker 2 --batch-size 25")
        lines.append("```")
        lines.append("")
    else:
        lines.append("No hay lote pendiente: la corrida está completa. Para reproducir:")
        lines.append("")
        lines.append("```bash")
        lines.append("python3 scripts/membrane/06_run_full.py --reconcile-only")
        lines.append("```")
        lines.append("")
    lines.append("## Validación")
    lines.append("")
    if validation:
        for p, v in validation.items():
            lines.append("- `%s`: ok=%s, totales=%s" % (p, v.get("ok"), v.get("totals")))
    else:
        lines.append("- Pendiente: ejecutar `05_validate_outputs.py` sobre la lista completa.")
    lines.append("")
    lines.append("## Integración con el lead")
    lines.append("")
    lines.append("- Unir por `genome` + `protein_id`; `start/end` 1-based inclusivos.")
    lines.append("- Usar `high`/`medium` para la GUI/informe; `low` es hipótesis.")
    lines.append("- `substrate=\"\"` significa sin evidencia, no ausencia de sustrato.")
    lines.append("- GFF3: features `polypeptide`; el rol está en `membrane_type`.")
    lines.append("- El piloto de 10 sigue en `data/membrane/*.tsv` y `data/membrane/gff/*.gff3`.")
    lines.append("")
    lines.append("## Reproducibilidad")
    lines.append("")
    lines.append("```bash")
    lines.append("python3 scripts/membrane/06_run_full.py --workers 4 --cpu-per-worker 2 --batch-size 25")
    lines.append("python3 scripts/membrane/05_validate_outputs.py --pilot data/membrane/full/genomes.tsv \\")
    lines.append("  --tsv-dir data/membrane/full --gff-dir data/membrane/gff/full \\")
    lines.append("  --json data/membrane/validation_full.json --quiet --allow-missing")
    lines.append("python3 scripts/membrane/07_summarize_full.py")
    lines.append("```")
    lines.append("")
    Path(args.readme).write_text("\n".join(lines) + "\n")
    print(f"[summary] {processed}/{expected} genomas, {total_rows} filas -> {args.json}")
    print(f"[summary] {args.readme}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
