"""Backend de la interfaz web: anota FASTA subidos por el usuario.

Sin dependencias de los datos del proyecto. Detecta qué herramientas hay y
ejecuta solo lo seleccionado. Pensado para distribuirse por GitHub.
"""
from __future__ import annotations

import csv
import gzip
import json
import io
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

# El motor de genes fantasma vive en scripts/ghost
_ROOT = Path(__file__).resolve().parents[1]
GHOST_PATH = None
for _cand in (_ROOT / "scripts" / "ghost", _ROOT):
    if (_cand / "ghost").is_dir():
        GHOST_PATH = _cand
        if str(_cand) not in sys.path:
            sys.path.insert(0, str(_cand))

from ghost.labels import SOURCES, label_row  # noqa: E402
from ghost.orfs import merge_intervals  # noqa: E402


# --------------------------------------------------------------------------
# Detección de herramientas
# --------------------------------------------------------------------------

def detect_tools() -> dict[str, bool]:
    def has(cmd: str) -> bool:
        return shutil.which(cmd) is not None

    try:
        import pyrodigal  # noqa: F401
        pyrodigal_ok = True
    except Exception:
        pyrodigal_ok = False
    try:
        import torch  # noqa: F401
        import esm  # noqa: F401
        esm_ok = True
    except Exception:
        esm_ok = False
    return {
        "pyrodigal": pyrodigal_ok,
        "trna": has("tRNAscan-SE"),
        "rrna": has("barrnap") or has("cmsearch"),
        "diamond": has("diamond"),
        "esm": esm_ok,
        "ghost": True,
    }


# --------------------------------------------------------------------------
# FASTA
# --------------------------------------------------------------------------

def parse_fasta(data: bytes) -> dict[str, str]:
    """Acepta FASTA plano o gzip; devuelve {contig: secuencia}."""
    if data[:2] == b"\x1f\x8b":
        data = gzip.decompress(data)
    text = data.decode("utf-8", errors="replace")
    contigs: dict[str, str] = {}
    name = None
    buf: list[str] = []
    for line in text.splitlines():
        if line.startswith(">"):
            if name is not None:
                contigs[name] = "".join(buf).upper()
            name = line[1:].split()[0]
            buf = []
        elif line.strip():
            buf.append(line.strip())
    if name is not None:
        contigs[name] = "".join(buf).upper()
    return contigs


def gc_content(contigs: dict[str, str]) -> float:
    g = c = tot = 0
    for s in contigs.values():
        g += s.count("G")
        c += s.count("C")
        tot += len(s)
    return (g + c) / tot if tot else 0.0


def safe_name(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]", "_", Path(name).stem) or "genome"


# --------------------------------------------------------------------------
# Anotadores
# --------------------------------------------------------------------------

def run_cds(contigs: dict[str, str]) -> tuple[list[dict], list[dict]]:
    import pyrodigal

    finder = pyrodigal.GeneFinder(meta=True)
    features: list[dict] = []
    proteins: list[dict] = []
    n = 0
    for contig, seq in contigs.items():
        for gene in finder.find_genes(seq.encode()):
            n += 1
            fid = f"cds_{n}"
            aa = gene.translate()
            features.append({
                "contig": contig, "type": "CDS", "start": int(gene.begin),
                "end": int(gene.end), "strand": "+" if gene.strand == 1 else "-",
                "phase": 0, "id": fid, "product": "hypothetical protein",
            })
            proteins.append({"id": fid, "aa": aa})
    return features, proteins


def run_trna(contigs: dict[str, str]) -> list[dict]:
    feats: list[dict] = []
    with tempfile.TemporaryDirectory() as td:
        fa = Path(td) / "in.fna"
        out = Path(td) / "trna.txt"
        with open(fa, "w") as fh:
            for name, seq in contigs.items():
                fh.write(f">{name}\n")
                for i in range(0, len(seq), 70):
                    fh.write(seq[i:i + 70] + "\n")
        subprocess.run(
            ["tRNAscan-SE", "-q", "-B", "-o", str(out), str(fa)],
            check=True, capture_output=True, timeout=1800,
        )
        for line in out.read_text().splitlines():
            p = line.split("\t")
            if len(p) < 9 or not p[0] or p[0].startswith("Sequence"):
                continue
            try:
                start, end = int(p[2]), int(p[3])
            except ValueError:
                continue
            feats.append({
                "contig": p[0], "type": "tRNA", "start": start, "end": end,
                "strand": "+" if end >= start else "-",
                "phase": 0, "id": f"trna_{len(feats)+1}",
                "product": f"tRNA-{p[4]} ({p[5]})",
            })
    return feats


def run_rrna(contigs: dict[str, str]) -> list[dict]:
    feats: list[dict] = []
    with tempfile.TemporaryDirectory() as td:
        fa = Path(td) / "in.fna"
        with open(fa, "w") as fh:
            for name, seq in contigs.items():
                fh.write(f">{name}\n")
                for i in range(0, len(seq), 70):
                    fh.write(seq[i:i + 70] + "\n")
        res = subprocess.run(
            ["barrnap", "--kingdom", "bac", "--threads", "1", str(fa)],
            check=True, capture_output=True, text=True, timeout=1800,
        )
        for line in res.stdout.splitlines():
            if line.startswith("#") or not line.strip():
                continue
            p = line.split("\t")
            if len(p) < 9 or p[2] != "rRNA":
                continue
            m = re.search(r"product=([^;]+)", p[8])
            feats.append({
                "contig": p[0], "type": "rRNA", "start": int(p[3]), "end": int(p[4]),
                "strand": p[6], "phase": 0, "id": f"rrna_{len(feats)+1}",
                "product": m.group(1) if m else "rRNA",
            })
    return feats


def run_ghost(contigs: dict[str, str], features: list[dict],
              min_region_len: int = 150, min_aa: int = 60,
              esm_top: int = 0) -> tuple[list[dict], dict]:
    from ghost.pipeline import build_candidates

    genome = "uploaded"
    with tempfile.TemporaryDirectory() as td:
        root = Path(td) / "data"
        (root / "raw").mkdir(parents=True)
        (root / "annotations_ab_initio" / genome).mkdir(parents=True)
        with gzip.open(root / "raw" / f"{genome}_genomic.fna.gz", "wt") as fh:
            for name, seq in contigs.items():
                fh.write(f">{name}\n")
                for i in range(0, len(seq), 70):
                    fh.write(seq[i:i + 70] + "\n")
        with open(root / "annotations_ab_initio" / genome / "annotation.gff3", "w") as fh:
            fh.write("##gff-version 3\n")
            for name, seq in contigs.items():
                fh.write(f"##sequence-region {name} 1 {len(seq)}\n")
            for f in features:
                fh.write(to_gff_line(f, "uploaded"))
        res = build_candidates(
            root, genome, min_region_len=min_region_len, min_aa=min_aa
        )
        if esm_top and res.candidates:
            try:
                from ghost.pipeline import add_esm_scores
                add_esm_scores(res.candidates, top_n=esm_top)
            except Exception:
                pass
        cands = sorted(res.candidates, key=lambda c: c["rank"], reverse=True)
        stats = res.stats()
        stats["tiers"] = {
            "high": sum(1 for c in cands if c["rank"] >= 0.75),
            "medium": sum(1 for c in cands if 0.55 <= c["rank"] < 0.75),
        }
    return cands, stats


# --------------------------------------------------------------------------
# GFF3
# --------------------------------------------------------------------------

def to_gff_line(f: dict, source: str = "uploaded") -> str:
    return "\t".join([
        f["contig"], source, f["type"], str(f["start"]), str(f["end"]),
        ".", f["strand"], str(f.get("phase", 0)),
        f"ID={f['id']};product={f.get('product', '')}",
    ]) + "\n"


def build_gff(contigs: dict[str, str], features: list[dict]) -> str:
    out = ["##gff-version 3\n"]
    for name, seq in contigs.items():
        out.append(f"##sequence-region {name} 1 {len(seq)}\n")
    for f in sorted(features, key=lambda x: (x["contig"], x["start"])):
        out.append(to_gff_line(f))
    return "".join(out)


def build_proteins_fasta(proteins: list[dict]) -> str:
    out = []
    for p in proteins:
        out.append(f">{p['id']}\n")
        aa = p["aa"]
        for i in range(0, len(aa), 60):
            out.append(aa[i:i + 60] + "\n")
    return "".join(out)


def candidates_tsv(cands: list[dict]) -> str:
    if not cands:
        return "contig\tstart\tend\tstrand\taa_len\tcodon_llr\trbs\trank\n"
    cols = ["contig", "start", "end", "strand", "aa_len", "start_codon",
            "codon_llr", "rbs", "esm_llr", "rank", "n_genomes", "label",
            "label_confidence", "aa"]
    lines = ["\t".join(cols)]
    for c in cands:
        lines.append("\t".join(str(c.get(k, "")) for k in cols))
    return "\n".join(lines) + "\n"


# --------------------------------------------------------------------------
# Desiertos de anotación
# --------------------------------------------------------------------------

def _gc_n(seq: str) -> tuple[float, float]:
    s = seq.upper()
    if not s:
        return 0.0, 0.0
    return (s.count("G") + s.count("C")) / len(s), s.count("N") / len(s)


def run_deserts(contigs: dict[str, str], features: list[dict],
                candidates: list[dict], min_len: int = 3000) -> list[dict]:
    """Tramos sin CDS/tRNA/rRNA. En subidas no hay oriC/móviles/defensa."""
    sig: dict[str, list[tuple[int, int]]] = {}
    for f in features:
        sig.setdefault(f["contig"], []).append((f["start"], f["end"]))
    rows: list[dict] = []
    for contig, seq in contigs.items():
        merged = merge_intervals(sig.get(contig, []))
        gaps: list[tuple[int, int]] = []
        cur = 1
        for s0, e0 in merged:
            if s0 - cur >= min_len:
                gaps.append((cur, s0 - 1))
            cur = max(cur, e0 + 1)
        if len(seq) - cur + 1 >= min_len:
            gaps.append((cur, len(seq)))
        for s0, e0 in gaps:
            inside = [c for c in candidates if c["contig"] == contig
                      and c["start"] >= s0 and c["end"] <= e0]
            g, n = _gc_n(seq[s0 - 1:e0])
            row = {
                "contig": contig, "start": s0, "end": e0, "length": e0 - s0 + 1,
                "gc": round(g, 3), "n_frac": round(n, 3),
                "ghost_orfs": len(inside),
                "max_rank": round(max((c["rank"] for c in inside), default=0.0), 3),
                "conserved_orfs": 0,
            }
            label_row(
                row, label="desert_region",
                source=("CDS/tRNA/rRNA (pyrodigal/tRNAscan/barrnap) + ORFs 09_ghost; "
                        "oriC, móviles, defensa y ARN no disponibles en genomas subidos"),
                evidence=(f"len={row['length']}; gc={row['gc']}; n={row['n_frac']}; "
                          f"orfs={row['ghost_orfs']}; max_rank={row['max_rank']}"),
                confidence="medium" if inside else "low",
                interpretation=("Hueco sin señales de anotación; sin datos de conservación "
                                "porque es un genoma único subido."),
            )
            rows.append(row)
    return rows


def compare_results(results: list[dict]) -> dict:
    """Deduplica candidatos entre los genomas subidos y etiqueta conservación."""
    if len(results) < 2:
        return {}
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        for r in results:
            doc = {"stats": {}, "candidates": r["candidates"]}
            (tmp / f"{r['genome']}_ghost_candidates.json").write_text(json.dumps(doc))
        out = tmp / "dedup"
        env = dict(os.environ)
        env["PYTHONPATH"] = str(GHOST_PATH) + os.pathsep + env.get("PYTHONPATH", "")
        subprocess.run(
            [sys.executable, "-m", "ghost.dedup", "--inputs", str(tmp / "*_ghost_candidates.json"),
             "--out", str(out), "--identity", "95", "--threads", "2"],
            check=True, capture_output=True, env=env,
        )
        cl = {r["cluster"]: int(r["n_genomes"])
              for r in csv.DictReader(open(out / "dedup_clusters.tsv"), delimiter="\t")}
        loc = {}
        for r in csv.DictReader(open(out / "dedup_members.tsv"), delimiter="\t"):
            loc[(r["genome"], r["locus"])] = r["cluster"]
    for r in results:
        for c in r["candidates"]:
            key = f"{c['contig']}:{c['start']}-{c['end']}{c['strand']}"
            cluster = loc.get((r["genome"], key))
            ng = cl.get(cluster, 1) if cluster else 1
            c["n_genomes"] = ng
            label_row(
                c, label="candidate_conserved" if ng >= 2 else "candidate_novel",
                source=("ORF de novo 09_ghost + " + SOURCES["conservation"]),
                evidence=(f"rank={c.get('rank')}; codon_llr={c.get('codon_llr')}; "
                          f"rbs={c.get('rbs')}; esm={c.get('esm_llr')}; ng={ng}"),
                confidence="high" if ng >= 2 else "medium",
                interpretation=("ORF conservado entre los genomas subidos: candidato a gen real."
                                if ng >= 2 else
                                "ORF de novo sin conservación entre los genomas subidos."),
            )
        for d in r.get("deserts", []):
            cons = sum(1 for c in r["candidates"]
                       if c["contig"] == d["contig"] and c["start"] >= d["start"]
                       and c["end"] <= d["end"] and c.get("n_genomes", 1) >= 2)
            d["conserved_orfs"] = cons
            if cons:
                d["label"] = "desert_candidate"
                d["label_confidence"] = "high"
                d["label_interpretation"] = ("Desierto con ORFs conservados entre los genomas "
                                             "subidos: candidato a gen/operón no anotado.")
    multi = sum(1 for v in cl.values() if v >= 2)
    return {"clusters": len(cl), "multi_genome": multi}


# --------------------------------------------------------------------------
# Orquestación
# --------------------------------------------------------------------------



def coverage_bp(contigs: dict[str, str], features: list[dict]) -> int:
    from collections import defaultdict
    per: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for f in features:
        per[f["contig"]].append((f["start"], f["end"]))
    total = 0
    for iv in per.values():
        for s, e in merge_intervals(iv):
            total += e - s + 1
    return total


def build_report(res: dict) -> str:
    from collections import Counter
    types = Counter(f["type"] for f in res.get("features", []))
    lines = [
        f"# Informe de anotación — {res.get('name','')}",
        "",
        f"- Genoma: **{res.get('genome','')}**",
        f"- Tamaño: {res.get('bp',0):,} pb · GC {res.get('gc',0)*100:.1f}%",
        f"- Cobertura anotada: **{res.get('coverage',0)*100:.1f}%** "
        f"({res.get('annotated_bp',0):,} pb)",
        f"- Features: CDS {types.get('CDS',0)}, tRNA {types.get('tRNA',0)}, "
        f"rRNA {types.get('rRNA',0)}",
        f"- Regiones no anotadas: {res.get('stats',{}).get('regions','n/d')}",
        f"- Desiertos ≥3 kb: {len(res.get('deserts',[]))}",
        f"- ORFs candidatos: {len(res.get('candidates',[]))}",
    ]
    lines += ["", "## Desiertos"]
    ds = res.get("deserts", [])
    if ds:
        lines += ["| contig | inicio | fin | pb | GC | ORFs | etiqueta | confianza |",
                  "|---|---|---|---|---|---|---|---|"]
        for d in ds:
            lines.append(f"| {d.get('contig')} | {d.get('start')} | {d.get('end')} | "
                         f"{d.get('length')} | {d.get('gc')} | {d.get('ghost_orfs')} | "
                         f"{d.get('label')} | {d.get('label_confidence')} |")
    else:
        lines.append("_Sin desiertos ≥3 kb._")
    lines += ["", "## Candidatos (top 20)", ""]
    cands = sorted(res.get("candidates", []), key=lambda c: c.get("rank", 0), reverse=True)[:20]
    if cands:
        lines += ["| contig | inicio | fin | aa | rank | etiqueta | confianza |",
                  "|---|---|---|---|---|---|---|"]
        for c in cands:
            lines.append(f"| {c.get('contig')} | {c.get('start')} | {c.get('end')} | "
                         f"{c.get('aa_len')} | {c.get('rank',0):.3f} | {c.get('label')} | "
                         f"{c.get('label_confidence')} |")
    else:
        lines.append("_Sin candidatos._")
    met = res.get("metabolic")
    if met:
        gaps = met.get("gaps", [])
        origins = Counter(g.get("origin", "") for g in gaps)
        lines += ["", "## Brechas metabólicas",
                  f"- Red: {met.get('source','')} · {met.get('reactions',0)} reacciones",
                  f"- Brechas detectadas: {len(gaps)}"]
        for k, v in origins.most_common():
            lines.append(f"- {k}: {v}")
        tp = Counter(g.get("transport_status", "") for g in gaps)
        dr = Counter(g.get("direction_status", "") for g in gaps)
        lines += ["",
                  "### Transporte",
                  f"- {tp.get('t1_bacterial',0)} importaciones T1 (homólogo bacteriano)",
                  f"- {tp.get('t1_review',0)} en revisión (eucariota/baja identidad)",
                  f"- {tp.get('transport_reaction_only',0)} solo con reacción de transporte",
                  f"- {tp.get('no_transport',0)} sin transporte",
                  "",
                  "### Dirección",
                  f"- {dr.get('direction_unspecified',0)} sin dirección (Rhea master UN)",
                  f"- {dr.get('thermo_forward',0)} favorables (ΔG)",
                  f"- {dr.get('thermo_reverse',0)} desfavorables (ΔG)",
                  "",
                  "### Evidencia por brecha (top 15)",
                  "| gap | soporte | rol | origen | transporte | dirección |",
                  "|---|---|---|---|---|---|"]
        from family_logic import explain_families
        fam_rows = explain_families(gaps)[:12]
        lines += ["", "### Familias de rutas (explicación)", "",
                  "| familia | brechas | homólogo | transporte | explicación | acción |",
                  "|---|---|---|---|---|---|"]
        for fr in fam_rows:
            lines.append(f"| {fr['family']} | {fr['gaps']} | {fr['homolog']} | "
                         f"{fr['transport']} | {fr['conclusion']} | {fr['action']} |")
        hom = Counter(g.get("homolog_status", "") for g in gaps)
        rc = Counter(g.get("review_class", "") for g in gaps)
        lines += ["", "### Clases de revisión",
                  f"- transport: {rc.get('transport',0)}",
                  f"- annotation_gap: {rc.get('annotation_gap',0)}",
                  f"- unknown: {rc.get('unknown',0)}"]
        ex_gaps = [g for g in gaps if g.get("known_exception")]
        if ex_gaps:
            lines += ["", "### Excepciones documentadas (aviso, no exención)"]
            for g in ex_gaps[:10]:
                lines.append(f"- **{g.get('gap')}**: {g.get('known_exception')} — {g.get('exception_note','')}")
        lines += ["", "### Enzima productora",
                  f"- homólogo fuerte/débil: {hom.get('homolog_strong',0)}/{hom.get('homolog_weak',0)}",
                  f"- sin homólogo: {hom.get('no_homolog',0)}",
                  f"- no comprobado (fuera de la base): {hom.get('not_checked',0)}"]
        sug = [g for g in gaps if g.get("homolog_rhea")]
        if sug:
            lines += ["", "### Sugerencias de mapeo (revisables)", "",
                      "| gap | ORF | homólogo | Rhea | EC | estado |",
                      "|---|---|---|---|---|---|"]
            for g in sug[:15]:
                lines.append(f"| {g.get('gap','')} | {g.get('homolog_query','') or g.get('homolog_acc','')} | "
                             f"{g.get('homolog_acc','')} | {g.get('homolog_rhea','')} | "
                             f"{g.get('homolog_ec','')} | {g.get('mapping_suggestion','suggested')} |")
        for g in sorted(gaps, key=lambda x: int(x.get("support", 0) or 0), reverse=True)[:15]:
            lines.append(f"| {g.get('gap','')} | {g.get('support',0)} | {g.get('role','')} | "
                         f"{g.get('origin','')} | {g.get('transport_status','')} | {g.get('direction_status','')} |")
        lines += ["",
                  "### Límites",
                  "- `direction_unspecified` no es una conclusión: falta ΔG o anotación fisiológica.",
                  "- Transporte se apoya en Rhea + homología; sin homólogo no se marca importación.",
                  "- Cierres por dominio no se calculan en subidos; requieren el catálogo de los 400.",
                  "- Fuentes: Rhea, 400 genomas, Pfam/HMM, eQuilibrator/TECRDB."]
    ec = res.get("esmfold_closure")
    if ec:
        lines += ["", "## Cierre por estructura (ESMFold)",
                  "- Predicción computacional; **no** es evidencia experimental ni cierre confirmado.",
                  "", "| gap | ORF | proveedor | job | estado | pLDDT | pTM | PDB |",
                  "|---|---|---|---|---|---|---|---|"]
        for row in ec[:15]:
            lines.append(f"| {row.get('gap','')} | {row.get('query','')} | {row.get('provider','')} | "
                         f"{row.get('job_id','')} | {row.get('status','')} | {row.get('plddt','')} | "
                         f"{row.get('ptm','')} | {row.get('pdb','')} |")
        if res.get("esmfold_closure_error"):
            lines.append(f"- Error: {res['esmfold_closure_error']}")
    hy = res.get("structural_hypotheses")
    if hy:
        lines += ["", "## Hipótesis estructurales revisables",
                  "- Capa separada de la anotación base; no confirma función.",
                  "", "| proteína | hipótesis | EC | cofactor | confianza | revisión |",
                  "|---|---|---|---|---|---|"]
        for row in hy[:20]:
            lines.append(f"| {row.get('protein_id','')} | {row.get('hypothesis','')} | "
                         f"{row.get('ec','') or '—'} | {row.get('cofactor','') or '—'} | "
                         f"{row.get('confidence','')} | {row.get('review_decision','pending')} |")
    ff = res.get("five_filters")
    if ff:
        lines += ["", "## 5 filtros estructurales",
                  "- Foldseek + EnzyMM + P2Rank + AlphaFill + sitio activo.",
                  "- Prioridad, no cierre.",
                  "", "| PDB | F1 | F2 | F3 | F4 | F5 | pasan | prioridad |",
                  "|---|---|---|---|---|---|---:|---|"]
        for row in ff[:15]:
            lines.append(f"| {row.get('pdb','')} | {row.get('f1_foldseek','')} | {row.get('f2_enzymm','')} | "
                         f"{row.get('f3_pocket','')} | {row.get('f4_alphafill','')} | {row.get('f5_active_site','')} | "
                         f"{row.get('n_pass','')} | {row.get('priority','')} |")
    enz = res.get("enzymm_scan")
    if enz:
        lines += ["", "## Anotación EnzyMM/PyJess (M-CSA)",
                  "- Match geométrico de plantillas catalíticas; no confirma función.",
                  "", "| query | template | M-CSA | EC | RMSD | matched |",
                  "|---|---|---|---|---:|---|"]
        for row in enz[:15]:
            lines.append(f"| {row.get('query_id','')} | {row.get('template_pdb_id','')} | "
                         f"{row.get('template_mcsa_id','')} | {row.get('template_ec','')} | "
                         f"{row.get('rmsd','')} | {row.get('matched_residues','')} |")
    if res.get("enzymm_error"):
        lines.append(f"- EnzyMM error: {res['enzymm_error']}")
    se = res.get("structural_evidence")
    if se:
        lines += ["", "## Evidencia estructural multi-fuente",
                  "- Foldseek PDB100 + AlphaFill + sitios activos locales.",
                  "- `multi_source_support` no es cierre confirmado.",
                  "", "| gap | ORF | Foldseek | AlphaFill | active cons | fuentes | veredicto |",
                  "|---|---|---|---|---:|---:|---|"]
        for row in se[:15]:
            lines.append(f"| {row.get('gap','')} | {row.get('query','')} | {row.get('foldseek_top','')} | "
                         f"{row.get('alphafill_ligands','') or '—'} | {row.get('active_site_conserved','')}/"
                         f"{row.get('active_site_rows','')} | {row.get('n_sources','')} | {row.get('verdict','')} |")
    lines += ["", "## Fuentes",
              "- CDS: pyrodigal; tRNA: tRNAscan-SE; rRNA: barrnap/cmsearch.",
              "- ORFs: motor 09_ghost; cada fila declara fuente, evidencia y confianza.",
              "- Rhea (reacciones), 400 genomas (soporte), Pfam/HMM (dominio)."]
    return "\n".join(lines) + "\n"


def annotate(name: str, data: bytes, options: dict) -> dict:
    """Anota un FASTA según las casillas seleccionadas."""
    contigs = parse_fasta(data)
    if not contigs:
        raise ValueError("El archivo no contiene secuencias FASTA válidas.")
    features: list[dict] = []
    proteins: list[dict] = []
    stats: dict = {}

    if options.get("cds"):
        f, p = run_cds(contigs)
        features += f
        proteins += p
        stats["cds"] = len(f)
    if options.get("trna"):
        f = run_trna(contigs)
        features += f
        stats["trna"] = len(f)
    if options.get("rrna"):
        f = run_rrna(contigs)
        features += f
        stats["rrna"] = len(f)
    if options.get("ghost"):
        cands, gstats = run_ghost(
            contigs, features,
            min_region_len=options.get("min_region_len", 150),
            min_aa=options.get("min_aa", 60),
            esm_top=options.get("esm_top", 0) if options.get("esm") else 0,
        )
        stats["ghost_orfs"] = gstats.get("n_orfs", 0)
        stats["ghost_high"] = gstats.get("tiers", {}).get("high", 0)
    else:
        cands = []

    gname = safe_name(name)
    for c in cands:
        c["genome"] = gname
    deserts = run_deserts(contigs, features, cands) if options.get("deserts") else []
    for c in cands:
        label_row(c, label="candidate_novel", source="ORF de novo 09_ghost",
                  evidence=(f"rank={c.get('rank')}; codon_llr={c.get('codon_llr')}; "
                            f"rbs={c.get('rbs')}; esm={c.get('esm_llr')}"),
                  confidence="medium",
                  interpretation="ORF de novo sin conservación evaluada aún (un solo genoma).")
    bp = sum(len(x) for x in contigs.values())
    ann = coverage_bp(contigs, features)
    res = {
        "name": name,
        "genome": gname,
        "contigs": len(contigs),
        "bp": sum(len(s) for s in contigs.values()),
        "gc": round(gc_content(contigs), 4),
        "features": features,
        "proteins": proteins,
        "candidates": cands,
        "deserts": deserts,
        "stats": stats,
        "gff": build_gff(contigs, features),
        "proteins_fasta": build_proteins_fasta(proteins),
        "candidates_tsv": candidates_tsv(cands),
        "annotated_bp": ann,
        "coverage": round(ann / bp, 4) if bp else 0.0,
    }
    res["report_md"] = build_report(res)
    return res


def deserts_tsv(deserts: list[dict]) -> str:
    if not deserts:
        return "contig\tstart\tend\tlength\tgc\tghost_orfs\tlabel\tlabel_confidence\n"
    cols = ["contig", "start", "end", "length", "gc", "n_frac", "ghost_orfs",
            "max_rank", "conserved_orfs", "label", "label_confidence",
            "label_interpretation"]
    lines = ["\t".join(cols)]
    for d in deserts:
        lines.append("\t".join(str(d.get(k, "")) for k in cols))
    return "\n".join(lines) + "\n"
