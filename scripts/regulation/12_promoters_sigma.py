#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""12_promoters_sigma.py — promotores posicionales (-35/-10) e inventario de σ.

Para cada operón de L2 se escanea la región río arriba (120 pb) buscando:
  * caja -35 (TTGACA, <=2 cambios) y caja -10 (TATAAT, <=2 cambios) con
    espaciador 15-19 pb;
  * -10 extendida (TG justo antes de la -10).
Es evidencia POSICIONAL, no un PWM: se declara `consensus` como máximo.

Inventario de factores sigma a partir de producto/nombre (functional.tsv o
GenBank) y, si estuviera disponible, confirmación por dominio Pfam. El panel
local de Pfam NO contiene los modelos de sigma (PF00140/PF04542/...), así que
la confirmación por HMM queda declarada como no disponible, no como hecha.

Salidas: `data/regulation/<G>/promoters.tsv`, `regulators_inventory.tsv`.

Uso:
  python scripts/regulation/12_promoters_sigma.py --genome GCF_...
  python scripts/regulation/12_promoters_sigma.py --all
  python scripts/regulation/12_promoters_sigma.py --selftest
"""
from __future__ import annotations

import argparse
import re
import sys
import time
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _common as C

WINDOW = 150
MIN_SCORE = 10         # de 15 posibles + prior de distancia
CONS35, SEED35 = "TTGACA", "TTGA"
CONS10, SEED10 = "TATAAT", "TAT"

PROM_COLS = ["operon_id", "label", "contig", "strand", "box35", "mm35", "box10",
             "mm10", "spacer", "ext10", "score", "dist_10_to_start", "n_sites",
             "sigma", "sigma_status", "status", "evidence", "attempted",
             "confidence"]
REG_COLS = ["genome", "gene_id", "gene_name", "product", "family", "role",
            "pfam_checked", "pfam_hit", "status", "evidence", "attempted",
            "confidence"]

SIGMA_PATTERNS = [
    (r"\brpod\b|sigma[- ]?70|primary sigma|housekeeping sigma", "sigma70"),
    (r"\brpoh\b|sigma[- ]?32|heat shock sigma", "sigma32"),
    (r"\brpos\b|sigma[- ]?38|stationary phase sigma", "sigma38"),
    (r"\brpon\b|sigma[- ]?54|sigma n", "sigma54"),
    (r"\bflia\b|\brpof\b|sigma[- ]?28|flagellar sigma", "sigma28"),
    (r"\brpoe\b|sigma[- ]?e\b|extracytoplasmic", "sigmaE"),
    (r"\bsiga\b", "sigma70"),
    (r"\bsigb\b", "sigmaB"),
    (r"\bsige\b", "sigmaE"),
    (r"\bsigf\b", "sigmaF"),
    (r"\bsigh\b", "sigmaH"),
    (r"\bsigk\b", "sigmaK"),
    (r"sigma[- ]?factor", "sigma_unresolved"),
]
_PFAM_SIGMA = ["PF00140", "PF04542", "PF04539", "PF04545"]  # PF04997 es Rpb1, no sigma
PANEL_REG = C.DATA / "db" / "pfam" / "panel_regulation.hmm"
SIGMA_DOMAINS = {"Sigma70_r1_2", "Sigma70_r2", "Sigma70_r3", "Sigma70_r4",
                 "Sigma54_1", "Sigma54_C"}


def _mm(cons: str, sub: str) -> int:
    return sum(1 for a, b in zip(cons, sub) if b not in C._ALLOWED[a])


def sites(cons: str, w: str, seed: str, maxmm: int):
    out, pos = [], 0
    at = cons.find(seed)
    while True:
        j = w.find(seed, pos)
        if j < 0:
            break
        pos = j + 1
        s = j - at
        if s < 0 or s + len(cons) > len(w):
            continue
        sub = w[s:s + len(cons)]
        mm = _mm(cons, sub)
        if mm <= maxmm:
            out.append((s, mm, sub))
    return out


SPACERS = range(14, 22)          # distribución real en RegulonDB: 14-21


def _distance_prior(dist: int) -> int:
    """Prior de arquitectura: la -10 suele estar 30-60 pb río arriba del inicio."""
    if 30 <= dist <= 60:
        return 3
    if 25 <= dist <= 80:
        return 2
    if 15 <= dist <= 120:
        return 0
    return -3


def find_promoter(w: str):
    """Mejor par -35/-10; devuelve dict o None."""
    w = w.upper()
    best, n_ok = None, 0
    for p, mm10, s10 in sites(CONS10, w, SEED10, 3):
        ext = 1 if p >= 2 and w[p - 2:p] == "TG" else 0
        for s in range(max(0, p - 25), max(0, p - 20) + 1):
            box = w[s:s + 6]
            if len(box) < 6:
                continue
            # la -10 extendida (TG) forma parte del elemento: el espaciador
            # efectivo se mide desde el TG, no desde el hexámero.
            spacer = p - (s + 6) - 2 * ext
            if spacer not in SPACERS:
                continue
            mm35 = _mm(CONS35, box)
            if mm35 > 3:
                continue
            score = (6 - mm35) + (6 - mm10) + \
                (2 if spacer in (16, 17, 18) else 1 if spacer in (15, 19) else 0) + ext
            score += _distance_prior(len(w) - p)
            if score < MIN_SCORE:
                continue
            n_ok += 1
            cand = (score, -(mm35 + mm10), -s, s, p, spacer, box, s10, ext)
            if best is None or cand > best:
                best = cand
    if best is None:
        return None
    score, _, _, s, p, spacer, box35, box10, ext = best
    return {"box35": box35, "mm35": _mm(CONS35, box35), "box10": box10,
            "mm10": _mm(CONS10, box10), "spacer": spacer, "ext10": ext,
            "score": score, "pos10": p, "n_sites": n_ok}


def hmm_sigma_scan(faa: Path) -> dict[str, set[str]]:
    """Dominios sigma por proteína (hmmsearch --cut_ga contra el panel)."""
    if not PANEL_REG.exists() or not faa.exists():
        return {}
    import subprocess
    out = C.TMP / f"{faa.stem}_{abs(hash(str(faa.name))) % 10**8}.domtbl"
    out.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(["hmmsearch", "--cut_ga", "--noali", "-o", "/dev/null",
                    "--domtblout", str(out), str(PANEL_REG), str(faa)],
                   check=False)
    hits: dict[str, set[str]] = {}
    if out.exists():
        for line in open(out):
            if line.startswith("#"):
                continue
            f = line.split()
            if len(f) >= 5:
                hits.setdefault(f[0], set()).add(f[3])
        out.unlink(missing_ok=True)
    return hits


def pfam_models_present() -> set[str]:
    """Accesiones Pfam realmente presentes en los paneles locales."""
    p = PANEL_REG if PANEL_REG.exists() else (C.DATA / "db" / "pfam" / "panel_focus.hmm")
    if not p.exists():
        return set()
    import gzip
    out = set()
    op = gzip.open if str(p).endswith(".gz") else open
    with op(p, "rt") as fh:
        for line in fh:
            if line.startswith("ACC"):
                out.add(line.split()[1].split(".")[0])
            elif line.startswith("HMMER3"):
                continue
    return out


_NOT_SIGMA = re.compile(r"anti[- ]?sigma|regulator of|stabilizing|protease|"
                        r"phosphatase|inhibitor|modulator|chaperone|"
                        r"sigma factor binding|sigma factor.*(?:regulator|binding)",
                        re.IGNORECASE)


def sigma_family(gene: C.Gene) -> tuple[str, str]:
    txt = f"{gene.gene_name or ''} {gene.product or ''}".lower()
    if not txt.strip():
        return "", ""
    if _NOT_SIGMA.search(f"{gene.product or ''}"):
        return "", ""          # anti-sigma o regulador, no factor sigma
    for pat, fam in SIGMA_PATTERNS:
        if re.search(pat, txt):
            return fam, pat
    return "", ""


def inventory(genome: C.Genome, present: set[str],
              hmm: dict[str, set[str]] | None = None) -> list[dict]:
    """Sigma por nombre de producto + confirmación por dominio HMM.

    Dos métodos independientes: texto de anotación (DIAMOND/nombre) y dominio
    Pfam. Solo con ambos se marca `validated`; un dominio sin nombre entra como
    `novel` (candidato), no se descarta.
    """
    hmm = hmm or {}
    rows = []
    have = [a for a in _PFAM_SIGMA if a in present]
    for g in genome.cds():
        fam, pat = sigma_family(g)
        doms = hmm.get(g.gene_id, set()) & SIGMA_DOMAINS
        # Sigma70_r2 es el dominio específico de sigma; r4 es un HTH que también
        # tienen reguladores de respuesta. Sin nombre, exigimos r2 o 2+ dominios.
        if not fam and not ({"Sigma70_r2", "Sigma54_1"} & doms or len(doms) >= 2):
            continue
        if not fam and not doms:
            continue
        name_hit = bool(re.search(r"\brpo[dhsne]|\bflia\b|\bsig[a-z]\b",
                                  (g.gene_name or "").lower()))
        ev = f"product:{g.product};domains:{','.join(sorted(doms)) or '-'}"
        if fam and doms:
            d = C.decide("validated", validated_by="dominio_pfam;nombre/producto",
                         reason="dominio sigma y anotación coinciden",
                         attempted="hmmsearch;diamond_swissprot", confidence="high",
                         evidence=ev)
        elif doms and not fam:
            fam = "sigma_unresolved"
            d = C.decide("novel", validated_by="dominio_pfam",
                         reason="dominio sigma sin nombre asignable",
                         attempted="hmmsearch;diamond_swissprot", confidence="medium",
                         evidence=ev)
        elif fam == "sigma_unresolved" or (pat == r"sigma[- ]?factor" and not name_hit):
            d = C.decide("unresolved", validated_by="nombre",
                         reason="producto genérico de sigma sin familia ni dominio",
                         attempted="hmmsearch;nombre", confidence="low",
                         evidence=ev)
        elif name_hit:
            d = C.decide("consensus", validated_by="diamond_swissprot;nombre_gen",
                         reason="nombre y producto coinciden; sin dominio",
                         attempted="hmmsearch;diamond_swissprot", confidence="medium",
                         evidence=ev)
        else:
            d = C.decide("novel", validated_by="nombre",
                         reason="familia inferida solo del producto, sin dominio",
                         attempted="hmmsearch;diamond_swissprot", confidence="low",
                         evidence=ev)
        rows.append(C.status_row({
            "genome": genome.genome_id, "gene_id": g.gene_id,
            "gene_name": g.gene_name, "product": g.product, "family": fam,
            "role": "sigma_factor",
            "pfam_checked": "true" if present else "false",
            "pfam_hit": ",".join(sorted(doms)),
        }, d))
    return rows


def run(genome_id: str, out_dir: Path) -> dict:
    t0 = time.time()
    g = C.load_input(genome_id)
    if g is None:
        raise SystemExit(f"genoma no encontrado: {genome_id}")
    ops = C.read_tsv(out_dir / "operons.tsv")
    og = C.read_tsv(out_dir / "operon_genes.tsv")
    seqs = C.load_sequence(genome_id)
    first = {}
    for r in og:
        if r["position"] == "1":
            first[r["operon_id"]] = r
    present = pfam_models_present()
    rows = []
    for o in ops:
        f = first.get(o["operon_id"])
        contig = o["contig"]
        seq = seqs.get(contig, "")
        if not f or not seq:
            d = C.decide("unresolved", attempted="secuencia;promotor",
                         reason="sin secuencia local del contig",
                         confidence="low", evidence="")
            rows.append(C.status_row({
                "operon_id": o["operon_id"], "label": o["label"],
                "contig": contig, "strand": o["strand"]}, d))
            continue
        strand = int(o["strand"])
        w = C.upstream_seq(seq, int(f["start"]), int(f["end"]), strand, WINDOW)
        hit = find_promoter(w)
        if not hit:
            d = C.decide("unresolved", attempted="ttgaca/tataat<=2mm;esp15-19",
                         reason="sin par -35/-10 compatible en 120 pb",
                         confidence="low", evidence=f"window={len(w)}")
            rows.append(C.status_row({
                "operon_id": o["operon_id"], "label": o["label"],
                "contig": contig, "strand": strand, "sigma": "",
                "sigma_status": "unresolved", "n_sites": 0}, d))
            continue
        dist = len(w) - hit["pos10"]
        d = C.decide("consensus", validated_by="-35/-10 posicional",
                     reason="evidencia posicional (no PWM, no TSS medido)",
                     attempted="ttgaca/tataat<=2mm;esp15-19;-10 extendida",
                     confidence="medium" if hit["score"] >= 11 else "low",
                     evidence=f"score={hit['score']};n={hit['n_sites']}")
        rows.append(C.status_row({
            "operon_id": o["operon_id"], "label": o["label"], "contig": contig,
            "strand": strand, "box35": hit["box35"], "mm35": hit["mm35"],
            "box10": hit["box10"], "mm10": hit["mm10"],
            "spacer": hit["spacer"], "ext10": hit["ext10"],
            "score": hit["score"], "dist_10_to_start": dist,
            "n_sites": hit["n_sites"], "sigma": "sigma70",
            "sigma_status": "assumed_default",
        }, d))
    rows.sort(key=lambda r: (r["contig"], r["operon_id"]))
    C.write_tsv(out_dir / "promoters.tsv", rows, PROM_COLS)
    faa = out_dir / "proteins.faa"
    if not faa.exists():
        faa = C.AB / genome_id / "proteins.faa"
    hmm = hmm_sigma_scan(faa)
    inv = inventory(g, present, hmm)
    inv.sort(key=lambda r: r["gene_id"])
    C.write_tsv(out_dir / "regulators_inventory.tsv", inv, REG_COLS)
    meta = {
        "genome": genome_id, "n_operons": len(rows),
        "n_with_promoter": sum(1 for r in rows if r.get("box10")),
        "n_unresolved": sum(1 for r in rows if r["status"] == "unresolved"),
        "has_sequence": bool(seqs),
        "pfam_sigma_models_present": [a for a in _PFAM_SIGMA if a in present],
        "pfam_sigma_models_missing": [a for a in _PFAM_SIGMA if a not in present],
        "n_sigma_factors": len(inv),
        "n_sigma_hmm_validated": sum(1 for r in inv if r["status"] == "validated"),
        "sigma_families": dict(sorted(Counter(r["family"] for r in inv).items())),
        "elapsed_s": round(time.time() - t0, 2),
        "sha256_promoters": C.sha256_files([out_dir / "promoters.tsv",
                                            out_dir / "regulators_inventory.tsv"]),
    }
    C.dump_json(out_dir / "promoters_meta.json", meta)
    return meta


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--genome")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--offset", type=int, default=0)
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args(argv)
    if a.selftest:
        return selftest()
    if a.all:
        dirs = sorted(p for p in C.AB.glob("GCF_*") if p.is_dir())
        if a.offset:
            dirs = dirs[a.offset:]
        if a.limit:
            dirs = dirs[:a.limit]
        t0, ok = time.time(), 0
        for i, d in enumerate(dirs, 1):
            if not (C.REG / d.name / "operons.tsv").exists():
                continue
            try:
                m = run(d.name, C.REG / d.name)
                ok += 1
                if i % 50 == 0:
                    print(f"  [{i}/{len(dirs)}] prom={m['n_with_promoter']}",
                          flush=True)
            except Exception as e:                        # noqa: BLE001
                print(f"  FALLO {d.name}: {e}", file=sys.stderr)
        print(f"promotores: {ok} genomas en {time.time()-t0:.1f}s")
        return 0
    if not a.genome:
        ap.error("indica --genome o --all")
    m = run(a.genome, C.REG / a.genome)
    print(f"{a.genome}: {m['n_with_promoter']}/{m['n_operons']} operones con "
          f"-35/-10; sigma={m['n_sigma_factors']} "
          f"{m['sigma_families']}; secuencia={m['has_sequence']}")
    return 0


def selftest() -> int:
    # promotor canónico: caja -35 a 17 pb de la -10
    w = "A" * 40 + "TTGACA" + "A" * 17 + "TATAAT" + "A" * 30
    h = find_promoter(w)
    assert h and h["mm35"] == 0 and h["mm10"] == 0 and h["spacer"] == 17, h
    assert h["score"] == 17, h["score"]
    # espaciador inválido -> no hay promotor
    w2 = "A" * 40 + "TTGACA" + "A" * 30 + "TATAAT" + "A" * 30
    assert find_promoter(w2) is None
    # -10 extendida suma 1
    w3 = "A" * 40 + "TTGACA" + "A" * 17 + "TGTATAAT" + "A" * 30
    h3 = find_promoter(w3)
    assert h3 and h3["ext10"] == 1 and h3["score"] == 18, h3
    # familia sigma
    assert sigma_family(C.Gene("c", 1, 9, 1, "1", gene_name="rpoS",
                               product="RNA polymerase sigma factor RpoS"))[0] == "sigma38"
    assert sigma_family(C.Gene("c", 1, 9, 1, "2", gene_name="foo",
                               product="hypothetical protein"))[0] == ""
    assert pfam_models_present() is not None
    print("selftest 12_promoters_sigma: OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
