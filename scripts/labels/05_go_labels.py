#!/usr/bin/env python3
"""Etiquetas de proceso desde GO OBO + GAF (K-12/B. subtilis)."""
from __future__ import annotations
import argparse, csv, gzip, json, re, sys
from collections import Counter, defaultdict
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))


def read_obo(path: Path):
    names, parents = {}, defaultdict(set)
    cur = None
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line == "[Term]":
            cur = None
        elif line.startswith("id: GO:"):
            cur = line.split("id: ", 1)[1].strip()
        elif cur and line.startswith("name: "):
            names[cur] = line.split("name: ", 1)[1].strip()
        elif cur and line.startswith("is_a: "):
            parents[cur].add(line.split("is_a: ", 1)[1].split("!")[0].strip())
    return names, parents


def ancestors(go, parents, cache={}):
    if go in cache: return cache[go]
    out = {go}; stack = [go]
    while stack:
        x = stack.pop()
        for p in parents.get(x, ()):
            if p not in out:
                out.add(p); stack.append(p)
    cache[go] = out
    return out


def read_gaf(path: Path, exclude_iea=True):
    out = defaultdict(set)
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!"): continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9: continue
            obj, go, ev = f[1], f[4], f[6]
            if exclude_iea and ev == "IEA": continue
            out[obj].add(go)
    return out


LABELS = {
 "heat_shock": r"response to heat|heat shock|protein folding|chaperon",
 "flagellum": r"flagell|chemotaxis",
 "nitrogen": r"nitrogen|nitrate|nitrite|glutamine",
 "stationary": r"stationary|starvation",
 "envelope": r"cell wall|outer membrane|envelope|peptidoglycan|lipopolysaccharide",
 "sporulation": r"sporulation|spore",
}


def label_go(names, go_set):
    labels = set()
    for go in go_set:
        n = names.get(go, "").lower()
        for lab, pat in LABELS.items():
            if n and re.search(pat, n):
                labels.add(lab)
    return labels


def operon_map_k12(data):
    out = defaultdict(set); g2o = {}
    for op in csv.DictReader(open(data/'regulation/ecoli_k12/operons.tsv'), delimiter='\t'):
        for gid in (op.get('gene_ids') or '').split(','):
            if gid: g2o[gid] = op['operon_id']
    return g2o


def operon_map_bsu(data):
    hits = json.loads((data.parent/'results/subtilis/subtilis.functional_annotation.json').read_text()).get('hits', {})
    # cds id -> uniprot; operons use cds ids directly
    return hits


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--data-root', type=Path, default=ROOT/'data')
    ap.add_argument('--ext', type=Path, default=ROOT/'data'/'labels'/'external')
    ap.add_argument('--out', type=Path, default=ROOT/'reports'/'labels_pilot')
    a = ap.parse_args(); a.out.mkdir(parents=True, exist_ok=True)
    names, parents = read_obo(a.ext/'go'/'go-basic.obo')
    cases = {
      'ecoli_k12': (a.data_root/'regulation/ecoli_k12/genes.tsv', 'acc', read_gaf(a.ext/'goa'/'ECOLI-uniprot.gaf.gz')),
      'subtilis': (a.data_root.parent/'results/subtilis/subtilis.functional_annotation.json', 'uniprot', read_gaf(a.ext/'goa'/'BACSU-uniprot.gaf.gz')),
    }
    summary = {}
    # K12
    genes = list(csv.DictReader(open(cases['ecoli_k12'][0]), delimiter='\t'))
    gaf = cases['ecoli_k12'][2]
    g2o = operon_map_k12(a.data_root)
    gene_labels = {}
    for g in genes:
        acc = (g.get('acc') or '').strip()
        gos = gaf.get(acc, set())
        labs = label_go(names, {x for go in gos for x in ancestors(go, parents)})
        if labs: gene_labels[g['gene_id']] = labs
    ops = list(csv.DictReader(open(a.data_root/'regulation/ecoli_k12/operons.tsv'), delimiter='\t'))
    op_labels = []
    for op in ops:
        labs = set()
        for gid in (op.get('gene_ids') or '').split(','):
            labs |= gene_labels.get(gid, set())
        op_labels.append({'operon_id': op['operon_id'], 'label': op['label'], 'labels': sorted(labs)})
    summary['ecoli_k12'] = {'n_genes': len(genes), 'n_gaf_ids': len(gaf),
                            'n_gene_labels': len(gene_labels),
                            'n_operons': len(ops),
                            'n_operons_with_label': sum(1 for x in op_labels if x['labels']),
                            'label_counts': dict(Counter(l for x in op_labels for l in x['labels']))}
    # BSU: hits cds -> uniprot; operon cds
    hits = cases['subtilis'][2] if False else json.loads(cases['subtilis'][0].read_text()).get('hits', {})
    gaf_b = cases['subtilis'][2]
    gene_labels_b = {}
    for cid, h in hits.items():
        acc = (h.get('uniprot') or '').strip()
        gos = gaf_b.get(acc, set())
        labs = label_go(names, {x for go in gos for x in ancestors(go, parents)})
        if labs: gene_labels_b[cid] = labs
    ops_b = list(csv.DictReader(open(a.data_root.parent/'results/subtilis/subtilis.operons_auto.tsv'), delimiter='\t'))
    op_labels_b = []
    for op in ops_b:
        labs = set()
        for cid in (op.get('genes') or '').split(';'):
            labs |= gene_labels_b.get(cid, set())
        op_labels_b.append({'operon_id': op['operon_id'], 'label': op.get('label',''), 'labels': sorted(labs)})
    summary['subtilis'] = {'n_genes': len(hits), 'n_gaf_ids': len(gaf_b),
                           'n_gene_labels': len(gene_labels_b),
                           'n_operons': len(ops_b),
                           'n_operons_with_label': sum(1 for x in op_labels_b if x['labels']),
                           'label_counts': dict(Counter(l for x in op_labels_b for l in x['labels']))}
    for g, opsx in (('ecoli_k12', op_labels), ('subtilis', op_labels_b)):
        with (a.out/f'{g}_go_operon_labels.tsv').open('w', newline='') as fh:
            w=csv.DictWriter(fh, fieldnames=['operon_id','label','labels'], delimiter='\t'); w.writeheader()
            for x in opsx: w.writerow({**x,'labels':','.join(x['labels'])})
    (a.out/'go_labels_summary.json').write_text(json.dumps(summary, indent=1), encoding='utf-8')
    print(json.dumps(summary, indent=1))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
