#!/usr/bin/env python3
"""Cruza KEGG gene→pathway a operones de K-12/B. subtilis y valida K-12."""
from __future__ import annotations
import argparse, csv, json, re, sys
from collections import Counter, defaultdict
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from Bio import SeqIO


def read(p: Path):
    if not p.exists(): return []
    with p.open(encoding='utf-8', errors='replace') as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def kegg_links(path: Path):
    out = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if "\t" not in line:
            continue
        a, b = line.split("\t", 1)
        gene = a.split(":", 1)[-1].strip()
        target = b.split(":", 1)[-1].strip()
        if gene and target:
            out.append((gene, target))
    return out


def bsu_locus_to_cds(data: Path) -> dict[str, str]:
    gbff = data/'dbtbs/subtilis.gbff'; gff = data.parent/'results/subtilis/subtilis.gff3'
    coord = {}
    for rec in SeqIO.parse(str(gbff), 'genbank'):
        for f in rec.features:
            if f.type != 'CDS': continue
            lt = (f.qualifiers.get('locus_tag') or [''])[0]
            if lt:
                coord[(int(f.location.start)+1, int(f.location.end), '+' if f.location.strand == 1 else '-')] = lt
    out = {}
    for line in gff.read_text().splitlines():
        if line.startswith('#'): continue
        p = line.split('\t')
        if len(p) != 9 or p[2] != 'CDS': continue
        m = re.search(r'ID=([^;]+)', p[8])
        if not m: continue
        lt = coord.get((int(p[3]), int(p[4]), p[6]))
        if lt:
            out[lt] = m.group(1)
            out[lt.replace("_", "")] = m.group(1)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--data-root', type=Path, default=ROOT/'data')
    ap.add_argument('--out', type=Path, default=ROOT/'reports/labels_pilot')
    a = ap.parse_args(); a.out.mkdir(parents=True, exist_ok=True)
    ext = a.data_root/'labels/external/kegg'
    summary = {}
    # K12
    genes = read(a.data_root/'regulation/ecoli_k12/genes.tsv')
    locus_to_gid = {g['locus_tag']: g['gene_id'] for g in genes if g.get('locus_tag')}
    eco_map = {}
    for gene, path in kegg_links(ext/'eco/pathway.tsv'):
        gid2 = locus_to_gid.get(gene, gene)
        eco_map.setdefault(gid2, set()).add(path)
    names = {r[0]: r[1] for r in (line.split('\t',1) for line in (ext/'eco/pathway_names.tsv').read_text().splitlines() if '\t' in line)}
    ops = read(a.data_root/'regulation/ecoli_k12/operons.tsv')
    rows = []
    for op in ops:
        ps = sorted({p for gid in (op.get('gene_ids') or '').split(',') for p in eco_map.get(gid, ())})
        rows.append({'operon_id': op['operon_id'], 'label': op['label'], 'pathways': ps})
    summary['ecoli_k12'] = {'n_operons': len(rows), 'n_with_pathway': sum(1 for r in rows if r['pathways']),
                            'n_pathways': len({p for r in rows for p in r['pathways']})}
    # BSU
    hits = json.loads((a.data_root.parent/'results/subtilis/subtilis.functional_annotation.json').read_text()).get('hits', {})
    lt2cds = bsu_locus_to_cds(a.data_root)
    bsu_map = {}
    for gene, path in kegg_links(ext/'bsu/pathway.tsv'):
        cds = lt2cds.get(gene)
        if cds: bsu_map.setdefault(cds, set()).add(path)
    opsb = read(a.data_root.parent/'results/subtilis/subtilis.operons_auto.tsv')
    rowsb = []
    for op in opsb:
        ps = sorted({p for cid in (op.get('genes') or '').split(';') for p in bsu_map.get(cid, ())})
        rowsb.append({'operon_id': op['operon_id'], 'label': op.get('label',''), 'pathways': ps})
    summary['subtilis'] = {'n_operons': len(rowsb), 'n_with_pathway': sum(1 for r in rowsb if r['pathways']),
                           'n_pathways': len({p for r in rowsb for p in r['pathways']})}
    for g, rr in (('ecoli_k12', rows), ('subtilis', rowsb)):
        with (a.out/f'{g}_kegg_pathway_operons.tsv').open('w', newline='') as fh:
            w=csv.DictWriter(fh, fieldnames=['operon_id','label','pathways'], delimiter='\t'); w.writeheader()
            for x in rr: w.writerow({**x,'pathways':';'.join(x['pathways'])})
    # K12 validation
    ops_by_label = {r['label']: r['operon_id'] for r in ops}
    op_path = {r['operon_id']: r['pathways'] for r in rows}
    net = read(a.data_root/'regulondb/NetworkSigmaTU.tsv')
    targets = defaultdict(set)
    for r in net:
        sig = (r.get('1)sigmaName') or '').strip()
        tu = re.split(r'[\[\s]', (r.get('2)regulatedTUName') or '').strip())[0]
        if sig and tu in ops_by_label: targets[sig].add(ops_by_label[tu])
    allp = Counter(p for ps in op_path.values() for p in ps)
    validation={}
    for sig,tids in sorted(targets.items()):
        cnt=Counter(p for oid in tids for p in op_path.get(oid,()))
        top=[]
        for p,n in cnt.most_common(8):
            bg=allp.get(p,0)-n
            if bg < 1: pass
            try:
                from scipy.stats import fisher_exact
                odds,pv=fisher_exact([[n,len(tids)-n],[bg,len(op_path)-len(tids)-bg]],alternative='greater')
            except Exception: odds,pv=None,None
            top.append({'pathway':names.get(p,p),'n':n,'odds':round(float(odds),2) if odds else None,'p':round(float(pv),5) if pv else None})
        validation[sig]={'n_targets_mapped':len(tids),'top':top}
    (a.out/'k12_kegg_validation.json').write_text(json.dumps(validation,indent=1),encoding='utf-8')
    (a.out/'kegg_summary.json').write_text(json.dumps(summary,indent=1),encoding='utf-8')
    print(json.dumps(summary,indent=1)); print(json.dumps(validation,indent=1)[:3000])
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
