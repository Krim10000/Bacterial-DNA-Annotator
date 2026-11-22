#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""09_make_gene_map.py — genera {genome}_gene_map.tsv a partir de anotación.

El 07 (portal) y el 08 (plegado) consumen un TSV con la procedencia génica de
cada reacción Rhea:

    rhea_id  gene_id  product  locus_tag  contig  start  end  strand
    gene_length_bp  hypothetical  [sequence]

Este script lo construye desde el resultado de una anotación estándar,
uniendo los números EC de los genes con el mapeo oficial EC->Rhea
(rhea2ec.tsv, mismo formato rhea2uniprot.tsv). Fuentes aceptadas:

  * --gff      GFF3 de NCBI/PGAP, prokka o bakta (.gff/.gff3, con o sin .gz):
               coords de las líneas CDS; EC extraído del product (regex
               "(EC x.y.z.w)") Y de los atributos ec=, ec_number=, Dbxref=EC:
               y note= (bakta escribe ec= en el atributo 9; PGAP lo lleva en
               Dbxref).
  * --gbk      GenBank/GBFF de NCBI (.gbk/.gbff, con o sin .gz), uno o varios
               registros (un contig por LOCUS): /EC_number, /db_xref="EC:...",
               /product y /translation (secuencia sin FAA).
  * --faa      FASTA de proteínas (.faa[.gz]): da el gene_id real y añade la
               columna `sequence` (el 08 la usa si no hay --proteins-dir).
  * --bakta    TSV de bakta (.bakta.tsv / functional.tsv, con o sin .gz).
  * --eggnog   .emapper.annotations de eggNOG-mapper (columna EC); NO trae
               coords: enriquece los genes de la fuente principal (join por
               id) sin duplicar filas.
  * --net-dir  opcional: filtra los rhea_id a los que existen en
               {genome}_network.json (los gene-maps quedan ajustados a la red).

Modo por defecto (auto-descubrimiento): con solo --genome busca en
--annotations-dir (data/annotations, data/annotations_ab_initio,
data/proteomes) los ficheros del genoma, incluido el layout por-genoma de
bakta (annotations_ab_initio/{genome}/annotation.gff3 + proteins.faa +
functional.tsv + annotation.gbk). Fuente de genes: gbk > bakta > gff.

Modo LOTE (sin argumentos): recorre los genomas con red en --net-dir
(por defecto data/metabolism) y procesa solo los PENDIENTES: los que aún
no tienen gene map y los que lo tienen vacío; los que ya tienen filas se
saltan (usa --force para rehacerlos también). --dry-run muestra el plan
sin generar nada y --limit N sirve para probar con unos pocos. Los fallos
de cada genoma no tumbar el lote y se listan en el resumen final.

Salida: {out_dir}/{genome}_gene_map.tsv ordenado por (contig, start).

Ejemplos:
  python scripts/metabolism/09_make_gene_map.py --net-dir data/metabolism
  python scripts/metabolism/09_make_gene_map.py \
      --genome GCF_000006985.1_ASM698v1 --net-dir data/metabolism
  python scripts/metabolism/09_make_gene_map.py --genome GENOMA_X \
      --gbk data/annotations/GENOMA_X.gbff.gz
"""
import argparse
import glob
import gzip
import io
import re
import sys
import time
import urllib.request
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import importlib
m04 = importlib.import_module("04_network_trace")
log = m04.log

PROJECT_ROOT = m04.PROJECT_ROOT
RHEA_EC_URL = "https://ftp.expasy.org/databases/rhea/tsv/rhea2ec.tsv"
EC_RE = re.compile(r"EC[:\s]*([0-9]+\.[0-9-]+(?:\.[0-9-]+)?(?:\.[0-9n-]+)?)",
                   re.IGNORECASE)
HYP_RE = re.compile(r"hypothetical|unknown function|uncharacterized|"
                    r"función desconocida|funcion desconocida", re.IGNORECASE)
EC_ATTR_KEYS = ("ec", "ec_number", "dbxref", "note")
DEFAULT_ANN_DIRS = ("data/annotations", "data/annotations_ab_initio",
                    "data/proteomes")
# layout por-genoma (bakta / NCBI datasets): annotations_ab_initio/{genome}/
SUBDIR_FILES = {"annotation.gff3": "gff", "annotation.gff": "gff",
                "annotation.gbk": "gbk", "annotation.gbff": "gbk",
                "proteins.faa": "faa", "functional.tsv": "bakta"}


def open_text(path):
    """Abre un fichero como texto UTF-8; .gz si lo es (sufijo o magic
    bytes) — NCBI/prokka/bakta/eggnog reparten en .gff.gz, .faa.gz..."""
    raw = open(path, "rb")
    if hasattr(raw, "peek"):
        gz = raw.peek(2)[:2] == b"\x1f\x8b"
    else:
        gz = raw.read(2) == b"\x1f\x8b"
        raw.seek(0)
    if gz:
        raw.close()
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return io.TextIOWrapper(raw, encoding="utf-8", errors="replace")


def parse_attrs(field):
    out = {}
    for kv in field.rstrip(";").split(";"):
        if "=" in kv:
            k, _, v = kv.partition("=")
            out[k.strip()] = v.strip()
    return out


def ec_list(text):
    """Todos los EC de un texto (columna EC, atributo GFF, qualifier GBFF o
    product). Entiende "1.1.1.1", "EC:1.1.1.1" y listas con , o ;."""
    if not text or text in ("-", ""):
        return set()
    if re.match(r"^\s*[0-9.]+\s*$", text):      # columna EC limpia
        return {text.strip()}
    if ";" in text or "," in text and not re.search(r"EC", text):
        return {t.strip() for t in re.split(r"[;,]", text) if EC_RE.match(t.strip())}
    return set(EC_RE.findall(text))


def attr_ecs(attrs):
    """EC de los atributos GFF relevantes (bakta: ec=; PGAP: Dbxref=EC:...)."""
    out = set()
    for k, v in attrs.items():
        if k.lower() in EC_ATTR_KEYS:
            out |= ec_list(v)
    return out


_RHEA2EC_MEMO = {}


def read_rhea2ec(rhea_dir, force_download=False):
    """EC -> {rhea_int}. Descarga rhea2ec.tsv si falta (como rhea2uniprot).
    Memoizado por (path, mtime, size): el lote de 400 lo pediría 400 veces."""
    path = rhea_dir / "rhea2ec.tsv"
    if not path.exists() or force_download:
        log.info(f"Descargando {RHEA_EC_URL} -> {path}")
        try:
            urllib.request.urlretrieve(RHEA_EC_URL, path)
        except Exception as e:  # noqa: BLE001
            log.error(f"No se pudo descargar ({e}). Hazlo a mano:")
            log.error(f"    curl -o {path} {RHEA_EC_URL}")
            raise SystemExit(2)
    st = path.stat()
    hit = _RHEA2EC_MEMO.get((str(path), st.st_mtime_ns, st.st_size))
    if hit is not None:
        return hit
    ec_of = defaultdict(set)
    with open_text(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 4:
                rid = m04.normalize_rhea_id(parts[0])
                if rid and parts[3] and parts[3] != "-":
                    ec_of[parts[3].strip()].add(rid[1])
    log.info(f"rhea2ec: {len(ec_of)} EC con reacción Rhea")
    _RHEA2EC_MEMO[(str(path), st.st_mtime_ns, st.st_size)] = ec_of
    return ec_of


def net_rheas(net_dir, genome):
    """Conjunto de ints Rhea presentes en {genome}_network.json (o None)."""
    if not net_dir:
        return None
    p = Path(net_dir) / f"{genome}_network.json"
    if not p.exists():
        log.warning(f"{p.name} no existe: no filtro rhea_id por la red")
        return None
    import json
    with open(p, encoding="utf-8") as f:
        data = json.load(f)
    out = set()
    for n in data.get("nodes", []):
        if n.get("type") == "reaction":
            norm = m04.normalize_rhea_id(n.get("id", ""))
            if norm:
                out.add(norm[1])
    log.info(f"Filtro por la red: {len(out)} reacciones en {p.name}")
    return out


def emit(rows, gene, ec_of, net_set, stats):
    for ec in gene["ecs"]:
        for rid in ec_of.get(ec, ()):
            if net_set is not None and rid not in net_set:
                continue
            rows.append({"rhea_id": f"RHEA:{rid}", "gene_id": gene["id"],
                         "product": gene["product"],
                         "locus_tag": gene["locus"],
                         "contig": gene["contig"], "start": gene["start"],
                         "end": gene["end"], "strand": gene["strand"],
                         "gene_length_bp": gene["len"],
                         "hypothetical": "true" if gene["hyp"] else "false",
                         "sequence": gene["seq"]})
            stats["rows"] += 1
            stats["rh_mapped"].add(rid)


def load_faa(path):
    """{(primer token del header): secuencia}."""
    out, key, seq = {}, None, []
    with open_text(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if key:
                    out[key] = "".join(seq)
                key, seq = line[1:].split()[0], []
            elif line:
                seq.append(line)
    if key:
        out[key] = "".join(seq)
    return out


def parse_gff(gff_path, faa):
    """Genes CDS de un GFF3 (NCBI/PGAP, prokka, bakta)."""
    cds = defaultdict(list)
    order = {}
    with open_text(gff_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[2] != "CDS":
                continue
            a = parse_attrs(p[8])
            key = (p[0], a.get("ID") or a.get("protein_id") or
                   a.get("locus_tag") or f"{p[0]}:{p[3]}-{p[4]}")
            cds[key].append((int(p[3]), int(p[4])))
            order.setdefault(key, a)
    genes = []
    for key, segs in cds.items():
        a = order[key]
        start, end = min(s for s, _ in segs), max(e for _, e in segs)
        product = a.get("product", "")
        # gene_id: preferimos el id que existe en el FAA (prokka usa
        # locus_tag como header; NCBI usa protein_id WP_...)
        cands = [c for c in (a.get("protein_id"), a.get("ID"),
                             a.get("locus_tag")) if c]
        if faa:
            gid = next((c for c in cands if c in faa), "")
        else:
            gid = cands[0] if cands else ""
        genes.append({
            "id": gid or (cands[0] if cands else key[1]),
            "_key": key,
            "locus": a.get("locus_tag", a.get("gene", "")),
            "product": product,
            "contig": key[0], "start": start, "end": end,
            "strand": "",   # lo completa main con parse_gff_strands
            "len": end - start + 1,
            "hyp": bool(HYP_RE.search(product)),
            "ecs": ec_list(product) | attr_ecs(a),
            "seq": faa.get(gid, "") if gid else "",
        })
    return genes


def parse_gff_strands(gff_path):
    """segunda pasada barata: strand y gene_id por key (los GFF parten CDS)."""
    strands = {}
    with open_text(gff_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[2] != "CDS":
                continue
            a = parse_attrs(p[8])
            key = (p[0], a.get("ID") or a.get("protein_id") or
                   a.get("locus_tag") or f"{p[0]}:{p[3]}-{p[4]}")
            strands.setdefault(key, p[6])
    return strands


def parse_gbk(gbk_path, faa=None):
    """Genes CDS de un GenBank/GBFF (PGAP): /EC_number, /db_xref="EC:...",
    /product, /translation. Cada LOCUS = un contig (replicón)."""
    genes = []
    acc = ""
    in_features = False
    feat = None
    qual = None   # [key, [lineas]]

    def finish_qual():
        if feat is not None and qual is not None:
            val = "".join(qual[1]).strip().strip('"')
            feat["q"].setdefault(qual[0].lower(), []).append(val)

    def finish_feat():
        nonlocal feat, qual
        finish_qual()
        qual = None
        if feat is None or feat.get("type") != "CDS":
            feat = None
            return
        q = feat["q"]
        loc_str = "".join(feat["loc"])
        nums = [int(x) for x in re.findall(r"\d+", loc_str)]
        product = (q.get("product") or [""])[0]
        locus = (q.get("locus_tag") or [""])[0]
        prot = (q.get("protein_id") or [""])[0]
        cands = [c for c in (prot, locus) if c]
        if faa:
            gid = next((c for c in cands if c in faa), "")
        else:
            gid = cands[0] if cands else ""
        ecs = {t for k in ("ec_number",) for t in q.get(k, [])}
        for v in q.get("db_xref", []):
            if v.upper().startswith("EC:"):
                ecs.add(v[3:].strip())
        ecs |= ec_list(product)
        for n in q.get("note", []):
            ecs |= ec_list(n)
        seq = "".join(q.get("translation", [])).strip()
        if nums:
            genes.append({
                "id": gid or locus or f"{acc}:{nums[0]}",
                "locus": locus, "product": product, "contig": acc,
                "start": min(nums), "end": max(nums),
                "strand": "-" if "complement" in loc_str.lower() else "+",
                "len": max(nums) - min(nums) + 1,
                "hyp": bool(HYP_RE.search(product)),
                "ecs": ecs,
                "seq": faa.get(gid, "") if (gid and faa) else seq,
            })
        feat = None

    with open_text(gbk_path) as f:
        for line in f:
            if line.startswith("LOCUS"):
                finish_feat()
                in_features = False
                acc = line.split()[1] if len(line.split()) > 1 else ""
                continue
            if line.startswith("VERSION") and len(line.split()) > 1:
                acc = line.split()[1]      # accesión con versión (NC_002695.1)
                continue
            if line.startswith("FEATURES"):
                in_features = True
                continue
            if line.startswith(("ORIGIN", "//")):
                finish_feat()
                in_features = False
                continue
            if not in_features:
                continue
            s = line.rstrip("\n")
            stripped = s.strip()
            m5 = re.match(r"^ {5}([A-Za-z_][^\s]*)\s", s)
            if m5:
                finish_feat()
                feat = {"type": m5.group(1), "loc": [s[5:].strip()], "q": {}}
            elif stripped.startswith("/"):
                mq = re.match(r"^/([^=\s]+)(?:\s*=\s*(.*))?$", stripped)
                finish_qual()
                qual = [mq.group(1), [mq.group(2) or ""]] if mq else None
            elif feat is not None:
                if qual is not None:
                    qual[1].append(stripped)
                else:
                    feat["loc"].append(stripped)
    finish_feat()
    return genes


def parse_bakta(tsv_path, faa):
    genes = []
    with open_text(tsv_path) as f:
        header = None
        for line in f:
            if not line.strip():
                continue
            if header is None:               # bakta: "# ID  Query ID  ..."
                header = [h.strip().lower().lstrip("#").strip()
                          for h in line.rstrip("\n").split("\t")]
                continue
            parts = line.rstrip("\n").split("\t")
            rec = dict(zip(header, parts))
            try:
                start, end = int(rec["start"]), int(rec["end"])
            except (KeyError, ValueError):
                continue
            product = rec.get("product", "")
            gid = rec.get("id", "")
            genes.append({
                "id": gid, "locus": rec.get("locus", rec.get("gene", "")),
                "product": product, "contig": rec.get("contig", ""),
                "start": start, "end": end,
                "strand": rec.get("strand", ""), "len": end - start + 1,
                "hyp": bool(HYP_RE.search(product)),
                "ecs": ec_list(rec.get("ec", "")) | ec_list(product),
                "seq": faa.get(gid, ""),
            })
    return genes


def parse_eggnog(emapper_path, gene_source, faa=None):
    """EC de eggNOG. Con genes ya parseados (gff/bakta/gbk), ENRIQUECE esos
    genes (union de ECs, secuencia si falta) y devuelve solo los genes sin
    equivalencia — nunca duplica la procedencia génica."""
    genes = {}
    col = {"query": 0, "ec": None, "product": None}
    with open_text(emapper_path) as f:
        for line in f:
            if line.startswith("#"):
                if col["ec"] is not None:
                    continue
                fields = [t for t in line.lstrip("#").split() if t]
                try:
                    col["query"] = fields.index("query")
                    col["ec"] = fields.index("EC")
                    col["product"] = fields.index("Description")
                except ValueError:
                    continue
                continue
            parts = line.rstrip("\n").split("\t")
            if col["ec"] is None or len(parts) <= col["ec"]:
                continue
            q = parts[col["query"]]
            ecs = {t.strip() for t in re.split(r"[;,]", parts[col["ec"]])
                   if t.strip() and t.strip() != "-"}
            product = (parts[col["product"]] if col["product"] is not None
                       else "")
            genes[q] = {"id": q, "locus": "", "product": product,
                        "contig": "", "start": 0, "end": 0, "strand": "",
                        "len": 0,
                        "hyp": bool(HYP_RE.search(product)) or not product,
                        "ecs": ecs, "seq": (faa or {}).get(q, "")}
    if not gene_source:
        return list(genes.values())
    by_id = {}
    for g in gene_source:
        if g["id"]:
            by_id[g["id"]] = g
        if g["locus"]:
            by_id.setdefault(g["locus"], g)
    matched = set()
    for q, eg in genes.items():
        src = by_id.get(q)
        if src is None:
            continue
        matched.add(q)
        src["ecs"] |= eg["ecs"]
        if not src["seq"] and eg["seq"]:
            src["seq"] = eg["seq"]
    log.info(f"eggNOG: {len(matched)}/{len(genes)} genes enriquecidos "
             f"({len(genes) - len(matched)} sin coords quedan aparte)")
    return [g for q, g in genes.items() if q not in matched]


def classify_annotation(name):
    """Rol de un fichero de anotación por sufijo: gff|faa|bakta|eggnog|gbk."""
    stem = name.lower()
    if stem.endswith(".gz"):
        stem = stem[:-3]
    if stem.endswith((".gff", ".gff3")):
        return "gff"
    if stem.endswith(".faa"):
        return "faa"
    if stem.endswith((".bakta.tsv", ".functional.tsv")):
        return "bakta"
    if stem.endswith(".emapper.annotations"):
        return "eggnog"
    if stem.endswith((".gbk", ".gbff", ".genbank", ".gb")):
        return "gbk"
    return None


def discover_annotations(ann_dirs, genome):
    """Busca la anotación de {genome} en una lista de dirs.

    Reconoce ficheros {genome}*.<ext> (con o sin .gz, nombres NCBI o
    renombrados) y el layout por-genoma {dir}/{genome}/annotation.gff3 +
    proteins.faa + functional.tsv + annotation.gbk.
    Devuelve dict {gff, faa, bakta, eggnog, gbk} con Paths o None."""
    out = {"gff": None, "faa": None, "bakta": None, "eggnog": None,
           "gbk": None}
    for d in ann_dirs:
        d = Path(d)
        if not d.is_dir():
            continue
        for p in sorted(d.glob(glob.escape(genome) + "*")):
            if p.is_dir():
                for q in sorted(p.iterdir()):
                    role = (SUBDIR_FILES.get(q.name.lower()) or
                            classify_annotation(q.name))
                    if role:
                        out[role] = out[role] or q
            else:
                role = classify_annotation(p.name)
                if role:
                    out[role] = out[role] or p
        if (out["gff"] or out["bakta"] or out["eggnog"] or out["gbk"]) \
                and out["faa"]:
            break
    return out


def build_gene_map(genome, gff=None, faa=None, bakta=None, eggnog=None,
                   gbk=None, annotations_dir=None, rhea_dir=None,
                   net_dir=None, out_dir=None):
    """Embudo completo: descubre la anotación si falta, une EC->Rhea y
    escribe {out_dir}/{genome}_gene_map.tsv. Devuelve el Path (también con
    el TSV vacío); SystemExit si no hay ninguna anotación."""
    rhea_dir = Path(rhea_dir) if rhea_dir else m04.RHEA_DIR
    out_dir = Path(out_dir) if out_dir else \
        PROJECT_ROOT / "data" / "metabolism" / "genemaps"
    out_dir.mkdir(parents=True, exist_ok=True)

    if not (gff or bakta or eggnog or gbk):
        dirs = ([annotations_dir] if annotations_dir else
                [PROJECT_ROOT / d for d in DEFAULT_ANN_DIRS])
        found = discover_annotations(dirs, genome)
        if not (found["gff"] or found["bakta"] or found["eggnog"] or
                found["gbk"]):
            raise SystemExit(
                f"sin anotación para {genome} en {', '.join(str(d) for d in dirs)}: "
                f"esperaba {{genome}}*.gff[3][.gz] (+{{genome}}*.faa[.gz]), "
                f"{{genome}}*.gbff[.gz], {{genome}}*.bakta.tsv[.gz] o "
                f"{{genome}}*.emapper.annotations[.gz]")
        log.info(f"auto-descubierta: gff={found['gff'].name if found['gff'] else '-'} "
                 f"faa={found['faa'].name if found['faa'] else '-'} "
                 f"bakta={found['bakta'].name if found['bakta'] else '-'} "
                 f"eggnog={found['eggnog'].name if found['eggnog'] else '-'} "
                 f"gbk={found['gbk'].name if found['gbk'] else '-'}")
        faa = faa or found["faa"]
        bakta = bakta or found["bakta"]
        eggnog = eggnog or found["eggnog"]
        gff = gff or found["gff"]
        gbk = gbk or found["gbk"]

    # una sola fuente de genes: gbk > bakta > gff (duplicar genes engaña)
    if gbk and (gff or bakta):
        log.info(f"gbk + {'gff' if gff else 'bakta'}: uso solo el GenBank "
                 "(los otros duplicarían genes)")
        gff = bakta = None
    if bakta and gff and not eggnog:
        log.info("bakta + gff: uso solo el TSV de bakta (el GFF duplicaría genes)")
        gff = None

    faa_map = load_faa(faa) if faa else {}
    if faa:
        log.info(f"FAA: {len(faa_map)} proteínas")

    genes = []
    if gbk:
        genes += parse_gbk(gbk, faa_map)
        log.info(f"GenBank: {len(genes)} CDS (contigs "
                 f"{len({g['contig'] for g in genes})})")
    elif gff:
        strands = parse_gff_strands(gff)
        for g in parse_gff(gff, faa_map):
            g["strand"] = strands.get(g["_key"], "")
            del g["_key"]
            genes.append(g)
        log.info(f"GFF: {len(genes)} CDS")
    elif bakta:
        genes += parse_bakta(bakta, faa_map)
        log.info(f"bakta: {len(genes)} genes")
    if eggnog:
        genes += parse_eggnog(eggnog, genes or None, faa_map)

    n_ec = sum(1 for g in genes if g["ecs"])
    log.info(f"genes con EC anotado: {n_ec}/{len(genes)}")
    if not n_ec:
        log.warning("Ningún gene tiene EC: el gene_map saldrá vacío.")
        samples = [g["product"][:60] for g in genes[:3] if g["product"]]
        if samples:
            log.warning(f"  productos de ejemplo: {samples}")
        log.warning("  La anotación no aporta números EC (ni en el product "
                    "ni en atributos ec=/Dbxref/EC_number). Los GFF3 de "
                    "NCBI no llevan EC: descarga los GenBank con "
                    "scripts/metabolism/10_fetch_gbff.py (llevan "
                    "/EC_number y /translation), o corrígelo con "
                    "eggNOG-mapper (--eggnog) o bakta.")

    ec_of = read_rhea2ec(rhea_dir)
    net_set = net_rheas(net_dir, genome)

    rows, stats = [], {"rows": 0, "rh_mapped": set()}
    for g in genes:
        emit(rows, g, ec_of, net_set, stats)
    rows.sort(key=lambda r: (r["contig"], r["start"], r["rhea_id"]))

    out = out_dir / f"{genome}_gene_map.tsv"
    cols = ["rhea_id", "gene_id", "product", "locus_tag", "contig", "start",
            "end", "strand", "gene_length_bp", "hypothetical"]
    if faa or gbk or any(g["seq"] for g in genes):
        cols.append("sequence")
    with open(out, "w", encoding="utf-8") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r[c]).replace("\t", " ") for c in cols)
                    + "\n")
    log.info(f"{out.name}: {stats['rows']} filas "
             f"({len(stats['rh_mapped'])} reacciones Rhea con gen) -> {out}")
    if not stats["rows"]:
        log.warning("gene_map vacío: el 07 mostrará 'sin procedencia génica' "
                    "y el 08 no tendrá hipotéticas que plegar")
    return out


def _gm_has_rows(path):
    """True si el gene map tiene al menos una fila de datos (los errores
    se propagan: un fichero ilegible cuenta como 'regenerar')."""
    if path.suffix == ".json":
        import json as _json
        return bool(_json.loads(path.read_text(encoding="utf-8")))
    with open(path, encoding="utf-8") as f:
        lines = [ln for ln in f if not ln.startswith("#")]
    return any(ln.strip() for ln in lines[1:])   # lines[0] = cabecera


def main():
    ap = argparse.ArgumentParser(
        description="Genera {genome}_gene_map.tsv desde anotación "
                    "(gbff/gff/faa, bakta o eggNOG) uniendola con "
                    "rhea2ec.tsv. SIN ARGUMENTOS procesa en lote los "
                    "pendientes de --net-dir (nuevos y vacíos).")
    ap.add_argument("--genome", default=None,
                    help="un genoma concreto; sin esta flag se procesan en "
                         "lote los pendientes de --net-dir")
    ap.add_argument("--gff", type=Path, default=None)
    ap.add_argument("--faa", type=Path, default=None)
    ap.add_argument("--bakta", type=Path, default=None)
    ap.add_argument("--eggnog", type=Path, default=None)
    ap.add_argument("--gbk", type=Path, default=None,
                    help="GenBank/GBFF (.gbk/.gbff[.gz]): /EC_number y "
                         "/translation; admite varios registros")
    ap.add_argument("--annotations-dir", type=Path, default=None,
                    help="dir de anotación para el auto-descubrimiento "
                         "(por defecto: data/annotations, "
                         "data/annotations_ab_initio y data/proteomes)")
    ap.add_argument("--rhea-dir", type=Path, default=None)
    ap.add_argument("--net-dir", type=Path, default=None,
                    help="dir de las redes ({genome}_network.json); filtra "
                         "los rhea_id por la red de cada genoma y en modo "
                         "lote define la lista de genomas (por defecto "
                         "data/metabolism)")
    ap.add_argument("--out-dir", type=Path, default=None,
                    help="destino de los gene maps (por defecto "
                         "data/metabolism/genemaps)")
    ap.add_argument("--force", action="store_true",
                    help="modo lote: rehacer TODOS los gene maps, incluso "
                         "los que ya tienen filas")
    ap.add_argument("--limit", type=int, default=None,
                    help="procesa como mucho N genomas (prueba a mano)")
    ap.add_argument("--dry-run", action="store_true",
                    help="muestra el plan (nuevos/vacíos/a rehacer) sin "
                         "generar nada")
    ap.add_argument("--regen-empty", action="store_true",
                    help="solo rehace los gene maps ya EXISTENTES que estén "
                         "vacíos (0 filas); el modo por defecto sin "
                         "--genome ya los incluye junto a los nuevos")
    args = ap.parse_args()

    out_dir = Path(args.out_dir) if args.out_dir else \
        PROJECT_ROOT / "data" / "metabolism" / "genemaps"
    net_dir = Path(args.net_dir) if args.net_dir else \
        PROJECT_ROOT / "data" / "metabolism"

    if args.regen_empty:
        empties = []
        for p in sorted(out_dir.glob("*_gene_map.*")):
            if p.suffix not in (".tsv", ".json"):
                continue
            try:
                if _gm_has_rows(p):
                    continue
            except Exception as e:  # noqa: BLE001 — ilegible: regenerar
                log.warning(f"{p.name}: ilegible ({e}), regenerando")
            empties.append(p.name.split("_gene_map.")[0])
        if not empties:
            log.info(f"sin gene maps vacíos en {out_dir}: nada que hacer")
            return
        log.info(f"regen-empty: {len(empties)} gene maps vacíos -> "
                 f"{out_dir}")
        ok = fail = 0
        for i, g in enumerate(empties, 1):
            log.info(f"[{i}/{len(empties)}] {g}")
            try:
                build_gene_map(g, gff=args.gff, faa=args.faa,
                               bakta=args.bakta, eggnog=args.eggnog,
                               gbk=args.gbk,
                               annotations_dir=args.annotations_dir,
                               rhea_dir=args.rhea_dir, net_dir=net_dir,
                               out_dir=out_dir)
                ok += 1
            except (Exception, SystemExit) as e:  # noqa: BLE001 — el 09 sale
                # con SystemExit si no hay anotación: no tumbar el lote
                log.error(f"  {g}: {e}")
                fail += 1
        log.info(f"regen-empty listo: {ok} regenerados, {fail} fallos")
        if fail:
            log.warning("los fallidos siguen vacíos: revisa el log de cada "
                        "uno (¿falta el GBFF de ese genoma?)")
        return

    if args.genome:
        build_gene_map(args.genome, gff=args.gff, faa=args.faa,
                       bakta=args.bakta, eggnog=args.eggnog, gbk=args.gbk,
                       annotations_dir=args.annotations_dir,
                       rhea_dir=args.rhea_dir, net_dir=args.net_dir,
                       out_dir=out_dir)
        return

    # ---- modo lote por defecto: los PENDIENTES de --net-dir
    nets = sorted(net_dir.glob("*_network.json")) if net_dir.is_dir() else []
    if not nets:
        ap.error(f"sin redes en {net_dir} (no hay *_network.json): pasa "
                 "--net-dir o genera las redes con el 04")
    genomes = [p.name[:-len("_network.json")] for p in nets]
    nuevos, vacios, rehacer, con_filas = [], [], [], []
    for g in genomes:
        p = out_dir / f"{g}_gene_map.tsv"
        if not p.exists():
            p = out_dir / f"{g}_gene_map.json"
        if not p.exists():
            nuevos.append(g)
            continue
        if args.force:
            rehacer.append(g)
            continue
        try:
            (con_filas if _gm_has_rows(p) else vacios).append(g)
        except Exception as e:  # noqa: BLE001 — ilegible: rehacer
            log.warning(f"{p.name}: ilegible ({e}), lo rehace")
            vacios.append(g)
    jobs = nuevos + vacios + rehacer
    log.info(f"modo lote: {len(genomes)} genomas con red en {net_dir}")
    log.info(f"  gene maps -> {out_dir}")
    log.info(f"  nuevos: {len(nuevos)}  vacíos: {len(vacios)}  "
             f"a rehacer (--force): {len(rehacer)}  "
             f"con filas (se saltan): {len(con_filas)}")
    if not jobs:
        log.info("nada pendiente: todos los gene maps tienen filas")
        return
    if args.dry_run:
        for g in jobs:
            print(f"  pendiente: {g}")
        print(f"DRY-RUN: {len(jobs)} gene maps por generar; nada escrito")
        return
    if args.limit:
        jobs = jobs[:args.limit]
        log.info(f"--limit: solo los primeros {len(jobs)}")

    t0 = time.perf_counter()
    ok = fail = 0
    fallidos = []
    for i, g in enumerate(jobs, 1):
        log.info(f"[{i}/{len(jobs)}] {g}")
        try:
            build_gene_map(g, gff=args.gff, faa=args.faa, bakta=args.bakta,
                           eggnog=args.eggnog, gbk=args.gbk,
                           annotations_dir=args.annotations_dir,
                           rhea_dir=args.rhea_dir, net_dir=net_dir,
                           out_dir=out_dir)
            ok += 1
        except (Exception, SystemExit) as e:  # noqa: BLE001 — seguir lote
            log.error(f"  {g}: {e}")
            fail += 1
            fallidos.append(g)
    dt = time.perf_counter() - t0
    log.info(f"lote listo en {dt/60:.1f} min: {ok} generados, {fail} "
             f"fallos, {len(con_filas)} se saltaron (ya tenían filas)")
    if fallidos:
        extra = (f" … y {len(fallidos) - 10} más"
                 if len(fallidos) > 10 else "")
        log.warning("fallidos (sin gene map usable): " +
                    ", ".join(fallidos[:10]) + extra)
        log.warning("re-lanza el mismo comando para reintentar (es "
                    "idempotente); --force para rehacer también los que "
                    "ya tienen filas")


if __name__ == "__main__":
    main()
