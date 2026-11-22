#!/usr/bin/env python3
"""05_classify.py — tipado de elementos móviles, defensa, CRISPR/Cas, AMR y virulencia.

Combina tres evidencias independientes y solo emite una fila cuando hay al menos
una evidencia específica (nombre de producto o dominio Pfam):

  1. **Nombre del producto** del CDS (producto NCBI mapeado al CDS ab initio, o
     descripción Swiss-Prot del pipeline ab initio).
  2. **Dominio Pfam** del panel curado (`_work/pfam_panel/panel.hmm`), con marca
     `passes_ga` para separar hits fiables de hits débiles.
  3. **Detección previa** (`data/mobile_elements.json`, `defense_systems.json`,
     `defense_virulence.json`) ya resuelta a coordenadas en `candidates.tsv`.

Las llamadas genéricas (p. ej. "methyltransferase", "ABC transporter",
"transposase" sin familia) se degradan o se descartan: la lista previa de AMR es
en gran parte transportadores no específicos y se documenta como descartada.

Además agrupa:
  - Profagos: CDS con marcadores de fago en regiones ≤10 kb (intact/partial/remnant).
  - Sistemas RM: endonucleasa + metilasa a ≤5 kb.
  - Sistemas TA: toxina + antitoxina a ≤3 kb.
  - Arrays CRISPR (de `04_crispr_arrays.py`) con contexto `cas` a ≤10 kb.

Salida: data/mobile_typing/_work/<genoma>/elements.tsv
        data/mobile_typing/_work/<genoma>/demoted.tsv
"""
from __future__ import annotations

import sys

import csv
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import paths as P  # noqa: E402

OUT = P.BASE
WORK = P.WORK

# ---------------------------------------------------------------- patrones ---
IS_RE = re.compile(r"\bIS(L?\d{1,4}[A-Z]?|NCY|200/IS605)\b", re.I)
CAS_RE = re.compile(
    r"\b(cas(?:1|2|3|4|5|6|7|8|9|10|11|12|13)[a-z]?|cse[12]|csm[2-6]|cmr[1-6]|"
    r"csy[1-4]|csa3|csn2|csf[1-5]|c2c[1-9])\b", re.I
)
PHAGE_NEG = re.compile(r"phage[- ]shock|\bpsp[a-f]?\b|phage shock", re.I)

AMR_RULES = [
    ("beta-lactamase", 2, r"beta[- ]?lactamase|penicillinase|cephalosporinase|carbapenemase|"
        r"\bndm\b|\bkpc\b|\bctx-m\b|\bshv\b|\btem-\d+|\bampc\b|\bmecA\b"),
    ("aminoglycoside", 2, r"aminoglycoside|\baac\(?\d|\baph\(?\d|\bant\(?\d|\badda\b|\bstr[ab]\b|"
        r"streptomycin|kanamycin|spectinomycin|gentamicin|tobramycin|neomycin"),
    ("tetracycline", 2, r"tetracycline|\btet(?:\([a-z0-9]{1,3}\)|[a-z])\b|\btetm\b|\bteto\b|\btetq\b|"
        r"\btetw\b|\btetl\b|\btetk\b|\btetb\b"),
    ("chloramphenicol", 2, r"chloramphenicol|\bcat\b|\bcmla\b|\bflor\b"),
    ("macrolide", 2, r"macrolide|\berm\b|\berm\(|\bmph\(|\bmef\(|\blnuc\b|\bmsr\b"),
    ("sulfonamide", 2, r"sulfonamide|\bsul[12]\b|sulfamethoxazole"),
    ("trimethoprim", 2, r"trimethoprim|\bdfr[ab]?\d*\b"),
    ("vancomycin", 2, r"vancomycin resistance|vancomycin|"
        r"(?:resistance|d-alanine|d-ala|peptidoglycan).{0,60}\bvan[abhrsxyz]\b|"
        r"\bvan[abhrsxyz]\b.{0,60}(?:resistance|d-alanine|d-ala|peptidoglycan)"),
    ("colistin", 2, r"colistin|\bmcr-?\d+|\bpmra\b|\bpmrb\b|\barna\b|\barnb\b|\barnc\b"),
    ("fosfomycin", 2, r"fosfomycin|\bfosa\b|\bfosb\b|\bfosx\b"),
    ("efflux", 2, r"multidrug|\bacr[abdef]\b|\bmex[abdf]\b|\bade[abcg]\b|\bcme[ab]\b|"
        r"\bmtr[cdé]\b|\bmdfa\b|\bemre\b|\bqaca?\b|acriflavin|chlorhexidine|"
        r"fluoroquinolone"),
    ("target_mutation_unverified", 1, r"(?:dna gyrase|topoisomerase|\brpob\b).{0,80}\b(?:gyra|gyrb|parc|pare|rpob)\b|"
        r"\b(?:gyra|gyrb|parc|pare|rpob)\b.{0,80}(?:dna gyrase|topoisomerase|rna polymerase|subunit)"),
]

VIR_RULES = [
    ("hemolysin_RTX", 2, r"hemolysin(?! iii)|haemolysin|\brtx\b|\bhly[abcd]?\b|cytolysin|leukotoxin"),
    ("toxin_CDT", 2, r"cytolethal|\bcdt[abc]\b|pertussis|cholera toxin|\bctx\b|shiga|verotoxin|"
        r"\bstx[12]\b|botulinum|tetanus|diphtheria"),
    ("toxin_colicin", 2, r"colicin|bacteriocin|microcin|pyocin"),
    ("adhesin", 2, r"adhesin|agglutinin|intimin|\beae\b|\bfim[a-z]?\b|\bpap[a-z]\b|\bcsga\b"),
    ("fimbriae", 2, r"fimbri|pili\b|pilin"),
    ("siderophore", 2, r"siderophore|enterobactin|aerobactin|yersiniabactin|pyochelin|pyoverdine|"
        r"\bfepa\b|\bfhua\b|\bfyua\b"),
    ("T3SS", 2, r"type iii secretion|\bt3ss\b|\bssa[abcfg]\b|\bhrp\b|\besca?\b|\bsepa?\b|"
        r"\binv[a-z]\b|\bipaa?[bcd]?\b|\byop[a-z]?\b"),
    ("T6SS", 2, r"type vi secretion|\bt6ss\b|\bhcp\b|\bvgrg\b|\bclpv\b|\bicmf\b"),
    ("capsule", 2, r"capsul|\bcps[a-z]?\b"),
    ("protease", 2, r"iga protease|igap|zinc protease"),
    ("phospholipase", 2, r"phospholipase c|lecithinase|sphingomyelinase"),
    ("virulence_other", 1, r"virulence"),
]

TA_FAMILIES = [
    ("RelE", r"\brele\b|\brelb\b"),
    ("MazF", r"\bmazf\b|\bmaze\b"),
    ("HicAB", r"\bhic[ab]\b"),
    ("VapBC", r"\bvap[bc]\b"),
    ("CcdAB", r"\bccd[ab]\b"),
    ("ParDE", r"\bpar[de]\b"),
    ("HigBA", r"\bhig[ab]\b"),
    ("YoeB", r"\byoeb\b|\byefm\b"),
    ("YafQ", r"\byafq\b|\bdinj\b"),
    ("MqsR", r"\bmqs[ra]\b"),
    ("DocPhd", r"\bdoc\b|\bphd\b"),
    ("CbtA", r"\bcbta\b"),
    ("BrnT", r"\bbrnt\b"),
]

PHAGE_MARKERS = [
    ("terminase", r"terminase"),
    ("portal", r"portal protein|phage portal"),
    ("capsid", r"capsid|major head|head protein"),
    ("tail", r"phage tail|tail fiber|tail fibre|tail tube|tail sheath|baseplate|"
             r"tail protein|tape measure|tail assembly|minor tail"),
    ("holin", r"holin"),
    ("lysin", r"lysin|lysozyme|endolysin"),
    ("integrase", r"integrase|recombinase"),
    ("phage_other", r"\bphage\b|prophage|bacteriophage"),
]

PFAM_IS = {
    "DDE_Tnp_IS1": "IS1", "HTH_Tnp_IS1": "IS1", "DDE_Tnp_IS66": "IS66",
    "DDE_Tnp_IS66_C": "IS66", "DDE_Tnp_ISL3": "ISL3", "HTH_Tnp_ISL3": "ISL3",
    "DDE_Tnp_Tn3": "Tn3", "DDE_Tnp_IS1595": "IS1595", "HTH_Tnp_IS1595": "IS1595",
    "DDE_Tnp_IS240": "IS240", "HTH_Tnp_IS630": "IS630",
}


def conf(score: int) -> str:
    return "high" if score >= 3 else "medium" if score == 2 else "low"


def parse_attrs_evidence(parts: dict) -> str:
    bits = []
    for key in ("prior", "name", "pfam", "map", "note", "system", "markers",
                "n_genes", "repeats", "spacers", "interrupted", "cas", "context"):
        if parts.get(key):
            bits.append(f"{key}={parts[key]}")
    return ";".join(bits)


def name_cas(product: str):
    if re.search(r"\bcascade\b", product, re.I):
        return "Cascade", "cascade", 2
    m = CAS_RE.search(product)
    if not m:
        if re.search(r"crispr[- ]associated", product, re.I):
            return "Cas_uncharacterized", "", 1
        return None
    raw = m.group(1).lower()
    if raw.startswith("cas"):
        subtype = raw.replace("cas", "Cas")
        subtype = re.sub(r"([0-9]+)([a-z])$", lambda x: x.group(1) + x.group(2).upper(), subtype)
        spec = 2
    elif raw in ("cse1", "cse2", "csm2", "csm3", "csm4", "csm5", "csm6", "cmr1",
                 "cmr2", "cmr3", "cmr4", "cmr5", "cmr6", "csy1", "csy2", "csy3",
                 "csy4", "csa3", "csn2", "csf1", "csf2", "csf3", "csf4", "csf5"):
        subtype = "Cascade"
        spec = 2
    else:
        subtype = raw.upper()
        spec = 2
    return subtype, m.group(1), spec


def name_amr(product: str):
    if re.search(r"transcriptional regulator|repressor|regulatory protein|\bbla[ir]\b", product, re.I):
        return None
    for subtype, spec, pattern in AMR_RULES:
        m = re.search(pattern, product, re.I)
        if not m:
            continue
        matched = m.group(0)
        # TetR es un represor, no un determinante de resistencia
        if subtype == "tetracycline" and re.fullmatch(r"tetr", matched, re.I):
            continue
        # PBPs y carboxipeptidasas son dianas, no determinantes de resistencia
        if subtype == "beta-lactamase" and re.search(
                r"carboxypeptidase|penicillin[- ]binding|d-alanyl|\bpbp\b|peptidase|transpeptidase",
                product, re.I) and not re.search(r"beta[- ]?lactamase", product, re.I):
            continue
        # ParE/ParC: solo cuentan como diana con contexto de topoisomerasa
        if subtype == "target_mutation_unverified" and re.search(
                r"toxin|antitoxin|plasmid|partition", product, re.I):
            continue
        return subtype, matched, spec
    return None


def name_virulence(product: str):
    if re.search(r"toxin[- ]antitoxin|antitoxin|protease modulator|modulator", product, re.I):
        return None
    if re.search(r"hemolysin", product, re.I) and re.search(
            r"hlyd|shlb|secretion|activation|expression", product, re.I):
        return None
    for subtype, spec, pattern in VIR_RULES:
        m = re.search(pattern, product, re.I)
        if m:
            return subtype, m.group(0), spec
    return None


def ta_family(text: str) -> str:
    for family, pattern in TA_FAMILIES:
        if re.search(pattern, text, re.I):
            return family
    return "unclassified"


def name_defense(product: str):
    """Devuelve (role, subtype, family, specificity) solo con evidencia de nombre.

    Los dominios Pfam de defensa se usan únicamente como refuerzo: DefenseFinder
    y PADLOC no están instalados y los modelos Pfam de restricción/TA aparecen en
    proteínas housekeeping (aminopeptidasas, ACC, kinasas), así que una llamada
    basada solo en dominio sería un falso positivo.
    """
    low = product.lower()
    # metilasas y subunidades S/M de sistemas RM (posibles parejas)
    if re.search(r"methyltransferase|methylase|subunit m\b|hsdm", product, re.I) and re.search(
            r"restriction|modification|hsd|type i", product, re.I):
        subtype = "RM_type_I" if re.search(r"hsd|type i", low) else "RM_type_II"
        return "RM_methylase", subtype, "", 2
    if re.search(r"restriction endonuclease|\bhsdr\b|\bmrr\b|\bmcrbc\b|\bmcrb\b|\bmcrc\b|"
                 r"\bre_[a-z0-9]+\b|type i restriction|type ii restriction|"
                 r"type iii restriction|type iv restriction", product, re.I):
        if re.search(r"subunit s\b|specificity subunit", low):
            return "RM_specificity", "RM_context_only", "", 2
        subtype = "RM_type_I" if re.search(r"hsdr|type i restriction", low) else (
            "RM_type_IV" if re.search(r"mrr|mcr|\btype iv\b", low) else "RM_type_II")
        return "RM_endonuclease", subtype, "", 2
    # TA estricto: exige contexto de toxina-antitoxina o símbolo génico conocido
    if re.search(r"toxin[- ]antitoxin|\bantitoxin\b", product, re.I):
        role = "TA_antitoxin" if re.search(r"antitoxin", product, re.I) else "TA_toxin"
        family = ta_family(product)
        return role, f"TA_type_II_{family}", family, 2 if family != "unclassified" else 1
    if re.search(r"\brele\b|\bmazf\b|\bhicA\b|\bvapc\b|\bccdb\b|\bhigb\b|\byoeb\b|"
                 r"\byafq\b|\bmqsr\b|\bparE\b.*toxin|toxin.*\bparE\b", product, re.I):
        family = ta_family(product)
        return "TA_toxin", f"TA_type_II_{family}", family, 2 if family != "unclassified" else 1
    if re.search(r"\brelb\b|\bmaze\b|\bhicb\b|\bvapb\b|\bccda\b|\bhiga\b|\bphd\b|"
                 r"\bmqsa\b|\bparD\b", product, re.I):
        family = ta_family(product)
        return "TA_antitoxin", f"TA_type_II_{family}", family, 1
    if re.search(r"abortive infection|\babi[a-z]{1,2}\b", product, re.I):
        return "other", "Abi", "", 2
    if re.search(r"\bbrex\b|\bbrx[a-z]\b", product, re.I):
        return "other", "BREX", "", 2
    if re.search(r"\bcbass\b|cgas|cd-ntase|cd_ntase", product, re.I):
        return "other", "CBASS", "", 2
    if re.search(r"thoeris", product, re.I):
        return "other", "Thoeris", "", 2
    if re.search(r"\bdnd[a-e]\b|dna sulfur|phosphorothio", product, re.I):
        return "other", "Dnd", "", 2
    if re.search(r"anti-?phage|antiphage|defense system|defence system", product, re.I):
        return "other", "defense_uncharacterized", "", 1
    return None


def name_mobile(product: str):
    """Devuelve (subtype, gene_hint, specificity, family)."""
    is_m = IS_RE.search(product)
    if is_m:
        family = is_m.group(1).upper()
        return (family if family.startswith("IS") else f"IS{family}"), is_m.group(0), 2, "IS"
    if re.search(r"\btn\d+\b", product, re.I):
        m = re.search(r"\btn(\d+)\b", product, re.I)
        return f"Tn{m.group(1)}", m.group(0), 2, "Tn"
    if re.search(r"transposon|transposase|insertion sequence", product, re.I):
        return "transposase_unspecified", "", 2 if re.search(r"transposase", product, re.I) else 1, "transposase"
    if re.search(r"conjugat|relaxase|\bmoba\b|\btrag\b|\btrad\b|type iv secretion", product, re.I):
        return "conjugative_element", "", 2, "conjugative"
    if re.search(r"integrase|recombinase", product, re.I):
        if re.search(r"phage|prophage", product, re.I):
            return "phage_integrase", "", 2, "phage_integrase"
        if re.search(r"serine|resolvase|tn916|xer[cd]", product, re.I):
            return "integrase_serine", "", 2, "integrase"
        return "integrase_tyrosine", "", 1, "integrase"
    if re.search(r"resolvase", product, re.I):
        return "resolvase_serine", "", 2, "integrase"
    return None


def phage_markers(product: str):
    if PHAGE_NEG.search(product):
        return []
    out = []
    for category, pattern in PHAGE_MARKERS:
        if re.search(pattern, product, re.I):
            out.append(category)
    return out


def pfam_best(hits, category_set, exclude=()):
    best = None
    for hit in hits:
        cats = set(hit["categories"].split(","))
        if not (cats & category_set):
            continue
        if hit["pfam_name"] in exclude:
            continue
        spec = 2 if hit["passes_ga"] == "yes" else 1
        key = (spec, float(hit["score"]))
        if best is None or key > best[0]:
            best = (key, hit, spec)
    return best


def pfam_cas(hits):
    best = pfam_best(hits, {"cas"})
    if not best:
        return None
    hit = best[1]
    m = re.search(r"cas(\d+)", hit["pfam_name"], re.I)
    if m:
        subtype = f"Cas{m.group(1)}"
    elif re.search(r"cse|csm|cmr|csy|csa|csn|csf|c2c", hit["pfam_name"], re.I):
        subtype = "Cascade"
    else:
        subtype = "Cas_uncharacterized"
    return subtype, hit, best[2]


def pfam_amr(hits):
    best = pfam_best(hits, {"amr_betalactam", "amr_aminoglycoside",
                            "amr_chloramphenicol", "amr_efflux", "amr_tetracycline"},
                     exclude={"MFS_1", "MFS_1_like"})
    if not best:
        return None
    hit = best[1]
    name = hit["pfam_name"]
    cats = set(hit["categories"].split(","))
    if "Beta-lactamase" in name or "Beta-lactamase2" in name or "Lactamase_B" in name:
        subtype = "beta-lactamase"
    elif "CAT" in name:
        subtype = "chloramphenicol"
    elif re.search(r"aminogly|AadA|APH|ANT|Phosphotransferase", name):
        subtype = "aminoglycoside"
    elif re.search(r"MFS|Multi_Drug|EmrE|SepA|AcrZ|BTP|efflux|FlqE|Rv2686", name):
        subtype = "efflux"
    else:
        subtype = "amr_uncharacterized"
    return subtype, hit, best[2]


def pfam_defense(hits):
    best = pfam_best(hits, {"restriction", "toxin_antitoxin", "defense_other"})
    if not best:
        return None
    hit = best[1]
    name = hit["pfam_name"]
    if re.search(r"^(HsdR|Mrr|McrBC|RE_|Restriction)", name):
        role = "RM_endonuclease"
        subtype = "RM_type_IV" if re.search(r"Mrr|McrBC", name) else "RM_type_II"
    elif re.search(r"^(HsdM|N6_N4_Mtase)", name):
        role = "RM_methylase"
        subtype = "RM_context_only"
    elif re.search(r"^(RelE|PemK|HicA|HigB|ParE|YafQ|YoeB|MqsR|CbtA|BrnT|VapC)", name):
        role = "TA_toxin"
        subtype = f"TA_type_II_{ta_family(name)}"
    elif re.search(r"^(MazE|HicB|Phd|CcdA|CopG|HigA)", name):
        role = "TA_antitoxin"
        subtype = f"TA_type_II_{ta_family(name)}"
    elif re.search(r"^Abi|^Brx|^BREX|CD_NTase|^TIR|PglZ|^Dnd", name):
        role = "other"
        subtype = "defense_uncharacterized"
    else:
        role = "other"
        subtype = "defense_uncharacterized"
    return role, subtype, hit, best[2]


def pfam_mobile(hits):
    best = pfam_best(hits, {"transposase", "integrase"})
    if not best:
        return None
    hit = best[1]
    name = hit["pfam_name"]
    if name in PFAM_IS:
        subtype = PFAM_IS[name]
    elif re.match(r"DDE_Tnp", name):
        subtype = "transposase_unspecified"
    elif name in ("Resolvase", "Recombinase"):
        subtype = "resolvase_serine"
    elif name == "Phage_integrase":
        subtype = "phage_integrase"
    else:
        subtype = "integrase_tyrosine"
    # los HTH de transposasa son dominios DNA-binding ubicuos: solo apoyan
    spec = 1 if re.match(r"HTH_Tnp", name) else best[2]
    return subtype, hit, spec


def pfam_virulence(hits):
    best = pfam_best(hits, {"virulence"})
    if not best:
        return None
    hit = best[1]
    name = hit["pfam_name"]
    if re.search(r"HemolysinCabind|HlyC|HlyE|HlyIII|Hemolysin", name):
        subtype = "hemolysin_RTX"
    elif re.search(r"Fimbrial", name):
        subtype = "fimbriae"
    elif re.search(r"Adhesin", name):
        subtype = "adhesin"
    elif re.search(r"Colicin", name):
        subtype = "toxin_colicin"
    elif re.search(r"CDtoxinA|Enterotoxin|CagA|Esp[ABFG]|Yop", name):
        subtype = "toxin_CDT"
    else:
        subtype = "virulence_uncharacterized"
    return subtype, hit, best[2]


def load_pilots():
    return P.genome_list()


def load_cds(genome):
    with (WORK / genome / "cds.tsv").open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def load_priors(genome):
    with (WORK / genome / "candidates.tsv").open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def load_pfam(genome):
    path = WORK / genome / "pfam_hits.tsv"
    by_protein = defaultdict(list)
    if path.exists():
        with path.open() as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                by_protein[row["protein_id"]].append(row)
    return by_protein


def load_arrays(genome):
    path = WORK / genome / "crispr_arrays.tsv"
    if not path.exists():
        return []
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def classify(genome):
    cds = load_cds(genome)
    priors = load_priors(genome)
    pfam = load_pfam(genome)
    arrays = load_arrays(genome)

    array_by_contig = defaultdict(list)
    for arr in arrays:
        array_by_contig[arr["contig"]].append((int(arr["start"]), int(arr["end"])))

    def near_array(rec, max_distance=10_000):
        for start, end in array_by_contig.get(rec["contig"], []):
            dist = max(0, max(int(rec["start"]), start) - min(int(rec["end"]), end))
            if dist <= max_distance:
                return dist
        return None

    prior_by_cds = defaultdict(list)
    for p in priors:
        if p["cds_id"]:
            prior_by_cds[p["cds_id"]].append(p)

    # ---- 1. clasificación por CDS -----------------------------------------
    info = {}
    for rec in cds:
        product = rec["product"] or ""
        hits = pfam.get(rec["cds_id"], [])
        priors_here = prior_by_cds.get(rec["cds_id"], [])
        prior_str = ",".join(sorted({f"{p['source_json']}/{p['source_section']}" for p in priors_here}))
        info[rec["cds_id"]] = {
            "rec": rec,
            "product": product,
            "hits": hits,
            "priors": priors_here,
            "prior_str": prior_str,
            "cas": None, "amr": None, "defense": None, "mobile": None,
            "vir": None, "phage": phage_markers(product),
        }
        cas_n = name_cas(product)
        cas_p = None
        if not cas_n:
            # Sin nombre de cas: se acepta un dominio Cas central con GA o, con
            # al menos dos dominios Cas distintos, un gen a <=10 kb de un array.
            cas_hits = [h for h in hits if "cas" in h["categories"].split(",")]
            core = [h for h in cas_hits if h["passes_ga"] == "yes" and re.match(
                r"(Cas\d|Cmr\d|Csm\d|Cse\d|Csy\d|Csa\d|Csf\d|Csn\d|CRISPR)", h["pfam_name"])]
            if core or (len({h["pfam_name"] for h in cas_hits}) >= 2 and near_array(
                    {"contig": rec["contig"], "start": rec["start"], "end": rec["end"]}) is not None):
                cas_p = pfam_cas(hits)
        if cas_n:
            subtype = cas_n[0]
            gene = cas_n[1]
            score = cas_n[2] + (1 if cas_p and cas_p[2] == 2 else 0)
        elif cas_p:
            subtype = cas_p[0]
            gene = ""
            score = 2 if cas_p[2] == 2 else 1
        if cas_n or cas_p:
            info[rec["cds_id"]]["cas"] = {
                "subtype": subtype, "gene": gene, "score": score,
                "name": cas_n, "pfam": cas_p,
            }
        # AMR: se exige evidencia de NOMBRE específica (gen/familia). Los dominios
        # Pfam amplios (Lactamase_B, APH, Penicillinase_R, MFS_1…) dan falsos
        # positivos sobre hidrolasas y transportadores housekeeping, así que solo
        # refuerzan la confianza de una llamada ya hecha por nombre.
        amr_n = name_amr(product)
        amr_p = pfam_amr(hits) if amr_n else None
        if amr_n:
            score = amr_n[2] + (2 if amr_p and amr_p[2] == 2 else 0)
            info[rec["cds_id"]]["amr"] = {
                "subtype": amr_n[0], "score": score, "name": amr_n, "pfm": amr_p,
                "pfam": amr_p,
            }
        dfn_n = name_defense(product)
        dfn_p = pfam_defense(hits) if dfn_n else None
        if dfn_n:
            role, subtype, family, spec = dfn_n
            score = spec + (1 if dfn_p and dfn_p[3] == 2 else 0)
            info[rec["cds_id"]]["defense"] = {
                "role": role, "subtype": subtype, "family": family,
                "score": score, "name": dfn_n, "pfam": dfn_p,
                "pfam_hit": dfn_p[2] if dfn_p else None,
            }
        mob_n = name_mobile(product)
        mob_p = pfam_mobile(hits)
        if mob_n or (mob_p and mob_p[2] == 2):
            subtype = mob_n[0] if mob_n else mob_p[0]
            score = (mob_n[2] if mob_n else 0) + (mob_p[2] if mob_p else 0)
            if not mob_n:
                score = max(score, 2)  # dominio DDE/rve con GA, sin nombre
            if re.search(r"domain-containing|putative|fragment|partial|truncat", product, re.I):
                score -= 1
            info[rec["cds_id"]]["mobile"] = {
                "subtype": subtype, "score": score, "name": mob_n, "pfam": mob_p,
            }
        # Virulencia: igual que AMR, se exige nombre específico; el dominio Pfam
        # solo refuerza (evita llamadas sobre "hypothetical protein").
        vir_n = name_virulence(product)
        vir_p = pfam_virulence(hits) if vir_n else None
        if vir_n:
            score = vir_n[2] + (2 if vir_p and vir_p[2] == 2 else 0)
            info[rec["cds_id"]]["vir"] = {
                "subtype": vir_n[0], "score": score, "name": vir_n, "pfam": vir_p,
            }

    rows = []
    covered = set()
    rows_by_cds = defaultdict(list)

    def add_row(**kw):
        row = {
            "genome": genome,
            "contig": kw["contig"],
            "start": int(kw["start"]),
            "end": int(kw["end"]),
            "strand": kw.get("strand", "."),
            "element_type": kw["element_type"],
            "subtype": kw["subtype"],
            "gene": kw.get("gene", ""),
            "evidence": parse_attrs_evidence(kw.get("evidence", {})),
            "confidence": kw["confidence"],
            "source_cds": ";".join(kw.get("source_cds", [])),
            "rule": kw.get("rule", ""),
            "source": kw.get("source", "de_novo"),
        }
        rows.append(row)
        for cid in kw.get("source_cds", []):
            covered.add(cid)
            rows_by_cds[cid].append(row)
        return row

    def gene_of(rec):
        return rec["gene"] or rec["locus_tag"] or rec["cds_id"]

    # ---- 2. profagos: regiones --------------------------------------------
    markers = []
    for cid, meta in info.items():
        if meta["phage"]:
            rec = meta["rec"]
            # un integrasa suelto se trata como elemento móvil, no como fago
            if meta["phage"] == ["integrase"] and not meta["mobile"]:
                continue
            markers.append((rec, meta["phage"]))
    markers.sort(key=lambda x: (x[0]["contig"], int(x[0]["start"])))
    regions = []
    for rec, cats in markers:
        placed = False
        for region in regions:
            if region["contig"] != rec["contig"]:
                continue
            gap = int(rec["start"]) - region["end"] - 1
            if gap <= 10_000:
                region["end"] = max(region["end"], int(rec["end"]))
                region["cds"].append(rec)
                region["markers"].update(cats)
                placed = True
                break
        if not placed:
            regions.append({"contig": rec["contig"], "start": int(rec["start"]),
                            "end": int(rec["end"]), "cds": [rec], "markers": set(cats)})

    core = {"terminase", "portal", "capsid", "tail"}
    emitted_regions = []
    for region in regions:
        core_hits = region["markers"] & core
        if len(region["cds"]) < 3 or not core_hits:
            continue
        if core_hits == core:
            subtype, confianza = "intact", "high"
        elif len(core_hits) >= 2:
            subtype, confianza = "partial", "medium"
        else:
            subtype, confianza = "remnant", "low"
        strands = Counter(r["strand"] for r in region["cds"])
        strand = strands.most_common(1)[0][0]
        mixed = "yes" if len(strands) > 1 else ""
        prior_here = sorted({f"{p['source_json']}/{p['source_section']}"
                             for r in region["cds"]
                             for p in prior_by_cds.get(r["cds_id"], [])
                             if p["source_section"] == "prophage"})
        add_row(
            contig=region["contig"], start=region["start"], end=region["end"],
            strand=strand or ".", element_type="prophage", subtype=subtype,
            gene="region", confidence=confianza,
            evidence={
                "n_genes": len(region["cds"]),
                "markers": ",".join(f"{k}:{sum(1 for r in region['cds'] if k in phage_markers(r['product'] or ''))}"
                                    for k in sorted(region["markers"])),
                "prior": ":".join(prior_here) if prior_here else "",
                "note": "mixed_strand=yes" if mixed else "",
                "map": "abi_cds",
            },
            source_cds=[r["cds_id"] for r in region["cds"]],
            rule="prophage_region",
            source="prior+de_novo" if prior_here else "de_novo",
        )
        emitted_regions.append(region)

    region_cds = {r["cds_id"] for region in emitted_regions for r in region["cds"]}

    # ---- 3. sistemas RM y TA (pares) --------------------------------------
    rm_endos, rm_mts = [], []
    for cid, meta in info.items():
        d = meta["defense"]
        if d and d["role"] == "RM_endonuclease":
            rm_endos.append(cid)
        # pareja potencial: metilasa anotada como RM o metilasa genérica vecina
        if (d and d["role"] in ("RM_methylase", "RM_specificity")) or re.search(
                r"methyltransferase|methylase", meta["product"], re.I):
            rm_mts.append(cid)

    def nearby(cid, candidates, max_gap, prefer_methylase=False):
        rec = info[cid]["rec"]
        best = None
        for other in candidates:
            if other == cid:
                continue
            oth = info[other]["rec"]
            if oth["contig"] != rec["contig"] or oth["strand"] != rec["strand"]:
                continue
            gap = max(int(oth["start"]), int(rec["start"])) - min(int(oth["end"]), int(rec["end"]))
            if gap > max_gap:
                continue
            d = info[other]["defense"]
            is_specificity = bool(d and d["role"] == "RM_specificity")
            # una metilasa M real vale más que una subunidad S
            key = (1 if (prefer_methylase and is_specificity) else 0, gap)
            if best is None or key < best[0]:
                best = (key, gap, other)
        if best is None:
            return None
        return best[1], best[2]

    consumed_defense = set()
    for cid in rm_endos:
        rec = info[cid]["rec"]
        d = info[cid]["defense"]
        pair = nearby(cid, rm_mts, 5_000, prefer_methylase=True)
        if pair:
            gap, mt_id = pair
            mt = info[mt_id]["rec"]
            partner = info[mt_id]["defense"]
            specificity_only = bool(partner and partner["role"] == "RM_specificity")
            consumed_defense.update({cid, mt_id})
            subtype = d["subtype"] if d["subtype"].startswith("RM_type") else "RM_type_II"
            add_row(
                contig=rec["contig"], start=min(int(rec["start"]), int(mt["start"])),
                end=max(int(rec["end"]), int(mt["end"])), strand=rec["strand"],
                element_type="defense", subtype=subtype,
                gene=f"{gene_of(rec)}/{gene_of(mt)}",
                confidence="medium" if specificity_only else ("high" if d["score"] >= 3 else "medium"),
                evidence={"system": "RM", "name": f"RE={rec['product']}|MTase={mt['product']}",
                          "pfam": f"{d['pfam_hit']['pfam_name']}({d['pfam_hit']['pfam_acc']})" if d["pfam_hit"] else "",
                          "note": "partner_specificity_subunit" if specificity_only else "",
                          "map": "abi_cds"},
                source_cds=[cid, mt_id], rule="RM_pair",
            )
        else:
            consumed_defense.add(cid)
            add_row(
                contig=rec["contig"], start=rec["start"], end=rec["end"], strand=rec["strand"],
                element_type="defense", subtype=f"{d['subtype']}_orphan", gene=gene_of(rec),
                confidence=conf(max(1, d["score"] - 1)),
                evidence={"system": "RM", "name": rec["product"], "note": "no_methylase_within_5kb",
                          "map": "abi_cds"},
                source_cds=[cid], rule="RM_orphan",
            )

    toxins, antitoxins, others = [], [], []
    for cid, meta in info.items():
        d = meta["defense"]
        if not d or cid in consumed_defense:
            continue
        if d["role"] == "TA_toxin":
            toxins.append(cid)
        elif d["role"] == "TA_antitoxin":
            antitoxins.append(cid)
        elif d["role"] == "other":
            others.append(cid)
    for cid in toxins:
        rec = info[cid]["rec"]
        d = info[cid]["defense"]
        pair = nearby(cid, antitoxins, 3_000)
        if pair:
            gap, at_id = pair
            at = info[at_id]["rec"]
            consumed_defense.update({cid, at_id})
            add_row(
                contig=rec["contig"], start=min(int(rec["start"]), int(at["start"])),
                end=max(int(rec["end"]), int(at["end"])), strand=rec["strand"],
                element_type="defense", subtype=d["subtype"],
                gene=f"{gene_of(rec)}/{gene_of(at)}",
                confidence="high" if d["score"] >= 3 else "medium",
                evidence={"system": "TA", "name": f"toxin={rec['product']}|antitoxin={at['product']}",
                          "map": "abi_cds"},
                source_cds=[cid, at_id], rule="TA_pair",
            )
        else:
            add_row(
                contig=rec["contig"], start=rec["start"], end=rec["end"], strand=rec["strand"],
                element_type="defense", subtype=d["subtype"], gene=gene_of(rec),
                confidence=conf(max(1, d["score"] - 1)),
                evidence={"system": "TA", "name": rec["product"], "note": "no_antitoxin_within_3kb",
                          "map": "abi_cds"},
                source_cds=[cid], rule="TA_single",
            )
            consumed_defense.add(cid)
    for cid in others:
        rec = info[cid]["rec"]
        d = info[cid]["defense"]
        add_row(
            contig=rec["contig"], start=rec["start"], end=rec["end"], strand=rec["strand"],
            element_type="defense", subtype=d["subtype"], gene=gene_of(rec),
            confidence=conf(d["score"]), evidence={"system": "other", "name": rec["product"],
                                                  "map": "abi_cds"},
            source_cds=[cid], rule="defense_single",
        )
        consumed_defense.add(cid)

    # ---- 4. Cas ------------------------------------------------------------
    for cid, meta in info.items():
        cas = meta["cas"]
        if not cas:
            continue
        rec = meta["rec"]
        add_row(
            contig=rec["contig"], start=rec["start"], end=rec["end"], strand=rec["strand"],
            element_type="cas", subtype=cas["subtype"], gene=cas["gene"] or gene_of(rec),
            confidence=conf(cas["score"]),
            evidence={"name": rec["product"],
                      "pfam": f"{cas['pfam'][1]['pfam_name']}({cas['pfam'][1]['pfam_acc']},GA={cas['pfam'][1]['passes_ga']})" if cas["pfam"] else "",
                      "prior": meta["prior_str"], "map": "abi_cds"},
            source_cds=[cid], rule="cas_gene",
            source="prior+de_novo" if meta["priors"] else "de_novo",
        )

    # ---- 5. AMR ------------------------------------------------------------
    for cid, meta in info.items():
        amr = meta["amr"]
        if not amr:
            continue
        rec = meta["rec"]
        add_row(
            contig=rec["contig"], start=rec["start"], end=rec["end"], strand=rec["strand"],
            element_type="amr", subtype=amr["subtype"], gene=gene_of(rec),
            confidence=conf(amr["score"]),
            evidence={"name": rec["product"],
                      "pfam": f"{amr['pfam'][1]['pfam_name']}({amr['pfam'][1]['pfam_acc']},GA={amr['pfam'][1]['passes_ga']})" if amr["pfam"] else "",
                      "prior": meta["prior_str"], "map": "abi_cds",
                      "note": "point_mutation_not_assessed" if amr["subtype"] == "target_mutation_unverified" else ""},
            source_cds=[cid], rule="amr_gene",
            source="prior+de_novo" if meta["priors"] else "de_novo",
        )

    # ---- 6. móviles (fuera de profagos) -----------------------------------
    for cid, meta in info.items():
        mob = meta["mobile"]
        if not mob or cid in region_cds or meta["cas"] or meta["amr"] or meta["defense"]:
            continue
        rec = meta["rec"]
        add_row(
            contig=rec["contig"], start=rec["start"], end=rec["end"], strand=rec["strand"],
            element_type="transposon", subtype=mob["subtype"], gene=gene_of(rec),
            confidence=conf(mob["score"]),
            evidence={"name": rec["product"],
                      "pfam": f"{mob['pfam'][1]['pfam_name']}({mob['pfam'][1]['pfam_acc']},GA={mob['pfam'][1]['passes_ga']})" if mob["pfam"] else "",
                      "prior": meta["prior_str"], "map": "abi_cds"},
            source_cds=[cid], rule="mobile_gene",
            source="prior+de_novo" if meta["priors"] else "de_novo",
        )

    # ---- 7. virulencia -----------------------------------------------------
    for cid, meta in info.items():
        vir = meta["vir"]
        if not vir or meta["cas"] or meta["amr"] or meta["defense"] or meta["mobile"]:
            continue
        rec = meta["rec"]
        add_row(
            contig=rec["contig"], start=rec["start"], end=rec["end"], strand=rec["strand"],
            element_type="virulence", subtype=vir["subtype"], gene=gene_of(rec),
            confidence=conf(vir["score"]),
            evidence={"name": rec["product"],
                      "pfam": f"{vir['pfam'][1]['pfam_name']}({vir['pfam'][1]['pfam_acc']})" if vir["pfam"] else "",
                      "prior": meta["prior_str"], "map": "abi_cds"},
            source_cds=[cid], rule="virulence_gene",
            source="prior+de_novo" if meta["priors"] else "de_novo",
        )

    # ---- 8. arrays CRISPR (con contexto cas) ------------------------------
    cas_cds = [meta["rec"] for meta in info.values() if meta["cas"]]
    for arr in arrays:
        start, end = int(arr["start"]), int(arr["end"])
        nearest = None
        for rec in cas_cds:
            if rec["contig"] != arr["contig"]:
                continue
            dist = max(0, max(start, int(rec["start"])) - min(end, int(rec["end"])))
            if nearest is None or dist < nearest[0]:
                nearest = (dist, rec)
        strand = "."
        cas_note = ""
        confianza = arr["confidence"]
        if nearest and nearest[0] <= 10_000:
            strand = nearest[1]["strand"]
            cas_note = f"{nearest[1]['gene'] or nearest[1]['locus_tag']}@{nearest[0]}bp"
            if confianza != "high":
                confianza = "high"
        add_row(
            contig=arr["contig"], start=start, end=end, strand=strand,
            element_type="crispr", subtype="array", gene=arr["array_id"],
            confidence=confianza,
            evidence={"repeats": arr["n_repeats"], "spacers": int(arr["n_repeats"]) - 1,
                      "interrupted": arr.get("interrupted", "no"), "cas": cas_note,
                      "note": f"repeat_len={arr['repeat_len']};identity={arr['mean_identity']};scan=crt_like"},
            source_cds=[], rule="crispr_array",
        )

    # ---- 9. previos sin CDS ab initio (ncbi_only) y demotions -------------
    demoted = []
    for p in priors:
        if p["resolution"] == "unplaced":
            demoted.append({**p, "reason": p.get("note") or "unplaced"})
            continue
        if p["cds_id"]:
            continue  # se evalúa abajo con cobertura
        # ncbi_only: solo nombre
        section_type = {"mobile": "transposon", "crispr": "cas", "prophage": "prophage",
                        "defense": "defense", "amr": "amr", "virulence": "virulence"}[p["source_section"]]
        product = p["product"] or p["prior_name"]
        found = None
        subtype, spec = "", 0
        if section_type == "prophage":
            # si el gen de fago previo ya cae dentro de una región emitida, la
            # región lo cubre y no hace falta una fila suelta
            inside = any(
                r["element_type"] == "prophage" and r["contig"] == p["contig"]
                and int(r["start"]) <= int(p["end"]) and int(p["start"]) <= int(r["end"])
                for r in rows
            )
            if inside:
                continue
        if section_type == "cas":
            found = name_cas(product)
            if found:
                subtype, spec = found[0], found[2]
        elif section_type == "amr":
            found = name_amr(product)
            if found:
                subtype, spec = found[0], found[2]
        elif section_type == "defense":
            found = name_defense(product)
            if found:
                subtype, spec = found[1], found[3]
        elif section_type == "transposon":
            found = name_mobile(product)
            if found:
                subtype, spec = found[0], found[2]
        elif section_type == "virulence":
            found = name_virulence(product)
            if found:
                subtype, spec = found[0], found[2]
        else:  # prophage
            found = phage_markers(product)
            if found:
                subtype, spec = "phage_gene", 2
        if not found or not subtype:
            demoted.append({**p, "reason": f"no_specific_evidence_for_{section_type}"})
            continue
        add_row(
            contig=p["contig"], start=p["start"], end=p["end"], strand=p["strand"] or ".",
            element_type=section_type, subtype=subtype, gene=p["gene"] or p["locus_tag"] or p["prior_name"],
            confidence=conf(spec),
            evidence={"prior": f"{p['source_json']}/{p['source_section']}", "name": product,
                      "map": "ncbi_only", "note": "no_ab_initio_cds"},
            source_cds=[], rule="prior_ncbi_only",
        )

    # cobertura de previos con CDS
    for p in priors:
        if not p["cds_id"] or p["resolution"] == "unplaced":
            continue
        rows_here = rows_by_cds.get(p["cds_id"], [])
        if not rows_here:
            demoted.append({**p, "reason": "cds_not_classified"})
            continue
        sections = {r["element_type"] for r in rows_here}
        expected = {"mobile": "transposon", "crispr": "cas", "prophage": "prophage",
                    "defense": "defense", "amr": "amr", "virulence": "virulence"}[p["source_section"]]
        if expected not in sections:
            demoted.append({**p, "reason": f"reclassified_as_{'/'.join(sorted(sections))}"})

    # ---- 10. salida --------------------------------------------------------
    with (WORK / genome / "elements.tsv").open("w", newline="") as fh:
        fields = ["genome", "contig", "start", "end", "strand", "element_type",
                  "subtype", "gene", "evidence", "confidence", "source_cds", "rule", "source"]
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    with (WORK / genome / "demoted.tsv").open("w", newline="") as fh:
        fields = ["genome", "source_json", "source_section", "prior_name", "product",
                  "contig", "start", "end", "resolution", "reason"]
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for p in demoted:
            writer.writerow({k: p.get(k, "") for k in fields})

    by_type = Counter(r["element_type"] for r in rows)
    summary = " ".join(f"{t}={by_type.get(t, 0)}" for t in
                      ["crispr", "cas", "transposon", "prophage", "defense", "amr", "virulence"])
    print(f"[{genome}] filas={len(rows)} ({summary}) descartados={len(demoted)}")
    return rows, demoted


def main() -> int:
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--shard", type=int, default=1)
    parser.add_argument("--shards", type=int, default=1)
    args, _ = parser.parse_known_args()
    genomes = load_pilots()
    if args.shards > 1:
        genomes = genomes[args.shard - 1 :: args.shards]
    failures = {}
    for genome in genomes:
        try:
            classify(genome)
        except Exception as exc:  # noqa: BLE001
            failures[genome] = f"{type(exc).__name__}: {exc}"
            print(f"[{genome}] FALLO: {failures[genome]}", file=sys.stderr)
    (WORK / "state").mkdir(parents=True, exist_ok=True)
    suffix = "" if args.shards == 1 else f"_{args.shard}"
    (WORK / "state" / f"classify_failures{suffix}.json").write_text(
        json.dumps(failures, indent=2, sort_keys=True) + "\n"
    )
    print(f"FIN shard {args.shard}/{args.shards}: {len(genomes) - len(failures)} ok, "
          f"{len(failures)} fallos")
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
