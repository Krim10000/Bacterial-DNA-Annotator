"""browser.py — navegador de genoma en Streamlit.

Dibuja una ventana del genoma con pistas apiladas (genes, operones, promotores,
horquillas) sobre un eje de coordenadas común, y permite buscar por coordenada,
gen, producto o reacción.

Principios de diseño
--------------------
- **Solo lectura**: no escribe datos del proyecto.
- **Sin dependencias nuevas**: el dibujo es SVG generado aquí; Streamlit solo lo
  incrusta. Nada de plotly ni componentes externos.
- **Un eje, varias pistas, una principal**: el usuario elige qué pista manda y
  las demás se atenúan, para no saturar.
- **Nada de decisiones automáticas**: aquí se muestra lo que otros pasos
  calcularon, con su clase y su confianza.

Uso desde `app.py`
------------------
    import browser
    browser.render(data_root=DATA_ROOT, genome=selected_genome)
"""
from __future__ import annotations

from pathlib import Path
import csv

import pandas as pd
import streamlit as st
from streamlit.components.v1 import declare_component

from i18n import tr

# ---------------------------------------------------------------- constantes
GENE_FWD = "#4a7ebb"
GENE_REV = "#b07a3a"
HAIRPIN = {"A_fuerte": "#b30000", "B_media": "#c58a00", "C_debil": "#9aa0ad"}
HAIRPIN_MEANING = {
    "A_fuerte": "estructura clara, tallo largo y GC por encima del fondo del genoma",
    "B_media": "intermedia: no se distingue del azar",
    "C_debil": "débil, probablemente sin función de terminación",
}
PROMOTER = "#177245"
TSS_C = "#177245"    # medido (sólido) o alto rendimiento (atenuado)
TTS_C = "#8a5a22"
OPERON = "#8e6cc9"

WINDOW_SIZES = {"5 kb": 5_000, "10 kb": 10_000, "25 kb": 25_000,
                "50 kb": 50_000, "100 kb": 100_000}

# Componente propio: convierte el SVG en algo clicable sin añadir dependencias.
_COMPONENT_DIR = Path(__file__).resolve().parent / "components" / "genome_track"
_genome_track = declare_component("genome_track", path=str(_COMPONENT_DIR))


def genome_track(svg_markup: str, key: str = "genome_track", *,
                 pos: int | None = None, contig: str | None = None,
                 window: int | None = None) -> dict | None:
    """Dibuja el SVG y devuelve el elemento pulsado (o None).

    ``pos``/``contig``/``window`` permiten que el componente construya una URL
    de respaldo si Streamlit no procesa el valor del clic.
    """
    value = _genome_track(svg=svg_markup, show_hint=True,
                          pos=pos, contig=contig, window=window, key=key)
    return value if isinstance(value, dict) else None


# ------------------------------------------------------------------ lectura
def _p(data_root: Path, *parts: str) -> Path:
    return Path(data_root).joinpath(*parts)


DEMO_DIR = Path(__file__).resolve().parent / "demo"


def load_tables(data_root: Path, genome: str) -> dict:
    """Carga las pistas de un genoma. Las que falten se devuelven vacías.

    El genoma ``demo`` se lee de ``webapp/demo/`` si no hay datos locales, para
    que el navegador funcione en un clon recién descargado.
    """
    reg = _p(data_root, "regulation", genome)
    if genome == "demo" and not (reg / "genes.tsv").exists():
        reg = DEMO_DIR
    out: dict[str, pd.DataFrame] = {}
    for name in ("genes", "operons", "promoters"):
        path = reg / f"{name}.tsv"
        out[name] = (pd.read_csv(path, sep="\t", dtype=str).fillna("")
                     if path.exists() else pd.DataFrame())
    hp = _p(data_root, "operon_finder", f"{genome}_hairpins.tsv")
    if not hp.exists() and genome == "demo":
        hp = DEMO_DIR / "hairpins.tsv"
    out["hairpins"] = (pd.read_csv(hp, sep="\t", dtype=str).fillna("")
                       if hp.exists() else pd.DataFrame())
    cx = _p(data_root, "complexes", f"{genome}.complexes.tsv")
    if not cx.exists() and genome == "demo":
        cx = DEMO_DIR / "complexes.tsv"
    out["complexes"] = (pd.read_csv(cx, sep="\t", dtype=str).fillna("")
                        if cx.exists() else pd.DataFrame())
    for key in ("genes", "operons", "complexes"):
        df = out[key]
        for col in ("start", "end"):
            if col in df:
                df[col] = pd.to_numeric(df[col], errors="coerce")
    df = out["hairpins"]
    for col in ("start", "end", "dG", "stem"):
        if col in df:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return out


TRANS_FILES = {
    "ecoli_k12": {
        "tss": ("regulondb", "tss_k12_union.tsv"),
        "tts": ("regulondb", "TerminatorSet.tsv"),
    },
    "ecoli": {"tss": ("regulondb", "tss_k12_union.tsv"),
              "tts": ("regulondb", "TerminatorSet.tsv")},
}


def load_transcription(data_root: Path, genome: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """TSS y TTS del genoma, con su procedencia (medido vs inferido).

    La distinción importa: un TSS curado es una medición; los de alto
    rendimiento son mediciones de otro estudio; y una predicción sería una
    inferencia. El navegador lo dice en cada marca.
    """
    root = Path(data_root)
    tss_rows, tts_rows = [], []
    # K-12: unión de curados + alto rendimiento
    tss_path = root / "regulondb" / "tss_k12_union.tsv"
    if genome.startswith("ecoli") and tss_path.exists():
        for r in csv.DictReader(open(tss_path), delimiter="\t"):
            if not (r.get("pos") or "").isdigit():
                continue
            src = "alto rendimiento" if (r.get("sigma") or "").startswith("HT") else "curado"
            tss_rows.append({"pos": int(r["pos"]), "strand": r.get("strand", "+"),
                             "source": src, "detail": r.get("sigma", "")})
    # B. subtilis: BSGatlas
    bsg = root / "dbtbs" / "gt_tss.tsv"
    if genome in ("subtilis", "bacillus") and bsg.exists():
        for r in csv.DictReader(open(bsg), delimiter="\t"):
            if (r.get("pos") or "").isdigit():
                tss_rows.append({"pos": int(r["pos"]), "strand": r.get("strand", "+"),
                                 "source": "BSGatlas", "detail": r.get("sigma", "")})
    # TTS
    term = root / "regulondb" / "TerminatorSet.tsv"
    if genome.startswith("ecoli") and term.exists():
        for r in csv.DictReader(open(term), delimiter="\t"):
            v = (r.get("2)leftEndPos") or "").strip()
            if v.isdigit():
                st = (r.get("4)strand") or "").strip().lower()
                tts_rows.append({"pos": int(v), "strand": "-" if st.startswith("r") else "+",
                                 "source": "curado", "detail": (r.get("7)type") or "")})
    bterm = root / "dbtbs" / "gt_terminators.tsv"
    if genome in ("subtilis", "bacillus") and bterm.exists():
        for r in csv.DictReader(open(bterm), delimiter="\t"):
            if (r.get("start") or "").isdigit():
                tts_rows.append({"pos": int(r["start"]), "strand": r.get("strand", "+"),
                                 "source": "BSGatlas", "detail": ""})
    return pd.DataFrame(tss_rows), pd.DataFrame(tts_rows)


def genome_gc(data_root: Path, genome: str) -> float | None:
    """GC del genoma a partir de los genes (suficiente para un rótulo)."""
    path = _p(data_root, "regulation", genome, "genes.tsv")
    if not path.exists():
        return None
    df = pd.read_csv(path, sep="\t", dtype=str, usecols=lambda c: c in ("gc",))
    if "gc" in df:
        v = pd.to_numeric(df["gc"], errors="coerce").dropna()
        if len(v):
            return float(v.iloc[0])
    return None


# ----------------------------------------------------------------- búsqueda
def search(data_root: Path, genome: str, query: str, tables: dict) -> list[dict]:
    """Busca por coordenadas, gen, producto, EC o reacción.

    Devuelve entradas {etiqueta, contig, pos, tipo} listas para saltar.
    """
    q = (query or "").strip()
    if not q:
        return []
    hits: list[dict] = []
    genes = tables.get("genes", pd.DataFrame())
    contigs = genes["contig"].unique().tolist() if "contig" in genes else []

    # 1) coordenadas "1000000-1050000" o "1000000"
    parts = q.replace("–", "-").split("-")
    if parts and all(p.strip().isdigit() for p in parts if p.strip()):
        if len(parts) == 2:
            a, b = int(parts[0]), int(parts[1])
            return [{"etiqueta": f"Región {a:,}–{b:,}", "contig": None,
                     "pos": a, "fin": b, "tipo": "coordenada"}]
        return [{"etiqueta": f"Posición {int(parts[0]):,}", "contig": None,
                 "pos": int(parts[0]), "fin": None, "tipo": "coordenada"}]

    # 2) reacción Rhea
    if q.upper().startswith("RHEA:"):
        gm = _p(data_root, "metabolism", "genemaps", f"{genome}_gene_map.tsv")
        if gm.exists():
            df = pd.read_csv(gm, sep="\t", dtype=str).fillna("")
            for _, r in df[df["rhea_id"].str.upper() == q.upper()].head(20).iterrows():
                hits.append({"etiqueta": f"{r.get('gene_id','')} · {r.get('product','')[:60]}",
                             "contig": r.get("contig"), "pos": int(r.get("start") or 0),
                             "fin": None, "tipo": "reacción"})
        return hits

    # 3) EC "2.7.1.39"
    if q.count(".") == 3 and q.replace(".", "").isdigit() and "ec" in genes:
        for _, r in genes[genes["ec"].astype(str).str.startswith(q)].head(20).iterrows():
            hits.append({"etiqueta": f"{r.get('locus_tag','')} · EC {r.get('ec','')}",
                         "contig": r.get("contig"), "pos": int(r.get("start") or 0),
                         "fin": None, "tipo": "EC"})
        return hits

    # 4) nombre de gen, locus o producto
    if not genes.empty:
        mask = pd.Series(False, index=genes.index)
        for col in ("locus_tag", "gene_name", "gene_id", "product"):
            if col in genes:
                mask |= genes[col].astype(str).str.contains(q, case=False, na=False)
        for _, r in genes[mask].head(25).iterrows():
            etiqueta = r.get("gene_name") or r.get("locus_tag") or r.get("gene_id")
            hits.append({"etiqueta": f"{etiqueta} · {str(r.get('product',''))[:60]}",
                         "contig": r.get("contig"), "pos": int(r.get("start") or 0),
                         "fin": None, "tipo": "gen"})
    return hits


# --------------------------------------------------------------------- svg
def _esc(s) -> str:
    return (str(s).replace("&", "&amp;").replace("<", "&lt;")
            .replace(">", "&gt;"))


def svg(tables: dict, contig: str, start: int, end: int,
        tracks: set[str] | None = None, focus: str = "",
        suspects: pd.DataFrame | None = None,
        tss: pd.DataFrame | None = None, tts: pd.DataFrame | None = None) -> str:
    """Dibuja la ventana como SVG ancho.

    Cada pista vive en su propia banda vertical, así que nunca se solapan: se
    pueden mostrar todas a la vez sin que se pisen. Las etiquetas de los genes
    se escalonan en dos filas cuando están muy juntas.
    """
    tracks = tracks or {"genes", "operones", "promotores", "horquillas",
                        "transcripcion", "complejos"}
    span = max(1, end - start)
    width, height = 1180, 420
    pad_l, pad_r = 96, 12
    plot_w = width - pad_l - pad_r

    def x(bp: float) -> float:
        return pad_l + (bp - start) / span * plot_w

    genes = tables.get("genes", pd.DataFrame())
    operons = tables.get("operons", pd.DataFrame())
    promoters = tables.get("promoters", pd.DataFrame())
    hairpins = tables.get("hairpins", pd.DataFrame())
    complexes = tables.get("complexes", pd.DataFrame())
    if "contig" in genes and contig:
        genes = genes[genes["contig"] == contig]
        if "contig" in hairpins:
            hairpins = hairpins[hairpins["contig"] == contig]
    if not genes.empty:
        genes = genes[(genes["end"] >= start) & (genes["start"] <= end)]
    if not hairpins.empty:
        hairpins = hairpins[(hairpins["end"] >= start) & (hairpins["start"] <= end)]
    if not operons.empty and "contig" in operons:
        operons = operons[(operons["contig"] == contig)
                          & (operons["end"] >= start) & (operons["start"] <= end)]
    if not complexes.empty and "contig" in complexes:
        complexes = complexes[(complexes["contig"] == contig)
                              & (complexes["end"] >= start)
                              & (complexes["start"] <= end)]

    Y = {"ruler": 30, "genes": 92, "operones": 166, "promotores": 214,
         "horquillas": 262, "transcripcion": 306, "complejos": 350,
         "sospecha": 392}
    o = [f'<svg xmlns="http://www.w3.org/2000/svg" width="100%" height="{height}" '
         f'viewBox="0 0 {width} {height}" font-family="sans-serif" '
         f'style="background:#fff;border-radius:8px">']

    # regla
    o.append(f'<line x1="{pad_l}" y1="{Y["ruler"]}" x2="{width-pad_r}" y2="{Y["ruler"]}" stroke="#c9cfda"/>')
    for k in range(6):
        bp = start + span * k / 5
        xx = x(bp)
        o.append(f'<line x1="{xx:.1f}" y1="{Y["ruler"]-5}" x2="{xx:.1f}" y2="{Y["ruler"]+5}" stroke="#c9cfda"/>')
        anchor = "start" if k == 0 else ("end" if k == 5 else "middle")
        o.append(f'<text x="{xx:.1f}" y="{Y["ruler"]-10}" font-size="11" fill="#8b90a0" '
                 f'text-anchor="{anchor}">{bp:,.0f}</text>')

    # rótulos de banda
    for key, label in (("genes", "GENES"), ("operones", "OPERONES"),
                       ("promotores", "PROMOTORES"), ("horquillas", "HORQUILLAS"),
                       ("transcripcion", "TSS / TTS"),
                       ("complejos", "COMPLEJOS")):
        active = key in tracks
        o.append(f'<text x="8" y="{Y[key]-6}" font-size="9.5" fill="#808495" '
                 f'font-weight="700" opacity="{1 if active else 0.25}">{label}</text>')

    yg = Y["genes"]
    # genes + etiquetas escalonadas
    if "genes" in tracks:
        labels: list[tuple[float, str]] = []
        for _, g in genes.iterrows():
            gx0, gx1 = x(g["start"]), x(g["end"])
            w = max(3.0, gx1 - gx0)
            rev = str(g.get("strand", "1")) in ("-1", "-")
            color = GENE_REV if rev else GENE_FWD
            name = g.get("gene_name") or g.get("locus_tag") or g.get("gene_id") or ""
            gid = (g.get("locus_tag") or g.get("gene_id") or g.get("gene_name")
                   or f"{int(g['start'])}-{int(g['end'])}")
            tip_txt = (f"{name or gid} · {int(g['start']):,}–{int(g['end']):,} pb · "
                       f"hebra {'−' if rev else '+'}")
            o.append(f'<g class="hit" data-kind="gene" data-id="{_esc(gid)}">'
                     f'<title>{_esc(tip_txt)}</title>')
            # zona de clic algo mayor que el rectángulo visible
            o.append(f'<rect x="{gx0-2:.1f}" y="{yg-16}" width="{max(6.0, w+4):.1f}" '
                     f'height="32" fill="transparent"/>')
            o.append(f'<rect x="{gx0:.1f}" y="{yg-11}" width="{w:.1f}" height="22" fill="{color}"/>')
            if w > 10:
                tip = min(8.0, w / 3)
                if rev:
                    o.append(f'<polygon points="{gx0:.1f},{yg-11} {gx0:.1f},{yg+11} '
                             f'{gx0-tip:.1f},{yg}" fill="{color}"/>')
                else:
                    o.append(f'<polygon points="{gx1:.1f},{yg-11} {gx1:.1f},{yg+11} '
                             f'{gx1+tip:.1f},{yg}" fill="{color}"/>')
            if focus and focus in (g.get("locus_tag"), g.get("gene_id"), g.get("gene_name")):
                o.append(f'<rect x="{gx0-3:.1f}" y="{yg-15}" width="{w+6:.1f}" height="30" '
                         f'fill="none" stroke="#ff4b4b" stroke-width="2"/>')
            o.append("</g>")
            if w > 24 and name:
                labels.append((gx0 + w / 2, str(name)))
        # dos filas para que no se pisen
        labels.sort()
        last_x, row = -999, 0
        for lx, name in labels:
            if lx - last_x < 52:
                row = 1 - row
            else:
                row = 0
            last_x = lx
            o.append(f'<text x="{lx:.1f}" y="{yg+30+row*15}" font-size="10" fill="#3b3f4a" '
                     f'text-anchor="middle">{_esc(name)[:16]}</text>')

    # operones
    if "operones" in tracks:
        for _, op in operons.iterrows():
            gx0, gx1 = x(op["start"]), x(op["end"])
            lab = op.get("label") or op.get("operon_id") or ""
            oid = op.get("operon_id") or lab or f"{int(op['start'])}-{int(op['end'])}"
            o.append(f'<g class="hit" data-kind="operon" data-id="{_esc(oid)}">'
                     f'<title>{_esc(lab or oid)} · operón</title>')
            o.append(f'<rect x="{gx0:.1f}" y="{Y["operones"]-8}" width="{max(4, gx1-gx0):.1f}" '
                     f'height="16" rx="8" fill="none" stroke="{OPERON}" stroke-width="1.8"/>')
            if gx1 - gx0 > 70 and lab:
                o.append(f'<text x="{gx0+5:.1f}" y="{Y["operones"]+3}" font-size="9" '
                         f'fill="{OPERON}">{_esc(lab)[:30]}</text>')
            o.append("</g>")

    # promotores
    if "promotores" in tracks and not promoters.empty and not operons.empty:
        by_id = {str(r.get("operon_id")): r for _, r in operons.iterrows()}
        for _, p in promoters.iterrows():
            op = by_id.get(str(p.get("operon_id")))
            if op is None:
                continue
            gx = x(op["start"])
            sigma = p.get("sigma", "")
            pid = p.get("operon_id") or op.get("operon_id") or f"{int(op['start'])}"
            o.append(f'<g class="hit" data-kind="promoter" data-id="{_esc(pid)}">'
                     f'<title>Promotor {_esc(str(pid))} · σ{_esc(str(sigma))}</title>')
            o.append(f'<line x1="{gx:.1f}" y1="{Y["promotores"]-9}" x2="{gx:.1f}" '
                     f'y2="{Y["promotores"]+9}" stroke="{PROMOTER}" stroke-width="2"/>')
            o.append(f'<circle cx="{gx:.1f}" cy="{Y["promotores"]-11}" r="3" fill="{PROMOTER}"/>')
            o.append("</g>")
            # el sigma no se repite en el dibujo: va en la leyenda y en el detalle

    # horquillas
    if "horquillas" in tracks:
        for hi, (_, h) in enumerate(hairpins.iterrows()):
            gx = x(h["start"])
            color = HAIRPIN.get(str(h.get("class", "C_debil")), HAIRPIN["C_debil"])
            hid = f"{contig}:{int(h['start'])}-{int(h['end'])}"
            o.append(f'<g class="hit" data-kind="hairpin" data-id="{_esc(hid)}">'
                     f'<title>Horquilla {_esc(str(h.get("class", "")))} · '
                     f'{int(h["start"]):,}–{int(h["end"]):,} pb</title>')
            o.append(f'<line x1="{gx:.1f}" y1="{Y["horquillas"]-10}" x2="{gx:.1f}" '
                     f'y2="{Y["horquillas"]+10}" stroke="{color}" stroke-width="1.6"/>')
            o.append(f'<circle cx="{gx:.1f}" cy="{Y["horquillas"]-12}" r="5" fill="none" '
                     f'stroke="{color}" stroke-width="2"/>')
            o.append("</g>")

    # TSS (triángulo hacia arriba) y TTS (hacia abajo): medidos vs inferidos
    if "transcripcion" in tracks:
        yt = Y["transcripcion"]
        for _, t in (tss if tss is not None else pd.DataFrame()).iterrows():
            if not (start <= t["pos"] <= end):
                continue
            gx = x(t["pos"])
            full = t.get("source") in ("curado", "BSGatlas")
            op = 1.0 if full else 0.55
            o.append(f'<g><title>TSS · {int(t["pos"]):,} pb · '
                     f'{_esc(str(t.get("source", "")))}</title>')
            o.append(f'<polygon points="{gx:.1f},{yt-16} {gx-4:.1f},{yt-9} '
                     f'{gx+4:.1f},{yt-9}" fill="{TSS_C}" opacity="{op}"/>')
            o.append(f'<line x1="{gx:.1f}" y1="{yt-9}" x2="{gx:.1f}" y2="{yt-3}" '
                     f'stroke="{TSS_C}" stroke-width="1" opacity="{op}"/>')
            o.append("</g>")
        for _, t in (tts if tts is not None else pd.DataFrame()).iterrows():
            if not (start <= t["pos"] <= end):
                continue
            gx = x(t["pos"])
            op = 1.0
            o.append(f'<g><title>TTS · {int(t["pos"]):,} pb</title>')
            o.append(f'<polygon points="{gx:.1f},{yt+16} {gx-4:.1f},{yt+9} '
                     f'{gx+4:.1f},{yt+9}" fill="{TTS_C}" opacity="{op}"/>')
            o.append(f'<line x1="{gx:.1f}" y1="{yt+3}" x2="{gx:.1f}" y2="{yt+9}" '
                     f'stroke="{TTS_C}" stroke-width="1" opacity="{op}"/>')
            o.append("</g>")

    # complejos por vecindad
    if "complejos" in tracks:
        for _, c in complexes.iterrows():
            gx0, gx1 = x(c["start"]), x(c["end"])
            status = str(c.get("status", ""))
            color = {"complex_complete": "#2e7d32",
                     "complex_incomplete": "#c58a00",
                     "complex_ambiguous": "#b3261e"}.get(status, "#8b90a0")
            cid = str(c.get("complex_id") or f"{int(c['start'])}-{int(c['end'])}")
            tip = (f"{c.get('name','complejo')} · {status} · "
                   f"{int(c['start']):,}–{int(c['end']):,} pb")
            o.append(f'<g class="hit" data-kind="complex" data-id="{_esc(cid)}">'
                     f'<title>{_esc(tip)}</title>')
            o.append(f'<rect x="{gx0:.1f}" y="{Y["complejos"]-8}" width="{max(6, gx1-gx0):.1f}" '
                     f'height="16" rx="8" fill="{color}" opacity="0.22" '
                     f'stroke="{color}" stroke-width="1.6"/>')
            if gx1 - gx0 > 70:
                o.append(f'<text x="{gx0+5:.1f}" y="{Y["complejos"]+3}" font-size="8.5" '
                         f'fill="{color}">{_esc(str(c.get("name", "")))[:32]}</text>')
            o.append("</g>")

    # sospecha por cadena metabólica, como banda propia
    if suspects is not None and not suspects.empty:
        o.append(f'<text x="8" y="{Y["sospecha"]-6}" font-size="9.5" fill="#c58a00" '
                 f'font-weight="700">SOSPECHA</text>')
        for si, (_, s) in enumerate(suspects.iterrows()):
            gx0, gx1 = x(s["start"]), x(s["end"])
            sid = f"{int(s['start'])}-{int(s['end'])}"
            o.append(f'<g class="hit" data-kind="suspect" data-id="{_esc(sid)}">'
                     f'<title>Sospecha metabólica · {_esc(str(s.get("n_genes")))} genes</title>')
            o.append(f'<rect x="{gx0:.1f}" y="{Y["sospecha"]-8}" width="{max(6, gx1-gx0):.1f}" '
                     f'height="16" rx="8" fill="#ffe9c7" stroke="#c58a00" stroke-width="1.6"/>')
            if gx1 - gx0 > 60:
                o.append(f'<text x="{gx0+5:.1f}" y="{Y["sospecha"]+3}" font-size="8.5" '
                         f'fill="#9a5b00">{_esc(str(s.get("n_genes")))} genes · cadena</text>')
            o.append("</g>")

    o.append("</svg>")
    return "".join(o)


# ------------------------------------------------------------------ página
def suspects_window(data_root: Path, genome: str, contig: str,
                    start: int, end: int) -> pd.DataFrame:
    """Sospechas por cadena metabólica que caen en la ventana."""
    path = Path(data_root) / "metabolism" / "operon_suspects.tsv"
    if not path.exists():
        return pd.DataFrame()
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    df = df[df["genome"] == genome]
    if "contig" in df and contig:
        df = df[df["contig"] == contig]
    for col in ("start", "end"):
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df[(df["end"] >= start) & (df["start"] <= end)]
    return df.sort_values("start").reset_index(drop=True)


# ------------------------------------------------------- ficha al hacer clic
def _match_any(df: pd.DataFrame, value, cols: tuple[str, ...]) -> pd.Series:
    """Máscara booleana tolerante a columnas ausentes."""
    if df.empty:
        return pd.Series(dtype=bool)
    mask = pd.Series(False, index=df.index)
    for col in cols:
        if col in df:
            mask |= df[col].astype(str) == str(value)
    return mask


def _row(df: pd.DataFrame, mask) -> dict | None:
    sub = df[mask] if not df.empty else df
    return sub.iloc[0].to_dict() if len(sub) else None


def _gene_detail(tables: dict, gid: str) -> None:
    g = _row(tables["genes"], _match_any(
        tables["genes"], gid, ("locus_tag", "gene_id", "gene_name")))
    if not g:
        st.warning(tr("Ese gen ya no está en la ventana actual."))
        return
    name = g.get("gene_name") or g.get("locus_tag") or g.get("gene_id")
    st.markdown(f"### 🧬 {name}")
    c1, c2, c3, c4 = st.columns(4)
    c1.metric(tr("Longitud"), f"{int(g['end']) - int(g['start']) + 1:,} pb")
    c2.metric(tr("Hebra"), "−" if str(g.get("strand")) in ("-1", "-") else "+")
    c3.metric(tr("Inicio"), f"{int(g['start']):,}")
    c4.metric(tr("Fin"), f"{int(g['end']):,}")
    st.markdown(f"**{tr('Producto')}:** {g.get('product') or tr('sin anotar')}")
    if g.get("locus_tag"):
        st.markdown(f"**Locus tag:** `{g['locus_tag']}`")
    if g.get("ec"):
        st.markdown(f"**EC:** `{g['ec']}`")
    if g.get("gene_type") or g.get("rna"):
        st.markdown(f"**{tr('Tipo')}:** {g.get('gene_type') or g.get('rna')}")
    st.caption(tr("Confianza: {c} · fuente: {s} · estado: {e}").format(
        c=g.get("confidence", "—"), s=g.get("source", "—"), e=g.get("status", "—")))


def _operon_detail(tables: dict, oid: str) -> None:
    o = _row(tables["operons"], _match_any(
        tables["operons"], oid, ("operon_id", "label")))
    if not o:
        st.warning(tr("Ese operón ya no está en la ventana actual."))
        return
    st.markdown(f"### 🧩 {o.get('label') or o.get('operon_id')}")
    c1, c2, c3 = st.columns(3)
    c1.metric(tr("Genes"), int(float(o.get("n_genes") or 0)))
    c2.metric(tr("Extensión"), f"{int(float(o.get('span_bp') or 0)):,} pb")
    c3.metric(tr("Hebra"), "−" if str(o.get("strand")) in ("-1", "-") else "+")
    st.markdown(f"**{tr('Genes')}:** {o.get('genes') or '—'}")
    st.markdown(f"**{tr('Nombres')}:** {o.get('gene_names') or '—'}")
    st.caption(tr("Hueco mínimo entre genes: {g} pb · evidencia: {e} · confianza: {c}").format(
        g=o.get("min_gap", "—"), e=o.get("evidence", "—"), c=o.get("confidence", "—")))


def _promoter_detail(tables: dict, oid: str) -> None:
    p = _row(tables["promoters"], _match_any(
        tables["promoters"], oid, ("operon_id", "label")))
    if not p:
        st.warning(tr("Ese promotor ya no está en la ventana actual."))
        return
    st.markdown(f"### 🎯 {tr('Promotor')} · {p.get('label') or p.get('operon_id')}")
    st.markdown(f"**σ:** {p.get('sigma') or '—'} · **score:** {p.get('score') or '—'}")
    st.markdown(f"**−35:** `{p.get('box35') or '—'}` ({tr('desajustes')} {p.get('mm35', '—')})")
    st.markdown(f"**−10:** `{p.get('box10') or '—'}` ({tr('desajustes')} {p.get('mm10', '—')})")
    st.caption(tr("Separación: {s} pb · distancia al inicio: {d} pb · sitios: {n}").format(
        s=p.get("spacer", "—"), d=p.get("dist_10_to_start", "—"), n=p.get("n_sites", "—")))


def _hairpin_detail(tables: dict, hid: str) -> None:
    hp = tables["hairpins"]
    mask = hp.apply(lambda r: f"{r.get('contig')}:{int(r['start'])}-{int(r['end'])}" == hid, axis=1)         if not hp.empty else pd.Series(dtype=bool)
    h = _row(hp, mask)
    if not h:
        st.warning(tr("Esa horquilla ya no está en la ventana actual."))
        return
    st.markdown(f"### 🌀 {tr('Horquilla')} {h.get('class', '')}")
    c1, c2, c3 = st.columns(3)
    c1.metric("ΔG", f"{h.get('dG', '—')} kcal/mol")
    c2.metric(tr("Tallo"), f"{h.get('stem', '—')} pb")
    c3.metric("GC", f"{h.get('gc', '—')}")
    st.markdown(f"**{tr('Bucle')}:** {h.get('loop', '—')} · **poli-U:** {h.get('polyU', '—')}")
    st.markdown(f"**{tr('Genes vecinos')}:** {h.get('left_gene') or '—'} ← → {h.get('right_gene') or '—'}")
    st.caption(HAIRPIN_MEANING.get(str(h.get("class")), ""))


def _tx_detail(tables: dict, kind: str, pos: str) -> None:
    label = "TSS" if kind == "tss" else "TTS"
    st.markdown(f"### 📡 {label} · {int(float(pos)):,} pb")
    st.caption(tr("Marca de {k}: {w}. Se muestra como triángulo en la pista de "
                  "transcripción.").format(
        k=label,
        w=tr("inicio de la transcripción (donde empieza el ARN)")
        if kind == "tss" else tr("fin de la transcripción (donde termina el ARN)")))


def _suspect_detail(data_root: Path, genome: str, sid: str) -> None:
    path = Path(data_root) / "metabolism" / "operon_suspects.tsv"
    if not path.exists():
        return
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    df = df[df["genome"] == genome]
    mask = df.apply(lambda r: f"{int(float(r['start']))}-{int(float(r['end']))}" == sid, axis=1)         if not df.empty else pd.Series(dtype=bool)
    r = _row(df, mask)
    if not r:
        st.warning(tr("Esa sospecha ya no está en la ventana actual."))
        return
    st.markdown(f"### 🔗 {tr('Sospecha metabólica')} · {r.get('n_genes')} genes")
    st.markdown(f"**{tr('Genes')}:** {r.get('genes') or '—'}")
    st.caption(tr("Genes contiguos que comparten metabolitos. Es una sospecha de "
                  "operón metabólico, no una prueba."))


def _complex_detail(tables: dict, ident: str) -> None:
    df = tables.get("complexes", pd.DataFrame())
    row = _row(df, _match_any(df, ident, ("complex_id",))) if not df.empty else None
    if not row:
        st.info(tr("Complejo no encontrado."))
        return
    st.markdown(f"### 🧩 {row.get('name', 'Complejo')}")
    st.caption(f"{row.get('status','')} · confianza {row.get('confidence','')}")
    st.write(f"**Subunidades:** {row.get('subunits','')}")
    st.write(f"**Esperadas:** {row.get('expected','')}")
    st.write(f"**Faltan:** {row.get('missing','') or '—'}")
    st.write(f"**Evidencia:** {row.get('evidence','')}")
    st.caption(row.get("limitation", ""))


def _render_pick(data_root: Path, genome: str, tables: dict, pick: dict) -> None:
    kind, ident = pick.get("kind"), str(pick.get("id", ""))
    with st.container(border=True):
        head, close = st.columns([6, 1])
        with head:
            st.caption(tr("Has hecho clic en el dibujo. Esta ficha explica ese elemento."))
        with close:
            if st.button("✖ " + tr("Cerrar"), key="br_close_pick"):
                st.session_state["br_dismissed"] = pick.get("ts")
                if st.query_params.get("sel"):
                    try:
                        del st.query_params["sel"]
                    except (KeyError, TypeError):
                        pass
                st.rerun()
        if kind == "gene":
            _gene_detail(tables, ident)
        elif kind == "operon":
            _operon_detail(tables, ident)
        elif kind == "promoter":
            _promoter_detail(tables, ident)
        elif kind == "hairpin":
            _hairpin_detail(tables, ident)
        elif kind in ("tss", "tts"):
            _tx_detail(tables, kind, ident)
        elif kind == "complex":
            _complex_detail(tables, ident)
        elif kind == "suspect":
            _suspect_detail(data_root, genome, ident)
        else:
            st.info(tr("Elemento sin ficha todavía."))


def render(data_root: Path, genome: str) -> None:
    """Dibuja la vista del navegador para un genoma."""
    st.subheader("🧬 Navegador de genoma")
    st.caption(
        "Ventana del genoma con sus pistas apiladas. El **eje es común**: cada "
        "posición en pb se alinea en todas las pistas, y cada pista tiene su "
        "banda para que no se pisen. Marca las que quieras ver."
    )

    tables = load_tables(data_root, genome)
    if genome == "demo":
        st.info(tr("Demo incluida: 150 kb de *Chlorobaculum tepidum* anotados "
                   "desde cero. Los genes son predicciones ab initio y los "
                   "operones, una agrupación por proximidad; sirve para probar "
                   "el navegador, no para conclusiones biológicas."))
    tss_df, tts_df = load_transcription(data_root, genome)
    genes = tables["genes"]
    if genes.empty:
        st.warning(f"No hay anotación para **{genome}**. Genera primero la capa "
                   "regulatoria con `run_pipeline.sh components`.")
        return
    contigs = genes["contig"].dropna().unique().tolist()
    if not contigs:
        st.warning("El genoma no tiene contigs anotados.")
        return

    if st.session_state.get("br_genome") != genome:
        ini = 50_000 if genome == "ecoli_k12" else int(genes["start"].min())
        st.session_state.update({"br_genome": genome, "br_contig": contigs[0],
                                 "br_start": ini, "br_size": 10_000})
    # enlace directo a una región: ?pos=123456&w=10+kb&sel=gene:b0048
    try:
        if st.query_params.get("pos"):
            st.session_state["br_start"] = int(st.query_params["pos"])
        wq = st.query_params.get("w")
        if wq in WINDOW_SIZES:
            st.session_state["br_size_sel"] = wq
        if st.query_params.get("contig") in contigs:
            st.session_state["br_contig"] = st.query_params["contig"]
    except (TypeError, ValueError):
        pass

    # --- búsqueda
    q = st.text_input(
        "🔎 Buscar por coordenada, gen, producto, EC o reacción (RHEA:…)",
        key="br_query", placeholder="p. ej. thrA, proteasa, 1000000-1050000, RHEA:12345",
    )
    if q:
        hits = search(data_root, genome, q, tables)
        if not hits:
            st.caption("Sin resultados en este genoma.")
        else:
            pick = st.selectbox(
                f"{len(hits)} resultado(s)", range(len(hits)),
                format_func=lambda i: f"{hits[i]['etiqueta']}  ·  {hits[i]['contig'] or ''} "
                                      f"{hits[i]['pos']:,}",
                key="br_hit",
            )
            if st.button("Saltar a este resultado", key="br_jump"):
                h = hits[pick]
                if h["contig"]:
                    st.session_state["br_contig"] = h["contig"]
                st.session_state["br_start"] = max(0, h["pos"] - 2_000)
                st.rerun()

    # --- controles
    c1, c2, c3, c4, c5 = st.columns([2, 2, 1.6, 1, 1])
    contig = c1.selectbox("Contig", contigs,
                          index=contigs.index(st.session_state["br_contig"])
                          if st.session_state["br_contig"] in contigs else 0,
                          key="br_contig_sel")
    st.session_state["br_contig"] = contig
    size_label = c2.selectbox("Ventana", list(WINDOW_SIZES), index=1, key="br_size_sel")
    size = WINDOW_SIZES[size_label]
    start = c3.number_input("Posición (pb)", 0, max(0, int(genes["end"].max())),
                            int(st.session_state["br_start"]), step=size, key="br_pos")
    st.session_state["br_start"] = int(start)
    if c4.button("◀", key="br_left"):
        st.session_state["br_start"] = max(0, int(start) - size)
        st.rerun()
    if c5.button("▶", key="br_right"):
        st.session_state["br_start"] = int(start) + size
        st.rerun()

    # --- pistas visibles: casillas, para poder superponer sin que se pisen
    t1, t2, t3, t4, t5, t6, t7 = st.columns(7)
    tracks: set[str] = set()
    if t1.checkbox("Genes", True, key="br_t_genes"):
        tracks.add("genes")
    if t2.checkbox("Operones", True, key="br_t_op"):
        tracks.add("operones")
    if t3.checkbox("Promotores", True, key="br_t_prom"):
        tracks.add("promotores")
    if t4.checkbox("Horquillas", True, key="br_t_hp"):
        tracks.add("horquillas")
    ver_sospecha = t5.checkbox("Sospecha metabólica", True, key="br_t_susp")
    if t6.checkbox("TSS / TTS", True, key="br_t_tx"):
        tracks.add("transcripcion")
    if t7.checkbox("Complejos", True, key="br_t_complex"):
        tracks.add("complejos")

    start, end = int(start), int(start) + size
    wg = genes[(genes["contig"] == contig)
               & (genes["end"] >= start) & (genes["start"] <= end)].sort_values("start")
    susp = suspects_window(data_root, genome, contig, start, end) if ver_sospecha else pd.DataFrame()
    st.caption(f"{genome} · {contig} · {start:,}–{end:,} pb · "
               f"{len(wg)} genes en la ventana")

    clicked = genome_track(
        svg(tables, contig, start, end, tracks=tracks, suspects=susp,
            tss=tss_df, tts=tts_df),
        key="genome_track",
    )
    if clicked and clicked.get("ts") != st.session_state.get("br_dismissed"):
        _render_pick(data_root, genome, tables, clicked)

    # --- distancias
    if len(wg) > 1:
        gaps = [wg.iloc[i]["start"] - wg.iloc[i - 1]["end"] - 1 for i in range(1, len(wg))]
        gaps = [g for g in gaps if g >= 0]
        if gaps:
            st.caption(f"Distancia entre genes contiguos: mediana "
                       f"{pd.Series(gaps).median():.0f} pb · mínimo {min(gaps)} pb · "
                       f"máximo {max(gaps)} pb")

    # --- leyenda
    st.markdown(
        "**Leyenda** · "
        '<span style="color:#b30000">●</span> horquilla clase A (fuerte) · '
        '<span style="color:#c58a00">●</span> clase B (intermedia) · '
        '<span style="color:#9aa0ad">●</span> clase C (débil) · '
        f'<span style="color:{PROMOTER}">│</span> promotor (la σ y las cajas se ven en el detalle) · '
        '<span style="color:#c58a00">▭</span> sospecha por cadena metabólica · '
        '<span style="color:#177245">▲</span> TSS medido · '
        '<span style="color:#177245;opacity:.55">▲</span> TSS alto rendimiento · '
        '<span style="color:#8a5a22">▼</span> TTS curado · '
        '<span style="color:#2e7d32">▭</span> complejo completo · '
        '<span style="color:#c58a00">▭</span> complejo incompleto',
        unsafe_allow_html=True,
    )

    # --- detalle por pestañas: aquí se "pincha" en cada elemento
    t_genes, t_ops, t_prom, t_hp, t_susp, t_tx = st.tabs(
        ["🧬 Genes", "🧩 Operones", "🎯 Promotores", "🌀 Horquillas",
         "🔗 Sospechas", "📡 Transcripción"]
    )

    def _elegir(df: pd.DataFrame, cols: list[str], key: str) -> int | None:
        """Tabla seleccionable; devuelve el índice de la fila elegida."""
        if df.empty:
            st.caption("Nada en esta ventana.")
            return None
        cols = [c for c in cols if c in df.columns]
        try:
            ev = st.dataframe(df[cols], width="stretch", hide_index=True,
                              on_select="rerun", selection_mode="single-row",
                              key=key)
            sel = ev["selection"]["rows"] if isinstance(ev, dict) else \
                  getattr(ev, "selection", {}).get("rows", [])
            if sel:
                return int(sel[0])
        except Exception:
            idx = st.selectbox("Elegir", range(len(df)), key=key + "_sb",
                               format_func=lambda i: " · ".join(
                                   str(df.iloc[i].get(c, ""))[:30] for c in cols[:3]))
            return int(idx)
        st.caption("Haz clic en una fila para ver el detalle.")
        return None

    with t_genes:
        i = _elegir(wg.reset_index(drop=True),
                    ["locus_tag", "gene_name", "product", "start", "end", "strand", "ec"],
                    "br_tab_gen")
        if i is not None:
            r = wg.reset_index(drop=True).iloc[i]
            st.markdown(f"### 🧬 {r.get('gene_name') or r.get('locus_tag')}")
            a, b, c = st.columns(3)
            a.markdown(f"**Producto**  \n{r.get('product','') or '—'}")
            b.markdown(f"**Coordenadas**  \n{r.get('contig')}:{r.get('start')}–{r.get('end')} · "
                       f"hebra {'−' if str(r.get('strand')) in ('-1','-') else '+'}")
            c.markdown(f"**EC / familia**  \n{r.get('ec','') or '—'} · {r.get('fam','') or '—'}")
            ops = tables["operons"]
            if not ops.empty and "genes" in ops:
                hit = ops[ops["genes"].astype(str).str.contains(
                    str(r.get("locus_tag")), regex=False)]
                if not hit.empty:
                    o = hit.iloc[0]
                    st.caption(f"Pertenece al operón **{o.get('label') or o.get('operon_id')}** "
                               f"({o.get('n_genes')} genes, {o.get('span_bp','')} pb).")
            # promotor del operón si lo hay
            prom = tables["promoters"]
            if not prom.empty and not ops.empty and "operon_id" in prom:
                hit = ops[ops["genes"].astype(str).str.contains(
                    str(r.get("locus_tag")), regex=False)]
                if not hit.empty:
                    p2 = prom[prom["operon_id"] == hit.iloc[0].get("operon_id")]
                    if not p2.empty:
                        pr = p2.iloc[0]
                        st.markdown(f"**Promotor**: −35 `{pr.get('box35','')}` · "
                                    f"−10 `{pr.get('box10','')}` · σ **{pr.get('sigma','')}** "
                                    f"(score {pr.get('score','')})")

    with t_ops:
        wops = tables["operons"]
        if not wops.empty and "contig" in wops:
            wops = wops[(wops["contig"] == contig)
                        & (wops["end"] >= start) & (wops["start"] <= end)]
        i = _elegir(wops.reset_index(drop=True),
                    ["label", "n_genes", "start", "end", "strand", "span_bp", "confidence"],
                    "br_tab_op")
        if i is not None:
            r = wops.reset_index(drop=True).iloc[i]
            st.markdown(f"### 🧩 {r.get('label') or r.get('operon_id')}")
            a, b, c = st.columns(3)
            a.markdown(f"**Genes**  \n{r.get('n_genes')} · {r.get('span_bp','')} pb")
            b.markdown(f"**Coordenadas**  \n{r.get('start')}–{r.get('end')} · hebra "
                       f"{'−' if str(r.get('strand')) in ('-1','-') else '+'}")
            c.markdown(f"**Confianza**  \n{r.get('confidence','') or '—'}")
            st.caption(f"Genes: {r.get('gene_names') or r.get('genes','')}")
            st.caption(f"Evidencia: {r.get('evidence','') or '—'}")

    with t_prom:
        prom = tables["promoters"]
        if not prom.empty and not tables["operons"].empty:
            ids = set(tables["operons"].loc[
                (tables["operons"]["contig"] == contig)
                & (tables["operons"]["end"] >= start)
                & (tables["operons"]["start"] <= end), "operon_id"].astype(str))
            prom = prom[prom["operon_id"].astype(str).isin(ids)]
        i = _elegir(prom.reset_index(drop=True),
                    ["operon_id", "sigma", "box35", "box10", "spacer", "score", "confidence"],
                    "br_tab_prom")
        if i is not None:
            r = prom.reset_index(drop=True).iloc[i]
            st.markdown(f"### 🎯 Promotor de {r.get('operon_id')}")
            st.markdown(f"**Caja −35** `{r.get('box35','')}` (desajustes {r.get('mm35','')})  \n"
                        f"**Caja −10** `{r.get('box10','')}` (desajustes {r.get('mm10','')})  \n"
                        f"**Espaciador** {r.get('spacer','')} pb · **σ** {r.get('sigma','')}  \n"
                        f"**Score** {r.get('score','')} · **Confianza** {r.get('confidence','')}")
            st.caption("Un promotor es la marca donde la ARN polimerasa empieza a "
                       "copiar. La caja −10 y la −35 son las dos secuencias que "
                       "reconoce el factor sigma.")

    with t_hp:
        hp = tables["hairpins"]
        if not hp.empty and "contig" in hp:
            hp = hp[(hp["contig"] == contig)
                    & (hp["end"] >= start) & (hp["start"] <= end)]
        i = _elegir(hp.reset_index(drop=True),
                    ["start", "end", "class", "dG", "stem", "gc", "polyU",
                     "left_gene", "right_gene"],
                    "br_tab_hp")
        if i is not None:
            r = hp.reset_index(drop=True).iloc[i]
            cls = str(r.get("class", ""))
            st.markdown(f"### 🌀 Horquilla en {r.get('start')}–{r.get('end')}")
            st.markdown(f"**Clase** {cls} — {HAIRPIN_MEANING.get(cls, '')}  \n"
                        f"**ΔG** {r.get('dG','')} kcal/mol · **tallo** {r.get('stem','')} nt · "
                        f"**GC** {r.get('gc','')} · **poli-U** {r.get('polyU','')}  \n"
                        f"**Entre** {r.get('left_gene','')} y {r.get('right_gene','')}")
            st.caption("La clase compara el tallo con el fondo de GC del genoma. "
                       "A: estructura clara; B: intermedia, no se distingue del azar; "
                       "C: débil.")

    with t_susp:
        if susp.empty:
            st.caption("Sin sospechas de cadena metabólica en esta ventana.")
        else:
            i = _elegir(susp.reset_index(drop=True),
                        ["n_genes", "products", "reactions", "gap_mean",
                         "shared_metabolites", "confidence"],
                        "br_tab_susp")
            if i is not None:
                r = susp.reset_index(drop=True).iloc[i]
                st.markdown(f"### 🔗 Sospecha: {r.get('n_genes')} genes contiguos")
                st.markdown(f"**Genes**: {r.get('genes','')}  \n"
                            f"**Enzimas**: {r.get('products','')}  \n"
                            f"**Reacciones**: {r.get('reactions','')}  \n"
                            f"**Metabolitos compartidos**: {r.get('shared_metabolites','')}")
                st.caption("Estas enzimas tocan los mismos metabolitos, así que "
                           "probablemente participan de la misma ruta. Al estar "
                           "contiguas en el ADN, es **sospecha** de operón "
                           "metabólico. No es prueba: hay rutas contiguas que no "
                           "son operones y operones sin ruta asociada.")

    with t_tx:
        wt = tss_df[(tss_df.pos >= start) & (tss_df.pos <= end)] if not tss_df.empty else tss_df
        wtt = tts_df[(tts_df.pos >= start) & (tts_df.pos <= end)] if not tts_df.empty else tts_df
        if wt.empty and wtt.empty:
            st.caption("Sin TSS ni TTS anotados en esta ventana.")
        else:
            st.caption(
                "Marcas de **transcripción**: dónde empieza (TSS) y dónde termina "
                "(TTS) el ARN. Las sólidas están curadas; las atenuadas vienen de "
                "alto rendimiento. Es el dato que decide las fronteras cuando existe."
            )
            if not wt.empty:
                st.markdown(f"**TSS en la ventana ({len(wt)})**")
                st.dataframe(wt[["pos", "strand", "source", "detail"]],
                             width="stretch", hide_index=True)
            if not wtt.empty:
                st.markdown(f"**TTS en la ventana ({len(wtt)})**")
                st.dataframe(wtt[["pos", "strand", "source", "detail"]],
                             width="stretch", hide_index=True)

    st.info(
        "**Terminador ≠ atenuador.** Una horquilla al final de un operón lo "
        "termina. Una horquilla *dentro* del operón puede ser un atenuador: "
        "corta el transcrito pero el operón continúa. Por eso una horquilla no "
        "basta para decidir dónde acaba un operón."
    )
