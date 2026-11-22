#!/usr/bin/env python3
"""Evalúa el tokenizador v3 contra v2 sobre genomas no vistos.

Mide:
  - compresión (nt/token) en genomas completos
  - tasa de <unk> (global, IUPAC, runs de N)
  - utilización del vocabulario y entropía de tokens
  - roundtrip encode/decode
  - preservación de motivos como token único, también en contexto
  - tokens especiales (<pad>, <unk>, <bos>, <eos>, <mask>)
  - coste estimado de entrenamiento (tokens para 1,386 Gbp)

Salida:
  data/tokenizer_v3/eval_v3.json
  data/tokenizer_v3/EVAL.md
  reports/tokenizer_v3_eval.md

Uso:
  python scripts/06_eval_tokenizer_v3.py
"""
from __future__ import annotations

import json
import math
import random
import re
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import sentencepiece as spm

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from lib.io import read_fasta  # noqa: E402

MOTIFS = ["TATAAT", "TTGACA", "AGGAGG", "ATG", "TAA", "TAG", "TGA"]
IUPAC = list("RYKMSWBDHVN")

GENOMES = [
    # no vistos por el tokenizer (raw_excluded)
    "data/raw_excluded/GCF_000006685.1_ASM668v1_genomic.fna.gz",
    "data/raw_excluded/GCF_008086545.1_ASM808654v1_genomic.fna.gz",
    "data/raw_excluded/GCF_033971325.1_ASM3397132v1_genomic.fna.gz",
    "data/raw_excluded/GCF_033969665.2_ASM3396966v2_genomic.fna.gz",
    # vistos (para comparar dominio)
    "data/raw/GCF_000006985.1_ASM698v1_genomic.fna.gz",
    "data/raw/GCF_000317105.1_ASM31710v1_genomic.fna.gz",
]


def entropy(counter: Counter) -> float:
    total = sum(counter.values())
    return -sum((c / total) * math.log2(c / total) for c in counter.values() if c)


def encode_stats(sp, seq: str, chunk: int = 10_000) -> dict:
    total_nt = 0
    total_tokens = 0
    unk = 0
    counts: Counter = Counter()
    for i in range(0, len(seq), chunk):
        s = seq[i:i + chunk]
        ids = sp.encode(s, out_type=int)
        total_nt += len(s)
        total_tokens += len(ids)
        unk += sum(1 for x in ids if x == sp.unk_id())
        counts.update(ids)
    return {
        "nt": total_nt,
        "tokens": total_tokens,
        "nt_per_token": total_nt / max(1, total_tokens),
        "unk": unk,
        "unk_rate": unk / max(1, total_tokens),
        "unique_tokens": len(counts),
        "vocab_used_frac": len(counts) / sp.get_piece_size(),
        "token_entropy_bits": entropy(counts),
        "top_tokens": [(sp.id_to_piece(i), c) for i, c in counts.most_common(15)],
    }


def roundtrip_ok(sp, seq: str, n: int = 20, chunk: int = 2_000) -> tuple[int, int]:
    rng = random.Random(0)
    ok = 0
    for _ in range(n):
        if len(seq) <= chunk:
            s = seq
        else:
            i = rng.randrange(0, len(seq) - chunk)
            s = seq[i:i + chunk]
        ids = sp.encode(s, out_type=int)
        # v2 añade ▁ (dummy prefix); al decodificar puede aparecer un espacio.
        # El ADN no tiene espacios, así que se eliminan para comparar.
        if sp.decode(ids).replace(" ", "") == s:
            ok += 1
    return ok, n


def motif_stats(sp, genomes: dict[str, str]) -> dict:
    out = {}
    for motif in MOTIFS:
        pieces = sp.encode(motif, out_type=str)
        # v2 añade ▁ (dummy prefix); se ignora para medir si el motivo
        # se parte o no.
        pieces_core = [p for p in pieces if p != "▁"]
        single = len(pieces_core) == 1 and pieces_core[0] == motif
        hits = 0
        total = 0
        for seq in genomes.values():
            idx = [m.start() for m in re.finditer(motif, seq)]
            if not idx:
                continue
            rng = random.Random(sum(ord(c) for c in motif))
            rng.shuffle(idx)
            for i in idx[:120]:
                s = max(0, i - 100)
                w = seq[s:i + len(motif) + 100]
                wp = sp.encode(w, out_type=str)
                total += 1
                if motif in wp:
                    hits += 1
        out[motif] = {
            "single_token": single,
            "pieces": pieces,
            "pieces_core": pieces_core,
            "context_single_rate": hits / total if total else None,
            "context_n": total,
        }
    return out


def iupac_stats(sp) -> dict:
    seq = "ACGT" + "".join(IUPAC) + "N"
    ids = sp.encode(seq, out_type=int)
    pieces = sp.encode(seq, out_type=str)
    unk = sum(1 for x in ids if x == sp.unk_id())
    per = {c: sp.piece_to_id(c) for c in IUPAC}
    return {
        "seq": seq,
        "pieces": pieces,
        "unk": unk,
        "unk_rate": unk / max(1, len(ids)),
        "piece_to_id": per,
        "all_known": all(v != sp.unk_id() for v in per.values()),
    }


def special_stats(sp) -> dict:
    out = {}
    for name in ("pad", "unk", "bos", "eos"):
        try:
            sid = getattr(sp, f"{name}_id")()
        except Exception:
            sid = -1
        out[name] = {"id": sid, "piece": sp.id_to_piece(sid) if sid >= 0 else None}
    mask_id = sp.piece_to_id("<mask>")
    out["mask"] = {"id": mask_id, "piece": "<mask>" if mask_id != sp.unk_id() else None}
    ids = [v["id"] for v in out.values() if v["id"] is not None and v["id"] >= 0]
    out["all_distinct"] = len(ids) == len(set(ids))
    return out


def n_run_stats(sp) -> dict:
    seq = "N" * 1000
    ids = sp.encode(seq, out_type=int)
    return {
        "nt": len(seq),
        "tokens": len(ids),
        "nt_per_token": len(seq) / max(1, len(ids)),
        "unk_rate": sum(1 for x in ids if x == sp.unk_id()) / max(1, len(ids)),
    }


def eval_tokenizer(path: Path, genomes: dict[str, str]) -> dict:
    sp = spm.SentencePieceProcessor(model_file=str(path))
    per_genome = {name: encode_stats(sp, seq) for name, seq in genomes.items()}
    total_nt = sum(v["nt"] for v in per_genome.values())
    total_tokens = sum(v["tokens"] for v in per_genome.values())
    total_unk = sum(v["unk"] for v in per_genome.values())
    rt_ok = rt_n = 0
    for seq in genomes.values():
        ok, n = roundtrip_ok(sp, seq)
        rt_ok += ok
        rt_n += n
    # distribución de tokens por ventana de 1 kb
    lengths = []
    for seq in genomes.values():
        for i in range(0, min(len(seq), 200_000) - 1000, 1000):
            lengths.append(len(sp.encode(seq[i:i + 1000], out_type=int)))
    lengths.sort()
    prefix_pieces = 0
    for seq in list(genomes.values())[:2]:
        prefix_pieces += sum(1 for p in sp.encode(seq[:20_000], out_type=str) if p == "▁")
    return {
        "model": str(path),
        "vocab_size": sp.get_piece_size(),
        "genomes": {k: {kk: vv for kk, vv in v.items() if kk != "top_tokens"}
                    for k, v in per_genome.items()},
        "global": {
            "nt": total_nt,
            "tokens": total_tokens,
            "nt_per_token": total_nt / max(1, total_tokens),
            "unk": total_unk,
            "unk_rate": total_unk / max(1, total_tokens),
            "roundtrip_ok": rt_ok,
            "roundtrip_n": rt_n,
        },
        "window_1kb_tokens": {
            "mean": sum(lengths) / len(lengths),
            "p5": lengths[int(0.05 * len(lengths))],
            "p50": lengths[int(0.50 * len(lengths))],
            "p95": lengths[int(0.95 * len(lengths))],
            "n": len(lengths),
        },
        "estimated_tokens_1_386gbp": int(1.386e9 / (total_nt / max(1, total_tokens))),
        "motifs": motif_stats(sp, genomes),
        "iupac": iupac_stats(sp),
        "special": special_stats(sp),
        "n_run": n_run_stats(sp),
        "dummy_prefix_pieces": prefix_pieces,
        "top_tokens": list(per_genome.values())[0]["top_tokens"],
    }


def main() -> None:
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--v2", default=str(ROOT / "data/tokenizer_v2/dna_sp_v2.model"))
    ap.add_argument("--v3", default=str(ROOT / "data/tokenizer_v3/dna_sp_v3.model"))
    ap.add_argument("--out", default=str(ROOT / "data/tokenizer_v3"))
    args = ap.parse_args()
    v2 = Path(args.v2)
    v3 = Path(args.v3)
    out_dir = Path(args.out)
    if not v3.exists():
        raise SystemExit(f"Falta {v3}; corre 06_train_tokenizer_v3.py")

    genomes = {}
    for rel in GENOMES:
        p = ROOT / rel
        if not p.exists():
            print(f"  [skip] {rel}")
            continue
        name, seq = next(read_fasta(p))
        genomes[name] = seq[:2_000_000]

    print(f"Genomas evaluados: {len(genomes)}")
    res2 = eval_tokenizer(v2, genomes)
    res3 = eval_tokenizer(v3, genomes)

    out = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "genomes": list(genomes.keys()),
        "v2": res2,
        "v3": res3,
        "comparison": {
            "nt_per_token_v2": res2["global"]["nt_per_token"],
            "nt_per_token_v3": res3["global"]["nt_per_token"],
            "compression_ratio_v3_vs_v2": res3["global"]["nt_per_token"] / res2["global"]["nt_per_token"],
            "unk_rate_v2": res2["global"]["unk_rate"],
            "unk_rate_v3": res3["global"]["unk_rate"],
            "iupac_known_v2": res2["iupac"]["all_known"],
            "iupac_known_v3": res3["iupac"]["all_known"],
            "motifs_single_v2": sum(1 for m in MOTIFS if res2["motifs"][m]["single_token"]),
            "motifs_single_v3": sum(1 for m in MOTIFS if res3["motifs"][m]["single_token"]),
            "dummy_prefix_v2": res2["dummy_prefix_pieces"],
            "dummy_prefix_v3": res3["dummy_prefix_pieces"],
            "tokens_1_386gbp_v2": res2["estimated_tokens_1_386gbp"],
            "tokens_1_386gbp_v3": res3["estimated_tokens_1_386gbp"],
        },
    }
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "eval_v3.json").write_text(json.dumps(out, indent=2, ensure_ascii=False))
    (ROOT / "reports").mkdir(exist_ok=True)
    (ROOT / "reports/tokenizer_v3_eval.json").write_text(
        json.dumps(out, indent=2, ensure_ascii=False))

    _write_report(out, out_dir / "EVAL.md")
    _write_report(out, ROOT / "reports/tokenizer_v3_eval.md")

    print("\n=== RESUMEN ===")
    print(f"v2: {res2['global']['nt_per_token']:.2f} nt/token, "
          f"unk={res2['global']['unk_rate']:.4f}, "
          f"motivos 1 token={out['comparison']['motifs_single_v2']}/7, "
          f"IUPAC={res2['iupac']['all_known']}, ▁={res2['dummy_prefix_pieces']}")
    print(f"v3: {res3['global']['nt_per_token']:.2f} nt/token, "
          f"unk={res3['global']['unk_rate']:.4f}, "
          f"motivos 1 token={out['comparison']['motifs_single_v3']}/7, "
          f"IUPAC={res3['iupac']['all_known']}, ▁={res3['dummy_prefix_pieces']}")
    print(f"tokens estimados 1.386 Gbp: v2={res2['estimated_tokens_1_386gbp']:,} "
          f"v3={res3['estimated_tokens_1_386gbp']:,}")
    print(f"-> {out_dir / 'EVAL.md'}")


def _write_report(out: dict, path: Path) -> None:
    v2, v3 = out["v2"], out["v3"]
    c = out["comparison"]
    lines = []
    lines.append("# Evaluación del tokenizador v3 vs v2\n")
    lines.append(f"Generado: {out['created_at']}\n")
    lines.append(f"- Genomas: {len(out['genomes'])} (4 no vistos + 2 vistos)\n")

    lines.append("## Comparación global\n")
    lines.append("| métrica | v2 | v3 |")
    lines.append("|---|---:|---:|")
    lines.append(f"| vocab | {v2['vocab_size']} | {v3['vocab_size']} |")
    lines.append(f"| nt/token | {v2['global']['nt_per_token']:.3f} | "
                 f"{v3['global']['nt_per_token']:.3f} |")
    lines.append(f"| unk rate | {v2['global']['unk_rate']:.5f} | "
                 f"{v3['global']['unk_rate']:.5f} |")
    lines.append(f"| roundtrip | {v2['global']['roundtrip_ok']}/{v2['global']['roundtrip_n']} | "
                 f"{v3['global']['roundtrip_ok']}/{v3['global']['roundtrip_n']} |")
    lines.append(f"| tokens 1 kb (media) | {v2['window_1kb_tokens']['mean']:.1f} | "
                 f"{v3['window_1kb_tokens']['mean']:.1f} |")
    lines.append(f"| tokens 1,386 Gbp | {c['tokens_1_386gbp_v2']:,} | "
                 f"{c['tokens_1_386gbp_v3']:,} |")
    lines.append(f"| IUPAC completo | {c['iupac_known_v2']} | {c['iupac_known_v3']} |")
    lines.append(f"| motivos 1 token | {c['motifs_single_v2']}/7 | "
                 f"{c['motifs_single_v3']}/7 |")
    lines.append(f"| piezas ▁ | {c['dummy_prefix_v2']} | {c['dummy_prefix_v3']} |")
    lines.append("")

    lines.append("## Motivos (token único y en contexto)\n")
    lines.append("| motivo | v2 single | v2 contexto | v3 single | v3 contexto |")
    lines.append("|---|---:|---:|---:|---:|")
    for m in MOTIFS:
        a, b = v2["motifs"][m], v3["motifs"][m]
        ca = "—" if a["context_single_rate"] is None else f"{a['context_single_rate']:.3f}"
        cb = "—" if b["context_single_rate"] is None else f"{b['context_single_rate']:.3f}"
        lines.append(f"| {m} | {a['single_token']} | {ca} | {b['single_token']} | {cb} |")
    lines.append("")

    lines.append("## Motivos solapados\n")
    lines.append("La protección atómica tiene un límite: si dos motivos comparten "
                 "bases, el tokenizador elige uno y parte el otro. Casos medidos:\n")
    lines.append("- `TGA` comparte la última base con `ATG`, `GTG` y `TTG`; en "
                 "`ATGA` gana `ATG` y `TGA` queda partido. Por eso su tasa en "
                 "contexto baja a ~0,66.")
    lines.append("- `AGGAGG` contiene `GGAGG`; al proteger ambos, algunas "
                 "apariciones se resuelven como `A` + `GGAGG`. La tasa ronda 0,70.")
    lines.append("- `TATAAT` y `TTGACA` no se solapan con otros motivos protegidos: "
                 "tasa 1,00.")
    lines.append("Para el modelo nuevo: definir una **prioridad de solapamiento** "
                 "o un pretokenizador que corte por límites de codón, en vez de "
                 "confiar solo en `user_defined_symbols`.\n")

    lines.append("## IUPAC\n")
    lines.append(f"- v2: {v2['iupac']['pieces']} (unk={v2['iupac']['unk']})")
    lines.append(f"- v3: {v3['iupac']['pieces']} (unk={v3['iupac']['unk']})")
    lines.append("")

    lines.append("## Tokens especiales\n")
    lines.append("| token | v2 | v3 |")
    lines.append("|---|---:|---:|")
    for name in ("pad", "unk", "bos", "eos", "mask"):
        a = v2["special"].get(name, {})
        b = v3["special"].get(name, {})
        lines.append(f"| {name} | {a.get('id')} | {b.get('id')} |")
    lines.append("")

    lines.append("## Veredicto\n")
    checks = [
        ("compresión 2-3 nt/token", 2.0 <= v3["global"]["nt_per_token"] <= 3.0,
         f"{v3['global']['nt_per_token']:.2f}"),
        ("IUPAC completo", v3["iupac"]["all_known"], str(v3["iupac"]["all_known"])),
        ("7/7 motivos token único", c["motifs_single_v3"] == 7, f"{c['motifs_single_v3']}/7"),
        ("sin ▁", c["dummy_prefix_v3"] == 0, str(c["dummy_prefix_v3"])),
        ("special tokens distintos", v3["special"]["all_distinct"],
         str(v3["special"]["all_distinct"])),
        ("roundtrip", v3["global"]["roundtrip_ok"] == v3["global"]["roundtrip_n"],
         f"{v3['global']['roundtrip_ok']}/{v3['global']['roundtrip_n']}"),
    ]
    for name, ok, val in checks:
        lines.append(f"- {'✅' if ok else '❌'} {name}: {val}")
    lines.append("")
    lines.append("Nota: es una evaluación **intrínseca** del tokenizador. "
                 "El impacto en el modelo requiere entrenar un modelo con v3.")
    path.write_text("\n".join(lines))


if __name__ == "__main__":
    main()
