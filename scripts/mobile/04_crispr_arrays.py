#!/usr/bin/env python3
"""04_crispr_arrays.py — detección de arrays CRISPR de novo sobre el FASTA.

No hay CRISPRCasFinder/CRT/minced instalados, así que se implementa un detector
determinista de repeticiones tipo CRISPR (mismo principio que CRT/minced):

  1. Ventanas solapadas de 12 kb (paso 6 kb) por contig.
  2. Se indexan los 12-mers de cada ventana y se conservan los que aparecen
     3–25 veces (semillas de repetición, no de bajo complejidad).
  3. Las posiciones de una semilla se agrupan en tandas con periodo 35–160 pb
     (repetición + espaciador).
  4. La repetición se extiende desde la semilla comparando las dos primeras
     copias; se exige longitud 23–50 pb y espaciador 20–80 pb.
  5. Se calcula identidad media entre copias (Hamming) y se descartan arrays
     que caen ≥80 % dentro de un CDS (repeticiones de proteína, no CRISPR) o de
     bajo complejidad.
  6. Se deduplican arrays detectados en ventanas contiguas.

Limitaciones (documentadas en el README): no detecta arrays con repeticiones
muy divergentes (<85 % de identidad) ni arrays de menos de 3 repeticiones; no
anota el tipo de sistema (eso lo hace 05_classify.py con los genes cas).

Salida: data/mobile_typing/_work/<genoma>/crispr_arrays.tsv
        data/mobile_typing/_work/<genoma>/crispr_spacers.tsv
        data/mobile_typing/_work/<genoma>/crispr_stats.json
"""
from __future__ import annotations

import csv
import gzip
import json
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data" / "raw"
OUT = ROOT / "data" / "mobile_typing"
WORK = OUT / "_work"

SEED = 12
WINDOW = 12_000
STEP = 6_000
MIN_REPEAT, MAX_REPEAT = 23, 50
MIN_SPACER, MAX_SPACER = 20, 80
MIN_IDENTITY = 0.85


def read_fasta(path: Path):
    name, chunks = None, []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    yield name, "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if name:
        yield name, "".join(chunks)


def load_cds(genome: str):
    """Intervalos CDS por contig (para descartar repeticiones dentro de genes)."""
    by_contig = defaultdict(list)
    path = WORK / genome / "cds.tsv"
    with path.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            by_contig[row["contig"]].append((int(row["start"]), int(row["end"])))
    for contig in by_contig:
        by_contig[contig].sort()
    return by_contig


def identity(a: str, b: str) -> float:
    n = min(len(a), len(b))
    if n == 0:
        return 0.0
    same = sum(1 for x, y in zip(a[:n], b[:n]) if x == y)
    return same / n


def span_repeat(seq: str, p0: int, p1: int):
    """Extiende la semilla a izquierda y derecha comparando las dos primeras copias.

    La semilla de 12 pb puede caer en cualquier posición dentro de la repetición;
    por eso se extiende en ambos sentidos hasta el primer mismatch (que marca el
    borde repetición/espaciador). Devuelve (offset_izquierdo, longitud_total,
    inicio_de_la_repetición).
    """
    left = 0
    while left < 30 and p0 - left - 1 >= 0 and seq[p0 - left - 1] == seq[p1 - left - 1]:
        left += 1
    right = 0
    while (right < 45 and p0 + SEED + right < len(seq)
           and p1 + SEED + right < len(seq)
           and seq[p0 + SEED + right] == seq[p1 + SEED + right]):
        right += 1
    length = left + SEED + right
    return left, length


def low_complexity(seq: str) -> bool:
    if not seq:
        return True
    counts = Counter(seq)
    if max(counts.values()) / len(seq) > 0.55:
        return True
    # homopolímero de 6 o más
    run = 1
    for a, b in zip(seq, seq[1:]):
        run = run + 1 if a == b else 1
        if run >= 6:
            return True
    return False


def find_arrays(contig: str, seq: str):
    found = []
    for w0 in range(0, len(seq), STEP):
        window = seq[w0 : w0 + WINDOW]
        if len(window) < 500:
            break
        positions = defaultdict(list)
        for i in range(len(window) - SEED + 1):
            positions[window[i : i + SEED]].append(i)
        consumed = []
        for seed, pos in sorted(positions.items(), key=lambda kv: kv[1][0]):
            if len(pos) < 3 or len(pos) > 500 or len(set(seed)) < 3 or low_complexity(seed):
                continue
            if any(pos[0] >= a and pos[0] <= b for a, b in consumed):
                continue
            # tandas con periodo 35–160 pb
            run = [pos[0]]
            for p in pos[1:]:
                if 35 <= p - run[-1] <= 160:
                    run.append(p)
                else:
                    if len(run) >= 3:
                        break
                    run = [p]
            if len(run) < 3:
                continue
            period = int(statistics.median([b - a for a, b in zip(run, run[1:])]))
            if not (35 <= period <= 160):
                continue
            p0, p1 = run[0], run[1]
            left, repeat_len = span_repeat(window, p0, p1)
            spacer_len = period - repeat_len
            if not (MIN_REPEAT <= repeat_len <= MAX_REPEAT and MIN_SPACER <= spacer_len <= MAX_SPACER):
                continue
            starts = [p - left for p in run]
            if any(x < 0 for x in starts):
                continue
            copies = [window[x : x + repeat_len] for x in starts]
            if low_complexity(copies[0]):
                continue
            ids = [identity(copies[0], c) for c in copies[1:]]
            mean_id = sum(ids) / len(ids)
            if mean_id < MIN_IDENTITY or min(ids) < 0.75:
                continue
            start = w0 + starts[0] + 1
            end = w0 + starts[-1] + repeat_len
            spacer_seqs = []
            spacer_coords = []
            for i in range(len(run) - 1):
                s0 = w0 + starts[i] + repeat_len + 1
                s1 = w0 + starts[i + 1]
                spacer_seqs.append(seq[s0 - 1 : s1])
                spacer_coords.append((s0, s1))
            if len(set(spacer_seqs)) < 2:
                continue
            # geometría consistente: cada espaciador debe tener longitud válida
            if any(s1 - s0 + 1 < MIN_SPACER for s0, s1 in spacer_coords):
                continue
            found.append(
                {
                    "contig": contig,
                    "start": start,
                    "end": end,
                    "n_repeats": len(run),
                    "interrupted": "no",
                    "repeat_len": repeat_len,
                    "spacer_len": spacer_len,
                    "period": period,
                    "mean_identity": round(mean_id, 3),
                    "repeat_unit": copies[0],
                    "spacers": spacer_seqs,
                    "spacer_coords": spacer_coords,
                }
            )
            consumed.append((starts[0], starts[-1] + repeat_len))
    return found


def overlap_frac(a_start, a_end, b_start, b_end) -> float:
    ov = min(a_end, b_end) - max(a_start, b_start) + 1
    return max(0.0, ov) / min(a_end - a_start + 1, b_end - b_start + 1)


def dedupe(arrays):
    arrays = sorted(arrays, key=lambda a: (a["start"], a["end"]))
    kept = []
    for arr in arrays:
        duplicate = False
        for old in kept:
            if overlap_frac(arr["start"], arr["end"], old["start"], old["end"]) >= 0.5:
                if arr["n_repeats"] > old["n_repeats"] or (
                    arr["n_repeats"] == old["n_repeats"] and arr["mean_identity"] > old["mean_identity"]
                ):
                    old.update(arr)
                duplicate = True
                break
        if not duplicate:
            kept.append(arr)
    return kept


def merge_adjacent(arrays, max_gap=2000, min_identity=0.9):
    """Une segmentos del mismo array separados por un hueco corto.

    Un array real puede quedar partido en dos por una repetición divergente o
    una inserción corta; ambos segmentos comparten la misma unidad de repetición
    y están a <2 kb. Se marca `interrupted=yes` cuando el hueco supera 200 pb,
    para que la evidencia diga explícitamente que el locus está interrumpido.
    """
    arrays = sorted(arrays, key=lambda a: (a["contig"], a["start"]))
    merged = []
    for arr in arrays:
        if merged:
            prev = merged[-1]
            gap = arr["start"] - prev["end"] - 1
            same_unit = identity(prev["repeat_unit"], arr["repeat_unit"]) >= min_identity
            if prev["contig"] == arr["contig"] and 0 <= gap <= max_gap and same_unit:
                if gap > 200:
                    prev["interrupted"] = "yes"
                prev["end"] = arr["end"]
                prev["n_repeats"] += arr["n_repeats"]
                prev["spacers"] += arr["spacers"]
                prev["spacer_coords"] += arr["spacer_coords"]
                prev["mean_identity"] = round(
                    (prev["mean_identity"] + arr["mean_identity"]) / 2, 3
                )
                prev["repeat_len"] = max(prev["repeat_len"], arr["repeat_len"])
                continue
        arr["interrupted"] = arr.get("interrupted", "no")
        merged.append(arr)
    return merged


def cds_overlap(intervals, start, end) -> float:
    best = 0.0
    for s, e in intervals:
        if s > end:
            break
        ov = min(end, e) - max(start, s) + 1
        if ov > 0:
            best = max(best, ov / (end - start + 1))
    return best


def detect_genome(genome, fna, cds_by_contig):
    """Devuelve (arrays, spacers_rows, stats) para un genoma.

    Reutilizado por el piloto (04) y por el barrido de 400 genomas (21).
    """
    arrays, spacers_rows = [], []
    filtered_cds = 0
    filtered_complexity = 0
    for contig, seq in read_fasta(fna):
        raw = find_arrays(contig, seq)
        raw = merge_adjacent(dedupe(raw))
        for arr in raw:
            if cds_overlap(cds_by_contig.get(contig, []), arr["start"], arr["end"]) >= 0.8:
                filtered_cds += 1
                continue
            if low_complexity(arr["repeat_unit"]):
                filtered_complexity += 1
                continue
            if arr["n_repeats"] >= 4 and arr["mean_identity"] >= 0.9:
                confidence = "high"
            elif arr["n_repeats"] == 3 or arr["mean_identity"] >= 0.85:
                confidence = "medium"
            else:
                confidence = "low"
            if arr["interrupted"] == "yes" and confidence == "high":
                confidence = "medium"
            arr["confidence"] = confidence
            arrays.append(arr)
    arrays.sort(key=lambda a: (a["contig"], a["start"]))
    for idx, arr in enumerate(arrays, start=1):
        arr["array_id"] = f"CRISPR{idx:03d}"
        for j, (spacer, (s0, s1)) in enumerate(zip(arr["spacers"], arr["spacer_coords"]), start=1):
            spacers_rows.append(
                {
                    "genome": genome,
                    "array_id": arr["array_id"],
                    "contig": arr["contig"],
                    "spacer_index": j,
                    "start": s0,
                    "end": s1,
                    "length": len(spacer),
                    "sequence": spacer,
                }
            )
    stats = {
        "arrays": len(arrays),
        "spacers": len(spacers_rows),
        "filtered_cds_overlap": filtered_cds,
        "filtered_low_complexity": filtered_complexity,
        "arrays_high": sum(1 for a in arrays if a["confidence"] == "high"),
        "arrays_medium": sum(1 for a in arrays if a["confidence"] == "medium"),
        "arrays_low": sum(1 for a in arrays if a["confidence"] == "low"),
    }
    return arrays, spacers_rows, stats


def write_detection(work, arrays, spacers_rows, stats):
    # escritura atómica: 05 puede estar leyendo mientras 21 produce
    tmp_arrays = work / "crispr_arrays.tsv.tmp"
    tmp_spacers = work / "crispr_spacers.tsv.tmp"
    with tmp_arrays.open("w", newline="") as fh:
        fields = ["array_id", "contig", "start", "end", "n_repeats", "repeat_len",
                  "spacer_len", "period", "mean_identity", "interrupted",
                  "confidence", "repeat_unit"]
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for arr in arrays:
            writer.writerow({k: arr[k] for k in fields})
    with tmp_spacers.open("w", newline="") as fh:
        fields = ["genome", "array_id", "contig", "spacer_index", "start", "end", "length", "sequence"]
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(spacers_rows)
    tmp_arrays.replace(work / "crispr_arrays.tsv")
    tmp_spacers.replace(work / "crispr_spacers.tsv")
    (work / "crispr_stats.json").write_text(json.dumps(stats, indent=2, sort_keys=True) + "\n")


def main() -> int:
    pilots = [row["genome"] for row in csv.DictReader((OUT / "pilot_genomes.tsv").open(), delimiter="\t")]
    for genome in pilots:
        work = WORK / genome
        cds = load_cds(genome)
        fna = RAW / f"{genome}_genomic.fna.gz"
        arrays, spacers_rows, stats = detect_genome(genome, fna, cds)
        write_detection(work, arrays, spacers_rows, stats)
        print(f"[{genome}] arrays={stats['arrays']} spacers={stats['spacers']} "
              f"(descartados: {stats['filtered_cds_overlap']} en CDS, "
              f"{stats['filtered_low_complexity']} baja complejidad)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
