"""Contexto de membrana/secuencia para el camino C (sin Swiss-Prot).

Solo aporta **hipótesis de contexto** calculadas con heurísticas declaradas:

- hidropatía de Kyte-Doolittle (1982), ventana de 19 aa, umbral de media 1.6,
  para contar segmentos hidrofóbicos tipo TM (``tm_segments``). **No es TMHMM**
  ni ningún predictor entrenado: es un conteo de ventanas hidrofóbicas.
- señal péptido N-terminal: primeros 30 aa, ventana hidrofóbica de 9 aa,
  motivo de corte Sec/SPI (residuos pequeños en -3/-1) y lipobox. **No es
  SignalP/Phobius**; los tres estados son ``yes``/``no``/``uncertain``.
- lipobox: regex ``[LVI][ASTM][AG]C`` en los primeros 40 aa (motivo de
  lipoproteína bacteriana).

No asigna función, no confirma localización ni complejo. No usa hmmscan,
DIAMOND ni red: es CPU pura sobre la secuencia de aminoácidos.
"""
from __future__ import annotations

import csv
import re
from pathlib import Path
from typing import Any

# Escala de hidropatía de Kyte & Doolittle (1982).
KD: dict[str, float] = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
    "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
    "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
    "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}

HYDRO_WINDOW = 19
HYDRO_THRESHOLD = 1.6
SIGNAL_HEAD = 30
SIGNAL_WINDOW = 9
SIGNAL_H_MAX = 1.8
SIGNAL_H_WEAK = 1.5
SIGNAL_H_START_MAX = 12
LIPOBOX_HEAD = 40
SMALL_RESIDUES = frozenset("ASTGCNQ")
RE_LIPOBOX = re.compile(r"[LVI][ASTM][AG]C")

CONTEXT_NOTE = (
    "Hipótesis de contexto por heurística de secuencia (Kyte-Doolittle "
    "ventana 19/umbral 1.6, N-terminal 30 aa, lipobox); no equivale a "
    "TMHMM/SignalP/Phobius, no asigna función ni confirma localización, "
    "topología o pertenencia a complejo."
)

CONTEXT_FIELDS = (
    "tm_segments", "signal_peptide", "lipobox",
    "context_hypothesis", "context_evidence", "context_note",
)


def _values(seq: str) -> list[float]:
    return [KD.get(aa, 0.0) for aa in seq.upper()]


def _max_window_mean(seq: str, window: int) -> float:
    vals = _values(seq)
    if len(vals) < window:
        return 0.0
    best = 0.0
    for i in range(len(vals) - window + 1):
        mean = sum(vals[i:i + window]) / window
        if mean > best:
            best = mean
    return best


def hydrophobic_segments(seq: str, window: int = HYDRO_WINDOW,
                         threshold: float = HYDRO_THRESHOLD
                         ) -> list[tuple[int, int, float]]:
    """Segmentos hidrofóbicos tipo Kyte-Doolittle (0-based, inclusivos).

    Desliza una ventana de ``window`` aa y marca las ventanas cuya media KD
    es >= ``threshold``. Ventanas contiguas/solapadas se fusionan en un solo
    segmento. Devuelve ``(inicio, fin, media_maxima)`` por segmento.

    Es una heurística de hidropatía, no un predictor de TM: no modela
    topología, carga del centro de la membrana ni longitud mínima de asa.
    """
    seq = seq.upper()
    vals = _values(seq)
    if len(vals) < window:
        return []
    runs: list[list[Any]] = []
    for i in range(len(vals) - window + 1):
        mean = sum(vals[i:i + window]) / window
        if mean < threshold:
            continue
        if runs and i <= runs[-1][1] + 1:
            runs[-1][1] = i + window - 1
            runs[-1][2] = max(runs[-1][2], mean)
        else:
            runs.append([i, i + window - 1, mean])
    return [(int(a), int(b), round(float(c), 3)) for a, b, c in runs]


def _sec_cleavage(seq: str, h_start: int, h_end: int) -> tuple[bool, int]:
    """Busca un motivo de corte Sec/SPI compatible con la región h.

    El motivo es un residuo pequeño en -3 y -1 (``[ASTGCNQ]``) dentro de los
    primeros 30 aa, con la ventana hidrofóbica inmediatamente antes.
    """
    upper = seq.upper()
    for cut in range(14, min(len(upper), SIGNAL_HEAD) + 1):
        motif = upper[cut - 3:cut]
        if len(motif) != 3:
            continue
        if (motif[0] in SMALL_RESIDUES and motif[2] in SMALL_RESIDUES
                and motif[1].isalpha()):
            if (h_start <= SIGNAL_H_START_MAX and h_start < cut
                    and h_end <= cut - 1 and h_end >= cut - 10):
                return True, cut
    return False, -1


def signal_peptide_call(seq: str) -> dict[str, Any]:
    """Heurística de señal péptido N-terminal (no es SignalP).

    Combina una ventana hidrofóbica de 9 aa en los primeros 30 aa con un
    motivo de corte Sec/SPI o con un lipobox (señal de lipoproteína).
    Devuelve un dict con ``status`` (yes/no/uncertain), ``hmax``, ``cleavage``,
    ``lipobox`` y ``why``.
    """
    head = seq.upper()[:SIGNAL_HEAD]
    lipobox = bool(RE_LIPOBOX.search(seq.upper()[:LIPOBOX_HEAD]))
    hmax = 0.0
    h_start = h_end = -1
    vals = _values(head)
    if len(vals) >= SIGNAL_WINDOW:
        last_start = min(len(vals) - SIGNAL_WINDOW, 21)
        for i in range(last_start + 1):
            mean = sum(vals[i:i + SIGNAL_WINDOW]) / SIGNAL_WINDOW
            if mean > hmax:
                hmax, h_start, h_end = mean, i, i + SIGNAL_WINDOW - 1
    cleavage, cut = _sec_cleavage(head, h_start, h_end)
    if lipobox and hmax >= HYDRO_THRESHOLD:
        status, why = "yes", "lipobox + region hidrofobica N-terminal"
    elif hmax >= SIGNAL_H_MAX and cleavage:
        status, why = "yes", "region h + motivo de corte Sec/SPI"
    elif hmax >= SIGNAL_H_MAX:
        status, why = "uncertain", "region h sin motivo de corte reconocido"
    elif hmax >= SIGNAL_H_WEAK and cleavage:
        status, why = "uncertain", "hidrofobicidad debil con posible corte"
    else:
        status, why = "no", "sin region h ni motivo de corte"
    return {
        "status": status,
        "hmax": round(hmax, 3),
        "cleavage": cleavage,
        "cut": cut,
        "lipobox": lipobox,
        "why": why,
    }


def context_call(seq: str) -> dict[str, Any]:
    """Calcula los campos de contexto de una secuencia (sin mutar nada).

    Una secuencia vacía no se interpreta como soluble: queda en ``unknown``.
    """
    seq = (seq or "").upper()
    if not seq:
        return {
            "tm_segments": 0,
            "signal_peptide": "uncertain",
            "lipobox": "no",
            "context_hypothesis": "unknown",
            "context_evidence": "sin_secuencia_para_la_proteina_C",
            "context_note": CONTEXT_NOTE,
        }
    segs = hydrophobic_segments(seq)
    sig = signal_peptide_call(seq)
    n_term_only = all(end <= 40 for _, end, _ in segs)
    kd19_max = _max_window_mean(seq, HYDRO_WINDOW)

    if sig["status"] == "yes" and sig["lipobox"]:
        hypothesis = "membrane_like"
        reason = "lipobox_ancla_lipoproteina"
    elif sig["status"] == "yes" and not n_term_only:
        hypothesis = "membrane_like"
        reason = "senal_exportacion_mas_segmento_hidrofobico_adicional"
    elif sig["status"] == "yes":
        hypothesis = "secreted_like"
        reason = "senal_exportacion_sin_segmento_transmembrana"
    elif segs and not n_term_only:
        hypothesis = "membrane_like"
        reason = "segmento_hidrofobico_no_N_terminal"
    elif segs:
        hypothesis = "unknown"
        reason = "solo_segmento_hidrofobico_N_terminal_ambiguo"
    elif sig["status"] == "uncertain":
        hypothesis = "unknown"
        reason = "senal_N_terminal_incierta"
    else:
        hypothesis = "soluble_like"
        reason = "sin_segmentos_hidrofobicos_ni_senal"

    seg_txt = "-" if not segs else ",".join(
        f"{start + 1}-{end + 1}" for start, end, _ in segs)
    evidence = (
        f"tm_segments={len(segs)} ({seg_txt}); "
        f"signal_peptide={sig['status']} (h9_max={sig['hmax']:.2f}, "
        f"cleavage={'si' if sig['cleavage'] else 'no'}); "
        f"lipobox={'si' if sig['lipobox'] else 'no'}; "
        f"kd19_max={kd19_max:.2f}; "
        f"hipotesis={hypothesis} ({reason})"
    )
    return {
        "tm_segments": len(segs),
        "signal_peptide": sig["status"],
        "lipobox": "yes" if sig["lipobox"] else "no",
        "context_hypothesis": hypothesis,
        "context_evidence": evidence,
        "context_note": CONTEXT_NOTE,
        "_reason": reason,
    }


def run(rows: list[dict], proteins: list[dict], outdir: Path, genome: str,
        cpu: int = 4) -> dict[str, Any]:
    """Añade contexto a las filas del camino C y escribe ``.context_c.tsv``.

    Solo procesa filas con ``path == "C"``; las filas A/B quedan intactas.
    Escribe un TSV con genome, protein_id y los campos de contexto.
    ``cpu`` se acepta por simetría con las otras capas; el cálculo es una
    heurística lineal por proteína y se ejecuta en serie (no lanza procesos).
    """
    try:
        outdir = Path(outdir)
        outdir.mkdir(parents=True, exist_ok=True)
    except Exception as exc:  # noqa: BLE001
        return {"ok": False, "reason": f"no se pudo crear {outdir}: {exc}",
                "n": 0, "rows": []}

    seqs = {str(p.get("id")): str(p.get("aa", ""))
            for p in (proteins or []) if p.get("id") is not None}
    out_rows: list[dict[str, Any]] = []
    for row in rows or []:
        if row.get("path") != "C":
            continue
        pid = str(row.get("protein_id", ""))
        seq = seqs.get(pid, "")
        fields = context_call(seq)
        fields.pop("_reason", None)
        row.update(fields)
        out_rows.append({"genome": genome, "protein_id": pid, **fields})

    path = outdir / f"{genome}.context_c.tsv"
    cols = ["genome", "protein_id", *CONTEXT_FIELDS]
    try:
        with path.open("w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, delimiter="\t", fieldnames=cols,
                                    extrasaction="ignore")
            writer.writeheader()
            writer.writerows(out_rows)
    except OSError as exc:
        return {"ok": False, "reason": f"no se pudo escribir {path}: {exc}",
                "n": len(out_rows), "rows": out_rows}
    return {"ok": True, "reason": "", "n": len(out_rows),
            "rows": out_rows, "tsv": str(path)}
