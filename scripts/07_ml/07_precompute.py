#!/usr/bin/env python3
"""Precomputa nt/bpe/spans/motif por secuencia. Un .npz por secuencia."""
import sys, re, json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # scripts/ (lib)
sys.path.insert(0, str(Path(__file__).resolve().parent))        # esta carpeta
from lib.io import (read_config, read_fasta, genome_id, setup_logging,
                    manifest_append, manifest_load)
import numpy as np
import sentencepiece as spm
from multiprocessing import Pool

_CFG = _SP = _NT = _USER = _PHY = None

def _init_worker(cfg, model_path, user_symbols, phylum_map):
    global _CFG, _SP, _NT, _USER, _PHY
    _CFG = cfg
    _SP = spm.SentencePieceProcessor(model_file=model_path)
    _NT = {"A": 0, "C": 1, "G": 2, "T": 3}
    _USER = user_symbols
    _PHY = phylum_map

def classify_sequence(header, seq):
    h = header.lower()
    if "plasmid" in h: return "plasmid"
    if "chromosome" in h or "complete genome" in h: return "chromosome"
    return "unknown"

def encode_chunked(sp, seq, chunk=2_000):
    """Trocea la secuencia en bloques de 2048 para respetar el limite de SP."""
    pieces, ids = [], []
    for i in range(0, len(seq), chunk):
        sub = seq[i:i+chunk]
        p = sp.encode(sub, out_type=str)
        # Ignorar el token de espacio inicial (▁) para que la longitud coincida
        if p and p[0] == "▁":
            p = p[1:]
        pieces.extend(p)
        # Convertir los tokens de texto a IDs numéricos
        for x in p:
            ids.append(sp.piece_to_id(x))
    return pieces, ids

def _process_one(fasta_path):
    gid = genome_id(fasta_path)
    out_dir = Path(_CFG["project_root"]) / _CFG["data"]["processed_dir"]
    min_len = _CFG["precompute"]["min_seq_length"]
    max_n = _CFG["precompute"]["max_n_fraction"]
    results = []

    for name, seq in read_fasta(fasta_path):
        if len(seq) < min_len: continue
        if seq.count("N") / max(len(seq),1) > max_n: continue
        
        seq_type = classify_sequence(name, seq)
        accession = name.split(".")[0] if "." in name else name
        seq_id = f"{gid}__{accession}"
        out_path = out_dir / f"{seq_id}.npz"
        
        if out_path.exists():
            results.append((seq_id, "skipped", {"seq_type": seq_type}))
            continue

        # 1. nt_ids
        nt_ids = np.fromiter((_NT.get(c,4) for c in seq), dtype=np.uint8, count=len(seq))
        
        # 2. BPE (usando la funcion troceada corregida)
        pieces, bpe_ids_list = encode_chunked(_SP, seq)
        bpe_ids = np.array(bpe_ids_list, dtype=np.uint16)
        bpe_lens = np.array([len(p) for p in pieces], dtype=np.uint8)
        bpe_starts = np.concatenate([[0], np.cumsum(bpe_lens)[:-1]]).astype(np.uint32)
            
        # 3. Motif mask
        motif_mask = np.zeros(len(seq), dtype=np.uint8)
        for m in _USER:
            for match in re.finditer(re.escape(m), seq):
                motif_mask[match.start():match.end()] = 1

        # 4. Guardar usando un file handle (Evita el bug de la extension .npz)
        tmp = out_path.with_suffix(".npz.tmp")
        try:
            with open(tmp, 'wb') as f:
                np.savez_compressed(f, nt=nt_ids, bpe=bpe_ids,
                                    bpe_lens=bpe_lens, bpe_starts=bpe_starts,
                                    motif=motif_mask)
            tmp.rename(out_path)
        except Exception as e:
            results.append((seq_id, "fail", {"reason": str(e)}))
            if tmp.exists(): tmp.unlink()
            continue

        results.append((seq_id, "ok", {
            "genome_id": gid, "accession": accession,
            "seq_type": seq_type, "phylum": _PHY.get(gid, "unknown"),
            "length_nt": len(seq), "n_bpe": len(bpe_ids),
            "compression_ratio": len(seq)/max(len(bpe_ids),1),
            "motif_fraction": float(motif_mask.mean()),
        }))
    return results

def main():
    cfg = read_config()
    log = setup_logging("07_precompute", cfg["logging"]["dir"])
    root = Path(cfg["project_root"])
    processed = root / cfg["data"]["processed_dir"]
    processed.mkdir(parents=True, exist_ok=True)
    manifest = processed / "manifest.jsonl"
    model = root / cfg["data"]["tokenizer_dir"] / "dna_sp.model"

    # phylum_map
    from lib.io import genome_id as gid_fn
    phylum_map = {}
    sel = root / "data" / "selected.tsv"
    if sel.exists():
        with open(sel) as f:
            f.readline()
            for line in f:
                p = line.rstrip("\n").split("\t")
                if len(p) >= 5:
                    fname = p[3].rstrip("/").split("/")[-1] + "_genomic.fna.gz"
                    phylum_map[gid_fn(fname)] = p[4]

    fastas = sorted((root / cfg["data"]["raw_dir"]).glob(cfg["data"]["fasta_pattern"]))
    done = manifest_load(manifest)
    # Filtramos los que ya existen fisicamente
    pending = []
    for f in fastas:
        gid = genome_id(f)
        expected_npz = list(processed.glob(f"{gid}__*.npz"))
        if not expected_npz:
            pending.append(f)

    log.info(f"Genomas: {len(fastas)} total, {len(fastas)-len(pending)} con .npz, {len(pending)} pendientes")
    print(f"Pendientes: {len(pending)}/{len(fastas)}")
    if not pending: return

    with Pool(cfg["precompute"]["n_workers"], initializer=_init_worker,
              initargs=(cfg, str(model), cfg["tokenizer"]["user_symbols"], phylum_map)) as pool:
        for i, seqs in enumerate(pool.imap_unordered(_process_one, pending, chunksize=1)):
            for seq_id, status, info in seqs:
                manifest_append(manifest, {"seq_id": seq_id, "status": status, "info": info})
            if (i+1) % 20 == 0:
                print(f"  [{i+1}/{len(pending)}] genomas procesados")
                log.info(f"[{i+1}/{len(pending)}]")

if __name__ == "__main__":
    main()
