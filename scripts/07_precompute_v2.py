#!/usr/bin/env python3
"""Precomputa nt/bpe/motif por secuencia usando el tokenizador v2."""
import sys, re, json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from lib.io import (read_config, read_fasta, genome_id, setup_logging,
                    manifest_append, manifest_load)
import numpy as np
import sentencepiece as spm
from multiprocessing import Pool

_CFG = _SP = _NT = _USER = _PHY = None

# Mapeo de IUPAC a IDs numéricos (para el array 'nt')
IUPAC_MAP = {
    "A": 0, "C": 1, "G": 2, "T": 3,
    "R": 4, "Y": 5, "K": 6, "M": 7,
    "S": 8, "W": 9, "B": 10, "D": 11,
    "H": 12, "V": 13, "N": 14
}

def _init_worker(cfg, model_path, user_symbols, phylum_map):
    global _CFG, _SP, _NT, _USER, _PHY
    _CFG = cfg
    _SP = spm.SentencePieceProcessor(model_file=model_path)
    _NT = IUPAC_MAP
    _USER = user_symbols
    _PHY = phylum_map

def classify_sequence(header, seq):
    h = header.lower()
    if "plasmid" in h: return "plasmid"
    if "chromosome" in h or "complete genome" in h: return "chromosome"
    return "unknown"

def encode_chunked(sp, seq, chunk=2_000):
    """Tokeniza en chunks para evitar el límite de SentencePiece."""
    pieces, ids = [], []
    for i in range(0, len(seq), chunk):
        sub = seq[i:i+chunk]
        p = sp.encode(sub, out_type=str)
        if p and p[0] == "▁":
            p = p[1:]
        pieces.extend(p)
        for x in p:
            ids.append(sp.piece_to_id(x))
    return pieces, ids

def _process_one(fasta_path):
    gid = genome_id(fasta_path)
    out_dir = Path(_CFG["project_root"]) / "data" / "processed_v2"
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

        # 1. nt_ids con soporte para IUPAC (usando el nuevo mapeo)
        nt_ids = np.fromiter((_NT.get(c, 14) for c in seq), dtype=np.uint8, count=len(seq))
        
        # 2. BPE v2
        pieces, bpe_ids_list = encode_chunked(_SP, seq)
        bpe_ids = np.array(bpe_ids_list, dtype=np.uint16)
        bpe_lens = np.array([len(p) for p in pieces], dtype=np.uint8)
        bpe_starts = np.concatenate([[0], np.cumsum(bpe_lens)[:-1]]).astype(np.uint32)
            
        # 3. Motif mask (buscando motivos protegidos)
        motif_mask = np.zeros(len(seq), dtype=np.uint8)
        for m in _USER:
            for match in re.finditer(re.escape(m), seq):
                motif_mask[match.start():match.end()] = 1

        # 4. Guardar
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
        }))
    return results

def main():
    cfg = read_config()
    log = setup_logging("07_precompute_v2", cfg["logging"]["dir"])
    root = Path(cfg["project_root"])
    processed = root / "data" / "processed_v2"
    processed.mkdir(parents=True, exist_ok=True)
    manifest = processed / "manifest.jsonl"
    model = root / "data" / "tokenizer_v2" / "dna_sp_v2.model"
    
    if not model.exists():
        print("❌ No se encontró el tokenizador v2. Corre 06_train_tokenizer_v2.py primero.")
        return

    # Cargar phylum_map desde el archivo existente
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

    fastas = sorted((root / "data" / "raw").glob(cfg["data"]["fasta_pattern"]))
    
    # Verificar cuáles ya tienen .npz v2
    pending = []
    for f in fastas:
        gid = genome_id(f)
        expected_npz = list(processed.glob(f"{gid}__*.npz"))
        if not expected_npz:
            pending.append(f)

    print(f"Genomas totales: {len(fastas)} | Pendientes de v2: {len(pending)}")
    if not pending: return

    # Obtener los user_symbols para el motif_mask
    meta_path = root / "data" / "tokenizer_v2" / "tokenizer_meta_v2.json"
    user_symbols = []
    if meta_path.exists():
        with open(meta_path) as f:
            user_symbols = json.load(f).get("user_symbols", [])

    with Pool(cfg["precompute"]["n_workers"], initializer=_init_worker,
              initargs=(cfg, str(model), user_symbols, phylum_map)) as pool:
        for i, seqs in enumerate(pool.imap_unordered(_process_one, pending, chunksize=1)):
            for seq_id, status, info in seqs:
                manifest_append(manifest, {"seq_id": seq_id, "status": status, "info": info})
            if (i+1) % 20 == 0:
                print(f"  [{i+1}/{len(pending)}] genomas procesados")
                log.info(f"[{i+1}/{len(pending)}]")

if __name__ == "__main__":
    main()
