"""Panel HMM de confirmación para cualquier genoma.

Se ejecuta después de detectar brechas: si hay brechas de familias dudosas
(cGAMP, hm-dCMP, sulfolípidos), confirma si la enzima está presente aunque no
esté mapeada a Rhea.
"""
from __future__ import annotations

import re
import subprocess
import tempfile
from pathlib import Path

MODEL_FUNCTION = {
    "DncV-like_NTFase": "cGAS/DncV (cGAMP)",
    "DncV_C": "cGAS/DncV (cGAMP)",
    "ThyX_C": "ThyX (timidilato sintasa)",
    "ThyX_C_Alphaprot": "ThyX",
    "Thymidylat_synt": "ThyA (timidilato sintasa)",
    "Thymidylate_kin": "timidilato quinasa",
    "Sulfatase": "sulfatasa",
    "Sulfatase_C": "sulfatasa",
    "Sulfatase_C_2": "sulfatasa",
    "Arylsulfotrans": "arilsulfotransferasa",
    "Arylsulfotran_2": "arilsulfotransferasa",
    "Arylsulfotran_N": "arilsulfotransferasa",
    "Sulfotransfer_1": "sulfotransferasa",
    "Sulfotransfer_2": "sulfotransferasa",
    "Sulfotransfer_3": "sulfotransferasa",
    "Sulfotransfer_4": "sulfotransferasa",
    "Sulfotransfer_5": "sulfotransferasa",
    "Sulfotransfer_6": "sulfotransferasa",
}
FOCUS = re.compile(r"cgamp|hydroxymethyl-dcmp|sulfo|sulfogalactosyl|thymidylate", re.I)


def scan(faa: Path, panel: Path, out: Path) -> list[dict]:
    subprocess.run(
        ["hmmscan", "--cut_ga", "--domtblout", str(out), "--noali",
         str(panel), str(faa)],
        check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
    )
    hits = []
    for line in open(out):
        if line.startswith("#"):
            continue
        f = line.split()
        if len(f) >= 14:
            hits.append({"model": f[0], "query": f[3],
                         "function": MODEL_FUNCTION.get(f[0], f[0]),
                         "evalue": float(f[6]), "score": float(f[7])})
    return hits


def confirm(gaps: list[dict], proteins: list[dict], panel: Path) -> list[dict]:
    """Escanea proteínas si hay brechas de familias cubiertas por el panel."""
    if not any(FOCUS.search(g.get("gap", "")) for g in gaps):
        return []
    if not Path(panel).exists():
        return []
    with tempfile.TemporaryDirectory() as td:
        faa = Path(td) / "proteins.faa"
        with open(faa, "w") as fh:
            for p in proteins:
                fh.write(f">{p['id']}\n")
                aa = p["aa"]
                for i in range(0, len(aa), 60):
                    fh.write(aa[i:i + 60] + "\n")
        return scan(faa, Path(panel), Path(td) / "hits.domtbl")
