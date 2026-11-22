"""Parche de compatibilidad para ejecutar ESMFold con versiones recientes de fair-esm."""
from pathlib import Path
p=Path("/usr/local/lib/python3.11/site-packages/esm/esmfold/v1/trunk.py")
s=p.read_text()
s=s.replace("from dataclasses import dataclass", "from dataclasses import dataclass, field")
s=s.replace("structure_module: StructureModuleConfig = StructureModuleConfig()",
            "structure_module: StructureModuleConfig = field(default_factory=StructureModuleConfig)")
p.write_text(s)
print("patched fair-esm")
