"""Descubrimiento automático de módulos.

Cada carpeta bajo ``modules/`` con un ``module.py`` que exporte ``MODULE`` es
un módulo. No hay lista central que mantener: añadir la carpeta y el módulo
aparece en la GUI y en ``python -m modules --list``.
"""
from __future__ import annotations

import importlib
import pkgutil
from functools import lru_cache
from pathlib import Path

from .base import Module, ModuleStatus

PACKAGE = __name__.rsplit(".", 1)[0]


@lru_cache(maxsize=1)
def discover() -> tuple[Module, ...]:
    """Importa y ordena todos los módulos disponibles."""
    found: dict[str, Module] = {}
    package = importlib.import_module(PACKAGE)
    for info in pkgutil.iter_modules(package.__path__):
        if info.name.startswith("_") or info.name in {"base", "registry"}:
            continue
        try:
            mod = importlib.import_module(f"{PACKAGE}.{info.name}.module")
        except Exception as exc:  # noqa: BLE001
            # Un módulo roto no debe tumbar la app: se reporta como no disponible.
            found[f"broken:{info.name}"] = Module(
                id=f"broken:{info.name}",
                name=f"⚠️ {info.name}",
                summary=f"No se pudo importar: {exc}",
                order=999,
                status_fn=lambda _root, msg=str(exc): ModuleStatus(False, (msg,)),
            )
            continue
        spec = getattr(mod, "MODULE", None)
        if not isinstance(spec, Module):
            continue
        if spec.id in found:
            raise RuntimeError(f"id de módulo duplicado: {spec.id}")
        found[spec.id] = spec
    return tuple(sorted(found.values(), key=lambda m: (m.order, m.name)))


def get(module_id: str) -> Module | None:
    return next((m for m in discover() if m.id == module_id), None)


def by_id() -> dict[str, Module]:
    return {m.id: m for m in discover()}


def list_ids() -> list[str]:
    return [m.id for m in discover()]


def status_all(data_root: Path) -> dict[str, ModuleStatus]:
    return {m.id: m.status(data_root) for m in discover()}
