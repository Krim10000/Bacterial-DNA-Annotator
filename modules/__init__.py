"""Módulos autodescubiertos del anotador.

Uso:
    from modules import discover, get
    from modules.base import RunContext

Cada módulo vive en ``modules/<id>/module.py`` y exporta ``MODULE``.
"""
from .base import Module, ModuleResult, ModuleStatus, Option, RunContext
from .registry import discover, get, list_ids

__all__ = ["Module", "ModuleResult", "ModuleStatus", "Option", "RunContext",
           "discover", "get", "list_ids"]
