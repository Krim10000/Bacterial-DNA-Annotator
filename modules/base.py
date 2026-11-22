"""Contrato de los módulos del anotador.

Un módulo declara **qué sabe hacer, qué opciones ofrece, qué necesita y cómo
se prueba**. La GUI descubre los módulos y construye sus controles sin que haya
que editar `app.py`.

Reglas de diseño
----------------
- El módulo no importa Streamlit: la interfaz se genera desde sus datos.
- Las opciones tienen tipo, valor por defecto, rango y dependencias.
- ``status()`` dice si se puede usar y por qué no.
- ``selftest()`` es rápido y no descarga datos: sirve para CI y para el botón
  «Probar módulo».
- ``run()`` recibe un contexto y devuelve un resultado con rutas y métricas.
"""
from __future__ import annotations

import shutil
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Sequence


@dataclass(frozen=True)
class Option:
    """Una opción de configuración que la GUI dibuja sola."""

    key: str
    label: str                     # texto en español; la GUI lo traduce
    kind: str                      # bool | int | float | text | select | multiselect
    default: Any
    help: str = ""
    choices: tuple[Any, ...] = ()
    min: float | None = None
    max: float | None = None
    step: float | None = None
    advanced: bool = False
    depends_on: str | None = None  # otra opción que debe ser verdadera
    unit: str = ""

    def __post_init__(self) -> None:
        if self.kind not in {"bool", "int", "float", "text", "select",
                             "multiselect"}:
            raise ValueError(f"tipo de opción no soportado: {self.kind}")
        if self.kind in {"select", "multiselect"} and not self.choices:
            raise ValueError(f"{self.key}: falta choices")
        if self.kind == "select" and self.default not in self.choices:
            raise ValueError(f"{self.key}: default fuera de choices")

    def clean(self, value: Any) -> Any:
        """Convierte y valida un valor; lanza ValueError si no sirve."""
        if self.kind == "bool":
            return bool(value)
        if self.kind in {"int", "float"}:
            num = int(value) if self.kind == "int" else float(value)
            if self.min is not None and num < self.min:
                raise ValueError(f"{self.key}: mínimo {self.min}")
            if self.max is not None and num > self.max:
                raise ValueError(f"{self.key}: máximo {self.max}")
            return num
        if self.kind == "select":
            if value not in self.choices:
                raise ValueError(f"{self.key}: opción inválida")
            return value
        if self.kind == "multiselect":
            vals = list(value)
            bad = [v for v in vals if v not in self.choices]
            if bad:
                raise ValueError(f"{self.key}: opciones inválidas {bad}")
            return vals
        return str(value)


@dataclass
class ModuleStatus:
    """Disponibilidad del módulo en esta máquina."""

    available: bool
    reasons: tuple[str, ...] = ()
    detail: dict[str, Any] = field(default_factory=dict)

    def badge(self) -> str:
        return "✅" if self.available else "⚠️"


@dataclass
class RunContext:
    """Todo lo que un módulo necesita para ejecutarse."""

    data_root: Path
    out_root: Path
    config: dict[str, Any]
    log: Callable[[str], None] = print
    results: dict[str, Any] = field(default_factory=dict)


@dataclass
class ModuleResult:
    """Salida uniforme de cualquier módulo."""

    module: str
    ok: bool
    elapsed: float = 0.0
    outputs: dict[str, Path] = field(default_factory=dict)
    metrics: dict[str, Any] = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)
    error: str | None = None


RunFn = Callable[[RunContext, Any], ModuleResult]
SelfTestFn = Callable[[], str]
EstimateFn = Callable[[dict[str, Any]], str]
StatusFn = Callable[[Path], ModuleStatus]


@dataclass
class Module:
    """Especificación completa de un módulo."""

    id: str
    name: str
    summary: str
    order: int = 100
    category: str = "general"
    version: str = "1.0"
    requires: tuple[str, ...] = ()       # binarios o módulos Python
    gpu: bool = False
    experimental: bool = False
    options: Sequence[Option] = ()
    docs: str = ""
    #: Qué entrada consume la GUI: fasta | proteins | genome | None.
    input_kind: str | None = None
    #: Texto del botón de ejecución (los módulos de solo lectura «resumen»).
    run_label: str = "Ejecutar módulo"
    run_fn: RunFn | None = None
    selftest_fn: SelfTestFn | None = None
    estimate_fn: EstimateFn | None = None
    status_fn: StatusFn | None = None

    def defaults(self) -> dict[str, Any]:
        return {o.key: o.default for o in self.options}

    def config(self, values: dict[str, Any] | None = None) -> dict[str, Any]:
        """Mezcla valores con defaults y valida cada opción."""
        raw = dict(self.defaults())
        if values:
            unknown = set(values) - set(raw)
            if unknown:
                raise ValueError(f"{self.id}: opciones desconocidas {sorted(unknown)}")
            raw.update(values)
        out = {}
        for opt in self.options:
            out[opt.key] = opt.clean(raw[opt.key])
        return out

    def status(self, data_root: Path) -> ModuleStatus:
        if self.status_fn is not None:
            return self.status_fn(data_root)
        missing = [r for r in self.requires if shutil.which(r) is None]
        if missing:
            return ModuleStatus(False, (f"Falta: {', '.join(missing)}",))
        return ModuleStatus(True)

    def estimate(self, config: dict[str, Any] | None = None) -> str:
        if self.estimate_fn is not None:
            return self.estimate_fn(config or self.defaults())
        return "depende del tamaño del genoma"

    def selftest(self) -> str:
        if self.selftest_fn is None:
            return "módulo registrado (sin autoprueba específica)"
        return self.selftest_fn()

    def run(self, ctx: RunContext, inputs: Any = None) -> ModuleResult:
        if self.run_fn is None:
            return ModuleResult(
                module=self.id, ok=False,
                error="Este módulo todavía no tiene ejecución desde la GUI; "
                      "consulta su pestaña de exploración.",
            )
        start = time.time()
        try:
            res = self.run_fn(ctx, inputs)
            res.elapsed = round(time.time() - start, 2)
            return res
        except Exception as exc:  # noqa: BLE001 - se reporta al usuario
            return ModuleResult(module=self.id, ok=False,
                                elapsed=round(time.time() - start, 2),
                                error=str(exc))
