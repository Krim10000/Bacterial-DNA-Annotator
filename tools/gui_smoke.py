"""Prueba de humo de la GUI con Playwright (manual, no CI).

Uso:
    pip install playwright
    playwright install chromium
    python tools/gui_smoke.py http://127.0.0.1:8501

Comprueba:
- que la app carga y muestra el título,
- que el selector de idioma cambia a inglés y chino,
- que el navegador de demostración dibuja genes,
- que al hacer clic en un gen aparece la ficha.

Guarda capturas en `.scratch/gui_smoke/`. Si Playwright no está instalado,
termina con un mensaje claro en vez de fallar.
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

try:
    from playwright.sync_api import sync_playwright
except ImportError:
    print("Falta Playwright: pip install playwright && playwright install chromium")
    raise SystemExit(0)

URL = sys.argv[1] if len(sys.argv) > 1 else "http://127.0.0.1:8501"
OUT = Path(__file__).resolve().parents[1] / ".scratch" / "gui_smoke"
OUT.mkdir(parents=True, exist_ok=True)


def main() -> int:
    errors: list[str] = []
    with sync_playwright() as p:
        browser = p.chromium.launch()
        page = browser.new_page(viewport={"width": 1500, "height": 1100})
        page.on("pageerror", lambda exc: errors.append(str(exc)))

        page.goto(f"{URL}/?lang=es", wait_until="domcontentloaded")
        page.wait_for_selector("text=Anotador de ADN bacteriano", timeout=120_000)
        page.wait_for_selector("text=Genomas", timeout=120_000)
        time.sleep(3)
        page.screenshot(path=str(OUT / "01_es.png"))

        page.goto(f"{URL}/?lang=en", wait_until="domcontentloaded")
        page.wait_for_selector("text=Bacterial DNA annotator", timeout=120_000)
        time.sleep(2)
        page.screenshot(path=str(OUT / "02_en.png"))

        page.goto(f"{URL}/?lang=zh", wait_until="domcontentloaded")
        page.wait_for_selector("text=细菌 DNA 注释器", timeout=120_000)
        time.sleep(2)
        page.screenshot(path=str(OUT / "03_zh.png"))

        # Vista modular: opciones generadas desde modules/
        page.goto(f"{URL}/?modo=modulos&lang=es", wait_until="domcontentloaded")
        page.wait_for_selector("text=Estado de todos los módulos", timeout=120_000)
        page.wait_for_selector("text=CDS con pyrodigal", timeout=120_000)
        page.screenshot(path=str(OUT / "05_modules.png"))

        # Navegador con la demo incluida
        page.goto(f"{URL}/?modo=navegador&genoma=demo", wait_until="domcontentloaded")
        page.wait_for_selector("iframe", timeout=120_000)
        time.sleep(4)
        frames = [f for f in page.frames if f.locator("[data-kind=gene]").count()]
        if frames:
            frame = frames[0]
            n = frame.locator("[data-kind=gene]").count()
            frame.locator("[data-kind=gene] rect").first.click(force=True)
            time.sleep(5)
            body = page.inner_text("body")
            page.screenshot(path=str(OUT / "04_demo_click.png"))
            print(f"demo: {n} genes; ficha={'Has hecho clic' in body}")
        else:
            errors.append("el navegador demo no dibujó genes")

        browser.close()

    if errors:
        print("ERRORES:")
        for e in errors[:10]:
            print(" -", e)
        return 1
    print(f"OK. Capturas en {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
