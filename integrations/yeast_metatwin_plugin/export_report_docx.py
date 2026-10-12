"""Export the plugin integration report to DOCX.

Reuses scripts/export_manuscript_docx.py's markdown converter so the Word output
keeps the same heading, table, list and inline-emphasis handling as the
manuscript exports; nothing about that converter is duplicated here.

Usage:
    python integrations/yeast_metatwin_plugin/export_report_docx.py
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CONVERTER = ROOT / "scripts/export_manuscript_docx.py"

SOURCE = HERE / "PLUGIN_INTEGRATION_REPORT_CN.md"
TARGET = HERE / "BioCandidateRanker插件接入报告_20261012.docx"


def _load_converter():
    """Import the manuscript exporter without copying its logic."""
    spec = importlib.util.spec_from_file_location("manuscript_docx_export", CONVERTER)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main() -> int:
    if not SOURCE.is_file():
        raise SystemExit(f"missing report source: {SOURCE}")
    if not CONVERTER.is_file():
        raise SystemExit(f"missing converter: {CONVERTER}")

    converter = _load_converter()
    converter.convert(SOURCE, TARGET)
    print(f"bytes: {TARGET.stat().st_size:,}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
