"""Export any markdown report in this repository to DOCX.

Delegates to the converter in scripts/export_manuscript_docx.py so report
exports keep the same heading, table, list, blockquote and inline-emphasis
handling as the manuscript exports; that converter is not duplicated here.

Reads a UTF-8 Markdown source and writes a .docx next to it (or at --target).
Governance role: presentation-only export wrapper; it never mutates the
source report or any artifact consumed by evaluations.

Usage:
    python scripts/export_docx_report.py --source docs/PROJECT_STATUS_REPORT_20261012_CN.md
    python scripts/export_docx_report.py --source a.md --target b.docx

中文：把仓库内任意 markdown 报告导出为 DOCX，复用手稿导出器的标题/表格/
列表/引用/行内强调处理，保证与手稿导出一致。
"""

from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CONVERTER = ROOT / "scripts/export_manuscript_docx.py"


def load_converter():
    """Load scripts/export_manuscript_docx.py as a module by file path.
中文：按文件路径把 export_manuscript_docx.py 作为模块加载。
"""
    spec = importlib.util.spec_from_file_location("manuscript_docx_export", CONVERTER)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True, help="markdown report to export")
    parser.add_argument("--target", type=Path, default=None,
                        help="output .docx (default: same path with .docx)")
    args = parser.parse_args()

    source = args.source if args.source.is_absolute() else ROOT / args.source
    if not source.is_file():
        raise SystemExit(f"missing source: {source}")
    if not CONVERTER.is_file():
        raise SystemExit(f"missing converter: {CONVERTER}")

    target = args.target or source.with_suffix(".docx")
    if not target.is_absolute():
        target = ROOT / target

    load_converter().convert(source, target)
    print(f"bytes: {target.stat().st_size:,}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
