"""Generate Figure 1: the data and task map.

Plots every corpus the project reads on two axes that decide whether a result
transfers: how novel the sequences are relative to training (homology-cold to
near-identical) and which endpoint the labels carry. This makes the scope of
each claim visible in one panel.

Writes a dependency-free SVG plus a markdown legend, so the figure can be
reviewed in git without a plotting toolchain.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

DEFAULT_OUTPUT = Path("manuscript/figures/figure1-data-task-map.svg")

# x: 0 = homology-cold, 3 = near-identical to training sequences
# y: 0 = absolute kinetic constant, 1 = relative fitness ranking, 2 = flux
TASKS = {
    "absolute": "Absolute kinetic constant",
    "ranking": "Relative engineering fitness",
    "flux": "Metabolic flux",
}

CORPORA = (
    {
        "id": "temporal-pool",
        "label": "Temporal absolute-kinetics pool",
        "novelty": 0.2,
        "task": "absolute",
        "records": 192,
        "use": "never scored",
        "license": "CC0 / CC-BY-4.0 per source",
        "note": "Below its 300-record and 30-family gates; predictions are not permitted.",
    },
    {
        "id": "enzengdb",
        "label": "EnzEngDB v1",
        "novelty": 0.1,
        "task": "ranking",
        "records": 6423,
        "use": "zero-shot external",
        "license": "CC-BY-4.0",
        "note": "51 homology-cold campaigns; macro-averaged Spearman 0.071.",
    },
    {
        "id": "lunzer",
        "label": "Lunzer IMDH landscape",
        "novelty": 2.6,
        "task": "absolute",
        "records": 512,
        "use": "zero-shot external",
        "license": "CC0-1.0",
        "note": "67-69% identity to training homologs: a mutation-sensitivity test, not homology-cold.",
    },
    {
        "id": "sabio-rk",
        "label": "SABIO-RK Km rows",
        "novelty": 0.2,
        "task": "absolute",
        "records": 0,
        "use": "blocked",
        "license": "not established",
        "note": "Permission request drafted but never sent.",
    },
    {
        "id": "internal-test",
        "label": "Frozen homology-cold test",
        "novelty": 0.0,
        "task": "absolute",
        "records": 1646,
        "use": "internal test, already observed",
        "license": "derived from training corpus",
        "note": "All headline RMSE numbers come from here. Not independent: it has been scored.",
    },
    {
        "id": "catapro-internal",
        "label": "CataPro retrained",
        "novelty": 0.0,
        "task": "absolute",
        "records": 1646,
        "use": "same-split comparator",
        "license": "code MIT; data unresolved",
        "note": "1,127 of 1,646 test pairs overlap CataPro's published ten-fold data.",
    },
    {
        "id": "unikp-mode-b",
        "label": "UniKP retrained (Mode B)",
        "novelty": 0.0,
        "task": "absolute",
        "records": 1646,
        "use": "same-split comparator",
        "license": "code MIT; data unresolved",
        "note": "Retrained on the frozen training partition, so it shares the test.",
    },
    {
        "id": "dlkcat-mode-b",
        "label": "DLKcat retrained (Mode B)",
        "novelty": 0.0,
        "task": "absolute",
        "records": 1646,
        "use": "same-split comparator",
        "license": "GPL-3.0-only",
        "note": "Retrained on the frozen training partition.",
    },
    {
        "id": "development",
        "label": "UniKP/DLKcat kcat corpus",
        "novelty": 3.0,
        "task": "absolute",
        "records": 16838,
        "use": "training",
        "license": "unknown",
        "note": "The training corpus itself; by definition not independent.",
    },
)

WIDTH = 1220
HEIGHT = 700
PLOT_LEFT = 392
PLOT_RIGHT = 840
LANE_TOP = 150
LANE_STEP = 48
LABEL_RIGHT = 322
DETAIL_LEFT = 872

TASK_COLORS = {
    "absolute": "#2c6fbb",
    "ranking": "#c8553d",
    "flux": "#6a994e",
}


def _x(value: float) -> float:
    return PLOT_LEFT + (value / 3.0) * (PLOT_RIGHT - PLOT_LEFT)


def _lane_y(index: int) -> float:
    return LANE_TOP + index * LANE_STEP


def _radius(records: int) -> float:
    if records <= 0:
        return 4.0
    # log scaling keeps a 16K corpus and a 192-record pool both legible
    scale = max(0.0, min(1.0, (len(str(records)) - 2.2) / 2.2))
    return 6.0 + scale * 10.0


def _escape(text: str) -> str:
    return text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def build_svg(corpora: tuple[dict[str, Any], ...]) -> str:
    """One lane per corpus: collisions are impossible, so every label stays legible."""
    parts: list[str] = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{WIDTH}" height="{HEIGHT}" '
        f'viewBox="0 0 {WIDTH} {HEIGHT}" font-family="Helvetica, Arial, sans-serif">',
        f'<rect width="{WIDTH}" height="{HEIGHT}" fill="white"/>',
        '<text x="40" y="42" font-size="20" font-weight="bold">'
        "Figure 1. Data and task map</text>",
        '<text x="40" y="68" font-size="13" fill="#555">'
        "Where each corpus sits, and what it can support. No corpus is both homology-cold "
        "and never scored.</text>",
    ]

    bottom = _lane_y(len(corpora) - 1) + 40

    # shaded target strip: independent homology-cold evidence would sit here
    parts.append(
        f'<rect x="{_x(0) - 26}" y="{LANE_TOP - 34}" width="52" height="{bottom - LANE_TOP - 12}" '
        'fill="#2c6fbb" fill-opacity="0.06"/>'
    )
    parts.append(
        f'<text x="{_x(0)}" y="{LANE_TOP - 44}" font-size="11" text-anchor="middle" '
        'fill="#2c6fbb">independent evidence goes here</text>'
    )

    for index, corpus in enumerate(corpora):
        y = _lane_y(index)
        color = TASK_COLORS[corpus["task"]]
        dashed = corpus["use"] == "never scored"

        parts.append(
            f'<line x1="40" y1="{y}" x2="{PLOT_RIGHT + 30}" y2="{y}" stroke="#f0f0f0"/>'
        )
        parts.append(
            f'<text x="{LABEL_RIGHT}" y="{y + 4}" font-size="12.5" text-anchor="end" '
            f'fill="#111">{_escape(corpus["label"])}</text>'
        )

        x = _x(corpus["novelty"])
        radius = _radius(corpus["records"])
        parts.append(
            f'<circle cx="{x:.1f}" cy="{y}" r="{radius:.1f}" fill="{color}" fill-opacity="0.35" '
            f'stroke="{color}" stroke-width="2"'
            + (' stroke-dasharray="4 3"' if dashed else "")
            + "/>"
        )
        if corpus["records"] == 0:
            parts.append(
                f'<text x="{x:.1f}" y="{y + 4}" font-size="11" text-anchor="middle" '
                'fill="#555">0</text>'
            )

        parts.append(
            f'<text x="{DETAIL_LEFT}" y="{y + 4}" font-size="12" fill="#333">'
            f'{corpus["records"]:,} records</text>'
        )
        parts.append(
            f'<text x="{DETAIL_LEFT + 130}" y="{y + 4}" font-size="12" fill="#666">'
            f'{_escape(corpus["use"])}</text>'
        )

    # x axis
    parts.append(
        f'<line x1="{PLOT_LEFT}" y1="{bottom}" x2="{PLOT_RIGHT}" y2="{bottom}" stroke="#999"/>'
    )
    tick_labels = {
        0: "homology-cold",
        1: "partly novel",
        2: "close homolog",
        3: "near-identical",
    }
    for value, label in tick_labels.items():
        x = _x(value)
        parts.append(
            f'<line x1="{x}" y1="{bottom}" x2="{x}" y2="{bottom + 6}" stroke="#999"/>'
        )
        parts.append(
            f'<text x="{x}" y="{bottom + 22}" font-size="11.5" text-anchor="middle" '
            f'fill="#333">{_escape(label)}</text>'
        )
    parts.append(
        f'<text x="{(PLOT_LEFT + PLOT_RIGHT) / 2}" y="{bottom + 48}" font-size="13" '
        'text-anchor="middle" fill="#111">Sequence novelty relative to training data</text>'
    )

    # legend
    legend_y = bottom + 82
    legend_x = 40
    for task, label in TASKS.items():
        parts.append(
            f'<circle cx="{legend_x}" cy="{legend_y}" r="6" fill="{TASK_COLORS[task]}" '
            f'fill-opacity="0.35" stroke="{TASK_COLORS[task]}" stroke-width="2"/>'
        )
        parts.append(
            f'<text x="{legend_x + 14}" y="{legend_y + 4}" font-size="12" fill="#333">'
            f"{_escape(label)}</text>"
        )
        legend_x += 40 + len(label) * 6.6
    parts.append(
        f'<circle cx="{legend_x}" cy="{legend_y}" r="6" fill="none" stroke="#666" '
        'stroke-width="2" stroke-dasharray="4 3"/>'
    )
    parts.append(
        f'<text x="{legend_x + 14}" y="{legend_y + 4}" font-size="12" fill="#333">'
        "never scored</text>"
    )
    parts.append(
        f'<text x="40" y="{legend_y + 26}" font-size="11.5" fill="#666">'
        "Marker area scales with record count. Lanes run from the most independent "
        "evidence at the top to the training corpus at the bottom.</text>"
    )

    parts.append("</svg>")
    return "\n".join(parts)


def build_legend(corpora: tuple[dict[str, Any], ...]) -> str:
    lines = ["# Figure 1 legend", ""]
    lines.append("| Corpus | Task | Records | Use | License | Note |")
    lines.append("|---|---|---:|---|---|---|")
    for corpus in corpora:
        lines.append(
            f"| {corpus['label']} | {TASKS[corpus['task']]} | {corpus['records']:,} | "
            f"{corpus['use']} | {corpus['license']} | {corpus['note']} |"
        )
    lines.append("")
    lines.append(
        "Read the panel as: no corpus sits in the top-left where an independent, "
        "homology-cold, never-scored benchmark would belong. The internal test has been "
        "scored, the comparators share it, EnzEngDB carries a different endpoint, and the "
        "temporal pool has never been scored but is below its own gate."
    )
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--legend", type=Path, default=None)
    parser.add_argument("--json-out", type=Path, default=None)
    args = parser.parse_args()

    root = args.root.resolve()
    output = args.output if args.output.is_absolute() else root / args.output
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(build_svg(CORPORA), encoding="utf-8")

    legend_path = args.legend or output.with_suffix(".legend.md")
    legend_path = legend_path if legend_path.is_absolute() else root / legend_path
    legend_path.write_text(build_legend(CORPORA), encoding="utf-8")

    summary = {
        "figure": output.relative_to(root).as_posix(),
        "legend": legend_path.relative_to(root).as_posix(),
        "corpus_count": len(CORPORA),
        "scored_corpora": [c["id"] for c in CORPORA if c["use"] != "never scored"],
        "never_scored": [c["id"] for c in CORPORA if c["use"] == "never scored"],
        "claim_boundary": (
            "The figure is generated from recorded project state, not measured here. "
            "Novelty positions are qualitative placements documented in the legend."
        ),
    }
    if args.json_out:
        json_path = args.json_out if args.json_out.is_absolute() else root / args.json_out
        json_path.parent.mkdir(parents=True, exist_ok=True)
        json_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
