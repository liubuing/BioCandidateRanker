"""Attach the yeast FBA context vector to BioCandidateRanker predictions.

Demonstrates the second fusion seam: the Yeast-MetaTwin model now registered
in biocandidate.fba_context feeds its governed 8-wide vector through the
workbench prediction path, alongside the kcat values used by the kcat-FBA
pilot.

Claim boundary (untrained input pathway): every shipped checkpoint was trained
without flux labels, so the FBA input pathway is untrained. Attaching the
context changes the output numbers only to prove the pipeline is connected; it
is NOT a scientifically validated input. Raw (context-free) predictions are
stored alongside for comparison.

Reads the governed FBA context from ``biocandidate.fba_context`` and the task-1
candidate list under ``artifacts/10h2da-yeast-prediction``; writes
``artifacts/fusion-yeast-pilot/context_attached_prediction_demo.json`` with raw
versus context-attached kcat values and their shifts.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

import serve_workbench as wb  # noqa: E402
from biocandidate import fba_context  # noqa: E402

MODEL_PATH = "artifacts/full-homology-baseline/best.pt"
TASK1 = ROOT / "artifacts/10h2da-yeast-prediction/task1_10h2da_pathway.json"


def main() -> int:
    """Run raw and FBA-attached predictions and record the demonstration artifact."""
    import torch

    device = "cuda" if torch.cuda.is_available() else "cpu"
    cache = wb.ModelCache(ROOT, device=device, capacity=1)

    context = fba_context.run_fba_context(model_id="yeast_metatwin")
    payload = {
        "metadata": context["metadata"],
        "features": context["features"],
    }
    candidates = [
        {k: row[k] for k in ("candidate_id", "sequence", "substrate_smiles",
                             "substrate_name", "organism", "ec",
                             "enzyme_type", "reaction")}
        for row in json.loads(TASK1.read_text(encoding="utf-8"))["candidates"]
    ]

    raw = wb.run_prediction(cache, MODEL_PATH, candidates, None)
    attached = wb.run_prediction(cache, MODEL_PATH, candidates, None, fba_payload=payload)

    def kcat(result):
        return {p["candidate_id"]: p["tasks"]["log10_kcat"]["mean"] for p in result["predictions"]}

    raw_kcat = kcat(raw)
    attached_kcat = kcat(attached)

    out = {
        "demo": "fusion_yeast_context_attach_demo",
        "claim_boundary": (
            "Untrained input pathway demonstration only: the checkpoint was trained "
            "without flux labels, so attached-context outputs have no validated meaning."
        ),
        "fba_context": {
            "model_id": context["model_id"],
            "condition_id": context["condition_id"],
            "medium": context["medium"],
            "growth": context["growth"],
            "feature_ids": context["feature_ids"],
            "features": context["features"],
            "provenance": context["provenance"],
        },
        "raw_log10_kcat": raw_kcat,
        "context_attached_log10_kcat": attached_kcat,
        "shift_log10": {
            k: attached_kcat[k] - raw_kcat[k] for k in raw_kcat
        },
        "input_identity": raw["input_identity"],
        "checkpoint": raw["checkpoint"],
    }
    out_path = ROOT / "artifacts/fusion-yeast-pilot/context_attached_prediction_demo.json"
    out_path.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    print(f"FBA context: growth={context['growth']:.4f} "
          f"condition_id={context['condition_id'][:19]}...")
    for k in raw_kcat:
        print(f"  {k:42s} raw={raw_kcat[k]:+.4f} attached={attached_kcat[k]:+.4f} "
              f"shift={attached_kcat[k]-raw_kcat[k]:+.4f}")
    print(f"-> {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
