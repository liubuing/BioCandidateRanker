from __future__ import annotations

import argparse
import json

from biocandidate.simulated_flux import run_pipeline


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Generate and train a governed simulated FBA flux smoke task")
    parser.add_argument("--model", required=True, help="Read-only Yeast-MetaTwin YAML model")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--seed", type=int, default=7)
    args = parser.parse_args()
    try:
        result = run_pipeline(args.model, args.output_dir, seed=args.seed)
    except RuntimeError as exc:
        print(json.dumps({"status": "blocked", "reason": str(exc)}, sort_keys=True))
        return 1
    print(json.dumps(result["metrics"], sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
