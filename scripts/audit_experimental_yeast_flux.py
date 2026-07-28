from __future__ import annotations

import argparse
import json

from biocandidate.experimental_flux_audit import audit


def main() -> None:
    parser = argparse.ArgumentParser(description="Audit targeted experimental yeast flux data")
    parser.add_argument("--model", required=True, help="Path to the Yeast-MetaTwin YAML asset")
    parser.add_argument(
        "--output",
        default="artifacts/experimental-yeast-flux/blocker.json",
        help="Machine-readable blocker output path",
    )
    args = parser.parse_args()
    print(json.dumps(audit(args.model, args.output), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
