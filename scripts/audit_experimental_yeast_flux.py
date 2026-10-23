"""CLI wrapper that audits targeted experimental yeast flux data, fail-closed.

Reads the Yeast-MetaTwin YAML model asset given by --model and delegates every
check to biocandidate.experimental_flux_audit.audit, which writes the
machine-readable blocker (default artifacts/experimental-yeast-flux/blocker.json).
The printed JSON is the audit result; any failed gate is recorded there rather
than raised.
"""

from __future__ import annotations

import argparse
import json

from biocandidate.experimental_flux_audit import audit


def main() -> None:
    """Parse CLI arguments, run the audit, and print its blocker JSON."""
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
