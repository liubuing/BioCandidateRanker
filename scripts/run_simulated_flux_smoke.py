"""Smoke-test the governed simulated FBA flux task end to end.

Calls biocandidate.simulated_flux.run_pipeline on a read-only Yeast-MetaTwin
YAML model: generates the simulated flux dataset in --output-dir, runs the CPU
training smoke with the given seed, and prints the resulting metrics as JSON.
Fail-closed: any pipeline failure clears partial outputs, writes blocker.json
carrying the claim boundary, and this wrapper prints a {"status": "blocked"}
record and exits non-zero.

中文：端到端冒烟受治理的模拟 FBA 通量任务：对只读 Yeast-MetaTwin YAML 模型
调用 run_pipeline，生成模拟通量数据集并跑 CPU 训练冒烟；受阻即非零退出。
"""

from __future__ import annotations

import argparse
import json

from biocandidate.simulated_flux import run_pipeline


def main() -> int:
    """Run the smoke pipeline, mapping a blocked run to a non-zero exit code.
中文：跑冒烟流水线；被阻塞的运行映射为非零退出码。
"""
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
