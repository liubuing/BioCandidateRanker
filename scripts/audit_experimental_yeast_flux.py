"""CLI wrapper that audits targeted experimental yeast flux data, fail-closed.

Reads the Yeast-MetaTwin YAML model asset given by --model and delegates every
check to biocandidate.experimental_flux_audit.audit, which writes the
machine-readable blocker (default artifacts/experimental-yeast-flux/blocker.json).
The printed JSON is the audit result; any failed gate is recorded there rather
than raised.

中文：本脚本是审计定向实验酵母通量数据的 CLI 包装，整体 fail-closed：读取 --model 指定
的 Yeast-MetaTwin YAML 模型资产，把全部检查委托给
biocandidate.experimental_flux_audit.audit，由后者写出机器可读阻断记录（默认
artifacts/experimental-yeast-flux/blocker.json）。打印的 JSON 即审计结果；任何失败的门
都会被记录在案而不是抛出异常。治理角色：fail-closed 阻断记录的命令行入口。
"""

from __future__ import annotations

import argparse
import json

from biocandidate.experimental_flux_audit import audit


def main() -> None:
    """Parse CLI arguments, run the audit, and print its blocker JSON.

    中文：解析 CLI 参数，运行审计并打印其阻断 JSON。
    """
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
