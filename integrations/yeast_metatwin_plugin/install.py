"""Install (or remove) the BioCandidateRanker plugin in the Yeast-MetaTwin Platform.

Writes only these paths inside the platform:

    plugin_runtime/biocandidate.py                 the adapter
    BioCandidateRanker/benchmark.py                manifest builder
    BioCandidateRanker/benchmark_manifest.json     generated
    BioCandidateRanker/readiness_manifest.json     generated
    BioCandidateRanker/ranker_root.txt             path to this checkout
    BioCandidateRanker/assets/esm2-t6-opt-seed7.best.pt
    BioCandidateRanker/LICENSE                     MIT, copied from this repo
    catalog.py                                     one delimited appended block

The catalog edit is a single marker-delimited block appended at end of file, so
--uninstall restores the file byte-for-byte and re-running install is a no-op.
Everything else is new files that --uninstall deletes.

Usage
    python install.py --platform "D:/biological/Metabolic model prediction/Integrated_Yeast_MetaTwin_Deployment"
    python install.py --platform ... --uninstall
    python install.py --platform ... --dry-run
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
RANKER_ROOT = HERE.parents[1]

BEGIN_MARKER = "# --- BEGIN BioCandidateRanker plugin (integrations/yeast_metatwin_plugin/install.py) ---"
END_MARKER = "# --- END BioCandidateRanker plugin ---"

CHECKPOINT_SOURCE = RANKER_ROOT / "artifacts/esm2-t6-opt-seed7/best.pt"
CHECKPOINT_TARGET_NAME = "esm2-t6-opt-seed7.best.pt"

PLUGIN_DIR_NAME = "BioCandidateRanker"
RANKER_ROOT_FILE = "ranker_root.txt"


def catalog_block(plugin_version: str, assets: tuple[tuple[str, int], ...],
                  packages: tuple[tuple[str, tuple[int, ...] | None, tuple[int, ...] | None], ...],
                  smoke: dict) -> str:
    """Render the appended catalog block. Pure data; no imports at render time."""
    asset_lines = "\n".join(
        f'    AssetSpec("{path}", min_size_bytes={size}),' for path, size in assets
    )
    package_lines = []
    for distribution, minimum, maximum in packages:
        args = [f'"{distribution}"']
        if minimum:
            args.append(str(minimum))
        if maximum:
            args.append(str(maximum))
        package_lines.append("        PackageSpec(" + ", ".join(args) + "),")
    packages_rendered = "\n".join(package_lines)

    return f"""
{BEGIN_MARKER}
# Appended by the BioCandidateRanker plugin installer. This block is the entire
# change to this file; removing it restores the original content exactly.
_BIOCANDIDATE_ASSETS = (
{asset_lines}
)

_BIOCANDIDATE_SPEC = PluginSpec(
    name="BioCandidateRanker",
    plugin_version="{plugin_version}",
    capabilities=("kcat_prediction",),
    required_assets=_BIOCANDIDATE_ASSETS,
    packages=(
{packages_rendered}
    ),
    source_url="https://github.com/liubuing/BioCandidateRanker",
)

PLUGIN_SPECS = PLUGIN_SPECS + (_BIOCANDIDATE_SPEC,)

FIXED_SMOKE_INPUTS = dict(FIXED_SMOKE_INPUTS)
FIXED_SMOKE_INPUTS["BioCandidateRanker"] = {json.dumps(smoke, indent=4).replace(chr(10), chr(10) + "    ")}
{END_MARKER}
"""


def locate_platform(root: Path) -> tuple[Path, Path, Path]:
    plugin_runtime = root / "04_prediction_plugins" / "plugin_runtime"
    plugin_dir = root / "04_prediction_plugins" / PLUGIN_DIR_NAME
    catalog = plugin_runtime / "catalog.py"
    if not catalog.is_file():
        raise SystemExit(f"not a Yeast-MetaTwin Platform deployment: {catalog} not found")
    return plugin_runtime, plugin_dir, catalog


def read_text_exact(path: Path) -> str:
    """Decode without newline translation, so the original bytes survive."""
    return path.read_bytes().decode("utf-8")


def write_text_exact(path: Path, text: str) -> None:
    """Write without newline translation.

    Path.write_text translates 
 to os.linesep, which silently converts a
    whole LF file to CRLF on Windows - that is how the first version of this
    installer broke byte-exact uninstall.
    """
    path.write_bytes(text.encode("utf-8"))


def dominant_newline(text: str) -> str:
    crlf = text.count("\r\n")
    lf_only = text.count("\n") - crlf
    return "\r\n" if crlf > lf_only else "\n"


def strip_block(text: str) -> tuple[str, bool]:
    """Remove the appended block and the blank separator line it added."""
    if BEGIN_MARKER not in text:
        return text, False
    start = text.index(BEGIN_MARKER)
    end = text.index(END_MARKER) + len(END_MARKER)
    remainder = text[end:]
    if remainder.startswith("\r\n"):
        remainder = remainder[2:]
    elif remainder.startswith("\n"):
        remainder = remainder[1:]
    newline = dominant_newline(text)
    head = text[:start].rstrip("\r\n") + newline
    return head + remainder.lstrip("\r\n"), True


def install(root: Path, dry_run: bool) -> int:
    plugin_runtime, plugin_dir, catalog = locate_platform(root)
    assets_dir = plugin_dir / "assets"

    if not CHECKPOINT_SOURCE.is_file():
        raise SystemExit(f"missing checkpoint to install: {CHECKPOINT_SOURCE}")

    actions = [
        f"write {plugin_runtime / 'biocandidate.py'}",
        f"write {plugin_dir / 'benchmark.py'}",
        f"write {assets_dir / CHECKPOINT_TARGET_NAME}",
        f"write {plugin_dir / 'LICENSE'}",
        f"write {plugin_dir / RANKER_ROOT_FILE}",
        f"append one delimited block to {catalog}",
        f"generate {plugin_dir / 'benchmark_manifest.json'} and readiness_manifest.json",
    ]
    if dry_run:
        print(json.dumps({"dry_run": True, "actions": actions}, indent=2))
        return 0

    plugin_runtime.mkdir(parents=True, exist_ok=True)
    assets_dir.mkdir(parents=True, exist_ok=True)

    shutil.copy2(HERE / "biocandidate.py", plugin_runtime / "biocandidate.py")
    shutil.copy2(HERE / "benchmark.py", plugin_dir / "benchmark.py")
    if (RANKER_ROOT / "LICENSE").is_file():
        shutil.copy2(RANKER_ROOT / "LICENSE", plugin_dir / "LICENSE")
    write_text_exact(
        plugin_dir / RANKER_ROOT_FILE,
        str(RANKER_ROOT).replace("\\", "/") + "\n",
    )

    target_checkpoint = assets_dir / CHECKPOINT_TARGET_NAME
    if not target_checkpoint.is_file() or (
        target_checkpoint.stat().st_size != CHECKPOINT_SOURCE.stat().st_size
    ):
        shutil.copy2(CHECKPOINT_SOURCE, target_checkpoint)

    # catalog: keep the pre-install content so --uninstall is byte-exact
    original = read_text_exact(catalog)
    backup = plugin_dir / "catalog.py.orig"
    stripped, had_block = strip_block(original)
    if not had_block and not backup.is_file():
        write_text_exact(backup, original)
    base = stripped
    newline = dominant_newline(original)

    from catalog_entry import ASSETS, PACKAGES, PLUGIN_VERSION, SMOKE_INPUT

    block = catalog_block(PLUGIN_VERSION, ASSETS, PACKAGES, SMOKE_INPUT)
    if newline != "\n":
        block = block.replace("\n", newline)
    write_text_exact(
        catalog,
        base.rstrip("\r\n") + newline + newline + block.lstrip("\r\n"),
    )

    # generate the manifests with the plugin's own builder
    completed = subprocess.run(
        [sys.executable, str(plugin_dir / "benchmark.py")],
        capture_output=True, text=True, cwd=str(plugin_dir), check=False,
    )
    report = {
        "installed": True,
        "platform": str(root),
        "checkpoint_sha256": hashlib.sha256(target_checkpoint.read_bytes()).hexdigest(),
        "ranker_root": str(RANKER_ROOT),
        "manifest_builder_exit": completed.returncode,
        "manifest_builder_stdout": completed.stdout.strip(),
        "manifest_builder_stderr": completed.stderr.strip()[-500:],
    }
    print(json.dumps(report, indent=2))
    return 0


def uninstall(root: Path) -> int:
    plugin_runtime, plugin_dir, catalog = locate_platform(root)

    text = read_text_exact(catalog)
    stripped, had_block = strip_block(text)
    if had_block:
        write_text_exact(catalog, stripped)

    removed = []
    for path in (plugin_runtime / "biocandidate.py",):
        if path.is_file():
            path.unlink()
            removed.append(str(path))
    if plugin_dir.exists():
        shutil.rmtree(plugin_dir)
        removed.append(str(plugin_dir))

    print(json.dumps({
        "uninstalled": True,
        "catalog_block_removed": had_block,
        "removed": removed,
    }, indent=2))
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--platform", type=Path, required=True,
                        help="path to Integrated_Yeast_MetaTwin_Deployment")
    parser.add_argument("--uninstall", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)

    root = args.platform.expanduser().resolve()
    return uninstall(root) if args.uninstall else install(root, args.dry_run)


if __name__ == "__main__":
    raise SystemExit(main())
