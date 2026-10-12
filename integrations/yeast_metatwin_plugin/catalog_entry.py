"""The catalog entry that registers BioCandidateRanker with the platform.

Kept as plain data so `install.py` can emit it without importing the
platform's package, and so the change to catalog.py is reviewable as data
rather than as an in-place edit. `to_plugin_spec()` converts the data using the
platform's own dataclasses when they are importable.
"""

from __future__ import annotations

from typing import Any


PLUGIN_NAME = "BioCandidateRanker"
PLUGIN_VERSION = "software-implementation-v6+esm2-t6-seed7"
CAPABILITIES = ("kcat_prediction",)
SOURCE_URL = "https://github.com/liubuing/BioCandidateRanker"

# (relative_path, min_size_bytes) tuples, mirroring the platform's AssetSpec.
ASSETS: tuple[tuple[str, int], ...] = (
    # The checkpoint is the whole model. The adapter is tracked too, so a
    # drifted adapter fails the same hash gate as a drifted weight file.
    ("BioCandidateRanker/assets/esm2-t6-opt-seed7.best.pt", 90_000_000),
    ("plugin_runtime/biocandidate.py", 1_000),
)

# (distribution, minimum, maximum_exclusive) triples, mirroring PackageSpec.
PACKAGES: tuple[tuple[str, tuple[int, ...] | None, tuple[int, ...] | None], ...] = (
    # torch 2.6+ is supported through the adapter's safe_globals shim, so no
    # upper bound is declared. fair-esm supplies the frozen ESM-2 backbone.
    ("torch", (2, 0), None),
    ("numpy", (1, 23), None),
    ("rdkit", None, None),
    ("fair-esm", (2, 0), None),
)

SMOKE_INPUT: dict[str, Any] = {
    "request_id": "smoke-biocandidate-001",
    "capability": "kcat_prediction",
    "sequence": "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQANL",
    "substrate_smiles": "CC(=O)O",
}


def to_asset_specs():
    """Build the platform's AssetSpec objects; requires the package context."""
    from .catalog import AssetSpec

    return tuple(AssetSpec(path, min_size_bytes=size) for path, size in ASSETS)


def to_package_specs():
    """Build the platform's PackageSpec objects; requires the package context."""
    from .catalog import PackageSpec

    specs = []
    for distribution, minimum, maximum in PACKAGES:
        if minimum and maximum:
            specs.append(PackageSpec(distribution, minimum, maximum))
        elif minimum:
            specs.append(PackageSpec(distribution, minimum))
        else:
            specs.append(PackageSpec(distribution))
    return tuple(specs)


def to_plugin_spec():
    """Build the platform's PluginSpec; requires the package context."""
    from .catalog import PluginSpec

    return PluginSpec(
        name=PLUGIN_NAME,
        plugin_version=PLUGIN_VERSION,
        capabilities=CAPABILITIES,
        required_assets=to_asset_specs(),
        packages=to_package_specs(),
        source_url=SOURCE_URL,
    )
