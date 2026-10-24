"""Contract tests for the Yeast-MetaTwin Platform plugin.

The point of these tests is that they exercise the adapter against the
platform's REAL schema and catalog modules, not a stub: the staging fixture
copies that project's plugin_runtime package, so a contract change there fails
here instead of at the platform.

Skipped when the platform checkout is absent.

中文：Yeast-MetaTwin 平台插件契约测试：对平台真实 schema 与 catalog 模块（非桩）执行，staging 夹具复制其 plugin_runtime 包。
"""

from __future__ import annotations

import importlib.util
import json
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
PLUGIN_DIR = ROOT / "integrations/yeast_metatwin_plugin"
PLATFORM = Path(r"D:\biological\Metabolic model prediction\Integrated_Yeast_MetaTwin_Deployment")
PLATFORM_RUNTIME = PLATFORM / "04_prediction_plugins/plugin_runtime"
CHECKPOINT = ROOT / "artifacts/esm2-t6-opt-seed7/best.pt"

requires_checkpoint = pytest.mark.skipif(
    not CHECKPOINT.is_file(),
    reason="the prepared workspace checkpoint is required",
)

requires_platform = pytest.mark.skipif(
    not PLATFORM_RUNTIME.is_dir(),
    reason="the Yeast-MetaTwin Platform checkout is required",
)
requires_checkpoint = pytest.mark.skipif(
    not CHECKPOINT.is_file(),
    reason="the local BCR checkpoint under gitignored artifacts/ is required",
)


def _load_plain(module_name: str, path: Path):
    """Import a module that has no platform-relative imports."""
    spec = importlib.util.spec_from_file_location(module_name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def staged_platform(tmp_path_factory):
    """A copy of the platform's plugin package with the adapter dropped in."""
    stage = tmp_path_factory.mktemp("platform")
    runtime = stage / "04_prediction_plugins" / "plugin_runtime"
    runtime.mkdir(parents=True)
    for name in ("__init__.py", "schema.py", "catalog.py"):
        shutil.copy(PLATFORM_RUNTIME / name, runtime / name)
    shutil.copy(PLUGIN_DIR / "biocandidate.py", runtime / "biocandidate.py")
    # catalog_entry uses relative imports, so it is validated from inside the package
    shutil.copy(PLUGIN_DIR / "catalog_entry.py", runtime / "catalog_entry.py")
    _pristine_catalog(runtime)

    plugin_root = stage / "04_prediction_plugins"
    assets = plugin_root / "BioCandidateRanker" / "assets"
    assets.mkdir(parents=True)
    if CHECKPOINT.is_file():
        shutil.copy(CHECKPOINT, assets / "esm2-t6-opt-seed7.best.pt")
    # the pointer the installer writes; without it the adapter falls back to an
    # importable biocandidate, which would hide a wrong plugin_root convention
    (plugin_root / "BioCandidateRanker" / "ranker_root.txt").write_text(
        str(ROOT).replace("\\", "/") + "\n", encoding="utf-8"
    )
    return stage, plugin_root


def _pristine_catalog(runtime: Path) -> None:
    """Normalise a copied catalog.py to its pre-plugin state.

    The live platform is expected to already carry the installer's block (that
    is the integration), so the tests must start from the block-free text to
    assert the real invariant: install then uninstall returns to pristine.
    """
    installer = _load_plain("plugin_installer_probe", PLUGIN_DIR / "install.py")
    catalog = runtime / "catalog.py"
    stripped, _ = installer.strip_block(installer.read_text_exact(catalog))
    installer.write_text_exact(catalog, stripped)


def _import_runtime(stage: Path):
    """Import the staged platform package, never the real checkout.

    Any other path that could provide plugin_runtime is removed first, so a
    test cannot accidentally validate against the live platform instead of the
    staged copy.
    """
    container = stage / "04_prediction_plugins"
    sys.path[:] = [
        entry for entry in sys.path
        if not (entry.endswith("04_prediction_plugins") and Path(entry) != container)
    ]
    if str(container) in sys.path:
        sys.path.remove(str(container))
    sys.path.insert(0, str(container))
    for name in list(sys.modules):
        if name.startswith("plugin_runtime"):
            del sys.modules[name]


def _adapter(stage: Path):
    _import_runtime(stage)
    from plugin_runtime import biocandidate, schema

    return biocandidate, schema


# --- catalog data ---------------------------------------------------------------


def test_catalog_entry_is_well_formed_without_the_platform():
    entry = _load_plain("catalog_entry_probe", PLUGIN_DIR / "catalog_entry.py")
    assert entry.CAPABILITIES == ("kcat_prediction",)
    assert entry.SMOKE_INPUT["capability"] == "kcat_prediction"
    paths = [path for path, _ in entry.ASSETS]
    assert "BioCandidateRanker/assets/esm2-t6-opt-seed7.best.pt" in paths
    assert "plugin_runtime/biocandidate.py" in paths
    distributions = [name for name, _, _ in entry.PACKAGES]
    assert "torch" in distributions and "fair-esm" in distributions


@requires_platform
def test_catalog_entry_converts_to_platform_dataclasses(staged_platform):
    stage, _ = staged_platform
    _import_runtime(stage)
    import plugin_runtime.catalog_entry as probe  # type: ignore[import-not-found]

    spec = probe.to_plugin_spec()
    assert spec.name == "BioCandidateRanker"
    assert spec.capabilities == ("kcat_prediction",)
    assert spec.required_assets and spec.packages


# --- input handling -------------------------------------------------------------


@requires_platform
def test_prepare_input_reports_transforms_and_context(staged_platform):
    stage, _ = staged_platform
    adapter, schema = _adapter(stage)

    canonical = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQANL"
    _, context, messages, transforms = adapter.prepare_input(
        schema.PluginInput(request_id="a", capability="kcat_prediction",
                           sequence=canonical, substrate_smiles="CC(=O)O")
    )
    assert [item.field for item in transforms] == ["sequence"]
    assert transforms[0].truncated is False
    assert context["enzyme_type"] == "unknown"
    assert any("defaults" in message for message in messages)

    _, _, _, with_substitutions = adapter.prepare_input(
        schema.PluginInput(request_id="b", capability="kcat_prediction",
                           sequence="MKTAY" + "XZ" + "IAK", substrate_smiles="CC(=O)O")
    )
    non_canonical = [item for item in with_substitutions if "non_canonical" in item.field]
    assert non_canonical and non_canonical[0].strategy == "mapped_to_unknown_token(2)"

    _, _, _, long_transforms = adapter.prepare_input(
        schema.PluginInput(request_id="c", capability="kcat_prediction",
                           sequence="ACDEFGHIKLMNPQRSTVWY" * 30, substrate_smiles="CC(=O)O")
    )
    assert long_transforms[0].truncated is True
    assert long_transforms[0].effective_length == adapter.SEQUENCE_LIMIT
    assert long_transforms[0].strategy == "head_first_512"


@requires_platform
def test_all_canonical_sequence_reports_no_substitutions(staged_platform):
    """Regression: counting distinct residues instead of positions flagged every sequence."""
    stage, _ = staged_platform
    adapter, schema = _adapter(stage)
    _, _, _, transforms = adapter.prepare_input(
        schema.PluginInput(request_id="d", capability="kcat_prediction",
                           sequence="MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQANL",
                           substrate_smiles="CC(=O)O")
    )
    assert not any("non_canonical" in item.field for item in transforms)


@requires_platform
def test_platform_schema_rejects_blank_sequence_before_the_adapter():
    """The platform contract guards blank input, so the adapter need not.

    Runs in a subprocess so importing the live platform here cannot shadow the
    staged copy that the other tests use.
    """
    script = "\n".join([
        f"import sys; sys.path.insert(0, r'{PLATFORM_RUNTIME.parent}')",
        "from plugin_runtime.schema import PluginInput",
        "try:",
        "    PluginInput(request_id='x', capability='kcat_prediction',"
        " sequence='', substrate_smiles='CC(=O)O')",
        "except ValueError as exc:",
        "    print('REJECTED:', exc)",
        "else:",
        "    print('ACCEPTED')",
        "",
    ])
    completed = subprocess.run([sys.executable, "-c", script],
                               capture_output=True, text=True)
    assert "REJECTED" in completed.stdout and "cannot be blank" in completed.stdout


@requires_platform
@requires_checkpoint
@requires_platform
@pytest.mark.parametrize(
    "capability,sequence,expected_status",
    [
        ("km_prediction", "ACDEFG", "unsupported"),      # capability not offered
        ("kcat_prediction", "ACDEF1G", "unsupported"),   # invalid residue
        ("kcat_prediction", "ACDEFG", "ready"),          # control: valid input
    ],
)
def test_bad_input_returns_unsupported_not_an_exception(
    staged_platform, capability, sequence, expected_status
):
    stage, plugin_root = staged_platform
    adapter, schema = _adapter(stage)
    result = adapter.predict(
        schema.PluginInput(request_id="bad", capability=capability,
                           sequence=sequence, substrate_smiles="CC(=O)O"),
        plugin_root=plugin_root,
    )
    assert result.status.value == expected_status
    # the contract forbids predictions on non-ready results; the control case
    # must produce them, which is what proves the rejection paths are specific
    if expected_status == "ready":
        assert result.predictions
    else:
        assert result.predictions == ()


# --- end to end -----------------------------------------------------------------


@requires_platform
@requires_checkpoint
def test_predict_satisfies_the_platform_contract(staged_platform):
    stage, plugin_root = staged_platform
    adapter, schema = _adapter(stage)
    # remove any already-importable biocandidate so the ranker_root.txt pointer
    # is genuinely exercised; the adapter must resolve the checkout itself
    source = str(ROOT / "src")
    saved = sys.path[:]
    sys.path[:] = [entry for entry in sys.path if entry != source]
    for name in list(sys.modules):
        if name == "biocandidate" or name.startswith("biocandidate."):
            del sys.modules[name]
    try:
        result = adapter.predict(
            schema.PluginInput(
                request_id="e2e", capability="kcat_prediction",
                sequence="MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQANL",
                substrate_smiles="CC(=O)O",
                context={"organism": "Saccharomyces cerevisiae", "ec": "1.1.1.1"},
            ),
            plugin_root=plugin_root,
        )
    finally:
        sys.path[:] = saved

    assert result.status.value == "ready"
    names = [prediction.name for prediction in result.predictions]
    assert names == ["log10_kcat", "kcat"]
    linear = next(p for p in result.predictions if p.name == "kcat")
    log10 = next(p for p in result.predictions if p.name == "log10_kcat")
    assert linear.value == pytest.approx(10 ** log10.value, rel=1e-9)
    assert linear.unit == "s^-1" and linear.scale == "linear"
    assert log10.scale == "log10"

    # Honest reporting the platform relies on
    assert result.applicability.assessment_status == "not_available"
    assert result.applicability.ood is None
    assert result.uncertainty.status == "available"
    assert result.uncertainty.calibrated is False
    assert result.uncertainty.tree_member_count is None
    assert result.benchmark.status == "not_run"
    assert result.benchmark.public_reference is False
    assert any("only log10_kcat is trained" in message.lower() for message in result.messages)
    # the serialised payload must be JSON-safe
    assert json.loads(json.dumps(result.to_dict()))["capability"] == "kcat_prediction"


# --- manifests ------------------------------------------------------------------


def test_benchmark_manifests_are_blocked_and_honest(tmp_path):
    module = _load_plain("plugin_benchmark_probe", PLUGIN_DIR / "benchmark.py")
    module.ROOT = tmp_path
    benchmark, readiness = module.build_manifests()
    module.validate_manifest(benchmark)
    module.validate_manifest(readiness)

    assert benchmark["status"] == "blocked"
    assert benchmark["independent"] is False
    assert benchmark["metrics_published"] is False
    assert benchmark["metrics"] == {}
    dataset = benchmark["datasets"][0]
    assert dataset["training_overlap"] is True
    assert dataset["independent_benchmark_eligible"] is False
    assert len(dataset["sha256"]) == 64

    assert readiness["status"] == "blocked"
    assert readiness["output_contract"]["ood_field"] == "applicability.ood"
    assert readiness["output_contract"]["ood_status"] == "not_implemented"
    assert readiness["output_contract"]["tree_member_count_field"] is None
    assert readiness["output_contract"]["uncertainty_calibrated"] is False
    assert readiness["trained_tasks"] == ["log10_kcat"]
    assert set(readiness["untrained_declared_tasks"]) == {
        "log10_km", "log10_activity", "log10_flux"
    }


# --- installer ------------------------------------------------------------------


@requires_checkpoint
@requires_platform
def test_installer_is_idempotent_and_uninstall_is_byte_exact(tmp_path):
    stage = tmp_path / "platform"
    runtime = stage / "04_prediction_plugins" / "plugin_runtime"
    runtime.mkdir(parents=True)
    for name in ("__init__.py", "schema.py", "catalog.py"):
        shutil.copy(PLATFORM_RUNTIME / name, runtime / name)
    _pristine_catalog(runtime)
    catalog = runtime / "catalog.py"
    before = catalog.read_bytes()

    installer = PLUGIN_DIR / "install.py"
    install = [sys.executable, str(installer), "--platform", str(stage)]

    first = subprocess.run(install, capture_output=True, text=True)
    assert first.returncode == 0, first.stderr[-800:]
    second = subprocess.run(install, capture_output=True, text=True)
    assert second.returncode == 0, second.stderr[-800:]

    plugin_dir = stage / "04_prediction_plugins" / "BioCandidateRanker"
    assert (runtime / "biocandidate.py").is_file()
    assert (plugin_dir / "readiness_manifest.json").is_file()
    # manifests must stay LF: Path.write_text would CRLF them on Windows
    assert b"\r\n" not in (plugin_dir / "readiness_manifest.json").read_bytes()
    assert b"\r\n" not in catalog.read_bytes(), "catalog line endings changed"

    removed = subprocess.run(install + ["--uninstall"], capture_output=True, text=True)
    assert removed.returncode == 0, removed.stderr[-800:]
    assert catalog.read_bytes() == before, "uninstall did not restore catalog.py byte-exactly"
    assert not plugin_dir.exists()
    assert not (runtime / "biocandidate.py").exists()


@requires_checkpoint
@requires_platform
def test_dry_run_writes_nothing(tmp_path):
    stage = tmp_path / "platform"
    runtime = stage / "04_prediction_plugins" / "plugin_runtime"
    runtime.mkdir(parents=True)
    for name in ("__init__.py", "schema.py", "catalog.py"):
        shutil.copy(PLATFORM_RUNTIME / name, runtime / name)
    before = (runtime / "catalog.py").read_bytes()

    completed = subprocess.run(
        [sys.executable, str(PLUGIN_DIR / "install.py"),
         "--platform", str(stage), "--dry-run"],
        capture_output=True, text=True,
    )
    assert completed.returncode == 0
    assert json.loads(completed.stdout)["dry_run"] is True
    assert (runtime / "catalog.py").read_bytes() == before
    assert not (stage / "04_prediction_plugins" / "BioCandidateRanker").exists()
