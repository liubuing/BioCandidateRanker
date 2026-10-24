"""Canonical record schema and label metadata contracts shared by all adapters.

:class:`EnzymeSubstrateRecord` is the single governed row type flowing into the
model: adapters normalize their sources into it and validation invariants run
in ``__post_init__``. ``TASK_NAMES`` fixes the task column order shared by the
collator, losses, and output heads; the metadata classes keep per-label
canonical units and provenance so predictions stay traceable to evidence.

中文：所有适配器共享的规范化记录 schema 与标签元数据契约。
EnzymeSubstrateRecord 是流入模型的唯一受治理行类型：各适配器把来源规范化为它，
校验不变量在 __post_init__ 中执行。TASK_NAMES 固定由 collator、损失与输出头
共享的任务列顺序；元数据类保存逐标签的规范化单位与出处，使预测可回溯到证据。
"""
from __future__ import annotations

import math
import json
from dataclasses import dataclass
from enum import Enum
from typing import Optional, Tuple


class EvidenceTier(str, Enum):
    """Label provenance tiers, from direct measurement down to unknown.
    中文：标签证据层级，从直接测量到未知。
    """

    DIRECT = "direct"
    CURATED = "curated"
    INFERRED = "inferred"
    UNKNOWN = "unknown"


FBA_FEATURE_SCHEMA_VERSION = 1
CANONICAL_LOG10_FLUX_UNIT = "log10(mmol gDW^-1 h^-1)"


@dataclass(frozen=True, slots=True)
class FBAFeatureMetadata:
    """Identity and ordered feature contract for one governed FBA context vector.
    中文：单个受治理 FBA 上下文向量的身份与有序特征契约。
    """

    schema_version: int
    feature_ids: Tuple[str, ...]
    model_id: str
    solver_id: str
    objective_id: str
    condition_id: str

    def __post_init__(self) -> None:
        if isinstance(self.schema_version, bool) or self.schema_version != FBA_FEATURE_SCHEMA_VERSION:
            raise ValueError(
                f"FBA feature schema_version must be {FBA_FEATURE_SCHEMA_VERSION}")
        feature_ids = tuple(str(value).strip() for value in self.feature_ids)
        if not feature_ids or any(not value for value in feature_ids):
            raise ValueError("FBA feature_ids must contain non-empty identifiers")
        if len({value.casefold() for value in feature_ids}) != len(feature_ids):
            raise ValueError("FBA feature_ids must be unique")
        object.__setattr__(self, "feature_ids", feature_ids)
        for name in ("model_id", "solver_id", "objective_id", "condition_id"):
            value = getattr(self, name)
            if not isinstance(value, str) or not value.strip():
                raise ValueError(f"FBA {name} must be a non-empty string")
            object.__setattr__(self, name, value.strip())


class FluxLabelType(str, Enum):
    """Origin of a flux label: simulated (FBA) or experimental.
    中文：通量标签来源：模拟（FBA）或实验。
    """

    SIMULATED = "simulated"
    EXPERIMENTAL = "experimental"


@dataclass(frozen=True, slots=True)
class FluxLabelMetadata:
    """Provenance and semantics for a supplied normalized flux target.
    中文：所供规范化通量目标的出处与语义。
    """

    label_type: FluxLabelType
    target_reaction_id: str
    canonical_unit: str
    provenance_json: str

    def __post_init__(self) -> None:
        if not isinstance(self.label_type, FluxLabelType):
            object.__setattr__(self, "label_type", FluxLabelType(self.label_type))
        for name in ("target_reaction_id", "canonical_unit"):
            value = getattr(self, name)
            if not isinstance(value, str) or not value.strip():
                raise ValueError(f"flux label {name} must be a non-empty string")
            object.__setattr__(self, name, value.strip())
        if self.canonical_unit != CANONICAL_LOG10_FLUX_UNIT:
            raise ValueError(
                f"flux label canonical_unit must be {CANONICAL_LOG10_FLUX_UNIT!r}")
        try:
            provenance = json.loads(self.provenance_json)
        except (TypeError, json.JSONDecodeError) as exc:
            raise ValueError("flux label provenance_json must contain valid JSON") from exc
        if not isinstance(provenance, dict) or not provenance:
            raise ValueError("flux label provenance must be a non-empty object")
        object.__setattr__(
            self, "provenance_json", json.dumps(provenance, sort_keys=True, separators=(",", ":")))

    @property
    def provenance(self) -> dict:
        return json.loads(self.provenance_json)


@dataclass(frozen=True, slots=True)
class KineticLabelMetadata:
    """Canonical unit, evidence tier, and provenance attached to one kinetic label.
    中文：附在单个动力学标签上的规范化单位、证据层级与出处。
    """

    task: str
    canonical_unit: str
    evidence_tier: EvidenceTier
    provenance_json: str

    def __post_init__(self) -> None:
        if not self.task.strip():
            raise ValueError("kinetic label task must not be empty")
        if not self.canonical_unit.strip():
            raise ValueError("kinetic label canonical_unit must not be empty")
        if not isinstance(self.evidence_tier, EvidenceTier):
            object.__setattr__(self, "evidence_tier", EvidenceTier(self.evidence_tier))
        try:
            provenance = json.loads(self.provenance_json)
        except (TypeError, json.JSONDecodeError) as exc:
            raise ValueError("kinetic label provenance_json must contain valid JSON") from exc
        if not isinstance(provenance, dict) or not provenance:
            raise ValueError("kinetic label provenance must be a non-empty object")
        object.__setattr__(self, "task", self.task.strip())
        object.__setattr__(self, "canonical_unit", self.canonical_unit.strip())
        object.__setattr__(
            self, "provenance_json", json.dumps(provenance, sort_keys=True, separators=(",", ":")))

    @property
    def provenance(self) -> dict:
        return json.loads(self.provenance_json)


@dataclass(frozen=True, slots=True)
class EnzymeSubstrateRecord:
    """One governed enzyme-substrate candidate row with sparse multi-task labels.

    ``None`` labels mean "not measured" and are masked out by the collator.
    Optional FBA context and flux/kinetic metadata each carry their own
    validation invariants enforced below.
    中文：一条带稀疏多任务标签的受治理酶-底物候选行。None 标签表示"未测量"，
    由 collator 掩掉。可选的 FBA 上下文与通量/动力学元数据各自携带在下方强制
    执行的校验不变量。
    """

    sequence: str
    substrate_smiles: str
    organism: str
    ec: str
    enzyme_type: str
    candidate_id: str = ""
    substrate_name: str = ""
    reaction: str = ""
    fba_context: Tuple[float, ...] = ()
    fba_feature_metadata: Optional[FBAFeatureMetadata] = None
    log10_kcat: Optional[float] = None
    log10_km: Optional[float] = None
    log10_kcat_per_km: Optional[float] = None
    log10_activity: Optional[float] = None
    log10_flux: Optional[float] = None
    evidence_tier: EvidenceTier = EvidenceTier.UNKNOWN
    source_dataset: str = ""
    source_row: int = 0
    source_rows: Tuple[int, ...] = ()
    replicate_count: int = 1
    label_stddev: Optional[float] = None
    split: Optional[str] = None
    activity_rank: Optional[float] = None
    campaign_group: str = ""
    kinetic_label_metadata: Tuple[KineticLabelMetadata, ...] = ()
    flux_label_metadata: Optional[FluxLabelMetadata] = None

    def __post_init__(self) -> None:
        sequence = "".join(self.sequence.split()).upper()
        if not sequence:
            raise ValueError("sequence must not be empty")
        if not self.substrate_smiles.strip():
            raise ValueError("substrate_smiles must not be empty")
        if self.source_row < 0:
            raise ValueError("source_row must be non-negative")
        if self.replicate_count < 1:
            raise ValueError("replicate_count must be positive")

        for name in (
            "log10_kcat", "log10_km", "log10_kcat_per_km", "log10_activity", "log10_flux",
        ):
            value = getattr(self, name)
            if value is not None and not math.isfinite(value):
                raise ValueError(f"{name} must be finite when provided")
        if self.activity_rank is not None and (
                not math.isfinite(self.activity_rank) or not 0.0 <= self.activity_rank <= 1.0):
            raise ValueError("activity_rank must be finite and between 0 and 1 when provided")

        object.__setattr__(self, "sequence", sequence)
        object.__setattr__(self, "substrate_smiles", self.substrate_smiles.strip())
        object.__setattr__(self, "organism", self.organism.strip())
        object.__setattr__(self, "ec", self.ec.strip())
        object.__setattr__(self, "enzyme_type", self.enzyme_type.strip())
        object.__setattr__(self, "candidate_id", self.candidate_id.strip())
        object.__setattr__(self, "substrate_name", self.substrate_name.strip())
        object.__setattr__(self, "reaction", self.reaction.strip())
        if any(not math.isfinite(value) for value in self.fba_context):
            raise ValueError("fba_context values must be finite")
        object.__setattr__(self, "fba_context", tuple(float(value) for value in self.fba_context))
        if self.fba_context:
            if not isinstance(self.fba_feature_metadata, FBAFeatureMetadata):
                raise ValueError("non-empty fba_context requires FBAFeatureMetadata")
            if len(self.fba_context) != len(self.fba_feature_metadata.feature_ids):
                raise ValueError("fba_context width must match FBA feature_ids exactly")
        elif self.fba_feature_metadata is not None:
            raise ValueError("FBAFeatureMetadata requires a non-empty fba_context")
        object.__setattr__(self, "source_dataset", self.source_dataset.strip())
        object.__setattr__(self, "campaign_group", self.campaign_group.strip())
        source_rows = self.source_rows or (self.source_row,)
        if any(row < 0 for row in source_rows):
            raise ValueError("source_rows must be non-negative")
        object.__setattr__(self, "source_rows", tuple(source_rows))
        if self.label_stddev is not None and (
                not math.isfinite(self.label_stddev) or self.label_stddev < 0):
            raise ValueError("label_stddev must be finite and non-negative")
        if not isinstance(self.evidence_tier, EvidenceTier):
            object.__setattr__(self, "evidence_tier", EvidenceTier(self.evidence_tier))
        metadata = tuple(self.kinetic_label_metadata)
        if any(not isinstance(item, KineticLabelMetadata) for item in metadata):
            raise ValueError("kinetic_label_metadata must contain KineticLabelMetadata values")
        metadata_tasks = [item.task for item in metadata]
        if len(metadata_tasks) != len(set(metadata_tasks)):
            raise ValueError("kinetic_label_metadata tasks must be unique")
        object.__setattr__(self, "kinetic_label_metadata", metadata)
        if self.log10_flux is not None:
            if not isinstance(self.flux_label_metadata, FluxLabelMetadata):
                raise ValueError("log10_flux requires FluxLabelMetadata")
            if self.fba_feature_metadata is not None:
                target_id = self.flux_label_metadata.target_reaction_id.casefold()
                input_ids = {value.casefold() for value in self.fba_feature_metadata.feature_ids}
                if target_id in input_ids:
                    raise ValueError("flux target reaction must not be included in FBA feature_ids")
        elif self.flux_label_metadata is not None:
            raise ValueError("FluxLabelMetadata requires a log10_flux value")


# Canonical task order; fixes the label column layout shared by the collator,
# the multitask losses, and the per-task output heads. Do not reorder.
TASK_NAMES = (
    "log10_kcat",
    "log10_km",
    "log10_kcat_per_km",
    "log10_activity",
    "log10_flux",
)
