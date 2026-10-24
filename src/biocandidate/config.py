"""Model architecture configuration shared by training, inference, and checkpoints.

``ModelConfig`` is the single source of truth for every architecture and
task-head knob. It is serialized verbatim into each checkpoint, and
``load_checkpoint`` rebuilds the model from that payload, so any field added
here needs a stable default to keep old checkpoints loadable.

中文：训练、推理和检查点共享的模型架构配置。``ModelConfig`` 是所有架构与任务头
参数的唯一事实来源。它会被逐字序列化进每个检查点，``load_checkpoint`` 依据该
载荷重建模型，因此此处新增的任何字段都必须提供稳定默认值，以保证旧检查点
仍可加载。
"""

from __future__ import annotations

from dataclasses import asdict, dataclass


@dataclass(slots=True)
class ModelConfig:
    """Architecture, modality, and task-head settings for :class:`BioCandidateRanker`.

    中文：:class:`BioCandidateRanker` 的架构、模态与任务头设置。
    """

    d_model: int = 128
    num_heads: int = 8
    protein_layers: int = 3
    protein_chunk_size: int = 128
    molecule_layers: int = 4
    fusion_layers: int = 2
    context_buckets: int = 4096
    context_fields: int = 4
    fba_context_dim: int = 8
    dropout: float = 0.1
    use_protein: bool = True
    use_molecule: bool = True
    use_context: bool = True
    shared_task_query: bool = False
    fusion_mode: str = "task_query"
    protein_encoder: str = "chunk_transformer"
    uncertainty_mode: str = "heteroscedastic"
    esm2_model_name: str = "esm2_t6_8M_UR50D"
    esm2_frozen: bool = True
    task_names: tuple[str, ...] = (
        "log10_kcat",
        "log10_km",
        "log10_kcat_per_km",
        "log10_activity",
        "log10_flux",
    )

    def __post_init__(self) -> None:
        if self.d_model % self.num_heads:
            raise ValueError("d_model must be divisible by num_heads")
        if self.protein_chunk_size < 1:
            raise ValueError("protein_chunk_size must be positive")
        if self.fba_context_dim < 1:
            raise ValueError("fba_context_dim must be positive")
        if not self.task_names:
            raise ValueError("at least one task is required")
        if not (self.use_protein or self.use_molecule or self.use_context):
            raise ValueError("at least one modality must be enabled")
        if self.fusion_mode not in {"task_query", "late_concat"}:
            raise ValueError("fusion_mode must be task_query or late_concat")
        if self.protein_encoder not in {"chunk_transformer", "global_mean", "esm2"}:
            raise ValueError(
                "protein_encoder must be chunk_transformer, global_mean, or esm2")
        if self.uncertainty_mode not in {"heteroscedastic", "fixed_variance_mse"}:
            raise ValueError(
                "uncertainty_mode must be heteroscedastic or fixed_variance_mse")

    def to_dict(self) -> dict:
        """Return the configuration as a plain JSON-serializable dict.

        中文：以可直接 JSON 序列化的普通 dict 形式返回本配置。
        """
        return asdict(self)

    @classmethod
    def from_dict(cls, values: dict) -> "ModelConfig":
        """Rebuild a config from :meth:`to_dict` output, coercing task_names back to a tuple.

        中文：依据 :meth:`to_dict` 的输出重建配置，并将 task_names 还原为元组。
        """
        values = dict(values)
        if "task_names" in values:
            values["task_names"] = tuple(values["task_names"])
        return cls(**values)
