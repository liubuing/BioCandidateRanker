from __future__ import annotations

from dataclasses import asdict, dataclass


@dataclass(slots=True)
class ModelConfig:
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
        return asdict(self)

    @classmethod
    def from_dict(cls, values: dict) -> "ModelConfig":
        values = dict(values)
        if "task_names" in values:
            values["task_names"] = tuple(values["task_names"])
        return cls(**values)
