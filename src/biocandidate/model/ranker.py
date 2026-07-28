from __future__ import annotations

import torch
import torch.nn as nn

from ..config import ModelConfig
from .encoders import (
    ContextEncoder,
    ESM2ProteinEncoder,
    GlobalMeanProteinEncoder,
    ProteinChunkEncoder,
    SparseMoleculeEncoder,
)


class BioCandidateRanker(nn.Module):
    """Task-query fusion over protein, molecular graph, and host context."""

    def __init__(self, config: ModelConfig) -> None:
        super().__init__()
        self.config = config
        if config.protein_encoder == "chunk_transformer":
            self.protein = ProteinChunkEncoder(
                config.d_model, config.num_heads, config.protein_layers,
                config.protein_chunk_size, config.dropout)
        elif config.protein_encoder == "esm2":
            self.protein = ESM2ProteinEncoder(
                config.d_model, config.protein_chunk_size, config.dropout,
                esm2_model_name=config.esm2_model_name,
                esm2_frozen=config.esm2_frozen)
        else:
            self.protein = GlobalMeanProteinEncoder(config.d_model)
        self.molecule = SparseMoleculeEncoder(
            config.d_model, config.molecule_layers, config.dropout)
        self.context = ContextEncoder(
            config.context_buckets, config.context_fields,
            config.fba_context_dim, config.d_model)
        if config.fusion_mode == "task_query":
            self.modality_embedding = nn.Embedding(3, config.d_model)
            query_count = 1 if config.shared_task_query else len(config.task_names)
            self.task_queries = nn.Parameter(torch.randn(query_count, config.d_model) * 0.02)
            self.cross_attention = nn.MultiheadAttention(
                config.d_model, config.num_heads, dropout=config.dropout, batch_first=True)
            query_layer = nn.TransformerEncoderLayer(
                d_model=config.d_model,
                nhead=config.num_heads,
                dim_feedforward=config.d_model * 4,
                dropout=config.dropout,
                activation="gelu",
                batch_first=True,
                norm_first=True,
            )
            self.query_fusion = nn.TransformerEncoder(
                query_layer, config.fusion_layers, enable_nested_tensor=False)
            self.output_norm = nn.LayerNorm(config.d_model)
        else:
            self.late_fusion = nn.Sequential(
                nn.Linear(3 * config.d_model, config.d_model),
                nn.SiLU(),
                nn.Dropout(config.dropout),
                nn.Linear(config.d_model, config.d_model),
                nn.LayerNorm(config.d_model),
            )
        output_dimension = 2 if config.uncertainty_mode == "heteroscedastic" else 1
        self.heads = nn.ModuleList(
            nn.Sequential(
                nn.Linear(config.d_model, config.d_model // 2),
                nn.SiLU(),
                nn.Dropout(config.dropout),
                nn.Linear(config.d_model // 2, output_dimension),
            )
            for _ in config.task_names
        )

    def forward(self, batch: dict[str, torch.Tensor]) -> dict[str, torch.Tensor]:
        protein_tokens, protein_mask = self.protein(
            batch["sequence_tokens"], batch["sequence_mask"])
        batch_size = protein_tokens.shape[0]
        molecule_tokens = self.molecule(
            batch["atom_features"], batch["edge_index"], batch["edge_features"],
            batch["graph_batch"], batch_size)
        context_tokens, context_mask = self.context(
            batch["context_ids"], batch["fba_context"], batch["fba_context_mask"])

        if not self.config.use_protein:
            protein_tokens = torch.zeros_like(protein_tokens)
            protein_mask = torch.zeros_like(protein_mask)
        if not self.config.use_molecule:
            molecule_tokens = torch.zeros_like(molecule_tokens)
        if not self.config.use_context:
            context_tokens = torch.zeros_like(context_tokens)
            context_mask = torch.zeros_like(context_mask)

        if self.config.fusion_mode == "late_concat":
            protein_denominator = protein_mask.sum(dim=1, keepdim=True).clamp_min(1)
            protein_pooled = (
                protein_tokens * protein_mask.unsqueeze(-1)).sum(dim=1) / protein_denominator
            context_denominator = context_mask.sum(dim=1, keepdim=True).clamp_min(1)
            context_pooled = (
                context_tokens * context_mask.unsqueeze(-1)).sum(dim=1) / context_denominator
            fused = self.late_fusion(torch.cat(
                (protein_pooled, molecule_tokens[:, 0], context_pooled), dim=-1))
            predictions = torch.stack([head(fused) for head in self.heads], dim=1)
            mean = predictions[..., 0]
            log_variance = (
                predictions[..., 1].clamp(-8.0, 8.0)
                if self.config.uncertainty_mode == "heteroscedastic"
                else torch.zeros_like(mean)
            )
            return {
                "mean": mean,
                "log_variance": log_variance,
                "standard_deviation": torch.exp(0.5 * log_variance),
            }

        protein_tokens = protein_tokens + self.modality_embedding.weight[0]
        molecule_tokens = molecule_tokens + self.modality_embedding.weight[1]
        context_tokens = context_tokens + self.modality_embedding.weight[2]
        modalities = torch.cat((protein_tokens, molecule_tokens, context_tokens), dim=1)
        modality_mask = torch.cat((
            protein_mask,
            torch.full(
                (batch_size, 1), self.config.use_molecule,
                device=protein_mask.device, dtype=torch.bool),
            context_mask,
        ), dim=1)

        queries = self.task_queries.unsqueeze(0).expand(batch_size, -1, -1)
        if self.config.shared_task_query:
            queries = queries.expand(-1, len(self.config.task_names), -1)
        attended, _ = self.cross_attention(
            queries, modalities, modalities, key_padding_mask=~modality_mask,
            need_weights=False)
        queries = self.output_norm(self.query_fusion(queries + attended))
        predictions = torch.stack(
            [head(queries[:, index]) for index, head in enumerate(self.heads)], dim=1)
        mean = predictions[..., 0]
        log_variance = (
            predictions[..., 1].clamp(-8.0, 8.0)
            if self.config.uncertainty_mode == "heteroscedastic"
            else torch.zeros_like(mean)
        )
        return {
            "mean": mean,
            "log_variance": log_variance,
            "standard_deviation": torch.exp(0.5 * log_variance),
        }
