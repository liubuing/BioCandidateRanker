from __future__ import annotations

import math

import torch
import torch.nn as nn
import torch.nn.functional as F


class ProteinChunkEncoder(nn.Module):
    """Local Transformer with O(L * chunk_size) sequence attention."""

    def __init__(self, d_model: int, num_heads: int, layers: int,
                 chunk_size: int, dropout: float) -> None:
        super().__init__()
        self.chunk_size = chunk_size
        self.token_embedding = nn.Embedding(22, d_model, padding_idx=0)
        self.position_embedding = nn.Embedding(chunk_size, d_model)
        layer = nn.TransformerEncoderLayer(
            d_model=d_model,
            nhead=num_heads,
            dim_feedforward=d_model * 4,
            dropout=dropout,
            activation="gelu",
            batch_first=True,
            norm_first=True,
        )
        self.encoder = nn.TransformerEncoder(layer, layers, enable_nested_tensor=False)
        self.output_norm = nn.LayerNorm(d_model)

    def forward(self, tokens: torch.Tensor, mask: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        batch_size, length = tokens.shape
        chunk_count = max(1, math.ceil(length / self.chunk_size))
        padded_length = chunk_count * self.chunk_size
        tokens = F.pad(tokens, (0, padded_length - length))
        mask = F.pad(mask, (0, padded_length - length))
        tokens = tokens.reshape(batch_size * chunk_count, self.chunk_size)
        token_mask = mask.reshape(batch_size * chunk_count, self.chunk_size)
        chunk_mask = token_mask.any(dim=1)

        # Transformer cannot consume rows where every key is padded.
        safe_mask = token_mask.clone()
        safe_mask[~chunk_mask, 0] = True
        positions = torch.arange(self.chunk_size, device=tokens.device)
        x = self.token_embedding(tokens) + self.position_embedding(positions).unsqueeze(0)
        x = self.encoder(x, src_key_padding_mask=~safe_mask)
        weights = token_mask.unsqueeze(-1).to(x.dtype)
        pooled = (x * weights).sum(dim=1) / weights.sum(dim=1).clamp_min(1)
        pooled = self.output_norm(pooled) * chunk_mask.unsqueeze(-1)
        return (
            pooled.reshape(batch_size, chunk_count, -1),
            chunk_mask.reshape(batch_size, chunk_count),
        )


class GlobalMeanProteinEncoder(nn.Module):
    """Position-free amino-acid embedding mean used as a low-capacity ablation."""

    def __init__(self, d_model: int) -> None:
        super().__init__()
        self.token_embedding = nn.Embedding(22, d_model, padding_idx=0)
        self.output_norm = nn.LayerNorm(d_model)

    def forward(self, tokens: torch.Tensor, mask: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        weights = mask.unsqueeze(-1).to(torch.float32)
        pooled = (self.token_embedding(tokens) * weights).sum(dim=1)
        pooled = pooled / weights.sum(dim=1).clamp_min(1)
        pooled = self.output_norm(pooled)
        return pooled.unsqueeze(1), mask.any(dim=1, keepdim=True)


class ESM2ProteinEncoder(nn.Module):
    """Pretrained ESM-2 backbone with chunk-pooled projection to d_model.

    Accepts the same (tokens, mask) interface as ProteinChunkEncoder where
    tokens use the project vocabulary (0=pad, 1=unknown, 2-21=ACDEFGHIKLMNPQRSTVWY).
    Internally remaps to ESM-2 vocabulary, runs the frozen backbone, and pools
    per-residue embeddings into chunk-level tokens.
    """

    # Project vocab order: A C D E F G H I K L M N P Q R S T V W Y (indices 2-21)
    # ESM-2 standard alphabet indices (from esm.constants.restype_order):
    _ESM2_RESIDUE_IDS = {
        "A": 5, "C": 23, "D": 13, "E": 9, "F": 18, "G": 6, "H": 21, "I": 12,
        "K": 15, "L": 4, "M": 20, "N": 17, "P": 14, "Q": 16, "R": 10, "S": 8,
        "T": 11, "V": 7, "W": 22, "Y": 19,
    }
    _PROJECT_AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"

    def __init__(self, d_model: int, chunk_size: int, dropout: float,
                 esm2_model_name: str = "esm2_t6_8M_UR50D",
                 esm2_frozen: bool = True, max_length: int = 512) -> None:
        super().__init__()
        import esm  # lazy import; requires `pip install fair-esm`

        self.chunk_size = chunk_size
        self.esm2_model_name = esm2_model_name
        self.esm2_frozen = esm2_frozen
        self.max_length = max_length
        self.truncation_count = 0

        loader = getattr(esm.pretrained, esm2_model_name, None)
        if loader is None:
            raise ValueError(
                f"unknown ESM-2 model {esm2_model_name!r}; "
                "available: esm2_t6_8M_UR50D, esm2_t12_35M_UR50D, "
                "esm2_t30_150M_UR50D, esm2_t33_650M_UR50D")
        model, alphabet = loader()
        self.backbone = model
        self.backbone.eval()

        esm2_dim = model.embed_dim
        # Build remapping table: project token ID -> ESM-2 token ID
        # Project: 0=pad, 1=unknown, 2+i = AMINO_ACIDS[i]
        # ESM-2: 1=pad, 3=unknown
        remap = torch.full((22,), 3, dtype=torch.long)  # default: unknown
        remap[0] = 1  # pad -> ESM-2 pad
        for i, amino_acid in enumerate(self._PROJECT_AMINO_ACIDS):
            remap[2 + i] = self._ESM2_RESIDUE_IDS[amino_acid]
        self.register_buffer("token_remap", remap, persistent=False)

        # ESM-2 special tokens: prepend <cls>=0, append <eos>=2
        self.cls_id = 0
        self.eos_id = 2

        self.projection = nn.Sequential(
            nn.Linear(esm2_dim, d_model),
            nn.SiLU(),
            nn.Dropout(dropout),
            nn.LayerNorm(d_model),
        )
        self.output_norm = nn.LayerNorm(d_model)

        if esm2_frozen:
            for param in self.backbone.parameters():
                param.requires_grad = False

    def train(self, mode: bool = True) -> "ESM2ProteinEncoder":
        super().train(mode)
        # Keep backbone in eval mode when frozen (dropout/batchnorm stability)
        if self.esm2_frozen:
            self.backbone.eval()
        return self

    def forward(self, tokens: torch.Tensor, mask: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        batch_size, length = tokens.shape
        device = tokens.device

        truncated = length > self.max_length
        if truncated:
            tokens = tokens[:, :self.max_length]
            mask = mask[:, :self.max_length]
            length = self.max_length
            self.truncation_count += batch_size
        else:
            self.truncation_count += 0

        # Remap project tokens to ESM-2 vocabulary
        esm_tokens = self.token_remap[tokens.clamp(0, 21)]  # [B, L]

        # Give every sequence its own <eos> immediately after its valid residues and pad
        # only after that EOS with the ESM-2 padding token (embedding is fixed at 0). ESM-2
        # uses positional embeddings, and a single shared trailing EOS after batch-wide
        # padding would shift the seen positions of short sequences and make scores depend
        # on batch composition; own-EOS-per-row + pad-token padding minimizes that coupling.
        valid_lengths = mask.sum(dim=1).clamp(max=length)  # [B]
        max_keep = int(valid_lengths.max().clamp(min=1))
        framed = torch.full(
            (batch_size, max_keep + 2),
            1,  # ESM-2 <pad> token
            dtype=torch.long,
            device=device,
        )
        framed[:, 0] = self.cls_id  # <cls>
        for row in range(batch_size):
            keep = int(valid_lengths[row])
            framed[row, 1:1 + keep] = esm_tokens[row, :keep]
            framed[row, keep + 1] = self.eos_id  # this row's own <eos>

        # Run backbone
        with torch.set_grad_enabled(not self.esm2_frozen and self.training):
            backbone_output = self.backbone(
                framed,
                repr_layers=[self.backbone.num_layers],
                return_contacts=False,
            )
        # representations shape: [B, max_keep+2, esm2_dim]
        # column 0 is <cls>; columns 1..keep are residues; keep+1 is this row's <eos>.
        per_row = backbone_output["representations"][self.backbone.num_layers]

        # Align each row's residue columns back onto the truncated token grid; anything
        # beyond the row's own length stays zero (padding).
        residue_repr = torch.zeros(
            (batch_size, length, per_row.shape[-1]),
            dtype=per_row.dtype,
            device=device,
        )
        for row in range(batch_size):
            keep = int(valid_lengths[row])
            residue_repr[row, :keep] = per_row[row, 1:1 + keep]

        # Project to d_model
        projected = self.projection(residue_repr)  # [B, L, d_model]

        # Chunk-pool to match ProteinChunkEncoder output interface
        chunk_count = max(1, math.ceil(length / self.chunk_size))
        padded_length = chunk_count * self.chunk_size
        if padded_length > length:
            projected = F.pad(projected, (0, 0, 0, padded_length - length))
            mask = F.pad(mask, (0, padded_length - length))

        projected = projected.reshape(batch_size, chunk_count, self.chunk_size, -1)
        chunk_mask = mask.reshape(batch_size, chunk_count, self.chunk_size)
        chunk_valid = chunk_mask.any(dim=2)  # [B, chunk_count]

        weights = chunk_mask.unsqueeze(-1).to(projected.dtype)
        pooled = (projected * weights).sum(dim=2) / weights.sum(dim=2).clamp_min(1)
        pooled = self.output_norm(pooled) * chunk_valid.unsqueeze(-1)

        return pooled, chunk_valid


class SparseMessageLayer(nn.Module):
    def __init__(self, d_model: int, dropout: float) -> None:
        super().__init__()
        self.message = nn.Sequential(
            nn.Linear(d_model * 2, d_model),
            nn.SiLU(),
            nn.Linear(d_model, d_model),
        )
        self.update = nn.GRUCell(d_model, d_model)
        self.norm = nn.LayerNorm(d_model)
        self.dropout = nn.Dropout(dropout)

    def forward(self, x: torch.Tensor, edge_index: torch.Tensor,
                edge_embedding: torch.Tensor) -> torch.Tensor:
        if edge_index.numel() == 0:
            return x
        source, target = edge_index
        messages = self.message(torch.cat((x[source], edge_embedding), dim=-1))
        aggregated = torch.zeros_like(x)
        aggregated.index_add_(0, target, messages)
        degree = torch.zeros(x.shape[0], device=x.device, dtype=x.dtype)
        degree.index_add_(0, target, torch.ones_like(target, dtype=x.dtype))
        aggregated = aggregated / degree.clamp_min(1).unsqueeze(-1)
        updated = self.update(aggregated, x)
        return self.norm(x + self.dropout(updated))


class SparseMoleculeEncoder(nn.Module):
    def __init__(self, d_model: int, layers: int, dropout: float) -> None:
        super().__init__()
        dims = (128, 12, 16, 2, 16, 12)
        self.atom_embeddings = nn.ModuleList(nn.Embedding(size, d_model) for size in dims)
        self.bond_type = nn.Embedding(8, d_model)
        self.bond_flags = nn.ModuleList(nn.Embedding(2, d_model) for _ in range(2))
        self.layers = nn.ModuleList(SparseMessageLayer(d_model, dropout) for _ in range(layers))
        self.output_norm = nn.LayerNorm(d_model)

    def forward(self, atom_features: torch.Tensor, edge_index: torch.Tensor,
                edge_features: torch.Tensor, graph_batch: torch.Tensor,
                batch_size: int) -> torch.Tensor:
        features = atom_features.long()
        limits = (127, 11, 15, 1, 15, 11)
        # Formal charge is shifted to keep embedding indices non-negative.
        features = features.clone()
        features[:, 2] += 7
        x = sum(
            embedding(features[:, index].clamp(0, limits[index]))
            for index, embedding in enumerate(self.atom_embeddings)
        )
        edge_features = edge_features.long()
        edge_embedding = self.bond_type(edge_features[:, 0].clamp(0, 7))
        edge_embedding = edge_embedding + self.bond_flags[0](edge_features[:, 1].clamp(0, 1))
        edge_embedding = edge_embedding + self.bond_flags[1](edge_features[:, 2].clamp(0, 1))
        for layer in self.layers:
            x = layer(x, edge_index, edge_embedding)

        pooled = torch.zeros(batch_size, x.shape[-1], device=x.device, dtype=x.dtype)
        pooled.index_add_(0, graph_batch, x)
        counts = torch.zeros(batch_size, device=x.device, dtype=x.dtype)
        counts.index_add_(0, graph_batch, torch.ones_like(graph_batch, dtype=x.dtype))
        return self.output_norm(pooled / counts.clamp_min(1).unsqueeze(-1)).unsqueeze(1)


class ContextEncoder(nn.Module):
    def __init__(self, buckets: int, fields: int, fba_context_dim: int,
                 d_model: int) -> None:
        super().__init__()
        self.value_embedding = nn.Embedding(buckets, d_model, padding_idx=0)
        self.field_embedding = nn.Embedding(fields, d_model)
        self.norm = nn.LayerNorm(d_model)
        self.fba_projection = nn.Sequential(
            nn.Linear(fba_context_dim, d_model),
            nn.SiLU(),
            nn.LayerNorm(d_model),
        )

    def forward(self, context_ids: torch.Tensor, fba_context: torch.Tensor,
                fba_context_mask: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        fields = torch.arange(context_ids.shape[1], device=context_ids.device)
        tokens = self.value_embedding(context_ids) + self.field_embedding(fields).unsqueeze(0)
        mask = context_ids.ne(0)
        tokens = self.norm(tokens) * mask.unsqueeze(-1)
        fba_token = self.fba_projection(fba_context).unsqueeze(1)
        fba_mask = fba_context_mask.unsqueeze(1)
        fba_token = fba_token * fba_mask.unsqueeze(-1)
        return torch.cat((tokens, fba_token), dim=1), torch.cat((mask, fba_mask), dim=1)
