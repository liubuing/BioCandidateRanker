"""Training loops, campaign-aware samplers, and governed checkpointing.

Epoch runners branch on the model's declared ``uncertainty_mode`` and normalize
the loss by observed-label evidence weight, so masked tasks contribute nothing.
Checkpoints embed the model config, RNG states, and the data-manifest
identities; downstream commands verify those identities before reusing a
checkpoint, so a checkpoint can only be evaluated or resumed against the exact
data it was trained on.

中文：训练循环、campaign 感知采样器与受治理的检查点。Epoch 执行器按模型的
uncertainty_mode 分支，损失按有标签的证据权重归一化；检查点内嵌模型配置、
RNG 状态与数据清单身份，下游命令复用前先校验身份。
"""

from __future__ import annotations

import json
import os
import random
import tempfile
from dataclasses import dataclass
from math import ceil
from pathlib import Path
from typing import Iterator, Sequence

import torch
from torch.utils.data import DataLoader, Dataset, Sampler

from .config import ModelConfig
from .data.schema import EnzymeSubstrateRecord
from .model.losses import (
    masked_multitask_gaussian_loss,
    masked_multitask_mse_loss,
    pairwise_logistic_ranking_loss,
)
from .model.ranker import BioCandidateRanker


def create_lr_scheduler(
    optimizer: torch.optim.Optimizer,
    *,
    warmup_epochs: int = 0,
    total_epochs: int = 1,
    min_lr_ratio: float = 0.01,
) -> torch.optim.lr_scheduler.LambdaLR:
    """Cosine annealing with optional linear warmup.
中文：带可选线性预热的余弦退火。
"""
    import math as _math

    def lr_lambda(epoch: int) -> float:
        if epoch < warmup_epochs:
            return (epoch + 1) / max(1, warmup_epochs)
        progress = (epoch - warmup_epochs) / max(1, total_epochs - warmup_epochs)
        cosine_decay = 0.5 * (1.0 + _math.cos(_math.pi * min(progress, 1.0)))
        return max(min_lr_ratio, min_lr_ratio + (1.0 - min_lr_ratio) * cosine_decay)

    return torch.optim.lr_scheduler.LambdaLR(optimizer, lr_lambda)


class RecordDataset(Dataset):
    """In-memory dataset of :class:`EnzymeSubstrateRecord` items.
中文：以 EnzymeSubstrateRecord 为元素的内存数据集。
"""

    def __init__(self, records: list[EnzymeSubstrateRecord] | tuple[EnzymeSubstrateRecord, ...]):
        self.records = tuple(records)

    def __len__(self) -> int:
        return len(self.records)

    def __getitem__(self, index: int) -> EnzymeSubstrateRecord:
        return self.records[index]


class CampaignBatchSampler(Sampler[list[int]]):
    """Fixed-size batches that never mix campaigns (ranking losses are within-campaign).

    Shuffling is deterministic per epoch: the permutation is seeded with
    ``seed + epoch``, so call :meth:`set_epoch` each epoch to vary it.
    
中文：定长且绝不混合 campaign 的批采样器（排序损失在 campaign 内）。
"""

    def __init__(
        self,
        records: Sequence[EnzymeSubstrateRecord] | Dataset,
        batch_size: int,
        *,
        shuffle: bool = False,
        seed: int = 0,
        drop_last: bool = False,
    ) -> None:
        if batch_size < 1:
            raise ValueError("batch_size must be positive")
        self.batch_size = batch_size
        self.shuffle = shuffle
        self.seed = seed
        self.drop_last = drop_last
        self.epoch = 0
        groups: dict[str, list[int]] = {}
        for index in range(len(records)):
            record = records[index]
            if not isinstance(record, EnzymeSubstrateRecord):
                raise TypeError("CampaignBatchSampler requires EnzymeSubstrateRecord items")
            groups.setdefault(record.campaign_group, []).append(index)
        self._groups = tuple((name, tuple(indices)) for name, indices in groups.items())

    def set_epoch(self, epoch: int) -> None:
        self.epoch = epoch

    def __iter__(self) -> Iterator[list[int]]:
        groups = [(name, list(indices)) for name, indices in self._groups]
        if self.shuffle:
            generator = random.Random(self.seed + self.epoch)
            generator.shuffle(groups)
            for _, indices in groups:
                generator.shuffle(indices)
        for _, indices in groups:
            for start in range(0, len(indices), self.batch_size):
                batch = indices[start:start + self.batch_size]
                if len(batch) == self.batch_size or not self.drop_last:
                    yield batch

    def __len__(self) -> int:
        if self.drop_last:
            return sum(len(indices) // self.batch_size for _, indices in self._groups)
        return sum(ceil(len(indices) / self.batch_size) for _, indices in self._groups)


class CampaignGroupedBatchSampler:
    """Yield every record of one campaign together as a single batch.

    Ranking pairs are then complete within each campaign instead of being restricted to
    the sub-chunk produced by a fixed ``batch_size``. Memory cost is quadratic in the
    largest campaign; consumers with very large groups should downsample pairs instead.
    
中文：把一个 campaign 的全部记录作为单批整体产出。
"""

    def __init__(self, records: Sequence[EnzymeSubstrateRecord], *, shuffle: bool = False,
                 seed: int = 42) -> None:
        self.shuffle = shuffle
        self.seed = seed
        self.epoch = 0
        groups: dict[str, list[int]] = {}
        for index, record in enumerate(records):
            if not isinstance(record, EnzymeSubstrateRecord):
                raise TypeError("CampaignGroupedBatchSampler requires EnzymeSubstrateRecord items")
            groups.setdefault(record.campaign_group, []).append(index)
        self._groups = tuple((name, tuple(indices)) for name, indices in groups.items())

    def set_epoch(self, epoch: int) -> None:
        self.epoch = epoch

    def __iter__(self) -> Iterator[list[int]]:
        groups = [(name, list(indices)) for name, indices in self._groups]
        if self.shuffle:
            generator = random.Random(self.seed + self.epoch)
            generator.shuffle(groups)
            for _, indices in groups:
                generator.shuffle(indices)
        for _, indices in groups:
            yield indices

    def __len__(self) -> int:
        return len(self._groups)


def move_batch(batch: dict[str, torch.Tensor], device: torch.device) -> dict[str, torch.Tensor]:
    return {name: value.to(device) for name, value in batch.items()}


@dataclass(slots=True)
class EpochMetrics:
    loss: float
    observed_labels: int


@dataclass(slots=True)
class RankingEpochMetrics:
    loss: float
    pair_count: int


def run_epoch(model: BioCandidateRanker, loader: DataLoader, device: torch.device,
              optimizer: torch.optim.Optimizer | None = None) -> EpochMetrics:
    """One pass over ``loader``; trains when ``optimizer`` is given, else evaluates.

    The reported loss is the evidence-weighted numerator summed over all
    observed labels divided by the summed denominator, matching the objective
    the masked loss functions normalize internally per batch.
    
中文：过一遍 loader；给 optimizer 则训练，否则评估。
"""
    training = optimizer is not None
    model.train(training)
    total_loss_numerator = 0.0
    total_loss_denominator = 0.0
    saw_batch = False
    observed = 0
    for batch in loader:
        saw_batch = True
        batch = move_batch(batch, device)
        with torch.set_grad_enabled(training):
            output = model(batch)
            if model.config.uncertainty_mode == "heteroscedastic":
                loss, _, loss_numerator, loss_denominator = masked_multitask_gaussian_loss(
                    output["mean"], output["log_variance"], batch["labels"],
                    batch["label_mask"], batch["evidence_weight"], return_components=True)
            else:
                loss, _, loss_numerator, loss_denominator = masked_multitask_mse_loss(
                    output["mean"], batch["labels"], batch["label_mask"],
                    batch["evidence_weight"], return_components=True)
            if training:
                optimizer.zero_grad(set_to_none=True)
                loss.backward()
                torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
                optimizer.step()
        total_loss_numerator += float(loss_numerator.detach())
        total_loss_denominator += float(loss_denominator.detach())
        observed += int(batch["label_mask"].sum())
    if not saw_batch:
        raise ValueError("data loader is empty")
    if total_loss_denominator <= 0:
        raise ValueError("data loader contains no observed labels with positive evidence weight")
    return EpochMetrics(total_loss_numerator / total_loss_denominator, observed)


def run_epoch_accum(model: BioCandidateRanker, loader: DataLoader, device: torch.device,
                    optimizer: torch.optim.Optimizer, *,
                    accum_steps: int = 1) -> EpochMetrics:
    """Training epoch with gradient accumulation for memory-constrained GPUs.

    Each accumulation window joins its micro-batches into the same objective as a single
    batch of the combined rows: the loss is ``sum(numerator)/sum(denominator)`` over the
    window with numerator and denominator accumulated *before* per-batch normalization.
    A micro-batch carrying few observed labels therefore contributes exactly as many
    gradient units as its share of the window's total observed-label weight. Each window
    backpropagates once with that window's total denominator; a trailing partial window is
    normalized by its own denominator, not the nominal ``accum_steps``.
    
中文：显存受限时带梯度累积的训练 epoch。
"""
    if accum_steps < 1:
        raise ValueError("accum_steps must be positive")
    model.train()
    total_loss_numerator = 0.0
    total_loss_denominator = 0.0
    saw_batch = False
    observed = 0
    optimizer.zero_grad(set_to_none=True)
    window_batches: list[tuple[dict[str, torch.Tensor], torch.Tensor, torch.Tensor]] = []
    window_denominator = 0.0
    for batch in loader:
        saw_batch = True
        batch = move_batch(batch, device)
        output = model(batch)
        if model.config.uncertainty_mode == "heteroscedastic":
            _, _, loss_numerator, loss_denominator = masked_multitask_gaussian_loss(
                output["mean"], output["log_variance"], batch["labels"],
                batch["label_mask"], batch["evidence_weight"], return_components=True)
        else:
            _, _, loss_numerator, loss_denominator = masked_multitask_mse_loss(
                output["mean"], batch["labels"], batch["label_mask"],
                batch["evidence_weight"], return_components=True)
        numerator = loss_numerator.detach()
        denominator = loss_denominator.detach()
        window_batches.append((batch, numerator, denominator))
        window_denominator += float(denominator)
        total_loss_numerator += float(numerator)
        total_loss_denominator += float(denominator)
        observed += int(batch["label_mask"].sum())
        if len(window_batches) >= accum_steps:
            _backprop_window(model, device, window_batches, window_denominator, optimizer)
            window_batches = []
            window_denominator = 0.0
    if window_batches:
        _backprop_window(model, device, window_batches, window_denominator, optimizer)
    if not saw_batch:
        raise ValueError("data loader is empty")
    if total_loss_denominator <= 0:
        raise ValueError("data loader contains no observed labels with positive evidence weight")
    return EpochMetrics(total_loss_numerator / total_loss_denominator, observed)


def _backprop_window(model: BioCandidateRanker, device: torch.device,
                     window_batches: list[tuple[dict[str, torch.Tensor], torch.Tensor, torch.Tensor]],
                     window_denominator: float,
                     optimizer: torch.optim.Optimizer) -> None:
    """Backpropagate one accumulated window as ``sum(numerator)/sum(denominator)``.
中文：以 sum(分子)/sum(分母) 回传一个累积窗口。
"""
    if window_denominator <= 0:
        return
    scale = 1.0 / window_denominator
    for batch, numerator, denominator in window_batches:
        if float(denominator) > 0:
            batch = move_batch(batch, device)
            output = model(batch)
            if model.config.uncertainty_mode == "heteroscedastic":
                loss, _, _, _ = masked_multitask_gaussian_loss(
                    output["mean"], output["log_variance"], batch["labels"],
                    batch["label_mask"], batch["evidence_weight"], return_components=True)
            else:
                loss, _, _, _ = masked_multitask_mse_loss(
                    output["mean"], batch["labels"], batch["label_mask"],
                    batch["evidence_weight"], return_components=True)
            (loss * (scale * float(denominator))).backward()
    torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
    optimizer.step()
    optimizer.zero_grad(set_to_none=True)


def run_ranking_epoch(
    model: BioCandidateRanker,
    loader: DataLoader,
    device: torch.device,
    optimizer: torch.optim.Optimizer | None = None,
    task_index: int = 0,
) -> RankingEpochMetrics:
    """One pass of within-campaign pairwise logistic ranking loss on ``task_index``.

    Trains when ``optimizer`` is given, else evaluates. The reported loss is the
    pair-count-weighted mean over the epoch.
    
中文：在 task_index 上以 campaign 内成对 logistic 排序损失跑一遍。
"""
    training = optimizer is not None
    model.train(training)
    weighted_loss = 0.0
    total_pairs = 0
    saw_batch = False
    for batch in loader:
        saw_batch = True
        batch = move_batch(batch, device)
        if batch["labels"].ndim != 2 or not 0 <= task_index < batch["labels"].shape[1]:
            raise ValueError("task_index is outside the collated label columns")
        observed = batch["label_mask"][:, task_index]
        with torch.set_grad_enabled(training):
            output = model(batch)
            scores = output["mean"][:, task_index][observed]
            labels = batch["labels"][:, task_index][observed]
            campaign_ids = batch["campaign_ids"][observed]
            loss, pair_count = pairwise_logistic_ranking_loss(scores, labels, campaign_ids)
            if training and pair_count:
                optimizer.zero_grad(set_to_none=True)
                loss.backward()
                torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
                optimizer.step()
        weighted_loss += float(loss.detach()) * pair_count
        total_pairs += pair_count
    if not saw_batch:
        raise ValueError("data loader is empty")
    if total_pairs == 0:
        raise ValueError("data loader contains no comparable ranking pairs")
    return RankingEpochMetrics(weighted_loss / total_pairs, total_pairs)


def save_checkpoint(path: str | Path, model: BioCandidateRanker,
                    optimizer: torch.optim.Optimizer, *, epoch: int,
                    data_manifest: dict, metrics: dict,
                    data_loader_generator_state: torch.Tensor | None = None) -> None:
    """Atomically save model/optimizer state plus the identities needed for exact resume.

    The payload records the model config, torch/CUDA/DataLoader RNG states, the
    data manifest, and metrics; ``train --resume`` restores all of them and
    refuses to continue when the embedded identities no longer match the run.
    
中文：原子保存模型/优化器状态及精确恢复所需身份。
"""
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "format_version": 1,
        "model_config": model.config.to_dict(),
        "model_state_dict": model.state_dict(),
        "optimizer_state_dict": optimizer.state_dict(),
        "epoch": epoch,
        "data_manifest": data_manifest,
        "metrics": metrics,
        "torch_version": str(torch.__version__),
        "torch_rng_state": torch.get_rng_state(),
        "cuda_rng_state_all": torch.cuda.get_rng_state_all() if torch.cuda.is_available() else None,
        "data_loader_generator_state": data_loader_generator_state,
    }
    with tempfile.NamedTemporaryFile(dir=target.parent, suffix=".tmp", delete=False) as handle:
        temporary = Path(handle.name)
    try:
        torch.save(payload, temporary)
        os.replace(temporary, target)
    finally:
        if temporary.exists():
            temporary.unlink()


def load_checkpoint(path: str | Path, device: torch.device) -> tuple[BioCandidateRanker, dict]:
    """Load a checkpoint and rebuild the model from its embedded config.

    Rejects payloads written by a different checkpoint format version.
    
中文：按内嵌配置重建并加载检查点，格式不符即拒绝。
"""
    payload = torch.load(path, map_location=device, weights_only=True)
    if payload.get("format_version") != 1:
        raise ValueError("unsupported checkpoint format")
    model = BioCandidateRanker(ModelConfig.from_dict(payload["model_config"]))
    model.load_state_dict(payload["model_state_dict"])
    return model.to(device), payload


def write_metrics(path: str | Path, values: dict) -> None:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(values, indent=2, sort_keys=True), encoding="utf-8")
