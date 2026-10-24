"""Feature-based MLP baseline for log10_kcat regression.

Trains a small multilayer perceptron on the classical features from
:mod:`biocandidate.baselines` under a governed protocol: seeded training,
standardization fitted on the train partition only, and the reported epoch
selected by validation RMSE. It exists to give the deep ranker a like-for-like
shallow reference on identical splits.

中文：基于特征的 log10_kcat 回归 MLP 基线。在治理协议下，用
:mod:`biocandidate.baselines` 的经典特征训练一个小型多层感知机：训练固定种子，
标准化统计只在训练分区上拟合，报告的 epoch 由验证集 RMSE 选出。它的存在是为了
在完全相同的切分上给深度排序器提供一个可比的浅层参照。
"""

from __future__ import annotations

import os
import tempfile
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, TensorDataset

from .evaluation import regression_metrics


class FeatureMLP(nn.Module):
    """Plain MLP regressor: Linear/SiLU/Dropout blocks ending in one scalar output.

    中文：普通 MLP 回归器：由 Linear/SiLU/Dropout 块组成，末尾输出一个标量。
    """

    def __init__(
        self, input_dimension: int, hidden_dimensions: tuple[int, ...], dropout: float,
    ) -> None:
        super().__init__()
        dimensions = (input_dimension,) + hidden_dimensions
        layers: list[nn.Module] = []
        for input_size, output_size in zip(dimensions, dimensions[1:]):
            layers.extend((nn.Linear(input_size, output_size), nn.SiLU(), nn.Dropout(dropout)))
        layers.append(nn.Linear(dimensions[-1], 1))
        self.network = nn.Sequential(*layers)

    def forward(self, features: torch.Tensor) -> torch.Tensor:
        return self.network(features).squeeze(-1)


def standardize_features(
    train: np.ndarray, validation: np.ndarray, test: np.ndarray,
) -> tuple[dict[str, np.ndarray], np.ndarray, np.ndarray]:
    """Z-score all partitions with train-only mean and standard deviation.

    Zero-variance columns are rescaled by 1.0 instead of dividing by zero, and
    the fitted statistics are returned so inference can reproduce them.

    中文：用仅在训练分区上计算的均值和标准差对所有分区做 z-score 标准化。
    零方差列改用 1.0 缩放以避免除零；同时返回拟合的统计量，供推理阶段复现。
    """
    mean = train.mean(axis=0, dtype=np.float64)
    standard_deviation = train.std(axis=0, dtype=np.float64)
    standard_deviation[standard_deviation == 0] = 1.0
    scaled = {
        "train": ((train - mean) / standard_deviation).astype(np.float32),
        "validation": ((validation - mean) / standard_deviation).astype(np.float32),
        "test": ((test - mean) / standard_deviation).astype(np.float32),
    }
    return scaled, mean, standard_deviation


def train_feature_mlp(
    features: dict[str, np.ndarray],
    labels: dict[str, np.ndarray],
    *,
    seed: int,
    epochs: int,
    batch_size: int,
    learning_rate: float,
    weight_decay: float,
    hidden_dimensions: tuple[int, ...],
    dropout: float,
    device: torch.device,
) -> tuple[FeatureMLP, dict, np.ndarray, np.ndarray]:
    """Train the MLP; return (model, metrics, feature mean, feature stddev).

    The seed pins torch (CPU and CUDA) and the DataLoader shuffling. The epoch
    with the lowest validation RMSE is restored before test metrics are
    computed, so test labels never influence model selection.

    中文：训练 MLP；返回（模型、指标、特征均值、特征标准差）。种子固定 torch
    （CPU 与 CUDA）及 DataLoader 洗牌。在计算测试指标前会恢复验证 RMSE 最低的
    epoch，因此测试标签绝不影响模型选择。
    """
    if set(features) != {"train", "validation", "test"} or set(labels) != set(features):
        raise ValueError("feature MLP requires train, validation, and test partitions")
    if epochs < 1 or batch_size < 1:
        raise ValueError("epochs and batch_size must be positive")
    scaled, mean, standard_deviation = standardize_features(
        features["train"], features["validation"], features["test"])
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
    model = FeatureMLP(
        scaled["train"].shape[1], hidden_dimensions, dropout).to(device)
    optimizer = torch.optim.AdamW(
        model.parameters(), lr=learning_rate, weight_decay=weight_decay)
    generator = torch.Generator().manual_seed(seed)
    train_dataset = TensorDataset(
        torch.from_numpy(scaled["train"]),
        torch.from_numpy(labels["train"].astype(np.float32)),
    )
    train_loader = DataLoader(
        train_dataset, batch_size=batch_size, shuffle=True, generator=generator)
    validation_features = torch.from_numpy(scaled["validation"]).to(device)
    validation_labels = torch.from_numpy(labels["validation"].astype(np.float32)).to(device)
    best_rmse = float("inf")
    best_epoch = 0
    best_state = None
    history = []
    for epoch in range(1, epochs + 1):
        model.train()
        losses = []
        for batch_features, batch_labels in train_loader:
            batch_features = batch_features.to(device)
            batch_labels = batch_labels.to(device)
            predictions = model(batch_features)
            loss = nn.functional.mse_loss(predictions, batch_labels)
            optimizer.zero_grad(set_to_none=True)
            loss.backward()
            optimizer.step()
            losses.append(float(loss.detach()))
        model.eval()
        with torch.no_grad():
            validation_predictions = model(validation_features)
            validation_rmse = float(torch.sqrt(nn.functional.mse_loss(
                validation_predictions, validation_labels)))
        history.append({
            "epoch": epoch,
            "train_mse": sum(losses) / len(losses),
            "validation_rmse": validation_rmse,
        })
        if validation_rmse < best_rmse:
            best_rmse = validation_rmse
            best_epoch = epoch
            best_state = {
                name: value.detach().cpu().clone() for name, value in model.state_dict().items()}
    if best_state is None:
        raise RuntimeError("feature MLP did not produce a checkpoint")
    model.load_state_dict(best_state)
    model.eval()
    with torch.no_grad():
        test_predictions = model(torch.from_numpy(scaled["test"]).to(device)).cpu()
    metrics = regression_metrics(
        test_predictions, torch.from_numpy(labels["test"].astype(np.float32)))
    return model, {
        "selected_epoch": best_epoch,
        "selected_validation_rmse": best_rmse,
        "history": history,
        "test_metrics": metrics,
    }, mean, standard_deviation


def save_feature_checkpoint(payload: dict, path: str | Path) -> None:
    """Write a torch payload atomically via a temp file plus ``os.replace``.

    中文：通过临时文件加 ``os.replace`` 原子地写入 torch 载荷。
    """
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    file_descriptor, temporary_name = tempfile.mkstemp(
        dir=destination.parent, prefix=destination.name, suffix=".tmp")
    os.close(file_descriptor)
    temporary = Path(temporary_name)
    try:
        torch.save(payload, temporary)
        os.replace(temporary, destination)
    finally:
        if temporary.exists():
            temporary.unlink()
