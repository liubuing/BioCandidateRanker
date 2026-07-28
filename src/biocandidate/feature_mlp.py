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
