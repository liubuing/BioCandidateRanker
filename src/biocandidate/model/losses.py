from __future__ import annotations

import torch
import torch.nn.functional as F


def pairwise_logistic_ranking_loss(
    scores: torch.Tensor,
    labels: torch.Tensor,
    campaign_ids: torch.Tensor,
) -> tuple[torch.Tensor, int]:
    """Logistic ranking loss over same-campaign unordered pairs with unequal labels."""
    if scores.ndim != 1 or labels.ndim != 1 or campaign_ids.ndim != 1:
        raise ValueError("scores, labels, and campaign_ids must be one-dimensional")
    if scores.shape != labels.shape or scores.shape != campaign_ids.shape:
        raise ValueError("scores, labels, and campaign_ids must have identical shapes")
    if not scores.is_floating_point() or not labels.is_floating_point():
        raise ValueError("scores and labels must be floating-point tensors")
    if campaign_ids.dtype != torch.int64:
        raise ValueError("campaign_ids must have dtype int64")
    if scores.device != labels.device or scores.device != campaign_ids.device:
        raise ValueError("scores, labels, and campaign_ids must be on the same device")
    if not torch.isfinite(scores).all() or not torch.isfinite(labels).all():
        raise ValueError("scores and labels must contain only finite values")

    upper_triangle = torch.triu(
        torch.ones((scores.numel(), scores.numel()), dtype=torch.bool, device=scores.device),
        diagonal=1,
    )
    comparable = (
        upper_triangle
        & campaign_ids[:, None].eq(campaign_ids[None, :])
        & labels[:, None].ne(labels[None, :])
    )
    pair_count = int(comparable.sum().item())
    if pair_count == 0:
        return scores.sum() * 0.0, 0
    score_differences = scores[:, None] - scores[None, :]
    directions = torch.sign(labels[:, None] - labels[None, :])
    return F.softplus(-directions[comparable] * score_differences[comparable]).mean(), pair_count


def masked_multitask_gaussian_loss(
    mean: torch.Tensor,
    log_variance: torch.Tensor,
    labels: torch.Tensor,
    label_mask: torch.Tensor,
    evidence_weight: torch.Tensor | None = None,
    *,
    return_components: bool = False,
) -> tuple[torch.Tensor, torch.Tensor] | tuple[
    torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor
]:
    """Heteroscedastic regression loss with missing-label and evidence masks."""
    if (
        mean.shape != labels.shape
        or log_variance.shape != labels.shape
        or label_mask.shape != labels.shape
    ):
        raise ValueError("predictions, labels, log variance, and mask must have identical shapes")
    mask = label_mask.to(mean.dtype)
    if evidence_weight is not None:
        mask = mask * evidence_weight.to(mean.dtype).unsqueeze(-1)
    per_entry = 0.5 * (torch.exp(-log_variance) * (labels - mean).square() + log_variance)
    per_entry = per_entry * mask
    numerator = per_entry.sum()
    denominator = mask.sum()
    total = numerator / denominator.clamp_min(1)
    per_task = per_entry.sum(dim=0) / mask.sum(dim=0).clamp_min(1)
    if return_components:
        return total, per_task, numerator, denominator
    return total, per_task


def masked_multitask_mse_loss(
    mean: torch.Tensor,
    labels: torch.Tensor,
    label_mask: torch.Tensor,
    evidence_weight: torch.Tensor | None = None,
    *,
    return_components: bool = False,
) -> tuple[torch.Tensor, torch.Tensor] | tuple[
    torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor
]:
    """Fixed-variance regression loss with missing-label and evidence masks."""
    if mean.shape != labels.shape or label_mask.shape != labels.shape:
        raise ValueError("predictions, labels, and mask must have identical shapes")
    mask = label_mask.to(mean.dtype)
    if evidence_weight is not None:
        mask = mask * evidence_weight.to(mean.dtype).unsqueeze(-1)
    per_entry = (labels - mean).square() * mask
    numerator = per_entry.sum()
    denominator = mask.sum()
    total = numerator / denominator.clamp_min(1)
    per_task = per_entry.sum(dim=0) / mask.sum(dim=0).clamp_min(1)
    if return_components:
        return total, per_task, numerator, denominator
    return total, per_task
