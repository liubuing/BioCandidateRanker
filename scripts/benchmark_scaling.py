from __future__ import annotations

import argparse
import json
import statistics
import time
from pathlib import Path

import torch

from biocandidate import BioCandidateRanker, ModelConfig


def synthetic_batch(batch_size: int, length: int, device: torch.device) -> dict[str, torch.Tensor]:
    atoms_per_graph = 24
    atom_count = batch_size * atoms_per_graph
    graph_batch = torch.arange(batch_size).repeat_interleave(atoms_per_graph)
    source = torch.arange(atom_count - 1)
    valid = (source + 1) % atoms_per_graph != 0
    source = source[valid]
    target = source + 1
    edge_index = torch.stack((torch.cat((source, target)), torch.cat((target, source))))
    return {
        "sequence_tokens": torch.randint(2, 22, (batch_size, length), device=device),
        "sequence_mask": torch.ones(batch_size, length, dtype=torch.bool, device=device),
        "atom_features": torch.randint(0, 2, (atom_count, 6), device=device),
        "edge_index": edge_index.to(device),
        "edge_features": torch.zeros(edge_index.shape[1], 3, dtype=torch.long, device=device),
        "graph_batch": graph_batch.to(device),
        "context_ids": torch.ones(batch_size, 4, dtype=torch.long, device=device),
        "fba_context": torch.zeros(batch_size, 8, device=device),
        "fba_context_mask": torch.zeros(batch_size, dtype=torch.bool, device=device),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default="artifacts/scaling.json")
    parser.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    parser.add_argument("--batch-size", type=int, default=4)
    args = parser.parse_args()
    device = torch.device(args.device)
    config = ModelConfig(
        d_model=64, num_heads=8, protein_layers=2, molecule_layers=2,
        fusion_layers=1, protein_chunk_size=128, dropout=0.0)
    model = BioCandidateRanker(config).to(device).eval()
    with torch.no_grad():
        model(synthetic_batch(args.batch_size, 256, device))
    if device.type == "cuda":
        torch.cuda.synchronize()
    results = []
    for length in (256, 512, 1024, 2048):
        batch = synthetic_batch(args.batch_size, length, device)
        if device.type == "cuda":
            torch.cuda.reset_peak_memory_stats()
            torch.cuda.synchronize()
        timings = []
        with torch.no_grad():
            for _ in range(10):
                if device.type == "cuda":
                    torch.cuda.synchronize()
                start = time.perf_counter()
                output = model(batch)
                if device.type == "cuda":
                    torch.cuda.synchronize()
                timings.append(time.perf_counter() - start)
        elapsed = statistics.median(timings)
        results.append({
            "length": length,
            "batch_size": args.batch_size,
            "seconds": elapsed,
            "samples_per_second": args.batch_size / elapsed,
            "peak_memory_mb": (
                torch.cuda.max_memory_allocated() / 1024 ** 2 if device.type == "cuda" else None),
            "finite": bool(torch.isfinite(output["mean"]).all()),
        })
    payload = {
        "device": str(device),
        "device_name": torch.cuda.get_device_name(device) if device.type == "cuda" else "CPU",
        "config": config.to_dict(),
        "results": results,
    }
    target = Path(args.output)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
