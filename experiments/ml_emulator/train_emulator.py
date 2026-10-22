#!/usr/bin/env python3
"""Train next-step emulators on PlaSiC output and compare architectures.

Baselines (persistence, climatology) plus four learned models: linear and
mlp (per-gridpoint), cnn and resnet (residual CNNs with spatial context).
Chronological train/val split, masked MSE, per-variable metrics, and
multi-step rollout evaluation.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

import torch
import torch.nn as nn
from torch.utils.data import DataLoader, Subset

from plasic_dataset import PlasicDataset


class Residual(nn.Module):
    """Predict the increment x(t+1) - x(t), not the full state."""

    def __init__(self, body: nn.Module) -> None:
        super().__init__()
        self.body = body

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return x + self.body(x)


class LinearEmulator(nn.Module):
    """Per-gridpoint linear map across channels (1x1 convolution)."""

    def __init__(self, channels: int, hidden: int = 32) -> None:
        super().__init__()
        self.net = nn.Conv2d(channels, channels, 1)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return x + self.net(x)


class MLPEmulator(nn.Module):
    """Per-gridpoint MLP: no spatial context, nonlinear in channels."""

    def __init__(self, channels: int, hidden: int = 32) -> None:
        super().__init__()
        self.net = nn.Sequential(
            nn.Conv2d(channels, hidden, 1), nn.ReLU(),
            nn.Conv2d(hidden, channels, 1),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return x + self.net(x)


class CNNEmulator(nn.Module):
    """Shallow residual CNN: 3x3 convolutions see spatial neighbours."""

    def __init__(self, channels: int, hidden: int = 32) -> None:
        super().__init__()
        self.net = nn.Sequential(
            nn.Conv2d(channels, hidden, 3, padding=1), nn.ReLU(),
            nn.Conv2d(hidden, hidden, 3, padding=1), nn.ReLU(),
            nn.Conv2d(hidden, channels, 3, padding=1),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return x + self.net(x)


class ResBlock(nn.Module):
    def __init__(self, hidden: int) -> None:
        super().__init__()
        self.net = nn.Sequential(
            nn.Conv2d(hidden, hidden, 3, padding=1), nn.ReLU(),
            nn.Conv2d(hidden, hidden, 3, padding=1),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return x + self.net(x)


class ResNetEmulator(nn.Module):
    """Deeper residual CNN: stacked 3x3 residual blocks."""

    def __init__(self, channels: int, hidden: int = 32,
                 n_blocks: int = 3) -> None:
        super().__init__()
        self.net = nn.Sequential(
            nn.Conv2d(channels, hidden, 3, padding=1), nn.ReLU(),
            *[ResBlock(hidden) for _ in range(n_blocks)],
            nn.Conv2d(hidden, channels, 3, padding=1),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return x + self.net(x)


MODELS: dict[str, type] = {
    "linear": LinearEmulator,
    "mlp": MLPEmulator,
    "cnn": CNNEmulator,
    "resnet": ResNetEmulator,
}


def masked_mse(pred: torch.Tensor, target: torch.Tensor,
               mask: torch.Tensor) -> torch.Tensor:
    se = (pred - target) ** 2
    return (se * mask).sum() / mask.sum().clamp(min=1)


@torch.no_grad()
def evaluate(predictor, loader: DataLoader,
             device: torch.device | None = None) -> float:
    """Validation MSE of a baseline or a model.

    predictor is 'persistence' / 'climatology' / nn.Module.
    'climatology' is the time mean over the validation set.
    """
    if predictor == "climatology":
        ys = torch.cat([y for _, y, _ in loader], dim=0)
        clim = ys.mean(dim=0, keepdim=True)  # (1, C, H, W)
    tot, n = 0.0, 0
    for x, y, m in loader:
        if device is not None:
            x, y, m = x.to(device), y.to(device), m.to(device)
        if predictor == "persistence":
            pred = x
        elif predictor == "climatology":
            pred = clim.to(x.device).expand_as(y)
        else:
            predictor.eval()
            pred = predictor(x)
        tot += masked_mse(pred, y, m).item() * x.shape[0]
        n += x.shape[0]
    return tot / max(n, 1)


@torch.no_grad()
def per_variable_mse(model: nn.Module, loader: DataLoader,
                     variables: list[str],
                     device: torch.device | None = None) -> dict[str, float]:
    """One-step MSE per variable."""
    model.eval()
    tot = np.zeros(len(variables))
    cnt = np.zeros(len(variables))
    for x, y, m in loader:
        if device is not None:
            x, y, m = x.to(device), y.to(device), m.to(device)
        se = (model(x) - y) ** 2  # (B, C, H, W)
        for c in range(len(variables)):
            mc = m[:, c]
            tot[c] += (se[:, c] * mc).sum().item()
            cnt[c] += mc.sum().item()
    return {v: tot[c] / max(cnt[c], 1) for c, v in enumerate(variables)}


@torch.no_grad()
def rollout_mse(model: nn.Module, dataset: PlasicDataset,
                i_start: int, i_end: int, steps: int = 4,
                device: torch.device | None = None) -> list[float]:
    """Autoregressive rollout MSE, one value per lead step.

    Feeds the model's own output back as the next input and records how
    the error grows with lead time.
    """
    model.eval()
    errs = np.zeros(steps)
    n = 0
    for i in range(i_start, min(i_end, len(dataset) - steps)):
        x, _, _ = dataset[i]
        cur = x.unsqueeze(0)
        if device is not None:
            cur = cur.to(device)
        for s in range(1, steps + 1):
            cur = model(cur)
            _, y_true, m_true = dataset[i + s - 1]
            if device is not None:
                y_true, m_true = y_true.to(device), m_true.to(device)
            errs[s - 1] += masked_mse(cur, y_true.unsqueeze(0),
                                      m_true.unsqueeze(0)).item()
        n += 1
    return (errs / max(n, 1)).tolist()


def train_one(model: nn.Module, train_loader: DataLoader,
              val_loader: DataLoader, device: torch.device,
              epochs: int, lr: float, name: str) -> tuple[nn.Module, list]:
    opt = torch.optim.Adam(model.parameters(), lr=lr)
    history = []
    for epoch in range(1, epochs + 1):
        model.train()
        tot, n = 0.0, 0
        for x, y, m in train_loader:
            x, y, m = x.to(device), y.to(device), m.to(device)
            opt.zero_grad()
            loss = masked_mse(model(x), y, m)
            loss.backward()
            opt.step()
            tot += loss.item() * x.shape[0]
            n += x.shape[0]
        val = evaluate(model, val_loader, device)
        history.append({"epoch": epoch, "train": tot / n, "val": val})
        print(f"  [{name}] epoch {epoch:3d}  train {tot / n:.6f}  val {val:.6f}")
    return model, history


def main() -> None:
    ap = argparse.ArgumentParser(description="Train tiny emulators on PlaSiC output")
    ap.add_argument("--data", required=True, help="dir with <var>.nc files")
    ap.add_argument("--vars", nargs="+", default=["tas", "ps", "pr"])
    ap.add_argument("--levels", type=int, default=3)
    ap.add_argument("--models", default="linear,mlp,cnn,resnet",
                    help="comma-separated subset of "
                         f"{sorted(MODELS)} or 'all'")
    ap.add_argument("--epochs", type=int, default=10)
    ap.add_argument("--hidden", type=int, default=32)
    ap.add_argument("--batch", type=int, default=8)
    ap.add_argument("--lr", type=float, default=1e-3)
    ap.add_argument("--val-frac", type=float, default=0.25)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--device", default="auto",
                    help="'auto', 'cpu', or 'cuda'")
    args = ap.parse_args()

    wanted = sorted(MODELS) if args.models == "all" \
        else [m.strip() for m in args.models.split(",")]
    unknown = [m for m in wanted if m not in MODELS]
    if unknown:
        raise SystemExit(f"unknown model(s): {unknown}; choose from {sorted(MODELS)}")

    torch.manual_seed(args.seed)
    np.random.seed(args.seed)
    device = torch.device(
        "cuda" if (args.device == "cuda" or
                   (args.device == "auto" and torch.cuda.is_available()))
        else "cpu")
    print(f"device: {device}")

    full = PlasicDataset(args.data, args.vars, levels=args.levels)
    n_val = max(1, int(len(full) * args.val_frac))
    n_train = len(full) - n_val
    # Split chronologically: validate on the most recent months, so no
    # future information leaks into training.
    train_loader = DataLoader(Subset(full, range(n_train)),
                              batch_size=args.batch, shuffle=True)
    val_loader = DataLoader(Subset(full, range(n_train, len(full))),
                            batch_size=args.batch)

    print(f"train pairs: {n_train}, val pairs: {n_val}, vars: {args.vars}")
    baselines = {
        "persistence": evaluate("persistence", val_loader),
        "climatology": evaluate("climatology", val_loader),
    }
    for k, v in baselines.items():
        print(f"baseline[{k}] val MSE: {v:.6f}")

    out_dir = Path(args.data)
    results: dict = {"vars": args.vars, "n_train": n_train, "n_val": n_val,
                     "baselines": baselines, "models": {}}
    summary = []
    for name in wanted:
        print(f"training [{name}] ...")
        model = MODELS[name](len(args.vars), hidden=args.hidden).to(device)
        n_params = sum(p.numel() for p in model.parameters())
        model, history = train_one(model, train_loader, val_loader,
                                   device, args.epochs, args.lr, name)
        per_var = per_variable_mse(model, val_loader, args.vars, device)
        ro = rollout_mse(model, full, n_train, len(full), steps=4,
                         device=device)
        ckpt = out_dir / f"emulator_{name}.pt"
        torch.save({"state_dict": model.state_dict(), "vars": args.vars,
                    "stats": full.stats}, ckpt)
        val_mse = history[-1]["val"]
        summary.append((name, n_params, val_mse))
        results["models"][name] = {
            "n_params": n_params,
            "history": history,
            "per_variable_val_mse": per_var,
            "rollout_val_mse": {f"+{s}m": e for s, e in enumerate(ro, 1)},
        }
        print(f"  [{name}] params: {n_params}, final val MSE: {val_mse:.6f}, "
              f"saved -> {ckpt}")

    print("\n=== val MSE comparison (lower is better) ===")
    rows = [("persistence", "-", baselines["persistence"]),
            ("climatology", "-", baselines["climatology"])]
    rows += [(n, str(p), v) for n, p, v in summary]
    for n, p, v in sorted(rows, key=lambda r: r[2]):
        print(f"  {n:12s} params={p:>8s}  val MSE={v:.6f}")

    mpath = out_dir / "metrics.json"
    mpath.write_text(json.dumps(results, indent=2))
    print(f"metrics -> {mpath}")


if __name__ == "__main__":
    main()
