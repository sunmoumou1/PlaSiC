#!/usr/bin/env python3
"""Figures for the ML emulator tutorial.

Reads metrics.json (+ checkpoints and <var>.nc files for the spatial
panel) from a finished run and writes PNGs:

- learning_curves.png : val MSE vs epoch, one line per model
- comparison.png       : final val MSE of models and baselines
- rollout.png          : rollout MSE vs lead month, one line per model
- spatial.png          : truth vs resnet vs persistence for one val sample

Usage:
    python plot_results.py --data data/synthetic --out figures
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import torch

from plasic_dataset import PlasicDataset


def load_run(data_dir: Path):
    metrics = json.loads((data_dir / "metrics.json").read_text())
    return metrics


def fig_learning_curves(metrics: dict, out: Path) -> None:
    plt.figure(figsize=(7, 4.5))
    for name, r in metrics["models"].items():
        epochs = [h["epoch"] for h in r["history"]]
        vals = [h["val"] for h in r["history"]]
        plt.plot(epochs, vals, marker="o", ms=3, label=name)
    plt.yscale("log")
    plt.xlabel("epoch")
    plt.ylabel("val MSE (log scale)")
    plt.title("Validation MSE during training")
    plt.legend()
    plt.tight_layout()
    plt.savefig(out / "learning_curves.png", dpi=120)
    plt.close()


def fig_comparison(metrics: dict, out: Path) -> None:
    rows = [(n, r["history"][-1]["val"]) for n, r in metrics["models"].items()]
    rows += [(k, v) for k, v in metrics["baselines"].items()]
    rows.sort(key=lambda r: r[1])
    names = [r[0] for r in rows]
    vals = [r[1] for r in rows]
    colors = ["tab:blue" if n in metrics["models"] else "tab:gray"
              for n in names]
    plt.figure(figsize=(7, 4.5))
    plt.barh(names, vals, color=colors)
    plt.xlabel("val MSE (lower is better)")
    plt.title("One-step validation MSE")
    for i, v in enumerate(vals):
        plt.text(v, i, f" {v:.4f}", va="center", fontsize=9)
    plt.tight_layout()
    plt.savefig(out / "comparison.png", dpi=120)
    plt.close()


def fig_rollout(metrics: dict, out: Path) -> None:
    plt.figure(figsize=(7, 4.5))
    for name, r in metrics["models"].items():
        ro = r["rollout_val_mse"]
        leads = sorted(ro, key=lambda k: int(k[1:-1]))
        plt.plot([int(k[1:-1]) for k in leads],
                 [ro[k] for k in leads], marker="o", ms=4, label=name)
    plt.xlabel("lead month")
    plt.ylabel("rollout val MSE")
    plt.title("Autoregressive rollout error growth")
    plt.legend()
    plt.tight_layout()
    plt.savefig(out / "rollout.png", dpi=120)
    plt.close()


def fig_spatial(data_dir: Path, metrics: dict, out: Path,
                model_name: str = "resnet") -> bool:
    """Truth vs model vs persistence maps for one validation sample.

    Returns False when the checkpoint or data is unavailable.
    """
    ckpt = data_dir / f"emulator_{model_name}.pt"
    if not ckpt.exists():
        return False
    from train_emulator import MODELS
    saved = torch.load(ckpt, map_location="cpu", weights_only=False)
    variables = metrics["vars"]
    state = saved["state_dict"]
    channels = len(variables)
    hidden = state["net.0.weight"].shape[0]
    # identify the architecture from the checkpoint keys
    if any(".net." in k for k in state):
        arch = "resnet"
        block_idx = {int(k.split(".")[1]) for k in state
                     if k.startswith("net.") and ".net." in k
                     and k.split(".")[1].isdigit()}
        model = MODELS[arch](channels, hidden=hidden,
                             n_blocks=len(block_idx))
    elif "net.4.weight" in state:
        arch = "cnn"
        model = MODELS[arch](channels, hidden=hidden)
    elif "net.2.weight" in state:
        arch = "mlp"
        model = MODELS[arch](channels, hidden=hidden)
    else:
        arch = "linear"
        model = MODELS[arch](channels, hidden=hidden)
    model.load_state_dict(state)
    model.eval()
    ds = PlasicDataset(data_dir, variables, stats=saved["stats"])
    n_val = metrics["n_val"]
    i = len(ds) - n_val // 2 - 1  # a mid-validation sample
    x, y, _ = ds[i]
    with torch.no_grad():
        pred = model(x.unsqueeze(0)).squeeze(0)
    # denormalize for physical units
    stats = saved["stats"]
    vi = 0  # first variable
    mu, sd = stats[variables[vi]]
    truth = y[vi].numpy() * sd + mu
    pm = pred[vi].numpy() * sd + mu
    pers = x[vi].numpy() * sd + mu
    vmin, vmax = truth.min(), truth.max()
    err = np.abs(pm - truth)
    fig, axes = plt.subplots(1, 4, figsize=(13, 3.5))
    panels = [("truth t+1", truth, vmin, vmax),
              (f"{model_name} pred", pm, vmin, vmax),
              ("persistence", pers, vmin, vmax),
              (f"|{model_name} err|", err, 0, err.max())]
    for ax, (title, fld, lo, hi) in zip(axes, panels):
        im = ax.imshow(fld, origin="lower", vmin=lo, vmax=hi,
                       cmap="RdYlBu_r")
        ax.set_title(title)
        ax.set_xticks([]); ax.set_yticks([])
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    fig.suptitle(f"{variables[vi]}: one-step prediction on a validation month")
    plt.tight_layout()
    plt.savefig(out / "spatial.png", dpi=120)
    plt.close()
    return True


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True, help="run dir with metrics.json")
    ap.add_argument("--out", required=True, help="output dir for PNGs")
    ap.add_argument("--prefix", default="",
                    help="filename prefix, e.g. 'era5_'")
    args = ap.parse_args()
    data_dir = Path(args.data)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    metrics = load_run(data_dir)
    fig_learning_curves(metrics, out)
    fig_comparison(metrics, out)
    fig_rollout(metrics, out)
    best = min(metrics["models"],
               key=lambda n: metrics["models"][n]["history"][-1]["val"])
    ok = fig_spatial(data_dir, metrics, out, model_name=best)
    # rename with prefix
    if args.prefix:
        for p in out.glob("*.png"):
            if not p.name.startswith(args.prefix):
                p.rename(out / f"{args.prefix}{p.name}")
    made = sorted(p.name for p in out.glob("*.png"))
    print("wrote:", made, "" if ok else "(spatial panel skipped)")


if __name__ == "__main__":
    main()
