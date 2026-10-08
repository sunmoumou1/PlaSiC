#!/usr/bin/env python3
"""PyTorch Dataset for PlaSiC monthly-mean NetCDF output.

Each variable lives in its own `<name>.nc` file. The dataset builds
(X_t, X_{t+1}) next-step pairs with per-variable standardization;
masked (_FillValue) points are excluded from the loss via the mask.
"""
from __future__ import annotations

from pathlib import Path

import netCDF4
import numpy as np
import torch
from torch.utils.data import Dataset


def load_variable(data_dir: str | Path, name: str,
                  levels: int | None = None) -> tuple[np.ndarray, np.ndarray]:
    """Load `<name>.nc` -> (data, valid_mask); data is (time, lat, lon)."""
    path = Path(data_dir) / f"{name}.nc"
    if not path.exists():
        raise FileNotFoundError(f"PlaSiC output not found: {path}")
    with netCDF4.Dataset(path, "r") as ds:
        var = ds.variables[name]
        arr = var[:]  # netCDF4 auto-masks _FillValue into a MaskedArray
        if isinstance(arr, np.ma.MaskedArray):
            valid = ~np.ma.getmaskarray(arr)
            data = arr.filled(np.nan).astype(np.float64)
        else:
            data = np.asarray(arr, dtype=np.float64)
            fill = getattr(var, "_FillValue", None)
            if fill is not None and not np.isnan(fill):
                valid = data != fill
                data = np.where(valid, data, np.nan)
            else:
                valid = np.ones(data.shape, dtype=bool)
        if data.ndim == 4:  # (time, level, lat, lon): average low levels
            nlev = data.shape[1] if levels is None else min(levels, data.shape[1])
            with np.errstate(invalid="ignore"):
                data = np.nanmean(data[:, :nlev], axis=1)
            valid = valid[:, :nlev].any(axis=1)
        elif data.ndim != 3:
            raise ValueError(f"unexpected dims for {name}: {var.dimensions}")
    return data.astype(np.float32), valid


def compute_stats(data: np.ndarray) -> tuple[float, float]:
    return float(np.nanmean(data)), float(np.nanstd(data)) + 1e-8


class PlasicDataset(Dataset):
    """(X_t, X_{t+1}) next-step pairs from PlaSiC monthly output.

    Each sample: x (C, H, W) standardized state at month t,
    y (C, H, W) standardized state at month t+1,
    m (C, H, W) bool mask (valid at both t and t+1).
    """

    def __init__(self, data_dir: str | Path, variables: list[str],
                 levels: int | None = 3,
                 stats: dict[str, tuple[float, float]] | None = None,
                 t_start: int = 0, t_end: int | None = None) -> None:
        self.variables = variables
        fields, masks = [], []
        for v in variables:
            data, valid = load_variable(data_dir, v, levels=levels)
            fields.append(data)
            masks.append(valid)
        self.n_time = fields[0].shape[0]
        assert all(f.shape[0] == self.n_time for f in fields), \
            "time dimension mismatch between variables"
        t_end = self.n_time - 1 if t_end is None else min(t_end, self.n_time - 1)
        self.t0, self.t1 = t_start, t_end
        self.stats = stats or {v: compute_stats(f)
                               for v, f in zip(variables, fields)}
        stacked = np.stack([(f - self.stats[v][0]) / self.stats[v][1]
                            for v, f in zip(variables, fields)])
        # Zero masked points so NaNs can't enter the loss; validity is
        # tracked by the mask, not the value.
        self.fields = np.nan_to_num(stacked, nan=0.0).astype(np.float32)
        self.masks = np.stack(masks)

    def __len__(self) -> int:
        return max(0, self.t1 - self.t0)

    def __getitem__(self, i: int):
        t = self.t0 + i
        x = torch.from_numpy(self.fields[:, t])
        y = torch.from_numpy(self.fields[:, t + 1])
        m = torch.from_numpy(self.masks[:, t + 1] & self.masks[:, t])
        return x, y, m
