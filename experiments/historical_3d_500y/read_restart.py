#!/usr/bin/env python3
"""Minimal reader for the PlaSiC sequential-unformatted restart format.

Each item is stored as

    uint32 name_bytes (== 16)
    16-byte blank-padded name
    uint32 name_bytes (== 16)
    uint32 payload_bytes
    payload_bytes of data
    uint32 payload_bytes

Array payloads are column-major: for each column the record stores
``stored_rows`` contiguous float32 values (``leading_dimension >= stored_rows``
is collapsed by the writer, so the stored data is simply
``columns * stored_rows`` floats).

This helper is used to inspect deep-ocean records (``deept``/``deeps``) saved by
the 500-year spin-up and to compute ocean-mean vertical profiles without having
to write monthly NetCDF diagnostics for the whole spin-up.
"""

from __future__ import annotations

import argparse
import struct
from pathlib import Path

import numpy as np


NAME_BYTES = 16


def read_records(path: Path) -> dict[str, np.ndarray]:
    """Return ``{name: flat float32 array}`` for every named record."""
    records: dict[str, np.ndarray] = {}
    with path.open("rb") as stream:
        while True:
            lead = stream.read(4)
            if lead == b"":
                break
            if len(lead) != 4:
                raise ValueError("truncated name marker")
            (name_bytes,) = struct.unpack("<I", lead)
            if name_bytes != NAME_BYTES:
                raise ValueError(f"unexpected name length {name_bytes}")
            name = stream.read(NAME_BYTES)
            trail = stream.read(4)
            if len(name) != NAME_BYTES or len(trail) != 4:
                raise ValueError("truncated name record")
            (trail_bytes,) = struct.unpack("<I", trail)
            if trail_bytes != name_bytes:
                raise ValueError("mismatched name markers")

            lead = stream.read(4)
            if len(lead) != 4:
                raise ValueError("truncated data marker")
            (payload_bytes,) = struct.unpack("<I", lead)
            payload = stream.read(payload_bytes)
            trail = stream.read(4)
            if len(payload) != payload_bytes or len(trail) != 4:
                raise ValueError("truncated data record")
            (trail_bytes,) = struct.unpack("<I", trail)
            if trail_bytes != payload_bytes:
                raise ValueError("mismatched data markers")

            key = name.decode("ascii", "replace").strip()
            records[key] = np.frombuffer(payload, dtype="<f4").copy()
    return records


def gaussian_weights(latitudes_deg: np.ndarray) -> np.ndarray:
    nodes, weights = np.polynomial.legendre.leggauss(latitudes_deg.size)
    sin_lat = np.sin(np.deg2rad(latitudes_deg))
    closest = np.argmin(np.abs(nodes[None, :] - sin_lat[:, None]), axis=1)
    return weights[closest].astype(np.float64)


def deep_profile(path: Path, nlat: int = 32, nlon: int = 64,
                 nlev: int = 16) -> dict[str, np.ndarray]:
    """Return global ocean-mean temperature/salinity profiles and weights."""
    records = read_records(path)
    if "deept" not in records:
        raise KeyError(f"{path} has no 'deept' record (not a 3-D-ocean restart)")
    land = records["yls"].reshape(nlat, nlon)
    land_mask = land >= 0.5
    nhor = nlat * nlon

    def profile(name: str) -> np.ndarray:
        field = records[name]
        if field.size != nhor * nlev:
            raise ValueError(
                f"{name} has {field.size} values, expected {nhor * nlev}")
        field = field.reshape(nlev, nlat, nlon)
        out = np.empty(nlev, dtype=np.float64)
        for level in range(nlev):
            layer = field[level]
            ocean = ~land_mask
            out[level] = layer[ocean].mean() if ocean.any() else np.nan
        return out

    # Area-weighted variant, using the Gaussian latitudes of the T21 grid.
    lat = np.array([
        85.76059, 80.268776, 74.74454, 69.212975, 63.678635, 58.142956,
        52.593815, 47.035643, 41.471423, 35.902974, 30.331470, 24.757678,
        19.182118, 13.605165, 8.027107, 2.448172, -2.448172, -8.027107,
        -13.605165, -19.182118, -24.757678, -30.331470, -35.902974,
        -41.471423, -47.035643, -52.593815, -58.142956, -63.678635,
        -69.212975, -74.74454, -80.268776, -85.76059,
    ])
    weights = gaussian_weights(lat)
    ocean_weights = np.where(land_mask, 0.0, weights[:, None])

    def weighted_profile(name: str) -> np.ndarray:
        field = records[name].reshape(nlev, nlat, nlon)
        out = np.empty(nlev, dtype=np.float64)
        for level in range(nlev):
            layer = np.where(land_mask, 0.0, field[level])
            out[level] = float((layer * ocean_weights).sum() /
                                ocean_weights.sum())
        return out

    return {
        "temperature": weighted_profile("deept"),
        "salinity": weighted_profile("deeps"),
        "temperature_unweighted": profile("deept"),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("restart", type=Path)
    parser.add_argument("--list", action="store_true",
                        help="list record names and sizes")
    parser.add_argument("--profile", action="store_true",
                        help="print the deep-ocean temperature profile")
    args = parser.parse_args()

    if args.list:
        for name, values in read_records(args.restart).items():
            print(f"{name:12s} {values.size}")
    if args.profile:
        result = deep_profile(args.restart)
        depths = np.array([
            25, 75, 145, 255, 405, 595, 825, 1095, 1405, 1755, 2145, 2575,
            3045, 3555, 4105, 4695], dtype=float)
        print(f"{'depth_m':>8s} {'T_K':>10s} {'S':>8s}")
        for depth, temp, salt in zip(
                depths, result["temperature"], result["salinity"]):
            print(f"{depth:8.0f} {temp:10.4f} {salt:8.3f}")
        heat = 0.0
        thickness = np.array([50.0] + [50.0 + 40.0 * (k - 1)
                                       for k in range(1, 16)])
        heat = float((result["temperature"] * thickness).sum())
        print(f"column-integrated temperature (K*m): {heat:.1f}")


if __name__ == "__main__":
    main()
