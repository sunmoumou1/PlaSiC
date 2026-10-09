#!/usr/bin/env python3
"""Generate tiny synthetic NetCDF files mimicking PlaSiC monthly output.

The signal is an eastward-propagating wave + seasonal cycle + noise,
so the task is learnable but not trivial.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

import netCDF4

FILL = -1.0e20


def write_var(path: Path, name: str, standard_name: str, long_name: str,
              units: str, data: np.ndarray, mask: np.ndarray) -> None:
    nt, nlat, nlon = data.shape
    lat = np.linspace(-90, 90, nlat, dtype=np.float32)
    lon = np.linspace(0, 360, nlon, dtype=np.float32, endpoint=False)
    with netCDF4.Dataset(path, "w") as ds:
        ds.createDimension("time", nt)
        ds.createDimension("lat", nlat)
        ds.createDimension("lon", nlon)
        ds.createDimension("bnds", 2)
        tv = ds.createVariable("time", "f8", ("time",))
        tv.standard_name = "time"
        tv.units = "days since 1850-01-01 00:00:00"
        tv.calendar = "noleap"
        tv.bounds = "time_bnds"
        tv[:] = 15.0 + 30.0 * np.arange(nt)  # mid-month
        bv = ds.createVariable("time_bnds", "f8", ("time", "bnds"))
        bv[:] = np.stack([30.0 * np.arange(nt), 30.0 * (np.arange(nt) + 1)], axis=1)
        lv = ds.createVariable("lat", "f4", ("lat",))
        lv.standard_name = "latitude"
        lv.units = "degrees_north"
        lv[:] = lat
        ov = ds.createVariable("lon", "f4", ("lon",))
        ov.standard_name = "longitude"
        ov.units = "degrees_east"
        ov[:] = lon
        vv = ds.createVariable(name, "f4", ("time", "lat", "lon"),
                               fill_value=FILL, zlib=True, complevel=2)
        vv.standard_name = standard_name
        vv.long_name = long_name
        vv.units = units
        vv.cell_methods = "time: mean"
        vv[:] = np.where(mask[None], data, FILL).astype(np.float32)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--months", type=int, default=24)
    ap.add_argument("--nlat", type=int, default=32)
    ap.add_argument("--nlon", type=int, default=64)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    nt, nlat, nlon = args.months, args.nlat, args.nlon
    lat = np.linspace(-90, 90, nlat)[:, None]
    lon = np.linspace(0, 360, nlon, endpoint=False)[None, :]
    t = np.arange(nt)[:, None, None]

    wave = np.sin(np.deg2rad(3 * lon - 25 * t)) * np.cos(np.deg2rad(2 * lat))
    seasonal = 12.0 * np.sin(2 * np.pi * t / 12.0) * np.cos(np.deg2rad(lat))
    noise = rng.normal(0, 0.8, (nt, nlat, nlon))

    tas = 273.0 + 28.0 * np.cos(np.deg2rad(lat)) + 0.4 * seasonal + 3.0 * wave + noise
    ps = 101325.0 - 8000.0 * (1 - np.cos(np.deg2rad(lat))) + 300.0 * wave + \
        rng.normal(0, 50, (nt, nlat, nlon))
    pr = np.maximum(0.0, 2.0 + 4.0 * wave + 0.05 * seasonal +
                    rng.normal(0, 1.5, (nt, nlat, nlon))) * 1e-5

    valid = np.ones((nlat, nlon), bool)
    valid[:, ::4] = False  # masked columns to cover the _FillValue path

    write_var(out / "tas.nc", "tas", "air_temperature",
              "Near-surface air temperature", "K", tas, valid)
    write_var(out / "ps.nc", "ps", "surface_air_pressure",
              "Surface pressure", "Pa", ps, valid)
    write_var(out / "pr.nc", "pr", "precipitation_flux",
              "Total precipitation flux", "kg m-2 s-1", pr, valid)
    print(f"wrote tas.nc ps.nc pr.nc -> {out} ({nt} months, {nlat}x{nlon})")


if __name__ == "__main__":
    main()
