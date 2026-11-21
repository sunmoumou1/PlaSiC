"""End-to-end tests: run via 'make -C src test-ozone' (Python + NumPy/netCDF4)."""
import argparse
import json
import subprocess
import tempfile
from pathlib import Path

import numpy as np
from netCDF4 import Dataset


def write_fixture(path, dates, amounts, problem=None):
    """Independent T21/L10 native coordinates, with a nonuniform spatial tracer."""
    s = np.arange(11) / 10
    half = 0.75 * s + 1.75 * s**3 - 1.5 * s**4
    axes = {
        "lev": ((half[1:] + half[:-1]) / 2, "1"),
        "lat": (np.rad2deg(np.arcsin(np.polynomial.legendre.leggauss(32)[0][::-1])),
                "degrees_north"),
        "lon": (np.arange(64) * 360 / 64, "degrees_east"),
    }
    with Dataset(path, "w", format="NETCDF3_CLASSIC") as ds:
        ds.createDimension("time", len(dates))
        date = ds.createVariable("yyyymm", "i4", ("time",))
        date[:] = dates
        if problem == "date_marker":
            date.missing_value = np.int32(dates[0])
        for name, (values, units) in axes.items():
            ds.createDimension(name, len(values))
            axis = ds.createVariable(name, "f4", (name,))
            axis.units = units
            axis[:] = values[::-1] if problem == name else values
            if problem == "axis_marker" and name == "lon":
                axis.missing_value = np.float32(values[0])

        dims = ("time", "lev", "lat", "lon")
        o3 = ds.createVariable("o3", "f8", dims if problem != "dims" else dims[::-1],
                               fill_value=1e-5)
        o3.units = "mol mol-1" if problem == "units" else "kg kg-1"
        if problem == "dims":
            return
        k, j, i = np.ogrid[:10, :32, :64]
        shape = (1 / (1 + k)) * (1 + j / 32) * (1 + i / 64)
        for t, amount in enumerate(amounts):
            o3[t] = amount * shape
        if problem == "fill":
            o3[0, 1, 2, 3] = o3._FillValue
        elif problem == "missing":
            o3.missing_value = 2e-7
            o3[0, 1, 2, 3] = o3.missing_value
        elif problem == "negative":
            o3[0, 1, 2, 3] = -1e-10
        elif problem == "nan":
            o3[0, 1, 2, 3] = np.nan
        elif problem == "packed":
            o3.scale_factor = 1e-6


def run_tests(executable, work):
    src = Path(__file__).resolve().parents[1] / "src"
    # Last step of Dec 1999; crossing New Year exercises both date components.
    clock = ["--start-year", "1999", "--start-yday", "365", "--start-hour",
             "23", "--start-minute", "15"]
    co2 = ["--co2-forcing", next((src / "data/forcing/CO2").glob("*.csv"))]
    fixed = ["--co2-ppm", "280"]
    index = 0

    def run(steps, ozone=None, *, time=None, gas=None, extra=(), failure=None):
        nonlocal index
        index += 1
        output = work / f"state-{index}"
        args = [executable, "--steps", str(steps), "--data-dir", src / "data/surface",
                "--output", output, *(clock if time is None else time),
                *(fixed if gas is None else gas), *extra]
        if ozone is not None:
            args += ["--ozone-forcing", ozone]
        result = subprocess.run(args, capture_output=True, text=True, timeout=60)
        if failure is not None:
            assert result.returncode and failure in result.stderr, result.stderr
            assert not output.exists()
            return
        assert result.returncode == 0, result.stderr
        return output

    a, b = 1e-7, 2e-6
    dates = [199912, 200001]
    varying, constant, invalid = (work / f"{name}.nc" for name in ("ab", "aa", "bad"))
    write_fixture(varying, dates, [a, b])
    write_fixture(constant, dates, [a, a])
    # A December-only run cannot see the January difference, even with CO2.
    dec = run(1, varying, gas=co2)
    assert dec.read_bytes() == run(1, constant, gas=co2).read_bytes()
    assert dec.read_bytes() != run(1, gas=co2).read_bytes()
    # A two-step integration must see the change in the January record.
    joined = run(2, varying, gas=co2)
    assert joined.read_bytes() != run(2, constant, gas=co2).read_bytes()
    # A fresh process resumes the forcing by the restart clock, not record 0.
    split = run(1, varying, time=["--restart", dec], gas=co2)
    assert split.read_bytes() == joined.read_bytes()
    # CLI and JSON entry points give identical states with fixed CO2 too.
    config = work / "config.json"
    config.write_text(json.dumps({"ozone_forcing": str(varying)}))
    assert run(2, varying).read_bytes() == run(2, extra=["--config", config]).read_bytes()
    config.write_text(json.dumps({"ozone_forcing": None}))
    assert run(1).read_bytes() == run(1, extra=["--config", config]).read_bytes()

    # Resetting the restart clock and selecting a 360-day cold-start clock
    # must use the effective model month (January 0), not Gregorian day counts
    # or the January 2000 date stored in dec.
    write_fixture(invalid, [1], [0.0])
    run(1, invalid, time=["--restart", dec, "--reset-nstep", "0"])
    run(1, invalid, time=["--calendar", "360day", "--start-yday", "30"])
    run(1, invalid, time=["--calendar", "360day", "--start-yday", "31"],
        failure="does not cover model date 0-02")
    write_fixture(invalid, dates[:1], [a])
    run(2, invalid, failure="does not cover model date 2000-01")
    write_fixture(invalid, [199912, 200002], [a, b])
    run(1, invalid, failure="consecutive months")
    for problem, message in [("lev", "'lev' coordinate"), ("lat", "'lat' coordinate"),
                             ("lon", "'lon' coordinate"), ("dims", "shape"),
                             ("units", "kg kg-1"), ("packed", "packing"),
                             ("date_marker", "invalid yyyymm"), ("axis_marker", "coordinate"),
                             ("fill", "invalid o3"), ("missing", "invalid o3"),
                             ("negative", "invalid o3"), ("nan", "invalid o3")]:
        write_fixture(invalid, dates, [a, b], problem)
        run(1, invalid, failure=message)
    run(1, work / "absent.nc", failure="absent.nc")
    print(f"Passed {index} T21/L10 ozone model invocations")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    with tempfile.TemporaryDirectory() as temp:
        run_tests(parser.parse_args().executable.resolve(), Path(temp))
