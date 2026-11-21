"""Check that a checkout contains the surface fields needed for cold starts."""

from pathlib import Path
import unittest


class SurfaceData(unittest.TestCase):
    def test_packaged_fields(self):
        surface = Path(__file__).resolve().parents[1] / "src/data/surface"
        fields = ("background_albedo", "glacier_fraction", "land_sea_mask", "orography",
                  "roughness_length", "sea_ice_cover", "soil_water_capacity", "surface_temperature")
        for grid, nlat in (("T21", 32), ("T31", 48), ("T42", 64), ("T85", 128)):
            for field in fields:
                path = surface / grid / f"N{nlat:03d}_surf_{field}.nc"
                with self.subTest(path=path), path.open("rb") as f:
                    magic = f.read(8)
                    self.assertTrue(magic.startswith(b"CDF") or magic == b"\x89HDF\r\n\x1a\n")
