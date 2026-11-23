"""Shared paths and input fingerprints for the historical ozone experiments."""
import hashlib
import os
from pathlib import Path


CASES = ('co2_ozone', 'ozone_only', 'co2_only')
RUNS = Path(os.environ.get('PLASIC_OZONE_RUNS', Path(__file__).resolve().parent / 'runs'))


def sha256(path):
    checksum = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            checksum.update(chunk)
    return checksum.hexdigest()
