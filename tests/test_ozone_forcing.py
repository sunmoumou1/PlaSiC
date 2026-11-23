"""Expose the selectable serial/MPI ozone driver to standard unittest discovery."""

if __package__:
    from .ozone_forcing import MonthlyOzone
else:
    from ozone_forcing import MonthlyOzone
