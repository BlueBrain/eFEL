"""Module to get units of efeatures."""
import json
import pkgutil


_units_raw = pkgutil.get_data(__name__, "units.json")

if _units_raw is None:
    raise ValueError("Failed to load units.json")

_units = json.loads(_units_raw)


def get_unit(feature_name: str) -> str:
    """Get the unit of a feature."""
    return _units[feature_name]
