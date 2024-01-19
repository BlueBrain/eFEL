"""Module containing access functions to C++ features for Python features."""
from __future__ import annotations
import warnings
import numpy as np
from efel import cppcore


def get_cpp_feature(feature_name: str, raise_warnings=False) -> np.ndarray | None:
    """Return value of feature implemented in cpp."""
    cppcoreFeatureValues: list[int | float] = list()
    exitCode = cppcore.getFeature(feature_name, cppcoreFeatureValues)
    if exitCode < 0:
        if raise_warnings:
            warnings.warn(
                f"Error while calculating {feature_name}, {cppcore.getgError()}",
                RuntimeWarning)
        return None
    return np.array(cppcoreFeatureValues)


def _get_cpp_data(data_name: str) -> float | int:
    """Get cpp data value."""
    try:
        return cppcore.getMapDoubleData(data_name)[0]
    except Exception:
        return cppcore.getMapIntData(data_name)[0]
