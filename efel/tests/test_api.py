"""Unit tests for api.py."""
from pathlib import Path

import numpy as np
import pandas as pd

import efel

testdata_dir = Path(__file__).resolve().parent / 'testdata'


def test_get_features_df():
    """Test returning the features as a dataframe."""
    trace_file = testdata_dir / 'basic' / 'min_AHP_values_single_peak.txt'
    trace_values = np.loadtxt(trace_file)

    trace = {
        "T": trace_values[:, 0],
        "V": trace_values[:, 1],
        "stim_start": [1950],
        "stim_end": [2050],
    }
    efeature_names = ["min_AHP_values", "min_AHP_indices", "peak_indices"]
    features_df = efel.get_features_df([trace], efeature_names)
    assert isinstance(features_df, pd.DataFrame)
    assert features_df["peak_indices"][0] == [237]
    assert features_df["min_AHP_indices"][0] is None
    assert features_df["min_AHP_values"][0] is None
