"""Unit tests for the validation module."""
from pathlib import Path

import numpy as np
import pytest

import efel
from efel.io import load_ascii_input
from efel.pyfeatures.validation import check_ais_initiation


testdata_dir = Path(__file__).parent / "testdata"
meanfrequency1_url = testdata_dir / "basic" / "mean_frequency_1.txt"


def test_check_ais_initiation():
    """Unit test for check_ais_initiation."""
    efel.reset()
    stim_start, stim_end = 1.0, 2.0  # dummy values not used in this feature

    time, voltage = load_ascii_input(meanfrequency1_url)
    soma_trace = ais_trace = {
        "T": time,
        "V": voltage,
        "stim_start": [stim_start],
        "stim_end": [stim_end],
    }
    assert check_ais_initiation(soma_trace, ais_trace)

    # test edge cases

    # 0. no spikes
    ais_trace = {
        "T": time[35000:],  # slice to remove spikes
        "V": voltage[35000:],
        "stim_start": [stim_start],
        "stim_end": [stim_end],
    }

    with pytest.raises(ValueError) as excinfo:
        check_ais_initiation(soma_trace, ais_trace)
    assert "AP_begin_time feature is not computed" in str(excinfo.value)

    # 1. different number of spikes
    ais_trace = {
        "T": time[15000:],  # slice to remove spikes
        "V": voltage[15000:],
        "stim_start": [stim_start],
        "stim_end": [stim_end],
    }
    with pytest.raises(ValueError) as excinfo:
        check_ais_initiation(soma_trace, ais_trace)
    assert "Number of APs in soma and AIS do not match" in str(excinfo.value)

    # 2. spike in soma starts earlier than in AIS
    ais_trace = {
        "T": time,
        "V": np.roll(voltage, 500),  # shift the spike by 500 samples,
        "stim_start": [stim_start],
        "stim_end": [stim_end],
    }

    with pytest.raises(ValueError) as excinfo:
        check_ais_initiation(soma_trace, ais_trace)
    assert "There is a spike that initiates in the soma before the axon" in str(
        excinfo.value
    )
