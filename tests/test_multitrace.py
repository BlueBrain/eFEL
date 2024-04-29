"""Unit tests for the multitrace module."""

from pathlib import Path
import efel
from efel.io import load_ascii_input
from efel.pyfeatures.multitrace import bpap_attenuation


testdata_dir = Path(__file__).parent / "testdata"
meanfrequency1_url = testdata_dir / "basic" / "mean_frequency_1.txt"


def test_bpap_attenuation():
    efel.reset()
    stim_start, stim_end = 1.0, 2.0  # dummy values not used in this feature

    time, voltage = load_ascii_input(meanfrequency1_url)
    soma_trace = dendrite_trace = {
        "T": time,
        "V": voltage,
        "stim_start": [stim_start],
        "stim_end": [stim_end],
    }
    assert bpap_attenuation(soma_trace, dendrite_trace) == 1.0

    # test voltage base subtraction
    soma_trace = {
        "T": time,
        "V": voltage,
        "stim_start": [stim_start],
        "stim_end": [stim_end],
    }
    # subtract 10 mv from V of soma_trace
    soma_trace["V"] = soma_trace["V"] - 10
    assert bpap_attenuation(soma_trace, dendrite_trace) == 1.0

    # divide by 2
    soma_trace["V"] = voltage / 2
    assert bpap_attenuation(soma_trace, dendrite_trace) == 0.5
