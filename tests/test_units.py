"""Unit tests for units module."""


from efel.units import get_unit


def test_get_unit():
    """Unit test for the get_unit function."""
    assert get_unit("time_to_last_spike") == "ms"
    assert get_unit("inv_second_ISI") == "Hz"
    assert get_unit("AP1_amp") != "wrong unit"
    assert get_unit("AP1_amp") == "mV"
    assert get_unit("ohmic_input_resistance") == "MΩ"
