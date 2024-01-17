"""Unit tests for units module."""

import importlib
import pytest
from efel.units import get_unit
from unittest.mock import patch


def test_get_unit():
    """Unit test for the get_unit function."""
    assert get_unit("time_to_last_spike") == "ms"
    assert get_unit("inv_second_ISI") == "Hz"
    assert get_unit("AP1_amp") != "wrong unit"
    assert get_unit("AP1_amp") == "mV"
    assert get_unit("ohmic_input_resistance") == "MΩ"


@patch('efel.units.pkgutil.get_data')
def test_get_data_failure(mock_get_data):
    """Test for handling failure in loading units.json."""
    mock_get_data.return_value = None

    with pytest.raises(ValueError) as excinfo:
        # Dynamically reload the module to simulate the import with mock
        importlib.reload(importlib.import_module('efel.units'))

    assert str(excinfo.value) == "Failed to load units.json"
