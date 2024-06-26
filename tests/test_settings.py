# pylint: disable=W0611, W0612, F0401, R0914, C0302

"""Settings tests of eFEL"""


"""
Copyright (c) 2015, Blue Brain Project/EPFL

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
  * Redistributions of source code must retain the above copyright
    notice, this list of conditions and the following disclaimer.
  * Redistributions in binary form must reproduce the above copyright
    notice, this list of conditions and the following disclaimer in the
    documentation and/or other materials provided with the distribution.
  * Neither the name of the copyright holder nor the
    names of its contributors may be used to endorse or promote products
    derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

import pytest

from efel.settings import Settings
from efel.api import set_setting, get_settings
import logging


def test_set_setting():
    """Test that the set_setting method correctly updates a setting."""
    settings = Settings()
    settings.set_setting("Threshold", -30.0)
    assert settings.Threshold == -30.0


@pytest.mark.parametrize("setting_name, new_value, converted_value, expected_type", [
    ("Threshold", "-30.0", -30.0, float),
    ("strict_stiminterval", 0, False, bool),
    ("initburst_freq_threshold", -50.9, -50, int),
    ("initburst_sahp_start", 5.5, 5, int)
])
def test_set_setting_conversion(caplog,
                                setting_name,
                                new_value,
                                converted_value,
                                expected_type):
    """Test that the set_setting method correctly updates a setting
    and logs a debug warning when converting types."""
    settings = Settings()

    if setting_name == "initburst_freq_threshold":
        logger_level = logging.WARNING
    else:
        logger_level = logging.DEBUG

    with caplog.at_level(logger_level):
        settings.set_setting(setting_name, new_value)

        expected_log_message = (
            "Value '%s' of type '%s' for setting '%s' "
            "has been converted to '%s' of type '%s'." % (
                new_value,
                type(new_value).__name__,
                setting_name,
                converted_value,
                expected_type.__name__
            )
        )
    assert any(record.message == expected_log_message for record in caplog.records)


def test_set_setting_new_setting(caplog):
    """Test that the set_setting method correctly adds a new setting
    when the setting is not present."""
    settings = Settings()
    setting_name = "stim_start"
    new_value = 100

    with caplog.at_level(logging.DEBUG):
        settings.set_setting(setting_name, new_value)

    assert getattr(settings, setting_name) == new_value
    expected_log_message = (
        "Setting '%s' not found in settings. "
        "Adding it as a new setting." % setting_name
    )
    assert any(record.message == expected_log_message for record in caplog.records)
    assert getattr(settings, setting_name) == new_value


def test_set_setting_invalid_type():
    """Test that the set_setting raises a ValueError when given an invalid type."""
    settings = Settings()
    with pytest.raises(ValueError):
        settings.set_setting("Threshold", [-30.0])


def test_set_setting_dependencyfile_path_not_found():
    """Test that the set_setting method raises a FileNotFoundError
    when given a nonexistent file."""
    settings = Settings()
    with pytest.raises(FileNotFoundError):
        settings.set_setting("dependencyfile_path", "nonexistent_file.txt")


def test_reset_to_default():
    """Test that the reset_to_default method correctly resets a setting to
    its default value."""
    settings = Settings()
    settings.Threshold = -30.0
    settings.reset_to_default()
    assert settings.Threshold == -20.0


def test_get_settings():
    """Test that the get_settings method returns an instance of efel.Settings."""
    settings = get_settings()
    assert isinstance(settings, Settings)


def test_str_method():
    """Test that the __str__ method returns the correct string representation."""
    settings = Settings()
    expected_output = (
        "Threshold: -20.0\n"
        "DerivativeThreshold: 10.0\n"
        "DownDerivativeThreshold: -12.0\n"
        f"dependencyfile_path: {settings.dependencyfile_path}\n"
        "spike_skipf: 0.1\n"
        "max_spike_skip: 2\n"
        "interp_step: 0.1\n"
        "burst_factor: 1.5\n"
        "strict_burst_factor: 2.0\n"
        "voltage_base_start_perc: 0.9\n"
        "voltage_base_end_perc: 1.0\n"
        "current_base_start_perc: 0.9\n"
        "current_base_end_perc: 1.0\n"
        "rise_start_perc: 0.0\n"
        "rise_end_perc: 1.0\n"
        "initial_perc: 0.1\n"
        "min_spike_height: 20.0\n"
        "strict_stiminterval: False\n"
        "initburst_freq_threshold: 50\n"
        "initburst_sahp_start: 5\n"
        "initburst_sahp_end: 100\n"
        "DerivativeWindow: 3\n"
        "voltage_base_mode: mean\n"
        "current_base_mode: mean\n"
        "precision_threshold: 1e-10\n"
        "sahp_start: 5.0\n"
        "ignore_first_ISI: True\n"
        "impedance_max_freq: 50.0\n"
        "AP_phaseslope_range: 2"
    )
    assert str(settings) == expected_output
