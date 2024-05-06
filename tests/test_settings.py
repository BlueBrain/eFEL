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


def test_set_setting():
    """Test that the set_setting method correctly updates a setting."""
    settings = Settings()
    settings.set_setting("Threshold", -30.0)
    assert settings.Threshold == -30.0


def test_set_setting_invalid_type():
    """Test that the set_setting method raises a ValueError
    when given an invalid type."""
    settings = Settings()
    with pytest.raises(ValueError):
        settings.set_setting("Threshold", "-30.0")


def test_set_setting_invalid_type2():
    """Test that the set_setting method raises a ValueError
    when given an invalid type."""
    settings = Settings()
    with pytest.raises(ValueError):
        settings.set_setting("voltage_base_mode", 1)


def test_set_setting_invalid_type3():
    """Test that the set_setting method raises a ValueError
    when given an invalid type."""
    settings = Settings()
    with pytest.raises(ValueError):
        settings.set_setting("max_spike_skip", 1.0)


def test_set_setting_convert_type():
    """Test that the set_setting method accepts int for
    float features"""
    settings = Settings()
    settings.set_setting("Threshold", -30)


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
        "impedance_max_freq: 50.0"
    )
    assert str(settings) == expected_output
