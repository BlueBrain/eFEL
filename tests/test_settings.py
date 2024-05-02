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

from pathlib import Path
import numpy

import pytest

from efel.settings import Settings
from efel.api import set_setting


def test_set_setting():
    """Test that the set_setting method correctly updates a setting."""
    settings = Settings()
    settings.set_setting("Threshold", -30.0)
    assert settings.Threshold == -30.0


def test_set_setting_invalid_type():
    """Test that the set_setting method raises a ValueError
    when given an invalid type."""
    settings = Settings()
    try:
        settings.set_setting("Threshold", "-30.0")
    except ValueError:
        assert True
    else:
        assert False


def test_set_setting_dependencyfile_path_not_found():
    """Test that the set_setting method raises a FileNotFoundError
    when given a nonexistent file."""
    settings = Settings()
    try:
        settings.set_setting("dependencyfile_path", "nonexistent_file.txt")
    except FileNotFoundError:
        assert True
    else:
        assert False


def test_reset_to_default():
    """Test that the reset_to_default method correctly resets a setting to
    its default value."""
    settings = Settings()
    settings.Threshold = -30.0
    settings.reset_to_default()
    assert settings.Threshold == -20.0
