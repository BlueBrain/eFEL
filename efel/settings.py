"""efel Settings class"""

"""
Copyright (c) 2015, EPFL/Blue Brain Project

 This file is part of eFEL <https://github.com/BlueBrain/eFEL>

 This library is free software; you can redistribute it and/or modify it under
 the terms of the GNU Lesser General Public License version 3.0 as published
 by the Free Software Foundation.

 This library is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
 details.

 You should have received a copy of the GNU Lesser General Public License
 along with this library; if not, write to the Free Software Foundation, Inc.,
 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
"""

from dataclasses import dataclass, fields, asdict
from pathlib import Path
from typing import Union
from typing_extensions import deprecated
import logging

logger = logging.getLogger(__name__)


@dataclass
class Settings:
    """eFEL settings class.

    Attributes:
        Threshold (float): Spike detection threshold (default: -20.0).
        DerivativeThreshold (float): Threshold value for derivative calculations
        (default: 10.0).
        DownDerivativeThreshold (float): Threshold value for downward derivative
        calculations (default: -12.0).
        dependencyfile_path (str): Path to the dependency file
        (default: 'DependencyV5.txt').
        spike_skipf (float): Fraction of spikes to skip (default: 0.1).
        max_spike_skip (int): Maximum number of spikes to skip (default: 2).
        interp_step (float): Interpolation step (default: 0.1).
        burst_factor (float): Burst factor (default: 1.5).
        strict_burst_factor (float): Strict burst factor (default: 2.0).
        voltage_base_start_perc (float): Voltage base start percentage (default: 0.9).
        voltage_base_end_perc (float): Voltage base end percentage (default: 1.0).
        current_base_start_perc (float): Current base start percentage (default: 0.9).
        current_base_end_perc (float): Current base end percentage (default: 1.0).
        rise_start_perc (float): Rise start percentage (default: 0.0).
        rise_end_perc (float): Rise end percentage (default: 1.0).
        initial_perc (float): Initial percentage (default: 0.1).
        min_spike_height (float): Minimum spike height (default: 20.0).
        strict_stiminterval (bool): Strict stimulus interval (default: False).
        initburst_freq_threshold (int): Initial burst frequency threshold
        (default: 50)
        initburst_sahp_start (int): Initial burst SAHP start (default: 5).
        initburst_sahp_end (int): Initial burst SAHP end (default: 100).
        DerivativeWindow (int): Derivative window (default: 3).
        voltage_base_mode (str): Voltage base mode (default: "mean").
        current_base_mode (str): Current base mode (default: "mean").
        precision_threshold (float): Precision threshold (default: 1e-10).
        sahp_start (float): SAHP start (default: 5.0).
        ignore_first_ISI (bool): Ignore first ISI (default: True).
        impedance_max_freq (float): Impedance maximum frequency (default: 50.0).
    """

    Threshold: float = -20.0
    DerivativeThreshold: float = 10.0
    DownDerivativeThreshold: float = -12.0
    dependencyfile_path: str = str(
        Path(__file__).parent.absolute() / "DependencyV5.txt"
    )

    spike_skipf: float = 0.1
    max_spike_skip: int = 2
    interp_step: float = 0.1
    burst_factor: float = 1.5
    strict_burst_factor: float = 2.0
    voltage_base_start_perc: float = 0.9
    voltage_base_end_perc: float = 1.0
    current_base_start_perc: float = 0.9
    current_base_end_perc: float = 1.0
    rise_start_perc: float = 0.0
    rise_end_perc: float = 1.0
    initial_perc: float = 0.1
    min_spike_height: float = 20.0
    strict_stiminterval: bool = False
    initburst_freq_threshold: int = 50
    initburst_sahp_start: int = 5
    initburst_sahp_end: int = 100
    DerivativeWindow: int = 3
    voltage_base_mode: str = "mean"
    current_base_mode: str = "mean"
    precision_threshold: float = 1e-10
    sahp_start: float = 5.0
    ignore_first_ISI: bool = True
    impedance_max_freq: float = 50.0
    AP_phaseslope_range: int = 2

    def set_setting(self,
                    setting_name: str,
                    new_value: Union[int, float, str, bool]) -> None:
        """Set a certain setting to a new value.

        Args:
            setting_name (str): Name of the setting to be modified.
            new_value (Union[int, float, str, bool]): New value for the setting.

        Raises:
            ValueError: If the value is of the wrong type.
            FileNotFoundError: If the path to the dependency file does not exist
            (for 'dependencyfile_path' setting).
        """
        if hasattr(self, setting_name):
            expected_types = {f.name: f.type for f in fields(self)}
            expected_type = expected_types.get(setting_name)
            if expected_type is None:
                raise TypeError(f"type not found for setting '{setting_name}'")
            try:
                converted_value = expected_type(new_value)
                if not isinstance(new_value, expected_type):
                    log_message = (
                        "Value '%s' of type '%s' for setting '%s' "
                        "has been converted to '%s' of type '%s'."
                    ) % (
                        new_value,
                        type(new_value).__name__,
                        setting_name,
                        converted_value,
                        expected_type.__name__
                    )
                    if expected_type is int and isinstance(new_value, float) and \
                       new_value != converted_value:
                        logger.warning(log_message)
                    else:
                        logger.debug(log_message)
            except (ValueError, TypeError):
                raise ValueError(f"Invalid value for setting '{setting_name}'. "
                                 f"Expected type: {expected_type.__name__}.")
        else:
            logger.debug("Setting '%s' not found in settings. "
                         "Adding it as a new setting.", setting_name)
            converted_value = new_value

        if setting_name == "dependencyfile_path":
            path = Path(str(converted_value))
            if not path.exists():
                raise FileNotFoundError(f"Path to dependency file {converted_value}"
                                        "doesn't exist")

        setattr(self, setting_name, converted_value)

    def reset_to_default(self):
        """Reset settings to their default values"""
        default_settings = Settings()
        for field in default_settings.__dataclass_fields__:
            setattr(self, field, getattr(default_settings, field))

    def __str__(self):
        attributes = asdict(self)
        return '\n'.join([f"{key}: {value}" for key, value in attributes.items()])
