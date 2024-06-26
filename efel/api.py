"""eFEL Python API functions.

This module provides the user-facing Python API of eFEL.
The convenience functions defined here call the underlying 'cppcore' library
to hide the lower level API from the user.


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
# pylint: disable=W0602,W0603,W0702, F0401, W0612, R0912
from __future__ import annotations

from pathlib import Path
from typing import Callable, Iterator, Literal, overload
from typing_extensions import deprecated
import numpy as np

import efel
import efel.cppcore as cppcore

import efel.pyfeatures as pyfeatures
from efel.pyfeatures.pyfeatures import get_cpp_feature

_settings = efel.Settings()


def set_setting(setting_name: str, new_value: int | float | str) -> None:
    """Set a certain setting to a new value.

    Args:
        setting_name: Name of the setting to change.
        new_value: New value for the setting.
    """
    _settings.set_setting(setting_name, new_value)


def get_settings() -> efel.Settings:
    """Returns the current settings of eFEL."""
    return _settings


def reset():
    """Resets the efel settings to their default values.
    see :py:func:`efel.Settings`
    """
    global _settings
    _settings = efel.Settings()
    _settings.reset_to_default()
    _initialise()


@deprecated("Use `set_setting('dependencyfile_path', location)` instead")
def set_dependency_file_location(location: str | Path) -> None:
    """Sets the location of the Dependency file.

    eFEL uses 'Dependency' files to let the user define versions of features to use.
    The installation directory of eFEL contains a default 'DependencyV5.txt' file.
    Unless users want to change this file, it is not necessary to call this function.
    Modifying the Dependency file can be useful in debugging.

    Args:
        location: Path to the location of a Dependency file.

    Raises:
        FileNotFoundError: If the path to the dependency file doesn't exist.
    """
    set_setting('dependencyfile_path', str(location))


def get_dependency_file_location() -> str:
    """Gets the location of the Dependency file.

    Returns:
        Path to the location of a Dependency file.
    """
    return _settings.dependencyfile_path


@deprecated("Use `set_setting('Threshold', new_threshold)` instead")
def set_threshold(new_threshold: float) -> None:
    """Set the spike detection threshold in the eFEL, default -20.0

    Args:
        new_threshold: The new spike detection threshold value (in the same units
                       as the traces, e.g. mV).
    """
    set_setting('Threshold', new_threshold)


@deprecated("Use `set_setting('DerivativeThreshold', "
            "new_derivative_threshold)` instead")
def set_derivative_threshold(new_derivative_threshold: float) -> None:
    """Set the threshold for the derivative for detecting the spike onset.

    Some features use a threshold on dV/dt to calculate the beginning of an
    action potential. This function allows you to set this threshold.

    Args:
        new_derivative_threshold: The new derivative threshold value (in the same units
                                  as the traces, e.g. mV/ms).
    """
    set_setting('DerivativeThreshold', new_derivative_threshold)


def get_feature_names() -> list[str]:
    """Return a list with the name of all the available features

    Returns:
        A list that contains all the feature names available in
        the eFEL. These names can be used in the feature_names
        argument of e.g. get_feature_values()
    """
    cppcore.Initialize(_settings.dependencyfile_path, "log")
    feature_names: list[str] = []
    cppcore.getFeatureNames(feature_names)

    feature_names += pyfeatures.all_pyfeatures

    return feature_names


def feature_name_exists(feature_name: str) -> bool:
    """Returns True if the feature name exists in eFEL, False otherwise."""
    return feature_name in get_feature_names()


def _get_feature(feature_name: str, raise_warnings=False) -> np.ndarray | None:
    """Get feature value, decide to use python or cpp"""
    if feature_name in pyfeatures.all_pyfeatures:
        return get_py_feature(feature_name)
    else:
        return get_cpp_feature(feature_name, raise_warnings=raise_warnings)


def get_distance(
        trace: dict,
        feature_name: str,
        mean: float,
        std: float,
        trace_check: bool = True,
        error_dist: float = 250) -> float:
    """Calculate distance value for a list of traces.

    Args:
        trace: Trace dict that represents one trace. The dict should have the
               following keys: 'T', 'V', 'stim_start', 'stim_end'
        feature_name: Name of the the features for which to calculate the distance
        mean: Mean to calculate the distance from
        std: Std to scale the distance with
        trace_check: Let the library check if there are spikes outside of stimulus
                     interval, default is True
        error_dist: Distance returned when error, default is 250

    Returns:
        The absolute number of standard deviation the feature is away
        from the mean. In case of anomalous results a value of
        'error_dist' standard deviations is returned.
        This can happen if: a feature generates an error, there are
        spikes outside of the stimulus interval, the feature returns
        a NaN, etc.
    """
    _initialise()

    # Next set time, voltage and the stimulus start and end
    for item in list(trace.keys()):
        cppcore.setFeatureDouble(item, [x for x in trace[item]])

    if trace_check:
        trace_check_success = get_feature_values(
            [trace], ['trace_check'], None, True, True)[0]

        if trace_check_success["trace_check"] is None:
            return error_dist

    feature_values = _get_feature(feature_name)

    distance = 0.0
    if feature_values is None or len(feature_values) < 1:
        return error_dist
    else:
        # Am not using anything more fancy to avoid breaking exact
        # reproducibility of legacy C++ code
        for feature_value in feature_values:
            distance += abs(feature_value - mean)

        distance = distance / std / len(feature_values)

        # Check for NaN
        if distance != distance:
            return error_dist

        return distance


def _initialise() -> None:
    """Set cppcore initial values."""
    cppcore.Initialize(_settings.dependencyfile_path, "log")
    # flush the GErrorString from previous runs by calling getgError()
    cppcore.getgError()

    # Set the settings in the cppcore
    settings_attrs = vars(_settings)
    for setting_name, setting_value in settings_attrs.items():
        if isinstance(setting_value, bool):
            setting_value = int(setting_value)
        if isinstance(setting_value, int):
            cppcore.setFeatureInt(setting_name, [setting_value])
        elif isinstance(setting_value, float):
            if isinstance(setting_value, list):
                cppcore.setFeatureDouble(setting_name, setting_value)
            else:
                cppcore.setFeatureDouble(setting_name, [setting_value])
        elif isinstance(setting_value, str):
            cppcore.setFeatureString(setting_name, setting_value)


@deprecated("Use `set_setting()` instead")
def set_int_setting(setting_name: str, new_value: int) -> None:
    """Set a certain integer setting to a new value."""
    set_setting(setting_name, new_value)


@deprecated("Use `set_setting()` instead")
def set_double_setting(setting_name: str, new_value: float) -> None:
    """Set a certain double setting to a new value."""
    set_setting(setting_name, new_value)


@deprecated("Use `set_setting()` instead")
def set_str_setting(setting_name: str, new_value: str) -> None:
    """Set a certain string setting to a new value."""
    set_setting(setting_name, new_value)


@overload
def get_feature_values(
    traces: list[dict],
    feature_names: list[str],
    parallel_map: Callable | None,
    return_list: Literal[True],
    raise_warnings: bool = True,
) -> list:
    ...


@overload
def get_feature_values(
    traces: list[dict],
    feature_names: list[str],
    parallel_map: Callable | None,
    return_list: Literal[False],
    raise_warnings: bool = True,
) -> Iterator:
    ...


def get_feature_values(
    traces: list[dict],
    feature_names: list[str],
    parallel_map: Callable | None = None,
    return_list: bool = True,
    raise_warnings: bool = True,
) -> list | Iterator:
    """Calculate feature values for a list of traces.

    This function is the core of eFEL API. A list of traces (in the form
    of dictionaries) is passed as argument, together with a list of feature
    names.

    The return value consists of a list of dictionaries, one for each input
    trace. The keys in the dictionaries are the names of the calculated
    features, the corresponding values are lists with the feature values.
    Beware that every feature returns an array of values. E.g. AP_amplitude
    will return a list with the amplitude of every action potential.

    Args:
        traces: Every trace dict represents one trace. The dict should have the
                following keys: 'T', 'V', 'stim_start', 'stim_end'
        feature_names: List with the names of the features to be calculated on all
                       the traces.
        parallel_map: Map function to parallelise over the traces. Default is the
                      serial map() function
        return_list: By default the function returns a list of dicts. This
                     optional argument can disable this, so that the result of the
                     parallel_map() is returned. Can be useful for performance
                     reasons when an iterator is preferred.
        raise_warnings: Raise warning when efel c++ returns an error

    Returns:
        For every input trace a feature value dict is returned (in
        the same order). The dict contains the keys of
        'feature_names', every key contains a numpy array with
        the feature values returned by the C++ efel code.
        The value is None if an error occured during the
        calculation of the feature.
    """

    if parallel_map is None:
        parallel_map = map

    traces_featurenames = ((trace, feature_names, raise_warnings) for trace in traces)
    map_result = parallel_map(_get_feature_values_serial, traces_featurenames)

    if return_list:
        return list(map_result)
    else:
        return map_result


def get_py_feature(feature_name: str) -> np.ndarray | None:
    """Return values of the given feature name."""
    return getattr(pyfeatures, feature_name)()


def _get_feature_values_serial(
    trace_featurenames: tuple[dict, list[str], bool]
) -> dict:
    """Single process of get_feature_values."""
    trace, feature_names, raise_warnings = trace_featurenames
    result = {}

    if 'stim_start' in trace and 'stim_end' in trace:
        try:
            len(trace['stim_start'])
            len(trace['stim_end'])
        except BaseException:
            raise Exception('Unable to determine length of stim_start or '
                            'stim_end, are you sure these are lists ?')

        if len(trace['stim_start']) == 1 and len(trace['stim_end']) == 1:
            if trace['stim_end'][0] <= trace['stim_start'][0]:
                raise Exception(
                    'stim_end needs to be larger than '
                    'stim_start:\nstim_start=%f stim_end=%f' %
                    (trace['stim_start'][0], trace['stim_end'][0]))
        else:
            raise Exception(
                'stim_start and stim_end in the trace '
                'dictionary need to be lists of exactly 1 element')

    else:
        raise Exception('stim_start or stim_end missing from trace')

    _initialise()

    # Next set time, voltage and the stimulus start and end
    for item in list(trace.keys()):
        cppcore.setFeatureDouble(item, [x for x in trace[item]])

    for feature_name in feature_names:
        result[feature_name] = _get_feature(
            feature_name, raise_warnings=raise_warnings)

    return result


def get_mean_feature_values(
        traces: list[dict],
        feature_names: list[str],
        raise_warnings: bool = True) -> list[dict]:
    """Convenience function that returns mean values from get_feature_values()

    Instead of return a list of values for every feature as get_feature_values()
    does, this function returns per trace one value for every feature, namely
    the mean value.

    Args:
        traces: Every trace dict represents one trace. The dict should have the
                following keys: 'T', 'V', 'stim_start', 'stim_end'
        feature_names: List with the names of the features to be calculated on all
                       the traces.
        raise_warnings: Raise warning when efel c++ returns an error

    Returns:
        For every input trace a feature value dict is returned (in
        the same order). The dict contains the keys of
        'feature_names', every key contains the mean of the array
        that is returned by get_feature_values()
        The value is None if an error occured during the
        calculation of the feature, or if the feature value array
        was empty.
    """
    featureDicts = get_feature_values(
        traces,
        feature_names,
        parallel_map=None,
        return_list=True,
        raise_warnings=raise_warnings)
    for featureDict in featureDicts:
        for (key, values) in list(featureDict.items()):
            if values is None or len(values) == 0:
                featureDict[key] = None
            else:
                featureDict[key] = np.mean(values)

    return featureDicts


def register_feature(feature_function: Callable):
    """Register a new feature."""
    pyfeatures.all_pyfeatures.append(feature_function.__name__)
    setattr(pyfeatures, feature_function.__name__, feature_function)


reset()


# Deprecated functions
@deprecated("Use set_threshold instead")
def setThreshold(newThreshold: float) -> None:
    set_threshold(newThreshold)


@deprecated("Use set_derivative_threshold instead")
def setDerivativeThreshold(newDerivativeThreshold: float) -> None:
    set_derivative_threshold(newDerivativeThreshold)


@deprecated("Use get_feature_names instead")
def getFeatureNames() -> list[str]:
    return get_feature_names()


@deprecated("Use feature_name_exists instead")
def FeatureNameExists(feature_name: str) -> bool:
    return feature_name_exists(feature_name)


@deprecated("Use get_distance instead")
def getDistance(
        trace,
        featureName,
        mean,
        std,
        trace_check=True,
        error_dist=250) -> float:
    return get_distance(trace, featureName, mean, std, trace_check, error_dist)


@deprecated("Use set_int_setting instead")
def setIntSetting(setting_name: str, new_value: int) -> None:
    set_int_setting(setting_name, new_value)


@deprecated("Use set_double_setting instead")
def setDoubleSetting(setting_name: str, new_value: float) -> None:
    set_double_setting(setting_name, new_value)


@deprecated("Use set_str_setting instead")
def setStrSetting(setting_name: str, new_value: str) -> None:
    set_str_setting(setting_name, new_value)


@deprecated("Use get_feature_values instead")
def getFeatureValues(
        traces,
        featureNames,
        parallel_map=None,
        return_list=True,
        raise_warnings=True):
    return get_feature_values(
        traces, featureNames, parallel_map, return_list, raise_warnings)


@deprecated("Use get_mean_feature_values instead")
def getMeanFeatureValues(
        traces,
        featureNames,
        raise_warnings=True):
    return get_mean_feature_values(traces, featureNames, raise_warnings)


@deprecated("Use get_dependency_file_location instead")
def getDependencyFileLocation() -> str:
    return get_dependency_file_location()


@deprecated("Use set_dependency_file_location instead")
def setDependencyFileLocation(location: str | Path) -> None:
    return set_dependency_file_location(location)
