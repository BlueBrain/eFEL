"""eFEL Python API functions.

This module provides the user-facing Python API of the eFEL.
The convenience functions defined here call the underlying 'cppcore' library
to hide the lower level API from the user.
All functions in this module can be called as efel.functionname, it is
not necessary to include 'api' as in efel.api.functionname.


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

import os
import numpy

import efel
import efel.cppcore as cppcore

"""
Disabling cppcore importerror override, it confuses users in case the error
is caused by something else
try:
except ImportError:
    six.reraise(ImportError, ImportError(
        '\n'
        'It looks like the efel.cppcore package could not be found.\n'
        'Could it be that you are running the \'import efel\' in a directory '
        'that has a subdirectory called \'efel\' '
        '(like e.g. the eFEL source directory) ?\n'
        'If this is the case, please try to import from another directory.\n'
        'If the issue persists, please create a ticket at '
        'github.com/BlueBrain/eFEL/issues.\n'), sys.exc_info()[2])
"""

_settings = efel.Settings()
_int_settings = {}
_double_settings = {}


def reset():
    """Resets the efel to its initial state

    The user can set certain values in the efel, like the spike threshold.
    These values are persisten. This function will reset these value to their
    original state.
    """

    global _settings
    _settings = efel.Settings()

    setDoubleSetting('spike_skipf', 0.1)
    setIntSetting('max_spike_skip', 2)
    setDoubleSetting('Threshold', _settings.threshold)
    setDoubleSetting('DerivativeThreshold', _settings.derivative_threshold)
    setDoubleSetting('interp_step', 0.1)
    setDoubleSetting('burst_factor', 1.5)
    setDoubleSetting('voltage_base_start_perc', 0.9)
    setDoubleSetting('voltage_base_end_perc', 1.0)
    setDoubleSetting("initial_perc", 0.1)
    setDoubleSetting("min_spike_height", 20.0)

    _initialise()


def setDependencyFileLocation(location):
    """Set the location of the Dependency file

    The eFEL uses 'Dependency' files to let the user define which versions
    of certain features are used to calculate.
    The installation directory of the eFEL contains a default
    'DependencyV5.txt' file. Unless the user wants to change this file,
    it is not necessary to call this function.

    Parameters
    ==========
    location : string
               path to the location of a Dependency file
    """

    global dependencyFileLocation
    if not os.path.exists(location):
        raise Exception(
            "Path to dependency file {%s} doesn't exist" %
            location)
    _settings.dependencyfile_path = location


def getDependencyFileLocation():
    """Get the location of the Dependency file

    The eFEL uses 'Dependency' files to let the user define which versions
    of certain features are used to calculate.
    The installation directory of the eFEL contains a default
    'DependencyV5.txt' file.

    Returns
    =======
    location : string
               path to the location of a Dependency file
    """

    return _settings.dependencyfile_path


def setThreshold(newThreshold):
    """Set the spike detection threshold in the eFEL

    Parameters
    ==========
    threshold : float
                The new spike detection threshold value (in the same units
                as the traces, e.g. mV).
    """
    _settings.threshold = newThreshold
    setDoubleSetting('Threshold', _settings.threshold)


def setDerivativeThreshold(newDerivativeThreshold):
    """Set the threshold for the derivate for detecting the spike onset

    Some featurea use a threshold on dV/dt to calculate the beginning of an
    action potential. This function allows you to set this threshold.

    Parameters
    ==========
    derivative_threshold : float
                The new derivative threshold value (in the same units
                as the traces, e.g. mV/ms).
    """
    _settings.derivative_threshold = newDerivativeThreshold
    setDoubleSetting('DerivativeThreshold', _settings.derivative_threshold)


def getFeatureNames():
    """Return a list with the name of all the available features

    Returns
    =======
    feature_names : list of strings
                    A list that contains all the feature names available in
                    the eFEL. These names can be used in the featureNames
                    argument of e.g. getFeatureValues()
    """

    cppcore.Initialize(_settings.dependencyfile_path, "log")
    feature_names = []
    cppcore.getFeatureNames(feature_names)

    return feature_names


def FeatureNameExists(feature_name):
    """Does a certain feature name exist ?

    Parameters
    ==========
    feature_name : string
                  Name of the feature to check

    Returns
    =======
    FeatureNameExists : bool
                    True if feature_name exists, otherwise False
    """

    return feature_name in getFeatureNames()


def getDistance(trace, featureName, mean, std, trace_check=None):
    """Calculate distance value for a list of traces.

    Parameters
    ==========
    trace : trace dicts
            Trace dict that represents one trace. The dict should have the
            following keys: 'T', 'V', 'stim_start', 'stim_end'
    featureName : string
                  Name of the the features for which to calculate the distance
    mean : float
           Mean to calculate the distance from
    std : float
          Std to scale the distance with
    trace_check : float
          Let the library check if there are spikes outside of stimulus interval

    Returns
    =======
    distance : float
               The absolute number of standard deviation the feature is away
               from the mean. In case of anomalous results a value of '250'
               standard deviations is returned. This can happen if: a feature
               generates an error, there are spikes outside of the stimulus
               interval, the feature returns a NaN, etc.
    """

    _initialise()

    # Next set time, voltage and the stimulus start and end
    for item in list(trace.keys()):
        cppcore.setFeatureDouble(item, [x for x in trace[item]])

    if trace_check is not None:
        return efel.cppcore.getDistance(
            featureName,
            mean,
            std,
            trace_check=1 if trace_check else 0)
    else:
        return efel.cppcore.getDistance(
            featureName,
            mean,
            std)


def _initialise():
    """Set cppcore initial values"""
    cppcore.Initialize(_settings.dependencyfile_path, "log")

    # First set some settings that are used by the feature extraction

    for setting_name, int_setting in list(_int_settings.items()):
        cppcore.setFeatureInt(setting_name, [int_setting])

    for setting_name, double_setting in list(_double_settings.items()):
        cppcore.setFeatureDouble(setting_name, [double_setting])


def setIntSetting(setting_name, new_value):
    """Set a certain integer setting to a new value"""

    _int_settings[setting_name] = new_value


def setDoubleSetting(setting_name, new_value):
    """Set a certain double setting to a new value"""

    _double_settings[setting_name] = new_value


def getFeatureValues(
        traces,
        featureNames,
        parallel_map=None,
        return_list=True,
        raise_warnings=True):
    """Calculate feature values for a list of traces.

    This function is the core of the eFEL API. A list of traces (in the form
    of dictionaries) is passed as argument, together with a list of feature
    names.

    The return value consists of a list of dictionaries, one for each input
    trace. The keys in the dictionaries are the names of the calculated
    features, the corresponding values are lists with the feature values.
    Beware that every feature returns an array of values. E.g. AP_amplitude
    will return a list with the amplitude of every action potential.

    Parameters
    ==========
    traces : list of trace dicts
             Every trace dict represent one trace. The dict should have the
             following keys: 'T', 'V', 'stim_start', 'stim_end'
    feature_names : list of string
                  List with the names of the features to be calculated on all
                  the traces.
    parallel_map : map function
                   Map function to parallelise over the traces. Default is the
                   serial map() function
    return_list: boolean
                 By default the function returns a list of dicts. This
                 optional argument can disable this, so that the result of the
                 parallel_map() is returned. Can be useful for performance
                 reasons when an iterator is preferred.
    raise_warnings: boolean
                    Raise warning when efel c++ returns an error

    Returns
    =======
    feature_values : list of dicts
                     For every input trace a feature value dict is return (in
                     the same order). The dict contains the keys of
                     'feature_names', every key contains a numpy array with
                     the feature values returned by the C++ efel code.
                     The value is None if an error occured during the
                     calculation of the feature.
    """

    if parallel_map is None:
        parallel_map = map

    traces_featurenames = (
        (trace, featureNames, raise_warnings)
        for trace in traces)
    map_result = parallel_map(_get_feature_values_serial, traces_featurenames)

    if return_list:
        return list(map_result)
    else:
        return map_result


def _get_feature_values_serial(trace_featurenames):
    """Single thread of getFeatureValues"""

    trace, featureNames, raise_warnings = trace_featurenames

    featureDict = {}

    if 'stim_start' in trace and 'stim_end' in trace:
        try:
            len(trace['stim_start'])
            len(trace['stim_end'])
        except:
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

    for featureName in featureNames:
        cppcoreFeatureValues = list()
        exitCode = cppcore.getFeature(featureName, cppcoreFeatureValues)

        if exitCode < 0:
            if raise_warnings:
                import warnings
                warnings.warn(
                    "Error while calculating feature %s: %s" %
                    (featureName, cppcore.getgError()),
                    RuntimeWarning)
            featureDict[featureName] = None
        else:
            featureDict[featureName] = numpy.array(cppcoreFeatureValues)

    return featureDict


def getMeanFeatureValues(traces, featureNames, raise_warnings=True):
    """Convenience function that returns the mean values from getFeatureValues()

    Instead of return a list of values for every feature as getFeatureValues()
    does, this function returns per trace one value for every feature, namely
    the mean value.

    Parameters
    ==========
    traces : list of trace dicts
             Every trace dict represent one trace. The dict should have the
             following keys: 'T', 'V', 'stim_start', 'stim_end'
    feature_names : list of string
                    List with the names of the features to be calculated on all
                    the traces.
    raise_warnings: boolean
                    Raise warning when efel c++ returns an error

    Returns
    =======
    feature_values : list of dicts
                     For every input trace a feature value dict is return (in
                     the same order). The dict contains the keys of
                     'feature_names', every key contains the mean of the array
                     that is returned by getFeatureValues()
                     The value is None if an error occured during the
                     calculation of the feature, or if the feature value array
                     was empty.
    """

    featureDicts = getFeatureValues(
        traces,
        featureNames,
        raise_warnings=raise_warnings)
    for featureDict in featureDicts:
        for (key, values) in list(featureDict.items()):
            if values is None or len(values) == 0:
                featureDict[key] = None
            else:
                featureDict[key] = numpy.mean(values)

    return featureDicts

reset()
