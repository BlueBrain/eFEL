"""General tests of eFEL"""

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

import nose.tools as nt

# pylint: disable=W0611, W0612, F0401


def test_import():
    """basic: Test importing of eFEL"""

    # pylint: disable=W0611
    import efel
    # pylint: enable=W0611


def test_version():
    """basic: Test if version number exists"""

    import efel

    nt.assert_true(efel.version is not None)


def test_empty_trace():
    """basic: Testing results for empty trace"""

    import efel
    import numpy

    max_time = 3000.0
    stim_start = 700.0
    stim_end = 2700.0
    dt = 0.02

    time = numpy.arange(0.0, max_time, dt)
    voltage = -80.0 * numpy.ones(len(time))

    trace = {}

    trace['T'] = time
    trace['V'] = voltage
    trace['stim_start'] = [stim_start]
    trace['stim_end'] = [stim_end]

    features = [
        'time_to_last_spike',
        'inv_time_to_first_spike',
        'inv_first_ISI',
        'inv_second_ISI',
        'inv_third_ISI',
        'inv_fourth_ISI',
        'inv_fifth_ISI',
        'inv_last_ISI']

    efel.getFeatureValues([trace], features)

    for feature, value in \
            efel.getFeatureValues([trace], features)[0].iteritems():

        nt.assert_equal(value[0], 0.0)


def test_ISI_log_slope_skip():
    """basic: Test ISI_log_slope_skip"""

    import efel
    import numpy

    stim_start = 31.2
    stim_end = 431.2

    data = numpy.loadtxt(
        'testdata/basic/zero_ISI_log_slope_skip95824004.abf.csv')

    time = data[:, 0]
    voltage = data[:, 1]

    trace = {}

    trace['T'] = time
    trace['V'] = voltage
    trace['stim_start'] = [stim_start]
    trace['stim_end'] = [stim_end]

    features = ['ISI_log_slope_skip']

    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        feature_values = \
            efel.getFeatureValues(
                [trace],
                features)
    nt.assert_equal(feature_values[0]['ISI_log_slope_skip'], None)


def test_mean_frequency1():
    """basic: Test mean_frequency 1"""

    import efel
    import numpy

    stim_start = 500.0
    stim_end = 900.0

    data = numpy.loadtxt('testdata/basic/mean_frequency_1.txt')

    time = data[:, 0]
    voltage = data[:, 1]

    trace = {}

    trace['T'] = time
    trace['V'] = voltage
    trace['stim_start'] = [stim_start]
    trace['stim_end'] = [stim_end]

    features = ['mean_frequency']

    feature_values = \
        efel.getFeatureValues(
            [trace],
            features)
    nt.assert_almost_equal(feature_values[0]['mean_frequency'], 15.2858453)


def test_getDistance1():
    """basic: Test getDistance 1"""

    import efel
    import numpy

    stim_start = 500.0
    stim_end = 900.0

    data = numpy.loadtxt('testdata/basic/mean_frequency_1.txt')

    time = data[:, 0]
    voltage = data[:, 1]

    trace = {}

    trace['T'] = time
    trace['V'] = voltage
    trace['stim_start'] = [stim_start]
    trace['stim_end'] = [stim_end]

    nt.assert_almost_equal(
        3.09045815935,
        efel.getDistance(
            trace,
            'AP_amplitude',
            50,
            10))


def test_APlast_amp():
    """basic: Test APlast_amp"""

    import efel
    import numpy

    stim_start = 500.0
    stim_end = 900.0

    data = numpy.loadtxt('testdata/basic/mean_frequency_1.txt')

    time = data[:, 0]
    voltage = data[:, 1]

    trace = {}

    trace['T'] = time
    trace['V'] = voltage
    trace['stim_start'] = [stim_start]
    trace['stim_end'] = [stim_end]

    features = ['AP_amplitude', 'APlast_amp']

    feature_values = \
        efel.getFeatureValues(
            [trace],
            features)

    APlast_amp = feature_values[0]['APlast_amp'][0]
    AP_amplitude = feature_values[0]['APlast_amp']
    nt.assert_equal(APlast_amp, AP_amplitude[-1])


def test_min_voltage_between_spikes1():
    """basic: Test min_voltage_between_spikes 1"""

    import efel
    import numpy

    stim_start = 500.0
    stim_end = 900.0

    data = numpy.loadtxt('testdata/basic/mean_frequency_1.txt')

    time = data[:, 0]
    voltage = data[:, 1]

    trace = {}

    trace['T'] = time
    trace['V'] = voltage
    trace['stim_start'] = [stim_start]
    trace['stim_end'] = [stim_end]

    features = ['min_voltage_between_spikes', 'peak_indices', 'voltage']

    feature_values = \
        efel.getFeatureValues(
            [trace],
            features)

    peak_indices = feature_values[0]['peak_indices']
    min_voltage_between_spikes = feature_values[
        0]['min_voltage_between_spikes']
    fel_voltage = feature_values[0]['voltage']

    for index, min_voltage_between_spikes_value in zip(
            range(len(peak_indices[:-1])),
            min_voltage_between_spikes):
        nt.assert_almost_equal(
            numpy.min(
                fel_voltage[
                    peak_indices[index]:peak_indices[
                        index +
                        1]]),
            min_voltage_between_spikes_value)


def test_getFeatureNames():
    """basic: Testing getting all feature names"""

    import efel
    import json
    with open('featurenames.json', 'r') as featurenames_json:
        expected_featurenames = json.load(featurenames_json)
    nt.assert_equal(efel.getFeatureNames(), expected_featurenames)


def test_steady_state_voltage_stimend():
    """basic: steady_state_voltage_stimend 1"""

    import efel
    import numpy

    stim_start = 500.0
    stim_end = 900.0

    data = numpy.loadtxt('testdata/basic/mean_frequency_1.txt')

    time = data[:, 0]
    voltage = data[:, 1]

    trace = {}

    trace['T'] = time
    trace['V'] = voltage
    trace['stim_start'] = [stim_start]
    trace['stim_end'] = [stim_end]

    features = ['steady_state_voltage_stimend']

    feature_values = \
        efel.getFeatureValues(
            [trace],
            features)[0]

    stim_duration = stim_end - stim_start
    begin_time = stim_end - 0.1 * stim_duration
    end_time = stim_end
    steady_state_voltage_stimend = [
        numpy.mean(voltage[numpy.where(
            (time <= end_time) & (time >= begin_time)
        )])]

    nt.assert_almost_equal(steady_state_voltage_stimend,
                           feature_values['steady_state_voltage_stimend'])
