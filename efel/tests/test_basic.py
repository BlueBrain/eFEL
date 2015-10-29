# pylint: disable=W0611, W0612, F0401, R0914, C0302

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

import os
from os.path import join as joinp
import nose.tools as nt


_multiprocess_can_split_ = True

testdata_dir = joinp(os.path.dirname(os.path.abspath(__file__)), 'testdata')


def test_import():
    """basic: Test importing of eFEL"""

    # pylint: disable=W0611
    import efel  # NOQA
    # pylint: enable=W0611


def test_version():
    """basic: Test if version number exists"""

    import efel
    efel.reset()

    nt.assert_true(efel.__version__ is not None)


def test_setDependencyFileLocation_wrongpath():
    """basic: Test if setDependencyFileLocation fails when path doesn't exist"""

    import efel
    efel.reset()
    nt.assert_raises(
        Exception,
        efel.setDependencyFileLocation, "thisfiledoesntexist")


def test_nonexisting_feature():
    """basic: Test nonexisting feature"""

    import efel
    efel.reset()

    import numpy
    trace = {}
    trace['T'] = numpy.arange(0, 100, 0.1)
    trace['V'] = numpy.ones(len(trace['T'])) * -80.0
    trace['stim_start'] = [25]
    trace['stim_end'] = [75]

    nt.assert_raises(
        TypeError,
        efel.getFeatureValues,
        [trace],
        ['nonexisting_feature'])


def test_failing_double_feature():
    """basic: Test failing double feature"""

    import efel
    efel.reset()

    import numpy
    trace = {}
    trace['T'] = numpy.arange(0, 100, 0.1)
    trace['V'] = numpy.ones(len(trace['T'])) * -80.0
    trace['stim_start'] = [25]
    trace['stim_end'] = [75]

    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        feature_value = efel.getFeatureValues(
            [trace],
            ['AP_amplitude'])[0]['AP_amplitude']

    nt.assert_equal(feature_value, None)


def test_failing_int_feature():
    """basic: Test failing int feature"""

    import efel
    efel.reset()

    import numpy
    trace = {}
    trace['T'] = numpy.arange(0, 100, 0.1)
    trace['V'] = numpy.ones(len(trace['T'])) * -80.0
    trace['stim_start'] = [25]
    trace['stim_end'] = [75]

    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        feature_value = efel.getFeatureValues(
            [trace],
            ['burst_number'])[0]['burst_number']

    nt.assert_equal(feature_value, None)


def test_empty_trace():
    """basic: Testing results for empty trace"""

    import efel
    efel.reset()

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

    # efel.getFeatureValues([trace], features)

    for feature, value in \
            efel.getFeatureValues([trace], features)[0].items():

        nt.assert_equal(value[0], 0.0)


def test_multiprocessing_traces():
    """basic: Test multiprocessing map"""

    import efel
    efel.reset()
    import numpy

    stim_start = 31.2
    stim_end = 431.2

    test_data_path = joinp(
        testdata_dir,
        'basic',
        'zero_ISI_log_slope_skip95824004.abf.csv')
    data1 = numpy.loadtxt(test_data_path)

    time1 = data1[:, 0]
    voltage1 = data1[:, 1]

    trace1 = {}

    trace1['T'] = time1
    trace1['V'] = voltage1
    trace1['stim_start'] = [stim_start]
    trace1['stim_end'] = [stim_end]

    feature_name = 'peak_time'

    test_data_path = joinp(
        testdata_dir,
        'basic',
        'AP_begin_indices_95810005.abf.csv')
    data2 = numpy.loadtxt(test_data_path)

    voltage2 = data2
    time2 = numpy.arange(len(voltage2)) * 0.1

    trace2 = {}

    trace2['T'] = time2
    trace2['V'] = voltage2
    trace2['stim_start'] = [stim_start]
    trace2['stim_end'] = [stim_end]

    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        feature_values_serial = efel.getFeatureValues(
            [trace1, trace2],
            [feature_name])

    import multiprocessing
    pool = multiprocessing.Pool()

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        feature_values_parallel = efel.getFeatureValues(
            [trace1, trace2],
            [feature_name], parallel_map=pool.map)

    nt.assert_equal(
        list(feature_values_serial[0]['peak_time']),
        list(feature_values_parallel[0]['peak_time']))
    nt.assert_equal(
        list(feature_values_serial[1]['peak_time']),
        list(feature_values_parallel[1]['peak_time']))

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        feature_values_async = efel.getFeatureValues(
            [trace1, trace2],
            [feature_name], parallel_map=pool.map_async, return_list=False)
        nt.assert_is_instance(
            feature_values_async,
            multiprocessing.pool.AsyncResult)


def test_consecutive_traces():
    """basic: Test if features from two different traces give other results"""

    import efel
    efel.reset()
    import numpy

    stim_start = 31.2
    stim_end = 431.2

    test_data_path = joinp(
        testdata_dir,
        'basic',
        'zero_ISI_log_slope_skip95824004.abf.csv')
    data1 = numpy.loadtxt(test_data_path)

    time1 = data1[:, 0]
    voltage1 = data1[:, 1]

    trace1 = {}

    trace1['T'] = time1
    trace1['V'] = voltage1
    trace1['stim_start'] = [stim_start]
    trace1['stim_end'] = [stim_end]

    feature_name = 'peak_time'

    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        feature_values1 = \
            efel.getFeatureValues(
                [trace1],
                [feature_name])

    test_data_path = joinp(
        testdata_dir,
        'basic',
        'AP_begin_indices_95810005.abf.csv')
    data2 = numpy.loadtxt(test_data_path)

    voltage2 = data2
    time2 = numpy.arange(len(voltage2)) * 0.1

    trace2 = {}

    trace2['T'] = time2
    trace2['V'] = voltage2
    trace2['stim_start'] = [stim_start]
    trace2['stim_end'] = [stim_end]

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        feature_values2 = \
            efel.getFeatureValues(
                [trace2],
                [feature_name])

    nt.assert_not_equal(
        len(feature_values1[0][feature_name]),
        len(feature_values2[0][feature_name]))


def test_stimstart_stimend():
    """basic: Test exception when stimstart or stimend are wrong"""

    import efel
    efel.reset()
    import numpy

    stim_start = 500.0
    stim_end = 900.0

    test_data_path = joinp(testdata_dir, 'basic', 'mean_frequency_1.txt')
    data = numpy.loadtxt(test_data_path)

    time = data[:, 0]
    voltage = data[:, 1]

    trace = {}

    trace['T'] = time
    trace['V'] = voltage
    trace['stim_start'] = stim_start
    trace['stim_end'] = stim_end

    features = ['AP_begin_voltage']

    nt.assert_raises(
        Exception,
        efel.getFeatureValues, [trace], features)

    trace['stim_start'] = [stim_end]
    trace['stim_end'] = [stim_start]

    nt.assert_raises(
        Exception,
        efel.getFeatureValues, [trace], features)

    trace['stim_start'] = [stim_start, stim_end]
    trace['stim_end'] = [stim_end]

    nt.assert_raises(
        Exception,
        efel.getFeatureValues, [trace], features)

    del trace['stim_start']

    nt.assert_raises(
        Exception,
        efel.getFeatureValues, [trace], features)


def test_setDerivativeThreshold():
    """basic: Test setDerivativeThreshold"""

    import efel
    efel.reset()
    import numpy

    stim_start = 500.0
    stim_end = 900.0

    test_data_path = joinp(testdata_dir, 'basic', 'mean_frequency_1.txt')
    data = numpy.loadtxt(test_data_path)

    time = data[:, 0]
    voltage = data[:, 1]

    trace = {}

    trace['T'] = time
    trace['V'] = voltage
    trace['stim_start'] = [stim_start]
    trace['stim_end'] = [stim_end]

    features = ['AP_begin_voltage']

    feature_values = \
        efel.getFeatureValues(
            [trace],
            features)
    AP_begin_voltage_orig = feature_values[0]['AP_begin_voltage'][1]

    efel.setDerivativeThreshold(5)
    feature_values = \
        efel.getFeatureValues(
            [trace],
            features)
    AP_begin_voltage = feature_values[0]['AP_begin_voltage'][1]
    nt.assert_almost_equal(AP_begin_voltage, -51.6400489995987)
    nt.assert_not_equal(AP_begin_voltage, AP_begin_voltage_orig)


def test_ISI_log_slope_skip():
    """basic: Test ISI_log_slope_skip"""

    import efel
    efel.reset()
    import numpy

    stim_start = 31.2
    stim_end = 431.2

    test_data_path = joinp(
        testdata_dir,
        'basic',
        'zero_ISI_log_slope_skip95824004.abf.csv')
    data = numpy.loadtxt(test_data_path)

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


def test_AP_begin_indices1():
    """basic: Test AP_begin_indices 1"""

    import efel
    efel.reset()
    import numpy

    stim_start = 31.2
    stim_end = 431.2

    test_data_path = joinp(
        testdata_dir,
        'basic',
        'AP_begin_indices_95810005.abf.csv')
    voltage = numpy.loadtxt(test_data_path)

    time = numpy.arange(len(voltage)) * 0.1

    trace = {}

    trace['V'] = voltage
    trace['T'] = time
    trace['stim_start'] = [stim_start]
    trace['stim_end'] = [stim_end]

    features = [
        'AP_begin_indices',
        'AP_amplitude',
        'peak_time',
        'AP_duration_half_width']

    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        feature_values = \
            efel.getFeatureValues(
                [trace],
                features)
    # Make sure the amount of peak_times, AP_begin_indices and AP_amplitude is
    # the same in this trace
    # There was originally an issue in this case due to the 'width' value
    # in AP_begin_indices, which caused a segmentation fault
    nt.assert_equal(
        len(feature_values[0]['AP_begin_indices']),
        len(feature_values[0]['AP_amplitude']))
    nt.assert_equal(
        len(feature_values[0]['AP_begin_indices']),
        len(feature_values[0]['peak_time']))
    nt.assert_equal(
        len(feature_values[0]['AP_begin_indices']),
        len(feature_values[0]['AP_duration_half_width']))


def test_mean_frequency1():
    """basic: Test mean_frequency 1"""

    import efel
    efel.reset()
    import numpy

    stim_start = 500.0
    stim_end = 900.0

    test_data_path = joinp(testdata_dir, 'basic', 'mean_frequency_1.txt')
    data = numpy.loadtxt(test_data_path)

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
    nt.assert_almost_equal(feature_values[0]['mean_frequency'][0], 15.2858453)


def test_ap_amplitude_outside_stim():
    """basic: Test AP amplitude with spike outside stim"""

    import efel
    efel.reset()
    import numpy

    stim_start = 700.0
    stim_end = 2700.0

    test_data_path = joinp(testdata_dir, 'basic', 'spike_outside_stim.txt')
    data = numpy.loadtxt(test_data_path)

    time = data[:, 0]
    voltage = data[:, 1]

    trace = {}

    trace['T'] = time
    trace['V'] = voltage
    trace['stim_start'] = [stim_start]
    trace['stim_end'] = [stim_end]

    features = ['AP_amplitude', 'peak_time']

    feature_values = \
        efel.getFeatureValues(
            [trace],
            features)

    # Make sure AP_amplitude doesn't pick up the two spikes outside of
    # the stimulus
    # (which are present in peak_time)
    nt.assert_equal(
        len(feature_values[0]['AP_amplitude']) + 2,
        len(feature_values[0]['peak_time']))


def test_ap_amplitude_from_voltagebase1():
    """basic: Test AP_amplitude_from_voltagebase 1"""

    import efel
    efel.reset()
    import numpy

    stim_start = 500.0
    stim_end = 900.0

    test_data_path = joinp(testdata_dir, 'basic', 'mean_frequency_1.txt')
    data = numpy.loadtxt(test_data_path)

    time = data[:, 0]
    voltage = data[:, 1]

    trace = {}

    trace['T'] = time
    trace['V'] = voltage
    trace['stim_start'] = [stim_start]
    trace['stim_end'] = [stim_end]

    features = ['AP_amplitude_from_voltagebase',
                'peak_voltage', 'voltage_base']

    feature_values = \
        efel.getFeatureValues(
            [trace],
            features)

    voltage_base = feature_values[0]['voltage_base'][0]
    for peak_voltage, ap_amplitude_from_voltagebase in zip(
            feature_values[0]['peak_voltage'],
            feature_values[0]['AP_amplitude_from_voltagebase']):
        nt.assert_almost_equal(peak_voltage - voltage_base,
                               ap_amplitude_from_voltagebase)


def test_voltagebase1():
    """basic: Test voltagebase 1"""

    import efel
    efel.reset()
    import numpy

    stim_start = 500.0
    stim_end = 900.0

    test_data_path = joinp(testdata_dir, 'basic', 'mean_frequency_1.txt')
    data = numpy.loadtxt(test_data_path)

    time = data[:, 0]
    voltage = data[:, 1]

    trace = {}

    trace['T'] = time
    trace['V'] = voltage
    trace['stim_start'] = [stim_start]
    trace['stim_end'] = [stim_end]

    features = ['voltage_base']

    feature_values = \
        efel.getFeatureValues(
            [trace],
            features)

    voltage_base = numpy.mean(voltage[numpy.where(
        (time >= 0.9 * stim_start) & (time <= stim_start))])

    nt.assert_almost_equal(voltage_base, feature_values[0]['voltage_base'][0])


def test_getDistance1():
    """basic: Test getDistance 1"""

    import efel
    efel.reset()
    import numpy

    stim_start = 500.0
    stim_end = 900.0

    test_data_path = joinp(testdata_dir, 'basic', 'mean_frequency_1.txt')
    data = numpy.loadtxt(test_data_path)

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
    efel.reset()
    import numpy

    stim_start = 500.0
    stim_end = 900.0

    test_data_path = joinp(testdata_dir, 'basic', 'mean_frequency_1.txt')
    data = numpy.loadtxt(test_data_path)

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


def test_spikecount1():
    """basic: Test Spikecount 1"""

    import efel
    efel.reset()
    import numpy

    stim_start = 500.0
    stim_end = 900.0

    test_data_path = joinp(testdata_dir, 'basic', 'mean_frequency_1.txt')
    data = numpy.loadtxt(test_data_path)

    time = data[:, 0]
    voltage = data[:, 1]

    trace = {}

    trace['T'] = time
    trace['V'] = voltage
    trace['stim_start'] = [stim_start]
    trace['stim_end'] = [stim_end]

    features = ['peak_indices', 'Spikecount']

    feature_values = \
        efel.getFeatureValues(
            [trace],
            features)

    peak_indices = feature_values[0]['peak_indices']
    spikecount = feature_values[0]['Spikecount'][0]
    nt.assert_equal(len(peak_indices), spikecount)


def test_spikecount_libv4peakindices():
    """basic: Test Spikecount in combination with LibV4 peak_indices"""

    import efel
    efel.reset()
    import numpy

    stim_start = 500.0
    stim_end = 900.0

    data = numpy.loadtxt(joinp(testdata_dir, 'basic', 'mean_frequency_1.txt'))

    time = data[:, 0]
    voltage = data[:, 1]

    trace = {}

    trace['T'] = time
    trace['V'] = voltage
    trace['stim_start'] = [stim_start]
    trace['stim_end'] = [stim_end]

    features = ['peak_indices', 'Spikecount']

    test_peak = joinp(os.path.dirname(os.path.abspath(__file__)),
                      'DependencyV5_LibV4peakindices.txt')
    efel.setDependencyFileLocation(test_peak)

    feature_values = \
        efel.getFeatureValues(
            [trace],
            features)

    peak_indices = feature_values[0]['peak_indices']
    spikecount = feature_values[0]['Spikecount'][0]
    nt.assert_equal(len(peak_indices), 5)
    nt.assert_equal(len(peak_indices), spikecount)


def test_spikecount2():
    """basic: Test Spikecount 2: test empty trace"""

    import efel
    efel.reset()
    import numpy

    stim_start = 500.0
    stim_end = 900.0

    time = numpy.arange(0, 1000.0, 0.1)
    voltage = numpy.ones(len(time)) * -80.0

    trace = {}

    trace['T'] = time
    trace['V'] = voltage
    trace['stim_start'] = [stim_start]
    trace['stim_end'] = [stim_end]

    features = ['Spikecount']

    feature_values = \
        efel.getFeatureValues(
            [trace],
            features)

    spikecount = feature_values[0]['Spikecount'][0]
    nt.assert_equal(spikecount, 0)


def test_min_voltage_between_spikes1():
    """basic: Test min_voltage_between_spikes 1"""

    import efel
    efel.reset()
    import numpy

    stim_start = 500.0
    stim_end = 900.0

    test_data_path = joinp(testdata_dir, 'basic', 'mean_frequency_1.txt')
    data = numpy.loadtxt(test_data_path)

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
            list(range(len(peak_indices[:-1]))),
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
    efel.reset()
    import json

    test_data_path = joinp(testdata_dir, '..', 'featurenames.json')
    with open(test_data_path, 'r') as featurenames_json:
        expected_featurenames = json.load(featurenames_json)
    nt.assert_equal(efel.getFeatureNames(), expected_featurenames)


def test_steady_state_voltage1():
    """basic: steady_state_voltage 1"""

    import efel
    efel.reset()
    import numpy

    stim_start = 500.0
    stim_end = 900.0

    test_data_path = joinp(testdata_dir, 'basic', 'mean_frequency_1.txt')
    data = numpy.loadtxt(test_data_path)

    time = data[:, 0]
    voltage = data[:, 1]

    trace = {}

    trace['T'] = time
    trace['V'] = voltage
    trace['stim_start'] = [stim_start]
    trace['stim_end'] = [stim_end]

    features = ['steady_state_voltage']

    feature_values = \
        efel.getFeatureValues(
            [trace],
            features)[0]

    begin_time = stim_end
    end_time = max(time)
    steady_state_voltage = numpy.mean(voltage[numpy.where(
        (time <= end_time) & (time > begin_time)
    )])

    nt.assert_almost_equal(steady_state_voltage,
                           feature_values['steady_state_voltage'][0])


def test_steady_state_voltage_stimend():
    """basic: steady_state_voltage_stimend 1"""

    import efel
    efel.reset()
    import numpy

    stim_start = 500.0
    stim_end = 900.0

    test_data_path = joinp(testdata_dir, 'basic', 'mean_frequency_1.txt')
    data = numpy.loadtxt(test_data_path)

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
    steady_state_voltage_stimend = numpy.mean(voltage[numpy.where(
        (time < end_time) & (time >= begin_time)
    )])

    nt.assert_almost_equal(steady_state_voltage_stimend,
                           feature_values['steady_state_voltage_stimend'][0])


def decay_time_constant_after_stim(time, voltage, interval_start,
                                   interval_end, stim_start, stim_end):
    '''numpy implementation'''
    import numpy

    def get_index(ts, t):
        """get_index"""
        return next(i for i in range(len(ts)) if ts[i] >= t)

    interval_indices = numpy.where(
        (time >= interval_start) & (time < interval_end))
    stim_start_index = get_index(time, stim_start)
    interval_time = time[interval_indices] - stim_end
    interval_voltage = abs(
        voltage[interval_indices] -
        voltage[stim_start_index])

    # fit
    log_interval_voltage = numpy.log(interval_voltage)
    slope, _ = numpy.polyfit(interval_time, log_interval_voltage, 1)

    tau = -1. / slope
    return abs(tau)


def test_decay_time_constant_after_stim1():
    """basic: decay_time_constant_after_stim 1"""

    import efel
    efel.reset()
    import numpy

    stim_start = 500.0
    stim_end = 900.0

    test_data_path = joinp(testdata_dir, 'basic', 'mean_frequency_1.txt')
    data = numpy.loadtxt(test_data_path)

    time = data[:, 0]
    voltage = data[:, 1]

    trace = {
        'T': time,
        'V': voltage,
        'stim_start': [stim_start],
        'stim_end': [stim_end],
    }

    features = ['decay_time_constant_after_stim']

    feature_values = efel.getFeatureValues([trace], features)[0]

    expected = decay_time_constant_after_stim(
        trace['T'],
        trace['V'],
        stim_end + 1.0,
        stim_end + 10.0,
        trace['stim_start'][0],
        trace['stim_end'][0])

    nt.assert_almost_equal(
        expected,
        feature_values['decay_time_constant_after_stim'][0])


def test_decay_time_constant_after_stim2():
    """basic: decay_time_constant_after_stim 2"""

    import efel
    efel.reset()
    import numpy

    stim_start = 100.0
    stim_end = 1000.0

    test_data_path = joinp(testdata_dir, 'basic', 'tau20.0.csv')
    data = numpy.loadtxt(test_data_path)

    time = data[:, 0]
    voltage = data[:, 1]

    trace = {
        'T': time,
        'V': voltage,
        'stim_start': [stim_start],
        'stim_end': [stim_end],
        'decay_start_after_stim': [1.0],
        'decay_end_after_stim': [10.0]
    }

    features = ['decay_time_constant_after_stim']

    feature_values = efel.getFeatureValues([trace], features)[0]

    nt.assert_almost_equal(
        20.0,
        feature_values['decay_time_constant_after_stim'][0], places=1)


def test_getmeanfeaturevalues():
    """basic: Test getMeanFeatureValues"""

    import efel
    efel.reset()
    import numpy

    stim_start = 500.0
    stim_end = 900.0

    test_data_path = joinp(testdata_dir, 'basic', 'mean_frequency_1.txt')
    data = numpy.loadtxt(test_data_path)

    time = data[:, 0]
    voltage = data[:, 1]

    trace = {}

    trace['T'] = time
    trace['V'] = voltage
    trace['stim_start'] = [stim_start]
    trace['stim_end'] = [stim_end]

    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        feature_values = \
            efel.getFeatureValues(
                [trace],
                ['AP_amplitude', 'BPAPHeightLoc1'])

        mean_feature_values = efel.getMeanFeatureValues(
            [trace], [
                'AP_amplitude', 'BPAPHeightLoc1'])

    nt.assert_equal(numpy.mean(feature_values[0]['AP_amplitude']),
                    mean_feature_values[0]['AP_amplitude'])
