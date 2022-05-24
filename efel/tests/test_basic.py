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
import numpy

import pytest

_multiprocess_can_split_ = True

testdata_dir = os.path.join(
    os.path.dirname(
        os.path.abspath(__file__)),
    'testdata')


meanfrequency1_url = 'file://%s' % os.path.join(os.path.abspath(testdata_dir),
                                                'basic',
                                                'mean_frequency_1.txt')

ahptest1_url = 'file://%s' % os.path.join(os.path.abspath(testdata_dir),
                                          'basic',
                                          'ahptest_1.txt')

tau20_0_url = 'file://%s' % os.path.join(os.path.abspath(testdata_dir),
                                         'basic',
                                         'tau20.0.csv')

spikeoutsidestim_url = 'file://%s' % os.path.join(
    os.path.abspath(testdata_dir),
    'basic',
    'spike_outside_stim.txt')

sagtrace1_url = 'file://%s' % os.path.join(os.path.abspath(testdata_dir),
                                           'basic',
                                           'sagtrace_1.txt')

zeroISIlog1_url = 'file://%s' % os.path.join(os.path.abspath(testdata_dir),
                                             'basic',
                                             'zero_ISI_log_slope_skip'
                                             '95824004.abf.csv')

derivwindow1_url = 'file://%s' % os.path.join(os.path.abspath(testdata_dir),
                                              'basic',
                                              'derivwindow.txt')

dendriticAP_url = 'file://%s' % os.path.join(os.path.abspath(testdata_dir),
                                             'basic',
                                             'dendritic_AP.txt')


def load_data(data_name, interp=False, interp_dt=0.1):
    """Load data file"""

    import efel

    trace = {}

    if data_name == 'mean_frequency1':
        stim_start = 500.0
        stim_end = 900.0

        time = efel.io.load_fragment('%s#col=1' % meanfrequency1_url)
        voltage = efel.io.load_fragment('%s#col=2' % meanfrequency1_url)
    elif data_name == 'tau20.0':
        stim_start = 100.0
        stim_end = 1000.0

        time = efel.io.load_fragment('%s#col=1' % tau20_0_url)
        voltage = efel.io.load_fragment('%s#col=2' % tau20_0_url)

        trace['decay_start_after_stim'] = [1.0]
        trace['decay_end_after_stim'] = [10.0]
    else:
        raise ValueError('Unknown data set name: %s' % data_name)

    trace['T'] = time
    trace['V'] = voltage
    trace['stim_start'] = [stim_start]
    trace['stim_end'] = [stim_end]

    if interp:
        time, voltage = interpolate(time, voltage, interp_dt)

    return trace, time, voltage, stim_start, stim_end


def test_import():
    """basic: Test importing of eFEL"""

    # pylint: disable=W0611
    import efel  # NOQA
    # pylint: enable=W0611


def test_version():
    """basic: Test if version number exists"""

    import efel
    efel.reset()

    assert efel.__version__ is not None


def test_setDependencyFileLocation_wrongpath():
    """basic: Test if setDependencyFileLocation fails if path doesn't exist"""

    import efel
    efel.reset()
    pytest.raises(
        Exception,
        efel.setDependencyFileLocation, "thisfiledoesntexist")


def test_nonexisting_feature():
    """basic: Test nonexisting feature"""

    import efel
    efel.reset()

    trace = {}
    trace['T'] = numpy.arange(0, 100, 0.1)
    trace['V'] = numpy.ones(len(trace['T'])) * -80.0
    trace['stim_start'] = [25]
    trace['stim_end'] = [75]

    pytest.raises(
        TypeError,
        efel.getFeatureValues,
        [trace],
        ['nonexisting_feature'])


def test_failing_double_feature():
    """basic: Test failing double feature"""

    import efel
    efel.reset()

    trace = {}
    trace['T'] = numpy.arange(0, 100, 0.1)
    trace['V'] = numpy.ones(len(trace['T'])) * -80.0
    trace['stim_start'] = [25]
    trace['stim_end'] = [75]

    feature_value = efel.getFeatureValues(
        [trace],
        ['AP_amplitude'], raise_warnings=False)[0]['AP_amplitude']

    assert feature_value is None


def test_raise_warnings():
    """basic: Test raise_warnings"""

    import efel
    efel.reset()

    trace = {}
    trace['T'] = numpy.arange(0, 100, 0.1)
    trace['V'] = numpy.ones(len(trace['T'])) * -80.0
    trace['stim_start'] = [25]
    trace['stim_end'] = [75]

    import warnings

    with warnings.catch_warnings(record=True) as warning:
        warnings.simplefilter("always")
        feature_value = efel.getFeatureValues(
            [trace],
            ['AP_amplitude'])[0]['AP_amplitude']

        assert feature_value is None
        assert len(warning) == 1
        assert ("Error while calculating feature AP_amplitude" in
                str(warning[0].message))

    with warnings.catch_warnings(record=True) as warning:
        warnings.simplefilter("always")
        feature_value = efel.getFeatureValues(
            [trace],
            ['AP_amplitude'], raise_warnings=False)[0]['AP_amplitude']

        assert feature_value is None
        assert len(warning) == 0


def test_failing_int_feature():
    """basic: Test failing int feature"""

    import efel
    efel.reset()

    trace = {}
    trace['T'] = numpy.arange(0, 100, 0.1)
    trace['V'] = numpy.ones(len(trace['T'])) * -80.0
    trace['stim_start'] = [25]
    trace['stim_end'] = [75]

    feature_value = efel.getFeatureValues(
        [trace],
        ['burst_number'], raise_warnings=False)[0]['burst_number']

    assert feature_value is None


def test_empty_trace():
    """basic: Test results for empty trace"""

    import efel
    efel.reset()

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

        assert value[0] == 0.0


def test_multiprocessing_traces():
    """basic: Test multiprocessing map"""
    import efel
    efel.reset()

    stim_start = 31.2
    stim_end = 431.2

    time1 = efel.io.load_fragment('%s#col=1' % zeroISIlog1_url)
    voltage1 = efel.io.load_fragment('%s#col=2' % zeroISIlog1_url)

    trace1 = {}

    trace1['T'] = time1
    trace1['V'] = voltage1
    trace1['stim_start'] = [stim_start]
    trace1['stim_end'] = [stim_end]

    feature_name = 'peak_time'

    test_data_path = os.path.join(
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

    feature_values_serial = efel.getFeatureValues(
        [trace1, trace2],
        [feature_name], raise_warnings=False)

    efel.reset()
    import multiprocessing
    pool = multiprocessing.Pool()

    feature_values_parallel = efel.getFeatureValues(
        [trace1, trace2],
        [feature_name], parallel_map=pool.map, raise_warnings=False)

    assert (
        list(feature_values_serial[0]['peak_time']) ==
        list(feature_values_parallel[0]['peak_time']))
    assert (
        list(feature_values_serial[1]['peak_time']) ==
        list(feature_values_parallel[1]['peak_time']))

    feature_values_async = efel.getFeatureValues(
        [trace1, trace2], [feature_name], parallel_map=pool.map_async,
        return_list=False, raise_warnings=False)
    assert isinstance(
        feature_values_async,
        multiprocessing.pool.MapResult)


def test_consecutive_traces():
    """basic: Test if features from two different traces give other results"""

    import efel
    efel.reset()

    stim_start = 31.2
    stim_end = 431.2

    time1 = efel.io.load_fragment('%s#col=1' % zeroISIlog1_url)
    voltage1 = efel.io.load_fragment('%s#col=2' % zeroISIlog1_url)

    trace1 = {}

    trace1['T'] = time1
    trace1['V'] = voltage1
    trace1['stim_start'] = [stim_start]
    trace1['stim_end'] = [stim_end]

    feature_name = 'peak_time'

    feature_values1 = \
        efel.getFeatureValues(
            [trace1],
            [feature_name], raise_warnings=False)

    test_data_path = os.path.join(
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

    feature_values2 = \
        efel.getFeatureValues(
            [trace2],
            [feature_name], raise_warnings=False)

    assert (
        len(feature_values1[0][feature_name]) !=
        len(feature_values2[0][feature_name]))


def test_stimstart_stimend():
    """basic: Test exception when stimstart or stimend are wrong"""

    import efel
    efel.reset()

    stim_start = 500.0
    stim_end = 900.0

    time = efel.io.load_fragment('%s#col=1' % meanfrequency1_url)
    voltage = efel.io.load_fragment('%s#col=2' % meanfrequency1_url)

    trace = {}

    trace['T'] = time
    trace['V'] = voltage
    trace['stim_start'] = stim_start
    trace['stim_end'] = stim_end

    features = ['AP_begin_voltage']

    pytest.raises(
        Exception,
        efel.getFeatureValues, [trace], features)

    trace['stim_start'] = [stim_end]
    trace['stim_end'] = [stim_start]

    pytest.raises(
        Exception,
        efel.getFeatureValues, [trace], features)

    trace['stim_start'] = [stim_start, stim_end]
    trace['stim_end'] = [stim_end]

    pytest.raises(
        Exception,
        efel.getFeatureValues, [trace], features)

    del trace['stim_start']

    pytest.raises(
        Exception,
        efel.getFeatureValues, [trace], features)


def test_setDerivativeThreshold():
    """basic: Test setDerivativeThreshold"""

    import efel
    efel.reset()

    stim_start = 500.0
    stim_end = 900.0

    time = efel.io.load_fragment('%s#col=1' % meanfrequency1_url)
    voltage = efel.io.load_fragment('%s#col=2' % meanfrequency1_url)

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
    numpy.testing.assert_allclose(AP_begin_voltage, -51.6400489995987)
    assert AP_begin_voltage != AP_begin_voltage_orig


def interpolate(time, voltage, new_dt):
    """Interpolate voltage to new dt"""

    interp_time = numpy.arange(time[0], time[-1] + new_dt, new_dt)
    interp_voltage = numpy.interp(interp_time, time, voltage)

    return interp_time, interp_voltage


def test_interpolate():
    """basic: Test interpolate"""

    import efel
    efel.reset()
    trace, time, voltage, stim_start, stim_end = load_data(
        'mean_frequency1', interp=True)

    features = ['time', 'voltage']

    feature_values = \
        efel.getFeatureValues(
            [trace],
            features, raise_warnings=False)
    interp_time = feature_values[0]['time']
    interp_voltage = feature_values[0]['voltage']
    assert len(interp_time) == len(time)
    assert len(interp_voltage) == len(voltage)
    assert numpy.allclose(interp_voltage, voltage)


def test_zero_ISI_log_slope_skip():
    """basic: Test zero ISI_log_slope_skip"""

    import efel
    efel.reset()

    stim_start = 31.2
    stim_end = 431.2

    time = efel.io.load_fragment('%s#col=1' % zeroISIlog1_url)
    voltage = efel.io.load_fragment('%s#col=2' % zeroISIlog1_url)
    trace = {}

    trace['T'] = time
    trace['V'] = voltage
    trace['stim_start'] = [stim_start]
    trace['stim_end'] = [stim_end]

    features = ['ISI_log_slope_skip']

    feature_values = \
        efel.getFeatureValues(
            [trace],
            features, raise_warnings=False)
    assert feature_values[0]['ISI_log_slope_skip'] is None


def test_peak_indices():
    """basic: Test peak_indices"""

    import efel
    efel.reset()

    stim_start = 650.0
    stim_end = 900.0

    time = efel.io.load_fragment('%s#col=1' % meanfrequency1_url)
    voltage = efel.io.load_fragment('%s#col=2' % meanfrequency1_url)
    trace = {}

    trace['T'] = time
    trace['V'] = voltage
    trace['stim_start'] = [stim_start]
    trace['stim_end'] = [stim_end]

    features = ['peak_indices']

    feature_values = \
        efel.getFeatureValues(
            [trace],
            features, raise_warnings=False)

    peak_indices = feature_values[0]['peak_indices']

    assert len(peak_indices) == 5


def test_min_AHP_indices():
    """basic: Test min_AHP_indices"""

    import efel
    efel.reset()

    stim_start = 650.0
    stim_end = 900.0

    time = efel.io.load_fragment('%s#col=1' % meanfrequency1_url)
    voltage = efel.io.load_fragment('%s#col=2' % meanfrequency1_url)
    trace = {}

    trace['T'] = time
    trace['V'] = voltage
    trace['stim_start'] = [stim_start]
    trace['stim_end'] = [stim_end]

    features = ['min_AHP_indices']

    feature_values = \
        efel.getFeatureValues(
            [trace],
            features, raise_warnings=False)

    min_AHP_indices = feature_values[0]['min_AHP_indices']

    assert len(min_AHP_indices) == 5


def test_min_AHP_indices_strict():
    """basic: Test min_AHP_indices with strict_stiminterval"""

    import efel

    for strict, n_of_ahp in [(False, 17), (True, 16)]:
        efel.reset()
        efel.setIntSetting('strict_stiminterval', strict)

        stim_start = 700.0
        stim_end = 2700.0

        time = efel.io.load_fragment('%s#col=1' % ahptest1_url)
        voltage = efel.io.load_fragment('%s#col=2' % ahptest1_url)

        trace = {}

        trace['T'] = time
        trace['V'] = voltage
        trace['stim_start'] = [stim_start]
        trace['stim_end'] = [stim_end]

        features = ['min_AHP_indices', 'AHP_time_from_peak', 'peak_time']

        feature_values = \
            efel.getFeatureValues(
                [trace],
                features, raise_warnings=False)

        min_AHP_indices = feature_values[0]['min_AHP_indices']
        AHP_time_from_peak = feature_values[0]['AHP_time_from_peak']

        assert len(min_AHP_indices) == n_of_ahp
        assert len(AHP_time_from_peak) == n_of_ahp


def test_min_AHP_indices_single_peak():
    """basic: Test min_AHP_indices with a single peak."""

    import efel

    trace_file = os.path.join(
        testdata_dir,
        'basic',
        'min_AHP_values_single_peak.txt')
    trace_values = numpy.loadtxt(trace_file)

    trace = {}
    trace["T"] = trace_values[:, 0]
    trace["V"] = trace_values[:, 1]
    trace["stim_start"] = [1950]
    trace["stim_end"] = [2050]

    feats = efel.getFeatureValues(
        [trace], ["min_AHP_values", "min_AHP_indices", "peak_indices"])

    assert len(feats[0]["peak_indices"]) == 1
    assert feats[0]["min_AHP_indices"] is None
    assert feats[0]["min_AHP_values"] is None


def test_strict_stiminterval():
    """basic: Test strict_stiminterval"""

    import efel

    for strict, n_of_spikes in [(False, 5), (True, 3)]:
        efel.reset()
        efel.setIntSetting("strict_stiminterval", strict)

        stim_start = 600.0
        stim_end = 750.0

        time = efel.io.load_fragment('%s#col=1' % meanfrequency1_url)
        voltage = efel.io.load_fragment('%s#col=2' % meanfrequency1_url)
        trace = {}

        trace['T'] = time
        trace['V'] = voltage
        trace['stim_start'] = [stim_start]
        trace['stim_end'] = [stim_end]

        features = ['peak_indices', 'peak_time', 'Spikecount']

        feature_values = \
            efel.getFeatureValues(
                [trace],
                features, raise_warnings=False)

        peak_indices = feature_values[0]['peak_indices']
        peak_time = feature_values[0]['peak_time']
        spikecount = feature_values[0]['Spikecount']

        assert len(peak_indices) == n_of_spikes
        assert len(peak_time) == n_of_spikes
        assert spikecount == n_of_spikes


def test_ISI_log_slope():
    """basic: Test ISI_log_slope"""

    import efel
    efel.reset()

    stim_start = 500.0
    stim_end = 900.0

    time = efel.io.load_fragment('%s#col=1' % meanfrequency1_url)
    voltage = efel.io.load_fragment('%s#col=2' % meanfrequency1_url)
    trace = {}

    trace['T'] = time
    trace['V'] = voltage
    trace['stim_start'] = [stim_start]
    trace['stim_end'] = [stim_end]

    features = ['ISI_values', 'ISI_log_slope']

    feature_values = \
        efel.getFeatureValues(
            [trace],
            features, raise_warnings=False)
    isi_values = feature_values[0]['ISI_values']
    x_values = numpy.arange(0, len(isi_values)) + 1.0

    # fit
    log_x_values = numpy.log(x_values)
    log_isi_values = numpy.log(isi_values)
    slope, _ = numpy.polyfit(log_x_values, log_isi_values, 1)

    numpy.testing.assert_allclose(feature_values[0]['ISI_log_slope'][0], slope)


def test_ISI_semilog_slope():
    """basic: Test ISI_semilog_slope"""

    import efel
    efel.reset()

    stim_start = 500.0
    stim_end = 900.0

    time = efel.io.load_fragment('%s#col=1' % meanfrequency1_url)
    voltage = efel.io.load_fragment('%s#col=2' % meanfrequency1_url)
    trace = {}

    trace['T'] = time
    trace['V'] = voltage
    trace['stim_start'] = [stim_start]
    trace['stim_end'] = [stim_end]

    features = ['ISI_values', 'ISI_semilog_slope']

    feature_values = \
        efel.getFeatureValues(
            [trace],
            features, raise_warnings=False)
    isi_values = feature_values[0]['ISI_values']
    x_values = numpy.arange(0, len(isi_values)) + 1.0

    # fit
    x_values = x_values
    log_isi_values = numpy.log(isi_values)
    slope, _ = numpy.polyfit(x_values, log_isi_values, 1)

    numpy.testing.assert_allclose(
        feature_values[0]['ISI_semilog_slope'][0], slope
    )


def test_AP_begin_indices1():
    """basic: Test AP_begin_indices 1"""

    import efel
    efel.reset()

    stim_start = 31.2
    stim_end = 431.2

    test_data_path = os.path.join(
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

    feature_values = \
        efel.getFeatureValues(
            [trace],
            features,
            raise_warnings=False)
    # Make sure the amount of peak_times, AP_begin_indices and AP_amplitude is
    # the same in this trace
    # There was originally an issue in this case due to the 'width' value
    # in AP_begin_indices, which caused a segmentation fault
    assert (
        len(feature_values[0]['AP_begin_indices']) ==
        len(feature_values[0]['AP_amplitude']))
    assert (
        len(feature_values[0]['AP_begin_indices']) ==
        len(feature_values[0]['peak_time']))
    assert (
        len(feature_values[0]['AP_begin_indices']) ==
        len(feature_values[0]['AP_duration_half_width']))


def test_AP_end_indices():
    """basic: Test AP end indices."""
    import efel
    efel.reset()

    stim_start = 31.2
    stim_end = 431.2

    test_data_path = os.path.join(
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
        'peak_indices',
        'AP_end_indices']

    feature_values = \
        efel.getFeatureValues(
            [trace],
            features,
            raise_warnings=False)

    begin_indices = feature_values[0]["AP_begin_indices"]
    peak_indices = feature_values[0]["peak_indices"]
    end_indices = feature_values[0]["AP_end_indices"]

    for begin, peak, end in zip(begin_indices, peak_indices, end_indices):
        # the voltage value for the end index should be closer than that of
        # begin index than the peak
        assert (abs(voltage[begin] - voltage[end])
                < abs(voltage[peak] - voltage[end]))
        assert end > begin

    efel.reset()

    efel.setDoubleSetting("DownDerivativeThreshold", -24)
    feature_values = \
        efel.getFeatureValues(
            [trace],
            features,
            raise_warnings=False)

    updated_end_indices = feature_values[0]["AP_end_indices"]

    for end_index, updated_end_index in zip(end_indices, updated_end_indices):
        assert end_index != updated_end_index


def test_mean_frequency1():
    """basic: Test mean_frequency 1"""

    import efel
    efel.reset()

    trace, time, voltage, stim_start, stim_end = load_data(
        'mean_frequency1', interp=True)

    features = ['mean_frequency', 'peak_time']

    feature_values = \
        efel.getFeatureValues(
            [trace],
            features)

    peak_times = feature_values[0]['peak_time']

    stim_spikes = peak_times[numpy.where((stim_start <= peak_times)
                                         & (peak_times <= stim_end))]
    n_of_spikes = len(stim_spikes)

    mean_frequency = float(n_of_spikes) * 1000 / \
        (stim_spikes[-1] - stim_start)

    numpy.testing.assert_allclose(
        feature_values[0]['mean_frequency'][0],
        mean_frequency)


def test_ap_amplitude_outside_stim():
    """basic: Test AP amplitude with spike outside stim"""

    import efel
    efel.reset()

    stim_start = 700.0
    stim_end = 2700.0

    test_data_path = os.path.join(
        testdata_dir,
        'basic',
        'spike_outside_stim.txt')
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
            features,
            raise_warnings=False)

    # Make sure AP_amplitude doesn't pick up the two spikes outside of
    # the stimulus
    # (which are present in peak_time)
    assert (
        len(feature_values[0]['AP_amplitude']) + 2 ==
        len(feature_values[0]['peak_time']))


def test_ap_amplitude_from_voltagebase1():
    """basic: Test AP_amplitude_from_voltagebase 1"""

    import efel
    efel.reset()

    stim_start = 500.0
    stim_end = 900.0

    time = efel.io.load_fragment('%s#col=1' % meanfrequency1_url)
    voltage = efel.io.load_fragment('%s#col=2' % meanfrequency1_url)

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
        numpy.testing.assert_allclose(peak_voltage - voltage_base,
                                      ap_amplitude_from_voltagebase)


def test_voltagebase1():
    """basic: Test voltagebase 1"""

    import efel
    efel.reset()

    trace, time, voltage, stim_start, stim_end = load_data(
        'mean_frequency1', interp=True)

    features = ['voltage_base']

    feature_values = \
        efel.getFeatureValues(
            [trace],
            features)

    interp_time, interp_voltage = interpolate(time, voltage, 0.1)

    voltage_base = numpy.mean(interp_voltage[numpy.where(
        (interp_time >= 0.9 * stim_start) & (interp_time <= stim_start))])

    numpy.testing.assert_allclose(voltage_base,
                                  feature_values[0]['voltage_base'][0],
                                  rtol=0, atol=1e-8)


def test_voltagebase_median():
    """basic: Test voltagebase computation with median option"""

    import efel
    efel.reset()
    efel.setStrSetting("voltage_base_mode", "median")

    trace, time, voltage, stim_start, stim_end = load_data(
        'mean_frequency1', interp=True)

    features = ['voltage_base']

    feature_values = \
        efel.getFeatureValues(
            [trace],
            features)

    interp_time, interp_voltage = interpolate(time, voltage, 0.1)

    voltage_base = numpy.median(interp_voltage[numpy.where(
        (interp_time >= 0.9 * stim_start) & (interp_time <= stim_start))])

    numpy.testing.assert_allclose(voltage_base,
                                  feature_values[0]['voltage_base'][0],
                                  rtol=0, atol=1e-8)


def test_currentbase():
    """basic: Test currentbase"""

    import efel
    efel.reset()

    data = numpy.loadtxt(os.path.join(os.path.abspath(testdata_dir),
                                      'basic',
                                      'current.txt'))
    current = data[:, 1]
    time = data[:, 0]
    stim_start = 2.0
    stim_end = 900.0  # not to be used

    trace = {'T': time, 'I': current,
             'stim_start': [stim_start], 'stim_end': [stim_end]}

    feature_values = efel.getFeatureValues([trace], ['current_base'])

    current_base = numpy.mean(current[numpy.where(
        (time >= 0.9 * stim_start) & (time <= stim_start))])

    # nt.set_trace()
    numpy.testing.assert_allclose(current_base,
                                  feature_values[0]['current_base'][0],
                                  rtol=0, atol=1e-8)


def test_currentbase_median():
    """basic: Test currentbase with median"""

    import efel
    efel.reset()
    efel.setStrSetting("current_base_mode", "median")

    data = numpy.loadtxt(os.path.join(os.path.abspath(testdata_dir),
                                      'basic',
                                      'current.txt'))
    current = data[:, 1]
    time = data[:, 0]
    stim_start = 2.0
    stim_end = 900.0  # not to be used

    trace = {'T': time, 'I': current,
             'stim_start': [stim_start], 'stim_end': [stim_end]}

    feature_values = efel.getFeatureValues([trace], ['current_base'])

    current_base = numpy.median(current[numpy.where(
        (time >= 0.9 * stim_start) & (time <= stim_start))])

    numpy.testing.assert_allclose(current_base,
                                  feature_values[0]['current_base'][0],
                                  rtol=0, atol=1e-8)


def test_getDistance1():
    """basic: Test getDistance 1"""

    import efel
    efel.reset()

    stim_start = 500.0
    stim_end = 900.0

    time = efel.io.load_fragment('%s#col=1' % meanfrequency1_url)
    voltage = efel.io.load_fragment('%s#col=2' % meanfrequency1_url)

    trace = {}

    trace['T'] = time
    trace['V'] = voltage
    trace['stim_start'] = [stim_start]
    trace['stim_end'] = [stim_end]

    numpy.testing.assert_allclose(
        3.09045815935,
        efel.getDistance(
            trace,
            'AP_amplitude',
            50,
            10))


def test_getDistance_error_dist():
    """basic: Test getDistance error_dist option"""

    import efel
    efel.reset()

    stim_start = 400.0
    stim_end = 500.0

    time = efel.io.load_fragment('%s#col=1' % meanfrequency1_url)
    voltage = efel.io.load_fragment('%s#col=2' % meanfrequency1_url)

    trace = {}

    trace['T'] = time
    trace['V'] = voltage
    trace['stim_start'] = [stim_start]
    trace['stim_end'] = [stim_end]

    score_normal = efel.getDistance(
        trace,
        'AP_amplitude',
        50,
        10)
    score_150 = efel.getDistance(
        trace,
        'AP_amplitude',
        50,
        10,
        error_dist=150)

    numpy.testing.assert_allclose(score_normal, 250)
    numpy.testing.assert_allclose(score_150, 150)


def test_getDistance_trace_check():
    """basic: Test getDistance trace_check option"""

    import efel
    efel.reset()

    dt = 0.1

    # voltage trace at constant -70 mV
    v = numpy.zeros(int(100 / dt)) - 70.0

    # create 'spikes' at 20, 40 and 60 ms
    v[int(20 / dt):int(25 / dt)] = 20.
    v[int(40 / dt):int(45 / dt)] = 20.
    v[int(60 / dt):int(65 / dt)] = 20.

    traces = []
    trace = {}
    trace['T'] = numpy.arange(len(v)) * dt
    trace['V'] = v
    trace['stim_start'] = [10]
    trace['stim_end'] = [70]
    traces.append(trace)

    numpy.testing.assert_allclose(
        efel.getDistance(trace, 'Spikecount', 0, 1), 3.0
    )

    trace['stim_end'] = [50]

    efel.reset()
    numpy.testing.assert_allclose(
        efel.getDistance(
            trace,
            'Spikecount',
            0,
            1,
            trace_check=False),
        3.0)

    efel.reset()
    numpy.testing.assert_allclose(
        efel.getDistance(trace, 'Spikecount', 0, 1), 250.0
    )


def test_APlast_amp():
    """basic: Test APlast_amp"""

    import efel
    efel.reset()

    stim_start = 500.0
    stim_end = 900.0

    time = efel.io.load_fragment('%s#col=1' % meanfrequency1_url)
    voltage = efel.io.load_fragment('%s#col=2' % meanfrequency1_url)

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
    assert APlast_amp == AP_amplitude[-1]


def test_APlast_width():
    """basic: Test APlast_width"""

    import efel
    efel.reset()

    stim_start = 500.0
    stim_end = 900.0

    time = efel.io.load_fragment('%s#col=1' % meanfrequency1_url)
    voltage = efel.io.load_fragment('%s#col=2' % meanfrequency1_url)

    trace = {}

    trace['T'] = time
    trace['V'] = voltage
    trace['stim_start'] = [stim_start]
    trace['stim_end'] = [stim_end]

    features = ['spike_half_width', 'APlast_width']

    feature_values = \
        efel.getFeatureValues(
            [trace],
            features)

    APlast_width = feature_values[0]['APlast_width'][0]
    spike_half_width = feature_values[0]['spike_half_width']
    assert APlast_width == spike_half_width[-1]


def test_derivwindow1():
    """basic: Test DerivativeWindow"""

    import efel
    efel.reset()

    stim_start = 100.0
    stim_end = 1000.0

    time = efel.io.load_fragment('%s#col=1' % derivwindow1_url)
    voltage = efel.io.load_fragment('%s#col=2' % derivwindow1_url)

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

    AP_begin_voltage = feature_values[0]['AP_begin_voltage'][0]
    numpy.testing.assert_allclose(AP_begin_voltage, -45.03627393790836)

    efel.reset()
    efel.setDoubleSetting('interp_step', 0.01)
    feature_values = \
        efel.getFeatureValues(
            [trace],
            features)

    AP_begin_voltage = feature_values[0]['AP_begin_voltage'][0]
    numpy.testing.assert_allclose(AP_begin_voltage, -83.57661997973835)

    efel.reset()
    efel.setDoubleSetting('interp_step', 0.01)
    efel.setIntSetting('DerivativeWindow', 30)
    feature_values = \
        efel.getFeatureValues(
            [trace],
            features)

    AP_begin_voltage = feature_values[0]['AP_begin_voltage'][0]
    numpy.testing.assert_allclose(AP_begin_voltage, -45.505521563640386)


def test_spikecount1():
    """basic: Test Spikecount 1"""

    import efel
    efel.reset()

    stim_start = 500.0
    stim_end = 900.0

    time = efel.io.load_fragment('%s#col=1' % meanfrequency1_url)
    voltage = efel.io.load_fragment('%s#col=2' % meanfrequency1_url)

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
    assert len(peak_indices) == spikecount


def test_spikecount_stimint1():
    """basic: Test Spikecount_stimint 1"""

    import efel
    efel.reset()

    stim_start = 700.0
    stim_end = 2700.0

    time = efel.io.load_fragment('%s#col=1' % spikeoutsidestim_url)
    voltage = efel.io.load_fragment('%s#col=2' % spikeoutsidestim_url)

    trace = {}

    trace['T'] = time
    trace['V'] = voltage
    trace['stim_start'] = [stim_start]
    trace['stim_end'] = [stim_end]

    features = ['peak_time', 'Spikecount_stimint', 'Spikecount']

    feature_values = \
        efel.getFeatureValues(
            [trace],
            features)

    peak_times = feature_values[0]['peak_time']
    spikecount = feature_values[0]['Spikecount'][0]
    spikecount_stimint = feature_values[0]['Spikecount_stimint'][0]

    interval_peaktimes, = \
        numpy.where((peak_times >= stim_start) & (peak_times <= stim_end))

    assert (
        len(interval_peaktimes) ==
        spikecount_stimint)

    assert (
        spikecount ==
        spikecount_stimint + 2)


def test_spikecount_libv4peakindices():
    """basic: Test Spikecount in combination with LibV4 peak_indices"""

    import efel
    efel.reset()

    stim_start = 500.0
    stim_end = 900.0

    time = efel.io.load_fragment('%s#col=1' % meanfrequency1_url)
    voltage = efel.io.load_fragment('%s#col=2' % meanfrequency1_url)

    trace = {}

    trace['T'] = time
    trace['V'] = voltage
    trace['stim_start'] = [stim_start]
    trace['stim_end'] = [stim_end]

    features = ['peak_indices', 'Spikecount']

    test_peak = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             'DependencyV5_LibV4peakindices.txt')
    efel.setDependencyFileLocation(test_peak)

    feature_values = \
        efel.getFeatureValues(
            [trace],
            features)

    peak_indices = feature_values[0]['peak_indices']
    spikecount = feature_values[0]['Spikecount'][0]
    assert len(peak_indices) == 5
    assert len(peak_indices) == spikecount


def test_ohmic_inputresistance():
    """basic: Test ohmic_input_resistance"""

    import efel
    efel.reset()

    stim_start = 500.0
    stim_end = 900.0

    time = efel.io.load_fragment('%s#col=1' % meanfrequency1_url)
    voltage = efel.io.load_fragment('%s#col=2' % meanfrequency1_url)

    trace = {}

    trace['T'] = time
    trace['V'] = voltage
    trace['stim_start'] = [stim_start]
    trace['stim_end'] = [stim_end]

    features = ['ohmic_input_resistance', 'voltage_deflection']

    stimulus_current = 10.0
    efel.setDoubleSetting('stimulus_current', stimulus_current)
    feature_values = \
        efel.getFeatureValues(
            [trace],
            features)

    voltage_deflection = feature_values[0]['voltage_deflection'][0]
    ohmic_input_resistance = feature_values[0]['ohmic_input_resistance'][0]
    assert (
        ohmic_input_resistance ==
        voltage_deflection /
        stimulus_current)


def test_sag_amplitude():
    """basic: Test sag_amplitude"""

    import efel
    efel.reset()

    stim_start = 800.0
    stim_end = 3800.0

    time = efel.io.load_fragment('%s#col=1' % sagtrace1_url)
    voltage = efel.io.load_fragment('%s#col=2' % sagtrace1_url)

    trace = {}

    trace['T'] = time
    trace['V'] = voltage
    trace['stim_start'] = [stim_start]
    trace['stim_end'] = [stim_end]

    features = [
        'sag_amplitude',
        'steady_state_voltage_stimend',
        'minimum_voltage']

    feature_values = \
        efel.getFeatureValues(
            [trace],
            features)

    steady_state_voltage_stimend = feature_values[
        0]['steady_state_voltage_stimend'][0]
    minimum_voltage = feature_values[0]['minimum_voltage'][0]
    sag_amplitude = feature_values[0]['sag_amplitude'][0]
    assert (
        sag_amplitude ==
        steady_state_voltage_stimend - minimum_voltage)


def test_sag_amplitude_pos_deflect():
    """basic: Test if sag_amplitude throws error on positive deflection"""

    import efel
    efel.reset()

    stim_start = 500.0
    stim_end = 900.0

    time = efel.io.load_fragment('%s#col=1' % meanfrequency1_url)
    voltage = efel.io.load_fragment('%s#col=2' % meanfrequency1_url)

    trace = {}

    trace['T'] = time
    trace['V'] = voltage
    trace['stim_start'] = [stim_start]
    trace['stim_end'] = [stim_end]

    features = ['sag_amplitude']

    feature_values = \
        efel.getFeatureValues(
            [trace],
            features,
            raise_warnings=False)

    assert (feature_values[0]['sag_amplitude'] is None)


def test_sag_ratio1():
    """basic: Test sag_ratio1"""

    import efel
    efel.reset()

    stim_start = 800.0
    stim_end = 3800.0

    time = efel.io.load_fragment('%s#col=1' % sagtrace1_url)
    voltage = efel.io.load_fragment('%s#col=2' % sagtrace1_url)

    trace = {}

    trace['T'] = time
    trace['V'] = voltage
    trace['stim_start'] = [stim_start]
    trace['stim_end'] = [stim_end]

    features = [
        'sag_ratio1',
        'sag_amplitude',
        'minimum_voltage',
        'voltage_base']

    feature_values = \
        efel.getFeatureValues(
            [trace],
            features)

    sag_amplitude = feature_values[0]['sag_amplitude'][0]
    minimum_voltage = feature_values[0]['minimum_voltage'][0]
    voltage_base = feature_values[0]['voltage_base'][0]
    sag_ratio1 = feature_values[0]['sag_ratio1'][0]
    assert (
        sag_ratio1 ==
        sag_amplitude / (voltage_base - minimum_voltage))


def test_sag_ratio1_empty():
    """basic: Test sag_ratio1 on empty trace"""

    import efel
    efel.reset()

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

    features = ['sag_ratio1']

    feature_values = \
        efel.getFeatureValues(
            [trace],
            features,
            raise_warnings=False)

    assert feature_values[0]['sag_ratio1'] is None


def test_sag_ratio2():
    """basic: Test sag_ratio2"""

    import efel
    efel.reset()

    stim_start = 800.0
    stim_end = 3800.0

    time = efel.io.load_fragment('%s#col=1' % sagtrace1_url)
    voltage = efel.io.load_fragment('%s#col=2' % sagtrace1_url)

    trace = {}

    trace['T'] = time
    trace['V'] = voltage
    trace['stim_start'] = [stim_start]
    trace['stim_end'] = [stim_end]

    features = [
        'sag_ratio2',
        'minimum_voltage',
        'steady_state_voltage_stimend',
        'voltage_base']

    feature_values = \
        efel.getFeatureValues(
            [trace],
            features)

    steady_state_voltage_stimend = \
        feature_values[0]['steady_state_voltage_stimend'][0]
    minimum_voltage = feature_values[0]['minimum_voltage'][0]
    voltage_base = feature_values[0]['voltage_base'][0]
    sag_ratio2 = feature_values[0]['sag_ratio2'][0]
    assert (
        sag_ratio2 ==
        (voltage_base - steady_state_voltage_stimend) /
        (voltage_base - minimum_voltage))


def test_ohmic_input_resistance_vb_ssse():
    """basic: Test ohmic_input_resistance_vb_ssse"""

    import efel
    efel.reset()

    stim_start = 500.0
    stim_end = 900.0

    time = efel.io.load_fragment('%s#col=1' % meanfrequency1_url)
    voltage = efel.io.load_fragment('%s#col=2' % meanfrequency1_url)

    trace = {}

    trace['T'] = time
    trace['V'] = voltage
    trace['stim_start'] = [stim_start]
    trace['stim_end'] = [stim_end]

    features = ['ohmic_input_resistance_vb_ssse', 'voltage_deflection_vb_ssse']

    stimulus_current = 10.0
    efel.setDoubleSetting('stimulus_current', stimulus_current)
    feature_values = \
        efel.getFeatureValues(
            [trace],
            features)

    voltage_deflection = feature_values[0]['voltage_deflection_vb_ssse'][0]
    ohmic_input_resistance = \
        feature_values[0]['ohmic_input_resistance_vb_ssse'][0]
    assert (
        ohmic_input_resistance ==
        voltage_deflection /
        stimulus_current)


def test_spikecount2():
    """basic: Test Spikecount 2: test empty trace"""

    import efel
    efel.reset()

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
    assert spikecount == 0


def test_min_voltage_between_spikes1():
    """basic: Test min_voltage_between_spikes 1"""

    import efel
    efel.reset()

    stim_start = 500.0
    stim_end = 900.0

    time = efel.io.load_fragment('%s#col=1' % meanfrequency1_url)
    voltage = efel.io.load_fragment('%s#col=2' % meanfrequency1_url)

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
        numpy.testing.assert_allclose(
            numpy.min(
                fel_voltage[
                    peak_indices[index]:peak_indices[
                        index +
                        1]]),
            min_voltage_between_spikes_value)


def test_getFeatureNames():
    """basic: Test getting all feature names"""

    import efel
    efel.reset()
    import json

    test_data_path = os.path.join(testdata_dir, '..', 'featurenames.json')
    with open(test_data_path, 'r') as featurenames_json:
        expected_featurenames = json.load(featurenames_json)
    assert set(efel.getFeatureNames()) == set(expected_featurenames)


def test_getFeatureNameExists():
    """basic: Test FeatureNameExists"""

    import efel
    efel.reset()
    assert efel.FeatureNameExists('voltage_base')
    assert not efel.FeatureNameExists('voltage_base_wrong')


def test_steady_state_voltage1():
    """basic: Test steady_state_voltage"""

    import efel
    efel.reset()

    trace, time, voltage, stim_start, stim_end = load_data(
        'mean_frequency1', interp=True)

    features = ['steady_state_voltage']

    feature_values = \
        efel.getFeatureValues(
            [trace],
            features)[0]

    steady_state_voltage = numpy.mean(voltage[numpy.where(time >= stim_end)])

    numpy.testing.assert_allclose(steady_state_voltage,
                                  feature_values['steady_state_voltage'][0])


def test_steady_state_voltage_stimend():
    """basic: Test steady_state_voltage_stimend"""

    import efel
    efel.reset()

    trace, time, voltage, stim_start, stim_end = load_data(
        'mean_frequency1', interp=True)

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

    numpy.testing.assert_allclose(
        steady_state_voltage_stimend,
        feature_values['steady_state_voltage_stimend'][0]
    )


def test_maximum_voltage_from_voltagebase():
    """basic: Test maximum_voltage_from_voltagebase"""

    import efel
    efel.reset()

    trace, time, voltage, stim_start, stim_end = load_data(
        'mean_frequency1', interp=True)

    features = ['maximum_voltage_from_voltagebase', 'voltage_base']

    feature_values = \
        efel.getFeatureValues(
            [trace],
            features)[0]

    maximum_voltage = numpy.max(voltage[numpy.where(
        (time <= stim_end) & (time >= stim_start)
    )])

    voltage_base = feature_values['voltage_base'][0]

    maximum_voltage_from_voltagebase = maximum_voltage - voltage_base

    numpy.testing.assert_allclose(
        maximum_voltage_from_voltagebase,
        feature_values['maximum_voltage_from_voltagebase'][0])


def decay_time_constant_after_stim(time, voltage, interval_start,
                                   interval_end, stim_start, stim_end):
    """decay_time_constant_after_stim numpy implementation"""

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
    """basic: Test decay_time_constant_after_stim 1"""

    import efel
    efel.reset()

    trace, time, voltage, stim_start, stim_end = load_data(
        'mean_frequency1', interp=True)

    features = ['decay_time_constant_after_stim']

    feature_values = efel.getFeatureValues([trace], features)[0]

    expected = decay_time_constant_after_stim(
        time,
        voltage,
        stim_end + 1.0,
        stim_end + 10.0,
        stim_start,
        stim_end)

    numpy.testing.assert_allclose(
        expected,
        feature_values['decay_time_constant_after_stim'][0])


def test_decay_time_constant_after_stim2():
    """basic: Test decay_time_constant_after_stim 2"""

    import efel
    efel.reset()

    trace, time, voltage, stim_start, stim_end = load_data(
        'tau20.0', interp=True)

    features = ['decay_time_constant_after_stim']

    feature_values = efel.getFeatureValues([trace], features)[0]

    numpy.testing.assert_allclose(
        19.9,
        feature_values['decay_time_constant_after_stim'][0], rtol=0, atol=1e-1)


def sag_time_constant(
        time, voltage, min_v, steady_state_v, sag_ampl, stim_start, stim_end):
    """sag_time_constant numpy implementation"""
    # select t, v in stimulus interval
    start_idx = numpy.where(time == stim_start)[0][0]
    end_idx = numpy.where(time == stim_end)[0][0]
    vinterval = voltage[start_idx:end_idx]
    tinterval = time[start_idx:end_idx]

    # get start decay
    start_decay = numpy.argmin(vinterval)

    # get end decay
    v90 = steady_state_v - 0.1 * sag_ampl
    end_decay = numpy.where(
        (tinterval > tinterval[start_decay]) & (vinterval >= v90)
    )[0][0]

    v_reference = vinterval[end_decay]

    # select t, v in decay interval
    interval_indices = numpy.arange(start_decay, end_decay)
    interval_time = tinterval[interval_indices]
    interval_voltage = abs(vinterval[interval_indices] - v_reference)

    # get tau
    log_interval_voltage = numpy.log(interval_voltage)
    slope, _ = numpy.polyfit(interval_time, log_interval_voltage, 1)
    tau = abs(1. / slope)

    return tau


def test_sag_time_constant():
    """basic: Test sag_time_constant"""

    import efel
    efel.reset()

    interp_dt = 0.1

    stim_start = 800.0
    stim_end = 3800.0
    time = efel.io.load_fragment('%s#col=1' % sagtrace1_url)
    voltage = efel.io.load_fragment('%s#col=2' % sagtrace1_url)
    time, voltage = interpolate(time, voltage, interp_dt)

    trace = {}
    trace['T'] = time
    trace['V'] = voltage
    trace['stim_start'] = [stim_start]
    trace['stim_end'] = [stim_end]

    features = [
        'minimum_voltage',
        'steady_state_voltage_stimend',
        'sag_time_constant',
        'sag_amplitude'
    ]
    feature_values = efel.getFeatureValues([trace], features)[0]

    min_v = feature_values['minimum_voltage'][0]
    steady_state_v = feature_values['steady_state_voltage_stimend'][0]
    sag_ampl = feature_values['sag_amplitude'][0]

    expected = sag_time_constant(
        time,
        voltage,
        min_v,
        steady_state_v,
        sag_ampl,
        stim_start,
        stim_end)

    numpy.testing.assert_allclose(
        expected,
        feature_values['sag_time_constant'][0])


def test_getmeanfeaturevalues():
    """basic: Test getMeanFeatureValues"""

    import efel
    efel.reset()

    stim_start = 500.0
    stim_end = 900.0

    time = efel.io.load_fragment('%s#col=1' % meanfrequency1_url)
    voltage = efel.io.load_fragment('%s#col=2' % meanfrequency1_url)

    trace = {}

    trace['T'] = time
    trace['V'] = voltage
    trace['stim_start'] = [stim_start]
    trace['stim_end'] = [stim_end]

    feature_values = \
        efel.getFeatureValues(
            [trace],
            ['AP_amplitude', 'BPAPHeightLoc1'], raise_warnings=False)

    mean_feature_values = efel.getMeanFeatureValues(
        [trace], [
            'AP_amplitude', 'BPAPHeightLoc1'], raise_warnings=False)

    assert (numpy.mean(feature_values[0]['AP_amplitude']) ==
            mean_feature_values[0]['AP_amplitude'])


def test_mean_AP_amplitude():
    """basic: Test mean_AP_amplitude"""

    import efel
    efel.reset()

    stim_start = 500.0
    stim_end = 900.0

    time = efel.io.load_fragment('%s#col=1' % meanfrequency1_url)
    voltage = efel.io.load_fragment('%s#col=2' % meanfrequency1_url)

    trace = {}

    trace['T'] = time
    trace['V'] = voltage
    trace['stim_start'] = [stim_start]
    trace['stim_end'] = [stim_end]

    feature_values = \
        efel.getFeatureValues(
            [trace],
            ['AP_amplitude', 'mean_AP_amplitude'], raise_warnings=False)

    assert (numpy.mean(feature_values[0]['AP_amplitude']) ==
            feature_values[0]['mean_AP_amplitude'])


def test_unfinished_peak():
    """basic: Test if unfinished peak doesn't break Spikecount"""

    import efel
    efel.setIntSetting('strict_stiminterval', True)

    dt = 0.1
    v = numpy.zeros(int(100 / dt)) - 70.0
    v[int(20 / dt):int(25 / dt)] = 20.
    v[int(40 / dt):int(45 / dt)] = 20.
    v[int(60 / dt):int(65 / dt)] = 20.

    trace = {}
    trace['T'] = numpy.arange(len(v)) * dt
    trace['V'] = v
    trace['stim_start'] = [10]
    trace['stim_end'] = [70]

    traces_results = efel.getFeatureValues([trace], ['Spikecount'])
    spikecount = traces_results[0]['Spikecount'][0]

    assert spikecount == 3

    # When the signal at the end of the trace is larger than the threshold,
    # Spikecount and possibly other features cannont be estimated.
    v[int(80 / dt):] = -19

    traces_results = efel.getFeatureValues([trace], ['Spikecount'])
    spikecount = traces_results[0]['Spikecount'][0]

    assert spikecount == 3


def rise_time_perc(
    time, voltage,
    AP_begin_indices,
    peak_indices,
    rise_start_perc,
    rise_end_perc
):
    """AP_rise_time numpy implementation with percentages"""
    rise_times = []
    AP_amp = voltage[peak_indices] - voltage[AP_begin_indices]
    begin_voltages = AP_amp * rise_start_perc + voltage[AP_begin_indices]
    end_voltages = AP_amp * rise_end_perc + voltage[AP_begin_indices]

    for AP_begin_indice, peak_indice, begin_v, end_v in zip(
        AP_begin_indices, peak_indices, begin_voltages, end_voltages
    ):
        voltage_window = voltage[AP_begin_indice:peak_indice]

        new_begin_indice = AP_begin_indice + numpy.min(
            numpy.where(voltage_window >= begin_v)[0]
        )
        new_end_indice = AP_begin_indice + numpy.max(
            numpy.where(voltage_window <= end_v)[0]
        )

        rise_times.append(time[new_end_indice] - time[new_begin_indice])

    return numpy.array(rise_times)


def test_rise_time_perc():
    """basic: Test AP rise time percentage"""

    import efel
    efel.reset()
    trace, time, voltage, stim_start, stim_end = load_data(
        'mean_frequency1', interp=True
    )

    trace['rise_start_perc'] = [0.2]
    trace['rise_end_perc'] = [0.8]

    features = ['AP_rise_time', 'AP_begin_indices', 'peak_indices']

    feature_values = efel.getFeatureValues(
        [trace], features, raise_warnings=False
    )
    ap_rise_time = feature_values[0]['AP_rise_time']
    AP_begin_indices = feature_values[0]['AP_begin_indices']
    peak_indices = feature_values[0]['peak_indices']

    expected = rise_time_perc(
        time, voltage, AP_begin_indices, peak_indices, 0.2, 0.8
    )

    for exp, rise_time in zip(expected, ap_rise_time):
        numpy.testing.assert_allclose(exp, rise_time)


def test_slow_ahp_start():
    """basic: Test AHP_depth_abs_slow with a custom after spike start time"""

    import efel
    efel.reset()
    trace, time, voltage, stim_start, stim_end = load_data(
        'mean_frequency1', interp=True
    )

    trace['sahp_start'] = [12.0]

    features = ['AHP_depth_abs_slow', 'peak_indices']

    feature_values = efel.getFeatureValues(
        [trace], features, raise_warnings=False
    )
    peak_indices = feature_values[0]['peak_indices']
    ahp_depth_abs_slow = feature_values[0]['AHP_depth_abs_slow']

    expected = []
    for i in range(1, len(peak_indices) - 1):
        new_start_time = time[peak_indices[i]] + trace['sahp_start'][0]
        new_idx = numpy.min(numpy.where(time >= new_start_time)[0])
        expected.append(numpy.min(voltage[new_idx:peak_indices[i + 1]]))

    for exp, ahp_slow in zip(expected, ahp_depth_abs_slow):
        numpy.testing.assert_allclose(exp, ahp_slow)


def test_AP_peak_upstroke():
    """basic: Test AP_peak_upstroke (maximum peak rise rate)"""

    import efel
    efel.reset()
    trace, time, voltage, stim_start, stim_end = load_data(
        'mean_frequency1', interp=True
    )

    features = ['AP_peak_upstroke', 'peak_indices', 'AP_begin_indices']

    feature_values = efel.getFeatureValues(
        [trace], features, raise_warnings=False
    )
    peak_indices = feature_values[0]['peak_indices']
    ap_begin_indices = feature_values[0]['AP_begin_indices']
    ap_peak_upstroke = feature_values[0]['AP_peak_upstroke']

    expected = []
    # compute dv/dt  omit dx and /2 that cancel out in division
    dv = (
        [voltage[1] - voltage[0]]
        + list(voltage[2:] - voltage[:-2])
        + [voltage[-1] - voltage[-2]]
    )
    dt = (
        [time[1] - time[0]]
        + list(time[2:] - time[:-2])
        + [time[-1] - time[-2]]
    )
    dvdt = numpy.array(dv) / numpy.array(dt)
    # compute ap peak upstroke
    for apbi, pi in zip(ap_begin_indices, peak_indices):
        expected.append(numpy.max(dvdt[apbi:pi]))

    for exp, pus in zip(expected, ap_peak_upstroke):
        numpy.testing.assert_allclose(exp, pus, rtol=0, atol=1e-6)


def test_AP_peak_downstroke():
    """basic: Test AP_peak_downstroke (minimum peak fall rate)"""

    import efel
    efel.reset()
    trace, time, voltage, stim_start, stim_end = load_data(
        'mean_frequency1', interp=True
    )

    features = ['AP_peak_downstroke', 'peak_indices', 'min_AHP_indices']

    feature_values = efel.getFeatureValues(
        [trace], features, raise_warnings=False
    )
    peak_indices = feature_values[0]['peak_indices']
    min_ahp_indices = feature_values[0]['min_AHP_indices']
    ap_peak_downstroke = feature_values[0]['AP_peak_downstroke']

    expected = []
    # compute dv/dt  omit dx and /2 that cancel out in division
    dv = (
        [voltage[1] - voltage[0]]
        + list(voltage[2:] - voltage[:-2])
        + [voltage[-1] - voltage[-2]]
    )
    dt = (
        [time[1] - time[0]]
        + list(time[2:] - time[:-2])
        + [time[-1] - time[-2]]
    )
    dvdt = numpy.array(dv) / numpy.array(dt)
    # compute ap peak downstroke
    for ahpi, pi in zip(min_ahp_indices, peak_indices):
        expected.append(numpy.min(dvdt[pi:ahpi]))

    for exp, pds in zip(expected, ap_peak_downstroke):
        numpy.testing.assert_allclose(exp, pds, rtol=0, atol=1e-6)


def test_min_between_peaks_indices():
    """basic: Test min_between_peaks_indices"""

    import efel
    efel.reset()

    stim_start = 200.0
    stim_end = 1200.0

    time = efel.io.load_fragment('%s#col=1' % dendriticAP_url)
    voltage = efel.io.load_fragment('%s#col=2' % dendriticAP_url)
    trace = {}

    trace['T'] = time
    trace['V'] = voltage
    trace['stim_start'] = [stim_start]
    trace['stim_end'] = [stim_end]

    features = ['min_AHP_indices', 'min_between_peaks_indices']

    feature_values = \
        efel.getFeatureValues(
            [trace],
            features, raise_warnings=False)

    min_AHP_indices = feature_values[0]['min_AHP_indices'][0]
    min_btw_peaks_indices = feature_values[0]['min_between_peaks_indices'][0]

    assert min_AHP_indices < min_btw_peaks_indices


def test_min_between_peaks_values():
    """basic: Test min_between_peaks_values"""

    import efel
    efel.reset()

    stim_start = 200.0
    stim_end = 1200.0

    time = efel.io.load_fragment('%s#col=1' % dendriticAP_url)
    voltage = efel.io.load_fragment('%s#col=2' % dendriticAP_url)
    time, voltage = interpolate(time, voltage, 0.1)
    trace = {}

    trace['T'] = time
    trace['V'] = voltage
    trace['stim_start'] = [stim_start]
    trace['stim_end'] = [stim_end]

    features = ['min_between_peaks_values', 'peak_indices']

    feature_values = \
        efel.getFeatureValues(
            [trace],
            features, raise_warnings=False)

    min_btw_peaks_value = feature_values[0]['min_between_peaks_values'][0]
    peak_idx = feature_values[0]['peak_indices'][0]

    expected = numpy.min(voltage[peak_idx:])

    numpy.testing.assert_allclose(min_btw_peaks_value, expected)


def test_AP_width_between_threshold():
    """basic: Test AP_width_between_threshold"""

    import efel
    efel.reset()

    threshold = -48
    efel.setDoubleSetting("Threshold", threshold)
    stim_start = 200.0
    stim_end = 1200.0

    time = efel.io.load_fragment('%s#col=1' % dendriticAP_url)
    voltage = efel.io.load_fragment('%s#col=2' % dendriticAP_url)
    time, voltage = interpolate(time, voltage, 0.1)
    trace = {}

    trace['T'] = time
    trace['V'] = voltage
    trace['stim_start'] = [stim_start]
    trace['stim_end'] = [stim_end]

    features = [
        'AP_width_between_threshold',
        'peak_indices',
        'min_between_peaks_indices'
    ]

    feature_values = \
        efel.getFeatureValues(
            [trace],
            features, raise_warnings=False)

    AP_width = feature_values[0]['AP_width_between_threshold'][0]
    peak_idx = feature_values[0]['peak_indices'][0]
    min_after_peak_idx = feature_values[0]['min_between_peaks_indices'][0]

    t0 = time[:peak_idx][voltage[:peak_idx] > threshold][0]
    t1 = time[peak_idx:min_after_peak_idx][
        voltage[peak_idx:min_after_peak_idx] < threshold
    ][0]

    numpy.testing.assert_allclose(AP_width, t1 - t0)


def test_AP_width_between_threshold_strict():
    """basic: Test AP_width_between_threshold with strict interval"""

    import efel
    efel.reset()
    efel.setIntSetting('strict_stiminterval', True)

    threshold = -48
    efel.setDoubleSetting("Threshold", threshold)
    stim_start = 200.0
    stim_end = 1200.0

    time = efel.io.load_fragment('%s#col=1' % dendriticAP_url)
    voltage = efel.io.load_fragment('%s#col=2' % dendriticAP_url)
    time, voltage = interpolate(time, voltage, 0.1)
    trace = {}

    trace['T'] = time
    trace['V'] = voltage
    trace['stim_start'] = [stim_start]
    trace['stim_end'] = [stim_end]

    features = ['AP_width_between_threshold']

    feature_values = \
        efel.getFeatureValues(
            [trace],
            features, raise_warnings=False)

    AP_width = feature_values[0]['AP_width_between_threshold']

    assert AP_width is None
