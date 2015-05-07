"""General tests of eFEL"""

import nose.tools as nt

# pylint: disable=W0611, W0612, F0401


def test_import():
    """Test importing of eFEL"""

    # pylint: disable=W0611
    import efel
    # pylint: enable=W0611


def test_empty_trace():
    """Testing results for empty trace"""

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

    print "Here"
    efel.getFeatureValues([trace], features)
    print "Here2"

    for feature, value in \
            efel.getFeatureValues([trace], features)[0].iteritems():

        nt.assert_equal(value[0], 0.0)


def test_ISI_log_slope_skip():
    """Test ISI_log_slope_skip"""

    import efel
    import numpy

    stim_start = 31.2
    stim_end = 431.2

    data = numpy.loadtxt('testdata/zero_ISI_log_slope_skip95824004.abf.csv')

    time = data[:, 0]
    voltage = data[:, 1]

    trace = {}

    trace['T'] = time
    trace['V'] = voltage
    trace['stim_start'] = [stim_start]
    trace['stim_end'] = [stim_end]

    features = ['ISI_log_slope_skip']

    feature_values = \
        efel.getFeatureValues(
            [trace],
            features)
    nt.assert_equal(feature_values[0]['ISI_log_slope_skip'], None)


def test_mean_frequency1():
    """Test mean_frequency 1"""

    import efel
    import numpy

    stim_start = 500.0
    stim_end = 900.0

    data = numpy.loadtxt('testdata/mean_frequency_1.txt')

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


def test_ap_amplitude_from_voltagebase1():
    """Test AP_amplitude_from_voltagebase 1"""

    import efel
    import numpy

    stim_start = 500.0
    stim_end = 900.0

    data = numpy.loadtxt('testdata/mean_frequency_1.txt')

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


def test_APlast_amp():
    """Test APlast_amp"""

    import efel
    import numpy

    stim_start = 500.0
    stim_end = 900.0

    data = numpy.loadtxt('testdata/mean_frequency_1.txt')

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
    """Test min_voltage_between_spikes 1"""

    import efel
    import numpy

    stim_start = 500.0
    stim_end = 900.0

    data = numpy.loadtxt('testdata/mean_frequency_1.txt')

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

    """
    for peak_voltage, ap_amplitude_from_voltagebase in zip(
            feature_values[0]['peak_voltage'],
            feature_values[0]['AP_amplitude_from_voltagebase']):
        nt.assert_almost_equal(peak_voltage - voltage_base,
                               ap_amplitude_from_voltagebase)
    """
