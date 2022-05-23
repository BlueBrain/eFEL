"""Tests for the python features of eFEL"""

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

import efel

testdata_dir = os.path.join(
    os.path.dirname(
        os.path.abspath(__file__)),
    'testdata')

traces_data = {
    'mean_frequency1': {
        'url': 'file://%s' % os.path.join(
            os.path.abspath(testdata_dir),
            'basic',
            'mean_frequency_1.txt'),
        't_col': 1,
        'v_col': 2,
        'stim_start': 500,
        'stim_end': 900},
    'init_burst1': {
        'url': 'file://%s' % os.path.join(
            os.path.abspath(testdata_dir),
            'basic',
            'init_burst1.txt'),
        't_col': 1,
        'v_col': 2,
        'stim_start': 250,
        'stim_end': 1600},
    'init_burst2': {
        'url': 'file://%s' % os.path.join(
            os.path.abspath(testdata_dir),
            'basic',
            'init_burst2.txt'),
        't_col': 1,
        'v_col': 2,
        'stim_start': 250,
        'stim_end': 1600},
    'init_burst3': {
        'url': 'file://%s' % os.path.join(
            os.path.abspath(testdata_dir),
            'basic',
            'init_burst3.txt'),
        't_col': 1,
        'v_col': 2,
        'stim_start': 250,
        'stim_end': 1600},
    'init_burst_sahp_error': {
        'url': 'file://%s' % os.path.join(
            os.path.abspath(testdata_dir),
            'basic',
            'initburst_sahp_error.txt'),
        't_col': 1,
        'v_col': 2,
        'stim_start': 800,
        'stim_end': 1600},
    'depol_block_subthresh': {
        'url': 'file://%s' % os.path.join(
            os.path.abspath(testdata_dir),
            'allfeatures',
            'testdb2data.txt'),
        't_col': 1,
        'v_col': 2,
        'stim_start': 419.995,
        'stim_end': 1419.995},
    'depol_block_subthresh_hyperpol': {
        'url': 'file://%s' % os.path.join(
            os.path.abspath(testdata_dir),
            'allfeatures',
            'testdb1data.txt'),
        't_col': 1,
        'v_col': 2,
        'stim_start': 419.995,
        'stim_end': 1419.995},
    'depol_block_spiking': {
        'url': 'file://%s' % os.path.join(
            os.path.abspath(testdata_dir),
            'allfeatures',
            'testdb3data.txt'),
        't_col': 1,
        'v_col': 2,
        'stim_start': 419.995,
        'stim_end': 1419.995},
    'depol_block_db': {
        'url': 'file://%s' % os.path.join(
            os.path.abspath(testdata_dir),
            'allfeatures',
            'testdbdata.txt'),
        't_col': 1,
        'v_col': 2,
        'stim_start': 419.995,
        'stim_end': 1419.995},
    'current': {
        'url': 'file://%s' % os.path.join(
            os.path.abspath(testdata_dir),
            'basic',
            'current.txt'),
        't_col': 1,
        'i_col': 2,
        'v_col': 3,
        'stim_start': 700.0,
        'stim_end': 2700.0},
}


def _load_trace(trace_name):
    """Load trace with a certain name"""

    trace_data = traces_data[trace_name]

    url = trace_data['url']

    trace = {
        'T': efel.io.load_fragment(
            '%s#col=%d' %
            (url, trace_data['t_col'])),
        'V': efel.io.load_fragment(
            '%s#col=%d' %
            (url, trace_data['v_col'])),
        'stim_start': [trace_data['stim_start']],
        'stim_end': [trace_data['stim_end']],
    }

    if 'i_col' in trace_data:
        trace['I'] = efel.io.load_fragment('%s#col=%d' %
                                           (url, trace_data['i_col']))

    return trace


def _test_expected_value(feature_name, expected_values):
    """Test expected values for feature"""

    for trace_name, expected_value in expected_values.items():
        trace = _load_trace(trace_name)

        feature_values = efel.getFeatureValues([trace], [feature_name])

        if expected_value is None:
            assert feature_values[0][feature_name] is None
        else:
            assert numpy.allclose(
                feature_values[0][feature_name],
                expected_value
            )


def test_initburst_sahp():
    """pyfeatures: Test initburst_sahp feature"""

    feature_name = 'initburst_sahp'
    expected_values = {
        'mean_frequency1': None, 'init_burst1': [-69.19999695],
        'init_burst2': [-74.34375],
        'init_burst3': [-76.72435972]}

    _test_expected_value(feature_name, expected_values)


def test_initburst_sahp_argmin_error():
    """pyfeatures: Test initburst_sahp argmin error"""

    feature_name = 'initburst_sahp'
    expected_values = {
        'init_burst_sahp_error': None}

    _test_expected_value(feature_name, expected_values)


def test_initburst_sahp_vb():
    """pyfeatures: Test initburst_sahp_vb feature"""

    feature_name = 'initburst_sahp_vb'
    expected_values = {
        'mean_frequency1': None, 'init_burst1': [13.80537756],
        'init_burst2': [11.43360025],
        'init_burst3': [6.08649806]}

    _test_expected_value(feature_name, expected_values)


def test_initburst_sahp_ssse():
    """pyfeatures: Test initburst_sahp_ssse feature"""

    feature_name = 'initburst_sahp_ssse'
    expected_values = {
        'mean_frequency1': None, 'init_burst1': [-11.33152931],
        'init_burst2': [2.02011574],
        'init_burst3': [-4.42756346]}

    _test_expected_value(feature_name, expected_values)


def test_ISIs():
    """pyfeatures: Test ISIs feature"""

    mf1_trace = _load_trace('mean_frequency1')

    feature_values = efel.getFeatureValues([mf1_trace], ['ISI_values', 'ISIs'])

    assert numpy.allclose(
        feature_values[0]['ISIs'][1:],
        feature_values[0]['ISI_values'])

    numpy.testing.assert_allclose(
        efel.getDistance(
            mf1_trace,
            'ISIs',
            1.0,
            1.0),
        64.25000000001484)


def test_depol_block():
    """pyfeatures: Test depolarization block feature"""
    feature_name = 'depol_block'
    expected_values = {
        'depol_block_subthresh': [1], 'depol_block_subthresh_hyperpol': [1],
        'depol_block_spiking': [1],
        'depol_block_db': None}

    _test_expected_value(feature_name, expected_values)


def test_pydistance():
    """pyfeatures: Test python distance against cpp version"""
    mf1_trace = _load_trace('mean_frequency1')

    feature_name = 'AP_height'
    mean = 1.0
    std = 1.0

    numpy.seterr(divide='ignore')

    # Check if cpp and python the same if:
    # - baseline
    # - std = 0.0
    # - trace_check is enabled
    # - trace_check is enabled on faulty trace
    # - trace_check is disabled on faulty trace
    for args, stim_end in [
        ((mf1_trace, feature_name, mean, std, None), 900),
        ((mf1_trace, feature_name, mean, 0.0, None), 900),
        ((mf1_trace, feature_name, mean, std, True), 900),
        ((mf1_trace, feature_name, mean, std, True), 600),
        ((mf1_trace, feature_name, mean, std, False), 600),
    ]:
        efel.reset()
        mf1_trace['stim_end'] = [stim_end]
        assert (
            efel.getDistance(*args) == efel.api._getDistance_cpp(*args))

    # Extra sanity checks for trace_check
    mf1_trace['stim_end'] = [600]

    efel.reset()
    numpy.testing.assert_allclose(efel.getDistance(
        mf1_trace,
        feature_name,
        mean,
        std,
        trace_check=False), 30.422218394481284)

    efel.reset()
    numpy.testing.assert_allclose(efel.api._getDistance_cpp(
        mf1_trace,
        feature_name,
        mean,
        std,
        trace_check=True), 250.0)


def test_pydistance_featurefail():
    """pyfeatures: Test failure of feature in getdistance"""

    mf1_trace = _load_trace('mean_frequency1')

    feature_name = 'initburst_sahp'
    mean = 1.0
    std = 1.0

    efel.reset()
    numpy.testing.assert_allclose(efel.getDistance(
        mf1_trace,
        feature_name,
        mean,
        std,
        trace_check=True), 250.0)


def test_current():
    """pyfeatures: Test current feature"""

    feature_name = 'current'
    data = numpy.loadtxt(os.path.join(os.path.abspath(testdata_dir),
                                      'basic',
                                      'current.txt'))
    current = data[:, 1]
    expected_values = {'current': current}
    _test_expected_value(feature_name, expected_values)


def test_interpolate_current():
    """pyfeatures: Test interpolation of current"""

    def interpolate(time, voltage, new_dt):
        """Interpolate voltage to new dt"""

        interp_time = numpy.arange(time[0], time[-1] + new_dt, new_dt)
        interp_voltage = numpy.interp(interp_time, time, voltage)

        return interp_time, interp_voltage

    data = numpy.loadtxt(os.path.join(os.path.abspath(testdata_dir),
                                      'basic',
                                      'current.txt'))
    time = data[:, 0]
    current = data[:, 1]
    voltage = data[:, 2]

    feature_name = ['time', 'current', 'voltage']
    trace = _load_trace('current')
    feature_values = efel.getFeatureValues([trace], ['current'])
    interp_time, interp_current = interpolate(time, current, new_dt=0.00025)

    assert len(interp_time) == len(time)
    assert len(interp_current) == len(current)
    assert len(voltage) == len(current)
    assert numpy.allclose(interp_current, current)
