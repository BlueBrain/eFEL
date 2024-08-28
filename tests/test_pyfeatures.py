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


import h5py
from pathlib import Path
import numpy

import efel
from efel.io import load_ascii_input

testdata_dir = Path(__file__).parent / 'testdata'

traces_data = {
    'mean_frequency1': {
        'url': testdata_dir / 'basic' / 'mean_frequency_1.txt',
        'stim_start': 500,
        'stim_end': 900},
    'init_burst1': {
        'url': testdata_dir / 'basic' / 'init_burst1.txt',
        'stim_start': 250,
        'stim_end': 1600},
    'init_burst2': {
        'url': testdata_dir / 'basic' / 'init_burst2.txt',
        'stim_start': 250,
        'stim_end': 1600},
    'init_burst3': {
        'url': testdata_dir / 'basic' / 'init_burst3.txt',
        'stim_start': 250,
        'stim_end': 1600},
    'init_burst_sahp_error': {
        'url': testdata_dir / 'basic' / 'initburst_sahp_error.txt',
        'stim_start': 800,
        'stim_end': 1600},
    'depol_block_subthresh': {
        'url': testdata_dir / 'allfeatures' / 'testdb2data.txt',
        'stim_start': 419.995,
        'stim_end': 1419.995},
    'depol_block_subthresh_hyperpol': {
        'url': testdata_dir / 'allfeatures' / 'testdb1data.txt',
        'stim_start': 419.995,
        'stim_end': 1419.995},
    'depol_block_spiking': {
        'url': testdata_dir / 'allfeatures' / 'testdb3data.txt',
        'stim_start': 419.995,
        'stim_end': 1419.995},
    'depol_block_db': {
        'url': testdata_dir / 'allfeatures' / 'testdbdata.txt',
        'stim_start': 419.995,
        'stim_end': 1419.995},
    'current': {
        'url': testdata_dir / 'basic' / 'current.txt',
        't_col': 1,
        'i_col': 2,
        'v_col': 3,
        'stim_start': 700.0,
        'stim_end': 2700.0},
    'impedance': {
        'url': testdata_dir / 'basic' / 'impedance.txt',
        't_col': 1,
        'v_col': 2,
        'i_col': 3,
        'stim_start': 100.0,
        'stim_end': 5100.0},
    'voltage_clamp': {
        'url': testdata_dir / 'basic' / 'rCell9430.nwb',
        'act_stim_start': 99.5,
        'act_stim_end': 600.0,
        'deact_stim_start': 401.0,
        'deact_stim_end': 599.0,
        'inact_stim_start': 100.0,
        'inact_stim_end': 1600.0
    }
}


def _load_trace(trace_name):
    """Load trace with a certain name"""
    trace_data = traces_data[trace_name]

    url = trace_data['url']

    trace = {
        'stim_start': [trace_data['stim_start']],
        'stim_end': [trace_data['stim_end']],
    }

    if 'i_col' in trace_data:
        data = numpy.loadtxt(url)
        if trace_data['i_col'] == 3:  # I is the last
            trace['T'] = data[:, 0]
            trace['V'] = data[:, 1]
            trace['I'] = data[:, 2]
        else:  # V is the last
            trace['T'] = data[:, 0]
            trace['V'] = data[:, 2]
            trace['I'] = data[:, 1]
        if trace_name == "current":
            trace['T'] = trace['T'] * 1000.0  # s -> ms
    else:
        trace['T'], trace['V'] = load_ascii_input(url)

    return trace


def _load_voltage_clamp_traces(exp_type, repetition_name="repetition1"):
    """Load the voltage clamp traces from nwb file.

    Args:
        exp_type (str): type of the experiment. Can be
            'Activation', 'Deactivation' or 'Inactivation'.
    """
    trace_data = traces_data["voltage_clamp"]
    url = trace_data['url']
    if exp_type == "Activation":
        stim_start = trace_data["act_stim_start"]
        stim_end = trace_data["act_stim_end"]
    elif exp_type == "Deactivation":
        stim_start = trace_data["deact_stim_start"]
        stim_end = trace_data["deact_stim_end"]
    elif exp_type == "Inactivation":
        stim_start = trace_data["inact_stim_start"]
        stim_end = trace_data["inact_stim_end"]
    else:
        raise ValueError(
            "exp_type should be 'Activation', 'Deactivation' or 'Inactivation'"
            f", got {exp_type}"
        )
    data = {}
    with h5py.File(url, "r") as content:
        reps = content["acquisition"]["timeseries"][exp_type]["repetitions"]
        for rep_name, rep in reps.items():
            data[rep_name] = {
                "dt": numpy.array(rep["x_interval"], dtype="float32"),
                "current": numpy.array(rep["data"], dtype="float32"),
            }
    traces = []
    for idx in range(len(data[repetition_name]["dt"])):
        trace = {}
        i = data[repetition_name]["current"][:, idx]
        t = numpy.arange(i.size) * data[repetition_name]["dt"][idx]
        # efel expects ms: s -> ms
        t = t * 1000.0
        trace["T"] = t
        trace["V"] = i  # trick: input current as if it was voltage
        trace["stim_start"] = [stim_start]
        trace["stim_end"] = [stim_end]
        traces.append(trace)

    return traces


def _test_expected_value(feature_name, expected_values):
    """Test expected values for feature"""
    for trace_name, expected_value in expected_values.items():
        trace = _load_trace(trace_name)

        feature_values = efel.get_feature_values([trace], [feature_name])

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
    efel.reset()

    mf1_trace = _load_trace('mean_frequency1')

    feature_values = efel.get_feature_values([mf1_trace], ['ISI_values', 'ISIs'])

    assert numpy.allclose(
        feature_values[0]['ISIs'][1:],
        feature_values[0]['ISI_values'])

    numpy.testing.assert_allclose(
        efel.get_distance(
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


def test_depol_block_bool():
    """pyfeatures: Test depolarization block bool feature"""
    feature_name = 'depol_block_bool'
    expected_values = {
        'depol_block_subthresh': [0], 'depol_block_subthresh_hyperpol': [0],
        'depol_block_spiking': [0],
        'depol_block_db': [1]}

    _test_expected_value(feature_name, expected_values)


def test_pydistance():
    """pyfeatures: Test python distance."""
    mf1_trace = _load_trace('mean_frequency1')

    feature_name = 'AP_height'
    mean = 1.0
    std = 1.0

    numpy.seterr(divide='ignore')

    mf1_trace['stim_start'] = [0]
    mf1_trace['stim_end'] = [900]
    # with successful trace_check
    numpy.testing.assert_allclose(efel.get_distance(
        mf1_trace,
        feature_name,
        mean,
        std,
        trace_check=True), 30.422218394481284)

    efel.reset()
    # with failed trace_check
    mf1_trace["stim_end"] = [600]
    error_value = 250.0
    res = efel.get_distance(
        mf1_trace,
        feature_name,
        mean,
        std,
        trace_check=True,
        error_dist=error_value)
    assert res == error_value


def test_pydistance_featurefail():
    """pyfeatures: Test failure of feature in get_distance"""
    mf1_trace = _load_trace('mean_frequency1')

    feature_name = 'initburst_sahp'
    mean = 1.0
    std = 1.0

    efel.reset()
    numpy.testing.assert_allclose(efel.get_distance(
        mf1_trace,
        feature_name,
        mean,
        std,
        trace_check=True), 250.0)


def test_interpolate_current():
    """pyfeatures: Test interpolation of current"""
    def interpolate(time, voltage, new_dt):
        """Interpolate voltage to new dt"""

        interp_time = numpy.arange(time[0], time[-1] + new_dt, new_dt)
        interp_voltage = numpy.interp(interp_time, time, voltage)

        return interp_time, interp_voltage

    data = numpy.loadtxt(testdata_dir / 'basic' / 'current.txt')
    time = data[:, 0] * 1000.0  # -> ms
    current = data[:, 1]
    voltage = data[:, 2]

    feature_name = ['time', 'current', 'voltage']
    trace = _load_trace('current')
    feature_values = efel.get_feature_values([trace], feature_name)
    feature_time = feature_values[0]["time"]
    feature_current = feature_values[0]["current"]
    feature_voltage = feature_values[0]["voltage"]
    interp_time, interp_current = interpolate(time, current, new_dt=0.1)
    _, interp_voltage = interpolate(time, voltage, new_dt=0.1)

    assert len(interp_time) == len(feature_time)
    assert len(interp_current) == len(feature_current)
    assert len(interp_voltage) == len(feature_voltage)
    assert len(feature_voltage) == len(feature_current)
    assert numpy.allclose(interp_current, feature_current, atol=1e-6)


def test_impedance():
    """pyfeatures: Test impedance feature"""
    feature_name = "impedance"

    expected_values = {feature_name: 4.80769231}
    _test_expected_value(feature_name, expected_values)


def test_trace_check():
    """Unit test for trace_check."""
    efel.reset()

    mf1_trace = _load_trace('mean_frequency1')
    feature_values = efel.get_feature_values([mf1_trace], ['trace_check'])
    assert feature_values[0]['trace_check'][0] == 0
    # failure
    mf1_trace["stim_end"] = [600]
    feature_values = efel.get_feature_values([mf1_trace], ['trace_check'])
    assert feature_values[0]['trace_check'] is None
    # no spikes
    mf1_trace["V"] = [0] * len(mf1_trace["V"])
    feature_values = efel.get_feature_values([mf1_trace], ['trace_check'])
    assert feature_values[0]['trace_check'][0] == 0


def test_phaseslope_max():
    """Unit test for phaseslope_max."""
    expected_values = {
        "mean_frequency1": 574.86769012,
        "depol_block_spiking": 93.51242208,
        "depol_block_db": 180.7325033,
    }
    _test_expected_value("phaseslope_max", expected_values)


def test_activation_time_constant():
    """Unit test for activation_time_constant."""
    efel.reset()

    traces = _load_voltage_clamp_traces("Activation")
    feats = efel.get_feature_values(traces, ["activation_time_constant"])
    act_tau = [feat_dict["activation_time_constant"][0] for feat_dict in feats]
    act_tau_ref = [
        4.42464151e-02, 6.36560480e-02, 1.35239568e+00, 6.788478e-08,
        2.00926636e-06, 1.08786116e+02, 3.51324930e+01, 1.01711405e+01,
        4.61677565e+00, 2.83639400e+00, 1.97298305e+00, 1.59010996e+00,
        1.21969111e+00, 1.07930254e+00, 8.55734307e-01, 7.39074539e-01,
        6.58848184e-01, 5.96009883e-01
    ]
    numpy.testing.assert_allclose(act_tau, act_tau_ref, rtol=1e-6, atol=1e-10)


def test_deactivation_time_constant():
    """Unit test for deactivation_time_constant."""
    efel.reset()

    traces = _load_voltage_clamp_traces("Deactivation")
    feats = efel.get_feature_values(traces, ["deactivation_time_constant"])
    deact_tau = [feat_dict["deactivation_time_constant"][0] for feat_dict in feats]
    deact_tau_ref = [
        1., 18.428159, 17.22834973, 27.27256786, 39.47847711,
        49.2782552, 65.81509027, 71.9253926, 81.54669955, 107.68102719,
        134.8814237, 120.34484955
    ]
    numpy.testing.assert_allclose(deact_tau, deact_tau_ref)


def test_inactivation_time_constant():
    """Unit test for inactivation_time_constant."""
    efel.reset()

    traces = _load_voltage_clamp_traces("Inactivation")
    feats = efel.get_feature_values(traces, ["inactivation_time_constant"])
    inact_tau = [feat_dict["inactivation_time_constant"][0] for feat_dict in feats]
    inact_tau_ref = [
        88.05764617, 397.85146878, 431.22382104, 239.30375528, 175.7584476,
        168.54772458, 143.78127593, 146.63424961, 145.0907385, 167.04060933,
        142.24530676, 137.62498461
    ]
    numpy.testing.assert_allclose(inact_tau, inact_tau_ref, rtol=1e-6)
