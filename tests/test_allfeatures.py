"""Test all features on an example trace"""

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

import json
import os
import warnings

import numpy as np

import efel

testdata_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            'testdata',
                            'allfeatures')


def get_allfeature_values():
    """Get back all the feature names and value"""
    efel.reset()

    all_featurenames = efel.get_feature_names()

    def load_data(filename):
        data = np.loadtxt(os.path.join(testdata_dir, filename))
        return data[:, 0], data[:, 1]

    soma_time, soma_voltage = load_data('testdata.txt')
    db_time, db_voltage = load_data('testdbdata.txt')
    bac_time, bac_voltage = load_data('testbacdata.txt')
    bap1_time, bap1_voltage = load_data('testbap1data.txt')
    bap2_time, bap2_voltage = load_data('testbap2data.txt')

    traces = [
        {
            'T': soma_time,
            'V': soma_voltage,
            'stim_start': [700],
            'stim_end': [2700]
        },
        {
            'T': soma_time,
            'V': soma_voltage,
            'stim_start': [700],
            'stim_end': [2700],
        },
        {
            'T': bac_time,
            'V': bac_voltage,
            'stim_start': [295],
            'stim_end': [600]
        },
        {
            'T': bap1_time,
            'V': bap1_voltage,
            'stim_start': [295],
            'stim_end': [500]
        },
        {
            'T': bap2_time,
            'V': bap2_voltage,
            'stim_start': [295],
            'stim_end': [500]
        }
    ]

    trace_db = {
        'T': db_time,
        'V': db_voltage,
        'stim_start': [419.995],
        'stim_end': [1419.995]
    }

    db_featurenames = ['depol_block', 'depol_block_bool']

    soma_featurenames = [
        x
        for x in all_featurenames
        if x
        not in db_featurenames
        + ["current", "current_base", "impedance", "steady_state_current_stimend"]
    ]

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        feature_values = efel.get_feature_values(traces, soma_featurenames)[0]

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        feature_values.update(efel.get_feature_values([trace_db], db_featurenames)[0])

    for feature_name in feature_values:
        if feature_values[feature_name] is not None:
            feature_values[feature_name] = list(feature_values[feature_name])

    return feature_values


def test_allfeatures():
    """allfeatures: Regression testing all features on a trace"""

    feature_values = get_allfeature_values()
    # drop Spikecount and Spikecount_stimint deprecated features
    feature_values.pop('Spikecount')
    feature_values.pop('Spikecount_stimint')
    # remove features that expects voltage clamp trace
    feature_values.pop('activation_time_constant')
    feature_values.pop('deactivation_time_constant')
    feature_values.pop('inactivation_time_constant')
    test_data_path = os.path.join(testdata_dir, 'expectedresults.json')
    with open(test_data_path, 'r') as expected_json:
        expected_results = json.load(expected_json)

    assert set(feature_values.keys()) == set(expected_results.keys())
    failed_feature = False
    for feature_name, feature_value in feature_values.items():
        expected_value = expected_results[feature_name]
        if feature_name is None:
            equal = (expected_value is None)
        if expected_value is None:
            equal = (feature_value is None)
        elif feature_value is None:
            equal = (expected_value is None)
        else:
            equal = (len(feature_value) == len(expected_value)) \
                and np.allclose(feature_value, expected_value)

        if not equal:
            print("Difference in feature %s: value=%s expected=%s" %
                  (feature_name, feature_value, expected_value))
            failed_feature = True

    assert not failed_feature


def test_allfeatures_on_constant_voltage():
    """Call all features on constant voltage input."""
    time = np.linspace(0, 999, 1000)
    voltage = np.full(1000, -80.0)

    efel.reset()
    traces = [{'T': time, 'V': voltage, 'stim_start': [100], 'stim_end': [999]}]
    all_featurenames = efel.get_feature_names()
    feature_values = efel.get_feature_values(traces, all_featurenames)[0]
    assert all(feature_values["voltage"] == -80.0)
    # Assert that each element in time is greater than or equal to the previous element
    assert np.all(feature_values["time"][1:] >= feature_values["time"][:-1])

    # Assert for array fields to be non-empty arrays
    # If you add a new feature in the future that results in an array, add here
    array_fields = [
        "minimum_voltage", "maximum_voltage", "maximum_voltage_from_voltagebase",
        "sag_amplitude", "voltage_after_stim", "steady_state_hyper",
        "steady_state_voltage", "steady_state_voltage_stimend",
        "voltage_deflection", "voltage_deflection_begin", "voltage_deflection_vb_ssse",
        "depol_block", "depol_block_bool", "voltage_base", "Spikecount",
        "Spikecount_stimint", "burst_number", "strict_burst_number", "trace_check",
        "spike_count", "spike_count_stimint", "phaseslope_max",
        "activation_time_constant", "deactivation_time_constant",
        "inactivation_time_constant"
    ]

    for field in array_fields:
        if field in feature_values:
            assert isinstance(feature_values[field], np.ndarray)
            assert len(feature_values[field]) > 0

    # Assert the rest of the fields are None
    excluded_fields = array_fields + ["voltage", "time"]
    for key in feature_values:
        if key not in excluded_fields:
            assert feature_values[key] is None, f"Field {key} is not None"
