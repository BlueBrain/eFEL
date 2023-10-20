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

import os
import warnings

# pylint: disable=R0914

testdata_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            'testdata',
                            'allfeatures')


def get_allfeature_values():
    """Get back all the feature names and value"""

    import efel
    import numpy
    efel.reset()

    all_featurenames = efel.getFeatureNames()

    def load_data(filename):
        data = numpy.loadtxt(os.path.join(testdata_dir, filename))
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

    bpap_featurenames = [
        'BPAPHeightLoc1',
        'BPAPHeightLoc2',
        'BPAPAmplitudeLoc1',
        'BPAPAmplitudeLoc2'
    ]

    bac_featurenames = ['BAC_width']
    db_featurenames = ['depol_block', 'depol_block_bool']

    soma_featurenames = [
        x
        for x in all_featurenames
        if x
        not in bpap_featurenames
        + bac_featurenames
        + db_featurenames
        + ["current", "current_base", "impedance"]
    ]

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        feature_values = efel.getFeatureValues(traces, soma_featurenames)[0]

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        feature_values.update(efel.getFeatureValues([trace_db], db_featurenames)[0])

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        efel.setThreshold(-30)
        feature_values.update(efel.getFeatureValues(traces, bpap_featurenames)[0])

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        efel.setThreshold(-55)
        feature_values.update(efel.getFeatureValues(traces, bac_featurenames)[0])

    for feature_name in feature_values:
        if feature_values[feature_name] is not None:
            feature_values[feature_name] = list(feature_values[feature_name])

    return feature_values


def test_allfeatures():
    """allfeatures: Regression testing all features on a trace"""

    feature_values = get_allfeature_values()

    import json
    test_data_path = os.path.join(testdata_dir, 'expectedresults.json')
    with open(test_data_path, 'r') as expected_json:
        expected_results = json.load(expected_json)

    import numpy
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
                and numpy.allclose(feature_value, expected_value)

        if not equal:
            print("Difference in feature %s: value=%s expected=%s" %
                  (feature_name, feature_value, expected_value))
            failed_feature = True

    assert not failed_feature
