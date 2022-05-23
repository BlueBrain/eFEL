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

# pylint: disable=R0914

testdata_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            'testdata',
                            'allfeatures')


def get_allfeature_values():
    """Get back all the feature names and value"""

    import efel
    efel.reset()
    import numpy

    all_featurenames = efel.getFeatureNames()

    soma_data = numpy.loadtxt(os.path.join(testdata_dir, 'testdata.txt'))
    soma_time = soma_data[:, 0]
    soma_voltage = soma_data[:, 1]

    db_data = numpy.loadtxt(os.path.join(testdata_dir, 'testdbdata.txt'))
    db_time = db_data[:, 0]
    db_voltage = db_data[:, 1]

    bac_data = numpy.loadtxt(os.path.join(testdata_dir, 'testbacdata.txt'))
    bac_time = bac_data[:, 0]
    bac_voltage = bac_data[:, 1]

    bap1_data = numpy.loadtxt(os.path.join(testdata_dir, 'testbap1data.txt'))
    bap1_time = bap1_data[:, 0]
    bap1_voltage = bap1_data[:, 1]

    bap2_data = numpy.loadtxt(os.path.join(testdata_dir, 'testbap2data.txt'))
    bap2_time = bap2_data[:, 0]
    bap2_voltage = bap2_data[:, 1]

    trace = {}
    trace_db = {}

    trace['T'] = soma_time
    trace['V'] = soma_voltage
    trace['stim_start'] = [700]
    trace['stim_end'] = [2700]
    trace['T;location_AIS'] = soma_time
    trace['V;location_AIS'] = soma_voltage
    trace['stim_start;location_AIS'] = [700]
    trace['stim_end;location_AIS'] = [2700]
    trace['T;location_epsp'] = bac_time
    trace['V;location_epsp'] = bac_voltage
    trace['stim_start;location_epsp'] = [295]
    trace['stim_end;location_epsp'] = [600]
    trace['T;location_dend1'] = bap1_time
    trace['V;location_dend1'] = bap1_voltage
    trace['stim_start;location_dend1'] = [295]
    trace['stim_end;location_dend1'] = [500]
    trace['T;location_dend2'] = bap2_time
    trace['V;location_dend2'] = bap2_voltage
    trace['stim_start;location_dend2'] = [295]
    trace['stim_end;location_dend2'] = [500]

    trace_db['T'] = db_time
    trace_db['V'] = db_voltage
    trace_db['stim_start'] = [419.995]
    trace_db['stim_end'] = [1419.995]

    bpap_featurenames = [
        'BPAPHeightLoc1',
        'BPAPHeightLoc2',
        'BPAPAmplitudeLoc1',
        'BPAPAmplitudeLoc2']

    bac_featurenames = [
        'BAC_width']

    db_featurenames = [
        'depol_block']

    soma_featurenames = all_featurenames[:]

    for feature_name in bpap_featurenames:
        soma_featurenames.remove(feature_name)

    for feature_name in bac_featurenames:
        soma_featurenames.remove(feature_name)

    for feature_name in db_featurenames:
        soma_featurenames.remove(feature_name)

    soma_featurenames = [x for x in soma_featurenames
                         if x not in ['current', 'current_base']]

    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        feature_values = efel.getFeatureValues([trace], soma_featurenames)[0]

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        feature_values = dict(
            list(feature_values.items()) +
            list(efel.getFeatureValues(
                [trace_db],
                db_featurenames)[0].items()))

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        efel.setThreshold(-30)
        feature_values = dict(
            list(feature_values.items()) +
            list(efel.getFeatureValues(
                [trace],
                bpap_featurenames)[0].items()))

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        efel.setThreshold(-55)
        feature_values = dict(
            list(feature_values.items()) +
            list(efel.getFeatureValues(
                [trace],
                bac_featurenames)[0].items()))

    for feature_name in feature_values:
        if feature_values[feature_name] is not None:
            feature_values[feature_name] = list(
                feature_values[feature_name])

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
