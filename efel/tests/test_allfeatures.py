"""Test all features on an example trace"""

import nose.tools as nt

# pylint: disable=R0914


def get_allfeature_values():
    """Get back all the feature names and value"""

    import efel
    import numpy

    all_featurenames = efel.getFeatureNames()

    soma_data = numpy.loadtxt('testdata/allfeatures/testdata.txt')
    soma_time = soma_data[:, 0]
    soma_voltage = soma_data[:, 1]

    bac_data = numpy.loadtxt('testdata/allfeatures/testbacdata.txt')
    bac_time = bac_data[:, 0]
    bac_voltage = bac_data[:, 1]

    bap1_data = numpy.loadtxt('testdata/allfeatures/testbap1data.txt')
    bap1_time = bap1_data[:, 0]
    bap1_voltage = bap1_data[:, 1]

    bap2_data = numpy.loadtxt('testdata/allfeatures/testbap2data.txt')
    bap2_time = bap2_data[:, 0]
    bap2_voltage = bap2_data[:, 1]

    trace = {}

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

    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        feature_values = efel.getFeatureValues([trace], all_featurenames)[0]

    for feature_name in feature_values:
        if feature_values[feature_name] is not None:
            feature_values[feature_name] = list(
                feature_values[feature_name])

    return feature_values


def test_allfeatures():
    """allfeatures: Regression testing all features on a trace"""

    feature_values = get_allfeature_values()

    import json
    with open('testdata/allfeatures/expectedresults.json', 'r') \
            as expected_json:
        expected_results = json.load(expected_json)

    all_equal = True
    for feature_name, feature_value in feature_values.items():
        # feature_value = feature_values[feature_name]
        expected_value = expected_results[feature_name]
        if feature_value != expected_value:
            print "Difference in feature %s: value=%s expected=%s" % \
                (feature_name, feature_value, expected_value)
            all_equal = False

    nt.assert_true(all_equal)
