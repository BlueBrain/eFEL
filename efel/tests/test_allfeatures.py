"""Test all features on an example trace"""

import nose.tools as nt


def test_allfeatures():
    """allfeatures: Regression testing all features on a trace"""
    import efel
    import numpy

    all_featurenames = efel.getFeatureNames()

    data = numpy.loadtxt('testdata/allfeatures/testdata.txt')

    time = data[:, 0]
    voltage = data[:, 1]

    trace = {}

    trace['T'] = time
    trace['V'] = voltage
    trace['stim_start'] = [700]
    trace['stim_end'] = [2700]

    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        feature_values = efel.getFeatureValues([trace], all_featurenames)[0]

    for feature_name in feature_values:
        if feature_values[feature_name] is not None:
            feature_values[feature_name] = list(
                feature_values[feature_name])
    import json
    with open('testdata/allfeatures/expectedresults.json', 'r') \
            as expected_json:
        expected_results = json.load(expected_json)

    for feature_name in all_featurenames:
        nt.assert_equal(
            feature_values[feature_name],
            expected_results[feature_name])
