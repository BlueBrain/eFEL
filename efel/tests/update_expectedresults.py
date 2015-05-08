"""Test all features on an example trace"""

import nose.tools as nt


if __name__ == '__main__':
    import test_allfeatures
    feature_values = test_allfeatures.get_allfeature_values()

    import json
    with open('testdata/allfeatures/expectedresults.json', 'w') \
            as expected_json:
        json.dump(
            feature_values,
            expected_json,
            indent=4,
            separators=(
                ',',
                ': '))
