"""Test all features on an example trace"""

import os

if __name__ == '__main__':
    import test_allfeatures
    feature_values = test_allfeatures.get_allfeature_values()

    testdata_dir = os.path.join(
        os.path.dirname(
            os.path.abspath(__file__)),
        'testdata/allfeatures')
    import json
    with open(os.path.join(testdata_dir, 'expectedresults.json'), 'w') \
            as expected_json:
        json.dump(
            feature_values,
            expected_json,
            indent=4,
            separators=(
                ',',
                ': '))
