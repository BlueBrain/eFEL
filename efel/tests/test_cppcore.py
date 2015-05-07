"""Test all features on an example trace"""

import nose.tools as nt


def test_import():
    """cppcore: Testing import of cppcore"""
    import efel.cppcore


class TestCppcore(object):

    def setup(self):
        """Setup"""
        import efel
        efel.cppcore.Initialize(efel.settings.dependencyfile_path, "log")

    def test_getFeatureNames(self):
        """cppcore: Testing getting all feature names"""
        import efel.cppcore
        feature_names = []
        efel.cppcore.getFeatureNames(feature_names)

        import json
        with open('featurenames.json', 'r') as featurenames_json:
            expected_featurenames = json.load(featurenames_json)
        nt.assert_equal(feature_names, expected_featurenames)

if __name__ == '__main__':
    test_getFeatureNames()
