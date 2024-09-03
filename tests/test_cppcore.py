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
from pathlib import Path
import tempfile

import numpy as np
import pytest

testdata_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            'testdata/')


def test_import():
    """cppcore: Testing import of cppcore"""
    import efel.cppcore  # NOQA


class TestCppcore:

    """Test cppcore"""

    def setup_method(self):  # pylint: disable=R0201
        """Setup"""
        import efel
        efel.cppcore.Initialize(efel.getDependencyFileLocation(), "log")

    def setup_data(self):  # pylint: disable=R0201
        """Set up data"""
        import efel
        stim_start = 500.0
        stim_end = 900.0

        test_data_path = os.path.join(
            testdata_dir,
            'basic/mean_frequency_1.txt')
        data = np.loadtxt(test_data_path)

        time = data[:, 0]
        voltage = data[:, 1]

        efel.cppcore.setFeatureDouble('T', [x for x in time])
        efel.cppcore.setFeatureDouble('V', [x for x in voltage])
        efel.cppcore.setFeatureDouble('stim_start', [stim_start])
        efel.cppcore.setFeatureDouble('stim_end', [stim_end])

        efel.cppcore.setFeatureDouble('spike_skipf', [0.1])
        efel.cppcore.setFeatureInt('max_spike_skip', [2])
        efel.cppcore.setFeatureDouble('Threshold',
                                      [-20.0])
        efel.cppcore.setFeatureDouble('DerivativeThreshold',
                                      [10.0])
        efel.cppcore.setFeatureDouble('interp_step', [0.1])
        efel.cppcore.setFeatureDouble('burst_factor', [1.5])
        efel.cppcore.setFeatureDouble("initial_perc", [0.1])

        efel.cppcore.setFeatureInt("DerivativeWindow", [3])

    def test_get_feature_names(self):  # pylint: disable=R0201
        """cppcore: Testing getting all feature names"""
        import efel
        feature_names = []
        efel.cppcore.getFeatureNames(feature_names)

        feature_names += efel.pyfeatures.all_pyfeatures

        import json
        test_data_path = os.path.join(testdata_dir, '../featurenames.json')
        with open(test_data_path, 'r') as featurenames_json:
            expected_featurenames = json.load(featurenames_json)
        # add the new names for the deprecated ones
        expected_featurenames += ["Spikecount", "Spikecount_stimint"]
        assert set(feature_names) == set(expected_featurenames)

    def test_getFeatureInt(self):
        """cppcore: Testing getFeatureInt"""
        import efel
        self.setup_data()
        # get int feature
        feature_values = list()
        efel.cppcore.getFeatureInt('AP_begin_indices', feature_values)
        assert isinstance(feature_values[0], int)

    def test_getFeatureDouble(self):
        """cppcore: Testing getFeatureDouble."""
        import efel
        self.setup_data()
        # get double feature
        feature_values = list()
        efel.cppcore.getFeatureDouble('AP_amplitude', feature_values)
        assert isinstance(feature_values[0], float)

    def test_getFeatureDouble_failure(self):  # pylint: disable=R0201
        """cppcore: Testing failure exit code in getFeatureDouble"""
        import efel
        feature_values = list()
        return_value = efel.cppcore.getFeatureDouble(
            "AP_amplitude",
            feature_values)
        assert return_value == -1

    @pytest.mark.xfail(raises=TypeError)
    def test_getFeatureDouble_wrong_type(self):  # pylint: disable=R0201
        """cppcore: Testing getFeatureDouble with wrong type"""
        import efel
        efel.cppcore.getFeatureDouble("AP_fall_indices", list())

    @pytest.mark.xfail(raises=TypeError)
    def test_getFeatureInt_wrong_type(self):  # pylint: disable=R0201
        """cppcore: Testing getFeatureInt with wrong type"""
        import efel
        efel.cppcore.getFeatureInt("AP_amplitude", list())

    def test_getFeature(self):
        """cppcore: Testing getFeature"""
        import efel
        self.setup_data()

        # get double feature
        feature_values = list()
        efel.cppcore.getFeature('AP_amplitude', feature_values)
        assert isinstance(feature_values[0], float)
        assert 5 == len(feature_values)
        assert np.allclose([80.45724099440199, 80.46320199354948,
                            80.73300299176428, 80.9965359926715,
                            81.87292599493423], feature_values)
        # get int feature
        feature_values = list()
        efel.cppcore.getFeature('AP_fall_indices', feature_values)
        assert isinstance(feature_values[0], int)
        assert 5 == len(feature_values)
        assert [5665, 6066, 6537, 7170, 8275] == feature_values

    def test_getFeature_failure(self):  # pylint: disable=R0201
        """cppcore: Testing failure exit code in getFeature"""
        import efel
        feature_values = list()
        return_value = efel.cppcore.getFeature("AP_amplitude", feature_values)
        assert return_value == -1

    @pytest.mark.xfail(raises=RuntimeError)
    def test_getFeature_non_existant(self):  # pylint: disable=R0201
        """cppcore: Testing failure exit code in getFeature"""
        import efel
        efel.cppcore.getFeature("does_not_exist", list())

    def test_getMapIntData(self):
        """cppcore: Testing getMapIntData."""
        import efel
        self.setup_data()
        # with non-existent key
        res = efel.cppcore.getMapIntData('AP_begin_indices')
        assert res == []
        feature_values = list()
        efel.cppcore.getFeatureInt('AP_begin_indices', feature_values)
        res = efel.cppcore.getMapIntData('AP_begin_indices')
        assert res == [5655, 6057, 6527, 7161, 8266]

    def test_getMapDoubleData(self):
        """cppcore: Testing getMapDoubleData."""
        import efel
        self.setup_data()
        # with non-existent key
        res = efel.cppcore.getMapDoubleData('AP_amplitude')
        assert res == []
        feature_values = list()
        efel.cppcore.getFeatureDouble('AP_amplitude', feature_values)
        res = efel.cppcore.getMapDoubleData('AP_amplitude')
        assert res[0] == 80.45724099440199
        assert len(res) == 5

    def test_featuretype(self):
        """cppcore: Testing featuretype."""
        import efel
        assert efel.cppcore.featuretype('AP_amplitude') == 'double'
        assert efel.cppcore.featuretype('AP_begin_indices') == 'int'

    def test_getgError(self):
        """cppcore: Testing getgError."""
        import efel
        efel.reset()
        assert efel.cppcore.getgError() == ''
        feature_values = list()
        efel.cppcore.getFeatureDouble('AP_amplitude', feature_values)
        assert "Feature T not found" in efel.cppcore.getgError()
        assert efel.cppcore.getgError() == ''

    def test_logging(self):  # pylint: disable=R0201
        """cppcore: Testing logging"""
        import efel
        with tempfile.TemporaryDirectory(prefix='efel_tests') as tempdir:
            efel.cppcore.Initialize(efel.getDependencyFileLocation(), tempdir)
            self.setup_data()
            with (Path(tempdir) / 'fllog.txt').open() as fd:
                contents = fd.read()
                assert 'Initializing' in contents
                # test vector working (if more than 10 elements, prints ...
                assert '...' in contents
            # re-call efel's Initialize with current dir
            # to remove pointer to tempdir.
            # this pointer was preventing the deletion of tempdir on windows.
            efel.cppcore.Initialize(efel.getDependencyFileLocation(), '.')

    @pytest.mark.parametrize(
        "feature_name", ['AP_amplitude', 'time_to_first_spike', 'peak_indices'])
    def test_caching(self, feature_name):
        """Test to make sure the caching mechanism works as intended.
        AP_amplitude: vector<double>
        time_to_first_spike: has a different implementation name
        peak_indices: vector<int>
        """
        with tempfile.TemporaryDirectory(prefix='efel_test_cache') as tempdir:
            import efel
            efel.cppcore.Initialize(efel.getDependencyFileLocation(), tempdir)
            self.setup_data()
            feature_values = list()
            efel.cppcore.getFeature(feature_name, feature_values)
            efel.cppcore.getFeature(feature_name, feature_values)
            efel.cppcore.getFeature(feature_name, feature_values)
            with (Path(tempdir) / 'fllog.txt').open() as fd:
                contents = fd.read()
            # re-call efel's Initialize with current dir to remove pointer to tempdir.
            # this pointer was preventing the deletion of tempdir on windows.
            efel.cppcore.Initialize(efel.getDependencyFileLocation(), '.')
        assert f"Calculated feature {feature_name}" in contents
        assert f"Reusing computed value of {feature_name}" in contents
        # make sure Calculated feature text occurs once
        assert contents.count(f"Calculated feature {feature_name}") == 1
        # make sure Reusing computed value of text occurs twice
        assert contents.count(f"Reusing computed value of {feature_name}") == 2

    def test_clear_function(self):
        """cppcore: Testing Clear function to reset state"""
        import efel
        self.setup_data()

        feature_values = list()
        efel.cppcore.getFeature('AP_amplitude', feature_values)
        assert len(feature_values) > 0  # Data should be present

        efel.cppcore.Clear()

        feature_values = list()
        return_value = efel.cppcore.getFeature('AP_amplitude', feature_values)
        assert return_value == -1  # Should return -1 since data is cleared


def test_efel_assertion_error():
    """Testing if C++ assertion error is propagated to Python correctly."""
    import efel
    efel.reset()
    trace = {
        "stim_start": [25],
        "stim_end": [75],
    }
    with pytest.raises(AssertionError):
        efel.get_feature_values([trace], ["__test_efel_assertion__"])
