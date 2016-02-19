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
import nose.tools as nt

import numpy as np

testdata_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            'testdata/')


def test_import():
    """cppcore: Testing import of cppcore"""
    import efel.cppcore


class TestCppcore(object):

    def setup(self):
        """Setup"""
        import efel
        efel.cppcore.Initialize(efel.getDependencyFileLocation(), "log")

    def setup_data(self):
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


    def test_getFeatureNames(self):
        """cppcore: Testing getting all feature names"""
        import efel.cppcore
        feature_names = []
        efel.cppcore.getFeatureNames(feature_names)

        import json
        test_data_path = os.path.join(testdata_dir, '../featurenames.json')
        with open(test_data_path, 'r') as featurenames_json:
            expected_featurenames = json.load(featurenames_json)
        nt.assert_equal(feature_names, expected_featurenames)

    def test_getFeatureDouble_failure(self):
        """cppcore: Testing failure exit code in getFeatureDouble"""
        import efel.cppcore
        feature_values = list()
        return_value = efel.cppcore.getFeatureDouble(
            "AP_amplitude",
            feature_values)
        nt.assert_equal(return_value, -1)

    @nt.raises(TypeError)
    def test_getFeatureDouble_wrong_type(self):
        import efel.cppcore
        efel.cppcore.getFeatureDouble("AP_fall_indices", list())

    @nt.raises(TypeError)
    def test_getFeatureInt_wrong_type(self):
        import efel.cppcore
        efel.cppcore.getFeatureInt("AP_amplitude", list())

    def test_getDistance(self):
        """cppcore: Testing getDistance()"""
        import efel.cppcore

        self.setup_data()
        nt.assert_almost_equal(
            3.09045815935,
            efel.cppcore.getDistance(
                'AP_amplitude',
                50.0,
                10.0))

    def test_getFeature(self):
        import efel.cppcore
        self.setup_data()

        #get double feature
        feature_values = list()
        efel.cppcore.getFeature('AP_amplitude', feature_values)
        nt.ok_(isinstance(feature_values[0], float))
        nt.eq_(5, len(feature_values))
        nt.ok_(np.allclose([80.45724099440199, 80.46320199354948, 80.73300299176428,
                            80.9965359926715, 81.87292599493423],
                           feature_values))

        #get int feature
        feature_values = list()
        efel.cppcore.getFeature('AP_fall_indices', feature_values)
        nt.ok_(isinstance(feature_values[0], int))
        nt.eq_(5, len(feature_values))
        nt.eq_([5665, 6066, 6537, 7170, 8275], feature_values)

    def test_getFeature_failure(self):
        """cppcore: Testing failure exit code in getFeature"""
        import efel.cppcore
        feature_values = list()
        return_value = efel.cppcore.getFeature("AP_amplitude", feature_values)
        nt.assert_equal(return_value, -1)

    @nt.raises(TypeError)
    def test_getFeature_non_existant(self):
        """cppcore: Testing failure exit code in getFeature"""
        import efel.cppcore
        efel.cppcore.getFeature("does_not_exist", list())

    def test_registerFeature(self):
        import efel.cppcore

        def getTwiceAP_Amp(output):
            feature_values = list()
            efel.cppcore.getFeature('AP_amplitude', feature_values)
            output.extend(2*x for x in feature_values)
            return len(feature_values)

        self.setup_data()
        efel.cppcore.registerFeature('getTwiceAP_Amp', getTwiceAP_Amp)
        features = list()
        ret = efel.cppcore.getFeature('getTwiceAP_Amp', features)
        nt.ok_(isinstance(features[0], float))
        nt.eq_(5, len(features))
        nt.ok_(np.allclose(2 * np.array([80.45724099440199, 80.46320199354948, 80.73300299176428,
                                         80.9965359926715, 81.87292599493423]),
                           features))
