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

    def test_getDistance(self):
        """cppcore: Testing getDistance()"""
        import efel.cppcore
        import numpy

        stim_start = 500.0
        stim_end = 900.0

        test_data_path = os.path.join(testdata_dir, 'basic/mean_frequency_1.txt')
        data = numpy.loadtxt(test_data_path)

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
        feature_values = list()
        efel.cppcore.getFeatureDouble('AP_amplitude', feature_values)

        nt.assert_almost_equal(
            3.09045815935,
            efel.cppcore.getDistance(
                'AP_amplitude',
                50.0,
                10.0))
