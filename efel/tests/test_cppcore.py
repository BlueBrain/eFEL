"""Test all features on an example trace"""

"""
Copyright (c) 2015, Blue Brain Project/EPFL and AUTHORS.txt

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
  * Redistributions of source code must retain the above copyright
    notice, this list of conditions and the following disclaimer.
  * Redistributions in binary form must reproduce the above copyright
    notice, this list of conditions and the following disclaimer in the
    documentation and/or other materials provided with the distribution.
  * Neither the name of the <organization> nor the
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
