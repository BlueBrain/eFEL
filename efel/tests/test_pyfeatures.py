"""Tests for the python features of eFEL"""

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
import numpy
import nose.tools as nt

# from nose.plugins.attrib import attr

import efel

testdata_dir = os.path.join(
    os.path.dirname(
        os.path.abspath(__file__)),
    'testdata')

traces_data = {
    'mean_frequency1': {
        'url': 'file://%s' % os.path.join(
            os.path.abspath(testdata_dir),
            'basic',
            'mean_frequency_1.txt'),
        't_col': 1,
        'v_col': 2,
        'stim_start': 500,
        'stim_end': 900},
    'init_burst1': {
        'url': 'file://%s' % os.path.join(
            os.path.abspath(testdata_dir),
            'basic',
            'init_burst1.txt'),
        't_col': 1,
        'v_col': 2,
        'stim_start': 250,
        'stim_end': 1600},
    'init_burst2': {
        'url': 'file://%s' % os.path.join(
            os.path.abspath(testdata_dir),
            'basic',
            'init_burst2.txt'),
        't_col': 1,
        'v_col': 2,
        'stim_start': 250,
        'stim_end': 1600},
    'init_burst3': {
        'url': 'file://%s' % os.path.join(
            os.path.abspath(testdata_dir),
            'basic',
            'init_burst3.txt'),
        't_col': 1,
        'v_col': 2,
        'stim_start': 250,
        'stim_end': 1600}}


def _load_trace(trace_name):
    """Load trace with a certain name"""

    trace_data = traces_data[trace_name]

    url = trace_data['url']

    trace = {
        'T': efel.io.load_fragment(
            '%s#col=%d' %
            (url, trace_data['t_col'])),
        'V': efel.io.load_fragment(
            '%s#col=%d' %
            (url, trace_data['v_col'])),
        'stim_start': [trace_data['stim_start']],
        'stim_end': [trace_data['stim_end']],
    }

    return trace


def _test_expected_value(feature_name, expected_values):
    """Test expected values for feature"""

    for trace_name, expected_value in expected_values.items():
        trace = _load_trace(trace_name)
        feature_values = efel.getFeatureValues([trace], [feature_name])

        nt.assert_almost_equal(feature_values[0][feature_name], expected_value)


def test_initburst_sahp():
    """pyfeatures: Test initburst_sahp feature"""

    feature_name = 'initburst_sahp'
    expected_values = {
        'mean_frequency1': None, 'init_burst1': [-69.19999695],
        'init_burst2': [-74.34375],
        'init_burst3': [-76.72435972]}

    _test_expected_value(feature_name, expected_values)


def test_initburst_sahp_vb():
    """pyfeatures: Test initburst_sahp_vb feature"""

    feature_name = 'initburst_sahp_vb'
    expected_values = {
        'mean_frequency1': None, 'init_burst1': [13.80537756],
        'init_burst2': [11.43360025],
        'init_burst3': [6.08649806]}

    _test_expected_value(feature_name, expected_values)


def test_initburst_sahp_ssse():
    """pyfeatures: Test initburst_sahp_ssse feature"""

    feature_name = 'initburst_sahp_ssse'
    expected_values = {
        'mean_frequency1': None, 'init_burst1': [-11.33152931],
        'init_burst2': [2.02011574],
        'init_burst3': [-4.42756346]}

    _test_expected_value(feature_name, expected_values)


def test_ISIs():
    """pyfeatures: Test ISIs feature"""

    mf1_trace = _load_trace('mean_frequency1')

    feature_values = efel.getFeatureValues([mf1_trace], ['ISI_values', 'ISIs'])

    nt.assert_true(numpy.allclose(
        feature_values[0]['ISIs'][1:],
        feature_values[0]['ISI_values']))
