"""Tests for extrafeats"""

import os

import numpy
import pytest

from efel.pyfeatures import extrafeats

testdata_dir = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), 'testdata', 'extrafeats'
)
waveforms_fpath = os.path.join(testdata_dir, 'mean_waveforms.dat')
waveforms = numpy.loadtxt(waveforms_fpath)
waveform = numpy.array([waveforms[0]])
sampling_freq = 10000


@pytest.mark.unit
def test_peak_to_valley():
    """extrafeats: Test peak_to_valley"""
    ptv = extrafeats.peak_to_valley(waveform, sampling_freq)
    assert len(ptv) == 1
    assert ptv[0] == pytest.approx(0.0013)


@pytest.mark.unit
def test_peak_trough_ratio():
    """extrafeats: Test peak_trough_ratio"""
    ptratio = extrafeats.peak_trough_ratio(waveform)
    assert len(ptratio) == 1
    assert ptratio[0] == pytest.approx(0.53804035)


@pytest.mark.unit
def test_halfwidth():
    """extrafeats: Test halfwidth"""
    ret = extrafeats.halfwidth(waveform, sampling_freq, True)
    assert len(ret) == 3

    hw = extrafeats.halfwidth(waveform, sampling_freq)
    assert len(hw) == 1
    assert hw[0] == pytest.approx(0.0015)


@pytest.mark.unit
def test_repolarization_slope():
    """extrafeats: Test repolarization_slope"""
    ret = extrafeats.repolarization_slope(
        waveform, sampling_freq, True
    )
    assert len(ret) == 2

    rslope = extrafeats.repolarization_slope(
        waveform, sampling_freq
    )
    assert len(rslope) == 1
    assert rslope[0] == pytest.approx(73.12572131)


@pytest.mark.unit
def test_recovery_slope():
    """extrafeats: Test recovery_slope"""
    window = 0.7
    rslope = extrafeats.recovery_slope(
        waveform, sampling_freq, window=window
    )
    assert len(rslope) == 1
    assert rslope[0] == pytest.approx(-3.63355521)


@pytest.mark.unit
def test_peak_image():
    """extrafeats: Test peak_image"""
    rel_peaks = extrafeats.peak_image(
        waveforms, sign="negative"
    )
    assert len(rel_peaks) == 209
    assert rel_peaks[0] == pytest.approx(0.06084468)

    rel_peaks = extrafeats.peak_image(
        waveforms, sign="positive"
    )
    assert len(rel_peaks) == 209
    assert rel_peaks[0] == pytest.approx(0.10850117)


@pytest.mark.unit
def test_relative_amplitude():
    """extrafeats: Test relative_amplitude"""
    rel_amp = extrafeats.relative_amplitude(
        waveforms, sign="negative"
    )
    assert len(rel_amp) == 209
    assert rel_amp[0] == pytest.approx(0.09513392)

    rel_amp = extrafeats.relative_amplitude(
        waveforms, sign="positive"
    )
    assert len(rel_amp) == 209
    assert rel_amp[0] == pytest.approx(0.2135929)


@pytest.mark.unit
def test_peak_time_diff():
    """extrafeats: Test peak_time_diff"""
    peak_t = extrafeats.peak_time_diff(
        waveforms, sampling_freq, sign="negative"
    )
    assert len(peak_t) == 209
    assert peak_t[0] == pytest.approx(0.0009)

    peak_t = extrafeats.peak_time_diff(
        waveforms, sampling_freq, sign="positive"
    )
    assert len(peak_t) == 209
    assert peak_t[0] == pytest.approx(0.0007)


@pytest.mark.unit
def test__get_trough_and_peak_idx():
    """extrafeats: Test _get_trough_and_peak_idx"""
    t_idx, p_idx = extrafeats._get_trough_and_peak_idx(
        waveform
    )
    assert t_idx == 102
    assert p_idx == 115


@pytest.mark.unit
def test_calculate_features():
    """extrafeats: Test calculate_features"""
    feats = extrafeats.calculate_features(
        waveforms, sampling_freq
    )
    for feature_name in extrafeats.all_1D_features:
        assert feature_name in feats
