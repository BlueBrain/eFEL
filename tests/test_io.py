"""Test eFEL io module"""

import os
from pathlib import Path
import numpy as np

import pytest

from efel.io import load_ascii_input


testdata_dir = Path(__file__).parent / 'testdata'
neo_test_files_dir = Path(__file__).parent / 'neo_test_files'
meanfrequency1_filename = testdata_dir / 'basic' / 'mean_frequency_1.txt'
meanfrequency1_url = str(meanfrequency1_filename)


def test_load_ascii_input():
    """Test loading of data from an ASCII file and splitting into time and voltage."""
    time, voltage = load_ascii_input(meanfrequency1_url)

    # Load data using numpy for comparison
    expected_data = np.loadtxt(meanfrequency1_url, delimiter=' ')
    expected_time, expected_voltage = expected_data[:, 0], expected_data[:, 1]

    # Assert that the arrays are equal
    np.testing.assert_array_equal(time, expected_time)
    np.testing.assert_array_equal(voltage, expected_voltage)


def test_load_neo_file_stim_time_arg():
    import efel
    file_name = os.path.join(neo_test_files_dir, "neo_test_file_no_times.mat")

    # test load_neo_file without stim time
    pytest.raises(ValueError, efel.io.load_neo_file, file_name)
    # test load_neo_file with stim time arguments
    result = efel.io.load_neo_file(file_name, stim_start=0, stim_end=20)
    # test load_neo_file with stim time incomplete arguments
    pytest.raises(
        ValueError,
        efel.io.load_neo_file,
        file_name,
        stim_start=0)
    assert (
        all(
            result[0][0][0]['T'] == [
                0., 1., 2., 3., 4., 5., 6., 7., 8., 9.]))
    assert (
        all(
            result[0][0][0]['V'] == [
                [0], [1], [2], [3], [4], [5], [6], [7], [8], [9]]))
    assert result[0][0][0]['stim_start'] == [0.0]
    assert result[0][0][0]['stim_end'] == [20.0]


def test_extract_stim_times_from_neo_data_two_epochs():
    import neo
    import efel
    import quantities as pq

    bl = neo.core.Block()
    seg = neo.core.Segment()
    data = range(10)
    rate = 1000 * pq.Hz
    signal = neo.core.AnalogSignal(data, sampling_rate=rate, units="mV")
    seg.analogsignals.append(signal)
    seg.epochs.append(
        neo.core.Epoch(
            times=pq.Quantity([0.0, 20.0],
                              units=pq.ms),
            name="stim"))
    seg.epochs.append(
        neo.core.Epoch(
            times=pq.Quantity([0.0, 20.0],
                              units=pq.ms),
            name="stim"))
    bl.segments.append(seg)

    pytest.raises(
        ValueError,
        efel.io.extract_stim_times_from_neo_data,
        [bl],
        None,
        None)


def test_extract_stim_times_from_neo_data_two_events_start():
    import neo
    import efel
    import quantities as pq

    bl = neo.core.Block()
    seg = neo.core.Segment()
    data = range(10)
    rate = 1000 * pq.Hz
    signal = neo.core.AnalogSignal(data, sampling_rate=rate, units="mV")
    seg.analogsignals.append(signal)
    seg.events.append(
        neo.core.Event(
            times=[0.0] * pq.ms,
            units=pq.ms,
            name="stim_start"))
    seg.events.append(
        neo.core.Event(
            times=[20.0] * pq.ms,
            units=pq.ms,
            name="stim_start"))
    bl.segments.append(seg)

    pytest.raises(
        ValueError,
        efel.io.extract_stim_times_from_neo_data,
        [bl],
        None,
        None)


def test_extract_stim_times_from_neo_data_two_events_end():
    import neo
    import efel
    import quantities as pq

    bl = neo.core.Block()
    seg = neo.core.Segment()
    data = range(10)
    rate = 1000 * pq.Hz
    signal = neo.core.AnalogSignal(data, sampling_rate=rate, units="mV")
    seg.analogsignals.append(signal)
    seg.events.append(
        neo.core.Event(
            times=[0.0] * pq.ms,
            units=pq.ms,
            name="stim_end"))
    seg.events.append(
        neo.core.Event(
            times=[20.0] * pq.ms,
            units=pq.ms,
            name="stim_end"))
    bl.segments.append(seg)

    pytest.raises(
        ValueError,
        efel.io.extract_stim_times_from_neo_data,
        [bl],
        None,
        None)


def test_extract_stim_times_from_neo_data_start_in_epoch_event():
    import neo
    import efel
    import quantities as pq

    bl = neo.core.Block()
    seg = neo.core.Segment()
    data = range(10)
    rate = 1000 * pq.Hz
    signal = neo.core.AnalogSignal(data, sampling_rate=rate, units="mV")
    seg.analogsignals.append(signal)
    seg.epochs.append(
        neo.core.Epoch(
            times=pq.Quantity([0.0, 20.0],
                              units=pq.ms),
            name="stim"))
    seg.events.append(
        neo.core.Event(
            times=[0.0] * pq.ms,
            units=pq.ms,
            name="stim_start"))
    bl.segments.append(seg)

    pytest.raises(
        ValueError,
        efel.io.extract_stim_times_from_neo_data,
        [bl],
        None,
        None)


def test_extract_stim_times_from_neo_data_end_in_epoch_event():
    import neo
    import efel
    import quantities as pq

    bl = neo.core.Block()
    seg = neo.core.Segment()
    data = range(10)
    rate = 1000 * pq.Hz
    signal = neo.core.AnalogSignal(data, sampling_rate=rate, units="mV")
    seg.analogsignals.append(signal)
    seg.epochs.append(
        neo.core.Epoch(
            times=pq.Quantity([0.0, 20.0],
                              units=pq.ms),
            name="stim"))
    seg.events.append(
        neo.core.Event(
            times=[0.0] * pq.ms,
            units=pq.ms,
            name="stim_end"))
    bl.segments.append(seg)

    pytest.raises(
        ValueError,
        efel.io.extract_stim_times_from_neo_data,
        [bl],
        None,
        None)


def test_extract_stim_times_from_neo_data_event_time_list():
    import neo
    import efel
    import quantities as pq

    bl = neo.core.Block()
    seg = neo.core.Segment()
    data = range(10)
    rate = 1000 * pq.Hz
    signal = neo.core.AnalogSignal(data, sampling_rate=rate, units="mV")
    seg.analogsignals.append(signal)
    seg.events.append(
        neo.core.Event(
            times=[
                0.0,
                1.0] * pq.ms,
            units=pq.ms,
            name="stim_start"))
    seg.events.append(
        neo.core.Event(
            times=[
                10.0,
                20.0] * pq.ms,
            units=pq.ms,
            name="stim_end"))
    bl.segments.append(seg)

    assert (
        (0.0, 20.0) == efel.io.extract_stim_times_from_neo_data(
            [bl], None, None))


def test_load_neo_file_stim_time_epoch():
    import efel
    file_name = os.path.join(
        neo_test_files_dir,
        "neo_test_file_epoch_times.mat")

    result = efel.io.load_neo_file(file_name)
    assert result[0][0][0]['stim_start'] == [0.0]
    assert result[0][0][0]['stim_end'] == [20.0]


def test_load_neo_file_stim_time_events():
    import efel
    file_name = os.path.join(
        neo_test_files_dir,
        "neo_test_file_events_time.mat")

    result = efel.io.load_neo_file(file_name)
    assert result[0][0][0]['stim_start'] == [0.0]
    assert result[0][0][0]['stim_end'] == [20.0]


def test_load_neo_file_stim_time_events_incomplete():
    import efel
    file_name = os.path.join(neo_test_files_dir,
                             "neo_test_file_events_time_incomplete.mat")

    pytest.raises(ValueError, efel.io.load_neo_file, file_name)
