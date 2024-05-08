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
nwb1_filename = testdata_dir / 'JY180308_A_1.nwb'


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


def test_load_neo_file_nwb():
    """Test loading of data from an NWB file."""
    import efel
    import neo

    efel_blocks = efel.io.load_neo_file(nwb1_filename, 250, 600)
    efel_trace = efel_blocks[0][0][0]

    reader = neo.io.NWBIO(filename=nwb1_filename)
    neo_blocks = reader.read()
    neo_trace = neo_blocks[0].segments[0].analogsignals[0]
    time_shifted = neo_trace.times[-1] - neo_trace.times[0]
    rescaled_time_shifted = time_shifted.rescale('ms').magnitude

    assert neo_trace[0].rescale('mV').magnitude == efel_trace['V'][0]
    assert rescaled_time_shifted == efel_trace['T'][-1]

    for efel_seg in efel_blocks:
        for traces in efel_seg:
            for trace in traces:
                assert trace['stim_start'] == [250]
                assert trace['stim_end'] == [600]


@pytest.fixture(params=['tmp.json', 'tmp.csv'])
def filename(request):
    yield request.param
    if os.path.exists(request.param):
        os.remove(request.param)


def load_data(filename):
    import csv
    import json
    if filename.endswith('.json'):
        with open(filename, 'r') as f:
            return json.load(f)
    elif filename.endswith('.csv'):
        loaded_data = {}
        with open(filename, 'r') as f:
            reader = csv.reader(f)
            header = next(reader)
            for col_name in header:
                loaded_data[col_name] = []

            for row in reader:
                for col_name, value in zip(header, row):
                    if value != '':
                        loaded_data[col_name].append(float(value))

            for loaded_key in loaded_data.keys():
                if not loaded_data[loaded_key]:
                    loaded_data[loaded_key] = None
        return loaded_data


@pytest.mark.parametrize("index", [2, 3, 5])
def test_save_feature(filename, index):
    """Test saving of the features to file."""
    import efel
    blocks = efel.io.load_neo_file(nwb1_filename, 250, 600)
    trace = blocks[0][0][index]
    features = ['peak_time', 'AP_height', 'peak_indices', 'spikes_per_burst']
    feature_values = efel.get_feature_values([trace], features)[0]

    if filename.endswith('.json'):
        efel.io.save_feature_to_json(feature_values, filename)
    elif filename.endswith('.csv'):
        efel.io.save_feature_to_csv(feature_values, filename)

    assert os.path.exists(filename)
    loaded_data = load_data(filename)
    for key in feature_values.keys():
        if feature_values[key] is not None and loaded_data[key]:
            assert np.allclose(feature_values[key],
                               loaded_data[key],
                               rtol=1e-05,
                               atol=1e-08)
        else:
            assert feature_values[key] == loaded_data[key]
