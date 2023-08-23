"""Test eFEL io module"""

# pylint: disable=F0401

import os

import pytest

testdata_dir = os.path.join(
    os.path.dirname(
        os.path.abspath(__file__)),
    'testdata')

neo_test_files_dir = os.path.join(
    os.path.dirname(
        os.path.abspath(__file__)),
    'neo_test_files')

meanfrequency1_filename = os.path.join(testdata_dir,
                                       'basic',
                                       'mean_frequency_1.txt')
meanfrequency1_url = 'file://%s' % meanfrequency1_filename


def test_import():
    """io: Testing import"""

    # pylint: disable=W0611
    import efel.io  # NOQA
    # pylint: enable=W0611


def test_import_without_urlparse():
    """io: Testing import without urlparse"""

    # The only purpose of this test is to get the code coverage to 100% :-)

    import sys
    del sys.modules['efel.io']

    python_version = sys.version_info[0]

    if python_version < 3:
        import __builtin__
    else:
        import builtins as __builtin__
    realimport = __builtin__.__import__

    def myimport(name, *args):  # global_s, local, fromlist, level):
        """Override import"""
        if name == 'urlparse':
            raise ImportError
        return realimport(name, *args)  # global_s, local, fromlist, level)
    __builtin__.__import__ = myimport

    try:
        import urllib.parse  # NOQA
        urllibparse_import_fails = False
    except ImportError:
        urllibparse_import_fails = True

    if urllibparse_import_fails:
        pytest.raises(ImportError, __builtin__.__import__, 'efel.io')
    else:
        import efel.io  # NOQA

    __builtin__.__import__ = realimport


def test_load_fragment_column_txt1():
    """io: Test loading of one column from txt file"""

    import efel
    import numpy

    time_io = efel.io.load_fragment('%s#col=1' % meanfrequency1_url)

    time_numpy = numpy.loadtxt(meanfrequency1_filename, usecols=[0])

    numpy.testing.assert_array_equal(time_io, time_numpy)


def test_load_fragment_strange_mimetype():
    """io: Test loading file with unresolvable mime type"""

    import efel

    pytest.raises(
        TypeError,
        efel.io.load_fragment, 'file://strange.mimetype')


def test_load_fragment_wrong_fragment_format():
    """io: Test loading file wrong fragment format"""

    import efel

    pytest.raises(
        TypeError,
        efel.io.load_fragment,
        '%s#co=1' %
        meanfrequency1_url)


def test_load_fragment_wrong_mimetype():
    """io: Test loading fragment wrong mime type"""

    import efel

    pytest.raises(
        TypeError,
        efel.io.load_fragment,
        '%s#col=1' % meanfrequency1_url, mime_type='application/json')


def test_load_fragment_allcolumns():
    """io: Test loading fragment without specifying columns"""

    import efel
    import numpy

    time_io = efel.io.load_fragment('%s' % meanfrequency1_url)

    time_numpy = numpy.loadtxt(meanfrequency1_filename)

    numpy.testing.assert_array_equal(time_io, time_numpy)


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
