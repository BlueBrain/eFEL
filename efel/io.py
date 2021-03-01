"""IO handler for eFEL"""

"""
Copyright (c) 2015, EPFL/Blue Brain Project

 This file is part of eFEL <https://github.com/BlueBrain/eFEL>

 This library is free software; you can redistribute it and/or modify it under
 the terms of the GNU Lesser General Public License version 3.0 as published
 by the Free Software Foundation.

 This library is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
 details.

 You should have received a copy of the GNU Lesser General Public License
 along with this library; if not, write to the Free Software Foundation, Inc.,
 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
"""

import os

# Python 2 has urlparse module, Python 3 has urllib.parse
try:
    import urlparse as up
    import urllib2 as ur
except ImportError:
    # pylint:disable=E0611, F0401
    import urllib.parse as up
    import urllib.request as ur
    # pylint:enable=E0611,F0401

import mimetypes


def windows_compatible(func):
    """Decorator allowing to use urlparse with windows."""

    def inner(fragment_url):
        """Change windows path part to url."""
        # if system is windows
        if os.name == "nt":
            fragment_url.replace("\\", "/")

        parsed_url = func(fragment_url)

        return parsed_url

    return inner


def load_fragment(fragment_url, mime_type=None):
    """Load fragment

    Load a fragment (e.g. time series data) from a given URL
    """
    parsed_url = windows_compatible(up.urlparse)(fragment_url)

    scheme = parsed_url.scheme
    server_loc = parsed_url.netloc
    path = parsed_url.path
    fragment_string = parsed_url.fragment

    # reform path for windows files
    if scheme == "file" and os.name == "nt":
        path = ur.url2pathname(r"\\" + server_loc + path)

    if mime_type is None:
        mimetypes.init()
        mime_type, _ = mimetypes.guess_type(path)
        if mime_type is None:
            raise TypeError(
                'load_fragment: impossible to guess MIME type from url, '
                'please specify the type manually as argument: %s' % path)

    if scheme == 'file' and os.name == 'nt':
        file_handle = open(path, 'r')
    elif scheme == 'file':
        file_handle = open(os.path.join(server_loc, path), 'r')

    if 'text/' in mime_type:
        import numpy

        if fragment_string == '':
            cols = None
        else:
            import re

            match = re.match("col=([0-9]+)", fragment_string)
            if match is None or len(match.groups()) != 1:
                raise TypeError(
                    "load_fragment: don't understand url fragment %s" %
                    fragment_string)
            else:
                cols = int(match.groups()[0]) - 1

        # Unfortunately we need this if statement
        # Setting usecols to None throws an error in the loadtxt call
        if cols is not None:
            fragment_content = numpy.loadtxt(file_handle, usecols=[cols])
        else:
            fragment_content = numpy.loadtxt(file_handle)

        return fragment_content
    else:
        raise TypeError('load_fragment: unknown mime type %s' % mime_type)


def extract_stim_times_from_neo_data(blocks, stim_start, stim_end):
    """
        Seeks for the stim_start and stim_end parameters inside the Neo data.

        Parameters
        ==========
        blocks : Neo object blocks
        stim_start : numerical value (ms) or None
        stim_end : numerical value (ms) or None

        Epoch.name should be one of "stim", "stimulus", "stimulation",
        "current_injection"
        First Event.name should be "stim_start", "stimulus_start",
        "stimulation_start", "current_injection_start"
        Second Event.name should be one of "stim_end",
        "stimulus_end", "stimulation_end", "current_injection_end"

        Returned objects
        ====================
        stim_start : numerical value (ms) or None
        stim_end : numerical value (ms) or None

    """
    # this part code aims to find informations about stimulations, if
    # stimulation time are None.
    if stim_start is None and stim_end is None:
        for bl in blocks:
            for seg in bl.segments:
                event_start_rescaled = False
                event_end_rescaled = False
                epoch_rescaled = False
                for epoch in seg.epochs:
                    if epoch.name in (
                            "stim",
                            "stimulus",
                            "stimulation",
                            "current_injection"):
                        if stim_start is None:
                            epoch = epoch.rescale('ms').magnitude
                            stim_start = epoch[0]
                            stim_end = epoch[-1]
                            epoch_rescaled = True
                        else:
                            raise ValueError(
                                'It seems that there are two epochs related '
                                'to stimulation, the program does not know '
                                'which one to chose')

                for event in seg.events:
                    # Test event stim_start
                    if event.name in (
                            "stim_start",
                            "stimulus_start",
                            "stimulation_start",
                            "current_injection_start"):
                        # tests if not already found once
                        if event_start_rescaled:
                            raise ValueError(
                                'It seems that stimulation start time is '
                                'defined more than one event.'
                                ' The program does not know which one '
                                'to choose')
                        if epoch_rescaled:
                            raise ValueError(
                                'It seems that stim_start is defined in both '
                                'epoch and an event.' +
                                ' The program does not know which one '
                                'to choose')

                        if stim_start is None:
                            event = event.rescale('ms').magnitude

                            try:
                                stim_start = event[0]
                            except BaseException:
                                stim_start = event

                            event_start_rescaled = True

                    # Test event stim_end
                    elif event.name in (
                        "stim_end",
                        "stimulus_end",
                        "stimulation_end",
                        "current_injection_end",
                    ):
                        # tests if not already found once
                        if event_end_rescaled:
                            raise ValueError(
                                'It seems that stimulation end time is '
                                'defined more than one event.'
                                ' The program does not know which one '
                                'to choose')
                        if epoch_rescaled:
                            raise ValueError(
                                'It seems that stim_end is defined in '
                                'both epoch and an event.'
                                ' The program does not know which one '
                                'to choose')

                        if stim_end is None:
                            event = event.rescale('ms').magnitude

                            try:
                                stim_end = event[-1]
                            except BaseException:
                                stim_end = event

                            event_end_rescaled = True

    return stim_start, stim_end


def load_neo_file(file_name, stim_start=None, stim_end=None, **kwargs):
    """
        Use neo to load a data file and convert it to be readable by eFEL.

        Parameters
        ==========
        file_name : string
                    path to the location of a Dependency file
        stim_start : numerical value (ms)
                    Optional if there is an Epoch or two Events in the file
        stim_end : numerical value (ms)
                Optional if there is an Epoch or two Events in the file
        kwargs : keyword arguments to be passed to the read() method of the
                Neo IO class

        Epoch.name should be one of "stim", "stimulus", "stimulation",
        "current_injection"
        First Event.name should be "stim_start", "stimulus_start",
        "stimulation_start", "current_injection_start"
        Second Event.name should be one of "stim_end", "stimulus_end",
        "stimulation_end", "current_injection_end"

        The returned object is presented like this :
            returned object : [Segments_1, Segments_2, ..., Segments_n]
            Segments_1 = [Traces_1, Traces_2, ..., Traces_n]
    """

    import neo

    reader = neo.io.get_io(file_name)
    blocks = reader.read(**kwargs)

    stim_start, stim_end = extract_stim_times_from_neo_data(
        blocks, stim_start, stim_end)
    if stim_start is None or stim_end is None:
        raise ValueError(
            'No stim_start or stim_end has been found inside epochs or events.'
            ' You can directly specify their value as argument "stim_start"'
            ' and "stim_end"')

    # this part of the code transforms the data format.
    efel_blocks = []
    for bl in blocks:
        efel_segments = []
        for seg in bl.segments:
            traces = []
            count_traces = 0
            analogsignals = seg.analogsignals

            for sig in analogsignals:
                traces.append({})
                traces[count_traces]['T'] = sig.times.rescale('ms').magnitude
                traces[count_traces]['V'] = sig.rescale('mV').magnitude
                traces[count_traces]['stim_start'] = [stim_start]
                traces[count_traces]['stim_end'] = [stim_end]
                count_traces += 1

            efel_segments.append(traces)
        efel_blocks.append(efel_segments)

    return efel_blocks
