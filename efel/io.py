from __future__ import annotations
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
from pathlib import Path
import neo
import numpy as np


def load_ascii_input(
    file_path: Path | str, delimiter: str = " "
) -> tuple[np.ndarray, np.ndarray]:
    """Loads electrophysiology data from an ASCII file.

    Returns: A tuple containing two numpy arrays, one for time and one for voltage.
    """
    file_path = Path(file_path)
    data = np.loadtxt(file_path, delimiter=delimiter)
    time, voltage = data[:, 0], data[:, 1]
    return time, voltage


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
    reader = neo.io.get_io(file_name)
    blocks = reader.read(**kwargs)

    stim_start, stim_end = extract_stim_times_from_neo_data(
        blocks, stim_start, stim_end)
    if stim_start is None or stim_end is None:
        raise ValueError(
            'No stim_start or stim_end found in epochs or events. '
            'Please specify "stim_start" and "stim_end" arguments.')

    # Convert data for eFEL
    efel_blocks = []
    for bl in blocks:
        efel_segments = []
        for seg in bl.segments:
            traces = []
            for sig in seg.analogsignals:
                trace = {
                    'T': sig.times.rescale('ms').magnitude,
                    'V': sig.rescale('mV').magnitude,
                    'stim_start': [stim_start],
                    'stim_end': [stim_end]
                }
                traces.append(trace)
            efel_segments.append(traces)
        efel_blocks.append(efel_segments)

    return efel_blocks
