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
import csv
import json
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


def extract_stim_times_from_neo_data(blocks, stim_start, stim_end) -> tuple:
    """
    Seeks for the stim_start and stim_end parameters inside the Neo data.

    Args:
        blocks (Neo object blocks): Description of what blocks represents.
        stim_start (numerical value or None): Start time of the stimulation in
            milliseconds. If not available, None should be used.
        stim_end (numerical value or None): End time of the stimulation in
            milliseconds. If not available, None should be used.

    Returns:
        tuple: A tuple containing:
            - stim_start (numerical value or None): Start time of the stimulation
              in milliseconds.
            - stim_end (numerical value or None): End time of the stimulation in
              milliseconds.

    Notes:
        - Epoch.name should be one of "stim", "stimulus", "stimulation",
          "current_injection".
        - First Event.name should be "stim_start", "stimulus_start",
          "stimulation_start", "current_injection_start".
        - Second Event.name should be one of "stim_end", "stimulus_end",
          "stimulation_end", "current_injection_end".

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


def load_neo_file(file_name: str, stim_start=None, stim_end=None, **kwargs) -> list:
    """
    Loads a data file using neo and converts it for eFEL readability.

    Args:
        file_name (string): Path to the Dependency file location.
        stim_start (numerical value, optional): Start time in ms. Optional if an Epoch
            or two Events are in the file.
        stim_end (numerical value, optional): End time in ms. Optional if an Epoch
            or two Events are in the file.
        **kwargs: Additional arguments for the read() method of Neo IO class.

    Returns:
        list of Segments: Segments containing traces, formatted as
            [Segments_1, Segments_2, ..., Segments_n], where each Segments_i is
            [Traces_1, Traces_2, ..., Traces_n].

    Notes:
        - Epoch.name should be "stim", "stimulus", "stimulation", "current_injection".
        - First Event.name: "stim_start", "stimulus_start", "stimulation_start",
          "current_injection_start".
        - Second Event.name: "stim_end", "stimulus_end", "stimulation_end",
          "current_injection_end".
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
                if 'V' in sig.units.dimensionality.string:
                    trace = {
                        'T': (sig.times - sig.times[0]).rescale('ms').magnitude,
                        'V': sig.rescale('mV').magnitude,
                        'stim_start': [stim_start],
                        'stim_end': [stim_end]
                    }
                    traces.append(trace)
            efel_segments.append(traces)
        efel_blocks.append(efel_segments)

    return efel_blocks


def save_feature_to_json(feature_values, filename):
    """Save feature values as a JSON file."""
    class NumpyEncoder(json.JSONEncoder):
        def default(self, o):
            if isinstance(o, np.integer):
                return int(o)
            elif isinstance(o, np.floating):
                return float(o)
            elif isinstance(o, np.ndarray):
                return o.tolist()
            else:
                return super().default(o)

    with open(filename, 'w') as f:
        json.dump(feature_values, f, cls=NumpyEncoder)


def save_feature_to_csv(feature_values, filename):
    """Save feature values as a CSV file."""
    with open(filename, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(feature_values.keys())
        max_list_length = max(len(value) if isinstance(value, (list, np.ndarray))
                              else 1 for value in feature_values.values())
        for i in range(max_list_length):
            row_data = [value[i] if isinstance(value, (list, np.ndarray))
                        and i < len(value) else ''
                        for value in feature_values.values()]
            writer.writerow(row_data)
