"""Python implementation of features"""

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


import numpy
import efel.cppcore


all_pyfeatures = [
    'voltage',
    'time',
    'ISIs',
    'initburst_sahp',
    'initburst_sahp_vb',
    'initburst_sahp_ssse']


def voltage():
    """Get voltage trace"""
    return _get_cpp_feature("voltage")


def time():
    """Get time trace"""
    return _get_cpp_feature("time")


def ISIs():
    """Get all ISIs"""

    peak_times = _get_cpp_feature("peak_time")

    return numpy.diff(peak_times)


def initburst_sahp_vb():
    """SlowAHP voltage from voltage base after initial burst"""

    # Required cpp features
    initburst_sahp_value = initburst_sahp()
    voltage_base = _get_cpp_feature("voltage_base")

    if initburst_sahp_value is None or voltage_base is None or \
            len(initburst_sahp_value) != 1 or len(voltage_base) != 1:
        return None
    else:
        return numpy.array([initburst_sahp_value[0] - voltage_base[0]])


def initburst_sahp_ssse():
    """SlowAHP voltage from steady_state_voltage_stimend after initial burst"""

    # Required cpp features
    initburst_sahp_value = initburst_sahp()
    ssse = _get_cpp_feature("steady_state_voltage_stimend")

    if initburst_sahp_value is None or ssse is None or \
            len(initburst_sahp_value) != 1 or len(ssse) != 1:
        return None
    else:
        return numpy.array([initburst_sahp_value[0] - ssse[0]])


def initburst_sahp():
    """SlowAHP voltage after initial burst"""

    # Required cpp features
    voltage = _get_cpp_feature("voltage")
    time = _get_cpp_feature("time")
    peak_times = _get_cpp_feature("peak_time")

    # Required python features
    all_isis = ISIs()

    # Required trace data
    stim_end = _get_cpp_data("stim_end")

    # Required settings
    initburst_freq_thresh = _get_cpp_data("initburst_freq_threshold")
    initburst_sahp_start = _get_cpp_data("initburst_sahp_start")
    initburst_sahp_end = _get_cpp_data("initburst_sahp_end")

    last_isi = None

    # Loop over ISIs until frequency higher than initburst_freq_threshold
    for isi_counter, isi in enumerate(all_isis):
        # Convert to Hz
        freq = 1000.0 / isi
        if freq < initburst_freq_thresh:
            # Threshold reached
            break
        else:
            # Add isi to initburst
            last_isi = isi_counter

    if last_isi is None:
        # No initburst found
        return None
    else:
        # Get index of second peak of last ISI
        last_peak = last_isi + 1

    # Get time of last peak
    last_peak_time = peak_times[last_peak]

    # Determine start of sahp interval
    sahp_interval_start = min(
        last_peak_time +
        initburst_sahp_start,
        stim_end)

    # Get next peak, we wont search beyond that
    next_peak = last_peak + 1

    # Determine end of sahp interval
    # Add initburst_slow_ahp_max to last peak time
    # If next peak or stim_end is earlier, use these
    # If no next peak, use stim end
    if next_peak < len(peak_times):
        next_peak_time = peak_times[next_peak]

        sahp_interval_end = min(
            last_peak_time + initburst_sahp_end, next_peak_time, stim_end)
    else:
        sahp_interval_end = min(
            last_peak_time + initburst_sahp_end, stim_end)

    if sahp_interval_end <= sahp_interval_start:
        return None
    else:
        sahp_interval = voltage[numpy.where(
            (time <= sahp_interval_end) &
            (time >= sahp_interval_start))]

        if len(sahp_interval) > 0:
            min_volt_index = numpy.argmin(sahp_interval)
        else:
            return None

        slow_ahp = sahp_interval[min_volt_index]

        return numpy.array([slow_ahp])


def _get_cpp_feature(feature_name):
    """Get cpp feature"""
    cppcoreFeatureValues = list()
    exitCode = efel.cppcore.getFeature(feature_name, cppcoreFeatureValues)

    if exitCode < 0:
        return None
    else:
        return numpy.array(cppcoreFeatureValues)


def _get_cpp_data(data_name):
    """Get cpp data value"""

    return efel.cppcore.getMapDoubleData(data_name)[0]
