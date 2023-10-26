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
from numpy.fft import *


all_pyfeatures = [
    'voltage',
    'time',
    'current',
    'ISIs',
    'initburst_sahp',
    'initburst_sahp_vb',
    'initburst_sahp_ssse',
    'depol_block',
    'depol_block_bool',
    'spikes_per_burst',
    'spikes_per_burst_diff',
    'spikes_in_burst1_burst2_diff',
    'spikes_in_burst1_burstlast_diff',
    'impedance',
]


def voltage():
    """Get voltage trace"""
    return _get_cpp_feature("voltage")


def time():
    """Get time trace"""
    return _get_cpp_feature("time")


def impedance():
    from scipy.ndimage.filters import gaussian_filter1d

    dt = _get_cpp_data("interp_step")
    Z_max_freq = _get_cpp_data("impedance_max_freq")
    voltage_trace = voltage()
    holding_voltage = _get_cpp_feature("voltage_base")
    normalized_voltage = voltage_trace - holding_voltage
    current_trace = current()
    if current_trace is not None:
        holding_current = _get_cpp_feature("current_base")
        normalized_current = current_trace - holding_current
        spike_count = _get_cpp_feature("Spikecount")
        if spike_count < 1:  # if there is no spikes in ZAP
            fft_volt = numpy.fft.fft(normalized_voltage)
            fft_cur = numpy.fft.fft(normalized_current)
            if any(fft_cur) == 0:
                return None
            # convert dt from ms to s to have freq in Hz
            freq = numpy.fft.fftfreq(len(normalized_voltage), d=dt / 1000.)
            Z = fft_volt / fft_cur
            norm_Z = abs(Z) / max(abs(Z))
            select_idxs = numpy.swapaxes(
                numpy.argwhere((freq > 0) & (freq <= Z_max_freq)), 0, 1
            )[0]
            smooth_Z = gaussian_filter1d(norm_Z[select_idxs], 10)
            ind_max = numpy.argmax(smooth_Z)
            return freq[ind_max]
        else:
            return None
    else:
        return None


def current():
    """Get current trace"""
    return _get_cpp_feature("current")


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
    time = time[:len(voltage)]
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


def depol_block():
    """Check for a depolarization block"""

    # if there is no depolarization block return 1
    # if there is a depolarization block return None
    # subthreshold traces will also return 1

    # Required trace data
    stim_start = _get_cpp_data("stim_start")
    stim_end = _get_cpp_data("stim_end")

    # Required cpp features
    voltage = _get_cpp_feature("voltage")
    time = _get_cpp_feature("time")
    AP_begin_voltage = _get_cpp_feature("AP_begin_voltage")
    stim_start_idx = numpy.flatnonzero(time >= stim_start)[0]
    stim_end_idx = numpy.flatnonzero(time >= stim_end)[0]

    if AP_begin_voltage is None:
        return numpy.array([1])  # if subthreshold no depolarization block
    elif AP_begin_voltage.size:
        depol_block_threshold = numpy.mean(AP_begin_voltage)  # mV
    else:
        depol_block_threshold = -50

    block_min_duration = 50.0  # ms
    long_hyperpol_threshold = -75.0  # mV

    bool_voltage = numpy.array(voltage > depol_block_threshold, dtype=int)
    up_indexes = numpy.flatnonzero(numpy.diff(bool_voltage) == 1)
    down_indexes = numpy.flatnonzero(numpy.diff(bool_voltage) == -1)
    if len(up_indexes) > len(down_indexes):
        down_indexes = numpy.append(down_indexes, [stim_end_idx])

    if len(up_indexes) == 0:
        # if it never gets high enough, that's not a good sign (meaning no
        # spikes)
        return None
    else:
        # if it stays in the depolarization block more than min_duration, flag
        # as depolarization block
        max_depol_duration = numpy.max(
            [time[down_indexes[k]] - time[up_idx] for k,
             up_idx in enumerate(up_indexes)])
        if max_depol_duration > block_min_duration:
            return None

    bool_voltage = numpy.array(voltage > long_hyperpol_threshold, dtype=int)
    up_indexes = numpy.flatnonzero(numpy.diff(bool_voltage) == 1)
    down_indexes = numpy.flatnonzero(numpy.diff(bool_voltage) == -1)
    down_indexes = down_indexes[(down_indexes > stim_start_idx) & (
        down_indexes < stim_end_idx)]
    if len(down_indexes) != 0:
        up_indexes = up_indexes[(up_indexes > stim_start_idx) & (
            up_indexes < stim_end_idx) & (up_indexes > down_indexes[0])]
        if len(up_indexes) < len(down_indexes):
            up_indexes = numpy.append(up_indexes, [stim_end_idx])
        max_hyperpol_duration = numpy.max(
            [time[up_indexes[k]] - time[down_idx] for k,
             down_idx in enumerate(down_indexes)])

        # if it stays in hyperpolarized stage for more than min_duration,
        # flag as depolarization block
        if max_hyperpol_duration > block_min_duration:
            return None

    return numpy.array([1])


def depol_block_bool():
    """Wrapper around the depol_block feature. Returns [1] if depol_block
    is None, [0] otherwise."""

    if depol_block() is None:
        return numpy.array([1])
    else:
        return numpy.array([0])


def spikes_per_burst():
    """Calculate the number of spikes per burst"""

    burst_begin_indices = _get_cpp_feature("burst_begin_indices")
    burst_end_indices = _get_cpp_feature("burst_end_indices")

    if burst_begin_indices is None:
        return None

    ap_per_bursts = []
    for idx_begin, idx_end in zip(burst_begin_indices, burst_end_indices):
        ap_per_bursts.append(idx_end - idx_begin + 1)

    return numpy.array(ap_per_bursts)


def spikes_per_burst_diff():
    """Calculate the diff between the spikes in each burst and the next one"""
    spikes_per_burst_values = spikes_per_burst()
    if spikes_per_burst_values is None:
        return None

    return spikes_per_burst_values[:-1] - spikes_per_burst_values[1:]


def spikes_in_burst1_burst2_diff():
    """Calculate the diff between the spikes in 1st and 2nd bursts"""
    spikes_per_burst_diff_values = spikes_per_burst_diff()
    if spikes_per_burst_diff_values is None or len(
        spikes_per_burst_diff_values
    ) < 1:
        return None

    return numpy.array([spikes_per_burst_diff_values[0]])


def spikes_in_burst1_burstlast_diff():
    """Calculate the diff between the spikes in 1st and last bursts"""
    spikes_per_burst_values = spikes_per_burst()
    if spikes_per_burst_values is None or len(spikes_per_burst_values) < 2:
        return None

    return numpy.array([
        spikes_per_burst_values[0] - spikes_per_burst_values[-1]
    ])


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
