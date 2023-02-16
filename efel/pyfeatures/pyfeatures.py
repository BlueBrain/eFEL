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
from scipy import signal
import efel.cppcore


all_pyfeatures = [
    'voltage',
    'time',
    'current',
    'ISIs',
    'initburst_sahp',
    'initburst_sahp_vb',
    'initburst_sahp_ssse',
    'depol_block',
    'subthreshold_oscillations',
    'subthreshold_frequencies',
]


def voltage():
    """Get voltage trace"""
    return _get_cpp_feature("voltage")


def time():
    """Get time trace"""
    return _get_cpp_feature("time")


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

        # if it stays in hyperpolarzed stage for more than min_duration,
        # flag as depolarization block
        if max_hyperpol_duration > block_min_duration:
            return None

    return numpy.array([1])


import matplotlib.pyplot as plt
def _subthreshold_frequencies():

    threshold = 500.
    freq_max = 0.1
    plot = False

    peak_times = _get_cpp_feature("peak_time")
    if list(peak_times):
        print('traces has spikes, so we cannot compute subthreshold features')
        return []

    stim_start = _get_cpp_data("stim_start")
    stim_end = _get_cpp_data("stim_end")
    stim_start = stim_start + 0.3*(stim_end-stim_start)
    raw_voltage = _get_cpp_feature("voltage")
    raw_time = _get_cpp_feature("time")
    raw_time = raw_time[:len(raw_voltage)]
    time = numpy.linspace(raw_time[0], raw_time[-1], 100001)
    voltage = numpy.interp(time, raw_time, raw_voltage)

    voltage = voltage[(time>stim_start)& (time<stim_end)]
    time = time[(time>stim_start)& (time<stim_end)]

    amps = numpy.fft.rfft(voltage)[1:]
    dt = time[1] - time[0]
    freqs = numpy.fft.rfftfreq(len(time) - 1, dt)[1:]
    power = abs(amps)**2 / dt
    from scipy.signal import savgol_filter
    _power = savgol_filter(power, 50, 3)
    base_power = numpy.median(_power[freqs < freq_max])

    peaks = signal.find_peaks(_power, height=threshold * base_power, width=30)#, )
    peaks = [freqs[peak] for peak in peaks[0] if freqs[peak] < freq_max]
    if len(peaks)>0:
        plt.figure()
        plt.plot(raw_time, raw_voltage)
        plt.savefig('test.pdf')
        plt.close()

        plt.figure()
        plt.axvline(freq_max)
        plt.axhline(base_power)
        plt.axhline(threshold*base_power, c='r')

        for peak in peaks:
            plt.axvline(peak, lw=0.5, c='k')
        plt.plot(freqs, power, lw=0.5)
        plt.loglog(freqs, _power, lw=0.5, c='b')
        #plt.gca().set_xlim(0, 0.2)
        plt.savefig('subthreshold_fft.pdf')
        plt.close()

        plt.figure()
        plt.plot(time, voltage)
        plt.savefig('voltage.pdf')

        plt.close()
    return peaks


def subthreshold_oscillations():
    peaks = _subthreshold_frequencies()
    return numpy.array([1 if len(peaks) > 0 else 0])


def subthreshold_frequencies():
    return _subthreshold_frequencies()


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
