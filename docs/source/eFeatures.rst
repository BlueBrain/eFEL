.. role:: red

=====================
eFeature descriptions
=====================

A pdf document describing the eFeatures is available
`here <http://bluebrain.github.io/eFEL/efeature-documentation.pdf>`_.

Not every eFeature has a description in this document yet,
the complete set will be available shortly.

Implemented eFeatures (to be continued)
=======================================

Spike event features
--------------------

.. image:: _static/figures/inv_ISI.png

LibV1 : time_to_first_spike
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Time from the start of the stimulus to the maximum of the first peak

- **Required features**: peak_time
- **Units**: ms
- **Pseudocode**: ::

    time_to_first_spike = peaktime[0] - stimstart


LibV5 : time_to_second_spike
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Time from the start of the stimulus to the maximum of the second peak

- **Required features**: peak_time
- **Units**: ms
- **Pseudocode**: ::

    time_to_second_spike = peaktime[1] - stimstart


LibV5 : inv_time_to_first_spike
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1.0 over time to first spike; returns 0 when no spike

- **Required features**: time_to_first_spike
- **Units**: Hz
- **Pseudocode**: ::

    if len(time_to_first_spike) > 0:
        inv_time_to_first_spike = 1.0 / time_to_first_spike[0]
    else:
        inv_time_to_first_spike = 0


LibV1 : ISI_values
~~~~~~~~~~~~~~~~~~

The interspike intervals (i.e. time intervals) between adjacent peaks, starting at the second peak.
The 1st spike is not taken into account, because some cells spike right after the stimulus onset and then stay silent for a while.

- **Required features**: peak_time (ms)
- **Units**: ms
- **Pseudocode**: ::

    isi_values = numpy.diff(peak_time)[1:]


LibV1 : doublet_ISI
~~~~~~~~~~~~~~~~~~~

The time interval between the first too peaks

- **Required features**: peak_time (ms)
- **Units**: ms
- **Pseudocode**: ::

    doublet_ISI = peak_time[1] - peak_time[0]


LibV5 : all_ISI_values, inv_first_ISI, inv_second_ISI, inv_third_ISI, inv_fourth_ISI, inv_fifth_ISI, inv_last_ISI
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1.0 over first/second/third/fourth/fith/last ISI; returns 0 when no ISI

- **Required features**: peak_time (ms)
- **Units**: Hz
- **Pseudocode**: ::

    all_isi_values_vec = numpy.diff(peak_time)

    if len(all_isi_values_vec) > 0:
        inv_first_ISI = 1000.0 / all_isi_values_vec[0]
    else:
        inv_first_ISI = 0

    if len(all_isi_values_vec) > 1:
        inv_second_ISI = 1000.0 / all_isi_values_vec[1]
    else:
        inv_second_ISI = 0

    if len(all_isi_values_vec) > 2:
        inv_third_ISI = 1000.0 / all_isi_values_vec[2]
    else:
        inv_third_ISI = 0

    if len(all_isi_values_vec) > 3:
        inv_fourth_ISI = 1000.0 / all_isi_values_vec[3]
    else:
        inv_fourth_ISI = 0

    if len(all_isi_values_vec) > 4:
        inv_fifth_ISI = 1000.0 / all_isi_values_vec[4]
    else:
        inv_fifth_ISI = 0

    if len(all_isi_values_vec) > 0:
        inv_last_ISI = 1000.0 / all_isi_values_vec[-1]
    else:
        inv_last_ISI = 0


LibV5 : time_to_last_spike
~~~~~~~~~~~~~~~~~~~~~~~~~~

time from stimulus start to last spike

- **Required features**: peak_time (ms), stimstart (ms)
- **Units**: ms
- **Pseudocode**: ::

    if len(peak_time) > 0:
        time_to_last_spike = peak_time[-1] - stimstart
    else:
        time_to_last_spike = 0

LibV1 : Spikecount
~~~~~~~~~~~~~~~~~~

number of spikes in the trace, including outside of stimulus interval

- **Required features**: LibV1:peak_indices
- **Units**: constant
- **Pseudocode**: ::

    Spikecount = len(peak_indices)

LibV5 : Spikecount_stimint
~~~~~~~~~~~~~~~~~~~~~~~~~~

number of spikes inside the stimulus interval

- **Required features**: LibV1:peak_time
- **Units**: constant
- **Pseudocode**: ::

    peaktimes_stimint = numpy.where((peak_time >= stim_start) & (peak_time <= stim_end)) 
    Spikecount_stimint = len(peaktimes_stimint)

LibV5 : number_initial_spikes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

number of spikes at the beginning of the stimulus

- **Required features**: LibV1:peak_time
- **Required parameters**: initial_perc (default=0.1)
- **Units**: constant
- **Pseudocode**: ::

    initial_length = (stimend - stimstart) * initial_perc
    number_initial_spikes = len(numpy.where( \
        (peak_time >= stimstart) & \
        (peak_time <= stimstart + initial_length)))

LibV1 : mean_frequency
~~~~~~~~~~~~~~~~~~~~~~

The mean frequency of the firing rate

- **Required features**: stim_start, stim_end, LibV1:peak_time
- **Units**: Hz
- **Pseudocode**: ::

    condition = np.all((stim_start < peak_time, peak_time < stim_end), axis=0)
    spikecount = len(peak_time[condition])
    last_spike_time = peak_time[peak_time < stim_end][-1]
    mean_frequency = 1000 * spikecount / (last_spike_time - stim_start)

LibV5 : ISI_semilog_slope
~~~~~~~~~~~~~~~~~~~~~~~~~

The slope of a linear fit to a semilog plot of the ISI values.

Attention: the 1st ISI is not taken into account. See LibV1: ISI_values feature for more details.

- **Required features**: t, V, stim_start, stim_end, ISI_values
- **Units**: ms
- **Pseudocode**: ::

    x = range(1, len(ISI_values)+1)
    log_ISI_values = numpy.log(ISI_values)
    slope, _ = numpy.polyfit(x, log_ISI_values, 1)

    ISI_semilog_slope = slope

LibV5 : ISI_log_slope
~~~~~~~~~~~~~~~~~~~~~

The slope of a linear fit to a loglog plot of the ISI values.

Attention: the 1st ISI is not taken into account. See LibV1: ISI_values feature for more details.

- **Required features**: t, V, stim_start, stim_end, ISI_values
- **Units**: ms
- **Pseudocode**: ::

    log_x = numpy.log(range(1, len(ISI_values)+1))
    log_ISI_values = numpy.log(ISI_values)
    slope, _ = numpy.polyfit(log_x, log_ISI_values, 1)

    ISI_log_slope = slope

LibV5 : ISI_log_slope_skip
~~~~~~~~~~~~~~~~~~~~~~~~~~

The slope of a linear fit to a loglog plot of the ISI values, but not taking into account the first ISI values.

The proportion of ISI values to be skipped is given by spike_skipf (between 0 and 1). 
However, if this number of ISI values to skip is higher than max_spike_skip, then max_spike_skip is taken instead.

- **Required features**: t, V, stim_start, stim_end, ISI_values
- **Parameters**: spike_skipf (default=0.1), max_spike_skip (default=2)
- **Units**: ms
- **Pseudocode**: ::

    start_idx = min([max_spike_skip, round((len(ISI_values) + 1) * spike_skipf)])
    ISI_values = ISI_values[start_idx:]
    log_x = numpy.log(range(1, len(ISI_values)+1))
    log_ISI_values = numpy.log(ISI_values)
    slope, _ = numpy.polyfit(log_x, log_ISI_values, 1)

    ISI_log_slope = slope

LibV1 : ISI_CV
~~~~~~~~~~~~~~

The coefficient of variation of the ISIs.

Attention: the 1st ISI is not taken into account. See LibV1: ISI_values feature for more details.

- **Required features**: ISI_values
- **Units**: constant
- **Pseudocode**: ::

    ISI_mean = numpy.mean(ISI_values)
    ISI_variance = numpy.sum(numpy.square(ISI_values-ISI_mean)) / (len(ISI_values)-1)
    ISI_std = math.sqrt(ISI_variance)
    ISI_CV = ISI_std / ISI_mean

LibV5 : irregularity_index
~~~~~~~~~~~~~~~~~~~~~~~~~~

Mean of the absolute difference of all ISIs, except the first one (see LibV1: ISI_values feature for more details.)

- **Required features**: ISI_values
- **Units**: ms
- **Pseudocode**: ::

    irregularity_index = numpy.mean(numpy.absolute(ISI_values[1:] - ISI_values[:-1]))


LibV5 : adaptation_index
~~~~~~~~~~~~~~~~~~~~~~~~

Normalized average difference of two consecutive ISIs, skipping the first ISIs

The proportion of ISI values to be skipped is given by spike_skipf (between 0 and 1). 
However, if this number of ISI values to skip is higher than max_spike_skip, then max_spike_skip is taken instead.

The adaptation index is zero for a constant firing rate and bigger than zero for a decreasing firing rate

- **Required features**: stim_start, stim_end, peak_time
- **Parameters**: offset (default=0), spike_skipf (default=0.1), max_spike_skip (default=2)
- **Units**: constant
- **Pseudocode**: ::

    # skip the first ISIs
    peak_selection = [peak_time >= stim_start - offset, peak_time <= stim_end - offset]
    spike_time = peak_time[numpy.all(peak_selection, axis=0)]

    start_idx = min([max_spike_skip, round(len(spike_time) * spike_skipf)])
    spike_time = spike_time[start_idx:]

    # compute the adaptation index
    ISI_values = spike_time[1:] - spike_time[:-1]
    ISI_sum = ISI_values[1:] + ISI_values[:-1]
    ISI_sub = ISI_values[1:] - ISI_values[:-1]
    adaptation_index = numpy.mean(ISI_sum / ISI_sub)


LibV5 : adaptation_index_2
~~~~~~~~~~~~~~~~~~~~~~~~~~

Normalized average difference of two consecutive ISIs, starting at the second ISI

The adaptation index is zero for a constant firing rate and bigger than zero for a decreasing firing rate

- **Required features**: stim_start, stim_end, peak_time
- **Parameters**: offset (default=0)
- **Units**: constant
- **Pseudocode**: ::

    # skip the first ISI
    peak_selection = [peak_time >= stim_start - offset, peak_time <= stim_end - offset]
    spike_time = peak_time[numpy.all(peak_selection, axis=0)]

    spike_time = spike_time[1:]

    # compute the adaptation index
    ISI_values = spike_time[1:] - spike_time[:-1]
    ISI_sum = ISI_values[1:] + ISI_values[:-1]
    ISI_sub = ISI_values[1:] - ISI_values[:-1]
    adaptation_index = numpy.mean(ISI_sum / ISI_sub)


LibV5 : check_AISInitiation
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Check initiation of AP in AIS

- **Required features**: t, V, stim_start, stim_end, AP_begin_time, AP_begin_time;location_AIS
- **Units**: constant
- **Pseudocode**: ::

    if len(AP_begin_time) != len(AP_begin_time;location_AIS):
        return None
    for soma_time, ais_time in zip(AP_begin_time, AP_begin_time;location_AIS):
        if soma_time < ais_time:
            return None
    return [1]

LibV1 : burst_mean_freq
~~~~~~~~~~~~~~~~~~~~~~~

The mean frequency during a burst for each burst

If burst_ISI_indices did not detect any burst beginning,
then the spikes are not considered to be part of any burst

- **Required features**: burst_ISI_indices, peak_time
- **Units**: Hz
- **Pseudocode**: ::

    if burst_ISI_indices is None:
        return None
    elif len(burst_ISI_indices) == 0:
        return []

    burst_mean_freq = []
    burst_index = numpy.insert(
        burst_index_tmp, burst_index_tmp.size, len(peak_time) - 1
    )

    # 1st burst
    span = peak_time[burst_index[0]] - peak_time[0]
    N_peaks = burst_index[0] + 1
    burst_mean_freq.append(N_peaks * 1000 / span)

    for i, burst_idx in enumerate(burst_index[:-1]):
        if burst_index[i + 1] - burst_idx != 1:
            span = peak_time[burst_index[i + 1]] - peak_time[burst_idx + 1]
            N_peaks = burst_index[i + 1] - burst_idx
            burst_mean_freq.append(N_peaks * 1000 / span)

    return burst_mean_freq

LibV1 : burst_number
~~~~~~~~~~~~~~~~~~~~

The number of bursts

- **Required features**: burst_mean_freq
- **Units**: constant
- **Pseudocode**: ::

    burst_number = len(burst_mean_freq)

LibV1 : interburst_voltage
~~~~~~~~~~~~~~~~~~~~~~~~~~

The voltage average in between two bursts

Iterating over the burst ISI indices determine the last peak before the burst. 
Starting 5 ms after that peak take the voltage average until 5 ms before the first peak of the subsequent burst.

- **Required features**: burst_ISI_indices, peak_indices
- **Units**: mV
- **Pseudocode**: ::

    interburst_voltage = []
    for idx in burst_ISI_idxs:
        ts_idx = peak_idxs[idx]
        t_start = time[ts_idx] + 5
        start_idx = np.argwhere(time < t_start)[-1][0]

        te_idx = peak_idxs[idx + 1]
        t_end = time[te_idx] - 5
        end_idx = np.argwhere(time > t_end)[0][0]

        interburst_voltage.append(np.mean(voltage[start_idx:end_idx + 1]))

LibV1 : single_burst_ratio
~~~~~~~~~~~~~~~~~~~~

Length of the second isi over the median of the rest of the isis. The first isi is not taken into account, because it could bias the feature.
See LibV1: ISI_values feature for more details.

- **Required features**: ISI_values
- **Units**: constant
- **Pseudocode**: ::

    single_burst_ratio = ISI_values[0] / numpy.mean(ISI_values)


Spike shape features
--------------------

.. image:: _static/figures/AP_Amplitude.png

LibV1 : peak_time
~~~~~~~~~~~~~~~~~

The times of the maxima of the peaks

- **Required features**: LibV5:peak_indices
- **Units**: ms
- **Pseudocode**: ::

    peak_time = time[peak_indices]


LibV1 : peak_voltage
~~~~~~~~~~~~~~~~~~~~

The voltages at the maxima of the peaks

- **Required features**: LibV5:peak_indices
- **Units**: mV
- **Pseudocode**: ::

    peak_voltage = voltage[peak_indices]


LibV1 : AP_Amplitude, AP1_amp, AP2_amp, APlast_amp
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The relative height of the action potential from spike onset

- **Required features**: LibV5:AP_begin_indices, LibV1:peak_voltage (mV)
- **Units**: mV
- **Pseudocode**: ::

    AP_Amplitude = peak_voltage - voltage[AP_begin_indices]
    AP1_amp = AP_Amplitude[0]
    AP2_amp = AP_Amplitude[1]
    APlast_amp = AP_Amplitude[-1]

LibV5 : AP_Amplitude_from_voltagebase
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The relative height of the action potential from voltage base

- **Required features**: LibV5:voltage_base, LibV1:peak_voltage (mV)
- **Units**: mV
- **Pseudocode**: ::

    AP_Amplitude_from_voltagebase = peak_voltage - voltage_base

LibV5 : AP1_peak, AP2_peak
~~~~~~~~~~~~~~~~~~~~~~~~~~

The peak voltage of the first and second action potentials

- **Required features**: LibV1:peak_voltage (mV)
- **Units**: mV
- **Pseudocode**: ::

    AP1_peak = peak_voltage[0]
    AP2_peak = peak_voltage[1]

LibV5 : AP2_AP1_diff
~~~~~~~~~~~~~~~~~~~~

Difference amplitude of the second to first spike

- **Required features**: LibV1:AP_amplitude (mV)
- **Units**: mV
- **Pseudocode**: ::

    AP2_AP1_diff = AP_amplitude[1] - AP_amplitude[0]

LibV5 : AP2_AP1_peak_diff
~~~~~~~~~~~~~~~~~~~~~~~~~

Difference peak voltage of the second to first spike

- **Required features**: LibV1:peak_voltage (mV)
- **Units**: mV
- **Pseudocode**: ::

    AP2_AP1_diff = peak_voltage[1] - peak_voltage[0]

.. image:: _static/figures/AHP.png

LibV5 : AHP_depth_abs
~~~~~~~~~~~~~~~~~~~~~

Absolute voltage values at the first after-hyperpolarization

- **Required features**: LibV5:min_AHP_values (mV)
- **Units**: mV

LibV1 : AHP_depth_abs_slow
~~~~~~~~~~~~~~~~~~~~~~~~~~

Absolute voltage values at the first after-hyperpolarization starting 
a given number of ms (default: 5) after the peak

- **Required features**: LibV1:peak_indices
- **Units**: mV

LibV1 : AHP_slow_time
~~~~~~~~~~~~~~~~~~~~~

Time difference between slow AHP (see AHP_depth_abs_slow) and peak, divided by
interspike interval 

- **Required features**: LibV1:AHP_depth_abs_slow
- **Units**: constant
  
LibV1 : AHP_depth
~~~~~~~~~~~~~~~~~

Relative voltage values at the first after-hyperpolarization

- **Required features**: LibV1:voltage_base (mV), LibV5:min_AHP_values (mV)
- **Units**: mV
- **Pseudocode**: ::

    min_AHP_values = first_min_element(voltage, peak_indices)
    AHP_depth = min_AHP_values[:] - voltage_base

LibV5 : AHP_depth_from_peak, AHP1_depth_from_peak, AHP2_depth_from_peak
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Voltage difference between AP peaks and first AHP depths

- **Required features**: LibV1:peak_indices, LibV5:min_AHP_indices
- **Units**: mV
- **Pseudocode**: ::

    AHP_depth_from_peak =  v[peak_indices] - v[min_AHP_indices]
    AHP1_depth_from_peak = AHP_depth_from_peak[0]
    AHP2_depth_from_peak = AHP_depth_from_peak[1]

LibV5 : AHP_time_from_peak
~~~~~~~~~~~~~~~~~~~~~~~~~~

Time between AP peaks and first AHP depths

- **Required features**: LibV1:peak_indices, LibV5:min_AHP_values (mV)
- **Units**: ms
- **Pseudocode**: ::

    min_AHP_indices = first_min_element(voltage, peak_indices)
    AHP_time_from_peak = t[min_AHP_indices[:]] - t[peak_indices[i]]

LibV5 : min_voltage_between_spikes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Minimal voltage between consecutive spikes

- **Required features**: LibV5:peak_indices
- **Units**: mV
- **Pseudocode**: ::

    min_voltage_between_spikes = []
    for peak1, peak2 in zip(peak_indices[:-1], peak_indices[1:]):
        min_voltage_between_spikes.append(numpy.min(voltage[peak1:peak2]))

LibV5 : min_between_peaks_values
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Minimal voltage between consecutive spikes

The last value is the minimum between last spike and stimulus end
if strict stiminterval is True, and minimum between last spike and last value
if strict stiminterval is False


- **Required features**: LibV5:min_between_peaks_indices
- **Units**: mV
- **Pseudocode**: ::

    min_between_peaks_values = v[min_between_peaks_indices]


.. image:: _static/figures/AP_duration_half_width.png


LibV2 : AP_duration_half_width
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Width of spike at half spike amplitude, with spike onset as described in LibV5: AP_begin_time

- **Required features**: LibV2: AP_rise_indices, LibV2: AP_fall_indices
- **Units**: ms
- **Pseudocode**: ::

    AP_rise_indices = index_before_peak((v(peak_indices) - v(AP_begin_indices)) / 2)
    AP_fall_indices = index_after_peak((v(peak_indices) - v(AP_begin_indices)) / 2)
    AP_duration_half_width = t(AP_fall_indices) - t(AP_rise_indices)

LibV1 : AP_width
~~~~~~~~~~~~~~~~

Width of spike at threshold, bounded by minimum AHP

Can use strict_stiminterval to not use minimum after stimulus end.

- **Required features**: LibV1: peak_indices, LibV5: min_AHP_indices, threshold
- **Units**: ms
- **Pseudocode**: ::

    min_AHP_indices = numpy.concatenate([[stim_start], min_AHP_indices])
    for i in range(len(min_AHP_indices)-1):
        onset_index = numpy.where(v[min_AHP_indices[i]:min_AHP_indices[i+1]] > threshold)[0]
        onset_time[i] = t[onset_index]
        offset_time[i] = t[numpy.where(v[onset_index:min_AHP_indices[i+1]] < threshold)[0]]
        AP_width[i] = t(offset_time[i]) - t(onset_time[i])

LibV5 : AP_width_between_threshold
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Width of spike at threshold, bounded by minimum between peaks

Can use strict_stiminterval to not use minimum after stimulus end.

- **Required features**: LibV1: peak_indices, LibV5: min_between_peaks_indices, threshold
- **Units**: ms
- **Pseudocode**: ::

    min_between_peaks_indices = numpy.concatenate([[stim_start], min_between_peaks_indices])
    for i in range(len(min_between_peaks_indices)-1):
        onset_index = numpy.where(v[min_between_peaks_indices[i]:min_between_peaks_indices[i+1]] > threshold)[0]
        onset_time[i] = t[onset_index]
        offset_time[i] = t[numpy.where(v[onset_index:min_between_peaks_indices[i+1]] < threshold)[0]]
        AP_width[i] = t(offset_time[i]) - t(onset_time[i])

LibV5 : spike_half_width, AP1_width, AP2_width, APlast_width
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Width of spike at half spike amplitude, 
with the spike amplitude taken as the difference between the minimum between two peaks and the next peak

- **Required features**: LibV5: peak_indices, LibV5: min_AHP_indices
- **Units**: ms
- **Pseudocode**: ::

    min_AHP_indices = numpy.concatenate([[stim_start], min_AHP_indices])
    for i in range(1, len(min_AHP_indices)):
        v_half_width = (v[peak_indices[i-1]] + v[min_AHP_indices[i]]) / 2.
        rise_idx = numpy.where(v[min_AHP_indices[i-1]:peak_indices[i-1]] > v_half_width)[0]
        v_dev = v_half_width - v[rise_idx]
        delta_v = v[rise_idx] - v[rise_idx - 1]
        delta_t = t[rise_idx] - t[rise_idx - 1]
        t_dev_rise = delta_t * v_dev / delta_v
        
        fall_idx = numpy.where(v[peak_indices[i-1]:min_AHP_indices[i]] < v_half_width)[0]
        v_dev = v_half_width - v[fall_idx]
        delta_v = v[fall_idx] - v[fall_idx - 1]
        delta_t = t[fall_idx] - t[fall_idx - 1]
        t_dev_fall = delta_t * v_dev / delta_v
        spike_half_width[i] = t[fall_idx] + t_dev_fall - t[rise_idx] - t_dev_rise

    AP1_width = spike_half_width[0]
    AP2_width = spike_half_width[1]
    APlast_width = spike_half_width[-1]


LibV1 : spike_width2
~~~~~~~~~~~~~~~~~~~~

Width of spike at half spike amplitude, with the spike onset taken as the maximum of the second derivative of the voltage in the range between
the minimum between two peaks and the next peak

- **Required features**: LibV5: peak_indices, LibV5: min_AHP_indices
- **Units**: ms
- **Pseudocode**: ::

    for i in range(len(min_AHP_indices)):
        dv2 = CentralDiffDerivative(CentralDiffDerivative(v[min_AHP_indices[i]:peak_indices[i + 1]]))
        peak_onset_idx = numpy.argmax(dv2) + min_AHP_indices[i]
        v_half_width = (v[peak_indices[i + 1]] + v[peak_onset_idx]) / 2.

        rise_idx = numpy.where(v[peak_onset_idx:peak_indices[i + 1]] > v_half_width)[0]
        v_dev = v_half_width - v[rise_idx]
        delta_v = v[rise_idx] - v[rise_idx - 1]
        delta_t = t[rise_idx] - t[rise_idx - 1]
        t_dev_rise = delta_t * v_dev / delta_v
        
        fall_idx = numpy.where(v[peak_indices[i + 1]:] < v_half_width)[0]
        v_dev = v_half_width - v[fall_idx]
        delta_v = v[fall_idx] - v[fall_idx - 1]
        delta_t = t[fall_idx] - t[fall_idx - 1]
        t_dev_fall = delta_t * v_dev / delta_v
        spike_width2[i] = t[fall_idx] + t_dev_fall - t[rise_idx] - t_dev_rise


LibV5 : AP_begin_width, AP1_begin_width, AP2_begin_width
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Width of spike at spike start

- **Required features**: LibV5: min_AHP_indices, LibV5: AP_begin_indices
- **Units**: ms
- **Pseudocode**: ::

    for i in range(len(min_AHP_indices)):
        rise_idx = AP_begin_indices[i]
        fall_idx = numpy.where(v[rise_idx + 1:min_AHP_indices[i]] < v[rise_idx])[0]
        AP_begin_width[i] = t[fall_idx] - t[rise_idx]

    AP1_begin_width = AP_begin_width[0]
    AP2_begin_width = AP_begin_width[1]

LibV5 : AP2_AP1_begin_width_diff
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Difference width of the second to first spike

- **Required features**: LibV5: AP_begin_width
- **Units**: ms
- **Pseudocode**: ::

    AP2_AP1_begin_width_diff = AP_begin_width[1] - AP_begin_width[0]

LibV5 : AP_begin_voltage, AP1_begin_voltage, AP2_begin_voltage
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Voltage at spike start

- **Required features**:  LibV5: AP_begin_indices
- **Units**: mV
- **Pseudocode**: ::

    AP_begin_voltage = v[AP_begin_indices]
    AP1_begin_voltage = AP_begin_voltage[0]
    AP2_begin_voltage = AP_begin_voltage[1]

LibV5 : AP_begin_time
~~~~~~~~~~~~~~~~~~~~~

Time at spike start. Spike start is defined as where the first derivative of the voltage trace is higher than 10 V/s , for at least 5 points

- **Required features**:  LibV5: AP_begin_indices
- **Units**: ms
- **Pseudocode**: ::

    AP_begin_time = t[AP_begin_indices]

LibV5 : AP_peak_upstroke
~~~~~~~~~~~~~~~~~~~~~~~~

Maximum of rise rate of spike

- **Required features**: LibV5: AP_begin_indices, LibV5: peak_indices
- **Units**: V/s
- **Pseudocode**: ::

    ap_peak_upstroke = []
    for apbi, pi in zip(ap_begin_indices, peak_indices):
        ap_peak_upstroke.append(numpy.max(dvdt[apbi:pi]))


LibV5 : AP_peak_downstroke
~~~~~~~~~~~~~~~~~~~~~~~~~~

Minimum of fall rate from spike

- **Required features**: LibV5: min_AHP_indices, LibV5: peak_indices
- **Units**: V/s
- **Pseudocode**: ::

    ap_peak_downstroke = []
    for ahpi, pi in zip(min_ahp_indices, peak_indices):
        ap_peak_downstroke.append(numpy.min(dvdt[pi:ahpi]))

LibV2 : AP_rise_time
~~~~~~~~~~~~~~~~~~~~

Time between the AP threshold and the peak, given a window
(default: from 0% to 100% of the AP amplitude)

- **Required features**: LibV5: AP_begin_indices, LibV5: peak_indices, LibV1: AP_amplitude
- **Units**: ms
- **Pseudocode**: ::

    rise_times = []
    begin_voltages = AP_amps * rise_start_perc + voltage[AP_begin_indices]
    end_voltages = AP_amps * rise_end_perc + voltage[AP_begin_indices]

    for AP_begin_indice, peak_indice, begin_v, end_v in zip(
        AP_begin_indices, peak_indices, begin_voltages, end_voltages
    ):
        voltage_window = voltage[AP_begin_indice:peak_indice]

        new_begin_indice = AP_begin_indice + numpy.min(
            numpy.where(voltage_window >= begin_v)[0]
        )
        new_end_indice = AP_begin_indice + numpy.max(
            numpy.where(voltage_window <= end_v)[0]
        )

        rise_times.append(time[new_end_indice] - time[new_begin_indice])

LibV5 : AP_phaseslope
~~~~~~~~~~~~~~~~~~~~~~

Slope of the V, dVdt phasespace plot at the beginning of every spike

(at the point where the derivative crosses the DerivativeThreshold)

- **Required features**: LibV5:AP_begin_indices
- **Parameters**: AP_phaseslope_range
- **Units**: 1/(ms)
- **Pseudocode**: ::

    range_max_idxs = AP_begin_indices + AP_phseslope_range
    range_min_idxs = AP_begin_indices - AP_phseslope_range
    AP_phaseslope = (dvdt[range_max_idxs] - dvdt[range_min_idxs]) / (v[range_max_idxs] - v[range_min_idxs])

LibV5 : AP_phaseslope_AIS
~~~~~~~~~~~~~~~~~~~~~~~~~

Same as AP_phaseslope, but for AIS location

Please, notice that you have to provide t, v, stim_start and stim_end for location.

- **Required features**: T;location_AIS, V;location_AIS, stim_start;location_AIS, stim_end;location_AIS, LibV5:AP_begin_indices;location_AIS
- **Parameters**: AP_phaseslope_range
- **Units**: 1/(ms)
- **Pseudocode**: ::

    range_max_idxs = AP_begin_indices + AP_phseslope_range
    range_min_idxs = AP_begin_indices - AP_phseslope_range
    AP_phaseslope_AIS = (dvdt[range_max_idxs] - dvdt[range_min_idxs]) / (v[range_max_idxs] - v[range_min_idxs])

LibV5 : BPAPHeightLoc1
~~~~~~~~~~~~~~~~~~~~~~

Voltage height (difference betwen peaks and voltage base) at dendrite location

Please, notice that you have to provide t, v, stim_start and stim_end for location.

- **Required features**: T;location_dend1, V;location_dend1, stim_start;location_dend1, stim_end;location_dend1, peak_voltage;location_dend1, voltage_base;location_dend1
- **Units**: mV
- **Pseudocode**: ::

    BPAPHeightLoc1 = peak_voltage - voltage_base

LibV5 : BPAPHeightLoc2
~~~~~~~~~~~~~~~~~~~~~~

Same as BPAPHeightLoc1, but for dend2 location

- **Required features**: T;location_dend2, V;location_dend2, stim_start;location_dend2, stim_end;location_dend2, peak_voltage;location_dend2, voltage_base;location_dend2
- **Units**: mV
- **Pseudocode**: ::

    BPAPHeightLoc2 = peak_voltage - voltage_base

LibV5 : BPAPAmplitudeLoc1
~~~~~~~~~~~~~~~~~~~~~~~~~

Amplitude at dendrite location

Please, notice that you have to provide t, v, stim_start and stim_end for location.

- **Required features**: T;location_dend1, V;location_dend1, stim_start;location_dend1, stim_end;location_dend1, peak_voltage;location_dend1, AP_begin_voltage;location_dend1
- **Units**: mV
- **Pseudocode**: ::

    BPAPAmplitudeLoc1 = peak_voltage - AP_begin_voltage

LibV5 : BPAPAmplitudeLoc2
~~~~~~~~~~~~~~~~~~~~~~~~~

Same as BPAPAmplitudeLoc1, but for dend2 location

- **Required features**: T;location_dend2, V;location_dend2, stim_start;location_dend2, stim_end;location_dend2, peak_voltage;location_dend2, AP_begin_voltage;location_dend2
- **Units**: mV
- **Pseudocode**: ::

    BPAPAmplitudeLoc2 = peak_voltage - AP_begin_voltage

LibV5 : BAC_width
~~~~~~~~~~~~~~~~~

AP width at epsp location

Please, notice that you have to provide t, v, stim_start and stim_end for location.

- **Required features**: T;location_epsp, V;location_epsp, stim_start;location_epsp, stim_end;location_epsp, AP_width;location_epsp
- **Units**: ms
- **Pseudocode**: ::

    BAC_width = AP_width

LibV5 : BAC_maximum_voltage
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Maximuum voltage at epsp location

Please, notice that you have to provide t, v, stim_start and stim_end for location.

- **Required features**: T;location_epsp, V;location_epsp, stim_start;location_epsp, stim_end;location_epsp, maximum_voltage;location_epsp
- **Units**: mV
- **Pseudocode**: ::

    BAC_maximum_voltage = maximum_voltage



Voltage features
----------------

.. image:: _static/figures/voltage_features.png


LibV5 : steady_state_voltage_stimend
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The average voltage during the last 10% of the stimulus duration.

- **Required features**: t, V, stim_start, stim_end
- **Units**: mV
- **Pseudocode**: ::

    stim_duration = stim_end - stim_start
    begin_time = stim_end - 0.1 * stim_duration
    end_time = stim_end
    steady_state_voltage_stimend = numpy.mean(voltage[numpy.where((t < end_time) & (t >= begin_time))])


LibV1 : steady_state_voltage
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The average voltage after the stimulus

- **Required features**: t, V, stim_end
- **Units**: mV
- **Pseudocode**: ::

    steady_state_voltage = numpy.mean(voltage[numpy.where((t <= max(t)) & (t > stim_end))])


LibV5 : voltage_base
~~~~~~~~~~~~~~~~~~~~

The average voltage during the last 10% of time before the stimulus.

- **Required features**: t, V, stim_start, stim_end
- **Parameters**: voltage_base_start_perc (default = 0.9), voltage_base_end_perc (default = 1.0)
- **Units**: mV
- **Pseudocode**: ::

    voltage_base = numpy.mean(voltage[numpy.where(
        (t >= voltage_base_start_perc * stim_start) &
        (t <= voltage_base_end_perc * stim_start))])

LibV5 : current_base
~~~~~~~~~~~~~~~~~~~~

The average current during the last 10% of time before the stimulus.

- **Required features**: t, I, stim_start, stim_end
- **Parameters**: current_base_start_perc (default = 0.9), current_base_end_perc (default = 1.0), precision_threshold (default = 1e-10), current_base_mode (can be "mean" or "median", default="mean")
- **Units**: mV
- **Pseudocode**: ::

    current_slice = I[numpy.where(
        (t >= current_base_start_perc * stim_start) &
        (t <= current_base_end_perc * stim_start))]
    if current_base_mode == "mean":
        current_base = numpy.mean(current_slice)
    elif current_base_mode == "median":
        current_base = numpy.median(current_slice)

LibV5 : decay_time_constant_after_stim
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The decay time constant of the voltage right after the stimulus

- **Required features**: t, V, stim_start, stim_end
- **Parameters**: decay_start_after_stim (default = 1.0 ms), decay_end_after_stim (default = 10.0 ms)
- **Units**: ms
- **Pseudocode**: ::

    time_interval = t[numpy.where(t => decay_start_after_stim &
                       t < decay_end_after_stim)] - t[numpy.where(t == stim_end)]
    voltage_interval = abs(voltages[numpy.where(t => decay_start_after_stim &
                                    t < decay_end_after_stim)]
                           - voltages[numpy.where(t == decay_start_after_stim)])

    log_voltage_interval = numpy.log(voltage_interval)
    slope, _ = numpy.polyfit(time_interval, log_voltage_interval, 1)

    decay_time_constant_after_stim = -1. / slope

LibV5 : sag_time_constant
~~~~~~~~~~~~~~~~~~~~~~~~~~

The decay time constant of the exponential voltage decay from the bottom of the sag to the steady-state.

The start of the decay is taken at the minimum voltage (the bottom of the sag).
The end of the decay is taken when the voltage crosses the steady state voltage minus 10% of the sag amplitude.
The time constant is the slope of the linear fit to the log of the voltage.
The golden search algorithm is not used, since the data is expected to be noisy and adding a parameter in the log
( log(voltage + x) ) is likely to increase errors on the fit.

- **Required features**: t, V, stim_start, stim_end, minimum_voltage, steady_state_voltage_stimend, sag_amplitude
- **Units**: ms
- **Pseudocode**: ::

    # get start decay
    start_decay = numpy.argmin(vinterval)

    # get end decay
    v90 = steady_state_v - 0.1 * sag_ampl
    end_decay = numpy.where((tinterval > tinterval[start_decay]) & (vinterval >= v90))[0][0]

    v_reference = vinterval[end_decay]

    # select t, v in decay interval
    interval_indices = numpy.arange(start_decay, end_decay)
    interval_time = tinterval[interval_indices]
    interval_voltage = abs(vinterval[interval_indices] - v_reference)

    # get tau
    log_interval_voltage = numpy.log(interval_voltage)
    slope, _ = numpy.polyfit(interval_time, log_interval_voltage, 1)
    tau = abs(1. / slope)

.. image:: _static/figures/sag.png

LibV5 : sag_amplitude
~~~~~~~~~~~~~~~~~~~~~

The difference between the minimal voltage and the steady state at stimend

- **Required features**: t, V, stim_start, stim_end, steady_state_voltage_stimend, minimum_voltage, voltage_deflection_stim_ssse
- **Parameters**: 
- **Units**: mV
- **Pseudocode**: ::

    if (voltage_deflection_stim_ssse <= 0):
        sag_amplitude = steady_state_voltage_stimend - minimum_voltage
    else:
        sag_amplitude = None


LibV5 : sag_ratio1
~~~~~~~~~~~~~~~~~~

The ratio between the sag amplitude and the maximal sag extend from voltage base

- **Required features**: t, V, stim_start, stim_end, sag_amplitude, voltage_base, minimum_voltage
- **Parameters**: 
- **Units**: constant
- **Pseudocode**: ::

    if voltage_base != minimum_voltage:
        sag_ratio1 = sag_amplitude / (voltage_base - minimum_voltage)
    else:
        sag_ratio1 = None

LibV5 : sag_ratio2
~~~~~~~~~~~~~~~~~~

The ratio between the maximal extends of sag from steady state and voltage base

- **Required features**: t, V, stim_start, stim_end, steady_state_voltage_stimend, voltage_base, minimum_voltage
- **Parameters**: 
- **Units**: constant
- **Pseudocode**: ::

    if voltage_base != minimum_voltage:
        sag_ratio2 = (voltage_base - steady_state_voltage_stimend) / (voltage_base - minimum_voltage)
    else:
        sag_ratio2 = None

LibV1 : ohmic_input_resistance
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ratio between the voltage deflection and stimulus current

- **Required features**: t, V, stim_start, stim_end, voltage_deflection
- **Parameters**: stimulus_current
- **Units**: mV/nA
- **Pseudocode**: ::

    ohmic_input_resistance = voltage_deflection / stimulus_current

LibV5 : ohmic_input_resistance_vb_ssse
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ratio between the voltage deflection (between voltage base and steady-state voltage at stimend) and stimulus current

- **Required features**: t, V, stim_start, stim_end, voltage_deflection_vb_ssse
- **Parameters**: stimulus_current
- **Units**: mV/nA
- **Pseudocode**: ::

    ohmic_input_resistance_vb_ssse = voltage_deflection_vb_ssse / stimulus_current

LibV5 : voltage_deflection_vb_ssse
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The voltage deflection between voltage base and steady-state voltage at stimend

The voltage base used is the average voltage during the last 10% of time before the stimulus
and the steady state voltage at stimend used is
the average voltage during the last 10% of the stimulus duration.

- **Required features**: t, V, stim_start, stim_end, voltage_base, steady_state_voltage_stimend
- **Units**: mV
- **Pseudocode**: ::

    voltage_deflection_vb_ssse = steady_state_voltage_stimend - voltage_base

LibV1 : voltage_deflection
~~~~~~~~~~~~~~~~~~~~~~~~~~
    
The voltage deflection between voltage base and steady-state voltage at stimend

The voltage base used is the average voltage during all of the time before the stimulus
and the steady state voltage at stimend used is
the average voltage of the five values before the last five values
before the end of the stimulus duration.

- **Required features**: t, V, stim_start, stim_end
- **Units**: mV
- **Pseudocode**: ::

    voltage_base = numpy.mean(V[t < stim_start])
    stim_end_idx = numpy.where(t >= stim_end)[0][0]
    steady_state_voltage_stimend = numpy.mean(V[stim_end_idx-10:stim_end_idx-5])
    voltage_deflection = steady_state_voltage_stimend - voltage_base

LibV5 : voltage_deflection_begin
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
The voltage deflection between voltage base and steady-state voltage soon after stimulation start

The voltage base used is the average voltage during all of the time before the stimulus
and the steady state voltage used is
the average voltage taken from 5% to 15% of the stimulus duration.

- **Required features**: t, V, stim_start, stim_end
- **Units**: mV
- **Pseudocode**: ::

    voltage_base = numpy.mean(V[t < stim_start])
    tstart = stim_start + 0.05 * (stim_end - stim_start)
    tend = stim_start + 0.15 * (stim_end - stim_start)
    condition = numpy.all((tstart < t, t < tend), axis=0)
    steady_state_voltage_stimend = numpy.mean(V[condition])
    voltage_deflection = steady_state_voltage_stimend - voltage_base

LibV5 : voltage_after_stim
~~~~~~~~~~~~~~~~~~~~~~~~~~
    
The mean voltage after the stimulus in
(stim_end + 25%*end_period, stim_end + 75%*end_period)

- **Required features**: t, V, stim_end
- **Units**: mV
- **Pseudocode**: ::

    tstart = stim_end + (t[-1] - stimEnd) * 0.25
    tend = stim_end + (t[-1] - stimEnd) * 0.75
    condition = numpy.all((tstart < t, t < tend), axis=0)
    voltage_after_stim = numpy.mean(V[condition])

LibV1: minimum_voltage
~~~~~~~~~~~~~~~~~~~~~~

The minimum of the voltage during the stimulus

- **Required features**: t, V, stim_start, stim_end
- **Units**: mV
- **Pseudocode**: ::

    minimum_voltage = min(voltage[numpy.where((t >= stim_start) & (t <= stim_end))])

LibV1: maximum_voltage
~~~~~~~~~~~~~~~~~~~~~~

The maximum of the voltage during the stimulus

- **Required features**: t, V, stim_start, stim_end
- **Units**: mV
- **Pseudocode**: ::

    maximum_voltage = max(voltage[numpy.where((t >= stim_start) & (t <= stim_end))])

LibV5: maximum_voltage_from_voltagebase
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Difference between maximum voltage during stimulus and voltage base

- **Required features**: maximum_voltage, voltage_base
- **Units**: mV
- **Pseudocode**: ::

    maximum_voltage_from_voltagebase = maximum_voltage - voltage_base



Requested eFeatures
===================

LibV1 : AHP_depth_last
~~~~~~~~~~~~~~~~~~~~~~

Relative voltage values at the last after-hyperpolarization

- **Required features**: LibV1:voltage_base (mV), LibV5:last_AHP_values (mV)
- **Units**: mV
- **Pseudocode**: ::

    last_AHP_values = last_min_element(voltage, peak_indices)
    AHP_depth = last_AHP_values[:] - voltage_base


LibV5 : AHP_time_from_peak_last
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Time between AP peaks and last AHP depths

- **Required features**: LibV1:peak_indices, LibV5:min_AHP_values (mV)
- **Units**: mV
- **Pseudocode**: ::

    last_AHP_indices = last_min_element(voltage, peak_indices)
    AHP_time_from_peak_last = t[last_AHP_indices[:]] - t[peak_indices[i]]


LibV5 : steady_state_voltage_stimend_from_voltage_base
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The average voltage during the last 90% of the stimulus duration realtive to voltage_base

- **Required features**: LibV5: steady_state_voltage_stimend (mV), LibV5: voltage_base (mV)
- **Units**: mV
- **Pseudocode**: ::

    steady_state_voltage_stimend_from_voltage_base = steady_state_voltage_stimend - voltage_base


LibV5 : min_duringstim_from_voltage_base
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The minimum voltage during stimulus

- **Required features**: LibV5: min_duringstim (mV), LibV5: voltage_base (mV)
- **Units**: mV
- **Pseudocode**: ::

    min_duringstim_from_voltage_base = minimum_voltage - voltage_base


LibV5 : max_duringstim_from_voltage_base
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The minimum voltage during stimulus

- **Required features**: LibV5: max_duringstim (mV), LibV5: voltage_base (mV)
- **Units**: mV
- **Pseudocode**: ::

    max_duringstim_from_voltage_base = maximum_voltage - voltage_base

LibV5 : diff_max_duringstim
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Difference between maximum and steady state during stimulation

- **Required features**: LibV5: max_duringstim (mV), LibV5: steady_state_voltage_stimend (mV)
- **Units**: mV
- **Pseudocode**: ::

    diff_max_duringstim: max_duringstim - steady_state_voltage_stimend

LibV5 : diff_min_duringstim
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Difference between minimum and steady state during stimulation

- **Required features**: LibV5: min_duringstim (mV), LibV5: steady_state_voltage_stimend (mV)
- **Units**: mV
- **Pseudocode**: ::

    diff_min_duringstim: min_duringstim - steady_state_voltage_stimend

