"""Features that are depending on the inter-spike intervals."""
from __future__ import annotations
import logging
from typing_extensions import deprecated
import warnings
import numpy as np
from efel.pyfeatures.cppfeature_access import (
    _get_cpp_data, get_cpp_feature
)


logger = logging.getLogger(__name__)


@deprecated("Use all_ISI_values instead")
def ISIs() -> np.ndarray | None:
    """Get all ISIs, inter-spike intervals."""
    return get_cpp_feature("all_ISI_values")


@deprecated("Use ISIs instead.")
def ISI_values() -> np.ndarray | None:
    """Get all ISIs, inter-spike intervals."""
    isi_values = get_cpp_feature("all_ISI_values")
    if isi_values is None:
        return None

    # Check "ignore_first_ISI" flag
    ignore_first_ISI = _get_cpp_data("ignore_first_ISI")
    if ignore_first_ISI:
        isi_values = isi_values[1:]
    return isi_values


def __ISI_CV(isi_values) -> float | None:
    if len(isi_values) < 2:
        return None

    # Calculate mean
    isi_mean = np.mean(isi_values)
    # Calculate coefficient of variation
    cv = np.std(isi_values, ddof=1) / isi_mean  # ddof 1 to replicate C++ impl
    return cv


def ISI_CV() -> np.ndarray | None:
    """Coefficient of variation of ISIs.

    If the ignore_first_ISI flag is set, the first ISI will be ignored.
    """
    isi_values = get_cpp_feature("all_ISI_values")
    if isi_values is None:
        return None

    # Check "ignore_first_ISI" flag
    ignore_first_ISI = _get_cpp_data("ignore_first_ISI")
    if ignore_first_ISI:
        isi_values = isi_values[1:]

    result = __ISI_CV(isi_values)
    if result is None:
        return None
    return np.array([result])


def single_burst_ratio() -> np.ndarray | None:
    """Calculates the single burst ratio.

    The ratio is the length of the first ISI over the average of the rest.
    If the ignore_first_ISI flag is set, the first ISI will be ignored.
    """
    isi_values = get_cpp_feature("all_ISI_values")
    if isi_values is None:
        return None

    # Check "ignore_first_ISI" flag
    ignore_first_ISI = _get_cpp_data("ignore_first_ISI")
    if ignore_first_ISI:
        isi_values = isi_values[1:]

    if len(isi_values) < 2:
        return None

    single_burst_ratio_value = isi_values[0] / np.mean(isi_values)
    return np.array([single_burst_ratio_value])


def irregularity_index() -> np.ndarray | None:
    """Calculate the irregularity index of ISI values.

    If the ignore_first_ISI flag is set, the first ISI will be ignored.
    """
    isi_values = get_cpp_feature("all_ISI_values")
    if isi_values is None:
        return None

    # Check "ignore_first_ISI" flag
    ignore_first_ISI = _get_cpp_data("ignore_first_ISI")
    if ignore_first_ISI:
        isi_values = isi_values[1:]

    # Calculate the absolute differences between consecutive ISI values
    isi_differences = np.abs(np.diff(isi_values))
    result = np.mean(isi_differences)

    return np.array([result])


def _isi_log_slope_core(
    isi_values, skip=False, spike_skipf=0.0, max_spike_skip=0, semilog=False
) -> np.ndarray | None:
    if isi_values is None or len(isi_values) == 0:
        return None

    if skip:
        isisToRemove = min(
            max_spike_skip, int((len(isi_values) + 1) * spike_skipf + 0.5)
        )
        isi_values = isi_values[isisToRemove:]

    log_isi_values = np.log(isi_values)
    x = np.arange(1, len(log_isi_values) + 1)
    if not semilog:
        x = np.log(x)

    try:
        slope, _ = np.polyfit(x, log_isi_values, 1)
    except (np.linalg.LinAlgError, SystemError) as e:
        warnings.warn(f"Error in polyfit: {e}")
        return None

    return np.array([slope])


def ISI_log_slope() -> np.ndarray | None:
    """The slope of a linear fit to a loglog plot of the ISI values.

    If the ignore_first_ISI flag is set, the first ISI will be ignored.
    """
    isi_values = get_cpp_feature("all_ISI_values")
    if isi_values is None:
        return None

    # Check "ignore_first_ISI" flag
    ignore_first_ISI = _get_cpp_data("ignore_first_ISI")
    if ignore_first_ISI:
        isi_values = isi_values[1:]

    return _isi_log_slope_core(isi_values, False, 0.0, 0, False)


def ISI_semilog_slope() -> np.ndarray | None:
    """The slope of a linear fit to a semilog plot of the ISI values.

    If the ignore_first_ISI flag is set, the first ISI will be ignored.
    """
    isi_values = get_cpp_feature("all_ISI_values")
    if isi_values is None:
        return None

    # Check "ignore_first_ISI" flag
    ignore_first_ISI = _get_cpp_data("ignore_first_ISI")
    if ignore_first_ISI:
        isi_values = isi_values[1:]

    return _isi_log_slope_core(isi_values, False, 0.0, 0, True)


def ISI_log_slope_skip() -> np.ndarray | None:
    """The slope of a linear fit to a loglog plot of the ISI values,
    but not taking into account the first ISI values.

    Uses the spike_skipf and max_spike_skip settings to determine how many
    ISIs to skip.
    ."""
    isi_values = get_cpp_feature("all_ISI_values")
    if isi_values is None:
        return None
    # Check "ignore_first_ISI" flag
    ignore_first_ISI = _get_cpp_data("ignore_first_ISI")
    if ignore_first_ISI:
        isi_values = isi_values[1:]

    spike_skipf = _get_cpp_data("spike_skipf")
    if spike_skipf < 0 or spike_skipf >= 1:
        raise ValueError("spike_skipf should lie between [0, 1).")
    max_spike_skip = _get_cpp_data("max_spike_skip")
    return _isi_log_slope_core(isi_values, True, spike_skipf, max_spike_skip, False)


def burst_ISI_indices() -> np.ndarray | None:
    """Calculate burst ISI indices based on burst factor and ISI values."""
    # Fetching necessary data
    isi_values = ISI_values()
    if isi_values is None:
        return None

    burst_factor = _get_cpp_data("burst_factor")

    if len(isi_values) < 4:
        logger.warning("4 or more spikes are needed for burst calculation.")
        return None

    burst_indices = []
    count = -1

    for i in range(1, len(isi_values) - 1):
        isi_p_copy = isi_values[count + 1: i]
        n = len(isi_p_copy)

        if n == 0:
            continue

        # Compute median
        d_median = np.median(isi_p_copy)

        # Check burst condition
        if isi_values[i] > (burst_factor * d_median) and isi_values[i + 1] < (
            isi_values[i] / burst_factor
        ):
            burst_indices.append(i + 1)
            count = i - 1

    if burst_indices == []:
        return None
    return np.array(burst_indices)


def burst_mean_freq() -> np.ndarray | None:
    """Calculate the mean frequency of bursts."""
    # Fetching required features
    peak_time = get_cpp_feature("peak_time")
    burst_isi_idx = burst_ISI_indices()
    if burst_isi_idx is None or peak_time is None:
        return None

    burst_mean_freq = []
    burst_index_tmp = burst_isi_idx
    burst_index = np.insert(
        burst_index_tmp, burst_index_tmp.size, len(peak_time) - 1
    )
    burst_index = burst_index.astype(int)

    # 1st burst
    span = peak_time[burst_index[0]] - peak_time[0]
    # + 1 because 1st ISI is ignored
    N_peaks = burst_index[0] + 1
    burst_mean_freq.append(N_peaks * 1000 / span)

    for i, burst_idx in enumerate(burst_index[:-1]):
        if burst_index[i + 1] - burst_idx != 1:
            span = peak_time[burst_index[i + 1]] - peak_time[burst_idx + 1]
            N_peaks = burst_index[i + 1] - burst_idx
            burst_mean_freq.append(N_peaks * 1000 / span)

    return np.array(burst_mean_freq)


def interburst_voltage() -> np.ndarray | None:
    """The voltage average in between two bursts."""
    peak_idx = get_cpp_feature("peak_indices")
    v = get_cpp_feature("voltage")
    t = get_cpp_feature("time")
    burst_isi_idx = burst_ISI_indices()
    if peak_idx is None or v is None or t is None or burst_isi_idx is None:
        return None

    interburst_voltage = []
    for idx in burst_isi_idx:
        ts_idx = peak_idx[idx]
        t_start = t[ts_idx] + 5
        start_idx = np.argwhere(t < t_start)[-1][0]

        te_idx = peak_idx[idx + 1]
        t_end = t[te_idx] - 5
        end_idx = np.argwhere(t > t_end)[0][0]

        interburst_voltage.append(np.mean(v[start_idx:end_idx + 1]))

    return np.array(interburst_voltage)


def initburst_sahp() -> np.ndarray | None:
    """SlowAHP voltage after initial burst."""
    # Required cpp features
    voltage = get_cpp_feature("voltage")
    time = get_cpp_feature("time")
    if voltage is None or time is None:
        return None
    time = time[:len(voltage)]
    peak_times = get_cpp_feature("peak_time")
    if peak_times is None:
        return None

    # Required python features
    all_isis = get_cpp_feature("all_ISI_values")

    # Required trace data
    stim_end = _get_cpp_data("stim_end")

    # Required settings
    initburst_freq_thresh = _get_cpp_data("initburst_freq_threshold")
    initburst_sahp_start = _get_cpp_data("initburst_sahp_start")
    initburst_sahp_end = _get_cpp_data("initburst_sahp_end")

    last_isi = None

    # Loop over ISIs until frequency higher than initburst_freq_threshold
    if all_isis is None:
        return None
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
        sahp_interval = voltage[np.where(
            (time <= sahp_interval_end) &
            (time >= sahp_interval_start))]

        if len(sahp_interval) > 0:
            min_volt_index = np.argmin(sahp_interval)
        else:
            return None

        slow_ahp = sahp_interval[min_volt_index]

        return np.array([slow_ahp])


def strict_burst_number() -> np.ndarray:
    """Calculate the strict burst number.

    This implementation does not assume that every spike belongs to a burst.
    The first spike is ignored by default. This can be changed by setting
    ignore_first_ISI to 0.

    The burst detection can be fine-tuned by changing the setting
    strict_burst_factor. Default value is 2.0."""
    burst_mean_freq = get_cpp_feature("strict_burst_mean_freq")
    if burst_mean_freq is None:
        return np.array([0])
    return np.array([burst_mean_freq.size])


def inv_ISI_values() -> np.ndarray | None:
    """Calculate the inverse of ISI values."""
    all_isi_values_vec = get_cpp_feature("all_ISI_values")
    return None if all_isi_values_vec is None else 1000.0 / all_isi_values_vec


def inv_first_ISI() -> np.ndarray | None:
    """Calculate the inverse of the first ISI."""
    inv_isi_values = inv_ISI_values()
    if inv_isi_values is None or len(inv_isi_values) == 0:
        return None
    return np.array([inv_isi_values[0]])


def inv_second_ISI() -> np.ndarray | None:
    """Calculate the inverse of the second ISI."""
    inv_isi_values = inv_ISI_values()
    if inv_isi_values is None or len(inv_isi_values) < 2:
        return None
    return np.array([inv_isi_values[1]])


def inv_third_ISI() -> np.ndarray | None:
    """Calculate the inverse of the third ISI."""
    inv_isi_values = inv_ISI_values()
    if inv_isi_values is None or len(inv_isi_values) < 3:
        return None
    return np.array([inv_isi_values[2]])


def inv_fourth_ISI() -> np.ndarray | None:
    """Calculate the inverse of the fourth ISI."""
    inv_isi_values = inv_ISI_values()
    if inv_isi_values is None or len(inv_isi_values) < 4:
        return None
    return np.array([inv_isi_values[3]])


def inv_fifth_ISI() -> np.ndarray | None:
    """Calculate the inverse of the fifth ISI."""
    inv_isi_values = inv_ISI_values()
    if inv_isi_values is None or len(inv_isi_values) < 5:
        return None
    return np.array([inv_isi_values[4]])


def inv_last_ISI() -> np.ndarray | None:
    """Calculate the inverse of the last ISI."""
    inv_isi_values = inv_ISI_values()
    if inv_isi_values is None or len(inv_isi_values) == 0:
        return None
    return np.array([inv_isi_values[-1]])
