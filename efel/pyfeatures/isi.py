"""Features that are depending on the inter-spike intervals."""
from __future__ import annotations
import numpy as np
from efel.pyfeatures.cppfeature_access import _get_cpp_data, get_cpp_feature


def ISIs() -> np.ndarray | None:
    """Get all ISIs, inter-spike intervals."""
    peak_times = get_cpp_feature("peak_time")
    if peak_times is None:
        return None
    else:
        return np.diff(peak_times)


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
    isi_values = ISIs()
    if isi_values is None:
        return None

    # Check "ignore_first_ISI" flag
    ignore_first_ISI = _get_cpp_data("ignore_first_ISI")
    if ignore_first_ISI:
        isi_values = isi_values[1:]

    result = __ISI_CV(isi_values)
    return np.array([result])


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
    all_isis = ISIs()

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
