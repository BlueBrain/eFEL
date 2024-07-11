Changelog
=========
All notable changes to this project will be documented in this file.

The format is based on `Keep a Changelog <https://keepachangelog.com/en/1.0.0/>`_,
and this project adheres to `Semantic Versioning <https://semver.org/spec/v2.0.0.html>`_.

5.7.1 - 2024-07
---------------

- When stimulus starts at 0 for impedance (which can happen often), use steady_state_voltage_stimend to get holding voltage, and equivalent for holding current.
  Since at the end of the stimulus, we expect the current and voltage to vary rapidly around the holding value, this is a good enough proxy to get it.
- implemented test_steady_state_current_stimend, along with documentation and test
- fixed current interpolation
- fixed stim_end value in impedance test data file

5.7.0 - 2024-07
----------------

- Adding extracellular features: `"peak_to_valley", "halfwidth", "peak_trough_ratio", "repolarization_slope", "recovery_slope","neg_peak_relative", "pos_peak_relative", "neg_peak_diff", "pos_peak_diff", "neg_image"` and `"pos_image"`
- Add example for extracellular features
- Documentation update to mention features that can be used fo phase analysis: `"AP_phaseslope", "phaseslope_max"`, `"AP_fall_rate", "AP_fall_rate_change", "AP_peak_downstroke", "AP_peak_upstroke", "AP_rise_rate"` and `"AP_rise_rate_change"``

5.6.29 - 2024-06
----------------

- Added AP_phaseslope_range to efel.settings as int with default value of 2

5.6.27 - 2024-05
----------------

- Fixed edge case (1 spike and no min_AHP_indices dependency) in AP_begin_indices

5.6.20 - 2024-05
----------------

- Refactored LibV1, LibV2, LibV3, LibV5 into more meaningful categories: SpikeEvent, SpikeShape and Subthreshold.

5.6.6 - 2024-04
---------------

- Adding AP_height to documentation

5.6.0 - 2024-02
----------------

- Reduce 3 alternative implementations to get ISIs into 1.
- "all_ISI_values" is recommended, "ISI_values" and "ISIs" are deprecated.
- The features depending on "ISI_values" are moved to Python and now they depend on "all_ISI_values".
- BUGFIX: single_burst_ratio, irregularity_index, burst_mean_freq, interburst_voltage features were ignoring the first two ISIs when the ignore_first_ISI was set.
- Added new feature: inv_ISI_values that computes and returns all of the inverse isi values.

5.5.5 - 2024-01
----------------
- Type annotate api.py's functions.
- Deprecate camel case function names in api.py.
- Start using same requirements_docs.txt in readthedocs and tox.
- Enable autodoc and typehints in the API documentation.
- Fix docstring errors in the io module.
- Add changelog to the documentation.

[5.5.4] - 2024-01
-----------------
- New feature: phaseslope_max

5.5.3 - 2024-01
----------------
- Add type stub for cppcore module to make Python recognise the C++ functions' arguments and return values.

5.5.0 - 2024-01
----------------
C++ changes
^^^^^^^^^^^
- AP_end_indices, AP_rise_time, AP_fall_time, AP_rise_rate, AP_fall_rate do not take into account peaks before stim_start anymore.
- New test and test data for spontaneous firing case. The data is provided by github user SzaBoglarka using cell `https://modeldb.science/114047 <https://modeldb.science/114047>`_.

5.4.0 - 2024-01
----------------
C++ changes
^^^^^^^^^^^
- New C++ function `getFeatures` replaced `getVec`.
- `getFeatures` automatically handles failures & distinguishes empty results from failures.
- Centralized error handling in `getFeatures` shortens the code by removing repetitions.
- C++ features' access is restricted. Read-only references are marked `const`.
- Removed wildcard features from C++ API. Use of Python is encouraged for that purpose.

Python changes
^^^^^^^^^^^^^^
- `bpap_attenuation` feature is added to the Python API.
- `Spikecount`, `Spikecount_stimint`, `burst_number`, `strict_burst_number` and `trace_check` features migrated to Python from C++.
- `check_ais_initiation` is added to the Python API.
