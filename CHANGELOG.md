
# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [5.5.3] - 2024-01

- Add type stub for cppcore module to make Python recognise the C++ functions' arguments and return values.

## [5.5.0] - 2024-01

### C++ changes
- AP_end_indices, AP_rise_time, AP_fall_time, AP_rise_rate, AP_fall_rate do not take into account peaks before stim_start anymore
- New test and test data for spontaneous firing case.
  The data is provided by github user SzaBoglarka using cell https://modeldb.science/114047


## [5.4.0] - 2024-01

### C++ changes
- New C++ function `getFeatures` replaced `getVec`.
- `getFeatures` automatically handles failures & distinguishes empty results from failures.
- Centralized error handling in `getFeatures` shortens the code by removing repetitions.
- C++ features' access is restricted. Read-only references are marked `const`.
- Removed wildcard features from C++ API. Use of Python is encouraged for that purpose.

### Python changes
- `bpap_attenuation` feature is added to the Python API.
- `Spikecount`, `Spikecount_stimint`, `burst_number`, `strict_burst_number` and `trace_check` features migrated to Python from C++.
- `check_ais_initiation` is added to the Python API.
