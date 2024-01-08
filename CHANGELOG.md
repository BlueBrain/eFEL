
# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).


## [5.4.0] - 2024-01

### C++ changes
- New C++ function `getFeatures` replaced `getVec`.
- `getFeatures` automatically handles failures & distinguishes empty results from failures.
- Centralized error handling in `getFeatures` shortens the code by removing repetitions.
- C++ features' access is restricted. Read-only references are marked `const`.
- Removed wildcard features from C++ API. Use of Python is encouraged for that purpose.

### Python changes
- `bpap_attenuation` feature is added to the Python API.
- `Spikecount`, `Spikecount_stimint`, `burst_number` and `strict_burst_number` features migrated to Python from C++.
- `check_ais_initiation` is added to the Python API.