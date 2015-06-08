eFeature descriptions
=====================

A pdf document describing the eFeatures is available 
`here <http://bluebrain.github.io/eFEL/efeature-documentation.pdf>`_. 

Not every eFeature has a description in this document yet, 
the complete set will be available shortly.

Requested features
==================

**LibV5 : steady_state_voltage_stimend**

The average voltage during the last 90% of the stimulus duration.

- **Required features**: t, V, stim_start, stim_end
- **Units**: mV
- **Pseudocode**: ::

    stim_duration = stim_end[0] - stim_start[0]
    begin_time = stim_end[0] - 0.9 * stim_duration
    end_time = stim_end[0]
    steady_state_voltage_stimend = [numpy.mean(voltage[numpy.where((t < end_time) & (t > end_time))])]


