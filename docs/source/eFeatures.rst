eFeature descriptions
=====================

A pdf document describing the eFeatures is available
`here <http://bluebrain.github.io/eFEL/efeature-documentation.pdf>`_.

Not every eFeature has a description in this document yet,
the complete set will be available shortly.

eFeatures (to be continued)
========================

    Spike event features
    --------------

        **LibV5 : inv_time_to_first_spike**

        1.0 over time to first spike; returns 0 when no spike

        - **Required features**: time_to_first_spike
        - **Units**: Hz
        - **Pseudocode**: ::

            if len(time_to_first_spike) > 0:
                inv_time_to_first_spike = 1.0 / time_to_first_spike[0]
            else:
                inv_time_to_first_spike = 0


        **LibV5 : inv_first_ISI, inv_second_ISI, inv_third_ISI, inv_fourth_ISI, inv_fith_ISI, inv_last_ISI**

        1.0 over first/second/third/fourth/fith/last ISI; returns 0 when no ISI

        - **Required features**: peak_time (ms)
        - **Units**: Hz
        - **Pseudocode**: ::

            all_isi_values_vec = numpy.diff(peak_time)
            if len(all_isi_values_vec) > 1:
                inv_first_ISI = 1000.0 / all_isi_values_vec[0]
            else:
                inv_first_ISI = 0

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


            **LibV5 : time_to_last_spike**

            time from stimulus start to last spike

            - **Required features**: peak_time (ms), stimstart (ms)
            - **Units**: ms
            - **Pseudocode**: ::

                if len(peak_time) > 0:
                    time_to_last_spike = peak_time[-1] - stimstart
                else:
                    time_to_last_spike = 0


        .. image:: figures/inv_ISI.png


    Spike shape features
    --------------


.. image:: figures/AHP.png
.. image:: figures/AP_Amplitude.png
.. image:: figures/AP_duration_half_width.png
.. image:: figures/voltage_features.png


Recently added eFeatures
========================

**LibV5 : steady_state_voltage_stimend**

The average voltage during the last 90% of the stimulus duration.

- **Required features**: t, V, stim_start, stim_end
- **Units**: mV
- **Pseudocode**: ::

    stim_duration = stim_end[0] - stim_start[0]
    begin_time = stim_end[0] - 0.1 * stim_duration
    end_time = stim_end[0]
    steady_state_voltage_stimend = [numpy.mean(voltage[numpy.where((t <= end_time) & (t >= begin_time))])]




Requested eFeatures
===================
