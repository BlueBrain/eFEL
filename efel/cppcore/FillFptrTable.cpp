/* Copyright (c) 2015, EPFL/Blue Brain Project
 *
 * This file is part of eFEL <https://github.com/BlueBrain/eFEL>
 *
 * This library is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License version 3.0 as published
 * by the Free Software Foundation.
 *
 * This library is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this library; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */
#include "FillFptrTable.h"

int FillFptrTable() {
  //****************** for FptrTableV1 *****************************
  FptrTableV1["interpolate"] = &LibV1::interpolate;
  FptrTableV1["ISI_values"] = &LibV1::ISI_values;
  FptrTableV1["peak_voltage"] = &LibV1::peak_voltage;
  FptrTableV1["mean_frequency"] = &LibV1::firing_rate;
  FptrTableV1["peak_time"] = &LibV1::peak_time;
  FptrTableV1["time_to_first_spike"] = &LibV1::first_spike_time;
  FptrTableV1["burst_ISI_indices"] = &LibV1::burst_ISI_indices;
  FptrTableV1["adaptation_index"] = &LibV1::adaptation_index;
  FptrTableV1["spike_width2"] = &LibV1::spike_width2;
  FptrTableV1["burst_mean_freq"] = &LibV1::burst_mean_freq;
  FptrTableV1["interburst_voltage"] = &LibV1::interburst_voltage;
  // passive properties
  FptrTableV1["time_constant"] = &LibV1::time_constant;
  FptrTableV1["voltage_deflection"] = &LibV1::voltage_deflection;
  FptrTableV1["ohmic_input_resistance"] = &LibV1::ohmic_input_resistance;
  FptrTableV1["maximum_voltage"] = &LibV1::maximum_voltage;
  FptrTableV1["minimum_voltage"] = &LibV1::minimum_voltage;
  FptrTableV1["steady_state_voltage"] = &LibV1::steady_state_voltage;

  FptrTableV1["AP_height"] = &LibV1::AP_height;
  FptrTableV1["AP_amplitude"] = &LibV1::AP_amplitude;
  FptrTableV1["single_burst_ratio"] = &LibV1::single_burst_ratio;
  FptrTableV1["AP_width"] = &LibV1::AP_width;
  FptrTableV1["doublet_ISI"] = &LibV1::doublet_ISI;
  FptrTableV1["adaptation_index2"] = &LibV1::adaptation_index2;
  FptrTableV1["ISI_CV"] = &LibV1::ISI_CV;
  FptrTableV1["AHP_depth_abs_slow"] = &LibV1::AHP_depth_abs_slow;
  FptrTableV1["AHP_slow_time"] = &LibV1::AHP_slow_time;
  FptrTableV1["AHP_depth"] = &LibV1::AHP_depth;
  FptrTableV1["AHP_depth_slow"] = &LibV1::AHP_depth_slow;
  FptrTableV1["AP_amplitude_diff"] = &LibV1::AP_amplitude_diff;
  FptrTableV1["AHP_depth_diff"] = &LibV1::AHP_depth_diff;

  //****************** for FptrTableV1 *****************************

  //****************** for FptrTableV2 *****************************
  /*
  FptrTableV2["AHP_min_indices"]      = &LibV2::AHP_min_indices;
  FptrTableV2["AHP_values"]           = &LibV2::AHP_values;
  FptrTableV2["burst_vector"]         = &LibV2::burst_vector;
  */

  // AP parameter
  FptrTableV2["AP_rise_indices"] = &LibV2::AP_rise_indices;
  FptrTableV2["AP_fall_indices"] = &LibV2::AP_fall_indices;
  // eFeatures
  FptrTableV2["AP_duration"] = &LibV2::AP_duration;
  FptrTableV2["AP_rise_time"] = &LibV2::AP_rise_time;
  FptrTableV2["AP_fall_time"] = &LibV2::AP_fall_time;
  FptrTableV2["AP_rise_rate"] = &LibV2::AP_rise_rate;
  FptrTableV2["AP_fall_rate"] = &LibV2::AP_fall_rate;
  FptrTableV2["fast_AHP"] = &LibV2::fast_AHP;
  FptrTableV2["AP_amplitude_change"] = &LibV2::AP_amplitude_change;
  FptrTableV2["AP_duration_change"] = &LibV2::AP_duration_change;
  FptrTableV2["AP_rise_rate_change"] = &LibV2::AP_rise_rate_change;
  FptrTableV2["AP_fall_rate_change"] = &LibV2::AP_fall_rate_change;
  FptrTableV2["fast_AHP_change"] = &LibV2::fast_AHP_change;
  FptrTableV2["AP_duration_half_width"] = &LibV2::AP_duration_half_width;
  FptrTableV2["AP_duration_half_width_change"] =
      &LibV2::AP_duration_half_width_change;
  FptrTableV2["steady_state_hyper"] = &LibV2::steady_state_hyper;
  FptrTableV2[string("amp_drop_first_second")] = &LibV2::amp_drop_first_second;
  FptrTableV2[string("amp_drop_first_last")] = &LibV2::amp_drop_first_last;
  FptrTableV2[string("amp_drop_second_last")] = &LibV2::amp_drop_second_last;
  FptrTableV2[string("max_amp_difference")] = &LibV2::max_amp_difference;
  // end of feature definition
  //****************** end of FptrTableV2 *****************************

  //******************  FptrTableV3 *****************************
  // eFeatures
  FptrTableV3["depolarized_base"] = &LibV3::depolarized_base;

  //****************** end of FptrTableV3 *****************************

  //******************  FptrTableV5 *****************************

  FptrTableV5["ISI_log_slope"] = &LibV5::ISI_log_slope;
  FptrTableV5["ISI_semilog_slope"] = &LibV5::ISI_semilog_slope;
  FptrTableV5["ISI_log_slope_skip"] = &LibV5::ISI_log_slope_skip;
  FptrTableV5["time_to_second_spike"] = &LibV5::time_to_second_spike;
  FptrTableV5["time_to_last_spike"] = &LibV5::time_to_last_spike;
  FptrTableV5["inv_first_ISI"] = [](mapStr2intVec& intData,
                                    mapStr2doubleVec& doubleData,
                                    mapStr2Str& strData) {
    return LibV5::inv_ISI_generic(intData, doubleData, strData, 0);
  };
  FptrTableV5["inv_second_ISI"] = [](mapStr2intVec& intData,
                                     mapStr2doubleVec& doubleData,
                                     mapStr2Str& strData) {
    return LibV5::inv_ISI_generic(intData, doubleData, strData, 1);
  };
  FptrTableV5["inv_third_ISI"] = [](mapStr2intVec& intData,
                                    mapStr2doubleVec& doubleData,
                                    mapStr2Str& strData) {
    return LibV5::inv_ISI_generic(intData, doubleData, strData, 2);
  };
  FptrTableV5["inv_fourth_ISI"] = [](mapStr2intVec& intData,
                                     mapStr2doubleVec& doubleData,
                                     mapStr2Str& strData) {
    return LibV5::inv_ISI_generic(intData, doubleData, strData, 3);
  };
  FptrTableV5["inv_fifth_ISI"] = [](mapStr2intVec& intData,
                                    mapStr2doubleVec& doubleData,
                                    mapStr2Str& strData) {
    return LibV5::inv_ISI_generic(intData, doubleData, strData, 4);
  };
  FptrTableV5["inv_last_ISI"] = &LibV5::inv_last_ISI;
  FptrTableV5["inv_time_to_first_spike"] = &LibV5::inv_time_to_first_spike;
  FptrTableV5["min_AHP_indices"] = &LibV5::min_AHP_indices;
  FptrTableV5["min_AHP_values"] = &LibV5::min_AHP_values;
  FptrTableV5["AHP_depth_abs"] = &LibV5::AHP_depth_abs;
  FptrTableV5["spike_half_width"] = &LibV5::spike_width1;
  FptrTableV5["AP_begin_indices"] = &LibV5::AP_begin_indices;
  FptrTableV5["AP_end_indices"] = &LibV5::AP_end_indices;
  FptrTableV5["irregularity_index"] = &LibV5::irregularity_index;
  FptrTableV5["number_initial_spikes"] = &LibV5::number_initial_spikes;
  FptrTableV5["AP1_amp"] = &LibV5::AP1_amp;
  FptrTableV5["APlast_amp"] = &LibV5::APlast_amp;
  FptrTableV5["AP1_peak"] = &LibV5::AP1_peak;
  FptrTableV5["AP1_width"] = &LibV5::AP1_width;
  FptrTableV5["AP2_amp"] = &LibV5::AP2_amp;
  FptrTableV5["AP2_peak"] = &LibV5::AP2_peak;
  FptrTableV5["AP2_AP1_diff"] = &LibV5::AP2_AP1_diff;
  FptrTableV5["AP2_AP1_peak_diff"] = &LibV5::AP2_AP1_peak_diff;
  FptrTableV5["AP2_width"] = &LibV5::AP2_width;
  FptrTableV5["APlast_width"] = &LibV5::APlast_width;
  FptrTableV5["AHP_depth_from_peak"] = &LibV5::AHP_depth_from_peak;
  FptrTableV5["AHP_time_from_peak"] = &LibV5::AHP_time_from_peak;
  FptrTableV5["AHP1_depth_from_peak"] = &LibV5::AHP1_depth_from_peak;
  FptrTableV5["AHP2_depth_from_peak"] = &LibV5::AHP2_depth_from_peak;
  FptrTableV5["AP_begin_width"] = &LibV5::AP_begin_width;
  FptrTableV5["AP_begin_time"] = &LibV5::AP_begin_time;
  FptrTableV5["AP_begin_voltage"] = &LibV5::AP_begin_voltage;

  FptrTableV5["AP1_begin_width"] = &LibV5::AP1_begin_width;
  FptrTableV5["AP2_begin_width"] = &LibV5::AP2_begin_width;

  FptrTableV5["AP1_begin_voltage"] = &LibV5::AP1_begin_voltage;
  FptrTableV5["AP2_begin_voltage"] = &LibV5::AP2_begin_voltage;

  FptrTableV5["voltage_deflection_begin"] = &LibV5::voltage_deflection_begin;
  FptrTableV5["is_not_stuck"] = &LibV5::is_not_stuck;
  FptrTableV5["mean_AP_amplitude"] = &LibV5::mean_AP_amplitude;
  FptrTableV5["voltage_after_stim"] = &LibV5::voltage_after_stim;

  FptrTableV5["AP2_AP1_begin_width_diff"] = &LibV5::AP2_AP1_begin_width_diff;

  FptrTableV5["AP_phaseslope"] = &LibV5::AP_phaseslope;

  FptrTableV5["all_ISI_values"] = &LibV5::all_ISI_values;

  FptrTableV5["AP_amplitude_from_voltagebase"] =
      &LibV5::AP_amplitude_from_voltagebase;
  FptrTableV5["min_voltage_between_spikes"] =
      &LibV5::min_voltage_between_spikes;
  FptrTableV5["voltage"] = &LibV5::voltage;
  FptrTableV5["current"] = &LibV5::current;
  FptrTableV5["time"] = &LibV5::time;
  FptrTableV5["steady_state_voltage_stimend"] =
      &LibV5::steady_state_voltage_stimend;
  FptrTableV5["voltage_base"] = &LibV5::voltage_base;
  FptrTableV5["current_base"] = &LibV5::current_base;
  FptrTableV5["decay_time_constant_after_stim"] =
      &LibV5::decay_time_constant_after_stim;
  FptrTableV5["multiple_decay_time_constant_after_stim"] =
      &LibV5::multiple_decay_time_constant_after_stim;
  FptrTableV5["sag_time_constant"] = &LibV5::sag_time_constant;

  FptrTableV5["ohmic_input_resistance_vb_ssse"] =
      &LibV5::ohmic_input_resistance_vb_ssse;
  FptrTableV5["voltage_deflection_vb_ssse"] =
      &LibV5::voltage_deflection_vb_ssse;
  FptrTableV5["maximum_voltage_from_voltagebase"] =
      &LibV5::maximum_voltage_from_voltagebase;

  FptrTableV5["peak_indices"] = &LibV5::peak_indices;
  FptrTableV5["sag_amplitude"] = &LibV5::sag_amplitude;
  FptrTableV5["sag_ratio1"] = &LibV5::sag_ratio1;
  FptrTableV5["sag_ratio2"] = &LibV5::sag_ratio2;
  FptrTableV5["AP_peak_upstroke"] = &LibV5::AP_peak_upstroke;
  FptrTableV5["AP_peak_downstroke"] = &LibV5::AP_peak_downstroke;
  FptrTableV5["min_between_peaks_indices"] = &LibV5::min_between_peaks_indices;
  FptrTableV5["min_between_peaks_values"] = &LibV5::min_between_peaks_values;
  FptrTableV5["AP_width_between_threshold"] =
      &LibV5::AP_width_between_threshold;

  FptrTableV5["burst_begin_indices"] = &LibV5::burst_begin_indices;
  FptrTableV5["burst_end_indices"] = &LibV5::burst_end_indices;
  FptrTableV5["strict_burst_mean_freq"] = &LibV5::strict_burst_mean_freq;
  FptrTableV5["strict_interburst_voltage"] = &LibV5::strict_interburst_voltage;

  FptrTableV5["ADP_peak_indices"] = &LibV5::ADP_peak_indices;
  FptrTableV5["ADP_peak_values"] = &LibV5::ADP_peak_values;
  FptrTableV5["ADP_peak_amplitude"] = &LibV5::ADP_peak_amplitude;

  FptrTableV5["interburst_min_indices"] = &LibV5::interburst_min_indices;
  FptrTableV5["interburst_min_values"] = &LibV5::interburst_min_values;
  FptrTableV5["postburst_min_indices"] = &LibV5::postburst_min_indices;
  FptrTableV5["postburst_min_values"] = &LibV5::postburst_min_values;
  FptrTableV5["time_to_interburst_min"] = &LibV5::time_to_interburst_min;

  //****************** end of FptrTableV5 *****************************

  return 1;
}
