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
  FptrTableV1["peak_indices"] = &LibV1::peak_indices;
  FptrTableV1["ISI_values"] = &LibV1::ISI_values;
  FptrTableV1["peak_voltage"] = &LibV1::peak_voltage;
  FptrTableV1["mean_frequency"] = &LibV1::firing_rate;
  FptrTableV1["peak_time"] = &LibV1::peak_time;
  FptrTableV1["time_to_first_spike"] = &LibV1::first_spike_time;
  FptrTableV1["min_AHP_indices"] = &LibV1::min_AHP_indices;
  FptrTableV1["min_AHP_values"] = &LibV1::min_AHP_values;
  FptrTableV1["voltage_base"] = &LibV1::rest_voltage_value;
  FptrTableV1["burst_ISI_indices"] = &LibV1::burst_ISI_indices;
  FptrTableV1["adaptation_index"] = &LibV1::adaptation_index;
  FptrTableV1["trace_check"] = &LibV1::trace_check;
  FptrTableV1["spike_half_width"] = &LibV1::spike_width1;
  FptrTableV1["spike_width2"] = &LibV1::spike_width2;
  FptrTableV1["burst_mean_freq"] = &LibV1::burst_mean_freq;
  FptrTableV1["interburst_voltage"] = &LibV1::interburst_voltage;
  FptrTableV1["AHP_depth_abs"] = &LibV1::AHP_depth_abs;
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
  FptrTableV1["Spikecount"] = &LibV1::Spikecount;
  FptrTableV1["AHP_depth"] = &LibV1::AHP_depth;
  FptrTableV1["burst_number"] = &LibV1::burst_number;
  FptrTableV1["AP_amplitude_diff"] = &LibV1::AP_amplitude_diff;
  FptrTableV1["AHP_depth_diff"] = &LibV1::AHP_depth_diff;
  FptrTableV1["threshold_current"] = &LibV1::threshold_current;

  //****************** for FptrTableV1 *****************************

  //****************** for FptrTableV2 *****************************
  /*
  FptrTableV2["peak_indices"]         = &LibV2::peak_indices;
  FptrTableV2["ISI_values"]           = &LibV2::ISI_values;
  FptrTableV2["peak_voltage"]         = &LibV2::peak_voltage;
  FptrTableV2["firing_rate"]          = &LibV2::firing_rate;
  FptrTableV2["peak_time"]    = &LibV2::peak_time;
  FptrTableV2["first_spike_time"]     = &LibV2::first_spike_time;
  FptrTableV2["AHP_min_indices"]      = &LibV2::AHP_min_indices;
  FptrTableV2["AHP_values"]           = &LibV2::AHP_values;
  FptrTableV2["rest_voltage_value"]   = &LibV2::rest_voltage_value;
  FptrTableV2["burst_vector"]         = &LibV2::burst_vector;
  FptrTableV2["adaptation_index"]     = &LibV2::adaptation_index;
  */

  // AP parameter
  FptrTableV2["AP_begin_indices"] = &LibV2::AP_begin_indices;
  FptrTableV2["AP_rise_indices"] = &LibV2::AP_rise_indices;
  FptrTableV2["AP_end_indices"] = &LibV2::AP_end_indices;
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
  FptrTableV2["E6"] = &LibV2::E6;
  FptrTableV2["E7"] = &LibV2::E7;

  FptrTableV2["BPAPatt2"] = &LibV2::BPAPatt2;
  FptrTableV2["BPAPatt3"] = &LibV2::BPAPatt3;
  FptrTableV2["E39"] = &LibV2::E39;
  FptrTableV2["E39_cod"] = &LibV2::E39_cod;
  FptrTableV2["E2"] = &LibV2::E2;
  FptrTableV2["E3"] = &LibV2::E3;
  FptrTableV2["E4"] = &LibV2::E4;
  FptrTableV2["E5"] = &LibV2::E5;
  FptrTableV2["E8"] = &LibV2::E8;
  FptrTableV2["E9"] = &LibV2::E9;
  FptrTableV2["E10"] = &LibV2::E10;
  FptrTableV2["E11"] = &LibV2::E11;
  FptrTableV2["E12"] = &LibV2::E12;
  FptrTableV2["E13"] = &LibV2::E13;
  FptrTableV2["E14"] = &LibV2::E14;
  FptrTableV2["E15"] = &LibV2::E15;
  FptrTableV2["E16"] = &LibV2::E16;
  FptrTableV2["E17"] = &LibV2::E17;
  FptrTableV2["E18"] = &LibV2::E18;
  FptrTableV2["E19"] = &LibV2::E19;
  FptrTableV2["E20"] = &LibV2::E20;
  FptrTableV2["E21"] = &LibV2::E21;
  FptrTableV2["E22"] = &LibV2::E22;
  FptrTableV2["E23"] = &LibV2::E23;
  FptrTableV2["E24"] = &LibV2::E24;
  FptrTableV2["E25"] = &LibV2::E25;
  FptrTableV2["E26"] = &LibV2::E26;
  FptrTableV2["E27"] = &LibV2::E27;
  FptrTableV2["E40"] = &LibV2::E40;
  FptrTableV2["steady_state_hyper"] = &LibV2::steady_state_hyper;
  FptrTableV2[string("amp_drop_first_second")] = &LibV2::amp_drop_first_second;
  FptrTableV2[string("amp_drop_first_last")] = &LibV2::amp_drop_first_last;
  FptrTableV2[string("amp_drop_second_last")] = &LibV2::amp_drop_second_last;
  FptrTableV2[string("max_amp_difference")] = &LibV2::max_amp_difference;
  // end of feature definition
  //****************** end of FptrTableV2 *****************************

  //******************  FptrTableV3 *****************************
  FptrTableV3["interpolate"] = &LibV3::interpolate;
  FptrTableV3["trace_check"] = &LibV3::trace_check;
  FptrTableV3["peak_indices"] = &LibV3::peak_indices;
  FptrTableV3["ISI_values"] = &LibV3::ISI_values;
  FptrTableV3["ISI_CV"] = &LibV3::ISI_CV;
  FptrTableV3["peak_voltage"] = &LibV3::peak_voltage;
  FptrTableV3["mean_frequency"] = &LibV3::firing_rate;
  FptrTableV3["peak_time"] = &LibV3::peak_time;
  FptrTableV3["time_to_first_spike"] = &LibV3::first_spike_time;
  FptrTableV3["spike_half_width"] = &LibV3::spike_width1;
  FptrTableV3["min_AHP_indices"] = &LibV3::min_AHP_indices;
  FptrTableV3["min_AHP_values"] = &LibV3::min_AHP_values;
  FptrTableV3["AHP_depth_abs"] = &LibV3::AHP_depth_abs;
  FptrTableV3["voltage_base"] = &LibV3::rest_voltage_value;
  FptrTableV3["adaptation_index2"] = &LibV3::adaptation_index2;
  FptrTableV3["AP_height"] = &LibV3::AP_height;
  FptrTableV3["AP_amplitude"] = &LibV3::AP_amplitude;
  FptrTableV3["AP_width"] = &LibV3::AP_width;
  FptrTableV3["doublet_ISI"] = &LibV3::doublet_ISI;

  FptrTableV3["AP_begin_indices"] = &LibV3::AP_begin_indices;
  FptrTableV3["AP_rise_indices"] = &LibV3::AP_rise_indices;
  FptrTableV3["AP_end_indices"] = &LibV3::AP_end_indices;
  FptrTableV3["AP_fall_indices"] = &LibV3::AP_fall_indices;
  // eFeatures
  FptrTableV3["AP_duration"] = &LibV3::AP_duration;

  FptrTableV3["depolarized_base"] = &LibV3::depolarized_base;

  //****************** end of FptrTableV3 *****************************

  FptrTableV4["peak_indices"] = &LibV4::peak_indices;

  //******************  FptrTableV5 *****************************

  FptrTableV5["ISI_log_slope"] = &LibV5::ISI_log_slope;
  FptrTableV5["ISI_semilog_slope"] = &LibV5::ISI_semilog_slope;
  FptrTableV5["ISI_log_slope_skip"] = &LibV5::ISI_log_slope_skip;
  FptrTableV5["time_to_second_spike"] = &LibV5::time_to_second_spike;
  FptrTableV5["time_to_last_spike"] = &LibV5::time_to_last_spike;
  FptrTableV5["inv_first_ISI"] = &LibV5::inv_first_ISI;
  FptrTableV5["inv_second_ISI"] = &LibV5::inv_second_ISI;
  FptrTableV5["inv_third_ISI"] = &LibV5::inv_third_ISI;
  FptrTableV5["inv_fourth_ISI"] = &LibV5::inv_fourth_ISI;
  FptrTableV5["inv_fifth_ISI"] = &LibV5::inv_fifth_ISI;
  FptrTableV5["inv_last_ISI"] = &LibV5::inv_last_ISI;
  FptrTableV5["inv_time_to_first_spike"] = &LibV5::inv_time_to_first_spike;
  FptrTableV5["min_AHP_indices"] = &LibV5::min_AHP_indices;
  FptrTableV5["min_AHP_values"] = &LibV5::min_AHP_values;
  FptrTableV5["AHP_depth_abs"] = &LibV5::AHP_depth_abs;
  FptrTableV5["spike_half_width"] = &LibV5::spike_width1;
  FptrTableV5["AP_begin_indices"] = &LibV5::AP_begin_indices;
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

  FptrTableV5["BPAPHeightLoc1"] = &LibV5::BPAPHeightLoc1;
  FptrTableV5["BPAPAmplitudeLoc1"] = &LibV5::BPAPAmplitudeLoc1;
  FptrTableV5["BPAPAmplitudeLoc2"] = &LibV5::BPAPAmplitudeLoc2;
  FptrTableV5["BPAPHeightLoc2"] = &LibV5::BPAPHeightLoc2;

  FptrTableV5["check_AISInitiation"] = &LibV5::check_AISInitiation;
  FptrTableV5["AP_phaseslope"] = &LibV5::AP_phaseslope;
  FptrTableV5["AP_phaseslope_AIS"] = &LibV5::AP_phaseslope_AIS;

  FptrTableV5["BAC_width"] = &LibV5::BAC_width;
  FptrTableV5["BAC_maximum_voltage"] = &LibV5::BAC_maximum_voltage;

  FptrTableV5["all_ISI_values"] = &LibV5::all_ISI_values;

  FptrTableV5["AP_amplitude_from_voltagebase"] =
      &LibV5::AP_amplitude_from_voltagebase;
  FptrTableV5["min_voltage_between_spikes"] =
      &LibV5::min_voltage_between_spikes;
  FptrTableV5["voltage"] = &LibV5::voltage;
  FptrTableV5["steady_state_voltage_stimend"] =
      &LibV5::steady_state_voltage_stimend;
  FptrTableV5["voltage_base"] = &LibV5::voltage_base;
  FptrTableV5["decay_time_constant_after_stim"] =
      &LibV5::decay_time_constant_after_stim;

  FptrTableV5["ohmic_input_resistance_vb_ssse"] =
      &LibV5::ohmic_input_resistance_vb_ssse;
  FptrTableV5["voltage_deflection_vb_ssse"] =
      &LibV5::voltage_deflection_vb_ssse;
  FptrTableV5["maximum_voltage_from_voltagebase"] =
      &LibV5::maximum_voltage_from_voltagebase;
  
  //****************** end of FptrTableV5 *****************************

  return 1;
}
