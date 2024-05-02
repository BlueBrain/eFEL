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
  //****************** for FptrTableBF *****************************
  FptrTableBF["interpolate"] = &BasicFeatures::interpolate;
  FptrTableBF["voltage"] = &BasicFeatures::voltage;
  FptrTableBF["current"] = &BasicFeatures::current;
  FptrTableBF["time"] = &BasicFeatures::time;

  //****************** for FptrTableSE *****************************

  FptrTableSE["peak_indices"] = &SpikeEvent::peak_indices;
  FptrTableSE["peak_time"] = &SpikeEvent::peak_time;
  FptrTableSE["time_to_first_spike"] = &SpikeEvent::first_spike_time;
  FptrTableSE["time_to_second_spike"] = &SpikeEvent::time_to_second_spike;
  FptrTableSE["inv_time_to_first_spike"] = &SpikeEvent::inv_time_to_first_spike;
  FptrTableSE["doublet_ISI"] = &SpikeEvent::doublet_ISI;
  FptrTableSE["all_ISI_values"] = &SpikeEvent::all_ISI_values;
  FptrTableSE["time_to_last_spike"] = &SpikeEvent::time_to_last_spike;
  FptrTableSE["number_initial_spikes"] = &SpikeEvent::number_initial_spikes;
  FptrTableSE["mean_frequency"] = &SpikeEvent::firing_rate;
  FptrTableSE["adaptation_index"] = &SpikeEvent::adaptation_index;
  FptrTableSE["adaptation_index2"] = &SpikeEvent::adaptation_index2;
  FptrTableSE["burst_begin_indices"] = &SpikeEvent::burst_begin_indices;
  FptrTableSE["burst_end_indices"] = &SpikeEvent::burst_end_indices;
  FptrTableSE["strict_burst_mean_freq"] = &SpikeEvent::strict_burst_mean_freq;
  FptrTableSE["strict_interburst_voltage"] = &SpikeEvent::strict_interburst_voltage;
  FptrTableSE["interburst_min_indices"] = &SpikeEvent::interburst_min_indices;
  FptrTableSE["interburst_min_values"] = &SpikeEvent::interburst_min_values;
  FptrTableSE["postburst_min_indices"] = &SpikeEvent::postburst_min_indices;
  FptrTableSE["postburst_min_values"] = &SpikeEvent::postburst_min_values;
  FptrTableSE["time_to_interburst_min"] = &SpikeEvent::time_to_interburst_min;
  FptrTableSE["postburst_slow_ahp_indices"] = &SpikeEvent::postburst_slow_ahp_indices;
  FptrTableSE["postburst_slow_ahp_values"] = &SpikeEvent::postburst_slow_ahp_values;
  FptrTableSE["time_to_postburst_slow_ahp"] = &SpikeEvent::time_to_postburst_slow_ahp;
  FptrTableSE["postburst_fast_ahp_indices"] = &SpikeEvent::postburst_fast_ahp_indices;
  FptrTableSE["postburst_fast_ahp_values"] = &SpikeEvent::postburst_fast_ahp_values;
  FptrTableSE["postburst_adp_peak_indices"] = &SpikeEvent::postburst_adp_peak_indices;
  FptrTableSE["postburst_adp_peak_values"] = &SpikeEvent::postburst_adp_peak_values;
  FptrTableSE["time_to_postburst_fast_ahp"] = &SpikeEvent::time_to_postburst_fast_ahp;
  FptrTableSE["time_to_postburst_adp_peak"] = &SpikeEvent::time_to_postburst_adp_peak;
  FptrTableSE["interburst_15percent_indices"] = [](mapStr2intVec& intData,
                                                   mapStr2doubleVec& doubleData,
                                                   mapStr2Str& strData) {
    return SpikeEvent::interburst_XXpercent_indices(intData, doubleData, strData, 15);
  };
  FptrTableSE["interburst_15percent_values"] = [](mapStr2intVec& intData,
                                                  mapStr2doubleVec& doubleData,
                                                  mapStr2Str& strData) {
    return SpikeEvent::interburst_XXpercent_indices(intData, doubleData, strData, 15);
  };
  FptrTableSE["interburst_20percent_indices"] = [](mapStr2intVec& intData,
                                                   mapStr2doubleVec& doubleData,
                                                   mapStr2Str& strData) {
    return SpikeEvent::interburst_XXpercent_indices(intData, doubleData, strData, 20);
  };
  FptrTableSE["interburst_20percent_values"] = [](mapStr2intVec& intData,
                                                  mapStr2doubleVec& doubleData,
                                                  mapStr2Str& strData) {
    return SpikeEvent::interburst_XXpercent_indices(intData, doubleData, strData, 20);
  };
  FptrTableSE["interburst_25percent_indices"] = [](mapStr2intVec& intData,
                                                   mapStr2doubleVec& doubleData,
                                                   mapStr2Str& strData) {
    return SpikeEvent::interburst_XXpercent_indices(intData, doubleData, strData, 25);
  };
  FptrTableSE["interburst_25percent_values"] = [](mapStr2intVec& intData,
                                                  mapStr2doubleVec& doubleData,
                                                  mapStr2Str& strData) {
    return SpikeEvent::interburst_XXpercent_indices(intData, doubleData, strData, 25);
  };
  FptrTableSE["interburst_30percent_indices"] = [](mapStr2intVec& intData,
                                                   mapStr2doubleVec& doubleData,
                                                   mapStr2Str& strData) {
    return SpikeEvent::interburst_XXpercent_indices(intData, doubleData, strData, 30);
  };
  FptrTableSE["interburst_30percent_values"] = [](mapStr2intVec& intData,
                                                  mapStr2doubleVec& doubleData,
                                                  mapStr2Str& strData) {
    return SpikeEvent::interburst_XXpercent_indices(intData, doubleData, strData, 30);
  };
  FptrTableSE["interburst_40percent_indices"] = [](mapStr2intVec& intData,
                                                   mapStr2doubleVec& doubleData,
                                                   mapStr2Str& strData) {
    return SpikeEvent::interburst_XXpercent_indices(intData, doubleData, strData, 40);
  };
  FptrTableSE["interburst_40percent_values"] = [](mapStr2intVec& intData,
                                                  mapStr2doubleVec& doubleData,
                                                  mapStr2Str& strData) {
    return SpikeEvent::interburst_XXpercent_indices(intData, doubleData, strData, 40);
  };
  FptrTableSE["interburst_60percent_indices"] = [](mapStr2intVec& intData,
                                                   mapStr2doubleVec& doubleData,
                                                   mapStr2Str& strData) {
    return SpikeEvent::interburst_XXpercent_indices(intData, doubleData, strData, 60);
  };
  FptrTableSE["interburst_60percent_values"] = [](mapStr2intVec& intData,
                                                  mapStr2doubleVec& doubleData,
                                                  mapStr2Str& strData) {
    return SpikeEvent::interburst_XXpercent_indices(intData, doubleData, strData, 60);
  };
  FptrTableSE["interburst_duration"] = &SpikeEvent::interburst_duration;
  FptrTableSE["is_not_stuck"] = &SpikeEvent::is_not_stuck;

  //******************  FptrTableSS *****************************
  // eFeatures
  FptrTableSS["peak_voltage"] = &SpikeShape::peak_voltage;
  FptrTableSS["AP_height"] = &SpikeShape::AP_height;
  FptrTableSS["AP_amplitude"] = &SpikeShape::AP_amplitude;
  FptrTableSS["AP1_amp"] = &SpikeShape::AP1_amp;
  FptrTableSS["AP2_amp"] = &SpikeShape::AP2_amp;
  FptrTableSS["mean_AP_amplitude"] = &SpikeShape::mean_AP_amplitude;
  FptrTableSS["APlast_amp"] = &SpikeShape::APlast_amp;
  FptrTableSS["AP_amplitude_change"] = &SpikeShape::AP_amplitude_change;
  FptrTableSS["AP_amplitude_from_voltagebase"] = &SpikeShape::AP_amplitude_from_voltagebase;
  FptrTableSS["AP1_peak"] = &SpikeShape::AP1_peak;
  FptrTableSS["AP2_peak"] = &SpikeShape::AP2_peak;
  FptrTableSS["AP2_AP1_diff"] = &SpikeShape::AP2_AP1_diff;
  FptrTableSS["AP2_AP1_peak_diff"] = &SpikeShape::AP2_AP1_peak_diff;
  FptrTableSS["amp_drop_first_second"] = &SpikeShape::amp_drop_first_second;
  FptrTableSS["amp_drop_first_last"] = &SpikeShape::amp_drop_first_last;
  FptrTableSS["amp_drop_second_last"] = &SpikeShape::amp_drop_second_last;
  FptrTableSS["max_amp_difference"] = &SpikeShape::max_amp_difference;
  FptrTableSS["AP_amplitude_diff"] = &SpikeShape::AP_amplitude_diff;
  FptrTableSS["min_AHP_indices"] = &SpikeShape::min_AHP_indices;
  FptrTableSS["min_AHP_values"] = &SpikeShape::min_AHP_values;
  FptrTableSS["AHP_depth_abs"] = &SpikeShape::AHP_depth_abs;
  FptrTableSS["AHP_depth_abs_slow"] = &SpikeShape::AHP_depth_abs_slow;
  FptrTableSS["AHP_slow_time"] = &SpikeShape::AHP_slow_time;
  FptrTableSS["AHP_depth_slow"] = &SpikeShape::AHP_depth_slow;
  FptrTableSS["AHP_depth"] = &SpikeShape::AHP_depth;
  FptrTableSS["AHP_depth_diff"] = &SpikeShape::AHP_depth_diff;
  FptrTableSS["fast_AHP"] = &SpikeShape::fast_AHP;
  FptrTableSS["fast_AHP_change"] = &SpikeShape::fast_AHP_change;
  FptrTableSS["AHP_depth_from_peak"] = &SpikeShape::AHP_depth_from_peak;
  FptrTableSS["AHP1_depth_from_peak"] = &SpikeShape::AHP1_depth_from_peak;
  FptrTableSS["AHP2_depth_from_peak"] = &SpikeShape::AHP2_depth_from_peak;
  FptrTableSS["AHP_time_from_peak"] = &SpikeShape::AHP_time_from_peak;
  FptrTableSS["ADP_peak_indices"] = &SpikeShape::ADP_peak_indices;
  FptrTableSS["ADP_peak_values"] = &SpikeShape::ADP_peak_values;
  FptrTableSS["ADP_peak_amplitude"] = &SpikeShape::ADP_peak_amplitude;
  FptrTableSS["depolarized_base"] = &SpikeShape::depolarized_base;
  FptrTableSS["min_voltage_between_spikes"] = &SpikeShape::min_voltage_between_spikes;
  FptrTableSS["min_between_peaks_indices"] = &SpikeShape::min_between_peaks_indices;
  FptrTableSS["min_between_peaks_values"] = &SpikeShape::min_between_peaks_values;
  FptrTableSS["AP_duration_half_width"] = &SpikeShape::AP_duration_half_width;
  FptrTableSS["AP_duration_half_width_change"] = &SpikeShape::AP_duration_half_width_change;
  FptrTableSS["AP_width"] = &SpikeShape::AP_width;
  FptrTableSS["AP_duration"] = &SpikeShape::AP_duration;
  FptrTableSS["AP_duration_change"] = &SpikeShape::AP_duration_change;
  FptrTableSS["AP_width_between_threshold"] = &SpikeShape::AP_width_between_threshold;
  FptrTableSS["spike_half_width"] = &SpikeShape::spike_width1;
  FptrTableSS["AP1_width"] = &SpikeShape::AP1_width;
  FptrTableSS["AP2_width"] = &SpikeShape::AP2_width;
  FptrTableSS["APlast_width"] = &SpikeShape::APlast_width;
  FptrTableSS["spike_width2"] = &SpikeShape::spike_width2;
  FptrTableSS["AP_begin_width"] = &SpikeShape::AP_begin_width;
  FptrTableSS["AP1_begin_width"] = &SpikeShape::AP1_begin_width;
  FptrTableSS["AP2_begin_width"] = &SpikeShape::AP2_begin_width;
  FptrTableSS["AP2_AP1_begin_width_diff"] = &SpikeShape::AP2_AP1_begin_width_diff;
  FptrTableSS["AP_begin_indices"] = &SpikeShape::AP_begin_indices;
  FptrTableSS["AP_end_indices"] = &SpikeShape::AP_end_indices;
  FptrTableSS["AP_begin_voltage"] = &SpikeShape::AP_begin_voltage;
  FptrTableSS["AP1_begin_voltage"] = &SpikeShape::AP1_begin_voltage;
  FptrTableSS["AP2_begin_voltage"] = &SpikeShape::AP2_begin_voltage;
  FptrTableSS["AP_begin_time"] = &SpikeShape::AP_begin_time;
  FptrTableSS["AP_peak_upstroke"] = &SpikeShape::AP_peak_upstroke;
  FptrTableSS["AP_peak_downstroke"] = &SpikeShape::AP_peak_downstroke;
  FptrTableSS["AP_rise_indices"] = &SpikeShape::AP_rise_indices;
  FptrTableSS["AP_fall_indices"] = &SpikeShape::AP_fall_indices;
  FptrTableSS["AP_rise_time"] = &SpikeShape::AP_rise_time;
  FptrTableSS["AP_fall_time"] = &SpikeShape::AP_fall_time;
  FptrTableSS["AP_rise_rate"] = &SpikeShape::AP_rise_rate;
  FptrTableSS["AP_fall_rate"] = &SpikeShape::AP_fall_rate;
  FptrTableSS["AP_rise_rate_change"] = &SpikeShape::AP_rise_rate_change;
  FptrTableSS["AP_fall_rate_change"] = &SpikeShape::AP_fall_rate_change;
  FptrTableSS["AP_phaseslope"] = &SpikeShape::AP_phaseslope;

  //******************  FptrTableST *****************************
  FptrTableST["steady_state_voltage_stimend"] = &Subthreshold::steady_state_voltage_stimend;
  FptrTableST["steady_state_hyper"] = &Subthreshold::steady_state_hyper;
  FptrTableST["steady_state_voltage"] = &Subthreshold::steady_state_voltage;
  FptrTableST["voltage_base"] = &Subthreshold::voltage_base;
  FptrTableST["current_base"] = &Subthreshold::current_base;
  FptrTableST["time_constant"] = &Subthreshold::time_constant;
  FptrTableST["decay_time_constant_after_stim"] = &Subthreshold::decay_time_constant_after_stim;
  FptrTableST["multiple_decay_time_constant_after_stim"] = &Subthreshold::multiple_decay_time_constant_after_stim;
  FptrTableST["sag_time_constant"] = &Subthreshold::sag_time_constant;
  FptrTableST["sag_amplitude"] = &Subthreshold::sag_amplitude;
  FptrTableST["sag_ratio1"] = &Subthreshold::sag_ratio1;
  FptrTableST["sag_ratio2"] = &Subthreshold::sag_ratio2;
  FptrTableST["ohmic_input_resistance"] = &Subthreshold::ohmic_input_resistance;
  FptrTableST["ohmic_input_resistance_vb_ssse"] = &Subthreshold::ohmic_input_resistance_vb_ssse;
  FptrTableST["voltage_deflection_vb_ssse"] = &Subthreshold::voltage_deflection_vb_ssse;
  FptrTableST["voltage_deflection"] = &Subthreshold::voltage_deflection;
  FptrTableST["voltage_deflection_begin"] = &Subthreshold::voltage_deflection_begin;
  FptrTableST["voltage_after_stim"] = &Subthreshold::voltage_after_stim;
  FptrTableST["maximum_voltage"] = &Subthreshold::maximum_voltage;
  FptrTableST["minimum_voltage"] = &Subthreshold::minimum_voltage;
  FptrTableST["maximum_voltage_from_voltagebase"] = &Subthreshold::maximum_voltage_from_voltagebase;
  
  return 1;
}
