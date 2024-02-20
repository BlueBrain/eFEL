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

#ifndef __LIBV5
#define __LIBV5

#include "mapoperations.h"
#include "Utils.h"

#include <vector>
#include <stdexcept>

using std::vector;

namespace LibV5 {
int time_to_second_spike(mapStr2intVec& IntFeatureData,
                         mapStr2doubleVec& DoubleFeatureData,
                         mapStr2Str& StringData);
int time_to_last_spike(mapStr2intVec& IntFeatureData,
                       mapStr2doubleVec& DoubleFeatureData,
                       mapStr2Str& StringData);
int inv_time_to_first_spike(mapStr2intVec& IntFeatureData,
                            mapStr2doubleVec& DoubleFeatureData,
                            mapStr2Str& StringData);
int inv_ISI_generic(mapStr2intVec& IntFeatureData,
                           mapStr2doubleVec& DoubleFeatureData,
                           mapStr2Str& StringData,
                           size_t index);
int inv_last_ISI(mapStr2intVec& IntFeatureData,
                 mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData);

int min_AHP_indices(mapStr2intVec& intfeaturedata,
                    mapStr2doubleVec& DoubleFeatureData,
                    mapStr2Str& StringData);
int min_AHP_values(mapStr2intVec& IntFeatureData,
                   mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData);

int AHP_depth_abs(mapStr2intVec& IntFeatureData,
                  mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData);

int spike_width1(mapStr2intVec& IntFeatureData,
                 mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData);

int AP_begin_indices(mapStr2intVec& IntFeatureData,
                     mapStr2doubleVec& DoubleFeatureData,
                     mapStr2Str& StringData);

int AP_end_indices(mapStr2intVec& IntFeatureData,
                   mapStr2doubleVec& DoubleFeatureData,
                   mapStr2Str& StringData);

int number_initial_spikes(mapStr2intVec& IntFeatureData,
                          mapStr2doubleVec& DoubleFeatureData,
                          mapStr2Str& StringData);

int AP1_amp(mapStr2intVec& IntFeatureData, mapStr2doubleVec& DoubleFeatureData,
            mapStr2Str& StringData);
int APlast_amp(mapStr2intVec& IntFeatureData, 
               mapStr2doubleVec& DoubleFeatureData,
               mapStr2Str& StringData);
int AP2_amp(mapStr2intVec& IntFeatureData, mapStr2doubleVec& DoubleFeatureData,
            mapStr2Str& StringData);

int AP1_peak(mapStr2intVec& IntFeatureData, mapStr2doubleVec& DoubleFeatureData,
             mapStr2Str& StringData);
int AP2_peak(mapStr2intVec& IntFeatureData, mapStr2doubleVec& DoubleFeatureData,
             mapStr2Str& StringData);

int AP2_AP1_diff(mapStr2intVec& IntFeatureData,
                 mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData);
int AP2_AP1_peak_diff(mapStr2intVec& IntFeatureData,
                      mapStr2doubleVec& DoubleFeatureData,
                      mapStr2Str& StringData);

int AP1_width(mapStr2intVec& IntFeatureData,
              mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData);
int AP2_width(mapStr2intVec& IntFeatureData,
              mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData);
int APlast_width(mapStr2intVec& IntFeatureData,
              mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData);

int AHP_depth_from_peak(mapStr2intVec& IntFeatureData,
                        mapStr2doubleVec& DoubleFeatureData,
                        mapStr2Str& StringData);

int AHP_time_from_peak(mapStr2intVec& IntFeatureData,
                       mapStr2doubleVec& DoubleFeatureData,
                       mapStr2Str& StringData);

int AHP1_depth_from_peak(mapStr2intVec& IntFeatureData,
                         mapStr2doubleVec& DoubleFeatureData,
                         mapStr2Str& StringData);
int AHP2_depth_from_peak(mapStr2intVec& IntFeatureData,
                         mapStr2doubleVec& DoubleFeatureData,
                         mapStr2Str& StringData);

int AP_begin_width(mapStr2intVec& IntFeatureData,
                   mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData);

int AP_begin_time(mapStr2intVec& IntFeatureData,
                  mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData);

int AP_begin_voltage(mapStr2intVec& IntFeatureData,
                     mapStr2doubleVec& DoubleFeatureData,
                     mapStr2Str& StringData);

int AP1_begin_voltage(mapStr2intVec& IntFeatureData,
                      mapStr2doubleVec& DoubleFeatureData,
                      mapStr2Str& StringData);
int AP2_begin_voltage(mapStr2intVec& IntFeatureData,
                      mapStr2doubleVec& DoubleFeatureData,
                      mapStr2Str& StringData);

int AP1_begin_width(mapStr2intVec& IntFeatureData,
                    mapStr2doubleVec& DoubleFeatureData,
                    mapStr2Str& StringData);
int AP2_begin_width(mapStr2intVec& IntFeatureData,
                    mapStr2doubleVec& DoubleFeatureData,
                    mapStr2Str& StringData);

int voltage_deflection_begin(mapStr2intVec& IntFeatureData,
                             mapStr2doubleVec& DoubleFeatureData,
                             mapStr2Str& StringData);

int mean_AP_amplitude(mapStr2intVec& intfeaturedata,
                      mapStr2doubleVec& doublefeaturedata,
                      mapStr2Str& StringData);

int is_not_stuck(mapStr2intVec& IntFeatureData,
                 mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData);

int voltage_after_stim(mapStr2intVec& IntFeatureData,
                       mapStr2doubleVec& DoubleFeatureData,
                       mapStr2Str& StringData);

int AP2_AP1_begin_width_diff(mapStr2intVec& IntFeatureData,
                             mapStr2doubleVec& DoubleFeatureData,
                             mapStr2Str& StringData);

int AP_phaseslope(mapStr2intVec&, mapStr2doubleVec&, mapStr2Str&);

int all_ISI_values(mapStr2intVec& IntFeatureData,
                   mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData);

int AP_amplitude_from_voltagebase(mapStr2intVec& intfeaturedata,
                                  mapStr2doubleVec& doublefeaturedata,
                                  mapStr2Str& StringData);
int min_voltage_between_spikes(mapStr2intVec& intfeaturedata,
                                  mapStr2doubleVec& doublefeaturedata,
                                  mapStr2Str& StringData);
int voltage(mapStr2intVec& intfeaturedata,
                                  mapStr2doubleVec& doublefeaturedata,
                                  mapStr2Str& StringData);
int current(mapStr2intVec& intfeaturedata,
                                  mapStr2doubleVec& doublefeaturedata,
                                  mapStr2Str& StringData);
int time(mapStr2intVec& intfeaturedata,
                                  mapStr2doubleVec& doublefeaturedata,
                                  mapStr2Str& StringData);
int steady_state_voltage_stimend(mapStr2intVec& IntFeatureData,           
                                 mapStr2doubleVec& DoubleFeatureData,             
                                 mapStr2Str& StringData);
int voltage_base(mapStr2intVec& IntFeatureData,           
                 mapStr2doubleVec& DoubleFeatureData,             
                 mapStr2Str& StringData);
int current_base(mapStr2intVec& IntFeatureData,           
                 mapStr2doubleVec& DoubleFeatureData,             
                 mapStr2Str& StringData);
int decay_time_constant_after_stim(mapStr2intVec& IntFeatureData,
                                   mapStr2doubleVec& DoubleFeatureData,
                                   mapStr2Str& StringData);
int multiple_decay_time_constant_after_stim(mapStr2intVec& IntFeatureData,
                                      mapStr2doubleVec& DoubleFeatureData,
                                      mapStr2Str& StringData);
int sag_time_constant(mapStr2intVec& IntFeatureData,
                            mapStr2doubleVec& DoubleFeatureData,
                            mapStr2Str& StringData);
int voltage_deflection_vb_ssse(mapStr2intVec& IntFeatureData,
                                   mapStr2doubleVec& DoubleFeatureData,
                                   mapStr2Str& StringData);
int ohmic_input_resistance_vb_ssse(mapStr2intVec& IntFeatureData,
                                   mapStr2doubleVec& DoubleFeatureData,
                                   mapStr2Str& StringData);
int maximum_voltage_from_voltagebase(mapStr2intVec& IntFeatureData,
                                   mapStr2doubleVec& DoubleFeatureData,
                                   mapStr2Str& StringData);
int peak_indices(mapStr2intVec& IntFeatureData,                                    
                       mapStr2doubleVec& DoubleFeatureData, 
                       mapStr2Str& StringData); 
int sag_amplitude(mapStr2intVec& IntFeatureData,                                    
                       mapStr2doubleVec& DoubleFeatureData, 
                       mapStr2Str& StringData); 
int sag_ratio1(mapStr2intVec& IntFeatureData,                                    
                       mapStr2doubleVec& DoubleFeatureData, 
                       mapStr2Str& StringData); 
int sag_ratio2(mapStr2intVec& IntFeatureData,                                    
                       mapStr2doubleVec& DoubleFeatureData, 
                       mapStr2Str& StringData);
int AP_peak_upstroke(mapStr2intVec& IntFeatureData,
                    mapStr2doubleVec& DoubleFeatureData,
                    mapStr2Str& StringData);
int AP_peak_downstroke(mapStr2intVec& IntFeatureData,
                            mapStr2doubleVec& DoubleFeatureData,
                            mapStr2Str& StringData);
int min_between_peaks_indices(mapStr2intVec& IntFeatureData,
                            mapStr2doubleVec& DoubleFeatureData,
                            mapStr2Str& StringData);
int min_between_peaks_values(mapStr2intVec& IntFeatureData,
                            mapStr2doubleVec& DoubleFeatureData,
                            mapStr2Str& StringData);
int AP_width_between_threshold(mapStr2intVec& IntFeatureData,
                            mapStr2doubleVec& DoubleFeatureData,
                            mapStr2Str& StringData);
int burst_begin_indices(mapStr2intVec& IntFeatureData,
                            mapStr2doubleVec& DoubleFeatureData,
                            mapStr2Str& StringData);
int burst_end_indices(mapStr2intVec& IntFeatureData,
                            mapStr2doubleVec& DoubleFeatureData,
                            mapStr2Str& StringData);
int strict_burst_mean_freq(mapStr2intVec& IntFeatureData,
                            mapStr2doubleVec& DoubleFeatureData,
                            mapStr2Str& StringData);
int strict_interburst_voltage(mapStr2intVec& IntFeatureData,
                            mapStr2doubleVec& DoubleFeatureData,
                            mapStr2Str& StringData);
int ADP_peak_indices(mapStr2intVec& IntFeatureData,
                            mapStr2doubleVec& DoubleFeatureData,
                            mapStr2Str& StringData);
int ADP_peak_values(mapStr2intVec& IntFeatureData,
                            mapStr2doubleVec& DoubleFeatureData,
                            mapStr2Str& StringData);
int ADP_peak_amplitude(mapStr2intVec& IntFeatureData,
                            mapStr2doubleVec& DoubleFeatureData,
                            mapStr2Str& StringData);
int interburst_min_indices(mapStr2intVec& IntFeatureData,
                            mapStr2doubleVec& DoubleFeatureData,
                            mapStr2Str& StringData);
int interburst_min_values(mapStr2intVec& IntFeatureData,
                            mapStr2doubleVec& DoubleFeatureData,
                            mapStr2Str& StringData);
int postburst_min_indices(mapStr2intVec& IntFeatureData,
                            mapStr2doubleVec& DoubleFeatureData,
                            mapStr2Str& StringData);
int postburst_min_values(mapStr2intVec& IntFeatureData,
                            mapStr2doubleVec& DoubleFeatureData,
                            mapStr2Str& StringData);
int time_to_interburst_min(mapStr2intVec& IntFeatureData,
                            mapStr2doubleVec& DoubleFeatureData,
                            mapStr2Str& StringData);
}
#endif
