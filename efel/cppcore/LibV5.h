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

using std::vector;

namespace LibV5 {
int ISI_log_slope(mapStr2intVec& IntFeatureData,
                  mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData);
int ISI_semilog_slope(mapStr2intVec& IntFeatureData,
                  mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData);
int ISI_log_slope_skip(mapStr2intVec& IntFeatureData,
                       mapStr2doubleVec& DoubleFeatureData,
                       mapStr2Str& StringData);

int time_to_second_spike(mapStr2intVec& IntFeatureData,
                         mapStr2doubleVec& DoubleFeatureData,
                         mapStr2Str& StringData);
int time_to_last_spike(mapStr2intVec& IntFeatureData,
                       mapStr2doubleVec& DoubleFeatureData,
                       mapStr2Str& StringData);
int inv_time_to_first_spike(mapStr2intVec& IntFeatureData,
                            mapStr2doubleVec& DoubleFeatureData,
                            mapStr2Str& StringData);
int inv_first_ISI(mapStr2intVec& IntFeatureData,
                  mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData);
int inv_second_ISI(mapStr2intVec& IntFeatureData,
                   mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData);
int inv_third_ISI(mapStr2intVec& IntFeatureData,
                  mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData);
int inv_fourth_ISI(mapStr2intVec& IntFeatureData,
                   mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData);
int inv_fifth_ISI(mapStr2intVec& IntFeatureData,
                  mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData);
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

int irregularity_index(mapStr2intVec& IntFeatureData,
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

int BPAPHeightLoc1(mapStr2intVec& IntFeatureData,
                   mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData);
int BPAPAmplitudeLoc1(mapStr2intVec& IntFeatureData,
                   mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData);
int BPAPAmplitudeLoc2(mapStr2intVec& IntFeatureData,
                   mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData);
int BPAPHeightLoc2(mapStr2intVec& IntFeatureData,
                   mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData);

int check_AISInitiation(mapStr2intVec&, mapStr2doubleVec&, mapStr2Str&);

int AP_phaseslope(mapStr2intVec&, mapStr2doubleVec&, mapStr2Str&);

int AP_phaseslope_AIS(mapStr2intVec&, mapStr2doubleVec&, mapStr2Str&);

int BAC_width(mapStr2intVec&, mapStr2doubleVec&, mapStr2Str&);
int BAC_maximum_voltage(mapStr2intVec&, mapStr2doubleVec&, mapStr2Str&);

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
int steady_state_voltage_stimend(mapStr2intVec& IntFeatureData,           
                                 mapStr2doubleVec& DoubleFeatureData,             
                                 mapStr2Str& StringData);
int voltage_base(mapStr2intVec& IntFeatureData,           
                 mapStr2doubleVec& DoubleFeatureData,             
                 mapStr2Str& StringData);
int decay_time_constant_after_stim(mapStr2intVec& IntFeatureData,
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
int Spikecount_stimint(mapStr2intVec& IntFeatureData,                                    
                       mapStr2doubleVec& DoubleFeatureData, 
                       mapStr2Str& StringData); 
int peak_indices(mapStr2intVec& IntFeatureData,                                    
                       mapStr2doubleVec& DoubleFeatureData, 
                       mapStr2Str& StringData); 
}
#endif
