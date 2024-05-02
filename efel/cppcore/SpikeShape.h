/* Copyright (c) 2015-2024, EPFL/Blue Brain Project                                   
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

#include "mapoperations.h"
#include "Utils.h"

#include <vector>
#include <stdexcept>

using std::vector;

namespace SpikeShape {
int peak_voltage(mapStr2intVec& IntFeatureData,
                 mapStr2doubleVec& DoubleFeatureData,
                 mapStr2Str& StringData);
int AP_height(mapStr2intVec& IntFeatureData,
              mapStr2doubleVec& DoubleFeatureData,
              mapStr2Str& StringData);
int AP_amplitude(mapStr2intVec& IntFeatureData,
                 mapStr2doubleVec& DoubleFeatureData,
                 mapStr2Str& StringData);
int AP1_amp(mapStr2intVec& IntFeatureData,
            mapStr2doubleVec& DoubleFeatureData,
            mapStr2Str& StringData);
int AP2_amp(mapStr2intVec& IntFeatureData,
            mapStr2doubleVec& DoubleFeatureData,
            mapStr2Str& StringData);
int mean_AP_amplitude(mapStr2intVec& IntFeatureData,
                      mapStr2doubleVec& DoubleFeatureData,
                      mapStr2Str& StringData);
int APlast_amp(mapStr2intVec& IntFeatureData,
               mapStr2doubleVec& DoubleFeatureData,
               mapStr2Str& StringData);
int AP_amplitude_change(mapStr2intVec& IntFeatureData,
                        mapStr2doubleVec& DoubleFeatureData,
                        mapStr2Str& StringData);
int AP_amplitude_from_voltagebase(mapStr2intVec& IntFeatureData,
                                  mapStr2doubleVec& DoubleFeatureData,
                                  mapStr2Str& StringData);
int AP1_peak(mapStr2intVec& IntFeatureData,
             mapStr2doubleVec& DoubleFeatureData,
             mapStr2Str& StringData);
int AP2_peak(mapStr2intVec& IntFeatureData,
             mapStr2doubleVec& DoubleFeatureData,
             mapStr2Str& StringData);
int AP2_AP1_diff(mapStr2intVec& IntFeatureData,
                 mapStr2doubleVec& DoubleFeatureData,
                 mapStr2Str& StringData);
int AP2_AP1_peak_diff(mapStr2intVec& IntFeatureData,
                      mapStr2doubleVec& DoubleFeatureData,
                      mapStr2Str& StringData);
int amp_drop_first_second(mapStr2intVec& IntFeatureData,
                          mapStr2doubleVec& DoubleFeatureData,
                          mapStr2Str& StringData);
int amp_drop_first_last(mapStr2intVec& IntFeatureData,
                        mapStr2doubleVec& DoubleFeatureData,
                        mapStr2Str& StringData);
int amp_drop_second_last(mapStr2intVec& IntFeatureData,
                         mapStr2doubleVec& DoubleFeatureData,
                         mapStr2Str& StringData);
int max_amp_difference(mapStr2intVec& IntFeatureData,
                       mapStr2doubleVec& DoubleFeatureData,
                       mapStr2Str& StringData);
int AP_amplitude_diff(mapStr2intVec& IntFeatureData,
                      mapStr2doubleVec& DoubleFeatureData,
                      mapStr2Str& StringData);
int min_AHP_indices(mapStr2intVec& IntFeatureData,
                    mapStr2doubleVec& DoubleFeatureData,
                    mapStr2Str& StringData);
int min_AHP_values(mapStr2intVec& IntFeatureData,
                   mapStr2doubleVec& DoubleFeatureData,
                   mapStr2Str& StringData);
int AHP_depth_abs(mapStr2intVec& IntFeatureData,
                  mapStr2doubleVec& DoubleFeatureData,
                  mapStr2Str& StringData);
int AHP_depth_abs_slow(mapStr2intVec& IntFeatureData,
                       mapStr2doubleVec& DoubleFeatureData,
                       mapStr2Str& StringData);
int AHP_slow_time(mapStr2intVec& IntFeatureData,
                  mapStr2doubleVec& DoubleFeatureData,
                  mapStr2Str& StringData);
int AHP_depth_slow(mapStr2intVec& IntFeatureData,
                   mapStr2doubleVec& DoubleFeatureData,
                   mapStr2Str& StringData);
int AHP_depth(mapStr2intVec& IntFeatureData,
              mapStr2doubleVec& DoubleFeatureData,
              mapStr2Str& StringData);
int AHP_depth_diff(mapStr2intVec& IntFeatureData,
                   mapStr2doubleVec& DoubleFeatureData,
                   mapStr2Str& StringData);
int fast_AHP(mapStr2intVec& IntFeatureData,
             mapStr2doubleVec& DoubleFeatureData,
             mapStr2Str& StringData);
int fast_AHP_change(mapStr2intVec& IntFeatureData,
                    mapStr2doubleVec& DoubleFeatureData,
                    mapStr2Str& StringData);
int AHP_depth_from_peak(mapStr2intVec& IntFeatureData,
                        mapStr2doubleVec& DoubleFeatureData,
                        mapStr2Str& StringData);
int AHP1_depth_from_peak(mapStr2intVec& IntFeatureData,
                         mapStr2doubleVec& DoubleFeatureData,
                         mapStr2Str& StringData);
int AHP2_depth_from_peak(mapStr2intVec& IntFeatureData,
                         mapStr2doubleVec& DoubleFeatureData,
                         mapStr2Str& StringData);
int AHP_time_from_peak(mapStr2intVec& IntFeatureData,
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
int depolarized_base(mapStr2intVec& IntFeatureData,
                     mapStr2doubleVec& DoubleFeatureData,
                     mapStr2Str& StringData);
int min_voltage_between_spikes(mapStr2intVec& IntFeatureData,
                               mapStr2doubleVec& DoubleFeatureData,
                               mapStr2Str& StringData);
int min_between_peaks_indices(mapStr2intVec& IntFeatureData,
                              mapStr2doubleVec& DoubleFeatureData,
                              mapStr2Str& StringData);
int min_between_peaks_values(mapStr2intVec& IntFeatureData,
                             mapStr2doubleVec& DoubleFeatureData,
                             mapStr2Str& StringData);
int AP_duration_half_width(mapStr2intVec& IntFeatureData,
                           mapStr2doubleVec& DoubleFeatureData,
                           mapStr2Str& StringData);
int AP_duration_half_width_change(mapStr2intVec& IntFeatureData,
                                  mapStr2doubleVec& DoubleFeatureData,
                                  mapStr2Str& StringData);
int AP_width(mapStr2intVec& IntFeatureData,
             mapStr2doubleVec& DoubleFeatureData,
             mapStr2Str& StringData);
int AP_duration(mapStr2intVec& IntFeatureData,
                mapStr2doubleVec& DoubleFeatureData,
                mapStr2Str& StringData);
int AP_duration_change(mapStr2intVec& IntFeatureData,
                       mapStr2doubleVec& DoubleFeatureData,
                       mapStr2Str& StringData);
int AP_width_between_threshold(mapStr2intVec& IntFeatureData,
                               mapStr2doubleVec& DoubleFeatureData,
                               mapStr2Str& StringData);
int spike_width1(mapStr2intVec& IntFeatureData,
                 mapStr2doubleVec& DoubleFeatureData,
                 mapStr2Str& StringData);
int AP1_width(mapStr2intVec& IntFeatureData,
              mapStr2doubleVec& DoubleFeatureData,
              mapStr2Str& StringData);
int AP2_width(mapStr2intVec& IntFeatureData,
              mapStr2doubleVec& DoubleFeatureData,
              mapStr2Str& StringData);
int APlast_width(mapStr2intVec& IntFeatureData,
                 mapStr2doubleVec& DoubleFeatureData,
                 mapStr2Str& StringData);
int spike_width2(mapStr2intVec& IntFeatureData,
                 mapStr2doubleVec& DoubleFeatureData,
                 mapStr2Str& StringData);
int AP_begin_width(mapStr2intVec& IntFeatureData,
                   mapStr2doubleVec& DoubleFeatureData,
                   mapStr2Str& StringData);
int AP1_begin_width(mapStr2intVec& IntFeatureData,
                    mapStr2doubleVec& DoubleFeatureData,
                    mapStr2Str& StringData);
int AP2_begin_width(mapStr2intVec& IntFeatureData,
                    mapStr2doubleVec& DoubleFeatureData,
                    mapStr2Str& StringData);
int AP2_AP1_begin_width_diff(mapStr2intVec& IntFeatureData,
                             mapStr2doubleVec& DoubleFeatureData,
                             mapStr2Str& StringData);
int AP_begin_indices(mapStr2intVec& IntFeatureData,
                     mapStr2doubleVec& DoubleFeatureData,
                     mapStr2Str& StringData);
int AP_end_indices(mapStr2intVec& IntFeatureData,
                   mapStr2doubleVec& DoubleFeatureData,
                   mapStr2Str& StringData);
int AP_begin_voltage(mapStr2intVec& IntFeatureData,
                     mapStr2doubleVec& DoubleFeatureData,
                     mapStr2Str& StringData);
int AP1_begin_voltage(mapStr2intVec& IntFeatureData,
                      mapStr2doubleVec& DoubleFeatureData,
                      mapStr2Str& StringData);
int AP2_begin_voltage(mapStr2intVec& IntFeatureData,
                      mapStr2doubleVec& DoubleFeatureData,
                      mapStr2Str& StringData);
int AP_begin_time(mapStr2intVec& IntFeatureData,
                  mapStr2doubleVec& DoubleFeatureData,
                  mapStr2Str& StringData);
int AP_peak_upstroke(mapStr2intVec& IntFeatureData,
                     mapStr2doubleVec& DoubleFeatureData,
                     mapStr2Str& StringData);
int AP_peak_downstroke(mapStr2intVec& IntFeatureData,
                       mapStr2doubleVec& DoubleFeatureData,
                       mapStr2Str& StringData);
int AP_rise_indices(mapStr2intVec& IntFeatureData,
                    mapStr2doubleVec& DoubleFeatureData,
                    mapStr2Str& StringData);
int AP_fall_indices(mapStr2intVec& IntFeatureData,
                    mapStr2doubleVec& DoubleFeatureData,
                    mapStr2Str& StringData);
int AP_rise_time(mapStr2intVec& IntFeatureData,
                 mapStr2doubleVec& DoubleFeatureData,
                 mapStr2Str& StringData);
int AP_fall_time(mapStr2intVec& IntFeatureData,
                 mapStr2doubleVec& DoubleFeatureData,
                 mapStr2Str& StringData);
int AP_rise_rate(mapStr2intVec& IntFeatureData,
                 mapStr2doubleVec& DoubleFeatureData,
                 mapStr2Str& StringData);
int AP_fall_rate(mapStr2intVec& IntFeatureData,
                 mapStr2doubleVec& DoubleFeatureData,
                 mapStr2Str& StringData);
int AP_rise_rate_change(mapStr2intVec& IntFeatureData,
                        mapStr2doubleVec& DoubleFeatureData,
                        mapStr2Str& StringData);
int AP_fall_rate_change(mapStr2intVec& IntFeatureData,
                        mapStr2doubleVec& DoubleFeatureData,
                        mapStr2Str& StringData);
int AP_phaseslope(mapStr2intVec& IntFeatureData,
                  mapStr2doubleVec& DoubleFeatureData,
                  mapStr2Str& StringData);
}