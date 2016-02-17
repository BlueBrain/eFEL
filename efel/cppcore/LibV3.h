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

#ifndef __LIBV3
#define __LIBV3
#include "mapoperations.h"
#include "Utils.h"

#include <vector>

using std::vector;

namespace LibV3 {
int interpolate(mapStr2intVec& IntFeatureData,
                mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData);

int trace_check(mapStr2intVec& IntFeatureData,
                mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData);

int peak_indices(mapStr2intVec& IntFeatureData,
                 mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData);

int ISI_values(mapStr2intVec& IntFeatureData,
               mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData);

int ISI_CV(mapStr2intVec& IntFeatureData, mapStr2doubleVec& DoubleFeatureData,
           mapStr2Str& StringData);

int peak_voltage(mapStr2intVec& IntFeatureData,
                 mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData);

int firing_rate(mapStr2intVec& IntFeatureData,
                mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData);

int peak_time(mapStr2intVec& IntFeatureData,
              mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData);

int first_spike_time(mapStr2intVec& IntFeatureData,
                     mapStr2doubleVec& DoubleFeatureData,
                     mapStr2Str& StringData);

int spike_width1(mapStr2intVec& IntFeatureData,
                 mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData);

int min_AHP_indices(mapStr2intVec& IntFeatureData,
                    mapStr2doubleVec& DoubleFeatureData,
                    mapStr2Str& StringData);
int min_AHP_values(mapStr2intVec& IntFeatureData,
                   mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData);

int AHP_depth_abs(mapStr2intVec& IntFeatureData,
                  mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData);

int rest_voltage_value(mapStr2intVec& IntFeatureData,
                       mapStr2doubleVec& DoubleFeatureData,
                       mapStr2Str& StringData);

int adaptation_index2(mapStr2intVec& IntFeatureData,
                      mapStr2doubleVec& DoubleFeatureData,
                      mapStr2Str& StringData);

int AP_height(mapStr2intVec& IntFeatureData,
              mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData);
int AP_amplitude(mapStr2intVec& IntFeatureData,
                 mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData);

int AP_width(mapStr2intVec& IntFeatureData, mapStr2doubleVec& DoubleFeatureData,
             mapStr2Str& StringData);

int doublet_ISI(mapStr2intVec& IntFeatureData,
                mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData);

int AP_begin_indices(mapStr2intVec& IntFeatureData,
                     mapStr2doubleVec& DoubleFeatureData,
                     mapStr2Str& StringData);

int AP_rise_indices(mapStr2intVec& IntFeatureData,
                    mapStr2doubleVec& DoubleFeatureData,
                    mapStr2Str& StringData);

int AP_fall_indices(mapStr2intVec& IntFeatureData,
                    mapStr2doubleVec& DoubleFeatureData,
                    mapStr2Str& StringData);

int AP_end_indices(mapStr2intVec& IntFeatureData,
                   mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData);

// eFeatures
int AP_duration(mapStr2intVec& IntFeatureData,
                mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData);

int depolarized_base(mapStr2intVec& IntFeatureData,
                     mapStr2doubleVec& DoubleFeatureData,
                     mapStr2Str& StringData);
}
#endif
