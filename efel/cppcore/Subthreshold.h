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

namespace Subthreshold {
int steady_state_voltage_stimend(mapStr2intVec& IntFeatureData,
                                 mapStr2doubleVec& DoubleFeatureData,
                                 mapStr2Str& StringData);
int steady_state_hyper(mapStr2intVec& IntFeatureData,
                       mapStr2doubleVec& DoubleFeatureData,
                       mapStr2Str& StringData);
int steady_state_voltage(mapStr2intVec& IntFeatureData,
                         mapStr2doubleVec& DoubleFeatureData,
                         mapStr2Str& StringData);
int voltage_base(mapStr2intVec& IntFeatureData,
                 mapStr2doubleVec& DoubleFeatureData,
                 mapStr2Str& StringData);
int current_base(mapStr2intVec& IntFeatureData,
                 mapStr2doubleVec& DoubleFeatureData,
                 mapStr2Str& StringData);
int time_constant(mapStr2intVec& IntFeatureData,
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
int sag_amplitude(mapStr2intVec& IntFeatureData,
                  mapStr2doubleVec& DoubleFeatureData,
                  mapStr2Str& StringData);
int sag_ratio1(mapStr2intVec& IntFeatureData,
               mapStr2doubleVec& DoubleFeatureData,
               mapStr2Str& StringData);
int sag_ratio2(mapStr2intVec& IntFeatureData,
               mapStr2doubleVec& DoubleFeatureData,
               mapStr2Str& StringData);
int ohmic_input_resistance(mapStr2intVec& IntFeatureData,
                           mapStr2doubleVec& DoubleFeatureData,
                           mapStr2Str& StringData);
int ohmic_input_resistance_vb_ssse(mapStr2intVec& IntFeatureData,
                                   mapStr2doubleVec& DoubleFeatureData,
                                   mapStr2Str& StringData);
int voltage_deflection_vb_ssse(mapStr2intVec& IntFeatureData,
                               mapStr2doubleVec& DoubleFeatureData,
                               mapStr2Str& StringData);
int voltage_deflection(mapStr2intVec& IntFeatureData,
                       mapStr2doubleVec& DoubleFeatureData,
                       mapStr2Str& StringData);
int voltage_deflection_begin(mapStr2intVec& IntFeatureData,
                             mapStr2doubleVec& DoubleFeatureData,
                             mapStr2Str& StringData);
int voltage_after_stim(mapStr2intVec& IntFeatureData,
                       mapStr2doubleVec& DoubleFeatureData,
                       mapStr2Str& StringData);
int maximum_voltage(mapStr2intVec& IntFeatureData,
                    mapStr2doubleVec& DoubleFeatureData,
                    mapStr2Str& StringData);
int minimum_voltage(mapStr2intVec& IntFeatureData,
                    mapStr2doubleVec& DoubleFeatureData,
                    mapStr2Str& StringData);
int maximum_voltage_from_voltagebase(mapStr2intVec& IntFeatureData,
                                     mapStr2doubleVec& DoubleFeatureData,
                                     mapStr2Str& StringData);
}