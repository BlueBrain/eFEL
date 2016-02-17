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

#ifndef __LIBV2
#define __LIBV2

#include "Utils.h"
#include "mapoperations.h"

#include <vector>

using std::vector;

namespace LibV2 {
// AP parameters of eCode Specification 1.04
// partly reimplemented Shaul's matlab code ap_points.m
//
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
int AP_rise_time(mapStr2intVec& IntFeatureData,
                 mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData);
int AP_fall_time(mapStr2intVec& IntFeatureData,
                 mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData);
int AP_rise_rate(mapStr2intVec& IntFeatureData,
                 mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData);
int AP_fall_rate(mapStr2intVec& IntFeatureData,
                 mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData);
int fast_AHP(mapStr2intVec& IntFeatureData, mapStr2doubleVec& DoubleFeatureData,
             mapStr2Str& StringData);
int AP_amplitude_change(mapStr2intVec& IntFeatureData,
                        mapStr2doubleVec& DoubleFeatureData,
                        mapStr2Str& StringData);
int AP_duration_change(mapStr2intVec& IntFeatureData,
                       mapStr2doubleVec& DoubleFeatureData,
                       mapStr2Str& StringData);
int AP_rise_rate_change(mapStr2intVec& IntFeatureData,
                        mapStr2doubleVec& DoubleFeatureData,
                        mapStr2Str& StringData);
int AP_fall_rate_change(mapStr2intVec& IntFeatureData,
                        mapStr2doubleVec& DoubleFeatureData,
                        mapStr2Str& StringData);
int fast_AHP_change(mapStr2intVec& IntFeatureData,
                    mapStr2doubleVec& DoubleFeatureData,
                    mapStr2Str& StringData);
int AP_duration_half_width(mapStr2intVec& IntFeatureData,
                           mapStr2doubleVec& DoubleFeatureData,
                           mapStr2Str& StringData);
int AP_duration_half_width_change(mapStr2intVec& IntFeatureData,
                                  mapStr2doubleVec& DoubleFeatureData,
                                  mapStr2Str& StringData);

int E6(mapStr2intVec& IntFeatureData, mapStr2doubleVec& DoubleFeatureData,
       mapStr2Str& StringData);
int E7(mapStr2intVec& IntFeatureData, mapStr2doubleVec& DoubleFeatureData,
       mapStr2Str& StringData);

int BPAPatt2(mapStr2intVec& IntFeatureData, mapStr2doubleVec& DoubleFeatureData,
             mapStr2Str& StringData);
int BPAPatt3(mapStr2intVec& IntFeatureData, mapStr2doubleVec& DoubleFeatureData,
             mapStr2Str& StringData);
int E39(mapStr2intVec& IntFeatureData, mapStr2doubleVec& DoubleFeatureData,
        mapStr2Str& StringData);
int E39_cod(mapStr2intVec& IntFeatureData, mapStr2doubleVec& DoubleFeatureData,
            mapStr2Str& StringData);
int E2(mapStr2intVec& IntFeatureData, mapStr2doubleVec& DoubleFeatureData,
       mapStr2Str& StringData);
int E3(mapStr2intVec& IntFeatureData, mapStr2doubleVec& DoubleFeatureData,
       mapStr2Str& StringData);
int E4(mapStr2intVec& IntFeatureData, mapStr2doubleVec& DoubleFeatureData,
       mapStr2Str& StringData);
int E5(mapStr2intVec& IntFeatureData, mapStr2doubleVec& DoubleFeatureData,
       mapStr2Str& StringData);
int E8(mapStr2intVec& IntFeatureData, mapStr2doubleVec& DoubleFeatureData,
       mapStr2Str& StringData);
int E9(mapStr2intVec& IntFeatureData, mapStr2doubleVec& DoubleFeatureData,
       mapStr2Str& StringData);
int E10(mapStr2intVec& IntFeatureData, mapStr2doubleVec& DoubleFeatureData,
        mapStr2Str& StringData);
int E11(mapStr2intVec& IntFeatureData, mapStr2doubleVec& DoubleFeatureData,
        mapStr2Str& StringData);
int E12(mapStr2intVec& IntFeatureData, mapStr2doubleVec& DoubleFeatureData,
        mapStr2Str& StringData);
int E13(mapStr2intVec& IntFeatureData, mapStr2doubleVec& DoubleFeatureData,
        mapStr2Str& StringData);
int E14(mapStr2intVec& IntFeatureData, mapStr2doubleVec& DoubleFeatureData,
        mapStr2Str& StringData);
int E15(mapStr2intVec& IntFeatureData, mapStr2doubleVec& DoubleFeatureData,
        mapStr2Str& StringData);
int E16(mapStr2intVec& IntFeatureData, mapStr2doubleVec& DoubleFeatureData,
        mapStr2Str& StringData);
int E17(mapStr2intVec& IntFeatureData, mapStr2doubleVec& DoubleFeatureData,
        mapStr2Str& StringData);
int E18(mapStr2intVec& IntFeatureData, mapStr2doubleVec& DoubleFeatureData,
        mapStr2Str& StringData);
int E19(mapStr2intVec& IntFeatureData, mapStr2doubleVec& DoubleFeatureData,
        mapStr2Str& StringData);
int E20(mapStr2intVec& IntFeatureData, mapStr2doubleVec& DoubleFeatureData,
        mapStr2Str& StringData);
int E21(mapStr2intVec& IntFeatureData, mapStr2doubleVec& DoubleFeatureData,
        mapStr2Str& StringData);
int E22(mapStr2intVec& IntFeatureData, mapStr2doubleVec& DoubleFeatureData,
        mapStr2Str& StringData);
int E23(mapStr2intVec& IntFeatureData, mapStr2doubleVec& DoubleFeatureData,
        mapStr2Str& StringData);
int E24(mapStr2intVec& IntFeatureData, mapStr2doubleVec& DoubleFeatureData,
        mapStr2Str& StringData);
int E25(mapStr2intVec& IntFeatureData, mapStr2doubleVec& DoubleFeatureData,
        mapStr2Str& StringData);
int E26(mapStr2intVec& IntFeatureData, mapStr2doubleVec& DoubleFeatureData,
        mapStr2Str& StringData);
int E27(mapStr2intVec& IntFeatureData, mapStr2doubleVec& DoubleFeatureData,
        mapStr2Str& StringData);
int steady_state_hyper(mapStr2intVec& IntFeatureData,
                       mapStr2doubleVec& DoubleFeatureData,
                       mapStr2Str& StringData);
int E40(mapStr2intVec& IntFeatureData, mapStr2doubleVec& DoubleFeatureData,
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
// end of feature definition
}
#endif
