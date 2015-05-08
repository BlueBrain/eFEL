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
#include <iterator>
#include <vector>
#include <algorithm>
#include <functional>
#include <math.h>
#include <map>
#include <string>
#include <fstream>

#include "Utils.h"
#include "mapoperations.h"
using namespace std;
typedef map<string, vector<int> > mapStr2intVec;
typedef map<string, vector<double> > mapStr2doubleVec;
typedef map<string, string> mapStr2Str;
namespace LibV2 {
// AP parameters of eCode Specification 1.04
// partly reimplemented Shaul's matlab code ap_points.m
//
int __AP_begin_indices(const vector<double>& t, const vector<double>& v,
                       double stimstart, double stimend,
                       const vector<int>& ahpi, vector<int>& apbi);
int AP_begin_indices(mapStr2intVec& IntFeatureData,
                     mapStr2doubleVec& DoubleFeatureData,
                     mapStr2Str& StringData);

int __AP_rise_indices(const vector<double>& v, const vector<int>& apbi,
                      const vector<int>& pi, vector<int>& apri);
int AP_rise_indices(mapStr2intVec& IntFeatureData,
                    mapStr2doubleVec& DoubleFeatureData,
                    mapStr2Str& StringData);

int __AP_fall_indices(const vector<double>& v, const vector<int>& apbi,
                      const vector<int>& apei, const vector<int>& pi,
                      vector<int>& apfi);
int AP_fall_indices(mapStr2intVec& IntFeatureData,
                    mapStr2doubleVec& DoubleFeatureData,
                    mapStr2Str& StringData);

int __AP_end_indices(const vector<double>& t, const vector<double>& v,
                     const vector<int>& pi, vector<int>& apei);
int AP_end_indices(mapStr2intVec& IntFeatureData,
                   mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData);

// eFeatures
int __AP_duration(const vector<double>& t, const vector<int>& apbeginindices,
                  const vector<int>& apendindices, vector<double>& apduration);
int AP_duration(mapStr2intVec& IntFeatureData,
                mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData);
int __AP_rise_time(const vector<double>& t, const vector<int>& apbeginindices,
                   const vector<int>& peakindices, vector<double>& aprisetime);
int AP_rise_time(mapStr2intVec& IntFeatureData,
                 mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData);
int __AP_fall_time(const vector<double>& t, const vector<int>& peakindices,
                   const vector<int>& apendindices, vector<double>& apfalltime);
int AP_fall_time(mapStr2intVec& IntFeatureData,
                 mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData);
int __AP_rise_rate(const vector<double>& t, const vector<double>& v,
                   const vector<int>& apbeginindices,
                   const vector<int>& peakindices, vector<double>& apriserate);
int AP_rise_rate(mapStr2intVec& IntFeatureData,
                 mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData);
int __AP_fall_rate(const vector<double>& t, const vector<double>& v,
                   const vector<int>& peakindices,
                   const vector<int>& apendindices, vector<double>& apfallrate);
int AP_fall_rate(mapStr2intVec& IntFeatureData,
                 mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData);
int __fast_AHP(const vector<double>& v, const vector<int>& apbeginindices,
               const vector<int>& minahpindices, vector<double>& fastahp);
int fast_AHP(mapStr2intVec& IntFeatureData, mapStr2doubleVec& DoubleFeatureData,
             mapStr2Str& StringData);
int __AP_amplitude_change(const vector<double>& apamplitude,
                          vector<double>& apamplitudechange);
int AP_amplitude_change(mapStr2intVec& IntFeatureData,
                        mapStr2doubleVec& DoubleFeatureData,
                        mapStr2Str& StringData);
int __AP_duration_change(const vector<double>& apduration,
                         vector<double>& apdurationchange);
int AP_duration_change(mapStr2intVec& IntFeatureData,
                       mapStr2doubleVec& DoubleFeatureData,
                       mapStr2Str& StringData);
int __AP_rise_rate_change(const vector<double>& apriserate,
                          vector<double>& apriseratechange);
int AP_rise_rate_change(mapStr2intVec& IntFeatureData,
                        mapStr2doubleVec& DoubleFeatureData,
                        mapStr2Str& StringData);
int __AP_fall_rate_change(const vector<double>& apfallrate,
                          vector<double>& apfallratechange);
int AP_fall_rate_change(mapStr2intVec& IntFeatureData,
                        mapStr2doubleVec& DoubleFeatureData,
                        mapStr2Str& StringData);
int __fast_AHP_change(const vector<double>& fastahp,
                      vector<double>& fastahpchange);
int fast_AHP_change(mapStr2intVec& IntFeatureData,
                    mapStr2doubleVec& DoubleFeatureData,
                    mapStr2Str& StringData);
int __AP_duration_half_width(const vector<double>& t,
                             const vector<int>& apriseindices,
                             const vector<int>& apfallindices,
                             vector<double>& apdurationhalfwidth);
int AP_duration_half_width(mapStr2intVec& IntFeatureData,
                           mapStr2doubleVec& DoubleFeatureData,
                           mapStr2Str& StringData);
int __AP_duration_half_width_change(const vector<double>& apdurationhalfwidth,
                                    vector<double>& apdurationhalfwidthchange);
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
int __steady_state_hyper(const vector<double>& v, const vector<double>& t,
                         double stimend, vector<double>& steady_state_hyper);
int steady_state_hyper(mapStr2intVec& IntFeatureData,
                       mapStr2doubleVec& DoubleFeatureData,
                       mapStr2Str& StringData);
int E40(mapStr2intVec& IntFeatureData, mapStr2doubleVec& DoubleFeatureData,
        mapStr2Str& StringData);
int __amp_drop_first_second(const vector<double>& peakvoltage,
                            vector<double>& ampdropfirstsecond);
int amp_drop_first_second(mapStr2intVec& IntFeatureData,
                          mapStr2doubleVec& DoubleFeatureData,
                          mapStr2Str& StringData);
int __amp_drop_first_last(const vector<double>& peakvoltage,
                          vector<double>& ampdropfirstlast);
int amp_drop_first_last(mapStr2intVec& IntFeatureData,
                        mapStr2doubleVec& DoubleFeatureData,
                        mapStr2Str& StringData);
int __amp_drop_second_last(const vector<double>& peakvoltage,
                           vector<double>& ampdropsecondlast);
int amp_drop_second_last(mapStr2intVec& IntFeatureData,
                         mapStr2doubleVec& DoubleFeatureData,
                         mapStr2Str& StringData);
int __max_amp_difference(const vector<double>& peakvoltage,
                         vector<double>& maxampdifference);
int max_amp_difference(mapStr2intVec& IntFeatureData,
                       mapStr2doubleVec& DoubleFeatureData,
                       mapStr2Str& StringData);
// end of feature definition
}
#endif
