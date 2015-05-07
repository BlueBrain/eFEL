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
#include <iterator>
#include <algorithm>
#include <map>
#include <string>
#include <vector>
#include <list>
#include "mapoperations.h"
#include "Utils.h"

using namespace std;

typedef map<string, vector<int> > mapStr2intVec;
typedef map<string, vector<double> > mapStr2doubleVec;
typedef map<string, string> mapStr2Str;

namespace LibV3 {
int interpolate(mapStr2intVec& IntFeatureData,
                mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData);

int trace_check(mapStr2intVec& IntFeatureData,
                mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData);

int __peak_indices(double dThreshold, vector<double>& V,
                   vector<int>& PeakIndex);
int peak_indices(mapStr2intVec& IntFeatureData,
                 mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData);

int ISI_values(mapStr2intVec& IntFeatureData,
               mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData);

int __ISI_CV(const vector<double>& isivalues, vector<double>& isicv);
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

int __spike_width1(const vector<double>& t, const vector<double>& V,
                   const vector<int>& PeakIndex, const vector<int>& minAHPIndex,
                   double stim_start, vector<double>& spike_width2);
int spike_width1(mapStr2intVec& IntFeatureData,
                 mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData);

int min_AHP_indices(mapStr2intVec& IntFeatureData,
                    mapStr2doubleVec& DoubleFeatureData,
                    mapStr2Str& StringData);
int min_AHP_values(mapStr2intVec& IntFeatureData,
                   mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData);

int AHP_depth_abs(mapStr2intVec& IntFeatureData,
                  mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData);

int __rest_voltage_value(vector<double>& V, vector<int>& PeakIndex,
                         vector<double>& PeakVoltage);
int rest_voltage_value(mapStr2intVec& IntFeatureData,
                       mapStr2doubleVec& DoubleFeatureData,
                       mapStr2Str& StringData);

int __adaptation_index2(double stimStart, double stimEnd, double Offset,
                        const vector<double>& peakvoltagetime,
                        vector<double>& adaptationindex2);
int adaptation_index2(mapStr2intVec& IntFeatureData,
                      mapStr2doubleVec& DoubleFeatureData,
                      mapStr2Str& StringData);

int AP_height(mapStr2intVec& IntFeatureData,
              mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData);
int AP_amplitude(mapStr2intVec& IntFeatureData,
                 mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData);

int __AP_width(const vector<double>& t, const vector<double>& v,
               double stimstart, double threshold,
               const vector<int>& peakindices, const vector<int>& minahpindices,
               vector<double>& apwidth);
int AP_width(mapStr2intVec& IntFeatureData, mapStr2doubleVec& DoubleFeatureData,
             mapStr2Str& StringData);

int doublet_ISI(mapStr2intVec& IntFeatureData,
                mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData);

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

int __depolarized_base(const vector<double>& t, const vector<double>& v,
                       double stimstart, double stimend,
                       const vector<int>& apbi, const vector<int>& apendi,
                       vector<double>& dep_base);

int depolarized_base(mapStr2intVec& IntFeatureData,
                     mapStr2doubleVec& DoubleFeatureData,
                     mapStr2Str& StringData);
}
#endif
