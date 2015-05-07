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

#ifndef __LIBV1
#define __LIBV1
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
namespace LibV1 {
int printVectorI(char* strName, vector<int> vec);
int printVectorD(char* strName, vector<double> vec);
// int __peak_indices(double dThreshold, vector<double> &V, vector<int>
// &PeakIndex);
int interpolate(mapStr2intVec& IntFeatureData,
                mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData);

int __peak_indices(double dThreshold, vector<double>& V,
                   vector<int>& PeakIndex);
int peak_indices(mapStr2intVec& IntFeatureData,
                 mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData);

int __ISI_values(vector<double>& V, vector<int>& PeakIndex,
                 vector<double>& PeakVoltage);
int ISI_values(mapStr2intVec& IntFeatureData,
               mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData);

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

int __spike_width2(vector<double>& t, vector<double>& V, vector<int>& PeakIndex,
                   vector<int>& minAHPIndex, vector<double>& spike_width2);
int spike_width2(mapStr2intVec& IntFeatureData,
                 mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData);

int min_AHP_indices(mapStr2intVec& IntFeatureData,
                    mapStr2doubleVec& DoubleFeatureData,
                    mapStr2Str& StringData);

int AHP_depth_abs(mapStr2intVec& IntFeatureData,
                  mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData);

int min_AHP_values(mapStr2intVec& IntFeatureData,
                   mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData);

int __rest_voltage_value(vector<double>& V, vector<int>& PeakIndex,
                         vector<double>& PeakVoltage);
int rest_voltage_value(mapStr2intVec& IntFeatureData,
                       mapStr2doubleVec& DoubleFeatureData,
                       mapStr2Str& StringData);

int __burst_ISI_indices(double BurstFactor, vector<int>& PeakIndex,
                        vector<double>& ISIValues, vector<int>& BurstIndex);
int burst_ISI_indices(mapStr2intVec& IntFeatureData,
                      mapStr2doubleVec& DoubleFeatureData,
                      mapStr2Str& StringData);

int __burst_mean_freq(vector<double>& PVTime, vector<int>& BurstIndex,
                      vector<double>& BurstMeanFreq);
int burst_mean_freq(mapStr2intVec& IntFeatureData,
                    mapStr2doubleVec& DoubleFeatureData,
                    mapStr2Str& StringData);

int __interburst_voltage(vector<int>& BrustIndex, vector<int>& PeakIndex,
                         vector<double>& T, vector<double>& V,
                         vector<double>& IBV);
int interburst_voltage(mapStr2intVec& IntFeatureData,
                       mapStr2doubleVec& DoubleFeatureData,
                       mapStr2Str& StringData);

int __adaptation_index(double spikeSkipf, int maxnSpike, double StimStart,
                       double StimEnd, double Offset, vector<double>& peakVTime,
                       vector<double>& adaptation_index);
int adaptation_index(mapStr2intVec& IntFeatureData,
                     mapStr2doubleVec& DoubleFeatureData,
                     mapStr2Str& StringData);

int trace_check(mapStr2intVec& IntFeatureData,
                mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData);

// passive properties
int time_constant(mapStr2intVec& IntFeatureData,
                  mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData);
int __time_constant(const vector<double>& V, const vector<double>& t,
                    double stimStart, double stimEnd, vector<double>& tc);

int voltage_deflection(mapStr2intVec& IntFeatureData,
                       mapStr2doubleVec& DoubleFeatureData,
                       mapStr2Str& StringData);
int __voltage_deflection(const vector<double>& V, const vector<double>& t,
                         double stimStart, double stimEnd, vector<double>& vd);

int ohmic_input_resistance(mapStr2intVec& IntFeatureData,
                           mapStr2doubleVec& DoubleFeatureData,
                           mapStr2Str& StringData);
int __ohmic_input_resistance(double voltage_deflection, double stimulus_current,
                             vector<double>& oir);

int maximum_voltage(mapStr2intVec& IntFeatureData,
                    mapStr2doubleVec& DoubleFeatureData,
                    mapStr2Str& StringData);
int minimum_voltage(mapStr2intVec& IntFeatureData,
                    mapStr2doubleVec& DoubleFeatureData,
                    mapStr2Str& StringData);
int __maxmin_voltage(const vector<double>& V, const vector<double>& t,
                     double stimStart, double stimEnd, vector<double>& maxV,
                     vector<double>& minV);

int steady_state_voltage(mapStr2intVec& IntFeatureData,
                         mapStr2doubleVec& DoubleFeatureData,
                         mapStr2Str& StringData);
int __steady_state_voltage(const vector<double>& V, const vector<double>& t,
                           double stimEnd, vector<double>& ssv);

int AP_height(mapStr2intVec& IntFeatureData,
              mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData);
int AP_amplitude(mapStr2intVec& IntFeatureData,
                 mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData);
int AHP_depth(mapStr2intVec& IntFeatureData,
              mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData);
int __single_burst_ratio(const vector<double>& isivalues,
                         vector<double>& singleburstratio);
int single_burst_ratio(mapStr2intVec& IntFeatureData,
                       mapStr2doubleVec& DoubleFeatureData,
                       mapStr2Str& StringData);
int threshold_current(mapStr2intVec& IntFeatureData,
                      mapStr2doubleVec& DoubleFeatureData,
                      mapStr2Str& StringData);
int __AP_width(const vector<double>& t, const vector<double>& v,
               double stimstart, double threshold,
               const vector<int>& peakindices, const vector<int>& minahpindices,
               vector<double>& apwidth);
int AP_width(mapStr2intVec& IntFeatureData, mapStr2doubleVec& DoubleFeatureData,
             mapStr2Str& StringData);
int doublet_ISI(mapStr2intVec& IntFeatureData,
                mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData);
int __adaptation_index2(double stimStart, double stimEnd, double Offset,
                        const vector<double>& peakvoltagetime,
                        vector<double>& adaptationindex2);
int adaptation_index2(mapStr2intVec& IntFeatureData,
                      mapStr2doubleVec& DoubleFeatureData,
                      mapStr2Str& StringData);
int __ISI_CV(const vector<double>& isivalues, vector<double>& isicv);
int ISI_CV(mapStr2intVec& IntFeatureData, mapStr2doubleVec& DoubleFeatureData,
           mapStr2Str& StringData);
int __AHP_depth_abs_slow_indices(const vector<double>& t,
                                 const vector<double>& v,
                                 const vector<int>& peakindices,
                                 vector<int>& adas_indices);
int AHP_depth_abs_slow(mapStr2intVec& IntFeatureData,
                       mapStr2doubleVec& DoubleFeatureData,
                       mapStr2Str& StringData);
int AHP_slow_time(mapStr2intVec& IntFeatureData,
                  mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData);
int Spikecount(mapStr2intVec& IntFeatureData,
               mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData);
int __AHP_depth(const vector<double>& voltagebase,
                const vector<double>& minahpvalues, vector<double>& ahpdepth);
int AHP_depth(mapStr2intVec& IntFeatureData,
              mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData);
int burst_number(mapStr2intVec& IntFeatureData,
                 mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData);
int __AP_amplitude_diff(const vector<double>& apamplitude,
                        vector<double>& apamplitudediff);
int AP_amplitude_diff(mapStr2intVec& IntFeatureData,
                      mapStr2doubleVec& DoubleFeatureData,
                      mapStr2Str& StringData);
int __AHP_depth_diff(const vector<double>& ahpdepth,
                     vector<double>& ahpdepthdiff);
int AHP_depth_diff(mapStr2intVec& IntFeatureData,
                   mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData);

// end of feature definition
}
#endif
