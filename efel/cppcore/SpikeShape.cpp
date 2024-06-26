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


#include "SpikeShape.h"

#include <math.h>

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <deque>
#include <functional>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <list>
#include <sstream>
#include <string>

#include "EfelExceptions.h"

using std::distance;
using std::find_if;
using std::max_element;
using std::min_element;
using std::transform;


int SpikeShape::peak_voltage(mapStr2intVec& IntFeatureData,
                        mapStr2doubleVec& DoubleFeatureData,
                        mapStr2Str& StringData) {
  const auto& doubleFeatures = getFeatures(DoubleFeatureData, {"V"});
  const auto& intFeatures = getFeatures(IntFeatureData, {"peak_indices"});
  vector<double> peakV;
  for (const auto& index : intFeatures.at("peak_indices")) {
    peakV.push_back(doubleFeatures.at("V")[index]);
  }
  setVec(DoubleFeatureData, StringData, "peak_voltage", peakV);
  return peakV.size();
}

int SpikeShape::AP_height(mapStr2intVec& IntFeatureData,
                     mapStr2doubleVec& DoubleFeatureData,
                     mapStr2Str& StringData) {
  const auto& doubleFeatures = getFeatures(DoubleFeatureData, {"peak_voltage"});
  setVec(DoubleFeatureData, StringData, "AP_height",
         doubleFeatures.at("peak_voltage"));
  return doubleFeatures.at("peak_voltage").size();
}

// spike amplitude: peak_voltage - v[AP_begin_indices]
int SpikeShape::AP_amplitude(mapStr2intVec& IntFeatureData,
                        mapStr2doubleVec& DoubleFeatureData,
                        mapStr2Str& StringData) {
  const auto& doubleFeatures =
      getFeatures(DoubleFeatureData,
                  {"V", "stim_start", "stim_end", "peak_voltage", "peak_time"});
  const auto& intFeatures = getFeatures(IntFeatureData, {"AP_begin_indices"});
  if (doubleFeatures.at("peak_voltage").size() != doubleFeatures.at("peak_time").size())
    throw FeatureComputationError("AP_amplitude: Not the same amount of peak_time and peak_voltage entries");
  vector<double> peakvoltage_duringstim;
  for (size_t i = 0; i < doubleFeatures.at("peak_time").size(); i++) {
    if (doubleFeatures.at("peak_time")[i] >=
            doubleFeatures.at("stim_start")[0] &&
        doubleFeatures.at("peak_time")[i] <= doubleFeatures.at("stim_end")[0]) {
      peakvoltage_duringstim.push_back(doubleFeatures.at("peak_voltage")[i]);
    }
  }
  if (peakvoltage_duringstim.size() > intFeatures.at("AP_begin_indices").size())
    throw FeatureComputationError("AP_amplitude: More peak_voltage entries during the stimulus than AP_begin_indices entries");
  vector<double> apamplitude;
  apamplitude.resize(peakvoltage_duringstim.size());
  for (size_t i = 0; i < apamplitude.size(); i++) {
    apamplitude[i] =
        peakvoltage_duringstim[i] -
        doubleFeatures.at("V")[intFeatures.at("AP_begin_indices")[i]];
  }
  setVec(DoubleFeatureData, StringData, "AP_amplitude", apamplitude);
  return apamplitude.size();
}

// Amplitude of the first spike
int SpikeShape::AP1_amp(mapStr2intVec& IntFeatureData,
                   mapStr2doubleVec& DoubleFeatureData,
                   mapStr2Str& StringData) {
  const vector<double>& AP_amplitudes =
      getFeature(DoubleFeatureData, "AP_amplitude");
  vector<double> AP1_amp;
  AP1_amp.push_back(AP_amplitudes[0]);
  setVec(DoubleFeatureData, StringData, "AP1_amp", AP1_amp);
  return 1;
}

// Amplitude of the second spike
int SpikeShape::AP2_amp(mapStr2intVec& IntFeatureData,
                   mapStr2doubleVec& DoubleFeatureData,
                   mapStr2Str& StringData) {
  const vector<double>& AP_amplitudes =
      getFeature(DoubleFeatureData, "AP_amplitude");
  vector<double> AP2_amp;
  if (AP_amplitudes.size() < 2) {
    throw FeatureComputationError(
        "Size of AP_amplitude should be >= 2 for AP2_amp");
  }
  AP2_amp.push_back(AP_amplitudes[1]);
  setVec(DoubleFeatureData, StringData, "AP2_amp", AP2_amp);
  return 1;
}

int SpikeShape::mean_AP_amplitude(mapStr2intVec& IntFeatureData,
                             mapStr2doubleVec& DoubleFeatureData,
                             mapStr2Str& StringData) {
  const vector<double>& AP_amplitude =
      getFeature(DoubleFeatureData, "AP_amplitude");
  double mean_amp = 0.0;
  for (const auto& amplitude : AP_amplitude) {
    mean_amp += amplitude;
  }

  mean_amp /= AP_amplitude.size();
  vector<double> mean_AP_amplitude = {mean_amp};

  setVec(DoubleFeatureData, StringData, "mean_AP_amplitude", mean_AP_amplitude);

  return mean_AP_amplitude.size();
}

// Amplitude of the first spike
int SpikeShape::APlast_amp(mapStr2intVec& IntFeatureData,
                      mapStr2doubleVec& DoubleFeatureData,
                      mapStr2Str& StringData) {
  const vector<double>& AP_amplitudes =
      getFeature(DoubleFeatureData, "AP_amplitude");
  vector<double> APlast_amp;
  APlast_amp.push_back(AP_amplitudes[AP_amplitudes.size() - 1]);
  setVec(DoubleFeatureData, StringData, "APlast_amp", APlast_amp);
  return 1;
}

// *** AP_amplitude_change according to E22 ***
static int __AP_amplitude_change(const vector<double>& apamplitude,
                                 vector<double>& apamplitudechange) {
  if (apamplitude.size() < 1) {
    return -1;
  }
  apamplitudechange.resize(apamplitude.size() - 1);
  for (size_t i = 0; i < apamplitudechange.size(); i++) {
    apamplitudechange[i] =
        (apamplitude[i + 1] - apamplitude[0]) / apamplitude[0];
  }
  return apamplitudechange.size();
}
int SpikeShape::AP_amplitude_change(mapStr2intVec& IntFeatureData,
                               mapStr2doubleVec& DoubleFeatureData,
                               mapStr2Str& StringData) {
  const auto& features = getFeatures(DoubleFeatureData, {"AP_amplitude"});
  const vector<double>& apamplitude = features.at("AP_amplitude");

  vector<double> apamplitudechange;
  int retval = __AP_amplitude_change(apamplitude, apamplitudechange);

  if (retval > 0) {
    setVec(DoubleFeatureData, StringData, "AP_amplitude_change",
           apamplitudechange);
  }
  return retval;
}

// spike amplitude: peak_voltage - voltage_base
int SpikeShape::AP_amplitude_from_voltagebase(mapStr2intVec& IntFeatureData,
                                         mapStr2doubleVec& DoubleFeatureData,
                                         mapStr2Str& StringData) {
  const auto& doubleFeatures =
      getFeatures(DoubleFeatureData, {"voltage_base", "peak_voltage"});
  vector<double> apamplitude;
  for (const auto& peak : doubleFeatures.at("peak_voltage")) {
    apamplitude.push_back(peak - doubleFeatures.at("voltage_base")[0]);
  }
  setVec(DoubleFeatureData, StringData, "AP_amplitude_from_voltagebase",
         apamplitude);
  return apamplitude.size();
}

// Peak voltage of the first spike
int SpikeShape::AP1_peak(mapStr2intVec& IntFeatureData,
                    mapStr2doubleVec& DoubleFeatureData,
                    mapStr2Str& StringData) {
  const vector<double>& peak_voltage =
      getFeature(DoubleFeatureData, "peak_voltage");
  vector<double> AP1_peak;
  AP1_peak.push_back(peak_voltage[0]);
  setVec(DoubleFeatureData, StringData, "AP1_peak", AP1_peak);
  return 1;
}

// Peak voltage of the second spike
int SpikeShape::AP2_peak(mapStr2intVec& IntFeatureData,
                    mapStr2doubleVec& DoubleFeatureData,
                    mapStr2Str& StringData) {
  const vector<double>& peak_voltage =
      getFeature(DoubleFeatureData, "peak_voltage");
  vector<double> AP2_peak;
  if (peak_voltage.size() < 2) {
    throw FeatureComputationError(
        "Size of peak_voltage should be >= 2 for AP2_peak");
  }
  AP2_peak.push_back(peak_voltage[1]);
  setVec(DoubleFeatureData, StringData, "AP2_peak", AP2_peak);
  return 1;
}

// Difference amplitude of the second to first spike
int SpikeShape::AP2_AP1_diff(mapStr2intVec& IntFeatureData,
                        mapStr2doubleVec& DoubleFeatureData,
                        mapStr2Str& StringData) {
  const vector<double>& AP_amplitudes =
      getFeature(DoubleFeatureData, "AP_amplitude");
  vector<double> AP2_AP1_diff;
  if (AP_amplitudes.size() < 2) {
    throw FeatureComputationError(
        "Size of AP_amplitude should be >= 2 for AP2_AP1_diff");
  }
  AP2_AP1_diff.push_back(AP_amplitudes[1] - AP_amplitudes[0]);
  setVec(DoubleFeatureData, StringData, "AP2_AP1_diff", AP2_AP1_diff);
  return 1;
}

// Difference peak_amp of the second to first spike
int SpikeShape::AP2_AP1_peak_diff(mapStr2intVec& IntFeatureData,
                             mapStr2doubleVec& DoubleFeatureData,
                             mapStr2Str& StringData) {
  const vector<double>& peak_voltage =
      getFeature(DoubleFeatureData, "peak_voltage");
  vector<double> AP2_AP1_peak_diff;
  if (peak_voltage.size() < 2) {
    throw FeatureComputationError(
        "Size of peak_voltage should be >= 2 for AP2_AP1_peak_diff");
  }
  AP2_AP1_peak_diff.push_back(peak_voltage[1] - peak_voltage[0]);
  setVec(DoubleFeatureData, StringData, "AP2_AP1_peak_diff",
          AP2_AP1_peak_diff);
  return 1;
}

// *** amp_drop_first_second ***
static int __amp_drop_first_second(const vector<double>& peakvoltage,
                                   vector<double>& ampdropfirstsecond) {
  ampdropfirstsecond.push_back(peakvoltage[0] - peakvoltage[1]);
  return ampdropfirstsecond.size();
}
int SpikeShape::amp_drop_first_second(mapStr2intVec& IntFeatureData,
                                 mapStr2doubleVec& DoubleFeatureData,
                                 mapStr2Str& StringData) {
  const auto& features = getFeatures(DoubleFeatureData, {"peak_voltage"});
  const vector<double> peakvoltage = features.at("peak_voltage");

  if (peakvoltage.size() < 2) {
    throw FeatureComputationError("At least 2 spikes needed for calculation of amp_drop_first_second.");
  }
  vector<double> ampdropfirstsecond;
  int retval = __amp_drop_first_second(peakvoltage, ampdropfirstsecond);
  if (retval > 0) {
    setVec(DoubleFeatureData, StringData, "amp_drop_first_second",
           ampdropfirstsecond);
  }
  return retval;
}

// *** amp_drop_first_last ***
static int __amp_drop_first_last(const vector<double>& peakvoltage,
                                 vector<double>& ampdropfirstlast) {
  ampdropfirstlast.push_back(peakvoltage[0] - peakvoltage.back());
  return ampdropfirstlast.size();
}
int SpikeShape::amp_drop_first_last(mapStr2intVec& IntFeatureData,
                               mapStr2doubleVec& DoubleFeatureData,
                               mapStr2Str& StringData) {
  const auto& peakVoltageFeature =
      getFeatures(DoubleFeatureData, {"peak_voltage"});
  const vector<double>& peakvoltage = peakVoltageFeature.at("peak_voltage");

  if (peakvoltage.size() < 2) {
    throw FeatureComputationError("At least 2 spikes needed for calculation of amp_drop_first_last.");
  }
  vector<double> ampdropfirstlast;
  int retval = __amp_drop_first_last(peakvoltage, ampdropfirstlast);
  if (retval > 0) {
    setVec(DoubleFeatureData, StringData, "amp_drop_first_last",
           ampdropfirstlast);
  }
  return retval;
}
// end of amp_drop_first_last

// *** amp_drop_second_last ***
static int __amp_drop_second_last(const vector<double>& peakvoltage,
                                  vector<double>& ampdropsecondlast) {
  ampdropsecondlast.push_back(peakvoltage[1] - peakvoltage.back());
  return ampdropsecondlast.size();
}
int SpikeShape::amp_drop_second_last(mapStr2intVec& IntFeatureData,
                                mapStr2doubleVec& DoubleFeatureData,
                                mapStr2Str& StringData) {
  const auto& peakVoltageFeatures =
      getFeatures(DoubleFeatureData, {"peak_voltage"});
  const vector<double>& peakvoltage = peakVoltageFeatures.at("peak_voltage");
  // Ensure there are at least 3 spikes for calculation
  if (peakvoltage.size() < 3) {
    throw FeatureComputationError("At least 3 spikes needed for calculation of amp_drop_second_last.");
  }
  vector<double> ampdropsecondlast;
  int retval = __amp_drop_second_last(peakvoltage, ampdropsecondlast);
  if (retval > 0) {
    setVec(DoubleFeatureData, StringData, "amp_drop_second_last",
           ampdropsecondlast);
  }
  return retval;
}

// *** max_amp_difference ***
static int __max_amp_difference(const vector<double>& peakvoltage,
                                vector<double>& maxampdifference) {
  vector<double> diff_peak_voltage;
  if (peakvoltage.size() < 1) {
    return -1;
  }
  diff_peak_voltage.resize(peakvoltage.size() - 1);
  for (size_t i = 0; i < diff_peak_voltage.size(); i++) {
    diff_peak_voltage[i] = peakvoltage[i] - peakvoltage[i + 1];
  }
  maxampdifference.push_back(
      *max_element(diff_peak_voltage.begin(), diff_peak_voltage.end()));
  return maxampdifference.size();
}

int SpikeShape::max_amp_difference(mapStr2intVec& IntFeatureData,
                              mapStr2doubleVec& DoubleFeatureData,
                              mapStr2Str& StringData) {
  const auto& features = getFeatures(DoubleFeatureData, {"peak_voltage"});

  // Ensure there are at least 2 spikes for calculation
  if (features.at("peak_voltage").size() < 2) {
    throw FeatureComputationError("At least 2 spikes needed for calculation of max_amp_difference.");
  }
  vector<double> maxampdifference;
  int retval =
      __max_amp_difference(features.at("peak_voltage"), maxampdifference);

  if (retval > 0) {
    setVec(DoubleFeatureData, StringData, "max_amp_difference",
           maxampdifference);
  }
  return retval;
}

// *** AP_amplitude_diff based on AP_amplitude_change but not normalized  ***
static int __AP_amplitude_diff(const vector<double>& apamplitude,
                               vector<double>& apamplitudediff) {
  if (apamplitude.size() <= 1) return -1;
  apamplitudediff.resize(apamplitude.size() - 1);
  for (size_t i = 0; i < apamplitudediff.size(); i++) {
    apamplitudediff[i] = (apamplitude[i + 1] - apamplitude[i]);
  }
  return apamplitudediff.size();
}

int SpikeShape::AP_amplitude_diff(mapStr2intVec& IntFeatureData,
                             mapStr2doubleVec& DoubleFeatureData,
                             mapStr2Str& StringData) {
  const auto& doubleFeatures = getFeatures(DoubleFeatureData, {"AP_amplitude"});
  vector<double> apamplitudediff;
  int retval =
      __AP_amplitude_diff(doubleFeatures.at("AP_amplitude"), apamplitudediff);
  if (retval > 0) {
    setVec(DoubleFeatureData, StringData, "AP_amplitude_diff", apamplitudediff);
  }
  return retval;
}

static int __min_AHP_indices(const vector<double>& t, const vector<double>& v,
                             const vector<int>& peak_indices,
                             const double stim_start, const double stim_end,
                             const bool strict_stiminterval,
                             vector<int>& min_ahp_indices,
                             vector<double>& min_ahp_values) {
  vector<int> peak_indices_plus = peak_indices;
  unsigned end_index = 0;

  if (strict_stiminterval) {
    end_index = distance(t.begin(),
                         find_if(t.begin(), t.end(), [stim_end](double t_val) {
                           return t_val >= stim_end;
                         }));
  } else {
    end_index = distance(t.begin(), t.end());
  }

  size_t ahpindex = 0;

  peak_indices_plus.push_back(end_index);

  for (size_t i = 0; i < peak_indices_plus.size() - 1; i++) {
    ahpindex = distance(
        v.begin(), first_min_element(v.begin() + peak_indices_plus[i],
                                     v.begin() + peak_indices_plus[i + 1]));

    if (ahpindex != end_index - 1) {
      min_ahp_indices.push_back(ahpindex);

      EFEL_ASSERT(ahpindex < v.size(),
                  "AHP index falls outside of voltage array");
      min_ahp_values.push_back(v[ahpindex]);
    }
  }

  return min_ahp_indices.size();
}

// min_AHP_indices
// find the first minimum between two spikes,
// and the first minimum between the last spike and the time the stimulus ends
int SpikeShape::min_AHP_indices(mapStr2intVec& IntFeatureData,
                           mapStr2doubleVec& DoubleFeatureData,
                           mapStr2Str& StringData) {
  int retVal;
  // Retrieve all required double features at once
  const auto& doubleFeatures =
      getFeatures(DoubleFeatureData, {"V", "T", "stim_start", "stim_end"});
  const auto& intFeatures = getFeatures(IntFeatureData, {"peak_indices"});

  vector<int> min_ahp_indices;
  vector<double> min_ahp_values;
  double stim_start = doubleFeatures.at("stim_start")[0];
  double stim_end = doubleFeatures.at("stim_end")[0];
  // Get strict_stiminterval
  vector<int> strict_stiminterval_vec;
  bool strict_stiminterval;
  retVal =
      getParam(IntFeatureData, "strict_stiminterval", strict_stiminterval_vec);
  if (retVal <= 0) {
    strict_stiminterval = false;
  } else {
    strict_stiminterval = bool(strict_stiminterval_vec[0]);
  }

  retVal =
      __min_AHP_indices(doubleFeatures.at("T"), doubleFeatures.at("V"),
                        intFeatures.at("peak_indices"), stim_start, stim_end,
                        strict_stiminterval, min_ahp_indices, min_ahp_values);

  if (retVal == 0) return -1;
  if (retVal > 0) {
    setVec(IntFeatureData, StringData, "min_AHP_indices", min_ahp_indices);
    setVec(DoubleFeatureData, StringData, "min_AHP_values", min_ahp_values);
  }
  return retVal;
}

int SpikeShape::min_AHP_values(mapStr2intVec& IntFeatureData,
                          mapStr2doubleVec& DoubleFeatureData,
                          mapStr2Str& StringData) {
  const auto& intFeatures = getFeatures(IntFeatureData, {"min_AHP_indices"});
  int retVal = intFeatures.at("min_AHP_indices").size();
  if (retVal <= 0) return -1;
  return retVal;
}

// Difference with SpikeShape is that this function doesn't return -1 if there are no
// min_AHP_values
int SpikeShape::AHP_depth_abs(mapStr2intVec& IntFeatureData,
                         mapStr2doubleVec& DoubleFeatureData,
                         mapStr2Str& StringData) {
  const auto& vAHP = getFeature(DoubleFeatureData, "min_AHP_values");
  setVec(DoubleFeatureData, StringData, "AHP_depth_abs", vAHP);
  return vAHP.size();
}

// *** AHP_depth_abs_slow ***
// same as AHP_depth_abs but the minimum search starts
// 5 ms (or custom duration) after the spike,
// first ISI is ignored
static int __AHP_depth_abs_slow_indices(const vector<double>& t,
                                        const vector<double>& v,
                                        const vector<int>& peakindices,
                                        double sahp_start,
                                        vector<int>& adas_indices) {
  adas_indices.resize(peakindices.size() - 2);
  for (size_t i = 0; i < adas_indices.size(); i++) {
    // start 5 ms (or custom duration) after last spike
    double t_start = t[peakindices[i + 1]] + sahp_start;
    adas_indices[i] = distance(
        v.begin(),
        min_element(v.begin() + distance(t.begin(),
                                         find_if(t.begin() + peakindices[i + 1],
                                                 t.begin() + peakindices[i + 2],
                                                 [t_start](double val) {
                                                   return val >= t_start;
                                                 })),
                    v.begin() + peakindices[i + 2]));
  }
  return adas_indices.size();
}

int SpikeShape::AHP_depth_abs_slow(mapStr2intVec& IntFeatureData,
                              mapStr2doubleVec& DoubleFeatureData,
                              mapStr2Str& StringData) {
  const auto& doubleFeatures =
      getFeatures(DoubleFeatureData, {"T", "V", "sahp_start"});
  const auto& intFeatures = getFeatures(IntFeatureData, {"peak_indices"});
  if (intFeatures.at("peak_indices").size() < 3)
    throw FeatureComputationError("At least 3 spikes needed for AHP_depth_abs_slow and AHP_slow_time.");
  double sahp_start = (doubleFeatures.at("sahp_start").empty())
                          ? 5
                          : doubleFeatures.at("sahp_start")[0];
  vector<int> adas_indices;
  int retval = __AHP_depth_abs_slow_indices(
      doubleFeatures.at("T"), doubleFeatures.at("V"),
      intFeatures.at("peak_indices"), sahp_start, adas_indices);
  vector<double> ahpdepthabsslow(adas_indices.size());
  vector<double> ahpslowtime(adas_indices.size());
  for (size_t i = 0; i < adas_indices.size(); i++) {
    ahpdepthabsslow[i] = doubleFeatures.at("V")[adas_indices[i]];
    ahpslowtime[i] =
        (doubleFeatures.at("T")[adas_indices[i]] -
         doubleFeatures.at("T")[intFeatures.at("peak_indices")[i + 1]]) /
        (doubleFeatures.at("T")[intFeatures.at("peak_indices")[i + 2]] -
         doubleFeatures.at("T")[intFeatures.at("peak_indices")[i + 1]]);
  }
  if (retval >= 0) {
    setVec(DoubleFeatureData, StringData, "AHP_depth_abs_slow",
           ahpdepthabsslow);
    setVec(DoubleFeatureData, StringData, "AHP_slow_time", ahpslowtime);
  }
  return retval;
}

int SpikeShape::AHP_slow_time(mapStr2intVec& IntFeatureData,
                         mapStr2doubleVec& DoubleFeatureData,
                         mapStr2Str& StringData) {
  const auto& doubleFeatures = getFeatures(DoubleFeatureData, {"AHP_depth_abs_slow"});
  int retVal = doubleFeatures.at("AHP_depth_abs_slow").size();
  if (retVal <= 0) return -1;
  return retVal;
}

// *** AHP_depth_slow ***
static int __AHP_depth_slow(const vector<double>& voltagebase,
                            const vector<double>& minahpvalues,
                            vector<double>& ahpdepth) {
  for (size_t i = 0; i < minahpvalues.size(); i++) {
    ahpdepth.push_back(minahpvalues[i] - voltagebase[0]);
  }
  return ahpdepth.size();
}

int SpikeShape::AHP_depth_slow(mapStr2intVec& IntFeatureData,
                          mapStr2doubleVec& DoubleFeatureData,
                          mapStr2Str& StringData) {
  const auto& doubleFeatures =
      getFeatures(DoubleFeatureData, {"voltage_base", "AHP_depth_abs_slow"});
  vector<double> ahpdepthslow;
  int retval =
      __AHP_depth_slow(doubleFeatures.at("voltage_base"),
                       doubleFeatures.at("AHP_depth_abs_slow"), ahpdepthslow);
  if (retval > 0) {
    setVec(DoubleFeatureData, StringData, "AHP_depth_slow", ahpdepthslow);
  }
  return retval;
}

// *** AHP_depth ***
static int __AHP_depth(const vector<double>& voltagebase,
                       const vector<double>& minahpvalues,
                       vector<double>& ahpdepth) {
  for (size_t i = 0; i < minahpvalues.size(); i++) {
    ahpdepth.push_back(minahpvalues[i] - voltagebase[0]);
  }
  return ahpdepth.size();
}
int SpikeShape::AHP_depth(mapStr2intVec& IntFeatureData,
                     mapStr2doubleVec& DoubleFeatureData,
                     mapStr2Str& StringData) {
  const auto& doubleFeatures =
      getFeatures(DoubleFeatureData, {"voltage_base", "min_AHP_values"});
  vector<double> ahpdepth;
  int retval = __AHP_depth(doubleFeatures.at("voltage_base"),
                           doubleFeatures.at("min_AHP_values"), ahpdepth);
  if (retval > 0) {
    setVec(DoubleFeatureData, StringData, "AHP_depth", ahpdepth);
  }
  return retval;
}

// *** AHP_depth_diff, returns AHP_depth[i+1] - AHP_depth[i]  ***
static int __AHP_depth_diff(const vector<double>& ahpdepth,
                            vector<double>& ahpdepthdiff) {
  if (ahpdepth.size() <= 1) return -1;
  ahpdepthdiff.resize(ahpdepth.size() - 1);
  for (size_t i = 0; i < ahpdepthdiff.size(); i++) {
    ahpdepthdiff[i] = (ahpdepth[i + 1] - ahpdepth[i]);
  }
  return ahpdepthdiff.size();
}
int SpikeShape::AHP_depth_diff(mapStr2intVec& IntFeatureData,
                          mapStr2doubleVec& DoubleFeatureData,
                          mapStr2Str& StringData) {
  const auto& doubleFeatures = getFeatures(DoubleFeatureData, {"AHP_depth"});
  vector<double> ahpdepthdiff;
  int retval = __AHP_depth_diff(doubleFeatures.at("AHP_depth"), ahpdepthdiff);
  if (retval > 0) {
    setVec(DoubleFeatureData, StringData, "AHP_depth_diff", ahpdepthdiff);
  }
  return retval;
}

// *** fast_AHP according to E13 and E21 ***
static int __fast_AHP(const vector<double>& v,
                      const vector<int>& apbeginindices,
                      const vector<int>& minahpindices,
                      vector<double>& fastahp) {
  if (apbeginindices.size() < 1) {
    return -1;
  }
  fastahp.resize(apbeginindices.size() - 1);
  for (size_t i = 0; i < fastahp.size(); i++) {
    fastahp[i] = v[apbeginindices[i]] - v[minahpindices[i]];
  }
  return fastahp.size();
}
int SpikeShape::fast_AHP(mapStr2intVec& IntFeatureData,
                    mapStr2doubleVec& DoubleFeatureData,
                    mapStr2Str& StringData) {
  const auto& doubleFeatures = getFeatures(DoubleFeatureData, {"V"});
  const auto& intFeatures =
      getFeatures(IntFeatureData, {"AP_begin_indices", "min_AHP_indices"});

  const vector<double>& v = doubleFeatures.at("V");
  const vector<int>& apbeginindices = intFeatures.at("AP_begin_indices");
  const vector<int>& minahpindices = intFeatures.at("min_AHP_indices");

  vector<double> fastahp;
  int retval = __fast_AHP(v, apbeginindices, minahpindices, fastahp);

  if (retval > 0) {
    setVec(DoubleFeatureData, StringData, "fast_AHP", fastahp);
  }
  return retval;
}

// *** fast_AHP_change according to E27 ***
static int __fast_AHP_change(const vector<double>& fastahp,
                             vector<double>& fastahpchange) {
  if (fastahp.size() < 1) {
    return -1;
  }
  fastahpchange.resize(fastahp.size() - 1);
  for (size_t i = 0; i < fastahpchange.size(); i++) {
    fastahpchange[i] = (fastahp[i + 1] - fastahp[0]) / fastahp[0];
  }
  return fastahpchange.size();
}
int SpikeShape::fast_AHP_change(mapStr2intVec& IntFeatureData,
                           mapStr2doubleVec& DoubleFeatureData,
                           mapStr2Str& StringData) {
  const auto& features = getFeatures(DoubleFeatureData, {"fast_AHP"});
  const vector<double>& fastahp = features.at("fast_AHP");
  vector<double> fastahpchange;
  int retval = __fast_AHP_change(fastahp, fastahpchange);
  if (retval > 0) {
    setVec(DoubleFeatureData, StringData, "fast_AHP_change", fastahpchange);
  }
  return retval;
}

static int __AHP_depth_from_peak(const vector<double>& v,
                                 const vector<int>& peakIndices,
                                 const vector<int>& minAHPIndices,
                                 vector<double>& ahpDepthFromPeak) {
  if (peakIndices.size() < minAHPIndices.size()) return -1;
  for (size_t i = 0; i < minAHPIndices.size(); i++) {
    ahpDepthFromPeak.push_back(v[peakIndices[i]] - v[minAHPIndices[i]]);
  }
  return ahpDepthFromPeak.size();
}

int SpikeShape::AHP_depth_from_peak(mapStr2intVec& IntFeatureData,
                               mapStr2doubleVec& DoubleFeatureData,
                               mapStr2Str& StringData) {
  const auto& doubleFeatures = getFeatures(DoubleFeatureData, {"V"});
  const auto& intFeatures =
      getFeatures(IntFeatureData, {"peak_indices", "min_AHP_indices"});
  vector<double> ahpDepthFromPeak;
  const vector<double>& V = doubleFeatures.at("V");
  const vector<int>& peakIndices = intFeatures.at("peak_indices");
  const vector<int>& minAHPIndices = intFeatures.at("min_AHP_indices");
  int retVal =
      __AHP_depth_from_peak(V, peakIndices, minAHPIndices, ahpDepthFromPeak);
  if (retVal > 0) {
    setVec(DoubleFeatureData, StringData, "AHP_depth_from_peak",
           ahpDepthFromPeak);
  }
  return retVal;
}

// AHP_depth_from_peak of first spike
int SpikeShape::AHP1_depth_from_peak(mapStr2intVec& IntFeatureData,
                                mapStr2doubleVec& DoubleFeatureData,
                                mapStr2Str& StringData) {
  const vector<double>& ahpDepthFromPeak =
      getFeature(DoubleFeatureData, "AHP_depth_from_peak");
  vector<double> ahp1DepthFromPeak;
  ahp1DepthFromPeak.push_back(ahpDepthFromPeak[0]);
  setVec(DoubleFeatureData, StringData, "AHP1_depth_from_peak",
         ahp1DepthFromPeak);
  return 1;
}

// AHP_depth_from_peak of second spike
int SpikeShape::AHP2_depth_from_peak(mapStr2intVec& IntFeatureData,
                                mapStr2doubleVec& DoubleFeatureData,
                                mapStr2Str& StringData) {
  const vector<double>& ahpDepthFromPeak =
      getFeature(DoubleFeatureData, "AHP_depth_from_peak");
  vector<double> ahp2DepthFromPeak;
  if (ahpDepthFromPeak.size() < 2) {
    throw FeatureComputationError(
        "Size of AHP_depth_from_peak should be >= 2 for AHP2_depth_from_peak");
  }
  ahp2DepthFromPeak.push_back(ahpDepthFromPeak[1]);
  setVec(DoubleFeatureData, StringData, "AHP2_depth_from_peak",
         ahp2DepthFromPeak);
  return 1;
}

static int __AHP_time_from_peak(const vector<double>& t,
                                const vector<int>& peakIndices,
                                const vector<int>& minAHPIndices,
                                vector<double>& ahpTimeFromPeak) {
  for (size_t i = 0; i < peakIndices.size() && i < minAHPIndices.size(); i++) {
    ahpTimeFromPeak.push_back(t[minAHPIndices[i]] - t[peakIndices[i]]);
  }
  return ahpTimeFromPeak.size();
}

int SpikeShape::AHP_time_from_peak(mapStr2intVec& IntFeatureData,
                              mapStr2doubleVec& DoubleFeatureData,
                              mapStr2Str& StringData) {
  const auto& doubleFeatures = getFeatures(DoubleFeatureData, {"T"});
  const auto& intFeatures =
      getFeatures(IntFeatureData, {"peak_indices", "min_AHP_indices"});
  vector<double> ahpTimeFromPeak;
  const vector<double>& T = doubleFeatures.at("T");
  const vector<int>& peakIndices = intFeatures.at("peak_indices");
  const vector<int>& minAHPIndices = intFeatures.at("min_AHP_indices");
  int retVal =
      __AHP_time_from_peak(T, peakIndices, minAHPIndices, ahpTimeFromPeak);
  if (retVal > 0) {
    setVec(DoubleFeatureData, StringData, "AHP_time_from_peak",
           ahpTimeFromPeak);
  }
  return retVal;
}

// strict_stiminterval should be True when using this feature
static int __ADP_peak_indices(const vector<double>& v,
                              const vector<int>& min_AHP_indices,
                              const vector<int>& min_between_peaks_indices,
                              vector<int>& ADP_peak_indices,
                              vector<double>& ADP_peak_values) {
  if (min_AHP_indices.size() > min_between_peaks_indices.size()) {
    throw FeatureComputationError("min_AHP_indices should not have less elements than min_between_peaks_indices");
  }

  unsigned adp_peak_index;
  for (size_t i = 0; i < min_AHP_indices.size(); i++) {
    adp_peak_index = max_element(v.begin() + min_AHP_indices[i],
                                 v.begin() + min_between_peaks_indices[i]) -
                     v.begin();

    ADP_peak_indices.push_back(adp_peak_index);
    ADP_peak_values.push_back(v[adp_peak_index]);
  }

  return ADP_peak_indices.size();
}

// strict_stiminterval should be True when using this feature
int SpikeShape::ADP_peak_indices(mapStr2intVec& IntFeatureData,
                            mapStr2doubleVec& DoubleFeatureData,
                            mapStr2Str& StringData) {
  const auto& doubleFeatures = getFeatures(DoubleFeatureData, {"V"});
  const auto& intFeatures = getFeatures(
      IntFeatureData, {"min_AHP_indices", "min_between_peaks_indices"});
  vector<int> ADP_peak_indices;
  vector<double> ADP_peak_values;
  int retVal = __ADP_peak_indices(doubleFeatures.at("V"),
                                  intFeatures.at("min_AHP_indices"),
                                  intFeatures.at("min_between_peaks_indices"),
                                  ADP_peak_indices, ADP_peak_values);
  if (retVal > 0) {
    setVec(IntFeatureData, StringData, "ADP_peak_indices", ADP_peak_indices);
    setVec(DoubleFeatureData, StringData, "ADP_peak_values", ADP_peak_values);
  }
  return retVal;
}

// strict_stiminterval should be True when using this feature
int SpikeShape::ADP_peak_values(mapStr2intVec& IntFeatureData,
                           mapStr2doubleVec& DoubleFeatureData,
                           mapStr2Str& StringData) {
  const auto& intFeatures = getFeatures(IntFeatureData, {"ADP_peak_indices"});
  int retVal = intFeatures.at("ADP_peak_indices").size();
  if (retVal <= 0) return -1;
  return retVal;
}

// strict_stiminterval should be True when using this feature
int SpikeShape::ADP_peak_amplitude(mapStr2intVec& IntFeatureData,
                              mapStr2doubleVec& DoubleFeatureData,
                              mapStr2Str& StringData) {
  const auto& doubleFeatures =
      getFeatures(DoubleFeatureData, {"min_AHP_values", "ADP_peak_values"});
  vector<double> ADP_peak_amplitude;
  const vector<double>& min_AHP_values = doubleFeatures.at("min_AHP_values");
  const vector<double>& ADP_peak_values = doubleFeatures.at("ADP_peak_values");

  if (min_AHP_values.size() != ADP_peak_values.size())
    throw FeatureComputationError("min_AHP_values and ADP_peak_values should have the same number of elements");
  for (size_t i = 0; i < ADP_peak_values.size(); i++) {
    ADP_peak_amplitude.push_back(ADP_peak_values[i] - min_AHP_values[i]);
  }
  setVec(DoubleFeatureData, StringData, "ADP_peak_amplitude",
         ADP_peak_amplitude);
  return ADP_peak_amplitude.size();
}

static int __depolarized_base(const vector<double>& t, const vector<double>& v,
                              double stimstart, double stimend,
                              const vector<int>& apbi,
                              const vector<int>& apendi,
                              vector<double>& dep_base) {
  int i, n, k, startIndex, endIndex, nPt;
  double baseValue;
  // to make sure it access minimum index of both length
  n = std::min(apendi.size(), apbi.size());

  if (n > 2) {
    dep_base.clear();
    for (i = 0; i < n - 1; i++) {
      nPt = 0;
      baseValue = 0;
      startIndex = apendi[i];
      endIndex = apbi[i + 1];
      for (k = startIndex; k < endIndex; k++) {
        if (k >= 0 && k < v.size()) {
          baseValue += v[k];
          ++nPt;
        }
      }
      if (nPt > 0) {
        baseValue /= nPt;
        dep_base.push_back(baseValue);
      }
    }
    return dep_base.size();
  }
  return -1;
}

int SpikeShape::depolarized_base(mapStr2intVec& IntFeatureData,
                            mapStr2doubleVec& DoubleFeatureData,
                            mapStr2Str& StringData) {
  // Retrieve all required double and int features at once
  const auto& doubleFeatures =
      getFeatures(DoubleFeatureData, {"T", "V", "stim_start", "stim_end"});
  const auto& intFeatures =
      getFeatures(IntFeatureData, {"AP_end_indices", "AP_begin_indices"});

  vector<double> dep_base;
  int retVal = __depolarized_base(
      doubleFeatures.at("T"), doubleFeatures.at("V"),
      doubleFeatures.at("stim_start").front(),
      doubleFeatures.at("stim_end").front(), intFeatures.at("AP_begin_indices"),
      intFeatures.at("AP_end_indices"), dep_base);

  if (retVal > 0) {
    setVec(DoubleFeatureData, StringData, "depolarized_base", dep_base);
  }
  return retVal;
}

// min_voltage_between_spikes: minimal voltage between consecutive spikes
int SpikeShape::min_voltage_between_spikes(mapStr2intVec& IntFeatureData,
                                      mapStr2doubleVec& DoubleFeatureData,
                                      mapStr2Str& StringData) {
  const auto& doubleFeatures = getFeatures(DoubleFeatureData, {"V"});
  const auto& intFeatures = getFeatures(IntFeatureData, {"peak_indices"});

  if (intFeatures.at("peak_indices").size() < 2) {
    throw FeatureComputationError(
        "Size of peak_indices should be >= 2 for min_voltage_between_spikes");
  }

  vector<double> min_voltage_between_spikes;
  for (size_t i = 0; i < intFeatures.at("peak_indices").size() - 1; i++) {
    min_voltage_between_spikes.push_back(*min_element(
        doubleFeatures.at("V").begin() + intFeatures.at("peak_indices")[i],
        doubleFeatures.at("V").begin() +
            intFeatures.at("peak_indices")[i + 1]));
  }

  setVec(DoubleFeatureData, StringData, "min_voltage_between_spikes",
         min_voltage_between_spikes);
  return min_voltage_between_spikes.size();
}

static int __min_between_peaks_indices(
    const vector<double>& t, const vector<double>& v,
    const vector<int>& peak_indices, const double stim_start,
    const double stim_end, const bool strict_stiminterval,
    vector<int>& min_btw_peaks_indices, vector<double>& min_btw_peaks_values) {
  vector<int> peak_indices_plus = peak_indices;
  unsigned end_index = 0;

  if (strict_stiminterval) {
    end_index = distance(t.begin(),
                         find_if(t.begin(), t.end(), [stim_end](double t_val) {
                           return t_val >= stim_end;
                         }));
  } else {
    end_index = distance(t.begin(), t.end());
  }

  size_t minindex = 0;

  peak_indices_plus.push_back(end_index);

  for (size_t i = 0; i < peak_indices_plus.size() - 1; i++) {
    minindex = distance(v.begin(),
                        std::min_element(v.begin() + peak_indices_plus[i],
                                         v.begin() + peak_indices_plus[i + 1]));

    min_btw_peaks_indices.push_back(minindex);

    EFEL_ASSERT(minindex < v.size(),
                "min index falls outside of voltage array");
    min_btw_peaks_values.push_back(v[minindex]);
  }

  return min_btw_peaks_indices.size();
}

// min_between_peaks_indices
// find the minimum between two spikes,
// and the minimum between the last spike and the time the stimulus ends.
// This is different from min_AHP_indices, for traces that have broad peaks
// with local minima and maxima in it (e.g. dendritic AP)
int SpikeShape::min_between_peaks_indices(mapStr2intVec& IntFeatureData,
                                     mapStr2doubleVec& DoubleFeatureData,
                                     mapStr2Str& StringData) {
  // Retrieve all required double and int features at once
  const auto& doubleFeatures =
      getFeatures(DoubleFeatureData, {"T", "V", "stim_start", "stim_end"});
  const auto& intFeatures =
      getFeatures(IntFeatureData, {"peak_indices", "strict_stiminterval"});

  vector<int> min_btw_peaks_indices;
  vector<double> min_btw_peaks_values;
  int retVal = __min_between_peaks_indices(
      doubleFeatures.at("T"), doubleFeatures.at("V"),
      intFeatures.at("peak_indices"), doubleFeatures.at("stim_start").front(),
      doubleFeatures.at("stim_end").front(),
      intFeatures.at("strict_stiminterval").empty()
          ? false
          : intFeatures.at("strict_stiminterval").front(),
      min_btw_peaks_indices, min_btw_peaks_values);

  if (retVal > 0) {
    setVec(IntFeatureData, StringData, "min_between_peaks_indices",
           min_btw_peaks_indices);
    setVec(DoubleFeatureData, StringData, "min_between_peaks_values",
           min_btw_peaks_values);
  }
  return retVal;
}

int SpikeShape::min_between_peaks_values(mapStr2intVec& IntFeatureData,
                                    mapStr2doubleVec& DoubleFeatureData,
                                    mapStr2Str& StringData) {
  const auto& intFeatures = getFeatures(IntFeatureData, {"min_between_peaks_indices"});
  int retVal = intFeatures.at("min_between_peaks_indices").size();
  if (retVal <= 0) return -1;
  return retVal;
}

// *** AP_duration_half_width according to E8 and E16 ***
static int __AP_duration_half_width(const vector<double>& t,
                                    const vector<int>& apriseindices,
                                    const vector<int>& apfallindices,
                                    vector<double>& apdurationhalfwidth) {
  apdurationhalfwidth.resize(apriseindices.size());
  for (size_t i = 0; i < apdurationhalfwidth.size(); i++) {
    apdurationhalfwidth[i] = t[apfallindices[i]] - t[apriseindices[i]];
  }
  return apdurationhalfwidth.size();
}
int SpikeShape::AP_duration_half_width(mapStr2intVec& IntFeatureData,
                                  mapStr2doubleVec& DoubleFeatureData,
                                  mapStr2Str& StringData) {
  // Fetching all required features in one go.
  const auto& doubleFeatures = getFeatures(DoubleFeatureData, {"T"});
  const auto& intFeatures =
      getFeatures(IntFeatureData, {"AP_rise_indices", "AP_fall_indices"});

  vector<double> apdurationhalfwidth;
  int retval = __AP_duration_half_width(
      doubleFeatures.at("T"), intFeatures.at("AP_rise_indices"),
      intFeatures.at("AP_fall_indices"), apdurationhalfwidth);
  if (retval > 0) {
    setVec(DoubleFeatureData, StringData, "AP_duration_half_width",
           apdurationhalfwidth);
  }
  return retval;
}

// *** AP_duration_half_width_change according to E24 ***
static int __AP_duration_half_width_change(
    const vector<double>& apdurationhalfwidth,
    vector<double>& apdurationhalfwidthchange) {
  if (apdurationhalfwidth.size() < 1) {
    return -1;
  }
  apdurationhalfwidthchange.resize(apdurationhalfwidth.size() - 1);
  for (size_t i = 0; i < apdurationhalfwidthchange.size(); i++) {
    apdurationhalfwidthchange[i] =
        (apdurationhalfwidth[i + 1] - apdurationhalfwidth[0]) /
        apdurationhalfwidth[0];
  }
  return apdurationhalfwidthchange.size();
}
int SpikeShape::AP_duration_half_width_change(mapStr2intVec& IntFeatureData,
                                         mapStr2doubleVec& DoubleFeatureData,
                                         mapStr2Str& StringData) {
  const auto& features =
      getFeatures(DoubleFeatureData, {"AP_duration_half_width"});
  const vector<double>& apdurationhalfwidth =
      features.at("AP_duration_half_width");
  vector<double> apdurationhalfwidthchange;
  int retval = __AP_duration_half_width_change(apdurationhalfwidth,
                                               apdurationhalfwidthchange);
  if (retval > 0) {
    setVec(DoubleFeatureData, StringData, "AP_duration_half_width_change",
           apdurationhalfwidthchange);
  }
  return retval;
}

// *** AP_width ***
//
// spike width calculation according to threshold value for the first spike
// unfortunately spike width means the width of the spike on onset not at half
// maximum
static int __AP_width(const vector<double>& t, const vector<double>& v,
                      double stimstart, double stimend, double threshold,
                      const vector<int>& peakindices,
                      const vector<int>& minahpindices,
                      const bool strict_stiminterval, vector<double>& apwidth) {
  //   printf("\n Inside AP_width...\n");
  //   printf("\nStimStart = %f , thereshold = %f ", stimstart, threshold);
  vector<int> indices;
  if (strict_stiminterval) {
    int start_index = distance(
        t.begin(), find_if(t.begin(), t.end(),
                           [stimstart](double x) { return x >= stimstart; }));
    int end_index = distance(
        t.begin(), find_if(t.begin(), t.end(),
                           [stimend](double x) { return x >= stimend; }));
    indices.push_back(start_index);
    for (size_t i = 0; i < minahpindices.size(); i++) {
      if (start_index < minahpindices[i] && minahpindices[i] < end_index) {
        indices.push_back(minahpindices[i]);
      }
    }
  } else {
    indices.push_back(0);
    for (size_t i = 0; i < minahpindices.size(); i++) {
      indices.push_back(minahpindices[i]);
    }
  }
  for (size_t i = 0; i < indices.size() - 1; i++) {
    /*
    // FWHM (not used):
    // half maximum
    double v_hm = (v[peakindices[i]] + threshold) / 2.;
    // half maximum indices
    int hm_index1 = distance(v.begin(), find_if(v.begin() + indices[i],
    v.begin() + indices[i + 1], bind2nd(greater_equal<double>(), v_hm)));
    int hm_index2 = distance(v.begin(), find_if(v.begin() + peakindices[i],
    v.begin() + indices[i + 1], bind2nd(less_equal<double>(), v_hm)));
    apwidth.push_back(t[hm_index2] - t[hm_index1]);
    */
    auto onset_index = distance(
        v.begin(), find_if(v.begin() + indices[i], v.begin() + indices[i + 1],
                           [threshold](double x) { return x >= threshold; }));
    // int end_index = distance(v.begin(), find_if(v.begin() + peakindices[i],
    // v.begin() + indices[i + 1], bind2nd(less_equal<double>(), threshold)));
    auto end_index = distance(
        v.begin(), find_if(v.begin() + onset_index, v.begin() + indices[i + 1],
                           [threshold](double x) { return x <= threshold; }));
    apwidth.push_back(t[end_index] - t[onset_index]);
  }
  return apwidth.size();
}

int SpikeShape::AP_width(mapStr2intVec& IntFeatureData,
                    mapStr2doubleVec& DoubleFeatureData,
                    mapStr2Str& StringData) {
  const auto& doubleFeatures = getFeatures(
      DoubleFeatureData, {"T", "V", "Threshold", "stim_start", "stim_end"});
  const auto& intFeatures =
      getFeatures(IntFeatureData,
                  {"peak_indices", "min_AHP_indices", "strict_stiminterval"});
  bool strict_stiminterval =
      intFeatures.at("strict_stiminterval").size() > 0
          ? bool(intFeatures.at("strict_stiminterval")[0])
          : false;
  vector<double> apwidth;
  int retval = __AP_width(
      doubleFeatures.at("T"), doubleFeatures.at("V"),
      doubleFeatures.at("stim_start")[0], doubleFeatures.at("stim_end")[0],
      doubleFeatures.at("Threshold")[0], intFeatures.at("peak_indices"),
      intFeatures.at("min_AHP_indices"), strict_stiminterval, apwidth);
  if (retval > 0) {
    setVec(DoubleFeatureData, StringData, "AP_width", apwidth);
  }
  return retval;
}

// *** AP_duration according to E7 and E15 ***
static int __AP_duration(const vector<double>& t,
                         const vector<int>& apbeginindices,
                         const vector<int>& endindices,
                         vector<double>& apduration) {
  apduration.resize(std::min(apbeginindices.size(), endindices.size()));
  for (size_t i = 0; i < apduration.size(); i++) {
    // printf("%d, %d, %d\n", t.size(), apbeginindices.size(),
    // endindices.size())
    apduration[i] = t[endindices[i]] - t[apbeginindices[i]];
  }
  return apduration.size();
}
int SpikeShape::AP_duration(mapStr2intVec& IntFeatureData,
                       mapStr2doubleVec& DoubleFeatureData,
                       mapStr2Str& StringData) {
  const auto& doubleFeatures = getFeatures(DoubleFeatureData, {"T"});
  const auto& intFeatures =
      getFeatures(IntFeatureData, {"AP_begin_indices", "AP_end_indices"});

  vector<double> apduration;
  int retval =
      __AP_duration(doubleFeatures.at("T"), intFeatures.at("AP_begin_indices"),
                    intFeatures.at("AP_end_indices"), apduration);

  if (retval > 0) {
    setVec(DoubleFeatureData, StringData, "AP_duration", apduration);
  }
  return retval;
}

// *** AP_duration_change according to E23 ***
static int __AP_duration_change(const vector<double>& apduration,
                                vector<double>& apdurationchange) {
  if (apduration.size() < 1) {
    return -1;
  }
  apdurationchange.resize(apduration.size() - 1);
  for (size_t i = 0; i < apdurationchange.size(); i++) {
    apdurationchange[i] = (apduration[i + 1] - apduration[0]) / apduration[0];
  }
  return apdurationchange.size();
}
int SpikeShape::AP_duration_change(mapStr2intVec& IntFeatureData,
                              mapStr2doubleVec& DoubleFeatureData,
                              mapStr2Str& StringData) {
  const auto& features = getFeatures(DoubleFeatureData, {"AP_duration"});
  const vector<double>& apduration = features.at("AP_duration");

  vector<double> apdurationchange;
  int retval = __AP_duration_change(apduration, apdurationchange);
  if (retval > 0) {
    setVec(DoubleFeatureData, StringData, "AP_duration_change",
           apdurationchange);
  }
  return retval;
}

// AP_width_between_threshold
// spike width calculation according to threshold value.
// min_between_peaks_indices are used to split the different APs
// do not add width if threshold crossing has not been found
static int __AP_width_between_threshold(
    const vector<double>& t, const vector<double>& v, double stimstart,
    double threshold, const vector<int>& min_between_peaks_indices,
    vector<double>& ap_width_threshold) {
  vector<int> indices(min_between_peaks_indices.size() + 1);
  int start_index = distance(
      t.begin(), find_if(t.begin(), t.end(),
                         [stimstart](double x) { return x >= stimstart; }));
  indices[0] = start_index;
  copy(min_between_peaks_indices.begin(), min_between_peaks_indices.end(),
       indices.begin() + 1);

  for (size_t i = 0; i < indices.size() - 1; i++) {
    int onset_index = distance(
        v.begin(), find_if(v.begin() + indices[i], v.begin() + indices[i + 1],
                           [threshold](double x) { return x >= threshold; }));
    int end_index = distance(
        v.begin(), find_if(v.begin() + onset_index, v.begin() + indices[i + 1],
                           [threshold](double x) { return x <= threshold; }));
    if (end_index != indices[i + 1]) {
      ap_width_threshold.push_back(t[end_index] - t[onset_index]);
    }
  }

  return ap_width_threshold.size();
}

int SpikeShape::AP_width_between_threshold(mapStr2intVec& IntFeatureData,
                                      mapStr2doubleVec& DoubleFeatureData,
                                      mapStr2Str& StringData) {
  // Retrieve all required double and int features at once
  const auto& doubleFeatures =
      getFeatures(DoubleFeatureData, {"T", "V", "Threshold", "stim_start"});
  const auto& intFeatures =
      getFeatures(IntFeatureData, {"min_between_peaks_indices"});

  vector<double> ap_width_threshold;
  int retval = __AP_width_between_threshold(
      doubleFeatures.at("T"), doubleFeatures.at("V"),
      doubleFeatures.at("stim_start").front(),
      doubleFeatures.at("Threshold").front(),
      intFeatures.at("min_between_peaks_indices"), ap_width_threshold);
  if (retval == 0)
    return -1;
  else if (retval > 0) {
    setVec(DoubleFeatureData, StringData, "AP_width_between_threshold",
           ap_width_threshold);
  }
  return retval;
}

// spike half width
// for spike amplitude = v_peak - v_AHP
static int __spike_width1(const vector<double>& t, const vector<double>& v,
                          const vector<int>& peak_indices,
                          const vector<int>& min_ahp_indices, double stim_start,
                          vector<double>& spike_width1) {
  int start_index = distance(
      t.begin(), find_if(t.begin(), t.end(), [stim_start](double t_val) {
        return t_val >= stim_start;
      }));
  vector<int> min_ahp_indices_plus(min_ahp_indices.size() + 1, start_index);
  copy(min_ahp_indices.begin(), min_ahp_indices.end(),
       min_ahp_indices_plus.begin() + 1);
  for (size_t i = 1; i < min_ahp_indices_plus.size(); i++) {
    double v_half = (v[peak_indices[i - 1]] + v[min_ahp_indices_plus[i]]) / 2.;
    // interpolate this one time step where the voltage is close to v_half in
    // the rising and in the falling edge
    double v_dev;
    double delta_v;
    double t_dev_rise;
    double t_dev_fall;
    double delta_t;
    int rise_index = distance(
        v.begin(), find_if(v.begin() + min_ahp_indices_plus[i - 1],
                           v.begin() + peak_indices[i - 1],
                           [v_half](double v_val) { return v_val >= v_half; }));
    v_dev = v_half - v[rise_index];
    delta_v = v[rise_index] - v[rise_index - 1];
    delta_t = t[rise_index] - t[rise_index - 1];
    t_dev_rise = delta_t * v_dev / delta_v;
    int fall_index = distance(
        v.begin(), find_if(v.begin() + peak_indices[i - 1],
                           v.begin() + min_ahp_indices_plus[i],
                           [v_half](double v_val) { return v_val <= v_half; }));
    v_dev = v_half - v[fall_index];
    delta_v = v[fall_index] - v[fall_index - 1];
    delta_t = t[fall_index] - t[fall_index - 1];
    t_dev_fall = delta_t * v_dev / delta_v;
    spike_width1.push_back(t[fall_index] - t_dev_rise - t[rise_index] +
                           t_dev_fall);
  }
  return spike_width1.size();
}

int SpikeShape::spike_width1(mapStr2intVec& IntFeatureData,
                        mapStr2doubleVec& DoubleFeatureData,
                        mapStr2Str& StringData) {
  const auto& doubleFeatures =
      getFeatures(DoubleFeatureData, {"V", "T", "stim_start"});
  const auto& intFeatures =
      getFeatures(IntFeatureData, {"min_AHP_indices", "peak_indices"});

  vector<double> spike_width1;
  // Calculate spike width
  int retVal = __spike_width1(doubleFeatures.at("T"), doubleFeatures.at("V"),
                              intFeatures.at("peak_indices"),
                              intFeatures.at("min_AHP_indices"),
                              doubleFeatures.at("stim_start")[0], spike_width1);

  if (retVal >= 0) {
    setVec(DoubleFeatureData, StringData, "spike_half_width", spike_width1);
  }
  return retVal;
}

// Width of the first spike
int SpikeShape::AP1_width(mapStr2intVec& IntFeatureData,
                     mapStr2doubleVec& DoubleFeatureData,
                     mapStr2Str& StringData) {
  const vector<double>& spike_half_width =
      getFeature(DoubleFeatureData, "spike_half_width");
  vector<double> AP1_width;
  AP1_width.push_back(spike_half_width[0]);
  setVec(DoubleFeatureData, StringData, "AP1_width", AP1_width);
  return 1;
}

// Width of the second spike
int SpikeShape::AP2_width(mapStr2intVec& IntFeatureData,
                     mapStr2doubleVec& DoubleFeatureData,
                     mapStr2Str& StringData) {
  const vector<double>& spike_half_width =
      getFeature(DoubleFeatureData, "spike_half_width");
  vector<double> AP2_width;
  if (spike_half_width.size() < 2) {
    throw FeatureComputationError(
        "Size of spike_half_width should be >= 2 for AP2_width");
  }
  AP2_width.push_back(spike_half_width[1]);
  setVec(DoubleFeatureData, StringData, "AP2_width", AP2_width);
  return 1;
}

// Width of the last spike
int SpikeShape::APlast_width(mapStr2intVec& IntFeatureData,
                        mapStr2doubleVec& DoubleFeatureData,
                        mapStr2Str& StringData) {
  const vector<double>& spike_half_width =
      getFeature(DoubleFeatureData, "spike_half_width");
  vector<double> APlast_width;
  APlast_width.push_back(spike_half_width[spike_half_width.size() - 1]);
  setVec(DoubleFeatureData, StringData, "APlast_width", APlast_width);
  return 1;
}

// To find spike width using Central difference derivative vec1[i] =
// ((vec[i+1]+vec[i-1])/2)/dx  and half width is between
// MinAHP and APThreshold
static int __spike_width2(const vector<double>& t, const vector<double>& V,
                          const vector<int>& PeakIndex,
                          const vector<int>& minAHPIndex,
                          vector<double>& spike_width2) {
  vector<double> v, dv1, dv2;
  double dx = t[1] - t[0];
  double VoltThreshold, VoltMax, HalfV, T0, V0, V1, fraction, TStart, TEnd;
  for (size_t i = 0; i < minAHPIndex.size() && i < PeakIndex.size() - 1; i++) {
    v.clear();
    dv1.clear();
    dv2.clear();

    for (int j = minAHPIndex[i]; j <= PeakIndex[i + 1]; j++) {
      if (j < 0) {
        throw FeatureComputationError("Invalid index");
      }
      v.push_back(V[j]);
    }

    getCentralDifferenceDerivative(dx, v, dv1);
    getCentralDifferenceDerivative(dx, dv1, dv2);
    double dMax = dv2[0];

    size_t index = 0;
    for (size_t j = 1; j < dv2.size(); ++j) {
      if (dMax <= dv2[j]) {
        dMax = dv2[j];
        index = j;
      }
    }

    // Take voltage at point where double derivative is maximum
    index += minAHPIndex[i];
    VoltThreshold = V[index];
    VoltMax = V[PeakIndex[i + 1]];
    HalfV = (VoltMax + VoltThreshold) / 2;

    // Find voltage where it crosses HalfV in the rising phase of action
    // potential
    for (size_t j = 0; j < v.size(); ++j) {
      if (v[j] > HalfV) {  // point is found where  v is crossing HalfV
        index = minAHPIndex[i] + j;
        break;
      }
    }

    // index is the index where v crossed HalfV now use linear interpolation to
    // find actual time at HalfV
    T0 = t[index - 1];
    V0 = V[index - 1];
    V1 = V[index];
    fraction = (HalfV - V0) / (V1 - V0);
    TStart = T0 + (fraction * dx);

    // Find voltage where it crosses HalfV in the falling phase of the action
    // potential
    for (size_t j = PeakIndex[i + 1]; j < V.size(); j++) {
      if (V[j] < HalfV) {
        index = j;
        break;
      }
    }

    if (index == V.size()) {
      throw FeatureComputationError("Falling phase of last spike is missing.");
    }

    T0 = t[index - 1];
    V0 = V[index - 1];
    V1 = V[index];
    fraction = (HalfV - V0) / (V1 - V0);
    TEnd = T0 + (fraction * dx);

    spike_width2.push_back(TEnd - TStart);
  }

  return spike_width2.size();
}

int SpikeShape::spike_width2(mapStr2intVec& IntFeatureData,
                        mapStr2doubleVec& DoubleFeatureData,
                        mapStr2Str& StringData) {
  const auto& doubleFeatures = getFeatures(DoubleFeatureData, {"V", "T"});
  const auto& intFeatures =
      getFeatures(IntFeatureData, {"min_AHP_indices", "peak_indices"});
  if (intFeatures.at("peak_indices").size() <= 1) {
    throw FeatureComputationError("More than one spike is needed for spikewidth2 calculation.");
  }
  vector<double> spike_width2;
  int retVal = __spike_width2(doubleFeatures.at("T"), doubleFeatures.at("V"),
                              intFeatures.at("peak_indices"),
                              intFeatures.at("min_AHP_indices"), spike_width2);
  if (retVal > 0) {
    setVec(DoubleFeatureData, StringData, "spike_width2", spike_width2);
  }
  return retVal;
}

// spike width at spike start
static int __AP_begin_width(const vector<double>& t, const vector<double>& v,
                            double stimstart,
                            const vector<int>& AP_begin_indices,
                            const vector<int>& min_ahp_indices,
                            vector<double>& AP_begin_width) {
  // int start_index = distance(t.begin(), find_if(t.begin(), t.end(),
  // std::bind2nd(std::greater_equal<double>(), stim_start)));
  /// vector<int> min_ahp_indices_plus(min_ahp_indices.size() + 1, start_index);
  // copy(min_ahp_indices.begin(), min_ahp_indices.end(),
  // min_ahp_indices_plus.begin());

  // keep only min_ahp_indices values that are after stim start
  // because AP_begin_indices are all after stim start
  // if not done, could cause cases where AP_begin_indices[i] > min_ahp_indices[i]
  // leading to segmentation faults
  auto it = find_if(t.begin(), t.end(),
                    [stimstart](double val) { return val >= stimstart; });
  int stimbeginindex = distance(t.begin(), it);
  vector<int> strict_min_ahp_indices;
  for (size_t i = 0; i < min_ahp_indices.size(); i++) {
    if (min_ahp_indices[i] > stimbeginindex) {
      strict_min_ahp_indices.push_back(min_ahp_indices[i]);
    }
  }

  if (AP_begin_indices.size() < strict_min_ahp_indices.size()) return -1;
  for (size_t i = 0; i < strict_min_ahp_indices.size(); i++) {
    double v_start = v[AP_begin_indices[i]];
    // interpolate this one time step where the voltage is close to v_start in
    // the falling edge
    int rise_index = AP_begin_indices[i];
    int fall_index = distance(
        v.begin(),
        find_if(v.begin() + rise_index + 1,
                v.begin() + strict_min_ahp_indices[i],
                [v_start](const auto& val) { return val <= v_start; }));
    // v_dev = v_start - v[fall_index];
    // delta_v = v[fall_index] - v[fall_index - 1];
    // delta_t = t[fall_index] - t[fall_index - 1];
    // t_dev_fall = delta_t * v_dev / delta_v;
    // printf("%f %f\n",delta_v, t_dev_fall);
    AP_begin_width.push_back(t[fall_index] - t[rise_index] /*+ t_dev_fall*/);
  }
  return AP_begin_width.size();
}

int SpikeShape::AP_begin_width(mapStr2intVec& IntFeatureData,
                          mapStr2doubleVec& DoubleFeatureData,
                          mapStr2Str& StringData) {
  const auto& doubleFeatures =
      getFeatures(DoubleFeatureData, {"V", "T", "stim_start"});
  const auto& intFeatures =
      getFeatures(IntFeatureData, {"min_AHP_indices", "AP_begin_indices"});
  vector<double> AP_begin_width;
  const vector<double>& V = doubleFeatures.at("V");
  const vector<double>& t = doubleFeatures.at("T");
  const vector<double>& stim_start = doubleFeatures.at("stim_start");
  const vector<int>& min_AHP_indices = intFeatures.at("min_AHP_indices");
  const vector<int>& AP_begin_indices = intFeatures.at("AP_begin_indices");
  int retVal = __AP_begin_width(t, V, stim_start[0], AP_begin_indices,
                                min_AHP_indices, AP_begin_width);
  if (retVal > 0) {
    setVec(DoubleFeatureData, StringData, "AP_begin_width", AP_begin_width);
  }
  return retVal;
}

int SpikeShape::AP1_begin_width(mapStr2intVec& IntFeatureData,
                           mapStr2doubleVec& DoubleFeatureData,
                           mapStr2Str& StringData) {
  const vector<double>& AP_begin_width =
      getFeature(DoubleFeatureData, "AP_begin_width");
  vector<double> AP1_begin_width;
  AP1_begin_width.push_back(AP_begin_width[0]);
  setVec(DoubleFeatureData, StringData, "AP1_begin_width", AP1_begin_width);
  return 1;
}

int SpikeShape::AP2_begin_width(mapStr2intVec& IntFeatureData,
                           mapStr2doubleVec& DoubleFeatureData,
                           mapStr2Str& StringData) {
  const vector<double>& AP_begin_width =
      getFeature(DoubleFeatureData, "AP_begin_width");
  vector<double> AP2_begin_width;
  if (AP_begin_width.size() < 2) {
    throw FeatureComputationError("There are less than 2 spikes in the trace.");
  }
  AP2_begin_width.push_back(AP_begin_width[1]);
  setVec(DoubleFeatureData, StringData, "AP2_begin_width", AP2_begin_width);
  return 1;
}

// Difference amplitude of the second to first spike
int SpikeShape::AP2_AP1_begin_width_diff(mapStr2intVec& IntFeatureData,
                                    mapStr2doubleVec& DoubleFeatureData,
                                    mapStr2Str& StringData) {
  const vector<double>& AP_begin_widths =
      getFeature(DoubleFeatureData, "AP_begin_width");
  vector<double> AP2_AP1_begin_width_diff;
  if (AP_begin_widths.size() < 2) {
    throw FeatureComputationError("There are less than 2 spikes in the trace.");
  }
  AP2_AP1_begin_width_diff.push_back(AP_begin_widths[1] - AP_begin_widths[0]);
  setVec(DoubleFeatureData, StringData, "AP2_AP1_begin_width_diff",
         AP2_AP1_begin_width_diff);
  return 1;
}

//
// *** AP begin indices ***
//
static int __AP_begin_indices(const vector<double>& t, const vector<double>& v,
                              double stimstart, double stimend,
                              const vector<int>& pi, const vector<int>& ahpi,
                              vector<int>& apbi, double dTh,
                              int derivative_window) {
  const double derivativethreshold = dTh;
  vector<double> dvdt(v.size());
  vector<double> dv;
  vector<double> dt;
  getCentralDifferenceDerivative(1., v, dv);
  getCentralDifferenceDerivative(1., t, dt);
  transform(dv.begin(), dv.end(), dt.begin(), dvdt.begin(),
            std::divides<double>());

  // restrict to time interval where stimulus is applied
  vector<int> minima, peak_indices;
  auto stimbeginindex =
      distance(t.begin(), find_if(t.begin(), t.end(), [stimstart](double time) {
                 return time >= stimstart;
               }));

  if (stimbeginindex > 1) {
    // to avoid skipping AP_begin when it is exactly at stimbeginindex
    // also because of float precision and interpolation, sometimes
    // stimbeginindex can be slightly above 'real' stim begin index
    minima.push_back(stimbeginindex - 2);
  } else if (stimbeginindex == 1) {
    minima.push_back(stimbeginindex - 1);
  } else {
    minima.push_back(stimbeginindex);
  }
  for (size_t i = 0; i < ahpi.size(); i++) {
    if (ahpi[i] > stimbeginindex) {
      minima.push_back(ahpi[i]);
    }
    // if(t[ahpi[i]] > stimend) {
    //    break;
    //}
  }
  for (size_t i = 0; i < pi.size(); i++) {
    if (pi[i] > stimbeginindex) {
      peak_indices.push_back(pi[i]);
    }
  }
  int endindex = distance(t.begin(), t.end());
  if (minima.size() < peak_indices.size())
    throw FeatureComputationError("More peaks than min_AHP in AP_begin_indices.");

  // printf("Found %d minima\n", minima.size());
  for (size_t i = 0; i < peak_indices.size(); i++) {
    // assure that the width of the slope is bigger than 4
    int newbegin = peak_indices[i];
    int begin = minima[i];
    int width = derivative_window;
    bool skip = false;

    // Detect where the derivate crosses derivativethreshold, and make sure
    // this happens in a window of 'width' sampling point
    do {
      begin = distance(dvdt.begin(),
                       find_if(
                           // use reverse iterator to get last occurence
                           // and avoid false positive long before the spike
                           dvdt.rbegin() + v.size() - newbegin,
                           dvdt.rbegin() + v.size() - minima[i],
                           [derivativethreshold](double val) {
                             return val <= derivativethreshold;
                           })
                           .base());
      // cover edge case to avoid going out of index in the while condition
      if (begin > endindex - width) {
        begin = endindex - width;
      }
      // printf("%d %d\n", newbegin, minima[i+1]);
      if (begin == minima[i]) {
        // printf("Skipping %d %d\n", newbegin, minima[i+1]);
        // could not find a spike in between these minima
        skip = true;
        break;
      }
      newbegin = begin - 1;
    } while (find_if(dvdt.begin() + begin, dvdt.begin() + begin + width,
                     [derivativethreshold](double val) {
                       return val < derivativethreshold;
                     }) != dvdt.begin() + begin + width);
    if (skip) {
      continue;
    }
    apbi.push_back(begin);
  }
  return apbi.size();
}

int SpikeShape::AP_begin_indices(mapStr2intVec& IntFeatureData,
                            mapStr2doubleVec& DoubleFeatureData,
                            mapStr2Str& StringData) {
  vector<int> min_ahp_idx;
  // Retrieve all required double and int features at once
  const auto& doubleFeatures =
      getFeatures(DoubleFeatureData, {"T", "V", "stim_start", "stim_end"});
  const auto& intFeatures =
      getFeatures(IntFeatureData, {"peak_indices"});
  try{
    const auto& intFeaturesSpecial = getFeatures(IntFeatureData,
                                                 {"min_AHP_indices"});
    min_ahp_idx = intFeaturesSpecial.at("min_AHP_indices");
  // Edge case for 1 spike and no min_AHP_indices:
  // continue with empty vector.
  } catch (const EmptyFeatureError& e) {}
    catch (const FeatureComputationError& e) {}

  vector<int> apbi;

  // Get DerivativeThreshold
  vector<double> dTh;
  int retVal = getParam(DoubleFeatureData, "DerivativeThreshold", dTh);
  if (retVal <= 0) {
    // derivative at peak start according to eCode specification 10mV/ms
    // according to Shaul 12mV/ms
    dTh.push_back(12.0);
  }

  // Get DerivativeWindow
  vector<int> derivative_window;
  retVal = getParam(IntFeatureData, "DerivativeWindow", derivative_window);
  if (retVal <= 0)
    throw FeatureComputationError("DerivativeWindow not set.");

  // Calculate feature
  retVal = __AP_begin_indices(
      doubleFeatures.at("T"), doubleFeatures.at("V"),
      doubleFeatures.at("stim_start")[0], doubleFeatures.at("stim_end")[0],
      intFeatures.at("peak_indices"), min_ahp_idx, apbi, dTh[0],
      derivative_window[0]);

  // Save feature value
  if (retVal >= 0) {
    setVec(IntFeatureData, StringData, "AP_begin_indices", apbi);
  }
  return retVal;
}

static int __AP_end_indices(const vector<double>& t, const vector<double>& v,
                            const double stimstart, const vector<int>& pi,
                            vector<int>& apei, double derivativethreshold) {

  vector<double> dvdt(v.size());
  vector<double> dv;
  vector<double> dt;
  vector<int> peak_indices;
  int max_slope;
  getCentralDifferenceDerivative(1., v, dv);
  getCentralDifferenceDerivative(1., t, dt);
  transform(dv.begin(), dv.end(), dt.begin(), dvdt.begin(),
            std::divides<double>());

  auto stimbeginindex = distance(t.begin(),
      find_if(t.begin(), t.end(),
          [stimstart](double time){ return time >= stimstart; }));

   for (size_t i = 0; i < pi.size(); i++) {
    if (pi[i] > stimbeginindex) {
      peak_indices.push_back(pi[i]);
    }
  }
  peak_indices.push_back(v.size() - 1);

  for (size_t i = 0; i < peak_indices.size() - 1; i++) {
    size_t start_index = peak_indices[i] + 1;
    size_t end_index = peak_indices[i + 1];

    if (start_index >= end_index || start_index >= dvdt.size() || end_index >= dvdt.size()) {
      continue;
    }

    auto min_element_it = std::min_element(dvdt.begin() + start_index, dvdt.begin() + end_index);
    auto max_slope = std::distance(dvdt.begin(), min_element_it);
    // assure that the width of the slope is bigger than 4
    auto threshold_it = std::find_if(dvdt.begin() + max_slope, dvdt.begin() + end_index,
                                    [derivativethreshold](double x) { return x >= derivativethreshold; });

    if (threshold_it != dvdt.begin() + end_index) {
      apei.push_back(std::distance(dvdt.begin(), threshold_it));
    }
  }
  return apei.size();
}

int SpikeShape::AP_end_indices(mapStr2intVec& IntFeatureData,
                          mapStr2doubleVec& DoubleFeatureData,
                          mapStr2Str& StringData) {
  const auto& T = getFeature(DoubleFeatureData, "T");
  const auto& V = getFeature(DoubleFeatureData, "V");
  const auto& stim_start = getFeature(DoubleFeatureData, "stim_start");
  const auto& peak_indices = getFeature(IntFeatureData, "peak_indices");

  vector<double> dTh;
  int retVal = getParam(DoubleFeatureData, "DownDerivativeThreshold", dTh);
  double downDerivativeThreshold = (retVal <= 0) ? -12.0 : dTh[0];

  vector<int> AP_end_indices;
  retVal = __AP_end_indices(T, V, stim_start[0], peak_indices, AP_end_indices,
                            downDerivativeThreshold);
  if (retVal >= 0) {
    setVec(IntFeatureData, StringData, "AP_end_indices", AP_end_indices);
  }
  return retVal;
}

static int __AP_begin_voltage(const vector<double>& t, const vector<double>& v,
                              const vector<int>& AP_begin_indices,
                              vector<double>& AP_begin_voltage) {
  for (size_t i = 0; i < AP_begin_indices.size(); i++) {
    AP_begin_voltage.push_back(v[AP_begin_indices[i]]);
  }
  return AP_begin_voltage.size();
}

int SpikeShape::AP_begin_voltage(mapStr2intVec& IntFeatureData,
                            mapStr2doubleVec& DoubleFeatureData,
                            mapStr2Str& StringData) {
  const auto& doubleFeatures = getFeatures(DoubleFeatureData, {"V", "T"});
  const auto& intFeatures = getFeatures(IntFeatureData, {"AP_begin_indices"});
  vector<double> AP_begin_voltage;
  const vector<double>& V = doubleFeatures.at("V");
  const vector<double>& t = doubleFeatures.at("T");
  const vector<int>& AP_begin_indices = intFeatures.at("AP_begin_indices");
  int retVal = __AP_begin_voltage(t, V, AP_begin_indices, AP_begin_voltage);
  if (retVal > 0) {
    setVec(DoubleFeatureData, StringData, "AP_begin_voltage", AP_begin_voltage);
  }
  return retVal;
}

int SpikeShape::AP1_begin_voltage(mapStr2intVec& IntFeatureData,
                             mapStr2doubleVec& DoubleFeatureData,
                             mapStr2Str& StringData) {
  const vector<double>& AP_begin_voltage =
      getFeature(DoubleFeatureData, "AP_begin_voltage");
  vector<double> AP1_begin_voltage;
  AP1_begin_voltage.push_back(AP_begin_voltage[0]);
  setVec(DoubleFeatureData, StringData, "AP1_begin_voltage", AP1_begin_voltage);
  return 1;
}

int SpikeShape::AP2_begin_voltage(mapStr2intVec& IntFeatureData,
                             mapStr2doubleVec& DoubleFeatureData,
                             mapStr2Str& StringData) {
  const vector<double>& AP_begin_voltage =
      getFeature(DoubleFeatureData, "AP_begin_voltage");
  vector<double> AP2_begin_voltage;

  if (AP_begin_voltage.size() < 2) {
    throw FeatureComputationError("There are less than 2 spikes in the trace.");
  }
  AP2_begin_voltage.push_back(AP_begin_voltage[1]);
  setVec(DoubleFeatureData, StringData, "AP2_begin_voltage", AP2_begin_voltage);
  return 1;
}

static int __AP_begin_time(const vector<double>& t, const vector<double>& v,
                           const vector<int>& AP_begin_indices,
                           vector<double>& AP_begin_time) {
  for (size_t i = 0; i < AP_begin_indices.size(); i++) {
    AP_begin_time.push_back(t[AP_begin_indices[i]]);
  }
  return AP_begin_time.size();
}

int SpikeShape::AP_begin_time(mapStr2intVec& IntFeatureData,
                         mapStr2doubleVec& DoubleFeatureData,
                         mapStr2Str& StringData) {
  const auto& doubleFeatures = getFeatures(DoubleFeatureData, {"V", "T"});
  const auto& intFeatures = getFeatures(IntFeatureData, {"AP_begin_indices"});
  vector<double> AP_begin_time;
  const vector<double>& V = doubleFeatures.at("V");
  const vector<double>& t = doubleFeatures.at("T");
  const vector<int>& AP_begin_indices = intFeatures.at("AP_begin_indices");
  int retVal = __AP_begin_time(t, V, AP_begin_indices, AP_begin_time);
  if (retVal > 0) {
    setVec(DoubleFeatureData, StringData, "AP_begin_time", AP_begin_time);
  }
  return retVal;
}

//
// *** Action potential peak upstroke ***
//
static int __AP_peak_upstroke(const vector<double>& t, const vector<double>& v,
                              const vector<int>& pi,    // peak indices
                              const vector<int>& apbi,  // AP begin indices
                              vector<double>& pus) {    // AP peak upstroke
  vector<double> dvdt(v.size());
  vector<double> dv;
  vector<double> dt;
  getCentralDifferenceDerivative(1., v, dv);
  getCentralDifferenceDerivative(1., t, dt);
  transform(dv.begin(), dv.end(), dt.begin(), dvdt.begin(),
            std::divides<double>());

  // Make sure that each value of pi is greater than its apbi counterpart
  vector<int> new_pi;
  size_t j = 0;
  for (size_t i = 0; i < apbi.size(); i++) {
    while (j < pi.size() && pi[j] < apbi[i]) {
      j++;
    }

    if (j < pi.size() && pi[j] >= apbi[i]) {
      new_pi.push_back(pi[j]);
      j++;
    }
  }

  for (size_t i = 0; i < std::min(apbi.size(), new_pi.size()); i++) {
    pus.push_back(
        *std::max_element(dvdt.begin() + apbi[i], dvdt.begin() + new_pi[i]));
  }

  return pus.size();
}

int SpikeShape::AP_peak_upstroke(mapStr2intVec& IntFeatureData,
                            mapStr2doubleVec& DoubleFeatureData,
                            mapStr2Str& StringData) {
  // Retrieve all required double and int features at once
  const auto& doubleFeatures = getFeatures(DoubleFeatureData, {"T", "V"});
  const auto& intFeatures =
      getFeatures(IntFeatureData, {"AP_begin_indices", "peak_indices"});

  vector<double> pus;
  int retVal = __AP_peak_upstroke(
      doubleFeatures.at("T"), doubleFeatures.at("V"),
      intFeatures.at("peak_indices"), intFeatures.at("AP_begin_indices"), pus);

  if (retVal >= 0) {
    setVec(DoubleFeatureData, StringData, "AP_peak_upstroke", pus);
  }
  return retVal;
}

//
// *** Action potential peak downstroke ***
//
static int __AP_peak_downstroke(const vector<double>& t,
                                const vector<double>& v,
                                const vector<int>& pi,    // peak indices
                                const vector<int>& ahpi,  // min AHP indices
                                vector<double>& pds) {    // AP peak downstroke
  vector<double> dvdt(v.size());
  vector<double> dv;
  vector<double> dt;
  getCentralDifferenceDerivative(1., v, dv);
  getCentralDifferenceDerivative(1., t, dt);
  transform(dv.begin(), dv.end(), dt.begin(), dvdt.begin(),
            std::divides<double>());

  for (size_t i = 0; i < std::min(ahpi.size(), pi.size()); i++) {
    pds.push_back(
        *std::min_element(dvdt.begin() + pi[i], dvdt.begin() + ahpi[i]));
  }

  return pds.size();
}

int SpikeShape::AP_peak_downstroke(mapStr2intVec& IntFeatureData,
                              mapStr2doubleVec& DoubleFeatureData,
                              mapStr2Str& StringData) {
  const auto& doubleFeatures = getFeatures(DoubleFeatureData, {"T", "V"});
  const auto& intFeatures =
      getFeatures(IntFeatureData, {"min_AHP_indices", "peak_indices"});

  vector<double> pds;
  int retVal = __AP_peak_downstroke(
      doubleFeatures.at("T"), doubleFeatures.at("V"),
      intFeatures.at("peak_indices"), intFeatures.at("min_AHP_indices"), pds);

  if (retVal >= 0) {
    setVec(DoubleFeatureData, StringData, "AP_peak_downstroke", pds);
  }
  return retVal;
}

// *** AP rise indices ***
//
static int __AP_rise_indices(const vector<double>& v, const vector<int>& apbi,
                             const vector<int>& pi, vector<int>& apri) {
  apri.resize(std::min(apbi.size(), pi.size()));
  for (size_t i = 0; i < apri.size(); i++) {
    double halfheight = (v[pi[i]] + v[apbi[i]]) / 2.;
    vector<double> vpeak;
    if (pi[i] < apbi[i]) {
      // For some reason the peak and begin indices are out of sync
      // Peak should always be later than begin index
      return -1;
    }
    vpeak.resize(pi[i] - apbi[i]);
    transform(v.begin() + apbi[i], v.begin() + pi[i], vpeak.begin(),
              [halfheight](double val) { return fabs(val - halfheight); });
    apri[i] = distance(vpeak.begin(), min_element(vpeak.begin(), vpeak.end())) +
              apbi[i];
  }
  return apri.size();
}
int SpikeShape::AP_rise_indices(mapStr2intVec& IntFeatureData,
                           mapStr2doubleVec& DoubleFeatureData,
                           mapStr2Str& StringData) {
  const auto& doubleFeatures = getFeatures(DoubleFeatureData, {"V"});
  const auto& intFeatures =
      getFeatures(IntFeatureData, {"AP_begin_indices", "peak_indices"});
  vector<int> apriseindices;
  int retval = __AP_rise_indices(doubleFeatures.at("V"),
                                 intFeatures.at("AP_begin_indices"),
                                 intFeatures.at("peak_indices"), apriseindices);
  if (retval > 0) {
    setVec(IntFeatureData, StringData, "AP_rise_indices", apriseindices);
  }
  return retval;
}

// *** AP fall indices ***
//
static int __AP_fall_indices(const vector<double>& v, const vector<int>& apbi,
                             const vector<int>& apei, const vector<int>& pi,
                             vector<int>& apfi) {
  apfi.resize(std::min(apbi.size(), pi.size()));
  for (size_t i = 0; i < apfi.size(); i++) {
    if (pi[i] >= v.size() || apbi[i] >= v.size() || apei[i] >= v.size() || pi[i] > apei[i]) {
        continue;
    }
    double halfheight = (v[pi[i]] + v[apbi[i]]) / 2.;
    vector<double> vpeak(&v[pi[i]], &v[apei[i]]);
    transform(vpeak.begin(), vpeak.end(), vpeak.begin(),
              [halfheight](double val) { return fabs(val - halfheight); });
    apfi[i] = distance(vpeak.begin(), min_element(vpeak.begin(), vpeak.end())) +
              pi[i];
  }
  return apfi.size();
}
int SpikeShape::AP_fall_indices(mapStr2intVec& IntFeatureData,
                           mapStr2doubleVec& DoubleFeatureData,
                           mapStr2Str& StringData) {
  const auto& doubleFeatures = getFeatures(DoubleFeatureData, {"V"});
  const auto& intFeatures = getFeatures(
      IntFeatureData, {"AP_begin_indices", "AP_end_indices", "peak_indices"});
  vector<int> apfallindices;
  int retval = __AP_fall_indices(doubleFeatures.at("V"),
                                 intFeatures.at("AP_begin_indices"),
                                 intFeatures.at("AP_end_indices"),
                                 intFeatures.at("peak_indices"), apfallindices);
  if (retval > 0) {
    setVec(IntFeatureData, StringData, "AP_fall_indices", apfallindices);
  }
  return retval;
}

// *** AP_rise_time according to E9 and E17 ***
static int __AP_rise_time(const vector<double>& t, const vector<double>& v,
                          const vector<int>& apbeginindices,
                          const vector<int>& peakindices,
                          double beginperc, double endperc,
                          vector<double>& aprisetime) {
  aprisetime.resize(std::min(apbeginindices.size(), peakindices.size()));
  // Make sure that we do not use peaks starting before the 1st AP_begin_index
  // Because AP_begin_indices only takes into account peaks after stimstart
  vector<int> newpeakindices;
  if (apbeginindices.size() > 0) {
    newpeakindices = peaks_after_stim_start(apbeginindices[0], peakindices);
  }
  double begin_v;
  double end_v;
  double begin_indice;
  double end_indice;
  double apamplitude;
  for (size_t i = 0; i < aprisetime.size(); i++) {
    // do not use AP_amplitude feature because it does not take into account
    // peaks after stim_end
    apamplitude = v[newpeakindices[i]] - v[apbeginindices[i]];
    begin_v = v[apbeginindices[i]] + beginperc * apamplitude;
    end_v = v[apbeginindices[i]] + endperc * apamplitude;

    // Get begin indice
    size_t j = apbeginindices[i];
    // change slightly begin_v for almost equal case
    // truncature error can change begin_v even when beginperc == 0.0
    while (j < newpeakindices[i] && v[j] < begin_v - 0.0000000000001){
      j++;
    }
    begin_indice = j;

    // Get end indice
    j = newpeakindices[i];
    // change slightly end_v for almost equal case
    // truncature error can change end_v even when beginperc == 0.0
    while (j > apbeginindices[i] && v[j] > end_v + 0.0000000000001) {
      j--;
    }
    end_indice = j;

    aprisetime[i] = t[end_indice] - t[begin_indice];
  }
  return aprisetime.size();
}
int SpikeShape::AP_rise_time(mapStr2intVec& IntFeatureData,
                        mapStr2doubleVec& DoubleFeatureData,
                        mapStr2Str& StringData) {
  // Fetching all required features in one go.
  const auto& doubleFeatures = getFeatures(
      DoubleFeatureData,
      {"T", "V", "rise_start_perc", "rise_end_perc"});
  const auto& intFeatures =
      getFeatures(IntFeatureData, {"AP_begin_indices", "peak_indices"});
  vector<double> aprisetime;
  int retval = __AP_rise_time(
      doubleFeatures.at("T"), doubleFeatures.at("V"),
      intFeatures.at("AP_begin_indices"), intFeatures.at("peak_indices"),
      doubleFeatures.at("rise_start_perc").empty()
          ? 0.0
          : doubleFeatures.at("rise_start_perc").front(),
      doubleFeatures.at("rise_end_perc").empty()
          ? 1.0
          : doubleFeatures.at("rise_end_perc").front(),
      aprisetime);
  if (retval > 0) {
    setVec(DoubleFeatureData, StringData, "AP_rise_time", aprisetime);
  }
  return retval;
}

// *** AP_fall_time according to E10 and E18 ***
static int __AP_fall_time(const vector<double>& t,
                          const double stimstart,
                          const vector<int>& peakindices,
                          const vector<int>& apendindices,
                          vector<double>& apfalltime) {
  apfalltime.resize(std::min(peakindices.size(), apendindices.size()));
  // Make sure that we do not use peaks starting before stim start
  // Because AP_end_indices only takes into account peaks after stim start
  vector<int> newpeakindices = peaks_after_stim_start(stimstart, peakindices, t);

  for (size_t i = 0; i < apfalltime.size(); i++) {
    apfalltime[i] = t[apendindices[i]] - t[newpeakindices[i]];
  }
  return apfalltime.size();
}
int SpikeShape::AP_fall_time(mapStr2intVec& IntFeatureData,
                        mapStr2doubleVec& DoubleFeatureData,
                        mapStr2Str& StringData) {
  const auto& doubleFeatures = getFeatures(DoubleFeatureData, {"T", "stim_start"});
  const auto& intFeatures =
      getFeatures(IntFeatureData, {"peak_indices", "AP_end_indices"});

  const vector<double>& t = doubleFeatures.at("T");
  const double stim_start = doubleFeatures.at("stim_start")[0];
  const vector<int>& peakindices = intFeatures.at("peak_indices");
  const vector<int>& apendindices = intFeatures.at("AP_end_indices");

  vector<double> apfalltime;
  int retval = __AP_fall_time(t, stim_start, peakindices, apendindices, apfalltime);

  if (retval > 0) {
    setVec(DoubleFeatureData, StringData, "AP_fall_time", apfalltime);
  }
  return retval;
}

// *** AP_rise_rate according to E11 and E19 ***
static int __AP_rise_rate(const vector<double>& t, const vector<double>& v,
                          const vector<int>& apbeginindices,
                          const vector<int>& peakindices,
                          vector<double>& apriserate) {
  apriserate.resize(std::min(peakindices.size(), apbeginindices.size()));
  vector<int> newpeakindices;
  if (apbeginindices.size() > 0) {
    newpeakindices = peaks_after_stim_start(apbeginindices[0], peakindices);
  }
  for (size_t i = 0; i < apriserate.size(); i++) {
    apriserate[i] = (v[newpeakindices[i]] - v[apbeginindices[i]]) /
                    (t[newpeakindices[i]] - t[apbeginindices[i]]);
  }
  return apriserate.size();
}
int SpikeShape::AP_rise_rate(mapStr2intVec& IntFeatureData,
                        mapStr2doubleVec& DoubleFeatureData,
                        mapStr2Str& StringData) {
  const auto& doubleFeatures = getFeatures(DoubleFeatureData, {"T", "V"});
  const auto& intFeatures =
      getFeatures(IntFeatureData, {"AP_begin_indices", "peak_indices"});

  const vector<double>& t = doubleFeatures.at("T");
  const vector<double>& v = doubleFeatures.at("V");
  const vector<int>& apbeginindices = intFeatures.at("AP_begin_indices");
  const vector<int>& peakindices = intFeatures.at("peak_indices");

  vector<double> apriserate;
  int retval = __AP_rise_rate(t, v, apbeginindices, peakindices, apriserate);
  if (retval > 0) {
    setVec(DoubleFeatureData, StringData, "AP_rise_rate", apriserate);
  }
  return retval;
}

// *** AP_fall_rate according to E12 and E20 ***
static int __AP_fall_rate(const vector<double>& t, const vector<double>& v,
                          const double stimstart,
                          const vector<int>& peakindices,
                          const vector<int>& apendindices,
                          vector<double>& apfallrate) {
  apfallrate.resize(std::min(apendindices.size(), peakindices.size()));
  vector<int> newpeakindices = peaks_after_stim_start(stimstart, peakindices, t);

  for (size_t i = 0; i < apfallrate.size(); i++) {
    apfallrate[i] = (v[apendindices[i]] - v[newpeakindices[i]]) /
                    (t[apendindices[i]] - t[newpeakindices[i]]);
  }
  return apfallrate.size();
}
int SpikeShape::AP_fall_rate(mapStr2intVec& IntFeatureData,
                        mapStr2doubleVec& DoubleFeatureData,
                        mapStr2Str& StringData) {
  const auto& doubleFeatures = getFeatures(DoubleFeatureData, {"T", "V", "stim_start"});
  const auto& intFeatures =
      getFeatures(IntFeatureData, {"peak_indices", "AP_end_indices"});

  const vector<double>& t = doubleFeatures.at("T");
  const vector<double>& v = doubleFeatures.at("V");
  const double stim_start = doubleFeatures.at("stim_start")[0];
  const vector<int>& peakindices = intFeatures.at("peak_indices");
  const vector<int>& apendindices = intFeatures.at("AP_end_indices");

  vector<double> apfallrate;
  int retval = __AP_fall_rate(t, v, stim_start, peakindices, apendindices, apfallrate);

  if (retval > 0) {
    setVec(DoubleFeatureData, StringData, "AP_fall_rate", apfallrate);
  }
  return retval;
}

// *** AP_rise_rate_change according to E25 ***
static int __AP_rise_rate_change(const vector<double>& apriserate,
                                 vector<double>& apriseratechange) {
  if (apriserate.size() < 1) {
    return -1;
  }
  apriseratechange.resize(apriserate.size() - 1);

  for (size_t i = 0; i < apriseratechange.size(); i++) {
    apriseratechange[i] = (apriserate[i + 1] - apriserate[0]) / apriserate[0];
  }
  return apriseratechange.size();
}
int SpikeShape::AP_rise_rate_change(mapStr2intVec& IntFeatureData,
                               mapStr2doubleVec& DoubleFeatureData,
                               mapStr2Str& StringData) {
  const auto& features = getFeatures(DoubleFeatureData, {"AP_rise_rate"});
  const vector<double>& apriserate = features.at("AP_rise_rate");
  vector<double> apriseratechange;
  int retval = __AP_rise_rate_change(apriserate, apriseratechange);
  if (retval > 0) {
    setVec(DoubleFeatureData, StringData, "AP_rise_rate_change",
           apriseratechange);
  }
  return retval;
}

// *** AP_fall_rate_change according to E26 ***
static int __AP_fall_rate_change(const vector<double>& apfallrate,
                                 vector<double>& apfallratechange) {
  if (apfallrate.size() < 1) {
    return -1;
  }
  apfallratechange.resize(apfallrate.size() - 1);
  for (size_t i = 0; i < apfallratechange.size(); i++) {
    apfallratechange[i] = (apfallrate[i + 1] - apfallrate[0]) / apfallrate[0];
  }
  return apfallratechange.size();
}

int SpikeShape::AP_fall_rate_change(mapStr2intVec& IntFeatureData,
                               mapStr2doubleVec& DoubleFeatureData,
                               mapStr2Str& StringData) {
  const auto& features = getFeatures(DoubleFeatureData, {"AP_fall_rate"});
  const vector<double>& apfallrate = features.at("AP_fall_rate");
  vector<double> apfallratechange;
  int retval = __AP_fall_rate_change(apfallrate, apfallratechange);
  if (retval > 0) {
    setVec(DoubleFeatureData, StringData, "AP_fall_rate_change",
           apfallratechange);
  }
  return retval;
}

static int __AP_phaseslope(const vector<double>& v, const vector<double>& t,
                           double stimStart, double stimEnd,
                           vector<double>& ap_phaseslopes, vector<int> apbi,
                           int range) {
  vector<double> dvdt(v.size());
  vector<double> dv;
  vector<double> dt;
  int apbegin_index, range_max_index, range_min_index;
  double ap_phaseslope;
  getCentralDifferenceDerivative(1., v, dv);
  getCentralDifferenceDerivative(1., t, dt);
  transform(dv.begin(), dv.end(), dt.begin(), dvdt.begin(),
            std::divides<double>());

  for (size_t i = 0; i < apbi.size(); i++) {
    apbegin_index = apbi[i];
    range_min_index = apbegin_index - range;
    range_max_index = apbegin_index + range;
    if (range_min_index < 0 || range_max_index < 0) return -1;
    if (range_min_index > (int)t.size() || range_max_index > (int)t.size())
      return -1;
    if (v[range_max_index] - v[range_min_index] == 0) return -1;
    ap_phaseslope = (dvdt[range_max_index] - dvdt[range_min_index]) /
                    (v[range_max_index] - v[range_min_index]);
    ap_phaseslopes.push_back(ap_phaseslope);
    // printf("slope %f, mint %f, minv %f, mindvdt %f\n", ap_phaseslope,
    // t[range_min_index], v[range_min_index], dvdt[range_min_index]);
    // printf("slope %f, maxt %f, maxv %f, maxdvdt %f\n", ap_phaseslope,
    // t[range_max_index], v[range_max_index], dvdt[range_max_index]);
  }

  return ap_phaseslopes.size();
}
/// Calculate the slope of the V, dVdt plot at the beginning of every spike
/// (at the point where the derivative crosses the DerivativeThreshold)
int SpikeShape::AP_phaseslope(mapStr2intVec& IntFeatureData,
                         mapStr2doubleVec& DoubleFeatureData,
                         mapStr2Str& StringData) {
  const auto& doubleFeatures =
      getFeatures(DoubleFeatureData,
                  {"V", "T", "stim_start", "stim_end"});
  const auto& intFeatures = 
      getFeatures(IntFeatureData,
                  {"AP_begin_indices", "AP_phaseslope_range"});
  vector<double> ap_phaseslopes;
  int retVal = __AP_phaseslope(doubleFeatures.at("V"), doubleFeatures.at("T"),
                               doubleFeatures.at("stim_start")[0],
                               doubleFeatures.at("stim_end")[0], ap_phaseslopes,
                               intFeatures.at("AP_begin_indices"),
                               intFeatures.at("AP_phaseslope_range")[0]);

  if (retVal > 0) {
    setVec(DoubleFeatureData, StringData, "AP_phaseslope", ap_phaseslopes);
  }
  return retVal;
}
