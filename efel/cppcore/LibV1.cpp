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

#include "LibV1.h"

#include <math.h>

#include <algorithm>
#include <cstdio>
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
using std::list;
using std::max_element;
using std::min_element;

template <typename T>
std::string to_string(const T& value) {
  std::ostringstream oss;
  oss << std::setprecision(17) << value;
  return oss.str();
}

int LibV1::interpolate(mapStr2intVec& IntFeatureData,
                       mapStr2doubleVec& DoubleFeatureData,
                       mapStr2Str& StringData) {
  vector<double> V, T, VIntrpol, TIntrpol, InterpStepVec;
  T = getFeature(DoubleFeatureData, "T");
  // interp_step is a stimulus independent parameter
  int retVal = getParam(DoubleFeatureData, "interp_step", InterpStepVec);
  double InterpStep = (retVal <= 0) ? 0.1 : InterpStepVec[0];

  try  // interpolate V if it's available
  {
    V = getFeature(DoubleFeatureData, "V");
    LinearInterpolation(InterpStep, T, V, TIntrpol, VIntrpol);
    setVec(DoubleFeatureData, StringData, "V", VIntrpol);
    setVec(DoubleFeatureData, StringData, "T", TIntrpol);
  } catch (...) {
    return -1;  // interpolation failed
  }

  // also interpolate current if present
  vector<double> I, IIntrpol, TIntrpolI;
  try {
    I = getFeature(DoubleFeatureData, "I");
    LinearInterpolation(InterpStep, T, I, TIntrpolI, IIntrpol);
    setVec(DoubleFeatureData, StringData, "I", IIntrpol);
    setVec(DoubleFeatureData, StringData, "T", TIntrpol);
  } catch (...) {
  }  // pass, it is optional
  return 1;
}

int LibV1::ISI_values(mapStr2intVec& IntFeatureData,
                      mapStr2doubleVec& DoubleFeatureData,
                      mapStr2Str& StringData) {
  const auto& doubleFeatures = getFeatures(DoubleFeatureData, {"peak_time"});
  if (doubleFeatures.at("peak_time").size() < 3)
    throw FeatureComputationError("Three spikes required for calculation of ISI_values.");

  const auto& intFeatures = getFeatures(IntFeatureData, {"ignore_first_ISI"});
  int IgnoreFirstISI = 1;
  if (intFeatures.at("ignore_first_ISI").size() > 0 &&
      intFeatures.at("ignore_first_ISI")[0] == 0) {
    IgnoreFirstISI = 0;
  }

  vector<double> VecISI;
  for (size_t i = IgnoreFirstISI + 1; i < doubleFeatures.at("peak_time").size();
       i++)
    VecISI.push_back(doubleFeatures.at("peak_time")[i] -
                     doubleFeatures.at("peak_time")[i - 1]);
  setVec(DoubleFeatureData, StringData, "ISI_values", VecISI);
  return VecISI.size();
}

// *** ISI_CV ***
// the coefficient of variation of the ISI
static int __ISI_CV(const vector<double>& isivalues, vector<double>& isicv) {
  // mean
  double isi_mean = 0.;
  for (size_t i = 0; i < isivalues.size(); i++) {
    isi_mean += isivalues[i];
  }
  isi_mean /= isivalues.size();

  // sigma^2
  double variance = 0.;
  for (size_t i = 0; i < isivalues.size(); i++) {
    double dev = isivalues[i] - isi_mean;
    variance += dev * dev;
  }
  // variation coefficient cv = sigma / mean
  isicv.push_back(sqrt(variance / (isivalues.size() - 1)) / isi_mean);
  return isicv.size();
}
int LibV1::ISI_CV(mapStr2intVec& IntFeatureData,
                  mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData) {
  const auto& doubleFeatures = getFeatures(DoubleFeatureData, {"ISI_values"});
  if (doubleFeatures.at("ISI_values").size() < 2) return -1;

  vector<double> isicv;
  int retval = __ISI_CV(doubleFeatures.at("ISI_values"), isicv);
  if (retval >= 0) {
    setVec(DoubleFeatureData, StringData, "ISI_CV", isicv);
  }
  return retval;
}

int LibV1::peak_voltage(mapStr2intVec& IntFeatureData,
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

int LibV1::firing_rate(mapStr2intVec& IntFeatureData,
                       mapStr2doubleVec& DoubleFeatureData,
                       mapStr2Str& StringData) {
  const auto& doubleFeatures =
      getFeatures(DoubleFeatureData, {"peak_time", "stim_start", "stim_end"});
  double lastAPTime = 0.;
  int nCount = 0;
  for (const auto& time : doubleFeatures.at("peak_time")) {
    if ((time >= doubleFeatures.at("stim_start")[0]) &&
        (time <= doubleFeatures.at("stim_end")[0])) {
      lastAPTime = time;
      nCount++;
    }
  }
  if (lastAPTime == doubleFeatures.at("stim_start")[0])
    throw FeatureComputationError("Prevent divide by zero.");
  vector<double> firing_rate;
  firing_rate.push_back(nCount * 1000 /
                        (lastAPTime - doubleFeatures.at("stim_start")[0]));
  setVec(DoubleFeatureData, StringData, "mean_frequency", firing_rate);
  return firing_rate.size();
}

int LibV1::peak_time(mapStr2intVec& IntFeatureData,
                     mapStr2doubleVec& DoubleFeatureData,
                     mapStr2Str& StringData) {
  const auto& doubleFeatures = getFeatures(DoubleFeatureData, {"T"});
  const auto& intFeatures = getFeatures(IntFeatureData, {"peak_indices"});
  vector<double> pvTime;
  for (const auto& index : intFeatures.at("peak_indices")) {
    pvTime.push_back(doubleFeatures.at("T")[index]);
  }
  setVec(DoubleFeatureData, StringData, "peak_time", pvTime);
  return pvTime.size();
}

// time from stimulus start to first threshold crossing
int LibV1::first_spike_time(mapStr2intVec& IntFeatureData,
                            mapStr2doubleVec& DoubleFeatureData,
                            mapStr2Str& StringData) {
  const auto& doubleFeatures =
      getFeatures(DoubleFeatureData, {"peak_time", "stim_start"});
  if (doubleFeatures.at("peak_time").size() < 1)
    throw FeatureComputationError("One spike required for time_to_first_spike.");
  vector<double> first_spike;
  first_spike.push_back(doubleFeatures.at("peak_time")[0] -
                        doubleFeatures.at("stim_start")[0]);
  setVec(DoubleFeatureData, StringData, "time_to_first_spike", first_spike);
  return first_spike.size();
}

int LibV1::AP_height(mapStr2intVec& IntFeatureData,
                     mapStr2doubleVec& DoubleFeatureData,
                     mapStr2Str& StringData) {
  const auto& doubleFeatures = getFeatures(DoubleFeatureData, {"peak_voltage"});
  setVec(DoubleFeatureData, StringData, "AP_height",
         doubleFeatures.at("peak_voltage"));
  return doubleFeatures.at("peak_voltage").size();
}

// spike amplitude: peak_voltage - v[AP_begin_indices]
int LibV1::AP_amplitude(mapStr2intVec& IntFeatureData,
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

int LibV1::AHP_depth_abs_slow(mapStr2intVec& IntFeatureData,
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
// end of AHP_depth_abs_slow

int LibV1::AHP_slow_time(mapStr2intVec& IntFeatureData,
                         mapStr2doubleVec& DoubleFeatureData,
                         mapStr2Str& StringData) {
  return -1;
}

static int __burst_ISI_indices(double BurstFactor, const vector<int>& PeakIndex,
                               const vector<double>& ISIValues,
                               vector<int>& BurstIndex) {
  vector<double> ISIpcopy;
  int n, count = -1;
  double dMedian;

  for (size_t i = 1; i < (ISIValues.size() - 1); i++) {
    ISIpcopy.clear();
    for (size_t j = count + 1; j < i; j++) ISIpcopy.push_back(ISIValues[j]);
    sort(ISIpcopy.begin(), ISIpcopy.end());
    n = ISIpcopy.size();
    if ((n % 2) == 0) {
      dMedian =
          (ISIpcopy[int((n - 1) / 2)] + ISIpcopy[int((n - 1) / 2) + 1]) / 2;
    } else {
      dMedian = ISIpcopy[int(n / 2)];
    }
    if (ISIValues[i] > (BurstFactor * dMedian) &&
        (ISIValues[i + 1] < ISIValues[i] / BurstFactor)) {
      BurstIndex.push_back(i + 1);
      count = i - 1;
    }
  }
  return BurstIndex.size();
}

int LibV1::burst_ISI_indices(mapStr2intVec& IntFeatureData,
                             mapStr2doubleVec& DoubleFeatureData,
                             mapStr2Str& StringData) {
  const auto& doubleFeatures =
      getFeatures(DoubleFeatureData, {"ISI_values", "burst_factor"});
  const auto& intFeatures = getFeatures(IntFeatureData, {"peak_indices"});
  if (intFeatures.at("peak_indices").size() < 5)
    throw FeatureComputationError("More than 5 spikes are needed for burst calculation.");
  double BurstFactor = (doubleFeatures.at("burst_factor").empty())
                           ? 2
                           : doubleFeatures.at("burst_factor")[0];
  vector<int> BurstIndex;
  int retVal = __burst_ISI_indices(BurstFactor, intFeatures.at("peak_indices"),
                                   doubleFeatures.at("ISI_values"), BurstIndex);
  if (retVal >= 0) {
    setVec(IntFeatureData, StringData, "burst_ISI_indices", BurstIndex);
  }
  return retVal;
}

// discrepancy betwen PVTime indices and BurstIndex indices because
// first ISI value is ignored
static int __burst_mean_freq(const vector<double>& PVTime,
                             const vector<int>& BurstIndex, int IgnoreFirstISI,
                             vector<double>& BurstMeanFreq) {
  // if no burst detected, do not consider all peaks in a single burst
  if (BurstIndex.size() == 0) return BurstMeanFreq.size();
  double span;
  size_t i;

  // 1st burst
  span = PVTime[BurstIndex[0] - 1 + IgnoreFirstISI] - PVTime[0];
  BurstMeanFreq.push_back((BurstIndex[0] + IgnoreFirstISI) * 1000 / span);

  for (i = 0; i < BurstIndex.size() - 1; i++) {
    if (BurstIndex[i + 1] - BurstIndex[i] > 1) {
      span = PVTime[BurstIndex[i + 1] - 1 + IgnoreFirstISI] -
             PVTime[BurstIndex[i] + IgnoreFirstISI];
      BurstMeanFreq.push_back((BurstIndex[i + 1] - BurstIndex[i]) * 1000 /
                              span);
    }
  }

  // last burst
  if (PVTime.size() - IgnoreFirstISI - BurstIndex[i] > 1) {
    span = PVTime[PVTime.size() - 1] - PVTime[BurstIndex[i] + IgnoreFirstISI];
    BurstMeanFreq.push_back((PVTime.size() - IgnoreFirstISI - BurstIndex[i]) *
                            1000 / span);
  }

  return BurstMeanFreq.size();
}

int LibV1::burst_mean_freq(mapStr2intVec& IntFeatureData,
                           mapStr2doubleVec& DoubleFeatureData,
                           mapStr2Str& StringData) {
  const auto& doubleFeatures = getFeatures(DoubleFeatureData, {"peak_time"});
  const auto& intFeatures =
      getFeatures(IntFeatureData, {"burst_ISI_indices"});
  int IgnoreFirstISI = 1;
  int retValIgnore;
  vector<int> retIgnore;
  retValIgnore = getParam(IntFeatureData, "ignore_first_ISI", retIgnore);
  if ((retValIgnore == 1) && (retIgnore.size() > 0) && (retIgnore[0] == 0)) {
    IgnoreFirstISI = 0;
  }
  vector<double> BurstMeanFreq;
  int retVal = __burst_mean_freq(doubleFeatures.at("peak_time"),
                                 intFeatures.at("burst_ISI_indices"),
                                 IgnoreFirstISI, BurstMeanFreq);
  if (retVal >= 0) {
    setVec(DoubleFeatureData, StringData, "burst_mean_freq", BurstMeanFreq);
  }
  return retVal;
}

// reminder: first ISI value is ignored in burst_ISI_indices if IgnoreFirstISI=1
static int __interburst_voltage(const vector<int>& BurstIndex,
                                const vector<int>& PeakIndex,
                                const vector<double>& T,
                                const vector<double>& V, int IgnoreFirstISI,
                                vector<double>& IBV) {
  if (BurstIndex.size() < 1) return 0;
  int j, pIndex, tsIndex, teIndex, cnt;
  double tStart, tEnd, vTotal = 0;
  for (size_t i = 0; i < BurstIndex.size(); i++) {
    pIndex = BurstIndex[i] + IgnoreFirstISI - 1;
    tsIndex = PeakIndex[pIndex];
    tStart = T[tsIndex] + 5;  // 5Millisecond after
    pIndex = BurstIndex[i] + IgnoreFirstISI;
    teIndex = PeakIndex[pIndex];
    tEnd = T[teIndex] - 5;  // 5Millisecond before

    for (j = tsIndex; j < teIndex; j++) {
      if (T[j] > tStart) break;
    }
    tsIndex = --j;

    for (j = teIndex; j > tsIndex; j--) {
      if (T[j] < tEnd) break;
    }
    teIndex = ++j;
    vTotal = 0;
    for (j = tsIndex, cnt = 1; j <= teIndex; j++, cnt++) vTotal = vTotal + V[j];
    if (cnt == 0) continue;
    IBV.push_back(vTotal / (cnt - 1));
  }
  return IBV.size();
}

int LibV1::interburst_voltage(mapStr2intVec& IntFeatureData,
                              mapStr2doubleVec& DoubleFeatureData,
                              mapStr2Str& StringData) {
  const auto& doubleFeatures = getFeatures(DoubleFeatureData, {"T", "V"});
  const auto& intFeatures =
      getFeatures(IntFeatureData, {"peak_indices", "burst_ISI_indices"});
  int IgnoreFirstISI;
  int retValIgnore;
  vector<int> retIgnore;
  retValIgnore = getParam(IntFeatureData, "ignore_first_ISI", retIgnore);
  if ((retValIgnore == 1) && (retIgnore.size() > 0) && (retIgnore[0] == 0)) {
    IgnoreFirstISI = 0;
  } else {
    IgnoreFirstISI = 1;
  }
  vector<double> IBV;
  int retVal = __interburst_voltage(
      intFeatures.at("burst_ISI_indices"), intFeatures.at("peak_indices"),
      doubleFeatures.at("T"), doubleFeatures.at("V"), IgnoreFirstISI, IBV);
  if (retVal >= 0) {
    setVec(DoubleFeatureData, StringData, "interburst_voltage", IBV);
  }
  return retVal;
}

static int __adaptation_index(double spikeSkipf, int maxnSpike,
                              double StimStart, double StimEnd, double Offset,
                              const vector<double>& peakVTime,
                              vector<double>& adaptation_index) {
  list<double> SpikeTime;
  vector<double> ISI;
  // Select spike time between given time scale (stim_start and stim_end )
  // consider Offset also if it is given as input
  for (size_t i = 0; i < peakVTime.size(); i++) {
    if ((peakVTime[i] >= (StimStart - Offset)) &&
        (peakVTime[i] <= (StimEnd + Offset))) {
      SpikeTime.push_back(peakVTime[i]);
    }
  }
  // Remove n spikes given by spike_skipf or max_spike_skip
  int spikeToRemove = (int)((SpikeTime.size() * spikeSkipf) + 0.5);
  // spike To remove is minimum of spike_skipf or max_spike_skip
  if (maxnSpike < spikeToRemove) {
    spikeToRemove = maxnSpike;
  }

  // Remove spikeToRemove spike from SpikeTime list
  for (int i = 0; i < spikeToRemove; ++i) {
    SpikeTime.pop_front();
  }

  // Adaptation index can not be calculated if nAPVec <4 or no of ISI is <3
  if (SpikeTime.size() < 4)
    throw FeatureComputationError("Minimum 4 spikes needed for feature [adaptation_index].");

  // Generate ISI vector
  list<double>::iterator lstItr = SpikeTime.begin();
  double lastValue = *lstItr;
  for (lstItr++; lstItr != SpikeTime.end(); lstItr++) {
    ISI.push_back(*lstItr - lastValue);
    lastValue = *lstItr;
  }

  // get addition and subtraction of ISIs
  double ISISum, ISISub, ADI;
  ADI = ISISum = ISISub = 0;
  for (size_t i = 1; i < ISI.size(); i++) {
    ISISum = ISI[i] + ISI[i - 1];
    ISISub = ISI[i] - ISI[i - 1];
    ADI = ADI + (ISISub / ISISum);
  }
  ADI = ADI / (ISI.size() - 1);
  adaptation_index.clear();
  adaptation_index.push_back(ADI);
  return 1;
}

int LibV1::adaptation_index(mapStr2intVec& IntFeatureData,
                            mapStr2doubleVec& DoubleFeatureData,
                            mapStr2Str& StringData) {
  const auto& doubleFeatures =
      getFeatures(DoubleFeatureData,
                  {"peak_time", "stim_start", "stim_end", "spike_skipf"});
  const auto& intFeatures = getFeatures(IntFeatureData, {"max_spike_skip"});

  if (doubleFeatures.at("spike_skipf")[0] < 0 ||
      doubleFeatures.at("spike_skipf")[0] >= 1) {
    throw FeatureComputationError("spike_skipf should lie between [0 1).");
  }
  vector<double> OffSetVec;
  double Offset;
  int retval = getParam(DoubleFeatureData, "offset", OffSetVec);
  if (retval < 0) {
    Offset = 0;  // Keep old behavior, set Offset to 0 if offset is not found
  } else {
    Offset = OffSetVec[0];  // Use the first element of OffSetVec if found
  }

  vector<double> adaptation_index;
  int retVal = __adaptation_index(
      doubleFeatures.at("spike_skipf")[0], intFeatures.at("max_spike_skip")[0],
      doubleFeatures.at("stim_start")[0], doubleFeatures.at("stim_end")[0],
      Offset, doubleFeatures.at("peak_time"), adaptation_index);
  if (retVal > 0) {
    setVec(DoubleFeatureData, StringData, "adaptation_index", adaptation_index);
  }
  return retVal;
}

// *** adaptation_index2 ***
// as adaptation_index, but start at the second ISI instead of the round(N *
// spikeskipf)
static int __adaptation_index2(double StimStart, double StimEnd, double Offset,
                               const vector<double>& peakVTime,
                               vector<double>& adaptation_index) {
  list<double> SpikeTime;
  vector<double> ISI;
  // Select spike time between given time scale (stim_start and stim_end )
  // considet Offset also if it is given as input
  for (size_t i = 0; i < peakVTime.size(); i++) {
    if ((peakVTime[i] >= (StimStart - Offset)) &&
        (peakVTime[i] <= (StimEnd + Offset))) {
      SpikeTime.push_back(peakVTime[i]);
    }
  }

  if (SpikeTime.size() < 4) {
    throw FeatureComputationError("At least 4 spikes within stimulus interval needed for adaptation_index2.");
  }
  // start at second ISI:
  SpikeTime.pop_front();

  // Generate ISI vector
  list<double>::iterator lstItr = SpikeTime.begin();
  double lastValue = *lstItr;
  for (++lstItr; lstItr != SpikeTime.end(); ++lstItr) {
    ISI.push_back(*lstItr - lastValue);
    lastValue = *lstItr;
  }

  // get addition and subtraction of ISIs
  double ISISum, ISISub, ADI;
  ADI = ISISum = ISISub = 0;
  for (size_t i = 1; i < ISI.size(); i++) {
    ISISum = ISI[i] + ISI[i - 1];
    ISISub = ISI[i] - ISI[i - 1];
    ADI = ADI + (ISISub / ISISum);
  }
  ADI = ADI / (ISI.size() - 1);
  adaptation_index.clear();
  adaptation_index.push_back(ADI);
  return 1;
}

int LibV1::adaptation_index2(mapStr2intVec& IntFeatureData,
                             mapStr2doubleVec& DoubleFeatureData,
                             mapStr2Str& StringData) {
  const auto& doubleFeatures =
      getFeatures(DoubleFeatureData, {"peak_time", "stim_start", "stim_end"});
  vector<double> OffSetVec;
  double Offset;
  int retval = getParam(DoubleFeatureData, "offset", OffSetVec);
  if (retval < 0)
    Offset = 0;
  else
    Offset = OffSetVec[0];

  if (doubleFeatures.at("peak_time").size() < 4) {
    throw FeatureComputationError("At least 4 spikes needed for adaptation_index2.");
  }

  vector<double> adaptationindex2;
  retval = __adaptation_index2(
      doubleFeatures.at("stim_start")[0], doubleFeatures.at("stim_end")[0],
      Offset, doubleFeatures.at("peak_time"), adaptationindex2);

  if (retval >= 0) {
    setVec(DoubleFeatureData, StringData, "adaptation_index2",
           adaptationindex2);
  }
  return retval;
}

// end of adaptation_index2

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

int LibV1::spike_width2(mapStr2intVec& IntFeatureData,
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

// passive properties implementation
//
// *** timeconstant ***
// requires hyperpolarizing stimulus
//
// the exponential fit works iteratively
//
static int __time_constant(const vector<double>& v, const vector<double>& t,
                           double stimStart, double stimEnd,
                           vector<double>& tc) {
  // value of the derivative near the minimum
  double min_derivative = 5e-3;
  // minimal required length of the decay (indices)
  size_t min_length = 10;
  // minimal required time length in ms
  int t_length = 70;
  size_t stimstartindex;
  for (stimstartindex = 0; t[stimstartindex] < stimStart; stimstartindex++)
    ;
  // increment stimstartindex to skip a possible transient
  stimstartindex += 10;
  // int stimendindex;
  // for(stimendindex = 0; t[stimendindex] < stimEnd; stimendindex++) ;
  // int stimmiddleindex = (stimstartindex + stimendindex) / 2;
  double mid_stim = (stimStart + stimEnd) / 2.;
  auto it_mid_stim =
      find_if(t.begin() + stimstartindex, t.end(),
              [mid_stim](double val) { return val >= mid_stim; });
  int stimmiddleindex = distance(t.begin(), it_mid_stim);

  if (stimstartindex >= v.size() || stimmiddleindex < 0 ||
      static_cast<size_t>(stimmiddleindex) >= v.size()) {
    return -1;
  }
  vector<double> part_v(&v[stimstartindex], &v[stimmiddleindex]);

  vector<double> part_t(&t[stimstartindex], &t[stimmiddleindex]);
  vector<double> dv;
  vector<double> dt;
  vector<double> dvdt(part_t.size());
  // calculate |dV/dt12
  // getCentralDifferenceDerivative(1.,part_v,dv);
  // getCentralDifferenceDerivative(1.,part_t,dt);
  getfivepointstencilderivative(part_v, dv);
  getfivepointstencilderivative(part_t, dt);
  transform(dv.begin(), dv.end(), dt.begin(), dvdt.begin(),
            std::divides<double>());
  // find start of the decay
  int i_start = 0;
  while (find_if(dvdt.begin() + i_start, dvdt.begin() + i_start + 5,
                 [min_derivative](double val) {
                   return val > -min_derivative;
                 }) != dvdt.begin() + i_start + 5) {
    if (dvdt.begin() + i_start + 5 == dvdt.end()) {
      throw FeatureComputationError("Could not find the decay.");
    }
    i_start++;
  }
  // find the flat
  // bool foundflat = false;
  int i_flat;
  for (i_flat = i_start;
       t[i_flat + stimstartindex] < t[stimmiddleindex] - t_length; i_flat++) {
    if (dvdt[i_flat] + min_derivative > 0.) {
      double sum = 0.;
      int length = 0;
      for (int i_mean = 0;
           t[i_flat + i_mean + stimstartindex] - t[i_flat + stimstartindex] <
           t_length;
           i_mean++) {
        sum += dvdt[i_flat + i_mean];
        length++;
      }
      double mean = sum / (double)length;
      if (mean + min_derivative > 0.) {
        // foundflat = true;
        break;
      }
    }
  }
  // if(!foundflat) {
  //    GErrorStr += "\nCould not locate plateau within range.\n";
  //    return -1;
  //}
  // containing the flat:
  vector<double> dvdt_decay(dvdt.begin() + i_start, dvdt.begin() + i_flat);
  vector<double> t_decay(part_t.begin() + i_start, part_t.begin() + i_flat);
  vector<double> v_decay(part_v.begin() + i_start, part_v.begin() + i_flat);
  if (dvdt_decay.size() < min_length) {
    throw FeatureComputationError("Trace fall time too short.");
  }

  // fit to exponential
  //
  vector<double> log_v(dvdt_decay.size(), 0.);

  // golden section search algorithm
  const double PHI = 1.618033988;
  vector<double> x(3, .0);
  // time_constant is searched in between 0 and 1000 ms
  x[2] = min_derivative * 1000.;
  x[1] = (x[0] * PHI + x[2]) / (1. + PHI);
  // calculate residuals at x[1]
  for (size_t i = 0; i < log_v.size(); i++) {
    log_v[i] = log(v_decay[i] - v_decay.back() + x[1]);
  }

  linear_fit_result fit;
  fit = slope_straight_line_fit(t_decay, log_v);
  double residuum = fit.normalized_std;
  bool right = true;
  double newx;
  while (x[2] - x[0] > .01) {
    // calculate new x value according to the golden section
    if (right) {
      newx = (x[1] * PHI + x[2]) / (1. + PHI);
    } else {
      newx = (x[0] + PHI * x[1]) / (1. + PHI);
    }
    // calculate residuals at newx
    for (size_t i = 0; i < log_v.size(); i++) {
      log_v[i] = log(v_decay[i] - v_decay.back() + newx);
    }
    fit = slope_straight_line_fit(t_decay, log_v);

    if (fit.normalized_std < residuum) {
      if (right) {
        x[0] = x[1];
        x[1] = newx;
      } else {
        x[2] = x[1];
        x[1] = newx;
      }
      residuum = fit.normalized_std;
    } else {
      if (right) {
        x[2] = newx;
      } else {
        x[0] = newx;
      }
      right = !right;
    }
  }
  tc.push_back(-1. / fit.slope);
  return 1;
}
int LibV1::time_constant(mapStr2intVec& IntFeatureData,
                         mapStr2doubleVec& DoubleFeatureData,
                         mapStr2Str& StringData) {
  const auto& doubleFeatures =
      getFeatures(DoubleFeatureData, {"V", "T", "stim_start", "stim_end"});
  vector<double> tc;
  int retVal = __time_constant(doubleFeatures.at("V"), doubleFeatures.at("T"),
                               doubleFeatures.at("stim_start")[0],
                               doubleFeatures.at("stim_end")[0], tc);
  if (retVal >= 0) {
    setVec(DoubleFeatureData, StringData, "time_constant", tc);
  }
  return retVal;
}

// *** voltage deflection ***

static int __voltage_deflection(const vector<double>& v,
                                const vector<double>& t, double stimStart,
                                double stimEnd, vector<double>& vd) {
  const size_t window_size = 5;

  size_t stimendindex = 0;
  double base = 0.;
  int base_size = 0;
  for (size_t i = 0; i < t.size(); i++) {
    if (t[i] < stimStart) {
      base += v[i];
      base_size++;
    }
    if (t[i] > stimEnd) {
      stimendindex = (int)i;
      break;
    }
  }
  if (base_size == 0) return -1;
  base /= base_size;
  double wind_mean = 0.;
  if (!(stimendindex >= 2 * window_size && v.size() > 0 &&
        stimendindex > window_size && stimendindex - window_size < v.size())) {
    return -1;
  }
  for (size_t i = stimendindex - 2 * window_size;
       i < stimendindex - window_size; i++) {
    wind_mean += v[i];
  }
  wind_mean /= window_size;
  vd.push_back(wind_mean - base);
  return 1;
}

int LibV1::voltage_deflection(mapStr2intVec& IntFeatureData,
                              mapStr2doubleVec& DoubleFeatureData,
                              mapStr2Str& StringData) {
  const auto& doubleFeatures =
      getFeatures(DoubleFeatureData, {"V", "T", "stim_start", "stim_end"});
  vector<double> vd;
  int retVal = __voltage_deflection(
      doubleFeatures.at("V"), doubleFeatures.at("T"),
      doubleFeatures.at("stim_start")[0], doubleFeatures.at("stim_end")[0], vd);
  if (retVal > 0) {
    setVec(DoubleFeatureData, StringData, "voltage_deflection", vd);
  }
  return retVal;
}

// *** ohmic input resistance ***

static int __ohmic_input_resistance(double voltage_deflection,
                                    double stimulus_current,
                                    vector<double>& oir) {
  if (stimulus_current == 0)
    throw FeatureComputationError("Stimulus current is zero which will result in division by zero.");
  oir.push_back(voltage_deflection / stimulus_current);
  return 1;
}

int LibV1::ohmic_input_resistance(mapStr2intVec& IntFeatureData,
                                  mapStr2doubleVec& DoubleFeatureData,
                                  mapStr2Str& StringData) {
  const auto& doubleFeatures = getFeatures(
      DoubleFeatureData, {"voltage_deflection", "stimulus_current"});
  vector<double> oir;
  int retVal =
      __ohmic_input_resistance(doubleFeatures.at("voltage_deflection")[0],
                               doubleFeatures.at("stimulus_current")[0], oir);
  if (retVal > 0) {
    setVec(DoubleFeatureData, StringData, "ohmic_input_resistance", oir);
  }
  return retVal;
}

static int __maxmin_voltage(const vector<double>& v, const vector<double>& t,
                            double stimStart, double stimEnd,
                            vector<double>& maxV, vector<double>& minV) {
  if (stimStart > t[t.size() - 1])
    throw FeatureComputationError("Stimulus start larger than max time in trace");

  if (stimEnd > t[t.size() - 1]) stimEnd = t.back();

  size_t stimstartindex = 0;
  while (t[stimstartindex] < stimStart && stimstartindex <= t.size())
    stimstartindex++;

  if (stimstartindex >= t.size()) {
    throw FeatureComputationError("Stimulus start index not found");
  }

  size_t stimendindex = 0;
  while (t[stimendindex] < stimEnd && stimstartindex <= t.size())
    stimendindex++;

  if (stimendindex >= t.size()) {
    throw FeatureComputationError("Stimulus end index not found");
  }

  maxV.push_back(*max_element(&v[stimstartindex], &v[stimendindex]));
  minV.push_back(*min_element(&v[stimstartindex], &v[stimendindex]));

  return 1;
}

int LibV1::maximum_voltage(mapStr2intVec& IntFeatureData,
                           mapStr2doubleVec& DoubleFeatureData,
                           mapStr2Str& StringData) {
  const auto& doubleFeatures =
      getFeatures(DoubleFeatureData, {"V", "T", "stim_start", "stim_end"});
  vector<double> maxV, minV;
  int retVal = __maxmin_voltage(doubleFeatures.at("V"), doubleFeatures.at("T"),
                                doubleFeatures.at("stim_start")[0],
                                doubleFeatures.at("stim_end")[0], maxV, minV);
  if (retVal > 0) {
    setVec(DoubleFeatureData, StringData, "maximum_voltage", maxV);
  }
  return retVal;
}

// *** maximum voltage ***
//
int LibV1::minimum_voltage(mapStr2intVec& IntFeatureData,
                           mapStr2doubleVec& DoubleFeatureData,
                           mapStr2Str& StringData) {
  const auto& doubleFeatures =
      getFeatures(DoubleFeatureData, {"V", "T", "stim_start", "stim_end"});
  vector<double> maxV, minV;
  int retVal = __maxmin_voltage(doubleFeatures.at("V"), doubleFeatures.at("T"),
                                doubleFeatures.at("stim_start")[0],
                                doubleFeatures.at("stim_end")[0], maxV, minV);
  if (retVal > 0) {
    setVec(DoubleFeatureData, StringData, "minimum_voltage", minV);
  }
  return retVal;
}

// *** steady state voltage ***
static int __steady_state_voltage(const vector<double>& v,
                                  const vector<double>& t, double stimEnd,
                                  vector<double>& ssv) {
  int mean_size = 0;
  double mean = 0;
  for (int i = t.size() - 1; t[i] > stimEnd; i--) {
    mean += v[i];
    mean_size++;
  }
  mean /= mean_size;
  ssv.push_back(mean);
  return 1;
}

int LibV1::steady_state_voltage(mapStr2intVec& IntFeatureData,
                                mapStr2doubleVec& DoubleFeatureData,
                                mapStr2Str& StringData) {
  const auto& doubleFeatures =
      getFeatures(DoubleFeatureData, {"V", "T", "stim_end"});
  if (doubleFeatures.at("stim_end").size() != 1) return -1;

  vector<double> ssv;
  int retVal =
      __steady_state_voltage(doubleFeatures.at("V"), doubleFeatures.at("T"),
                             doubleFeatures.at("stim_end")[0], ssv);
  if (retVal > 0) {
    setVec(DoubleFeatureData, StringData, "steady_state_voltage", ssv);
  }
  return retVal;
}

// *** single_burst_ratio ***
// according to Shaul: measures the length of the first isi over the median of
// the rest of the isis
static int __single_burst_ratio(const vector<double>& isivalues,
                                vector<double>& singleburstratio) {
  if (isivalues.size() < 2) {
    return 0;
  }
  // calculate the average instead of the median
  double average = 0.;
  for (size_t i = 1; i < isivalues.size(); i++) {
    average += isivalues[i];
  }
  average /= isivalues.size() - 1;
  singleburstratio.push_back(isivalues[0] / average);
  return singleburstratio.size();
}

int LibV1::single_burst_ratio(mapStr2intVec& IntFeatureData,
                              mapStr2doubleVec& DoubleFeatureData,
                              mapStr2Str& StringData) {
  const auto& doubleFeatures = getFeatures(DoubleFeatureData, {"ISI_values"});
  vector<double> singleburstratio;
  int retval =
      __single_burst_ratio(doubleFeatures.at("ISI_values"), singleburstratio);
  if (retval > 0) {
    setVec(DoubleFeatureData, StringData, "single_burst_ratio",
           singleburstratio);
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

int LibV1::AP_width(mapStr2intVec& IntFeatureData,
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
// end of AP_width

// *** doublet_ISI ***
// value of the first ISI
int LibV1::doublet_ISI(mapStr2intVec& IntFeatureData,
                       mapStr2doubleVec& DoubleFeatureData,
                       mapStr2Str& StringData) {
  const auto& doubleFeatures = getFeatures(DoubleFeatureData, {"peak_time"});
  if (doubleFeatures.at("peak_time").size() < 2) {
    throw FeatureComputationError("Need at least two spikes for doublet_ISI.");
  }
  vector<double> doubletisi(
      1, doubleFeatures.at("peak_time")[1] - doubleFeatures.at("peak_time")[0]);
  setVec(DoubleFeatureData, StringData, "doublet_ISI", doubletisi);
  return doubleFeatures.at("peak_time").size();
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

int LibV1::AHP_depth_slow(mapStr2intVec& IntFeatureData,
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
int LibV1::AHP_depth(mapStr2intVec& IntFeatureData,
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

int LibV1::AP_amplitude_diff(mapStr2intVec& IntFeatureData,
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
int LibV1::AHP_depth_diff(mapStr2intVec& IntFeatureData,
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
// end of feature definition
