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


#include "SpikeEvent.h"

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
using std::list;
using std::max_element;
using std::min_element;


static int __peak_indices(double threshold, const vector<double>& V,
                          const vector<double>& t, vector<int>& PeakIndex,
                          bool strict_stiminterval, double stim_start,
                          double stim_end) {
  vector<int> upVec, dnVec;
  double dtmp;
  size_t itmp = 0;
  bool itmp_set = false;

  for (size_t i = 1; i < V.size(); i++) {
    if (V[i] > threshold && V[i - 1] < threshold) {
      upVec.push_back(i);
    } else if (V[i] < threshold && V[i - 1] > threshold) {
      dnVec.push_back(i);
    }
  }
  if (dnVec.size() == 0)
    throw FeatureComputationError("Voltage never goes below or above threshold in spike detection.");
  if (upVec.size() == 0)
    throw FeatureComputationError("Voltage never goes above threshold in spike detection.");

  // case where voltage starts above threshold: remove 1st dnVec
  while (dnVec.size() > 0 && dnVec[0] < upVec[0]) {
    dnVec.erase(dnVec.begin());
  }

  if (upVec.size() > dnVec.size()) {
    size_t size_diff = upVec.size() - dnVec.size();
    for (size_t i = 0; i < size_diff; i++) {
      upVec.pop_back();
    }
  }

  PeakIndex.clear();
  int j = 0;
  for (size_t i = 0; i < upVec.size(); i++) {
    dtmp = -1e9;
    itmp = 0;
    itmp_set = false;
    EFEL_ASSERT(i < dnVec.size(), "dnVec array too small");
    for (j = upVec[i]; j <= dnVec[i]; j++) {
      if (dtmp < V[j]) {
        dtmp = V[j];
        itmp = j;
        itmp_set = true;
      }
    }
    if (itmp_set) {
      if (strict_stiminterval) {
        EFEL_ASSERT(itmp < t.size(), "peak_time falls outside of time array");
        if (t[itmp] >= stim_start && t[itmp] <= stim_end) {
          PeakIndex.push_back(itmp);
        }
      } else {
        PeakIndex.push_back(itmp);
      }
    }
  }
  return PeakIndex.size();
}

int SpikeEvent::peak_indices(mapStr2intVec& IntFeatureData,
                        mapStr2doubleVec& DoubleFeatureData,
                        mapStr2Str& StringData) {
  const auto& doubleFeatures = getFeatures(
      DoubleFeatureData, {"V", "T", "Threshold", "stim_start", "stim_end"});
  bool strict_stiminterval;
  try {
    const auto& intFeatures =
        getFeatures(IntFeatureData, {"strict_stiminterval"});
    strict_stiminterval = bool(intFeatures.at("strict_stiminterval")[0]);
  } catch (const std::runtime_error& e) {
    strict_stiminterval = false;
  }
  vector<int> PeakIndex;
  int retVal = __peak_indices(
      doubleFeatures.at("Threshold")[0], doubleFeatures.at("V"),
      doubleFeatures.at("T"), PeakIndex, strict_stiminterval,
      doubleFeatures.at("stim_start")[0], doubleFeatures.at("stim_end")[0]);

  if (retVal > 0) {
    setVec(IntFeatureData, StringData, "peak_indices", PeakIndex);
  }
  return retVal;
}

int SpikeEvent::peak_time(mapStr2intVec& IntFeatureData,
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
int SpikeEvent::first_spike_time(mapStr2intVec& IntFeatureData,
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

// time from stimulus start to second threshold crossing
int SpikeEvent::time_to_second_spike(mapStr2intVec& IntFeatureData,
                                mapStr2doubleVec& DoubleFeatureData,
                                mapStr2Str& StringData) {
  const auto doubleFeatures =
      getFeatures(DoubleFeatureData, {"peak_time", "stim_start"});
  const auto& peaktime = doubleFeatures.at("peak_time");
  const auto& stimstart = doubleFeatures.at("stim_start");
  if (peaktime.size() < 2)
    throw FeatureComputationError("Two spikes required for time_to_second_spike.");

  vector<double> second_spike = {peaktime[1] - stimstart[0]};
  setVec(DoubleFeatureData, StringData, "time_to_second_spike", second_spike);
  return 1;
}

// 1.0 over time to first spike (in Hz); returns 0 when no spike
int SpikeEvent::inv_time_to_first_spike(mapStr2intVec& IntFeatureData,
                                   mapStr2doubleVec& DoubleFeatureData,
                                   mapStr2Str& StringData) {
  vector<double> time_to_first_spike_vec =
      getFeature(DoubleFeatureData, "time_to_first_spike");
  vector<double> inv_time_to_first_spike_vec;

  double inv_time_to_first_spike = 1000.0 / time_to_first_spike_vec[0];
  inv_time_to_first_spike_vec.push_back(inv_time_to_first_spike);

  setVec(DoubleFeatureData, StringData, "inv_time_to_first_spike",
         inv_time_to_first_spike_vec);
  return 1;
}

// *** doublet_ISI ***
// value of the first ISI
int SpikeEvent::doublet_ISI(mapStr2intVec& IntFeatureData,
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

int SpikeEvent::all_ISI_values(mapStr2intVec& IntFeatureData,
                          mapStr2doubleVec& DoubleFeatureData,
                          mapStr2Str& StringData) {
  const vector<double>& peak_time = getFeature(DoubleFeatureData, "peak_time");
  if (peak_time.size() < 2)
    throw FeatureComputationError("Two spikes required for calculation of all_ISI_values.");

  vector<double> VecISI;
  for (size_t i = 1; i < peak_time.size(); i++) {
    VecISI.push_back(peak_time[i] - peak_time[i - 1]);
  }
  setVec(DoubleFeatureData, StringData, "all_ISI_values", VecISI);
  return VecISI.size();
}

// time from stimulus start to last spike
int SpikeEvent::time_to_last_spike(mapStr2intVec& IntFeatureData,
                              mapStr2doubleVec& DoubleFeatureData,
                              mapStr2Str& StringData) {
  const auto doubleFeatures =
      getFeatures(DoubleFeatureData, {"peak_time", "stim_start"});
  const auto& peaktime = doubleFeatures.at("peak_time");
  const auto& stimstart = doubleFeatures.at("stim_start");

  vector<double> last_spike = {peaktime.back() - stimstart[0]};

  setVec(DoubleFeatureData, StringData, "time_to_last_spike", last_spike);
  return 1;
}

static int __number_initial_spikes(const vector<double>& peak_times,
                                   double stimstart, double stimend,
                                   double initial_perc,
                                   vector<int>& number_initial_spikes) {
  double initialLength = (stimend - stimstart) * initial_perc;

  int startIndex =
      distance(peak_times.begin(),
               find_if(peak_times.begin(), peak_times.end(),
                       [stimstart](double t) { return t >= stimstart; }));
  int endIndex = distance(peak_times.begin(),
                          find_if(peak_times.begin(), peak_times.end(),
                                  [stimstart, initialLength](double t) {
                                    return t >= stimstart + initialLength;
                                  }));

  number_initial_spikes.push_back(endIndex - startIndex);

  return 1;
}

// Number of spikes in the initial_perc interval
int SpikeEvent::number_initial_spikes(mapStr2intVec& IntFeatureData,
                                 mapStr2doubleVec& DoubleFeatureData,
                                 mapStr2Str& StringData) {
  const auto& doubleFeatures =
      getFeatures(DoubleFeatureData,
                  {"peak_time", "initial_perc", "stim_start", "stim_end"});
  vector<int> number_initial_spikes;

  const vector<double>& peak_times = doubleFeatures.at("peak_time");
  const vector<double>& initial_perc = doubleFeatures.at("initial_perc");
  const vector<double>& stimstart = doubleFeatures.at("stim_start");
  const vector<double>& stimend = doubleFeatures.at("stim_end");

  if ((initial_perc[0] < 0) || (initial_perc[0] >= 1)) {
    throw FeatureComputationError("initial_perc should lie between [0 1).");
  }

  int retVal = __number_initial_spikes(peak_times, stimstart[0], stimend[0],
                                       initial_perc[0], number_initial_spikes);
  if (retVal >= 0) {
    setVec(IntFeatureData, StringData, "number_initial_spikes",
           number_initial_spikes);
  }
  return retVal;
}

int SpikeEvent::firing_rate(mapStr2intVec& IntFeatureData,
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

int SpikeEvent::adaptation_index(mapStr2intVec& IntFeatureData,
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

int SpikeEvent::adaptation_index2(mapStr2intVec& IntFeatureData,
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

// the recorded indices correspond to the peak indices
// does not skip the first ISI by default
static int __burst_indices(double burst_factor, int IgnoreFirstISI,
                           const vector<double> ISI_values,
                           vector<int>& burst_begin_indices,
                           vector<int>& burst_end_indices) {
  vector<double> ISIpcopy;
  vector<double>::iterator it1, it2;
  int n;
  double dMedian;
  bool in_burst;
  int first_ISI = IgnoreFirstISI, count = IgnoreFirstISI;

  burst_begin_indices.push_back(first_ISI);

  for (size_t i = first_ISI + 1; i < (ISI_values.size()); i++) {
    // get median
    ISIpcopy.clear();
    for (size_t j = count; j < i; j++) ISIpcopy.push_back(ISI_values[j]);
    sort(ISIpcopy.begin(), ISIpcopy.end());
    n = ISIpcopy.size();
    if ((n % 2) == 0) {
      dMedian =
          (ISIpcopy[int((n - 1) / 2)] + ISIpcopy[int((n - 1) / 2) + 1]) / 2;
    } else {
      dMedian = ISIpcopy[int(n / 2)];
    }

    in_burst = (burst_end_indices.size() == 0 ||
                burst_begin_indices.back() > burst_end_indices.back());

    // look for end burst
    if (in_burst && ISI_values[i] > (burst_factor * dMedian)) {
      burst_end_indices.push_back(i);
      count = i;
    }

    if (ISI_values[i] < ISI_values[i - 1] / burst_factor) {
      if (in_burst) {
        burst_begin_indices.back() = i;
      } else {
        burst_begin_indices.push_back(i);
      }
      count = i;
    }
  }

  in_burst = (burst_end_indices.size() == 0 ||
              burst_begin_indices.back() > burst_end_indices.back());
  if (in_burst) {
    burst_end_indices.push_back(ISI_values.size());
  }

  return burst_begin_indices.size();
}

int SpikeEvent::burst_begin_indices(mapStr2intVec& IntFeatureData,
                               mapStr2doubleVec& DoubleFeatureData,
                               mapStr2Str& StringData) {
  int retVal;
  vector<int> burst_begin_indices, burst_end_indices, retIgnore;
  vector<double> ISI_values, tVec;
  int IgnoreFirstISI;
  double burst_factor = 0;
  ISI_values = getFeature(DoubleFeatureData, "all_ISI_values");
  if (ISI_values.size() < 2) {
    GErrorStr +=
        "\nError: At least than 3 spikes are needed for burst calculation.\n";
    return -1;
  }
  retVal = getParam(DoubleFeatureData, "strict_burst_factor", tVec);
  if (retVal < 0)
    burst_factor = 2;
  else
    burst_factor = tVec[0];

  retVal = getParam(IntFeatureData, "ignore_first_ISI", retIgnore);
  if ((retVal == 1) && (retIgnore.size() > 0) && (retIgnore[0] == 0))
    IgnoreFirstISI = 0;
  else
    IgnoreFirstISI = 1;

  retVal = __burst_indices(burst_factor, IgnoreFirstISI, ISI_values,
                           burst_begin_indices, burst_end_indices);
  if (retVal >= 0) {
    setVec(IntFeatureData, StringData, "burst_begin_indices",
           burst_begin_indices);
    setVec(IntFeatureData, StringData, "burst_end_indices", burst_end_indices);
  }
  return retVal;
}

int SpikeEvent::burst_end_indices(mapStr2intVec& IntFeatureData,
                             mapStr2doubleVec& DoubleFeatureData,
                             mapStr2Str& StringData) {
  const auto& intFeatures = getFeatures(IntFeatureData, {"burst_begin_indices"});
  int retVal = intFeatures.at("burst_begin_indices").size();
  if (retVal <= 0) return -1;
  return retVal;
}

static int __strict_burst_mean_freq(const vector<double>& PVTime,
                                    const vector<int>& burst_begin_indices,
                                    const vector<int>& burst_end_indices,
                                    vector<double>& BurstMeanFreq) {
  if (burst_begin_indices.size() == 0) return BurstMeanFreq.size();
  double span;
  size_t i;

  for (i = 0; i < burst_begin_indices.size(); i++) {
    if (burst_end_indices[i] - burst_begin_indices[i] > 0) {
      span = PVTime[burst_end_indices[i]] - PVTime[burst_begin_indices[i]];
      BurstMeanFreq.push_back(
          (burst_end_indices[i] - burst_begin_indices[i] + 1) * 1000 / span);
    }
  }

  return BurstMeanFreq.size();
}

int SpikeEvent::strict_burst_mean_freq(mapStr2intVec& IntFeatureData,
                                  mapStr2doubleVec& DoubleFeatureData,
                                  mapStr2Str& StringData) {
  // Retrieve all required double and int features at once
  const auto& doubleFeatures = getFeatures(DoubleFeatureData, {"peak_time"});
  const auto& intFeatures =
      getFeatures(IntFeatureData, {"burst_begin_indices", "burst_end_indices"});

  vector<double> BurstMeanFreq;
  int retVal = __strict_burst_mean_freq(
      doubleFeatures.at("peak_time"), intFeatures.at("burst_begin_indices"),
      intFeatures.at("burst_end_indices"), BurstMeanFreq);

  if (retVal >= 0) {
    setVec(DoubleFeatureData, StringData, "strict_burst_mean_freq",
           BurstMeanFreq);
  }
  return retVal;
}

static int __strict_interburst_voltage(const vector<int>& burst_begin_indices,
                                       const vector<int>& PeakIndex,
                                       const vector<double>& T,
                                       const vector<double>& V,
                                       vector<double>& IBV) {
  if (burst_begin_indices.size() < 1) return 0;
  int j, pIndex, tsIndex, teIndex, cnt;
  double tStart, tEnd, vTotal = 0;
  for (size_t i = 1; i < burst_begin_indices.size(); i++) {
    pIndex = burst_begin_indices[i] - 1;
    tsIndex = PeakIndex[pIndex];
    tStart = T[tsIndex] + 5;  // 5 millisecond after
    pIndex = burst_begin_indices[i];
    teIndex = PeakIndex[pIndex];
    tEnd = T[teIndex] - 5;  // 5 millisecond before

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

int SpikeEvent::strict_interburst_voltage(mapStr2intVec& IntFeatureData,
                                     mapStr2doubleVec& DoubleFeatureData,
                                     mapStr2Str& StringData) {
  const auto& doubleFeatures = getFeatures(DoubleFeatureData, {"T", "V"});
  const auto& intFeatures =
      getFeatures(IntFeatureData, {"peak_indices", "burst_begin_indices"});
  vector<double> IBV;
  int retVal = __strict_interburst_voltage(
      intFeatures.at("burst_begin_indices"), intFeatures.at("peak_indices"),
      doubleFeatures.at("T"), doubleFeatures.at("V"), IBV);

  if (retVal > 0) {
    setVec(DoubleFeatureData, StringData, "strict_interburst_voltage", IBV);
  }
  return retVal;
}

static int __interburst_min_indices(const vector<double>& v,
                                    const vector<int>& peak_indices,
                                    const vector<int>& burst_end_indices,
                                    vector<int>& interburst_min_indices,
                                    vector<double>& interburst_min_values) {
  unsigned interburst_min_index;
  for (size_t i = 0; i < burst_end_indices.size() &&
                     burst_end_indices[i] + 1 < peak_indices.size();
       i++) {
    interburst_min_index =
        min_element(v.begin() + peak_indices[burst_end_indices[i]],
                    v.begin() + peak_indices[burst_end_indices[i] + 1]) -
        v.begin();

    interburst_min_indices.push_back(interburst_min_index);
    interburst_min_values.push_back(v[interburst_min_index]);
  }
  return interburst_min_indices.size();
}

int SpikeEvent::interburst_min_indices(mapStr2intVec& IntFeatureData,
                                  mapStr2doubleVec& DoubleFeatureData,
                                  mapStr2Str& StringData) {
  const auto& doubleFeatures = getFeatures(DoubleFeatureData, {"V"});
  const auto& intFeatures =
      getFeatures(IntFeatureData, {"peak_indices", "burst_end_indices"});
  vector<int> interburst_min_indices;
  vector<double> interburst_min_values;
  const vector<double>& v = doubleFeatures.at("V");
  const vector<int>& peak_indices = intFeatures.at("peak_indices");
  const vector<int>& burst_end_indices = intFeatures.at("burst_end_indices");
  int retVal =
      __interburst_min_indices(v, peak_indices, burst_end_indices,
                               interburst_min_indices, interburst_min_values);
  if (retVal > 0) {
    setVec(IntFeatureData, StringData, "interburst_min_indices",
           interburst_min_indices);
    setVec(DoubleFeatureData, StringData, "interburst_min_values",
           interburst_min_values);
  }
  return retVal;
}

int SpikeEvent::interburst_min_values(mapStr2intVec& IntFeatureData,
                                 mapStr2doubleVec& DoubleFeatureData,
                                 mapStr2Str& StringData) {
  const auto& intFeatures = getFeatures(IntFeatureData, {"interburst_min_indices"});
  int retVal = intFeatures.at("interburst_min_indices").size();
  if (retVal <= 0) return -1;
  return retVal;
}

static int __postburst_min_indices(const vector<double>& t,
                                   const vector<double>& v,
                                   const vector<int>& peak_indices,
                                   const vector<int>& burst_end_indices,
                                   vector<int>& postburst_min_indices,
                                   vector<double>& postburst_min_values,
                                   const double stim_end) {
  unsigned postburst_min_index, stim_end_index, end_index;
  stim_end_index = distance(
      t.begin(), find_if(t.begin(), t.end(),
                         [stim_end](double x) { return x >= stim_end; }));
  end_index = distance(t.begin(), t.end());
  for (size_t i = 0; i < burst_end_indices.size(); i++) {
    if (burst_end_indices[i] + 1 < peak_indices.size()){
      postburst_min_index = min_element(
        v.begin() + peak_indices[burst_end_indices[i]],
        v.begin() + peak_indices[burst_end_indices[i] + 1]
      ) - v.begin();
    } else if (peak_indices[burst_end_indices[i]] < stim_end_index){
      postburst_min_index = min_element(
        v.begin() + peak_indices[burst_end_indices[i]],
        v.begin() + stim_end_index
      ) - v.begin();
      if (postburst_min_index == stim_end_index){
        continue;
      }
    } else {
      postburst_min_index = min_element(
        v.begin() + peak_indices[burst_end_indices[i]],
        v.begin() + end_index
      ) - v.begin();
      if (postburst_min_index == end_index){
        continue;
      }
    }

    postburst_min_indices.push_back(postburst_min_index);
    postburst_min_values.push_back(v[postburst_min_index]);
  }

  return postburst_min_indices.size();
}

int SpikeEvent::postburst_min_indices(mapStr2intVec& IntFeatureData,
                                 mapStr2doubleVec& DoubleFeatureData,
                                 mapStr2Str& StringData) {
  const auto& doubleFeatures =
      getFeatures(DoubleFeatureData, {"T", "V", "stim_end"});
  const auto& intFeatures =
      getFeatures(IntFeatureData, {"peak_indices", "burst_end_indices"});
  vector<int> postburst_min_indices;
  vector<double> postburst_min_values;
  double stim_end = doubleFeatures.at("stim_end").front();
  int retVal = __postburst_min_indices(
      doubleFeatures.at("T"), doubleFeatures.at("V"),
      intFeatures.at("peak_indices"), intFeatures.at("burst_end_indices"),
      postburst_min_indices, postburst_min_values, stim_end);
  if (retVal > 0) {
    setVec(IntFeatureData, StringData, "postburst_min_indices",
           postburst_min_indices);
    setVec(DoubleFeatureData, StringData, "postburst_min_values",
           postburst_min_values);
  }
  return retVal;
}

int SpikeEvent::postburst_min_values(mapStr2intVec& IntFeatureData,
                                mapStr2doubleVec& DoubleFeatureData,
                                mapStr2Str& StringData) {
  const auto& intFeatures = getFeatures(IntFeatureData, {"postburst_min_indices"});
  int retVal = intFeatures.at("postburst_min_indices").size();
  if (retVal <= 0) return -1;
  return retVal;
}

int SpikeEvent::time_to_interburst_min(mapStr2intVec& IntFeatureData,
                                  mapStr2doubleVec& DoubleFeatureData,
                                  mapStr2Str& StringData) {
  // Retrieve all required double and int features at once
  const auto& doubleFeatures =
      getFeatures(DoubleFeatureData, {"T", "peak_time"});
  const auto& intFeatures = getFeatures(
      IntFeatureData, {"burst_end_indices", "interburst_min_indices"});

  vector<double> time_to_interburst_min;
  const vector<double>& time = doubleFeatures.at("T");
  const vector<double>& peak_time = doubleFeatures.at("peak_time");
  const vector<int>& burst_end_indices = intFeatures.at("burst_end_indices");
  const vector<int>& interburst_min_indices =
      intFeatures.at("interburst_min_indices");

  if (burst_end_indices.size() < interburst_min_indices.size()) {
    throw FeatureComputationError("burst_end_indices should not have less elements than interburst_min_indices");
  }

  for (size_t i = 0; i < interburst_min_indices.size(); i++) {
    time_to_interburst_min.push_back(time[interburst_min_indices[i]] -
                                     peak_time[burst_end_indices[i]]);
  }

  setVec(DoubleFeatureData, StringData, "time_to_interburst_min",
         time_to_interburst_min);
  return time_to_interburst_min.size();
}

static int __postburst_slow_ahp_indices(const vector<double>& t,
                                   const vector<double>& v,
                                   const vector<int>& peak_indices,
                                   const vector<int>& burst_end_indices,
                                   vector<int>& postburst_slow_ahp_indices,
                                   vector<double>& postburst_slow_ahp_values,
                                   const double stim_end,
                                   const double sahp_start) {
  unsigned postburst_slow_ahp_index, stim_end_index, end_index, t_start_index;
  stim_end_index =
      distance(t.begin(),
                find_if(t.begin(), t.end(),
                        [stim_end](double x) { return x >= stim_end; }));
  end_index = distance(t.begin(), t.end());
  for (size_t i = 0; i < burst_end_indices.size(); i++) {
    double t_start = t[peak_indices[burst_end_indices[i]]] + sahp_start;

    if (burst_end_indices[i] + 1 < peak_indices.size()){
      t_start_index = find_if(
        t.begin() + peak_indices[burst_end_indices[i]],
        t.begin() + peak_indices[burst_end_indices[i] + 1],
        [t_start](double x) { return x >= t_start; }
      ) - t.begin();
      postburst_slow_ahp_index = min_element(
        v.begin() + t_start_index,
        v.begin() + peak_indices[burst_end_indices[i] + 1]
      ) - v.begin();
    } else if (peak_indices[burst_end_indices[i]] < stim_end_index){
      t_start_index = find_if(
        t.begin() + peak_indices[burst_end_indices[i]],
        t.begin() + stim_end_index,
        [t_start](double x) { return x >= t_start; }
      ) - t.begin();
      if (t_start_index < stim_end_index){
        postburst_slow_ahp_index = min_element(
          v.begin() + t_start_index, v.begin() + stim_end_index
        ) - v.begin();
      } else {
        // edge case: stim_end_index is 1 index after stim_end
        continue;
      }
    } else {
      t_start_index = find_if(
        t.begin() + peak_indices[burst_end_indices[i]],
        t.begin() + end_index,
        [t_start](double x) { return x >= t_start; }
      ) - t.begin();
      if (t_start_index < end_index){
        postburst_slow_ahp_index = min_element(
          v.begin() + t_start_index, v.begin() + end_index
        ) - v.begin();
      } else{
        // edge case: end_index is 1 index after end
        continue;
      }
    }
    
    postburst_slow_ahp_indices.push_back(postburst_slow_ahp_index);
    postburst_slow_ahp_values.push_back(v[postburst_slow_ahp_index]);
  }

  return postburst_slow_ahp_indices.size();
}

int SpikeEvent::postburst_slow_ahp_indices(mapStr2intVec& IntFeatureData,
                                  mapStr2doubleVec& DoubleFeatureData,
                                  mapStr2Str& StringData) {
  const auto& doubleFeatures =
      getFeatures(DoubleFeatureData, {"T", "V", "stim_end", "sahp_start"});
  const auto& intFeatures =
      getFeatures(IntFeatureData, {"peak_indices", "burst_end_indices"});

  vector<double> postburst_slow_ahp_values;
  vector<int> postburst_slow_ahp_indices;
  int retVal = __postburst_slow_ahp_indices(
      doubleFeatures.at("T"),
      doubleFeatures.at("V"),
      intFeatures.at("peak_indices"),
      intFeatures.at("burst_end_indices"),
      postburst_slow_ahp_indices,
      postburst_slow_ahp_values,
      doubleFeatures.at("stim_end").front(),
      // time after the spike in ms after which to start searching for minimum
      doubleFeatures.at("sahp_start").empty() ? 5.0 : doubleFeatures.at("sahp_start").front()
  );
  if (retVal >= 0) {
    setVec(IntFeatureData, StringData, "postburst_slow_ahp_indices",
           postburst_slow_ahp_indices);
    setVec(DoubleFeatureData, StringData, "postburst_slow_ahp_values",
           postburst_slow_ahp_values);
  }
  return retVal;
}

int SpikeEvent::postburst_slow_ahp_values(mapStr2intVec& IntFeatureData,
                                 mapStr2doubleVec& DoubleFeatureData,
                                 mapStr2Str& StringData) {
  const auto& intFeatures = getFeatures(IntFeatureData, {"postburst_slow_ahp_indices"});
  int retVal = intFeatures.at("postburst_slow_ahp_indices").size();
  if (retVal <= 0) return -1;
  return retVal;
}

int SpikeEvent::time_to_postburst_slow_ahp(mapStr2intVec& IntFeatureData,
                                      mapStr2doubleVec& DoubleFeatureData,
                                      mapStr2Str& StringData) {
  const auto& doubleFeatures =
      getFeatures(DoubleFeatureData, {"T", "peak_time"});
  const auto& intFeatures =
      getFeatures(IntFeatureData, {"postburst_slow_ahp_indices", "burst_end_indices"});

  vector<double> time_to_postburst_slow_ahp;
  const vector<double>& time = doubleFeatures.at("T");
  const vector<double>& peak_time = doubleFeatures.at("peak_time");
  const vector<int>& burst_end_indices = intFeatures.at("burst_end_indices");
  const vector<int>& postburst_slow_ahp_indices = intFeatures.at("postburst_slow_ahp_indices");

  if (burst_end_indices.size() < postburst_slow_ahp_indices.size()){
    GErrorStr +=
        "\nburst_end_indices should not have less elements than postburst_slow_ahp_indices\n";
    return -1;
  }

  for (size_t i = 0; i < postburst_slow_ahp_indices.size(); i++) {
    time_to_postburst_slow_ahp.push_back(time[postburst_slow_ahp_indices[i]] -
                                         peak_time[burst_end_indices[i]]);
  }
  setVec(DoubleFeatureData, StringData, "time_to_postburst_slow_ahp",
         time_to_postburst_slow_ahp);
  return (time_to_postburst_slow_ahp.size());
}

static int __postburst_fast_ahp_indices(const vector<double>& t, const vector<double>& v,
                                        const vector<int>& peak_indices,
                                        const vector<int>& burst_end_indices,
                                        const double stim_end,
                                        vector<int>& postburst_fast_ahp_indices,
                                        vector<double>& postburst_fast_ahp_values) {
  vector<int> start_indices, end_indices;
  for (size_t i = 0; i < burst_end_indices.size(); i++) {
    start_indices.push_back(peak_indices[burst_end_indices[i]]);
    if (burst_end_indices[i] + 1 < peak_indices.size()){
      end_indices.push_back(peak_indices[burst_end_indices[i] + 1]);
    }
  }

  unsigned end_index = 0;
  if (t[start_indices.back()] < stim_end) {
    end_index =
        distance(t.begin(),
                 find_if(t.begin(), t.end(),
                         [stim_end](double x) { return x >= stim_end; }));
  } else {
    end_index = distance(t.begin(), t.end());
  }

  if (end_indices.size() < start_indices.size()){
    end_indices.push_back(end_index);
  }

  size_t fahpindex = 0;
  for (size_t i = 0; i < start_indices.size(); i++) {
    // can use first_min_element because dv/dt is very steep before fash ahp
    // and noise is very unlikely to make voltage go up before reaching fast ahp
    fahpindex = distance(
        v.begin(), first_min_element(v.begin() + start_indices[i],
                                     v.begin() + end_indices[i]));

    if (fahpindex != end_index - 1) {
      postburst_fast_ahp_indices.push_back(fahpindex);

      EFEL_ASSERT(fahpindex < v.size(),
                  "fast AHP index falls outside of voltage array");
      postburst_fast_ahp_values.push_back(v[fahpindex]);
    }
  }

  return postburst_fast_ahp_indices.size();
}

int SpikeEvent::postburst_fast_ahp_indices(mapStr2intVec& IntFeatureData,
                           mapStr2doubleVec& DoubleFeatureData,
                           mapStr2Str& StringData) {
  const auto& doubleFeatures =
      getFeatures(DoubleFeatureData, {"T", "V", "stim_end"});
  const auto& intFeatures =
      getFeatures(IntFeatureData, {"peak_indices", "burst_end_indices"});

  vector<int> postburst_fast_ahp_indices;
  vector<double> postburst_fast_ahp_values;
  int retVal =
      __postburst_fast_ahp_indices(
        doubleFeatures.at("T"),
        doubleFeatures.at("V"),
        intFeatures.at("peak_indices"),
        intFeatures.at("burst_end_indices"),
        doubleFeatures.at("stim_end").front(),
        postburst_fast_ahp_indices,
        postburst_fast_ahp_values
      );

  if (retVal > 0) {
    setVec(IntFeatureData, StringData, "postburst_fast_ahp_indices", postburst_fast_ahp_indices);
    setVec(DoubleFeatureData, StringData, "postburst_fast_ahp_values",
                 postburst_fast_ahp_values);
    return postburst_fast_ahp_indices.size();
  }
  return -1;
}

int SpikeEvent::postburst_fast_ahp_values(mapStr2intVec& IntFeatureData,
                                 mapStr2doubleVec& DoubleFeatureData,
                                 mapStr2Str& StringData) {
  const auto& intFeatures = getFeatures(IntFeatureData, {"postburst_fast_ahp_indices"});
  int retVal = intFeatures.at("postburst_fast_ahp_indices").size();
  if (retVal <= 0) return -1;
  return retVal;
}

static int __postburst_adp_peak_indices(const vector<double>& t, const vector<double>& v,
                                        const vector<int>& postburst_fast_ahp_indices,
                                        const vector<int>& postburst_slow_ahp_indices,
                                        vector<int>& postburst_adp_peak_indices,
                                        vector<double>& postburst_adp_peak_values) {

  if (postburst_slow_ahp_indices.size() > postburst_fast_ahp_indices.size()){
    GErrorStr +=
        "\n postburst_slow_ahp should not have more elements than "
        "postburst_fast_ahp for postburst_adp_peak_indices calculation.\n";
    return -1;
  }
  size_t adppeakindex = 0;
  for (size_t i = 0; i < postburst_slow_ahp_indices.size(); i++) {
    if (postburst_slow_ahp_indices[i] < postburst_fast_ahp_indices[i]){
      continue;
    }
    adppeakindex = distance(
        v.begin(), max_element(v.begin() + postburst_fast_ahp_indices[i],
                               v.begin() + postburst_slow_ahp_indices[i]));

    if (adppeakindex < postburst_slow_ahp_indices[i] - 1) {
      postburst_adp_peak_indices.push_back(adppeakindex);

      EFEL_ASSERT(adppeakindex < v.size(),
                  "ADP peak index falls outside of voltage array");
      postburst_adp_peak_values.push_back(v[adppeakindex]);
    }
  }

  return postburst_adp_peak_indices.size();
}

int SpikeEvent::postburst_adp_peak_indices(mapStr2intVec& IntFeatureData,
                           mapStr2doubleVec& DoubleFeatureData,
                           mapStr2Str& StringData) {
  const auto& doubleFeatures =
      getFeatures(DoubleFeatureData, {"T", "V"});
  const auto& intFeatures =
      getFeatures(IntFeatureData, {"postburst_fast_ahp_indices", "postburst_slow_ahp_indices"});
  vector<int> postburst_adp_peak_indices;
  vector<double> postburst_adp_peak_values;
  int retVal =
      __postburst_adp_peak_indices(
        doubleFeatures.at("T"),
        doubleFeatures.at("V"),
        intFeatures.at("postburst_fast_ahp_indices"),
        intFeatures.at("postburst_slow_ahp_indices"),
        postburst_adp_peak_indices,
        postburst_adp_peak_values
      );

  if (retVal > 0) {
    setVec(IntFeatureData, StringData, "postburst_adp_peak_indices", postburst_adp_peak_indices);
    setVec(DoubleFeatureData, StringData, "postburst_adp_peak_values",
                 postburst_adp_peak_values);
    return retVal;
  }
  return -1;
}

int SpikeEvent::postburst_adp_peak_values(mapStr2intVec& IntFeatureData,
                                 mapStr2doubleVec& DoubleFeatureData,
                                 mapStr2Str& StringData) {
  const auto& intFeatures = getFeatures(IntFeatureData, {"postburst_adp_peak_indices"});
  int retVal = intFeatures.at("postburst_adp_peak_indices").size();
  if (retVal <= 0) return -1;
  return retVal;
}

int SpikeEvent::time_to_postburst_fast_ahp(mapStr2intVec& IntFeatureData,
                                      mapStr2doubleVec& DoubleFeatureData,
                                      mapStr2Str& StringData) {
  const auto& doubleFeatures =
      getFeatures(DoubleFeatureData, {"T", "peak_time"});
  const auto& intFeatures =
      getFeatures(IntFeatureData, {"postburst_fast_ahp_indices", "burst_end_indices"});

  vector<double> time_to_postburst_fast_ahp;
  const vector<double>& time = doubleFeatures.at("T");
  const vector<double>& peak_time = doubleFeatures.at("peak_time");
  const vector<int>& burst_end_indices = intFeatures.at("burst_end_indices");
  const vector<int>& postburst_fast_ahp_indices = intFeatures.at("postburst_fast_ahp_indices");

  if (burst_end_indices.size() < postburst_fast_ahp_indices.size()){
    GErrorStr +=
        "\nburst_end_indices should not have less elements than postburst_fast_ahp_indices\n";
    return -1;
  }

  for (size_t i = 0; i < postburst_fast_ahp_indices.size(); i++) {
    time_to_postburst_fast_ahp.push_back(time[postburst_fast_ahp_indices[i]] -
                                         peak_time[burst_end_indices[i]]);
  }
  setVec(DoubleFeatureData, StringData, "time_to_postburst_fast_ahp",
         time_to_postburst_fast_ahp);
  return (time_to_postburst_fast_ahp.size());
}

int SpikeEvent::time_to_postburst_adp_peak(mapStr2intVec& IntFeatureData,
                                      mapStr2doubleVec& DoubleFeatureData,
                                      mapStr2Str& StringData) {
  const auto& doubleFeatures =
      getFeatures(DoubleFeatureData, {"T", "peak_time"});
  const auto& intFeatures =
      getFeatures(IntFeatureData, {"postburst_adp_peak_indices", "burst_end_indices"});

  vector<double> time_to_postburst_adp_peak;
  const vector<double>& time = doubleFeatures.at("T");
  const vector<double>& peak_time = doubleFeatures.at("peak_time");
  const vector<int>& burst_end_indices = intFeatures.at("burst_end_indices");
  const vector<int>& postburst_adp_peak_indices = intFeatures.at("postburst_adp_peak_indices");

  if (burst_end_indices.size() < postburst_adp_peak_indices.size()){
    GErrorStr +=
        "\nburst_end_indices should not have less elements than postburst_adp_peak_indices\n";
    return -1;
  }

  for (size_t i = 0; i < postburst_adp_peak_indices.size(); i++) {
    // there are not always an adp peak after each burst
    // so make sure that the burst and adp peak indices are consistent
    size_t k = 0;
    while (burst_end_indices[i] + k + 1 < peak_time.size() &&
           peak_time[burst_end_indices[i] + k + 1] < time[postburst_adp_peak_indices[i]]){
      k++;
    }
    time_to_postburst_adp_peak.push_back(time[postburst_adp_peak_indices[i]] -
                                        peak_time[burst_end_indices[i] + k]);
  }
  setVec(DoubleFeatureData, StringData, "time_to_postburst_adp_peak",
         time_to_postburst_adp_peak);
  return (time_to_postburst_adp_peak.size());
}

// index and voltage value at a given percentage of the duration of the interburst after fast AHP
int __interburst_percent_indices(const vector<double>& t, const vector<double>& v,
                                        const vector<int>& postburst_fast_ahp_indices,
                                        const vector<int>& peak_indices,
                                        const vector<int>& burst_end_indices,
                                        vector<int>& interburst_percent_indices,
                                        vector<double>& interburst_percent_values,
                                        // percentage should be a value between 0 and 1
                                        double fraction) {

  double time_interval, time_at_fraction;
  size_t index_at_fraction;
  for (size_t i = 0; i < postburst_fast_ahp_indices.size(); i++) {
    if (i < burst_end_indices.size()){
      if (burst_end_indices[i] + 1 < peak_indices.size()){
        time_interval = t[peak_indices[burst_end_indices[i] + 1]] - t[postburst_fast_ahp_indices[i]];
        time_at_fraction = t[postburst_fast_ahp_indices[i]] + time_interval * fraction;
        index_at_fraction =
          distance(t.begin(),
                    find_if(t.begin(), t.end(),
                            [time_at_fraction](double x){ return x >= time_at_fraction; }));
        interburst_percent_indices.push_back(index_at_fraction);
        interburst_percent_values.push_back(v[index_at_fraction]);
      }
    }
  }
  return interburst_percent_indices.size();
}

int SpikeEvent::interburst_XXpercent_indices(mapStr2intVec& IntFeatureData,
                           mapStr2doubleVec& DoubleFeatureData,
                           mapStr2Str& StringData, int percent) {
  const auto& doubleFeatures =
      getFeatures(DoubleFeatureData, {"T", "V"});
  const auto& intFeatures =
      getFeatures(IntFeatureData, {"peak_indices", "burst_end_indices", "postburst_fast_ahp_indices"});

  vector<int> interburst_XXpercent_indices;
  vector<double> interburst_XXpercent_values;
  char featureNameIndices[30], featureNameValues[30];
  sprintf(featureNameIndices, "interburst_%dpercent_indices", percent);
  sprintf(featureNameValues, "interburst_%dpercent_values", percent);

  int retVal =
      __interburst_percent_indices(
        doubleFeatures.at("T"),
        doubleFeatures.at("V"),
        intFeatures.at("postburst_fast_ahp_indices"),
        intFeatures.at("peak_indices"),
        intFeatures.at("burst_end_indices"),
        interburst_XXpercent_indices,
        interburst_XXpercent_values,
        percent / 100.
      );

  if (retVal > 0) {
    setVec(IntFeatureData, StringData, featureNameIndices, interburst_XXpercent_indices);
    setVec(DoubleFeatureData, StringData, featureNameValues,
                 interburst_XXpercent_values);
    return interburst_XXpercent_indices.size();
  }
  return -1;
}

static int __interburst_duration(const vector<double>& peak_time,
                                const vector<int>& burst_end_indices,
                                vector<double>& interburst_duration) {

  double duration;
  for (size_t i = 0; i < burst_end_indices.size(); i++) {
    if (burst_end_indices[i] + 1 < peak_time.size()){
      duration = peak_time[burst_end_indices[i] + 1] - peak_time[burst_end_indices[i]];
      interburst_duration.push_back(duration);
    }
  }
  return interburst_duration.size();
}

int SpikeEvent::interburst_duration(mapStr2intVec& IntFeatureData,
                              mapStr2doubleVec& DoubleFeatureData,
                              mapStr2Str& StringData) {
  const auto& doubleFeatures =
      getFeatures(DoubleFeatureData, {"peak_time"});
  const auto& intFeatures =
      getFeatures(IntFeatureData, {"burst_end_indices"});

  vector<double> interburst_duration;
  int retVal =
      __interburst_duration(
        doubleFeatures.at("peak_time"), intFeatures.at("burst_end_indices"), interburst_duration);

  if (retVal > 0) {
    setVec(DoubleFeatureData, StringData, "interburst_duration", interburst_duration);
    return interburst_duration.size();
  }
  return -1;
}

// Check if a cell is transiently stuck (i.e. not firing any spikes) at the end
// of
// retval will be -1 if the cell gets stuck, retval will be 1 if the cell
// doesn't get stuck
int SpikeEvent::is_not_stuck(mapStr2intVec& IntFeatureData,
                        mapStr2doubleVec& DoubleFeatureData,
                        mapStr2Str& StringData) {
  const vector<double>& peak_time = getFeature(DoubleFeatureData, "peak_time");
  const vector<double>& stim_start =
      getFeature(DoubleFeatureData, "stim_start");
  const vector<double>& stim_end = getFeature(DoubleFeatureData, "stim_end");
  bool stuck = true;
  for (const auto& pt : peak_time) {
    if (pt > stim_end[0] * 0.5 && pt < stim_end[0]) {
      stuck = false;
      break;
    }
  }
  if (!stuck) {
    vector<int> tc = {1};
    setVec(IntFeatureData, StringData, "is_not_stuck", tc);
    return tc.size();
  } else {
    return -1;
  }
}
