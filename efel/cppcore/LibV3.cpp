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

#include "LibV3.h"

#include <algorithm>
#include <functional>
#include <list>
#include <math.h>

using std::bind2nd;
using std::find_if;
using std::greater_equal;
using std::less_equal;
using std::list;
using std::min_element;

int LibV3::interpolate(mapStr2intVec& IntFeatureData,
                       mapStr2doubleVec& DoubleFeatureData,
                       mapStr2Str& StringData) {
  int retVal, nSize;
  retVal = CheckInMap(IntFeatureData, StringData, "interpolate", nSize);
  if (retVal)
    return nSize;

  vector<double> V, T, VIntrpol, TIntrpol, InterpStepVec;
  vector<int> intrpolte;
  double InterpStep;
  // getVec takes care of stimulus suffix
  retVal = getVec(DoubleFeatureData, StringData, "V", V);
  if (retVal <= 0) return -1;
  retVal = getVec(DoubleFeatureData, StringData, "T", T);
  if (retVal <= 0) return -1;
  //
  // interp_step is a stimulus independent parameter
  retVal = getDoubleParam(DoubleFeatureData, "interp_step", InterpStepVec);
  if (retVal <= 0)
    InterpStep = 0.1;
  else
    InterpStep = InterpStepVec[0];

  LinearInterpolation(InterpStep, T, V, TIntrpol, VIntrpol);

  setVec(DoubleFeatureData, StringData, "V", VIntrpol);
  setVec(DoubleFeatureData, StringData, "T", TIntrpol);
  setVec(IntFeatureData, StringData, "interpolate", intrpolte);
  return retVal;
}

int LibV3::trace_check(mapStr2intVec& IntFeatureData,
                       mapStr2doubleVec& DoubleFeatureData,
                       mapStr2Str& StringData) {
  int retval;
  int size;
  retval = CheckInMap(IntFeatureData, StringData, "trace_check", size);
  if (retval)
    return size;

  vector<double> peak_time;
  vector<double> stim_start;
  vector<double> stim_end;
  vector<int> tc;
  retval = getVec(DoubleFeatureData, StringData, "peak_time", peak_time);
  if (retval < 0) return -1;
  retval = getVec(DoubleFeatureData, StringData, "stim_start", stim_start);
  if (retval < 0) return -1;
  retval = getVec(DoubleFeatureData, StringData, "stim_end", stim_end);
  if (retval < 0) return -1;

  bool sane = true;
  for (size_t i = 0; i < peak_time.size(); i++) {
    if (peak_time[i] < stim_start[0] || peak_time[i] > stim_end[0] * 1.05) {
      sane = false;
      break;
    }
  }

  if (sane) {
    tc.push_back(0);
    setVec(IntFeatureData, StringData, "trace_check", tc);
    return tc.size();
  } else {
    GErrorStr +=
        "Trace sanity check failed, there were spike outside the stimulus "
        "interval.\n";
    return -1;
  }
}

static int __peak_indices(double dThreshold, vector<double>& V,
                          vector<int>& PeakIndex) {
  vector<int> upVec, dnVec;
  double dtmp;
  int itmp;
  for (size_t i = 1; i < V.size(); i++) {
    if (V[i] > dThreshold && V[i - 1] < dThreshold) {
      upVec.push_back(i);
    } else if (V[i] < dThreshold && V[i - 1] > dThreshold) {
        dnVec.push_back(i);
    }
  }
  if ((dnVec.size() != upVec.size()) || (dnVec.size() == 0)) {
    GErrorStr += "\nBad Trace Shape.\n";
    return 0;
  }
  PeakIndex.clear();
  int j = 0;
  for (size_t i = 0; i < upVec.size(); i++) {
    dtmp = -1e9;
    itmp = -1;
    for (j = upVec[i]; j <= dnVec[i]; j++) {
      if (dtmp < V[j]) {
        dtmp = V[j];
        itmp = j;
      }
    }
    if (itmp != -1) PeakIndex.push_back(itmp);
  }
  return PeakIndex.size();
}

int LibV3::peak_indices(mapStr2intVec& IntFeatureData,
                        mapStr2doubleVec& DoubleFeatureData,
                        mapStr2Str& StringData) {
  // printf("\n  LibV1  This is inside peak_indices()");
  int retVal, nSize;
  retVal =
      CheckInMap(IntFeatureData, StringData, "peak_indices", nSize);
  if (retVal)
    return nSize;

  vector<int> PeakIndex;
  vector<double> v, Th;
  retVal = getVec(DoubleFeatureData, StringData, "V", v);
  if (retVal <= 0) return -1;
  retVal = getDoubleParam(DoubleFeatureData, "Threshold", Th);
  if (retVal <= 0) return -1;

  int retval = __peak_indices(Th[0], v, PeakIndex);
  if (retval >= 0)
    setVec(IntFeatureData, StringData, "peak_indices", PeakIndex);
  return retval;
}

int LibV3::ISI_values(mapStr2intVec& IntFeatureData,
                      mapStr2doubleVec& DoubleFeatureData,
                      mapStr2Str& StringData) {
  //  printf("\n  LibV1  This is inside ISI_values()");
  int retVal, nSize;
  retVal = CheckInMap(DoubleFeatureData, StringData, "ISI_values",
                            nSize);
  if (retVal)
    return nSize;

  vector<double> VecISI, pvTime;
  retVal = getVec(DoubleFeatureData, StringData, "peak_time", pvTime);
  if (retVal < 3) {
    GErrorStr += "\n Three spikes required for calculation of ISI_values.\n";
    return -1;
  }
  for (size_t i = 2; i < pvTime.size(); i++) {
    VecISI.push_back(pvTime[i] - pvTime[i - 1]);
  }
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
int LibV3::ISI_CV(mapStr2intVec& IntFeatureData,
                  mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData) {
  int retval;
  int nsize;
  retval =
      CheckInMap(DoubleFeatureData, StringData, "ISI_CV", nsize);
  if (retval) {
    return nsize;
  }
  vector<double> isivalues;
  retval = getVec(DoubleFeatureData, StringData, "ISI_values",
                        isivalues);
  if (retval < 2) return -1;
  vector<double> isicv;
  retval = __ISI_CV(isivalues, isicv);
  if (retval >= 0) {
    setVec(DoubleFeatureData, StringData, "ISI_CV", isicv);
  }
  return retval;
}

int LibV3::peak_voltage(mapStr2intVec& IntFeatureData,
                        mapStr2doubleVec& DoubleFeatureData,
                        mapStr2Str& StringData) {
  // printf("\n  LibV1  Inside PeakVoltage ... LibV1");
  int retVal, nSize;
  retVal = CheckInMap(DoubleFeatureData, StringData,
                            "peak_voltage", nSize);
  if (retVal)
    return nSize;

  // vector<int> PeakI = getVec(IntFeatureData, StringData,
  // string("peak_indices"));
  vector<int> PeakI;
  vector<double> V, peakV;
  retVal = getVec(IntFeatureData, StringData, "peak_indices", PeakI);
  if (retVal <= 0) return -1;
  retVal = getVec(DoubleFeatureData, StringData, "V", V);
  if (retVal <= 0) return -1;

  for (size_t i = 0; i < PeakI.size(); i++) {
    peakV.push_back(V[PeakI[i]]);
  }
  setVec(DoubleFeatureData, StringData, "peak_voltage", peakV);
  return peakV.size();
}

int LibV3::firing_rate(mapStr2intVec& IntFeatureData,
                       mapStr2doubleVec& DoubleFeatureData,
                       mapStr2Str& StringData) {
  // printf("\n  LibV1  This is inside firing_rate()");
  int retVal, nSize;
  retVal = CheckInMap(DoubleFeatureData, StringData,
                            "mean_frequency", nSize);
  if (retVal)
    return nSize;

  vector<double> stimStart, stimEnd, peakVTime, firing_rate;
  double lastAPTime = 0.;
  retVal = getVec(DoubleFeatureData, StringData, "peak_time", peakVTime);
  if (retVal <= 0) return -1;
  retVal = getVec(DoubleFeatureData, StringData, "stim_start", stimStart);
  if (retVal <= 0) return -1;
  retVal = getVec(DoubleFeatureData, StringData, "stim_end", stimEnd);
  if (retVal <= 0) return -1;

  int nCount = 0;
  for (size_t i = 0; i < peakVTime.size(); i++) {
    if ((peakVTime[i] >= stimStart[0]) && (peakVTime[i] <= stimEnd[0])) {
      lastAPTime = peakVTime[i];
      nCount++;
    }
  }
  if (lastAPTime == stimStart[0]) {
    GErrorStr += "\nPrevent divide by zero.\n";
    return -1;
  }
  firing_rate.push_back(nCount * 1000 / (lastAPTime - stimStart[0]));
  firing_rate.push_back(nCount * 1000 / (lastAPTime - stimStart[0]));
  setVec(DoubleFeatureData, StringData, "mean_frequency", firing_rate);
  return firing_rate.size();
}

int LibV3::peak_time(mapStr2intVec& IntFeatureData,
                     mapStr2doubleVec& DoubleFeatureData,
                     mapStr2Str& StringData) {
  // printf("\n  LibV1  This is inside peak_time()");
  int retVal, nSize;
  retVal = CheckInMap(DoubleFeatureData, StringData, "peak_time", nSize);
  if (retVal)
    return nSize;

  vector<int> PeakI;
  vector<double> T, pvTime;
  retVal = getVec(IntFeatureData, StringData, "peak_indices", PeakI);
  if (retVal <= 0) return -1;
  retVal = getVec(DoubleFeatureData, StringData, "T", T);
  if (retVal <= 0) return -1;

  for (size_t i = 0; i < PeakI.size(); i++) {
    pvTime.push_back(T[PeakI[i]]);
  }
  setVec(DoubleFeatureData, StringData, "peak_time", pvTime);
  return pvTime.size();
}

// time from stimulus start to first threshold crossing
int LibV3::first_spike_time(mapStr2intVec& IntFeatureData,
                            mapStr2doubleVec& DoubleFeatureData,
                            mapStr2Str& StringData) {
  int retVal, nSize;
  retVal = CheckInMap(DoubleFeatureData, StringData,
                            "time_to_first_spike", nSize);
  if (retVal)
    return nSize;

  vector<double> first_spike;
  vector<double> peaktime;
  vector<double> stimstart;
  retVal = getVec(DoubleFeatureData, StringData, "peak_time", peaktime);
  if (retVal < 1) {
    GErrorStr += "\n One spike required for time_to_first_spike.\n";
    return -1;
  }
  retVal = getVec(DoubleFeatureData, StringData, "stim_start", stimstart);
  if (retVal <= 0) return -1;

  first_spike.push_back(peaktime[0] - stimstart[0]);
  setVec(DoubleFeatureData, StringData, "time_to_first_spike",
               first_spike);
  return first_spike.size();
}

// spike half width
// for spike amplitude = v_peak - v_AHP
static int __spike_width1(const vector<double>& t, const vector<double>& v,
                          const vector<int>& peak_indices,
                          const vector<int>& min_ahp_indices, double stim_start,
                          vector<double>& spike_width1) {
  int start_index = distance(
      t.begin(), find_if(t.begin(), t.end(),
                         bind2nd(greater_equal<double>(), stim_start)));
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
    int rise_index =
        distance(v.begin(), find_if(v.begin() + min_ahp_indices_plus[i - 1],
                                    v.begin() + peak_indices[i - 1],
                                    bind2nd(greater_equal<double>(), v_half)));
    v_dev = v_half - v[rise_index];
    delta_v = v[rise_index] - v[rise_index - 1];
    delta_t = t[rise_index] - t[rise_index - 1];
    t_dev_rise = delta_t * v_dev / delta_v;
    int fall_index =
        distance(v.begin(), find_if(v.begin() + peak_indices[i - 1],
                                    v.begin() + min_ahp_indices_plus[i],
                                    bind2nd(less_equal<double>(), v_half)));
    v_dev = v_half - v[fall_index];
    delta_v = v[fall_index] - v[fall_index - 1];
    delta_t = t[fall_index] - t[fall_index - 1];
    t_dev_fall = delta_t * v_dev / delta_v;
    spike_width1.push_back(t[fall_index] + t_dev_rise - t[rise_index] +
                           t_dev_fall);
  }
  return spike_width1.size();
}

int LibV3::spike_width1(mapStr2intVec& IntFeatureData,
                        mapStr2doubleVec& DoubleFeatureData,
                        mapStr2Str& StringData) {
  int retVal, nSize;
  retVal = CheckInMap(DoubleFeatureData, StringData,
                            "spike_half_width", nSize);
  if (retVal)
    return nSize;

  vector<int> PeakIndex, minAHPIndex;
  vector<double> V, t, dv1, dv2, spike_width1;
  vector<double> stim_start;
  retVal = getVec(DoubleFeatureData, StringData, "V", V);
  if (retVal < 0) return -1;
  retVal = getVec(DoubleFeatureData, StringData, "T", t);
  if (retVal < 0) return -1;
  retVal = getVec(DoubleFeatureData, StringData, "stim_start", stim_start);
  if (retVal < 0) return -1;
  retVal = getVec(IntFeatureData, StringData, "min_AHP_indices", minAHPIndex);
  if (retVal < 0) return -1;
  retVal = getVec(IntFeatureData, StringData, "peak_indices", PeakIndex);
  if (retVal < 0) return -1;

  if (PeakIndex.size() <= 1) {
    GErrorStr += "\nError: More than one spike is needed for spikewidth "
        "calculation.\n";
    return -1;
  }

  // Take derivative of voltage from 1st AHPmin to the peak of the spike
  // Using Central difference derivative vec1[i] = ((vec[i+1]+vec[i-1])/2)/dx
  retVal = __spike_width1(t, V, PeakIndex, minAHPIndex, stim_start[0],
                          spike_width1);
  if (retVal >= 0) {
    setVec(DoubleFeatureData, StringData, "spike_half_width",
                 spike_width1);
  }
  return retVal;
}

// min_AHP_indices
// find the minimum between two spikes,
// and the minimum between the last spike and the time the stimulus ends
int LibV3::min_AHP_indices(mapStr2intVec& IntFeatureData,
                           mapStr2doubleVec& DoubleFeatureData,
                           mapStr2Str& StringData) {
  int retVal, nSize;
  retVal = CheckInMap(IntFeatureData, StringData, "min_AHP_indices",
                         nSize);
  if (retVal)
    return nSize;

  vector<int> peak_indices_plus;
  vector<int> min_ahp_indices;
  vector<double> v;
  vector<double> min_ahp_values;
  vector<double> stim_end;
  vector<double> t;
  retVal = getVec(DoubleFeatureData, StringData, "V", v);
  if (retVal <= 0) return -1;
  retVal = getVec(IntFeatureData, StringData, "peak_indices",
                     peak_indices_plus);
  if (retVal < 1) {
    GErrorStr +=
        "\n At least one spike required for calculation of "
        "min_AHP_indices.\n";
    return -1;
  }
  retVal = getVec(DoubleFeatureData, StringData, "stim_end", stim_end);
  if (retVal <= 0) return -1;
  retVal = getVec(DoubleFeatureData, StringData, "T", t);
  if (retVal <= 0) return -1;

  int end_index = distance(
      t.begin(), find_if(t.begin(), t.end(),
                         bind2nd(greater_equal<double>(), stim_end[0])));
  // if the last spike happens to be close to the end of the stimulus
  // there will not be a proper AHP, this case is not properly dealt with here
  if (end_index > peak_indices_plus.back() + 5) {
    peak_indices_plus.push_back(end_index);
  }
  for (size_t i = 0; i < peak_indices_plus.size() - 1; i++) {
    int ahpindex = distance(
        v.begin(), min_element(v.begin() + peak_indices_plus[i],
                               v.begin() + peak_indices_plus[i + 1]));
    min_ahp_indices.push_back(ahpindex);
    min_ahp_values.push_back(v[ahpindex]);
  }
  setVec(IntFeatureData, StringData, "min_AHP_indices", min_ahp_indices);
  setVec(DoubleFeatureData, StringData, "min_AHP_values",
               min_ahp_values);
  return min_ahp_indices.size();
}

int LibV3::min_AHP_values(mapStr2intVec& IntFeatureData,
                          mapStr2doubleVec& DoubleFeatureData,
                          mapStr2Str& StringData) {
  int retVal, nSize;
  retVal = CheckInMap(DoubleFeatureData, StringData,
                            "min_AHP_values", nSize);
  if (retVal) return nSize;
  return -1;
}

// AHP_depth_abs
// naming conflict here AHP_depth_abs does the same as min_AHP_values.
// In my opinion the should not be a feature called 'AHP_depth_abs',
// use min_AHP_values instead.
// A more interesting feature would be 'AHP_depth' anyways, which calculates the
// depth of the AHP relative to the voltage base
int LibV3::AHP_depth_abs(mapStr2intVec& IntFeatureData,
                         mapStr2doubleVec& DoubleFeatureData,
                         mapStr2Str& StringData) {
  int retVal, nSize;
  retVal = CheckInMap(DoubleFeatureData, StringData,
                            "AHP_depth_abs", nSize);
  if (retVal)
    return nSize;

  vector<double> vAHP;
  retVal = getVec(DoubleFeatureData, StringData, "min_AHP_values", vAHP);
  if (retVal <= 0) return -1;

  setVec(DoubleFeatureData, StringData, "AHP_depth_abs", vAHP);
  return vAHP.size();
}

int LibV3::rest_voltage_value(mapStr2intVec& IntFeatureData,
                              mapStr2doubleVec& DoubleFeatureData,
                              mapStr2Str& StringData) {
  int retVal, nSize;
  retVal = CheckInMap(DoubleFeatureData, StringData,
                            "voltage_base", nSize);
  if (retVal)
    return nSize;

  vector<double> v, t, stimStart, vRest;
  double startTime, endTime;
  retVal = getVec(DoubleFeatureData, StringData, "V", v);
  if (retVal < 0) return -1;
  retVal = getVec(DoubleFeatureData, StringData, "T", t);
  if (retVal < 0) return -1;
  retVal = getVec(DoubleFeatureData, StringData, "stim_start", stimStart);
  if (retVal < 0) return -1;
  startTime = stimStart[0] * .25;  // It is 25% from start (0), so if stimulus
                                   // starts at 100ms then StartTime will be
                                   // 25mS
  // as in above case end time will be 25% less than startTime
  endTime = stimStart[0] * .75;
  int nCount = 0;
  double vSum = 0;
  // calculte the mean of voltage between startTime and endTime
  for (size_t i = 0; i < t.size(); i++) {
    if (t[i] >= startTime) {
      vSum = vSum + v[i];
      nCount++;
    }
    if (t[i] > endTime) {
      break;
    }
  }
  vRest.push_back(vSum / nCount);
  setVec(DoubleFeatureData, StringData, "voltage_base", vRest);
  return 1;
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
    GErrorStr +=
        "\n At least 4 spikes within stimulus interval needed for "
        "adaptation_index2.\n";
    return -1;
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

int LibV3::adaptation_index2(mapStr2intVec& IntFeatureData,
                             mapStr2doubleVec& DoubleFeatureData,
                             mapStr2Str& StringData) {
  int retval;
  int nsize;
  retval = CheckInMap(DoubleFeatureData, StringData,
                            "adaptation_index2", nsize);
  if (retval) {
    return nsize;
  }
  vector<double> peakvoltagetime;
  retval = getVec(DoubleFeatureData, StringData, "peak_time",
                        peakvoltagetime);
  if (retval < 4) {
    GErrorStr += "\n At least 4 spikes needed for adaptation_index2.\n";
    return -1;
  }
  vector<double> stimStart;
  retval = getVec(DoubleFeatureData, StringData, "stim_start", stimStart);
  {
    if (retval < 0) return -1;
  };
  vector<double> stimEnd;
  retval = getVec(DoubleFeatureData, StringData, "stim_end", stimEnd);
  {
    if (retval < 0) return -1;
  };
  vector<double> OffSetVec;
  double Offset;
  retval = getDoubleParam(DoubleFeatureData, "offset", OffSetVec);
  if (retval < 0)
    Offset = 0;
  else
    Offset = OffSetVec[0];
  vector<double> adaptationindex2;
  retval = __adaptation_index2(stimStart[0], stimEnd[0], Offset,
                               peakvoltagetime, adaptationindex2);
  if (retval >= 0) {
    setVec(DoubleFeatureData, StringData, "adaptation_index2",
                 adaptationindex2);
  }
  return retval;
}
// end of adaptation_index2

int LibV3::AP_height(mapStr2intVec& IntFeatureData,
                     mapStr2doubleVec& DoubleFeatureData,
                     mapStr2Str& StringData) {
  int retVal, nSize;
  retVal = CheckInMap(DoubleFeatureData, StringData, "AP_height",
                            nSize);
  if (retVal)
    return nSize;

  vector<double> vPeak;
  retVal = getVec(DoubleFeatureData, StringData, "peak_voltage", vPeak);
  if (retVal <= 0) return -1;
  setVec(DoubleFeatureData, StringData, "AP_height", vPeak);
  return vPeak.size();
}

// spike amplitude: peak_voltage - v[AP_begin_indices]
int LibV3::AP_amplitude(mapStr2intVec& IntFeatureData,
                        mapStr2doubleVec& DoubleFeatureData,
                        mapStr2Str& StringData) {
  int retVal, nSize;
  retVal =
      CheckInMap(DoubleFeatureData, StringData, "AP_amplitude", nSize);
  if (retVal > 0)
    return nSize;

  vector<double> peakvoltage;
  vector<int> apbeginindices;
  vector<double> v;
  retVal = getVec(DoubleFeatureData, StringData, "V", v);
  if (retVal <= 0) { return -1; }
  retVal = getVec(DoubleFeatureData, StringData, "peak_voltage",
                        peakvoltage);
  if (retVal <= 0) { return -1; }
  retVal = getVec(IntFeatureData, StringData, "AP_begin_indices",
                     apbeginindices);
  if (retVal <= 0) { return -1; }

  vector<double> apamplitude;
  apamplitude.resize(peakvoltage.size());
  for (size_t i = 0; i < apamplitude.size(); i++) {
    apamplitude[i] = peakvoltage[i] - v[apbeginindices[i]];
  }
  setVec(DoubleFeatureData, StringData, "AP_amplitude", apamplitude);
  return apamplitude.size();
}

// *** AP_width ***
//
// spike width calculation according to threshold value for the first spike
// unfortunately spike width means the width of the spike on onset not at half
// maximum
static int __AP_width(const vector<double>& t, const vector<double>& v,
                      double stimstart, double threshold,
                      const vector<int>& peakindices,
                      const vector<int>& minahpindices,
                      vector<double>& apwidth) {
  //   printf("\n Inside AP_width...\n");
  //   printVectorD("t", t);
  //   printVectorD("v", v);
  //   printVectorI("peakindices", peakindices);
  //   printVectorI("minahpindices", minahpindices);
  //   printf("\nStimStart = %f , thereshold = %f ", stimstart, threshold);
  vector<int> indices(minahpindices.size() + 1);
  int start_index = distance(
      t.begin(),
      find_if(t.begin(), t.end(), bind2nd(greater_equal<double>(), stimstart)));
  indices[0] = start_index;
  copy(minahpindices.begin(), minahpindices.end(), indices.begin() + 1);
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
    int onset_index = distance(
        v.begin(), find_if(v.begin() + indices[i], v.begin() + indices[i + 1],
                           bind2nd(greater_equal<double>(), threshold)));
    // int end_index = distance(v.begin(), find_if(v.begin() + peakindices[i],
    // v.begin() + indices[i + 1], bind2nd(less_equal<double>(), threshold)));
    int end_index = distance(
        v.begin(), find_if(v.begin() + onset_index, v.begin() + indices[i + 1],
                           bind2nd(less_equal<double>(), threshold)));
    apwidth.push_back(t[end_index] - t[onset_index]);
  }
  return apwidth.size();
}

int LibV3::AP_width(mapStr2intVec& IntFeatureData,
                    mapStr2doubleVec& DoubleFeatureData,
                    mapStr2Str& StringData) {
  int retval;
  int nsize;
  retval = CheckInMap(DoubleFeatureData, StringData, "AP_width",
                            nsize);
  if (retval) {
    return nsize;
  }
  vector<double> t;
  retval = getVec(DoubleFeatureData, StringData, "T", t);
  if (retval < 0) return -1;
  vector<double> v;
  retval = getVec(DoubleFeatureData, StringData, "V", v);
  if (retval < 0) return -1;
  vector<double> threshold;
  retval = getDoubleParam(DoubleFeatureData, "Threshold", threshold);
  if (retval < 0) return -1;
  vector<double> stimstart;
  retval = getVec(DoubleFeatureData, StringData, "stim_start", stimstart);
  if (retval < 0) return -1;
  vector<int> peakindices;
  retval = getVec(IntFeatureData, StringData, "peak_indices", peakindices);
  if (retval <= 0) {
    GErrorStr += "\nNo spike in trace.\n";
    return -1;
  }
  vector<int> minahpindices;
  retval = getVec(IntFeatureData, StringData, "min_AHP_indices", minahpindices);
  if (retval < 0) return -1;
  vector<double> apwidth;
  retval = __AP_width(t, v, stimstart[0], threshold[0], peakindices,
                      minahpindices, apwidth);
  if (retval >= 0) {
    setVec(DoubleFeatureData, StringData, "AP_width", apwidth);
  }
  return retval;
}
// end of AP_width

// *** doublet_ISI ***
// value of the first ISI
int LibV3::doublet_ISI(mapStr2intVec& IntFeatureData,
                       mapStr2doubleVec& DoubleFeatureData,
                       mapStr2Str& StringData) {
  int retval;
  int nsize;
  retval = CheckInMap(DoubleFeatureData, StringData,
                            "doublet_ISI", nsize);
  if (retval) {
    return nsize;
  }
  vector<double> pvt;
  retval = getVec(DoubleFeatureData, StringData, "peak_time", pvt);
  if (retval < 2) {
    GErrorStr += "\nNeed at least two spikes for doublet_ISI.\n";
    return -1;
  }
  vector<double> doubletisi(1, pvt[1] - pvt[0]);
  setVec(DoubleFeatureData, StringData, "doublet_ISI", doubletisi);
  return retval;
}
// end of doublet_ISI

// AP parameters
//
// *** AP begin indices ***
//
static int __AP_begin_indices(const vector<double>& t, const vector<double>& v,
                              double stimstart, double stimend,
                              const vector<int>& ahpi, vector<int>& apbi) {
  // derivative at peak start according to eCode specification 10mV/ms
  // according to Shaul 12mV/ms
  const double derivativethreshold = 12.;
  vector<double> dvdt(v.size());
  vector<double> dv;
  vector<double> dt;
  getCentralDifferenceDerivative(1., v, dv);
  getCentralDifferenceDerivative(1., t, dt);
  transform(dv.begin(), dv.end(), dt.begin(), dvdt.begin(), std::divides<double>());

  // restrict to time interval where stimulus is applied
  vector<int> minima;
  int stimbeginindex = distance(
      t.begin(),
      find_if(t.begin(), t.end(), bind2nd(greater_equal<double>(), stimstart)));
  minima.push_back(stimbeginindex);
  for (size_t i = 0; i < ahpi.size(); i++) {
    if (ahpi[i] > stimbeginindex) {
      minima.push_back(ahpi[i]);
    }
    if (t[ahpi[i]] > stimend) {
      break;
    }
  }
  // if the AHP_indices are already restricted make sure that we do not miss
  // the last spike
  if (t[minima.back()] < stimend) {
    int stimendindex =
        distance(t.begin(), find_if(t.begin() + minima.back(), t.end(),
                                    bind2nd(greater_equal<double>(), stimend)));
    minima.push_back(stimendindex);
  }
  for (size_t i = 0; i < minima.size() - 1; i++) {
    // assure that the width of the slope is bigger than 4
    int newbegin = minima[i];
    int begin = minima[i];
    int width = 5;
    bool skip = false;
    do {
      begin = distance(
          dvdt.begin(),
          find_if(dvdt.begin() + newbegin, dvdt.begin() + minima[i + 1],
                  bind2nd(greater_equal<double>(), derivativethreshold)));
      if (begin == minima[i + 1]) {
        // could not find a spike in between these minima
        skip = true;
        break;
      }
      newbegin = begin + 1;
    } while (find_if(dvdt.begin() + begin, dvdt.begin() + begin + width,
                     bind2nd(std::less<double>(), derivativethreshold)) !=
             dvdt.begin() + begin + width);
    if (skip) {
      continue;
    }
    apbi.push_back(begin);
  }
  return apbi.size();
}
int LibV3::AP_begin_indices(mapStr2intVec& IntFeatureData,
                            mapStr2doubleVec& DoubleFeatureData,
                            mapStr2Str& StringData) {
  int retVal;
  int nSize;
  retVal = CheckInMap(IntFeatureData, StringData, "AP_begin_indices",
                         nSize);
  if (retVal) {
    return nSize;
  }
  vector<double> t;
  retVal = getVec(DoubleFeatureData, StringData, "T", t);
  if (retVal < 0) return -1;
  vector<double> v;
  retVal = getVec(DoubleFeatureData, StringData, "V", v);
  if (retVal < 0) return -1;
  vector<double> stimstart;
  retVal = getVec(DoubleFeatureData, StringData, "stim_start", stimstart);
  if (retVal < 0) return -1;
  vector<double> stimend;
  retVal = getVec(DoubleFeatureData, StringData, "stim_end", stimend);
  if (retVal < 0) return -1;
  vector<int> ahpi;
  retVal = getVec(IntFeatureData, StringData, "min_AHP_indices", ahpi);
  if (retVal < 0) return -1;
  vector<int> apbi;
  retVal = __AP_begin_indices(t, v, stimstart[0], stimend[0], ahpi, apbi);
  if (retVal >= 0) {
    setVec(IntFeatureData, StringData, "AP_begin_indices", apbi);
  }
  return retVal;
}

// *** AP end indices ***
//
static int __AP_end_indices(const vector<double>& t, const vector<double>& v,
                            const vector<int>& pi, vector<int>& apei) {
  // derivative at peak end according to eCode specification -10mV/ms
  // according to Shaul -12mV/ms
  const double derivativethreshold = -12.;
  // assume constant time steps
  double timestep = t[1] - t[0];
  vector<double> dvdt;
  getCentralDifferenceDerivative(timestep, v, dvdt);

  apei.resize(pi.size());
  vector<int> picopy(pi.begin(), pi.end());
  picopy.push_back(v.size() - 1);
  for (size_t i = 0; i < apei.size(); i++) {
    // assure that the width of the slope is bigger than 4
    apei[i] = distance(
        dvdt.begin(),
        find_if(dvdt.begin() + picopy[i] + 1, dvdt.begin() + picopy[i + 1],
                bind2nd(greater_equal<double>(), derivativethreshold)));
  }
  return apei.size();
}
int LibV3::AP_end_indices(mapStr2intVec& IntFeatureData,
                          mapStr2doubleVec& DoubleFeatureData,
                          mapStr2Str& StringData) {
  int retVal;
  int nSize;
  retVal = CheckInMap(IntFeatureData, StringData, "AP_end_indices",
                         nSize);
  if (retVal) {
    return nSize;
  }

  vector<double> t;
  retVal = getVec(DoubleFeatureData, StringData, "T", t);
  if (retVal < 0) return -1;
  vector<double> v;
  retVal = getVec(DoubleFeatureData, StringData, "V", v);
  if (retVal < 0) return -1;
  vector<int> pi;
  retVal = getVec(IntFeatureData, StringData, "peak_indices", pi);
  if (retVal < 0) return -1;
  vector<int> apei;
  retVal = __AP_end_indices(t, v, pi, apei);
  if (retVal >= 0) {
    setVec(IntFeatureData, StringData, "AP_end_indices", apei);
  }
  return retVal;
}

// *** AP rise indices ***
//
static int __AP_rise_indices(const vector<double>& v, const vector<int>& apbi,
                             const vector<int>& pi, vector<int>& apri) {
  apri.resize(apbi.size());
  for (size_t i = 0; i < apri.size(); i++) {
    double halfheight = (v[pi[i]] + v[apbi[i]]) / 2.;
    vector<double> vpeak;
    vpeak.resize(pi[i] - apbi[i]);
    transform(v.begin() + apbi[i], v.begin() + pi[i], vpeak.begin(),
              bind2nd(std::minus<double>(), halfheight));
    transform(vpeak.begin(), vpeak.end(), vpeak.begin(), 
              static_cast<double(*)(double)>(fabs));
    apri[i] = distance(vpeak.begin(), min_element(vpeak.begin(), vpeak.end())) +
              apbi[i];
  }
  return apri.size();
}
int LibV3::AP_rise_indices(mapStr2intVec& IntFeatureData,
                           mapStr2doubleVec& DoubleFeatureData,
                           mapStr2Str& StringData) {
  int retVal;
  int nSize;
  retVal = CheckInMap(IntFeatureData, StringData, "AP_rise_indices",
                         nSize);
  if (retVal) {
    return nSize;
  }

  vector<double> v;
  retVal = getVec(DoubleFeatureData, StringData, "V", v);
  if (retVal < 0) return -1;
  vector<int> apbi;
  retVal = getVec(IntFeatureData, StringData, "AP_begin_indices", apbi);
  if (retVal < 0) return -1;
  vector<int> pi;
  retVal = getVec(IntFeatureData, StringData, "peak_indices", pi);
  if (retVal < 0) return -1;
  vector<int> apri;
  retVal = __AP_rise_indices(v, apbi, pi, apri);
  if (retVal >= 0) {
    setVec(IntFeatureData, StringData, "AP_rise_indices", apri);
  }
  return retVal;
}

// *** AP fall indices ***
//
static int __AP_fall_indices(const vector<double>& v, const vector<int>& apbi,
                             const vector<int>& apei, const vector<int>& pi,
                             vector<int>& apfi) {
  apfi.resize(apbi.size());
  for (size_t i = 0; i < apfi.size(); i++) {
    double halfheight = (v[pi[i]] + v[apbi[i]]) / 2.;
    vector<double> vpeak(&v[pi[i]], &v[apei[i]]);
    transform(vpeak.begin(), vpeak.end(), vpeak.begin(),
              bind2nd(std::minus<double>(), halfheight));
    transform(vpeak.begin(), vpeak.end(), vpeak.begin(), 
              static_cast<double(*)(double)>(fabs));
    apfi[i] = distance(vpeak.begin(), min_element(vpeak.begin(), vpeak.end())) +
              pi[i];
  }
  return apfi.size();
}
int LibV3::AP_fall_indices(mapStr2intVec& IntFeatureData,
                           mapStr2doubleVec& DoubleFeatureData,
                           mapStr2Str& StringData) {
  int retVal;
  int nSize;
  retVal = CheckInMap(IntFeatureData, StringData, "AP_fall_indices", nSize);
  if (retVal) {
    return nSize;
  }

  vector<double> v;
  retVal = getVec(DoubleFeatureData, StringData, "V", v);
  if (retVal < 0) return -1;
  vector<int> apbi;
  retVal = getVec(IntFeatureData, StringData, "AP_begin_indices", apbi);
  if (retVal < 0) return -1;
  vector<int> apei;
  retVal = getVec(IntFeatureData, StringData, "AP_end_indices", apei);
  if (retVal < 0) return -1;
  vector<int> pi;
  retVal = getVec(IntFeatureData, StringData, "peak_indices", pi);
  if (retVal < 0) return -1;
  vector<int> apfi;
  retVal = __AP_fall_indices(v, apbi, apei, pi, apfi);
  if (retVal >= 0) {
    setVec(IntFeatureData, StringData, "AP_fall_indices", apfi);
  }
  return retVal;
}

// eFeatures
// *** AP_duration according to E7 and E15 ***
static int __AP_duration(const vector<double>& t,
                         const vector<int>& apbeginindices,
                         const vector<int>& endindices,
                         vector<double>& apduration) {
  apduration.resize(apbeginindices.size());
  for (size_t i = 0; i < apduration.size(); i++) {
    apduration[i] = t[endindices[i]] - t[apbeginindices[i]];
  }
  return apduration.size();
}
int LibV3::AP_duration(mapStr2intVec& IntFeatureData,
                       mapStr2doubleVec& DoubleFeatureData,
                       mapStr2Str& StringData) {
  int retval;
  int nsize;
  retval = CheckInMap(DoubleFeatureData, StringData,
                            "AP_duration", nsize);
  if (retval) {
    return nsize;
  }
  vector<double> t;
  retval = getVec(DoubleFeatureData, StringData, "T", t);
  if (retval < 0) return -1;
  vector<int> apbeginindices;
  retval = getVec(IntFeatureData, StringData, "AP_begin_indices",
                     apbeginindices);
  if (retval < 0) return -1;
  vector<int> endindices;
  retval = getVec(IntFeatureData, StringData, "AP_end_indices",
                     endindices);
  if (retval < 0) return -1;
  vector<double> apduration;
  retval = __AP_duration(t, apbeginindices, endindices, apduration);
  if (retval >= 0) {
    setVec(DoubleFeatureData, StringData, "AP_duration", apduration);
  }
  return retval;
}
// end of AP_duration

static int __depolarized_base(const vector<double>& t, const vector<double>& v,
                              double stimstart, double stimend,
                              const vector<int>& apbi,
                              const vector<int>& apendi,
                              vector<double>& dep_base) {
  int i, n, k, startIndex, endIndex, nPt;
  double baseValue;
  // to make sure it access minimum index of both length
  if (apendi.size() < apbi.size())
    n = apendi.size();
  else
    n = apbi.size();

  if (apendi.size() == apbi.size()) n = apendi.size() - 1;

  if (n > 2) {
    dep_base.clear();
    for (i = 0; i < n; i++) {
      nPt = 0;
      baseValue = 0;
      startIndex = apendi[i];
      endIndex = apbi[i + 1];
      for (k = startIndex; k < endIndex; k++) {
        baseValue += v[k];
        nPt++;
      }
      baseValue = baseValue / nPt;
      dep_base.push_back(baseValue);
    }
    return dep_base.size();
  }
  return -1;
}

int LibV3::depolarized_base(mapStr2intVec& IntFeatureData,
                            mapStr2doubleVec& DoubleFeatureData,
                            mapStr2Str& StringData) {
  int retVal;
  int nSize;
  retVal = CheckInMap(IntFeatureData, StringData, "depolarized_base",
                         nSize);
  if (retVal) {
    return nSize;
  }
  vector<double> t;
  retVal = getVec(DoubleFeatureData, StringData, "T", t);
  if (retVal < 0) return -1;
  vector<double> v;
  retVal = getVec(DoubleFeatureData, StringData, "V", v);
  if (retVal < 0) return -1;
  vector<double> stimstart;
  retVal = getVec(DoubleFeatureData, StringData, "stim_start", stimstart);
  if (retVal < 0) return -1;
  vector<double> stimend;
  retVal = getVec(DoubleFeatureData, StringData, "stim_end", stimend);
  if (retVal < 0) return -1;
  vector<int> apendi;
  retVal = getVec(IntFeatureData, StringData, "AP_end_indices", apendi);
  if (retVal < 0) return -1;
  vector<int> apbi;
  retVal = getVec(IntFeatureData, StringData, "AP_begin_indices", apbi);
  if (retVal < 0) return -1;

  vector<double> dep_base;
  retVal = __depolarized_base(t, v, stimstart[0], stimend[0], apbi, apendi,
                              dep_base);
  if (retVal >= 0) {
    setVec(DoubleFeatureData, StringData, "depolarized_base", dep_base);
  }
  return retVal;
}
