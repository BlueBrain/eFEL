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
#include <algorithm>
#include <iterator>
#include <math.h>
#include <string>

int LibV1::interpolate(mapStr2intVec& IntFeatureData,
                       mapStr2doubleVec& DoubleFeatureData,
                       mapStr2Str& StringData) {
  int retVal, nSize;
  retVal =
      CheckInIntmap(IntFeatureData, StringData, string("interpolate"), nSize);
  if (retVal)
    return nSize;
  else {
    vector<double> V, T, VIntrpol, TIntrpol, InterpStepVec;
    vector<int> intrpolte;
    double InterpStep;
    // getDoubleVec takes care of stimulus suffix
    retVal = getDoubleVec(DoubleFeatureData, StringData, string("V"), V);
    if (retVal <= 0) return -1;
    retVal = getDoubleVec(DoubleFeatureData, StringData, string("T"), T);
    if (retVal <= 0) return -1;
    // interp_step is a stimulus independent parameter
    retVal =
        getDoubleParam(DoubleFeatureData, string("interp_step"), InterpStepVec);
    if (retVal <= 0)
      InterpStep = 0.1;
    else
      InterpStep = InterpStepVec[0];

    LinearInterpolation(InterpStep, T, V, TIntrpol, VIntrpol);

    setDoubleVec(DoubleFeatureData, StringData, "V", VIntrpol);
    setDoubleVec(DoubleFeatureData, StringData, "T", TIntrpol);
    setIntVec(IntFeatureData, StringData, "interpolate", intrpolte);
    return retVal;
  }
  return -1;
}

int LibV1::__peak_indices(double dThreshold, vector<double>& V,
                          vector<int>& PeakIndex) {
  vector<int> upVec, dnVec;
  double dtmp;
  int itmp;

  for (unsigned i = 1; i < V.size(); i++) {
    if (V[i] > dThreshold && V[i - 1] < dThreshold) {
      upVec.push_back(i);
    } else if (V[i] < dThreshold && V[i - 1] > dThreshold) {
      dnVec.push_back(i);
    }
  }
  if (dnVec.size() == 0) {
    GErrorStr =
        GErrorStr +
        "\nVoltage never goes below or above threshold in spike detection.\n";
    return 0;
  }

  if (dnVec.size() != upVec.size()) {
    GErrorStr =
        GErrorStr + "\nVoltage never goes below threshold after last spike.\n";
    return 0;
  }

  PeakIndex.clear();
  int j = 0;
  for (unsigned i = 0; i < upVec.size(); i++) {
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
int LibV1::peak_indices(mapStr2intVec& IntFeatureData,
                        mapStr2doubleVec& DoubleFeatureData,
                        mapStr2Str& StringData) {
  // printf("\n  LibV1  This is inside peak_indices()");
  int retVal, nSize;
  retVal =
      CheckInIntmap(IntFeatureData, StringData, string("peak_indices"), nSize);
  if (retVal)
    return nSize;
  else {

    vector<int> PeakIndex;
    vector<double> v, Th;
    retVal = getDoubleVec(DoubleFeatureData, StringData, string("V"), v);
    if (retVal <= 0) return -1;
    retVal = getDoubleParam(DoubleFeatureData, string("Threshold"), Th);
    if (retVal <= 0) return -1;
    int retval = __peak_indices(Th[0], v, PeakIndex);
    if (retval >= 0)
      setIntVec(IntFeatureData, StringData, "peak_indices", PeakIndex);
    return retval;
  }
}

// *** Spikecount ***
int LibV1::Spikecount(mapStr2intVec& IntFeatureData,
                      mapStr2doubleVec& DoubleFeatureData,
                      mapStr2Str& StringData) {
  int retval;
  int nsize;
  retval =
      CheckInIntmap(IntFeatureData, StringData, string("Spikecount"), nsize);
  if (retval) {
    return nsize;
  }
  vector<int> peakindices;
  retval = getIntVec(IntFeatureData, StringData, string("peak_indices"),
                     peakindices);
  if (retval <= 0) return -1;
  vector<int> spikecount(1, peakindices.size());
  if (retval >= 0) {
    setIntVec(IntFeatureData, StringData, "Spikecount", spikecount);
  }
  return retval;
}
// end of Spikecount

// ISI starting at the 2nd spike, the first ISI is the feature doublet_ISI
int LibV1::__ISI_values(vector<double>& V, vector<int>& PeakIndex,
                        vector<double>& PeakVoltage) {
  //    printf("\n  LibV1  This is inside __ISI_values()");
  return 1;
}
int LibV1::ISI_values(mapStr2intVec& IntFeatureData,
                      mapStr2doubleVec& DoubleFeatureData,
                      mapStr2Str& StringData) {
  //  printf("\n  LibV1  This is inside ISI_values()");
  int retVal, nSize;
  retVal = CheckInDoublemap(DoubleFeatureData, StringData, string("ISI_values"),
                            nSize);
  if (retVal)
    return nSize;
  else {
    vector<double> VecISI, pvTime;
    int retVal = 0;
    retVal = getDoubleVec(DoubleFeatureData, StringData, string("peak_time"),
                          pvTime);
    if (retVal < 3) {
      GErrorStr = GErrorStr +
                  "\n Three spikes required for calculation of ISI_values.\n";
      return -1;
    }
    for (unsigned i = 2; i < pvTime.size(); i++)
      VecISI.push_back(pvTime[i] - pvTime[i - 1]);
    setDoubleVec(DoubleFeatureData, StringData, "ISI_values", VecISI);
    return VecISI.size();
  }
  return -1;
}

// *** ISI_CV ***
// the coefficient of variation of the ISI
int LibV1::__ISI_CV(const vector<double>& isivalues, vector<double>& isicv) {
  // mean
  double isi_mean = 0.;
  for (unsigned i = 0; i < isivalues.size(); i++) {
    isi_mean += isivalues[i];
  }
  isi_mean /= isivalues.size();

  // sigma^2
  double variance = 0.;
  for (unsigned i = 0; i < isivalues.size(); i++) {
    double dev = isivalues[i] - isi_mean;
    variance += dev * dev;
  }
  // variation coefficient cv = sigma / mean
  isicv.push_back(sqrt(variance / (isivalues.size() - 1)) / isi_mean);
  return isicv.size();
}
int LibV1::ISI_CV(mapStr2intVec& IntFeatureData,
                  mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData) {
  int retval;
  int nsize;
  retval =
      CheckInDoublemap(DoubleFeatureData, StringData, string("ISI_CV"), nsize);
  if (retval) {
    return nsize;
  }
  vector<double> isivalues;
  retval = getDoubleVec(DoubleFeatureData, StringData, string("ISI_values"),
                        isivalues);
  if (retval < 2) return -1;
  vector<double> isicv;
  retval = __ISI_CV(isivalues, isicv);
  if (retval >= 0) {
    setDoubleVec(DoubleFeatureData, StringData, "ISI_CV", isicv);
  }
  return retval;
}
// end of ISI_CV

int LibV1::peak_voltage(mapStr2intVec& IntFeatureData,
                        mapStr2doubleVec& DoubleFeatureData,
                        mapStr2Str& StringData) {
  // printf("\n  LibV1  Inside PeakVoltage ... LibV1");
  int retVal, nSize;
  retVal = CheckInDoublemap(DoubleFeatureData, StringData,
                            string("peak_voltage"), nSize);
  if (retVal)
    return nSize;
  else {
    // vector<int> PeakI = getIntVec(IntFeatureData, StringData,
    // string("peak_indices"));
    vector<int> PeakI;
    vector<double> V, peakV;
    int retVal;
    retVal =
        getIntVec(IntFeatureData, StringData, string("peak_indices"), PeakI);
    {
      if (retVal <= 0) return -1;
    };
    retVal = getDoubleVec(DoubleFeatureData, StringData, string("V"), V);
    {
      if (retVal <= 0) return -1;
    };
    for (unsigned i = 0; i < PeakI.size(); i++) {
      peakV.push_back(V[PeakI[i]]);
    }
    setDoubleVec(DoubleFeatureData, StringData, "peak_voltage", peakV);
    return peakV.size();
  }
  return 1;
}

int LibV1::firing_rate(mapStr2intVec& IntFeatureData,
                       mapStr2doubleVec& DoubleFeatureData,
                       mapStr2Str& StringData) {
  // printf("\n  LibV1  This is inside firing_rate()");
  int retVal, nSize;
  retVal = CheckInDoublemap(DoubleFeatureData, StringData,
                            string("mean_frequency"), nSize);
  if (retVal)
    return nSize;
  else {
    vector<double> stimStart, stimEnd, peakVTime, firing_rate;
    double lastAPTime;
    int retVal;
    retVal = getDoubleVec(DoubleFeatureData, StringData, string("peak_time"),
                          peakVTime);
    {
      if (retVal <= 0) return -1;
    };
    retVal =
        getDoubleVec(DoubleFeatureData, StringData, "stim_start", stimStart);
    {
      if (retVal <= 0) return -1;
    };
    retVal = getDoubleVec(DoubleFeatureData, StringData, "stim_end", stimEnd);
    {
      if (retVal <= 0) return -1;
    };
    int nCount = 0;
    for (unsigned i = 0; i < peakVTime.size(); i++) {
      if ((peakVTime[i] >= stimStart[0]) && (peakVTime[i] <= stimEnd[0])) {
        lastAPTime = peakVTime[i];
        nCount++;
      }
    }
    firing_rate.push_back(nCount * 1000 / (lastAPTime - stimStart[0]));
    setDoubleVec(DoubleFeatureData, StringData, "mean_frequency", firing_rate);
    return firing_rate.size();
  }
  return -1;
}

int LibV1::peak_time(mapStr2intVec& IntFeatureData,
                     mapStr2doubleVec& DoubleFeatureData,
                     mapStr2Str& StringData) {
  // printf("\n  LibV1  This is inside peak_time()");
  int retVal, nSize;
  retVal = CheckInDoublemap(DoubleFeatureData, StringData, string("peak_time"),
                            nSize);
  if (retVal)
    return nSize;
  else {
    vector<int> PeakI;
    vector<double> T, pvTime;
    int retVal = 0;
    retVal =
        getIntVec(IntFeatureData, StringData, string("peak_indices"), PeakI);
    if (retVal < 0) return -1;
    retVal = getDoubleVec(DoubleFeatureData, StringData, string("T"), T);
    if (retVal < 0) return -1;
    for (unsigned i = 0; i < PeakI.size(); i++) {
      pvTime.push_back(T[PeakI[i]]);
    }
    setDoubleVec(DoubleFeatureData, StringData, "peak_time", pvTime);
    return pvTime.size();
  }
  return -1;
}

// time from stimulus start to first threshold crossing
int LibV1::first_spike_time(mapStr2intVec& IntFeatureData,
                            mapStr2doubleVec& DoubleFeatureData,
                            mapStr2Str& StringData) {
  int retVal, nSize;
  retVal = CheckInDoublemap(DoubleFeatureData, StringData,
                            "time_to_first_spike", nSize);
  if (retVal)
    return nSize;
  else {
    vector<double> first_spike;
    vector<double> peaktime;
    vector<double> stimstart;
    int retVal;
    retVal = getDoubleVec(DoubleFeatureData, StringData, "peak_time", peaktime);
    if (retVal < 1) {
      GErrorStr += "\n One spike required for time_to_first_spike.\n";
      return -1;
    }
    retVal =
        getDoubleVec(DoubleFeatureData, StringData, "stim_start", stimstart);
    if (retVal <= 0) return -1;
    first_spike.push_back(peaktime[0] - stimstart[0]);
    setDoubleVec(DoubleFeatureData, StringData, "time_to_first_spike",
                 first_spike);
    return first_spike.size();
  }
  return -1;
}

// min_AHP_indices
// find the minimum between two spikes,
// and the minimum between the last spike and the time the stimulus ends
int LibV1::min_AHP_indices(mapStr2intVec& IntFeatureData,
                           mapStr2doubleVec& DoubleFeatureData,
                           mapStr2Str& StringData) {
  int retVal, nSize;
  retVal = CheckInIntmap(IntFeatureData, StringData, string("min_AHP_indices"),
                         nSize);
  if (retVal)
    return nSize;
  else {
    vector<int> peak_indices_plus;
    vector<int> min_ahp_indices;
    vector<double> v;
    vector<double> min_ahp_values;
    vector<double> stim_end;
    vector<double> t;
    retVal = getDoubleVec(DoubleFeatureData, StringData, "V", v);
    if (retVal <= 0) return -1;
    retVal = getIntVec(IntFeatureData, StringData, "peak_indices",
                       peak_indices_plus);
    if (retVal < 1) {
      GErrorStr +=
          "\n At least one spike required for calculation of "
          "min_AHP_indices.\n";
      return -1;
    }
    retVal = getDoubleVec(DoubleFeatureData, StringData, "stim_end", stim_end);
    retVal = getDoubleVec(DoubleFeatureData, StringData, "T", t);

    int end_index = distance(
        t.begin(), find_if(t.begin(), t.end(),
                           bind2nd(greater_equal<double>(), stim_end[0])));
    // if the last spike happens to be close to the end of the stimulus
    // there will not be a proper AHP, this case is not properly dealt with here
    if (end_index > peak_indices_plus.back() + 5) {
      peak_indices_plus.push_back(end_index);
    }
    for (unsigned i = 0; i < peak_indices_plus.size() - 1; i++) {
      int ahpindex = distance(
          v.begin(), min_element(v.begin() + peak_indices_plus[i],
                                 v.begin() + peak_indices_plus[i + 1]));
      min_ahp_indices.push_back(ahpindex);
      min_ahp_values.push_back(v[ahpindex]);
    }
    setIntVec(IntFeatureData, StringData, "min_AHP_indices", min_ahp_indices);
    setDoubleVec(DoubleFeatureData, StringData, "min_AHP_values",
                 min_ahp_values);
    return min_ahp_indices.size();
  }
  return -1;
}

int LibV1::min_AHP_values(mapStr2intVec& IntFeatureData,
                          mapStr2doubleVec& DoubleFeatureData,
                          mapStr2Str& StringData) {
  int retVal, nSize;
  retVal = CheckInDoublemap(DoubleFeatureData, StringData,
                            string("min_AHP_values"), nSize);
  if (retVal) return nSize;
  return -1;
}

int LibV1::AP_height(mapStr2intVec& IntFeatureData,
                     mapStr2doubleVec& DoubleFeatureData,
                     mapStr2Str& StringData) {
  int retVal, nSize;
  retVal = CheckInDoublemap(DoubleFeatureData, StringData, string("AP_height"),
                            nSize);
  if (retVal)
    return nSize;
  else {
    vector<double> vPeak;
    retVal = getDoubleVec(DoubleFeatureData, StringData, string("peak_voltage"),
                          vPeak);
    if (retVal <= 0) return -1;
    setDoubleVec(DoubleFeatureData, StringData, "AP_height", vPeak);
    return vPeak.size();
  }
  return -1;
}

// spike amplitude: peak_voltage - v[AP_begin_indices]
int LibV1::AP_amplitude(mapStr2intVec& IntFeatureData,
                        mapStr2doubleVec& DoubleFeatureData,
                        mapStr2Str& StringData) {
  int retVal, nSize;
  retVal =
      CheckInDoublemap(DoubleFeatureData, StringData, "AP_amplitude", nSize);
  if (retVal > 0) {
    return nSize;
  } else {
    vector<double> peakvoltage;
    vector<int> apbeginindices;
    vector<double> v;
    retVal = getDoubleVec(DoubleFeatureData, StringData, "V", v);
    if (retVal <= 0) {
      return -1;
    }
    retVal = getDoubleVec(DoubleFeatureData, StringData, "peak_voltage",
                          peakvoltage);
    if (retVal <= 0) {
      GErrorStr += "Error calculating peak_voltage for AP_amplitude";
      return -1;
    }
    retVal = getIntVec(IntFeatureData, StringData, "AP_begin_indices",
                       apbeginindices);
    if (retVal <= 0) {
      GErrorStr += "Error calculating AP_begin_indicies for AP_amplitude";
      return -1;
    }

    if (peakvoltage.size() > apbeginindices.size()) {
      GErrorStr += "More peakvoltage entries than AP begin entries";
      return -1;
    }

    vector<double> apamplitude;
    apamplitude.resize(peakvoltage.size());
    for (unsigned i = 0; i < apamplitude.size(); i++) {
      apamplitude[i] = peakvoltage[i] - v[apbeginindices[i]];
    }
    setDoubleVec(DoubleFeatureData, StringData, "AP_amplitude", apamplitude);
    return apamplitude.size();
  }
  return -1;
}

// AHP_depth_abs
// naming conflict here AHP_depth_abs does the same as min_AHP_values.
// In my opinion the should not be a feature called 'AHP_depth_abs',
// use min_AHP_values instead.
// A more interesting feature would be 'AHP_depth' anyways, which calculates the
// depth of the AHP relative to the voltage base
int LibV1::AHP_depth_abs(mapStr2intVec& IntFeatureData,
                         mapStr2doubleVec& DoubleFeatureData,
                         mapStr2Str& StringData) {
  int retVal, nSize;
  retVal = CheckInDoublemap(DoubleFeatureData, StringData,
                            string("AHP_depth_abs"), nSize);
  if (retVal)
    return nSize;
  else {
    vector<double> vAHP;
    retVal = getDoubleVec(DoubleFeatureData, StringData,
                          string("min_AHP_values"), vAHP);
    if (retVal <= 0) return -1;
    setDoubleVec(DoubleFeatureData, StringData, "AHP_depth_abs", vAHP);
    return vAHP.size();
  }
  return -1;
}

// *** AHP_depth_abs_slow ***
// same as AHP_depth_abs but the minimum search starts 5 ms after the spike,
// first ISI is ignored
int LibV1::__AHP_depth_abs_slow_indices(const vector<double>& t,
                                        const vector<double>& v,
                                        const vector<int>& peakindices,
                                        vector<int>& adas_indices) {
  adas_indices.resize(peakindices.size() - 2);
  for (unsigned i = 0; i < adas_indices.size(); i++) {
    // start 5 ms after last spike
    double t_start = t[peakindices[i + 1]] + 5;
    adas_indices[i] = distance(
        v.begin(),
        min_element(
            v.begin() +
                distance(t.begin(),
                         find_if(t.begin() + peakindices[i + 1],
                                 t.begin() + peakindices[i + 2],
                                 bind2nd(greater_equal<double>(), t_start))),
            v.begin() + peakindices[i + 2]));
  }
  return adas_indices.size();
}

int LibV1::AHP_depth_abs_slow(mapStr2intVec& IntFeatureData,
                              mapStr2doubleVec& DoubleFeatureData,
                              mapStr2Str& StringData) {
  int retval;
  int nsize;
  retval = CheckInDoublemap(DoubleFeatureData, StringData,
                            string("AHP_depth_abs_slow"), nsize);
  if (retval) {
    return nsize;
  }
  vector<double> t;
  retval = getDoubleVec(DoubleFeatureData, StringData, "T", t);
  if (retval < 0) return -1;
  vector<double> v;
  retval = getDoubleVec(DoubleFeatureData, StringData, "V", v);
  if (retval < 0) return -1;
  vector<int> peakindices;
  retval = getIntVec(IntFeatureData, StringData, "peak_indices", peakindices);
  if (retval < 3) {
    GErrorStr +=
        "\n At least 3 spikes needed for AHP_depth_abs_slow and "
        "AHP_slow_time.\n";
    return -1;
  }
  vector<int> adas_indices;
  retval = __AHP_depth_abs_slow_indices(t, v, peakindices, adas_indices);
  vector<double> ahpdepthabsslow(adas_indices.size());
  vector<double> ahpslowtime(adas_indices.size());
  for (unsigned i = 0; i < adas_indices.size(); i++) {
    ahpdepthabsslow[i] = v[adas_indices[i]];
    ahpslowtime[i] = (t[adas_indices[i]] - t[peakindices[i + 1]]) /
                     (t[peakindices[i + 2]] - t[peakindices[i + 1]]);
  }
  if (retval >= 0) {
    setDoubleVec(DoubleFeatureData, StringData, "AHP_depth_abs_slow",
                 ahpdepthabsslow);
    setDoubleVec(DoubleFeatureData, StringData, "AHP_slow_time", ahpslowtime);
  }
  return retval;
}
// end of AHP_depth_abs_slow

int LibV1::AHP_slow_time(mapStr2intVec& IntFeatureData,
                         mapStr2doubleVec& DoubleFeatureData,
                         mapStr2Str& StringData) {
  int retVal, nSize;
  retVal =
      CheckInDoublemap(DoubleFeatureData, StringData, "AHP_slow_time", nSize);
  if (retVal) return nSize;
  return -1;
}

int LibV1::__rest_voltage_value(vector<double>& V, vector<int>& PeakIndex,
                                vector<double>& PeakVoltage) {
  // printf("\n  LibV1  This is inside  __rest_voltage_value()");
  return 1;
}
int LibV1::rest_voltage_value(mapStr2intVec& IntFeatureData,
                              mapStr2doubleVec& DoubleFeatureData,
                              mapStr2Str& StringData) {
  int retVal, nSize;
  retVal = CheckInDoublemap(DoubleFeatureData, StringData,
                            string("voltage_base"), nSize);
  if (retVal)
    return nSize;
  else {
    vector<double> v, t, stimStart, vRest;
    double startTime, endTime;
    retVal = getDoubleVec(DoubleFeatureData, StringData, string("V"), v);
    if (retVal < 0) return -1;
    retVal = getDoubleVec(DoubleFeatureData, StringData, string("T"), t);
    if (retVal < 0) return -1;
    retVal =
        getDoubleVec(DoubleFeatureData, StringData, "stim_start", stimStart);
    if (retVal < 0) return -1;
    startTime = stimStart[0] * .25;  // It is 25% from start (0), so if stimulus
                                     // starts at 100ms then StartTime will be
                                     // 25mS
    endTime = stimStart[0] *
              .75;  // as in above case end time will be 25% less than startTime
    int nCount = 0;
    double vSum = 0;
    // calculte the mean of voltage between startTime and endTime
    for (unsigned i = 0; i < t.size(); i++) {
      if (t[i] >= startTime) {
        vSum = vSum + v[i];
        nCount++;
      }
      if (t[i] > endTime) break;
    }
    vRest.push_back(vSum / nCount);
    setDoubleVec(DoubleFeatureData, StringData, "voltage_base", vRest);
    return 1;
  }
  return -1;
}

int LibV1::__burst_ISI_indices(double BurstFactor, vector<int>& PeakIndex,
                               vector<double>& ISIValues,
                               vector<int>& BurstIndex) {
  vector<double> ISIpcopy;
  vector<double>::iterator it1, it2;
  int n, count = -1;
  double dMedian;

  for (unsigned i = 1; i < (ISIValues.size() - 1); i++) {
    ISIpcopy.clear();
    for (unsigned j = count + 1; j < i; j++) ISIpcopy.push_back(ISIValues[j]);
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
  int retVal, nSize;
  retVal = CheckInDoublemap(DoubleFeatureData, StringData,
                            string("burst_ISI_indices"), nSize);
  if (retVal)
    return nSize;
  else {
    vector<int> PeakIndex, BurstIndex;
    vector<double> ISIValues, tVec;
    double BurstFactor = 0;
    retVal = getIntVec(IntFeatureData, StringData, string("peak_indices"),
                       PeakIndex);
    if (retVal < 0) return -1;
    if (PeakIndex.size() < 5) {
      GErrorStr =
          GErrorStr +
          "\nError: More than 5 spike is needed for burst calculation.\n";
      return -1;
    }
    retVal = getDoubleVec(DoubleFeatureData, StringData, string("ISI_values"),
                          ISIValues);
    if (retVal < 0) return -1;
    retVal = getDoubleParam(DoubleFeatureData, string("burst_factor"), tVec);
    {
      if (retVal < 0)
        BurstFactor = 2;
      else
        BurstFactor = tVec[0];
    }
    retVal = __burst_ISI_indices(BurstFactor, PeakIndex, ISIValues, BurstIndex);
    if (retVal >= 0) {
      setIntVec(IntFeatureData, StringData, "burst_ISI_indices", BurstIndex);
    }
    return retVal;
  }
  return -1;
}

int LibV1::__burst_mean_freq(vector<double>& PVTime, vector<int>& BurstIndex,
                             vector<double>& BurstMeanFreq) {
  vector<double> tmpVec;
  BurstIndex.insert(BurstIndex.begin(), 0);
  for (unsigned i = 0; i < BurstIndex.size(); i++) tmpVec.push_back(0);
  double span;
  unsigned i;
  for (i = 0; i < BurstIndex.size() - 1; i++) {
    if (BurstIndex[i + 1] - BurstIndex[i] == 1) {
      tmpVec.push_back(0);
    } else {
      span = PVTime[BurstIndex[i + 1] - 1] - PVTime[BurstIndex[i]];
      tmpVec.push_back((BurstIndex[i + 1] - BurstIndex[i] + 1) * 1000 / span);
    }
  }
  // NOW LAST BURST
  span = PVTime[PVTime.size() - 1] - PVTime[BurstIndex[i]];
  tmpVec.push_back((PVTime.size() - 1 - BurstIndex[i] + 1) * 1000 / span);

  for (i = 0; i < tmpVec.size(); i++) {
    if (tmpVec[i] != 0) BurstMeanFreq.push_back(tmpVec[i]);
  }
  return BurstMeanFreq.size();
}

int LibV1::burst_mean_freq(mapStr2intVec& IntFeatureData,
                           mapStr2doubleVec& DoubleFeatureData,
                           mapStr2Str& StringData) {
  int retVal, nSize;
  retVal = CheckInDoublemap(DoubleFeatureData, StringData,
                            string("burst_mean_freq"), nSize);
  if (retVal)
    return nSize;
  else {
    vector<int> BurstIndex;
    vector<double> BurstMeanFreq, PVTime;
    retVal = getDoubleVec(DoubleFeatureData, StringData, string("peak_time"),
                          PVTime);
    if (retVal < 0) return -1;
    retVal = getIntVec(IntFeatureData, StringData, string("burst_ISI_indices"),
                       BurstIndex);
    if (retVal < 0) return -1;
    retVal = __burst_mean_freq(PVTime, BurstIndex, BurstMeanFreq);
    if (retVal >= 0) {
      setDoubleVec(DoubleFeatureData, StringData, "burst_mean_freq",
                   BurstMeanFreq);
    }
    return retVal;
  }
  return -1;
}

int LibV1::burst_number(mapStr2intVec& IntFeatureData,
                        mapStr2doubleVec& DoubleFeatureData,
                        mapStr2Str& StringData) {
  int retVal, nSize;
  retVal =
      CheckInIntmap(IntFeatureData, StringData, string("burst_number"), nSize);
  if (retVal)
    return nSize;
  else {
    vector<double> BurstMeanFreq;
    vector<int> BurstNum;
    retVal = getDoubleVec(DoubleFeatureData, StringData,
                          string("burst_mean_freq"), BurstMeanFreq);
    if (retVal < 0) return -1;
    if (retVal >= 0) {
      BurstNum.push_back(BurstMeanFreq.size());
      setIntVec(IntFeatureData, StringData, "burst_number", BurstNum);
      return (BurstNum.size());
    }
    return retVal;
  }
  return -1;
}

int LibV1::__interburst_voltage(vector<int>& BurstIndex, vector<int>& PeakIndex,
                                vector<double>& T, vector<double>& V,
                                vector<double>& IBV) {
  if (BurstIndex.size() < 2) return 0;
  int j, pIndex, tsIndex, teIndex, cnt;
  double tStart, tEnd, vTotal = 0;
  for (unsigned i = 0; i < BurstIndex.size(); i++) {
    pIndex = BurstIndex[i] - 1;
    tsIndex = PeakIndex[pIndex];
    tStart = T[tsIndex] + 5;  // 5Millisecond after
    pIndex = BurstIndex[i];
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
  int retVal, nSize;
  retVal = CheckInDoublemap(DoubleFeatureData, StringData,
                            string("interburst_voltage"), nSize);
  if (retVal)
    return nSize;
  else {
    vector<int> BurstIndex, PeakIndex;
    vector<double> V, T, IBV;
    retVal = getIntVec(IntFeatureData, StringData, string("peak_indices"),
                       PeakIndex);
    if (retVal < 0) return -1;
    retVal = getDoubleVec(DoubleFeatureData, StringData, string("T"), T);
    if (retVal < 0) return -1;
    retVal = getIntVec(IntFeatureData, StringData, string("burst_ISI_indices"),
                       BurstIndex);
    if (retVal < 0) return -1;
    retVal = getDoubleVec(DoubleFeatureData, StringData, string("V"), V);
    if (retVal < 0) return -1;
    retVal = __interburst_voltage(BurstIndex, PeakIndex, T, V, IBV);
    if (retVal >= 0) {
      setDoubleVec(DoubleFeatureData, StringData, "interburst_voltage", IBV);
    }
    return retVal;
  }
  return -1;
}

int LibV1::__adaptation_index(double spikeSkipf, int maxnSpike,
                              double StimStart, double StimEnd, double Offset,
                              vector<double>& peakVTime,
                              vector<double>& adaptation_index) {
  list<double> SpikeTime;
  vector<double> ISI;
  // Select spike time between given time scale (stim_start and stim_end )
  // considet Offset also if it is given as input
  for (unsigned i = 0; i < peakVTime.size(); i++) {
    if ((peakVTime[i] >= (StimStart - Offset)) &&
        (peakVTime[i] <= (StimEnd + Offset))) {
      SpikeTime.push_back(peakVTime[i]);
    }
  }
  // Remove n spikes given by spike_skipf or max_spike_skip
  unsigned spikeToRemove = (unsigned)((SpikeTime.size() * spikeSkipf) + 0.5);
  // spike To remove is minimum of spike_skipf or max_spike_skip
  if (maxnSpike < spikeToRemove) spikeToRemove = maxnSpike;

  // Remove spikeToRemove spike from SpikeTime list
  for (unsigned i = 0; i < spikeToRemove; i++) {
    SpikeTime.pop_front();
  }

  // Adaptation index can not be calculated if nAPVec <4 or no of ISI is <3
  if (SpikeTime.size() < 4) {
    GErrorStr = GErrorStr +
                "\nMinimum 4 spike needed for feature [adaptation_index].\n";
    return -1;
  }

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
  for (unsigned i = 1; i < ISI.size(); i++) {
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
  int retVal, nSize;
  retVal = CheckInDoublemap(DoubleFeatureData, StringData,
                            string("adaptation_index"), nSize);
  if (retVal)
    return nSize;
  else {
    vector<double> peakVTime, stimStart, stimEnd, OffSetVec, spikeSkipf,
        adaptation_index;
    vector<int> maxnSpike;
    int retVal;
    double Offset;
    retVal = getDoubleVec(DoubleFeatureData, StringData, string("peak_time"),
                          peakVTime);
    {
      if (retVal < 0) return -1;
    };
    retVal =
        getDoubleVec(DoubleFeatureData, StringData, "stim_start", stimStart);
    {
      if (retVal < 0) return -1;
    };
    retVal = getDoubleVec(DoubleFeatureData, StringData, "stim_end", stimEnd);
    {
      if (retVal < 0) return -1;
    };
    retVal =
        getDoubleParam(DoubleFeatureData, string("spike_skipf"), spikeSkipf);
    {
      if (retVal < 0) return -1;
    };
    // spikeSkipf is a fraction hence value should lie between >=0 and <1. [0 1)
    if ((spikeSkipf[0] < 0) || (spikeSkipf[0] >= 1)) {
      GErrorStr = GErrorStr + "\nspike_skipf should lie between [0 1).\n";
      return -1;
    }
    retVal = getIntParam(IntFeatureData, string("max_spike_skip"), maxnSpike);
    {
      if (retVal < 0) return -1;
    };
    retVal = getDoubleParam(DoubleFeatureData, string("offset"), OffSetVec);
    if (retVal < 0)
      Offset = 0;
    else
      Offset = OffSetVec[0];
    retVal =
        __adaptation_index(spikeSkipf[0], maxnSpike[0], stimStart[0],
                           stimEnd[0], Offset, peakVTime, adaptation_index);
    if (retVal >= 0)
      setDoubleVec(DoubleFeatureData, StringData, "adaptation_index",
                   adaptation_index);
    return retVal;
  }
  return -1;
}

// *** adaptation_index2 ***
// as adaptation_index, but start at the second ISI instead of the round(N *
// spikeskipf)
int LibV1::__adaptation_index2(double StimStart, double StimEnd, double Offset,
                               const vector<double>& peakVTime,
                               vector<double>& adaptation_index) {
  list<double> SpikeTime;
  vector<double> ISI;
  // Select spike time between given time scale (stim_start and stim_end )
  // considet Offset also if it is given as input
  for (unsigned i = 0; i < peakVTime.size(); i++) {
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
  for (unsigned i = 1; i < ISI.size(); i++) {
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
  int retval;
  int nsize;
  retval = CheckInDoublemap(DoubleFeatureData, StringData,
                            string("adaptation_index2"), nsize);
  if (retval) {
    return nsize;
  }
  vector<double> peakvoltagetime;
  retval = getDoubleVec(DoubleFeatureData, StringData, string("peak_time"),
                        peakvoltagetime);
  if (retval < 4) {
    GErrorStr += "\n At least 4 spikes needed for adaptation_index2.\n";
    return -1;
  }
  vector<double> stimStart;
  retval = getDoubleVec(DoubleFeatureData, StringData, "stim_start", stimStart);
  {
    if (retval < 0) return -1;
  };
  vector<double> stimEnd;
  retval = getDoubleVec(DoubleFeatureData, StringData, "stim_end", stimEnd);
  {
    if (retval < 0) return -1;
  };
  vector<double> OffSetVec;
  double Offset;
  retval = getDoubleParam(DoubleFeatureData, string("offset"), OffSetVec);
  if (retval < 0)
    Offset = 0;
  else
    Offset = OffSetVec[0];
  vector<double> adaptationindex2;
  retval = __adaptation_index2(stimStart[0], stimEnd[0], Offset,
                               peakvoltagetime, adaptationindex2);
  if (retval >= 0) {
    setDoubleVec(DoubleFeatureData, StringData, "adaptation_index2",
                 adaptationindex2);
  }
  return retval;
}
// end of adaptation_index2

int LibV1::trace_check(mapStr2intVec& IntFeatureData,
                       mapStr2doubleVec& DoubleFeatureData,
                       mapStr2Str& StringData) {
  int retval;
  int size;
  retval = CheckInIntmap(IntFeatureData, StringData, "trace_check", size);
  if (retval) {
    return size;
  } else {
    vector<double> peak_time;
    vector<double> stim_start;
    vector<double> stim_end;
    vector<int> tc;
    retval =
        getDoubleVec(DoubleFeatureData, StringData, "peak_time", peak_time);
    if (retval < 0) return -1;
    retval =
        getDoubleVec(DoubleFeatureData, StringData, "stim_start", stim_start);
    if (retval < 0) return -1;
    retval = getDoubleVec(DoubleFeatureData, StringData, "stim_end", stim_end);
    if (retval < 0) return -1;
    bool sane = true;
    for (unsigned i = 0; i < peak_time.size(); i++) {
      if (peak_time[i] < stim_start[0] || peak_time[i] > stim_end[0] * 1.05) {
        sane = false;
        break;
      }
    }
    if (sane) {
      tc.push_back(0);
      setIntVec(IntFeatureData, StringData, "trace_check", tc);
      return tc.size();
    } else {
      GErrorStr +=
          "Trace sanity check failed, there were spike outside the stimulus "
          "interval.\n";
      return -1;
    }
  }
}

// To find spike width using Central difference derivative vec1[i] =
// ((vec[i+1]+vec[i-1])/2)/dx  and half width is between
// MinAHP and APThreshold
int LibV1::__spike_width2(vector<double>& t, vector<double>& V,
                          vector<int>& PeakIndex, vector<int>& minAHPIndex,
                          vector<double>& spike_width2) {
  vector<double> v, dv1, dv2;
  double dx = t[1] - t[0];
  double VoltThreshold, VoltMax, HalfV, T0, V0, V1, fraction, TStart, TEnd;
  int j, index;
  for (unsigned i = 0; i < minAHPIndex.size() && i < PeakIndex.size() - 1; i++) {
    v.clear();
    dv1.clear();
    dv2.clear();
    for (unsigned j = minAHPIndex[i]; j <= PeakIndex[i + 1]; j++) {
      v.push_back(V[j]);
    }
    index = v.size();  // tbr
    getCentralDifferenceDerivative(dx, v, dv1);
    getCentralDifferenceDerivative(dx, dv1, dv2);
    double dMax = dv2[0];
    index = 0;
    for (j = 1; j < dv2.size(); j++) {
      if (dMax <= dv2[j]) {
        dMax = dv2[j];
        index = j;
      }
    }
    // Take voltage at point where double derivative is maximum
    index = minAHPIndex[i] + index;
    VoltThreshold = V[index];
    VoltMax = V[PeakIndex[i + 1]];
    HalfV = (VoltMax + VoltThreshold) / 2;

    // Find voltage where it crosses HalfV in the rising phase of action
    // potential
    for (j = 0; j < v.size(); j++) {
      if (v[j] > HalfV) break;  // point is found where  v is crossing HalfV
    }
    index = minAHPIndex[i] + j;
    // index is the index where v crossed HalfV now use linear interpolation to
    // find actual time at HalfV
    T0 = t[index - 1];
    V0 = V[index - 1];
    V1 = V[index];
    fraction = (HalfV - V0) / (V1 - V0);
    TStart = T0 + (fraction * dx);

    // Find voltage where it crosses HalfV in the falling phase of the action
    // potential
    for (j = PeakIndex[i + 1]; j < V.size(); j++) {
      if (V[j] < HalfV) break;
    }
    if (j == V.size()) {
      GErrorStr = GErrorStr + "\nFalling phase of last spike is missing.\n";
      return -1;
    }
    T0 = t[j - 1];
    V0 = V[j - 1];
    V1 = V[j];
    fraction = (HalfV - V0) / (V1 - V0);
    TEnd = T0 + (fraction * dx);
    spike_width2.push_back(TEnd - TStart);
  }
  return spike_width2.size();
}

int LibV1::spike_width2(mapStr2intVec& IntFeatureData,
                        mapStr2doubleVec& DoubleFeatureData,
                        mapStr2Str& StringData) {
  int retVal, nSize;
  retVal = CheckInDoublemap(DoubleFeatureData, StringData,
                            string("spike_width2"), nSize);
  if (retVal)
    return nSize;
  else {
    vector<int> PeakIndex, minAHPIndex;
    vector<double> V, t, dv1, dv2, spike_width2;
    retVal = getDoubleVec(DoubleFeatureData, StringData, string("V"), V);
    if (retVal < 0) return -1;
    retVal = getDoubleVec(DoubleFeatureData, StringData, string("T"), t);
    if (retVal < 0) return -1;
    retVal = getIntVec(IntFeatureData, StringData, string("min_AHP_indices"),
                       minAHPIndex);
    if (retVal < 0) return -1;
    retVal = getIntVec(IntFeatureData, StringData, string("peak_indices"),
                       PeakIndex);
    if (retVal < 0) return -1;
    if (PeakIndex.size() <= 1) {
      GErrorStr = GErrorStr +
                  "\nError: More than one spike is needed for spikewidth "
                  "calculation.\n";
      return -1;
    }
    // Take derivative of voltage from 1st AHPmin to the peak of the spike
    // Using Central difference derivative vec1[i] = ((vec[i+1]+vec[i-1])/2)/dx
    retVal = __spike_width2(t, V, PeakIndex, minAHPIndex, spike_width2);
    if (retVal >= 0) {
      setDoubleVec(DoubleFeatureData, StringData, "spike_width2", spike_width2);
    }
    return retVal;
  }
  return -1;
}

// spike half width
// for spike amplitude = v_peak - v_AHP
int LibV1::__spike_width1(const vector<double>& t, const vector<double>& v,
                          const vector<int>& peak_indices,
                          const vector<int>& min_ahp_indices, double stim_start,
                          vector<double>& spike_width1) {
  int start_index = distance(
      t.begin(), find_if(t.begin(), t.end(),
                         bind2nd(greater_equal<double>(), stim_start)));
  vector<int> min_ahp_indices_plus(min_ahp_indices.size() + 1, start_index);
  copy(min_ahp_indices.begin(), min_ahp_indices.end(),
       min_ahp_indices_plus.begin() + 1);
  for (unsigned i = 1; i < min_ahp_indices_plus.size(); i++) {
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

int LibV1::spike_width1(mapStr2intVec& IntFeatureData,
                        mapStr2doubleVec& DoubleFeatureData,
                        mapStr2Str& StringData) {
  int retVal, nSize;
  retVal = CheckInDoublemap(DoubleFeatureData, StringData,
                            string("spike_half_width"), nSize);
  if (retVal)
    return nSize;
  else {
    vector<int> PeakIndex, minAHPIndex;
    vector<double> V, t, dv1, dv2, spike_width1;
    vector<double> stim_start;
    retVal = getDoubleVec(DoubleFeatureData, StringData, string("V"), V);
    if (retVal < 0) return -1;
    retVal = getDoubleVec(DoubleFeatureData, StringData, string("T"), t);
    if (retVal < 0) return -1;
    retVal =
        getDoubleVec(DoubleFeatureData, StringData, "stim_start", stim_start);
    if (retVal < 0) return -1;
    retVal = getIntVec(IntFeatureData, StringData, string("min_AHP_indices"),
                       minAHPIndex);
    if (retVal < 0) return -1;
    retVal = getIntVec(IntFeatureData, StringData, string("peak_indices"),
                       PeakIndex);
    if (retVal < 0) return -1;
    if (PeakIndex.size() <= 1) {
      GErrorStr = GErrorStr +
                  "\nError: More than one spike is needed for spikewidth "
                  "calculation.\n";
      return -1;
    }
    // Take derivative of voltage from 1st AHPmin to the peak of the spike
    // Using Central difference derivative vec1[i] = ((vec[i+1]+vec[i-1])/2)/dx
    retVal = __spike_width1(t, V, PeakIndex, minAHPIndex, stim_start[0],
                            spike_width1);
    if (retVal >= 0) {
      setDoubleVec(DoubleFeatureData, StringData, "spike_half_width",
                   spike_width1);
    }
    return retVal;
  }
  return -1;
}

// passive properties implementation
//
// *** timeconstant ***
// requires hyperpolarizing stimulus
//
// the exponential fit works iteratively
//
int LibV1::__time_constant(const vector<double>& v, const vector<double>& t,
                           double stimStart, double stimEnd,
                           vector<double>& tc) {
  // value of the derivative near the minimum
  double min_derivative = 5e-3;
  // minimal required length of the decay (indices)
  int min_length = 10;
  // minimal required time length in ms
  int t_length = 70;
  int stimstartindex;
  for (stimstartindex = 0; t[stimstartindex] < stimStart; stimstartindex++)
    ;
  stimstartindex += 10;
  // int stimendindex;
  // for(stimendindex = 0; t[stimendindex] < stimEnd; stimendindex++) ;
  // int stimmiddleindex = (stimstartindex + stimendindex) / 2;
  int stimmiddleindex = distance(
      t.begin(),
      find_if(t.begin() + stimstartindex, t.end(),
              bind2nd(greater_equal<double>(), (stimStart + stimEnd) / 2.)));
  if (stimstartindex < 0 || stimstartindex >= v.size() || stimmiddleindex < 0 ||
      stimmiddleindex >= v.size()) {
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
  transform(dv.begin(), dv.end(), dt.begin(), dvdt.begin(), divides<double>());
  // find start of the decay
  int i_start = 0;
  while (find_if(dvdt.begin() + i_start, dvdt.begin() + i_start + 5,
                 bind2nd(greater<double>(), -min_derivative)) !=
         dvdt.begin() + i_start + 5) {
    if (dvdt.begin() + i_start + 5 == dvdt.end()) {
      GErrorStr += "Could not find the decay.\n";
      return -1;
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
    GErrorStr += "\nTrace fall time too short.\n";
    return -1;
  }

  // fit to exponential
  //
  vector<double> log_v(dvdt_decay.size(), 0.);

  // stores slope and relative deviation
  vector<double> slope;

  // golden section search algorithm
  const double PHI = 1.618033988;
  vector<double> x(3, .0);
  // time_constant is searched in between 0 and 200 ms
  x[2] = min_derivative * 200.;
  x[1] = (x[0] * PHI + x[2]) / (1. + PHI);
  // calculate residuals at x[1]
  for (unsigned i = 0; i < log_v.size(); i++) {
    log_v[i] = log(v_decay[i] - v_decay.back() + x[1]);
  }
  slope_straight_line_fit(t_decay, log_v, slope);
  double residuum = slope[1];
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
    for (unsigned i = 0; i < log_v.size(); i++) {
      log_v[i] = log(v_decay[i] - v_decay.back() + newx);
    }
    slope_straight_line_fit(t_decay, log_v, slope);

    if (slope[1] < residuum) {
      if (right) {
        x[0] = x[1];
        x[1] = newx;
      } else {
        x[2] = x[1];
        x[1] = newx;
      }
      residuum = slope[1];
    } else {
      if (right) {
        x[2] = newx;
      } else {
        x[0] = newx;
      }
      right = !right;
    }
  }
  tc.push_back(-1. / slope[0]);
  return 1;
}
int LibV1::time_constant(mapStr2intVec& IntFeatureData,
                         mapStr2doubleVec& DoubleFeatureData,
                         mapStr2Str& StringData) {
  int retVal;
  int nSize;
  retVal = CheckInDoublemap(DoubleFeatureData, StringData,
                            string("time_constant"), nSize);
  if (retVal) {
    return nSize;
  } else {
    vector<double> v;
    vector<double> t;
    vector<double> stimStart;
    vector<double> stimEnd;
    retVal = getDoubleVec(DoubleFeatureData, StringData, string("V"), v);
    if (retVal < 0) return -1;
    retVal = getDoubleVec(DoubleFeatureData, StringData, string("T"), t);
    if (retVal < 0) return -1;
    retVal =
        getDoubleVec(DoubleFeatureData, StringData, "stim_start", stimStart);
    if (retVal < 0) return -1;
    retVal = getDoubleVec(DoubleFeatureData, StringData, "stim_end", stimEnd);
    if (retVal < 0) return -1;
    vector<double> tc;
    retVal = __time_constant(v, t, stimStart[0], stimEnd[0], tc);
    if (retVal >= 0) {
      setDoubleVec(DoubleFeatureData, StringData, "time_constant", tc);
    }
    return retVal;
  }
}

// *** voltage deflection ***
//
int LibV1::voltage_deflection(mapStr2intVec& IntFeatureData,
                              mapStr2doubleVec& DoubleFeatureData,
                              mapStr2Str& StringData) {
  int retVal;
  int nSize;
  retVal = CheckInDoublemap(DoubleFeatureData, StringData,
                            string("voltage_deflection"), nSize);
  if (retVal) {
    return nSize;
  } else {
    vector<double> v;
    vector<double> t;
    vector<double> stimStart;
    vector<double> stimEnd;
    retVal = getDoubleVec(DoubleFeatureData, StringData, string("V"), v);
    if (retVal < 0) return -1;
    retVal = getDoubleVec(DoubleFeatureData, StringData, string("T"), t);
    if (retVal < 0) return -1;
    retVal =
        getDoubleVec(DoubleFeatureData, StringData, "stim_start", stimStart);
    if (retVal < 0) return -1;
    retVal = getDoubleVec(DoubleFeatureData, StringData, "stim_end", stimEnd);
    if (retVal < 0) return -1;
    vector<double> vd;
    retVal = __voltage_deflection(v, t, stimStart[0], stimEnd[0], vd);
    if (retVal >= 0) {
      setDoubleVec(DoubleFeatureData, StringData, "voltage_deflection", vd);
    }
    return retVal;
  }
}
int LibV1::__voltage_deflection(const vector<double>& v,
                                const vector<double>& t, double stimStart,
                                double stimEnd, vector<double>& vd) {
  int wind = 5;
  int stimendindex = 0;
  double base = 0.;
  int base_size = 0;
  for (unsigned i = 0; i < t.size(); i++) {
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
  if (stimendindex - 2 * wind >= 0 && v.size() > 0 &&
      stimendindex - wind < v.size()) {
    for (int i = stimendindex - 2 * wind; i < stimendindex - wind; i++) {
      wind_mean += v[i];
    }
  } else {
    return -1;
  }
  wind_mean /= wind;
  vd.push_back(wind_mean - base);
  return 1;
}

// *** ohmic input resistance ***
//
int LibV1::ohmic_input_resistance(mapStr2intVec& IntFeatureData,
                                  mapStr2doubleVec& DoubleFeatureData,
                                  mapStr2Str& StringData) {
  int retVal;
  int nSize;
  retVal = CheckInDoublemap(DoubleFeatureData, StringData,
                            string("ohmic_input_resistance"), nSize);
  if (retVal) {
    return nSize;
  } else {
    vector<double> voltage_deflection;
    retVal = getDoubleVec(DoubleFeatureData, StringData,
                          string("voltage_deflection"), voltage_deflection);
    if (retVal < 0) {
      return -1;
    }
    vector<double> stimulus_current;
    retVal = getDoubleVec(DoubleFeatureData, StringData, "stimulus_current",
                          stimulus_current);
    if (retVal < 0) {
      return -1;
    }
    vector<double> oir;
    retVal = __ohmic_input_resistance(voltage_deflection[0],
                                      stimulus_current[0], oir);
    if (retVal >= 0) {
      setDoubleVec(DoubleFeatureData, StringData, "ohmic_input_resistance",
                   oir);
    }
    return retVal;
  }
}
int LibV1::__ohmic_input_resistance(double voltage_deflection,
                                    double stimulus_current,
                                    vector<double>& oir) {
  oir.push_back(voltage_deflection / stimulus_current);
  return 1;
}

// *** maximum voltage ***
//
int LibV1::maximum_voltage(mapStr2intVec& IntFeatureData,
                           mapStr2doubleVec& DoubleFeatureData,
                           mapStr2Str& StringData) {
  int retVal;
  int nSize;
  retVal = CheckInDoublemap(DoubleFeatureData, StringData,
                            string("maximum_voltage"), nSize);
  if (retVal) {
    return nSize;
  } else {
    vector<double> v;
    vector<double> t;
    vector<double> stimStart;
    vector<double> stimEnd;
    retVal = getDoubleVec(DoubleFeatureData, StringData, string("V"), v);
    if (retVal < 0) return -1;
    retVal = getDoubleVec(DoubleFeatureData, StringData, string("T"), t);
    if (retVal < 0) return -1;
    retVal =
        getDoubleVec(DoubleFeatureData, StringData, "stim_start", stimStart);
    if (retVal < 0) return -1;
    retVal = getDoubleVec(DoubleFeatureData, StringData, "stim_end", stimEnd);
    if (retVal < 0) return -1;
    vector<double> maxV, minV;
    retVal = __maxmin_voltage(v, t, stimStart[0], stimEnd[0], maxV, minV);
    if (retVal >= 0) {
      setDoubleVec(DoubleFeatureData, StringData, "maximum_voltage", maxV);
    }
    return retVal;
  }
}

// *** maximum voltage ***
//
int LibV1::minimum_voltage(mapStr2intVec& IntFeatureData,
                           mapStr2doubleVec& DoubleFeatureData,
                           mapStr2Str& StringData) {
  int retVal;
  int nSize;
  retVal = CheckInDoublemap(DoubleFeatureData, StringData,
                            string("minimum_voltage"), nSize);
  if (retVal) {
    return nSize;
  } else {
    vector<double> v;
    vector<double> t;
    vector<double> stimStart;
    vector<double> stimEnd;
    retVal = getDoubleVec(DoubleFeatureData, StringData, string("V"), v);
    if (retVal < 0) return -1;
    retVal = getDoubleVec(DoubleFeatureData, StringData, string("T"), t);
    if (retVal < 0) return -1;
    retVal =
        getDoubleVec(DoubleFeatureData, StringData, "stim_start", stimStart);
    if (retVal < 0) return -1;
    retVal = getDoubleVec(DoubleFeatureData, StringData, "stim_end", stimEnd);
    if (retVal < 0) return -1;
    vector<double> maxV, minV;
    retVal = __maxmin_voltage(v, t, stimStart[0], stimEnd[0], maxV, minV);
    if (retVal >= 0) {
      setDoubleVec(DoubleFeatureData, StringData, "minimum_voltage", minV);
    }
    return retVal;
  }
}

int LibV1::__maxmin_voltage(const vector<double>& v, const vector<double>& t,
                            double stimStart, double stimEnd,
                            vector<double>& maxV, vector<double>& minV) {
  int stimstartindex;
  if (stimStart >= t[t.size() - 1]) {
    GErrorStr += "\nStimulus start larger than max time in trace\n";
    return -1;
  }
  if (stimEnd >= t[t.size() - 1]) {
    GErrorStr += "\nStimulus end larger than max time in trace\n";
    return -1;
  }

  for (stimstartindex = 0; t[stimstartindex] < stimStart; stimstartindex++) {
  };
  int stimendindex;
  for (stimendindex = 0; t[stimendindex] < stimEnd; stimendindex++) {
  };
  maxV.push_back(*max_element(&v[stimstartindex], &v[stimendindex]));
  minV.push_back(*min_element(&v[stimstartindex], &v[stimendindex]));
  return 1;
}

// *** steady state voltage ***
//
int LibV1::steady_state_voltage(mapStr2intVec& IntFeatureData,
                                mapStr2doubleVec& DoubleFeatureData,
                                mapStr2Str& StringData) {
  int retVal;
  int nSize;
  retVal = CheckInDoublemap(DoubleFeatureData, StringData,
                            string("steady_state_voltage"), nSize);
  if (retVal) {
    return nSize;
  } else {
    vector<double> v;
    retVal = getDoubleVec(DoubleFeatureData, StringData, string("V"), v);
    if (retVal < 0) return -1;
    vector<double> t;
    retVal = getDoubleVec(DoubleFeatureData, StringData, string("T"), t);
    if (retVal < 0) return -1;
    vector<double> stimEnd;
    retVal = getDoubleVec(DoubleFeatureData, StringData, "stim_end", stimEnd);
    if (retVal < 0) return -1;
    vector<double> ssv;
    retVal = __steady_state_voltage(v, t, stimEnd[0], ssv);
    if (retVal >= 0) {
      setDoubleVec(DoubleFeatureData, StringData, "steady_state_voltage", ssv);
    }
    return retVal;
  }
}
int LibV1::__steady_state_voltage(const vector<double>& v,
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

// *** single_burst_ratio ***
// according to Shaul: measures the length of the first isi over the median of
// the rest of the isis
int LibV1::__single_burst_ratio(const vector<double>& isivalues,
                                vector<double>& singleburstratio) {
  if (isivalues.size() < 2) {
    return 0;
  }
  // calculate the average instead of the median
  double average = 0.;
  for (unsigned i = 1; i < isivalues.size(); i++) {
    average += isivalues[i];
  }
  average /= isivalues.size() - 1;
  singleburstratio.push_back(isivalues[0] / average);
  return singleburstratio.size();
}
int LibV1::single_burst_ratio(mapStr2intVec& IntFeatureData,
                              mapStr2doubleVec& DoubleFeatureData,
                              mapStr2Str& StringData) {
  int retval;
  int nsize;
  retval = CheckInDoublemap(DoubleFeatureData, StringData,
                            string("single_burst_ratio"), nsize);
  if (retval) {
    return nsize;
  }
  vector<double> isivalues;
  retval = getDoubleVec(DoubleFeatureData, StringData, string("ISI_values"),
                        isivalues);
  if (retval < 0) return -1;
  vector<double> singleburstratio;
  retval = __single_burst_ratio(isivalues, singleburstratio);
  if (retval >= 0) {
    setDoubleVec(DoubleFeatureData, StringData, "single_burst_ratio",
                 singleburstratio);
  }
  return retval;
}
// end of single_burst_ratio

// *** threshold_current ***
// just return the value for threshold_current which has been inserted
int LibV1::threshold_current(mapStr2intVec& IntFeatureData,
                             mapStr2doubleVec& DoubleFeatureData,
                             mapStr2Str& StringData) {
  int retval;
  int nsize;
  retval = CheckInDoublemap(DoubleFeatureData, StringData,
                            string("threshold_current"), nsize);
  if (retval) {
    return nsize;
  }
  return retval;
}
// end of threshold_current

int LibV1::printVectorI(char* strName, vector<int> vec) {
  int nSize = 0;
  vector<int>::iterator pos1, pos2;
  nSize = vec.size();
  printf("\nName = [%s] size = [%d] values = [", strName, nSize);
  if (nSize > 0) {
    if (nSize < 100.0) {
      for (unsigned i = 0; i < nSize; i++) {
        printf("%d  ", vec[i]);
      }
    }
    pos1 = max_element(vec.begin(), vec.end());
    pos2 = min_element(vec.begin(), vec.end());
    cout << "max :" << *pos1 << " min :" << *pos2;
  }
  printf("]\n");
  return 0;
}
int LibV1::printVectorD(char* strName, vector<double> vec) {
  int nSize = 0;
  vector<double>::iterator pos1, pos2;
  nSize = vec.size();
  printf("\nName = [%s] size = [%d] values = [", strName, nSize);
  if (nSize > 0) {
    if (nSize < 100.0) {
      for (unsigned i = 0; i < nSize; i++) {
        printf("%f  ", vec[i]);
      }
    }
    pos1 = max_element(vec.begin(), vec.end());
    pos2 = min_element(vec.begin(), vec.end());
    cout << "max :" << *pos1 << " min :" << *pos2;
  }
  printf("]\n");
  return 0;
}

// *** AP_width ***
//
// spike width calculation according to threshold value for the first spike
// unfortunately spike width means the width of the spike on onset not at half
// maximum
int LibV1::__AP_width(const vector<double>& t, const vector<double>& v,
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
  for (unsigned i = 0; i < indices.size() - 1; i++) {
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

int LibV1::AP_width(mapStr2intVec& IntFeatureData,
                    mapStr2doubleVec& DoubleFeatureData,
                    mapStr2Str& StringData) {
  int retval;
  int nsize;
  retval = CheckInDoublemap(DoubleFeatureData, StringData, string("AP_width"),
                            nsize);
  if (retval) {
    return nsize;
  }
  vector<double> t;
  retval = getDoubleVec(DoubleFeatureData, StringData, "T", t);
  if (retval < 0) return -1;
  vector<double> v;
  retval = getDoubleVec(DoubleFeatureData, StringData, "V", v);
  if (retval < 0) return -1;
  vector<double> threshold;
  retval = getDoubleParam(DoubleFeatureData, "Threshold", threshold);
  if (retval < 0) return -1;
  vector<double> stimstart;
  retval = getDoubleVec(DoubleFeatureData, StringData, "stim_start", stimstart);
  if (retval < 0) return -1;
  vector<int> peakindices;
  retval = getIntVec(IntFeatureData, StringData, "peak_indices", peakindices);
  if (retval <= 0) {
    GErrorStr += "\nNo spike in trace.\n";
    return -1;
  }
  vector<int> minahpindices;
  retval =
      getIntVec(IntFeatureData, StringData, "min_AHP_indices", minahpindices);
  if (retval < 0) return -1;
  vector<double> apwidth;
  retval = __AP_width(t, v, stimstart[0], threshold[0], peakindices,
                      minahpindices, apwidth);
  if (retval >= 0) {
    setDoubleVec(DoubleFeatureData, StringData, "AP_width", apwidth);
  }
  return retval;
}
// end of AP_width

// *** doublet_ISI ***
// value of the first ISI
int LibV1::doublet_ISI(mapStr2intVec& IntFeatureData,
                       mapStr2doubleVec& DoubleFeatureData,
                       mapStr2Str& StringData) {
  int retval;
  int nsize;
  retval = CheckInDoublemap(DoubleFeatureData, StringData,
                            string("doublet_ISI"), nsize);
  if (retval) {
    return nsize;
  }
  vector<double> pvt;
  retval =
      getDoubleVec(DoubleFeatureData, StringData, string("peak_time"), pvt);
  if (retval < 2) {
    GErrorStr += "\nNeed at least two spikes for doublet_ISI.\n";
    return -1;
  }
  vector<double> doubletisi(1, pvt[1] - pvt[0]);
  setDoubleVec(DoubleFeatureData, StringData, "doublet_ISI", doubletisi);
  return retval;
}
// end of doublet_ISI

// *** AHP_depth ***
int LibV1::__AHP_depth(const vector<double>& voltagebase,
                       const vector<double>& minahpvalues,
                       vector<double>& ahpdepth) {
  for (unsigned i = 0; i < minahpvalues.size(); i++) {
    ahpdepth.push_back(minahpvalues[i] - voltagebase[0]);
  }
  return ahpdepth.size();
}
int LibV1::AHP_depth(mapStr2intVec& IntFeatureData,
                     mapStr2doubleVec& DoubleFeatureData,
                     mapStr2Str& StringData) {
  int retval;
  int nsize;
  retval = CheckInDoublemap(DoubleFeatureData, StringData, string("AHP_depth"),
                            nsize);
  if (retval) {
    return nsize;
  }
  vector<double> voltagebase;
  retval = getDoubleVec(DoubleFeatureData, StringData, string("voltage_base"),
                        voltagebase);
  if (retval < 0) return -1;
  vector<double> minahpvalues;
  retval = getDoubleVec(DoubleFeatureData, StringData, string("min_AHP_values"),
                        minahpvalues);
  if (retval < 0) return -1;
  vector<double> ahpdepth;
  retval = __AHP_depth(voltagebase, minahpvalues, ahpdepth);
  if (retval >= 0) {
    setDoubleVec(DoubleFeatureData, StringData, "AHP_depth", ahpdepth);
  }
  return retval;
}
// end of AHP_depth

// *** AP_amplitude_diff based on AP_amplitude_change but not normalized  ***
int LibV1::__AP_amplitude_diff(const vector<double>& apamplitude,
                               vector<double>& apamplitudediff) {
  if (apamplitude.size() <= 1) return -1;
  apamplitudediff.resize(apamplitude.size() - 1);
  for (unsigned i = 0; i < apamplitudediff.size(); i++) {
    apamplitudediff[i] = (apamplitude[i + 1] - apamplitude[i]);
  }
  return apamplitudediff.size();
}
int LibV1::AP_amplitude_diff(mapStr2intVec& IntFeatureData,
                             mapStr2doubleVec& DoubleFeatureData,
                             mapStr2Str& StringData) {
  int retval;
  int nsize;
  retval = CheckInDoublemap(DoubleFeatureData, StringData,
                            string("AP_amplitude_diff"), nsize);
  if (retval) {
    return nsize;
  }
  vector<double> apamplitude;
  retval = getDoubleVec(DoubleFeatureData, StringData, string("AP_amplitude"),
                        apamplitude);
  if (retval < 0) return -1;
  vector<double> apamplitudediff;
  retval = __AP_amplitude_diff(apamplitude, apamplitudediff);
  if (retval >= 0) {
    setDoubleVec(DoubleFeatureData, StringData, "AP_amplitude_diff",
                 apamplitudediff);
  }
  return retval;
}
// end of AP_amplitude_diff

// *** AHP_depth_diff, returns AHP_depth[i+1] - AHP_depth[i]  ***
int LibV1::__AHP_depth_diff(const vector<double>& ahpdepth,
                            vector<double>& ahpdepthdiff) {
  if (ahpdepth.size() <= 1) return -1;
  ahpdepthdiff.resize(ahpdepth.size() - 1);
  for (unsigned i = 0; i < ahpdepthdiff.size(); i++) {
    ahpdepthdiff[i] = (ahpdepth[i + 1] - ahpdepth[i]);
  }
  return ahpdepthdiff.size();
}
int LibV1::AHP_depth_diff(mapStr2intVec& IntFeatureData,
                          mapStr2doubleVec& DoubleFeatureData,
                          mapStr2Str& StringData) {
  int retval;
  int nsize;
  retval = CheckInDoublemap(DoubleFeatureData, StringData,
                            string("AHP_depth_diff"), nsize);
  if (retval) {
    return nsize;
  }
  vector<double> ahpdepth;
  retval = getDoubleVec(DoubleFeatureData, StringData, string("AHP_depth"),
                        ahpdepth);
  if (retval < 0) return -1;
  vector<double> ahpdepthdiff;
  retval = __AHP_depth_diff(ahpdepth, ahpdepthdiff);
  if (retval >= 0) {
    setDoubleVec(DoubleFeatureData, StringData, "AHP_depth_diff", ahpdepthdiff);
  }
  return retval;
}
// end of AHP_depth_diff

// end of feature definition
