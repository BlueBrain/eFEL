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
