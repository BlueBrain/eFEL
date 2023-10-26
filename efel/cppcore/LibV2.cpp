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

#include "LibV2.h"

#include <algorithm>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <math.h>

using std::find_if;
using std::min_element;
using std::max_element;
using std::transform;

// AP parameters
//
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
int LibV2::AP_rise_indices(mapStr2intVec& IntFeatureData,
                           mapStr2doubleVec& DoubleFeatureData,
                           mapStr2Str& StringData) {
  int retVal;
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
  apfi.resize(std::min(apbi.size(), pi.size()));
  for (size_t i = 0; i < apfi.size(); i++) {
    double halfheight = (v[pi[i]] + v[apbi[i]]) / 2.;
    vector<double> vpeak(&v[pi[i]], &v[apei[i]]);
    transform(vpeak.begin(), vpeak.end(), vpeak.begin(),
              [halfheight](double val) { return fabs(val - halfheight); });
    apfi[i] = distance(vpeak.begin(), min_element(vpeak.begin(), vpeak.end())) +
              pi[i];
  }
  return apfi.size();
}
int LibV2::AP_fall_indices(mapStr2intVec& IntFeatureData,
                           mapStr2doubleVec& DoubleFeatureData,
                           mapStr2Str& StringData) {
  int retVal;
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
  apduration.resize(std::min(apbeginindices.size(), endindices.size()));
  for (size_t i = 0; i < apduration.size(); i++) {
    // printf("%d, %d, %d\n", t.size(), apbeginindices.size(),
    // endindices.size())
    apduration[i] = t[endindices[i]] - t[apbeginindices[i]];
  }
  return apduration.size();
}
int LibV2::AP_duration(mapStr2intVec& IntFeatureData,
                       mapStr2doubleVec& DoubleFeatureData,
                       mapStr2Str& StringData) {
  int retval;
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
int LibV2::AP_duration_half_width(mapStr2intVec& IntFeatureData,
                                  mapStr2doubleVec& DoubleFeatureData,
                                  mapStr2Str& StringData) {
  int retval;
  vector<double> t;
  retval = getVec(DoubleFeatureData, StringData, "T", t);
  if (retval < 0) return -1;
  vector<int> apriseindices;
  retval = getVec(IntFeatureData, StringData, "AP_rise_indices",
                     apriseindices);
  if (retval < 0) return -1;
  vector<int> apfallindices;
  retval = getVec(IntFeatureData, StringData, "AP_fall_indices",
                     apfallindices);
  if (retval < 0) return -1;
  vector<double> apdurationhalfwidth;
  retval = __AP_duration_half_width(t, apriseindices, apfallindices,
                                    apdurationhalfwidth);
  if (retval >= 0) {
    setVec(DoubleFeatureData, StringData, "AP_duration_half_width",
                 apdurationhalfwidth);
  }
  return retval;
}
// end of AP_duration_half_width

// *** AP_rise_time according to E9 and E17 ***
static int __AP_rise_time(const vector<double>& t, 
                          const vector<double>& v,
                          const vector<int>& apbeginindices,
                          const vector<int>& peakindices,
                          const vector<double>& apamplitude,
                          double beginperc,
                          double endperc,
                          vector<double>& aprisetime) {
  aprisetime.resize(std::min(apbeginindices.size(), peakindices.size()));
  double begin_v;
  double end_v;
  double begin_indice;
  double end_indice;
  for (size_t i = 0; i < aprisetime.size(); i++) {
    begin_v = v[apbeginindices[i]] + beginperc * apamplitude[i];
    end_v = v[apbeginindices[i]] + endperc * apamplitude[i];

    // Get begin indice
    size_t j=apbeginindices[i];
    // change slightly begin_v for almost equal case
    // truncature error can change begin_v even when beginperc == 0.0
    while (j<peakindices[i] && v[j] < begin_v - 0.0000000000001){
      j++;
    }
    begin_indice = j;

    // Get end indice
    j=peakindices[i];
    // change slightly end_v for almost equal case
    // truncature error can change end_v even when beginperc == 0.0
    while (j>apbeginindices[i] && v[j] > end_v + 0.0000000000001){
      j--;
    }
    end_indice = j;
    
    aprisetime[i] = t[end_indice] - t[begin_indice];
  }
  return aprisetime.size();
}
int LibV2::AP_rise_time(mapStr2intVec& IntFeatureData,
                        mapStr2doubleVec& DoubleFeatureData,
                        mapStr2Str& StringData) {
  int retval;
  vector<double> t;
  retval = getVec(DoubleFeatureData, StringData, "T", t);
  if (retval < 0) return -1;
  vector<int> apbeginindices;
  retval = getVec(IntFeatureData, StringData, "AP_begin_indices",
                     apbeginindices);
  if (retval < 0) return -1;
  vector<int> peakindices;
  retval = getVec(IntFeatureData, StringData, "peak_indices",
                     peakindices);
  if (retval < 0) return -1;
  vector<double> v;
  retval = getVec(DoubleFeatureData, StringData, "V", v);
  if (retval < 0) return -1;
  vector<double> AP_amplitude;
  retval =
      getVec(DoubleFeatureData, StringData, "AP_amplitude", AP_amplitude);
  if (retval < 0) {
    GErrorStr += "Error calculating AP_amplitude for mean_AP_amplitude";
    return -1;
  } else if (retval == 0) {
    GErrorStr += "No spikes found when calculating mean_AP_amplitude";
    return -1;
  } else if (AP_amplitude.size() == 0) {
    GErrorStr += "No spikes found when calculating mean_AP_amplitude";
    return -1;
  }
  // Get rise begin percentage
  vector<double> risebeginperc;
  retval = getVec(DoubleFeatureData, StringData, "rise_start_perc", risebeginperc);
  if (retval <= 0) {
    risebeginperc.push_back(0.0);
  }
  // Get rise end percentage
  vector<double> riseendperc;
  retval = getVec(DoubleFeatureData, StringData, "rise_end_perc", riseendperc);
  if (retval <= 0) {
    riseendperc.push_back(1.0);
  }
  vector<double> aprisetime;
  retval = __AP_rise_time(t, v, apbeginindices, peakindices, AP_amplitude, risebeginperc[0], riseendperc[0], aprisetime);
  if (retval >= 0) {
    setVec(DoubleFeatureData, StringData, "AP_rise_time", aprisetime);
  }
  return retval;
}
// end of AP_rise_time

// *** AP_fall_time according to E10 and E18 ***
static int __AP_fall_time(const vector<double>& t,
                          const vector<int>& peakindices,
                          const vector<int>& apendindices,
                          vector<double>& apfalltime) {
  apfalltime.resize(std::min(peakindices.size(), apendindices.size()));
  for (size_t i = 0; i < apfalltime.size(); i++) {
    apfalltime[i] = t[apendindices[i]] - t[peakindices[i]];
  }
  return apfalltime.size();
}
int LibV2::AP_fall_time(mapStr2intVec& IntFeatureData,
                        mapStr2doubleVec& DoubleFeatureData,
                        mapStr2Str& StringData) {
  int retval;
  vector<double> t;
  retval = getVec(DoubleFeatureData, StringData, "T", t);
  if (retval < 0) return -1;
  vector<int> peakindices;
  retval = getVec(IntFeatureData, StringData, "peak_indices",
                     peakindices);
  if (retval < 0) return -1;
  vector<int> apendindices;
  retval = getVec(IntFeatureData, StringData, "AP_end_indices",
                     apendindices);
  if (retval < 0) return -1;
  vector<double> apfalltime;
  retval = __AP_fall_time(t, peakindices, apendindices, apfalltime);
  if (retval >= 0) {
    setVec(DoubleFeatureData, StringData, "AP_fall_time", apfalltime);
  }
  return retval;
}
// end of AP_fall_time

// *** AP_rise_rate according to E11 and E19 ***
static int __AP_rise_rate(const vector<double>& t, const vector<double>& v,
                          const vector<int>& apbeginindices,
                          const vector<int>& peakindices,
                          vector<double>& apriserate) {
  apriserate.resize(std::min(peakindices.size(), apbeginindices.size()));
  for (size_t i = 0; i < apriserate.size(); i++) {
    apriserate[i] = (v[peakindices[i]] - v[apbeginindices[i]]) /
                    (t[peakindices[i]] - t[apbeginindices[i]]);
  }
  return apriserate.size();
}
int LibV2::AP_rise_rate(mapStr2intVec& IntFeatureData,
                        mapStr2doubleVec& DoubleFeatureData,
                        mapStr2Str& StringData) {
  int retval;
  vector<double> t;
  retval = getVec(DoubleFeatureData, StringData, "T", t);
  if (retval < 0) return -1;
  vector<double> v;
  retval = getVec(DoubleFeatureData, StringData, "V", v);
  if (retval < 0) return -1;
  vector<int> apbeginindices;
  retval = getVec(IntFeatureData, StringData, "AP_begin_indices",
                     apbeginindices);
  if (retval < 0) return -1;
  vector<int> peakindices;
  retval = getVec(IntFeatureData, StringData, "peak_indices",
                     peakindices);
  if (retval < 0) return -1;
  vector<double> apriserate;
  retval = __AP_rise_rate(t, v, apbeginindices, peakindices, apriserate);
  if (retval >= 0) {
    setVec(DoubleFeatureData, StringData, "AP_rise_rate", apriserate);
  }
  return retval;
}
// end of AP_rise_rate

// *** AP_fall_rate according to E12 and E20 ***
static int __AP_fall_rate(const vector<double>& t, const vector<double>& v,
                          const vector<int>& peakindices,
                          const vector<int>& apendindices,
                          vector<double>& apfallrate) {
  apfallrate.resize(std::min(apendindices.size(), peakindices.size()));
  for (size_t i = 0; i < apfallrate.size(); i++) {
    apfallrate[i] = (v[apendindices[i]] - v[peakindices[i]]) /
                    (t[apendindices[i]] - t[peakindices[i]]);
  }
  return apfallrate.size();
}
int LibV2::AP_fall_rate(mapStr2intVec& IntFeatureData,
                        mapStr2doubleVec& DoubleFeatureData,
                        mapStr2Str& StringData) {
  int retval;
  vector<double> t;
  retval = getVec(DoubleFeatureData, StringData, "T", t);
  if (retval < 0) return -1;
  vector<double> v;
  retval = getVec(DoubleFeatureData, StringData, "V", v);
  if (retval < 0) return -1;
  vector<int> peakindices;
  retval = getVec(IntFeatureData, StringData, "peak_indices", peakindices);
  if (retval < 0) return -1;
  vector<int> apendindices;
  retval = getVec(IntFeatureData, StringData, "AP_end_indices",
                     apendindices);
  if (retval < 0) return -1;
  vector<double> apfallrate;
  retval = __AP_fall_rate(t, v, peakindices, apendindices, apfallrate);
  if (retval >= 0) {
    setVec(DoubleFeatureData, StringData, "AP_fall_rate", apfallrate);
  }
  return retval;
}
// end of AP_fall_rate

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
int LibV2::fast_AHP(mapStr2intVec& IntFeatureData,
                    mapStr2doubleVec& DoubleFeatureData,
                    mapStr2Str& StringData) {
  int retval;
  vector<double> v;
  retval = getVec(DoubleFeatureData, StringData, "V", v);
  if (retval < 0) return -1;
  vector<int> apbeginindices;
  retval = getVec(IntFeatureData, StringData, "AP_begin_indices",
                     apbeginindices);
  if (retval < 0) return -1;
  vector<int> minahpindices;
  retval = getVec(IntFeatureData, StringData, "min_AHP_indices",
                     minahpindices);
  if (retval < 0) return -1;
  vector<double> fastahp;
  retval = __fast_AHP(v, apbeginindices, minahpindices, fastahp);
  if (retval >= 0) {
    setVec(DoubleFeatureData, StringData, "fast_AHP", fastahp);
  }
  return retval;
}
// end of fast_AHP

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
int LibV2::AP_amplitude_change(mapStr2intVec& IntFeatureData,
                               mapStr2doubleVec& DoubleFeatureData,
                               mapStr2Str& StringData) {
  int retval;
  vector<double> apamplitude;
  retval = getVec(DoubleFeatureData, StringData, "AP_amplitude",
                        apamplitude);
  if (retval < 0) return -1;
  vector<double> apamplitudechange;
  retval = __AP_amplitude_change(apamplitude, apamplitudechange);
  if (retval >= 0) {
    setVec(DoubleFeatureData, StringData, "AP_amplitude_change",
                 apamplitudechange);
  }
  return retval;
}
// end of AP_amplitude_change

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
int LibV2::AP_duration_change(mapStr2intVec& IntFeatureData,
                              mapStr2doubleVec& DoubleFeatureData,
                              mapStr2Str& StringData) {
  int retval;
  vector<double> apduration;
  retval = getVec(DoubleFeatureData, StringData, "AP_duration",
                        apduration);
  if (retval < 0) return -1;
  vector<double> apdurationchange;
  retval = __AP_duration_change(apduration, apdurationchange);
  if (retval >= 0) {
    setVec(DoubleFeatureData, StringData, "AP_duration_change",
                 apdurationchange);
  }
  return retval;
}
// end of AP_duration_change

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
int LibV2::AP_duration_half_width_change(mapStr2intVec& IntFeatureData,
                                         mapStr2doubleVec& DoubleFeatureData,
                                         mapStr2Str& StringData) {
  int retval;
  vector<double> apdurationhalfwidth;
  retval = getVec(DoubleFeatureData, StringData,
                        "AP_duration_half_width", apdurationhalfwidth);
  if (retval < 0) return -1;
  vector<double> apdurationhalfwidthchange;
  retval = __AP_duration_half_width_change(apdurationhalfwidth,
                                           apdurationhalfwidthchange);
  if (retval >= 0) {
    setVec(DoubleFeatureData, StringData, "AP_duration_half_width_change",
                 apdurationhalfwidthchange);
  }
  return retval;
}
// end of AP_duration_half_width_change

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
int LibV2::AP_rise_rate_change(mapStr2intVec& IntFeatureData,
                               mapStr2doubleVec& DoubleFeatureData,
                               mapStr2Str& StringData) {
  int retval;
  vector<double> apriserate;
  retval = getVec(DoubleFeatureData, StringData, "AP_rise_rate",
                        apriserate);
  if (retval < 0) return -1;
  vector<double> apriseratechange;
  retval = __AP_rise_rate_change(apriserate, apriseratechange);
  if (retval >= 0) {
    setVec(DoubleFeatureData, StringData, "AP_rise_rate_change",
                 apriseratechange);
  }
  return retval;
}
// end of AP_rise_rate_change

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
int LibV2::AP_fall_rate_change(mapStr2intVec& IntFeatureData,
                               mapStr2doubleVec& DoubleFeatureData,
                               mapStr2Str& StringData) {
  int retval;
  vector<double> apfallrate;
  retval = getVec(DoubleFeatureData, StringData, "AP_fall_rate",
                        apfallrate);
  if (retval < 0) return -1;
  vector<double> apfallratechange;
  retval = __AP_fall_rate_change(apfallrate, apfallratechange);
  if (retval >= 0) {
    setVec(DoubleFeatureData, StringData, "AP_fall_rate_change",
                 apfallratechange);
  }
  return retval;
}
// end of AP_fall_rate_change

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
int LibV2::fast_AHP_change(mapStr2intVec& IntFeatureData,
                           mapStr2doubleVec& DoubleFeatureData,
                           mapStr2Str& StringData) {
  int retval;
  vector<double> fastahp;
  retval = getVec(DoubleFeatureData, StringData, "fast_AHP", fastahp);
  if (retval < 0) return -1;
  vector<double> fastahpchange;
  retval = __fast_AHP_change(fastahp, fastahpchange);
  if (retval >= 0) {
    setVec(DoubleFeatureData, StringData, "fast_AHP_change",
                 fastahpchange);
  }
  return retval;
}
// end of fast_AHP_change

// *** E6 ***
int LibV2::E6(mapStr2intVec& IntFeatureData,
              mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData) {
  int retval;
  vector<double> e6;
  retval = mean_traces_double(DoubleFeatureData, "AP_amplitude", "APWaveForm",
                              0, e6);
  if (retval >= 0) {
    e6.resize(1);
    setVec(DoubleFeatureData, StringData, "E6", e6);
  }
  return retval;
}
// end of E6

// *** E7 ***
int LibV2::E7(mapStr2intVec& IntFeatureData,
              mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData) {
  int retval;
  vector<double> e7;
  retval =
      mean_traces_double(DoubleFeatureData, "AP_duration", "APWaveForm", 0, e7);
  if (retval >= 0) {
    e7.resize(1);
    setVec(DoubleFeatureData, StringData, "E7", e7);
  }
  return retval;
}
// end of E7

// *** BPAPatt2 ***
int LibV2::BPAPatt2(mapStr2intVec& IntFeatureData,
                    mapStr2doubleVec& DoubleFeatureData,
                    mapStr2Str& StringData) {
  int retval;
  vector<double> peakvoltage;
  retval = getDoubleParam(DoubleFeatureData, "peak_voltage;location_soma",
                          peakvoltage);
  // one spike required
  if (retval <= 0) return -1;
  // voltage base
  vector<double> vb_dend;
  retval = getDoubleParam(DoubleFeatureData, "voltage_base;location_dend620",
                          vb_dend);
  if (retval <= 0) return -1;
  vector<double> v_dend;
  retval = getDoubleParam(DoubleFeatureData, "V;location_dend620", v_dend);
  if (retval <= 0) return -1;
  vector<double> vb_soma;
  retval =
      getDoubleParam(DoubleFeatureData, "voltage_base;location_soma", vb_soma);
  if (retval <= 0) return -1;
  vector<double> bpapatt;
  // this is according to Etay's definition of the backpropagating action
  // potential:
  // the ratio of the height of somatic and dendritic spike
  // bpapatt.push_back((peakvoltage[0] - vb_soma[0]) /
  // (*max_element(v_dend.begin(), v_dend.end()) - vb_dend[0]));
  //
  // this is according to the definition of the error of the backpropagating
  // action potential:
  // the height of the dendritic spike
  bpapatt.push_back(*max_element(v_dend.begin(), v_dend.end()) - vb_dend[0]);
  setVec(DoubleFeatureData, StringData, "BPAPatt2", bpapatt);
  return retval;
}
// end of BPAPatt2

// *** BPAPatt3 ***
int LibV2::BPAPatt3(mapStr2intVec& IntFeatureData,
                    mapStr2doubleVec& DoubleFeatureData,
                    mapStr2Str& StringData) {
  int retval;
  vector<double> peakvoltage;
  retval = getDoubleParam(DoubleFeatureData, "peak_voltage;location_soma",
                          peakvoltage);
  // one spike required
  if (retval <= 0) return -1;
  // voltage base
  vector<double> vb_dend;
  retval = getDoubleParam(DoubleFeatureData, "voltage_base;location_dend800",
                          vb_dend);
  if (retval <= 0) return -1;
  vector<double> v_dend;
  retval = getDoubleParam(DoubleFeatureData, "V;location_dend800", v_dend);
  if (retval <= 0) return -1;
  vector<double> vb_soma;
  retval =
      getDoubleParam(DoubleFeatureData, "voltage_base;location_soma", vb_soma);
  if (retval <= 0) return -1;
  vector<double> bpapatt;
  // this is according to Etay's definition of the backpropagating action
  // potential:
  // the ratio of the height of somatic and dendritic spike
  // bpapatt.push_back((peakvoltage[0] - vb_soma[0]) /
  // (*max_element(v_dend.begin(), v_dend.end()) - vb_dend[0]));
  //
  // this is according to the definition of the error of the backpropagating
  // action potential:
  // the height of the dendritic spike
  bpapatt.push_back(*max_element(v_dend.begin(), v_dend.end()) - vb_dend[0]);
  setVec(DoubleFeatureData, StringData, "BPAPatt3", bpapatt);
  return retval;
}
// end of BPAPatt3

// *** E39 ***
// unit: Hz / nA
int LibV2::E39(mapStr2intVec& IntFeatureData,
               mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData) {
  vector<string> stim_params;
  // retrieve the complete suffixes of all traces where the suffix matches
  // "IDthreshold":
  getTraces(DoubleFeatureData, "IDthreshold", stim_params);
  if (stim_params.size() > 1) {
    vector<double> current(stim_params.size());
    vector<double> frequency(stim_params.size());
    // iterate over these traces:
    for (size_t i = 0; i < stim_params.size(); i++) {
      vector<double> stimulus_current;
      // retrieve the trace data
      // note that we call getDoubleParam with suffix appended,
      // this is the direct access to the global map
      getDoubleParam(DoubleFeatureData, "stimulus_current" + stim_params[i],
                     stimulus_current);
      current[i] = stimulus_current[0];
      vector<double> freq;
      getDoubleParam(DoubleFeatureData, "mean_frequency" + stim_params[i],
                     freq);
      frequency[i] = freq[0];
    }
    linear_fit_result fit;
    fit = slope_straight_line_fit(current, frequency);
    vector<double> e39(1, fit.slope);
    vector<double> e39_cod(1, fit.r_square);
    setVec(DoubleFeatureData, StringData, "E39", e39);
    setVec(DoubleFeatureData, StringData, "E39_cod", e39_cod);
    return e39.size();
  }
  GErrorStr += "\nMore than 1 trace required for calculation of E39";
  return -1;
}
// end of E39

// *** E39_cod ***
// coefficient of determination for the slope of the IF curve
int LibV2::E39_cod(mapStr2intVec& IntFeatureData,
                   mapStr2doubleVec& DoubleFeatureData,
                   mapStr2Str& StringData) {
  return 1;
}
// end of E39_cod

// *** amp_drop_first_second ***
static int __amp_drop_first_second(const vector<double>& peakvoltage,
                                   vector<double>& ampdropfirstsecond) {
  ampdropfirstsecond.push_back(peakvoltage[0] - peakvoltage[1]);
  return ampdropfirstsecond.size();
}
int LibV2::amp_drop_first_second(mapStr2intVec& IntFeatureData,
                                 mapStr2doubleVec& DoubleFeatureData,
                                 mapStr2Str& StringData) {
  int retval;
  vector<double> peakvoltage;
  retval = getVec(DoubleFeatureData, StringData, "peak_voltage",
                        peakvoltage);
  if (retval < 2) {
    GErrorStr +=
        "At least 2 spikes needed for calculation of amp_drop_first_second.\n";
    return -1;
  }
  vector<double> ampdropfirstsecond;
  retval = __amp_drop_first_second(peakvoltage, ampdropfirstsecond);
  if (retval >= 0) {
    setVec(DoubleFeatureData, StringData, "amp_drop_first_second",
                 ampdropfirstsecond);
  }
  return retval;
}
// end of amp_drop_first_second

// *** E2 ***
int LibV2::E2(mapStr2intVec& IntFeatureData,
              mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData) {
  int retval;
  vector<double> e2;
  retval = mean_traces_double(DoubleFeatureData, "amp_drop_first_second",
                              "APDrop", 0, e2);
  if (retval > 0) {
    setVec(DoubleFeatureData, StringData, "E2", e2);
    return 1;
  }
  return retval;
}
// end of E2

// *** amp_drop_first_last ***
static int __amp_drop_first_last(const vector<double>& peakvoltage,
                                 vector<double>& ampdropfirstlast) {
  ampdropfirstlast.push_back(peakvoltage[0] - peakvoltage.back());
  return ampdropfirstlast.size();
}
int LibV2::amp_drop_first_last(mapStr2intVec& IntFeatureData,
                               mapStr2doubleVec& DoubleFeatureData,
                               mapStr2Str& StringData) {
  int retval;
  vector<double> peakvoltage;
  retval = getVec(DoubleFeatureData, StringData, "peak_voltage",
                        peakvoltage);
  if (retval < 2) {
    GErrorStr +=
        "At least 2 spikes needed for calculation of amp_drop_first_last.\n";
    return -1;
  }
  vector<double> ampdropfirstlast;
  retval = __amp_drop_first_last(peakvoltage, ampdropfirstlast);
  if (retval >= 0) {
    setVec(DoubleFeatureData, StringData, "amp_drop_first_last",
                 ampdropfirstlast);
  }
  return retval;
}
// end of amp_drop_first_last

// *** E3 ***
int LibV2::E3(mapStr2intVec& IntFeatureData,
              mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData) {
  int retval;
  vector<double> e3;
  retval = mean_traces_double(DoubleFeatureData, "amp_drop_first_last",
                              "APDrop", 0, e3);
  if (retval > 0) {
    setVec(DoubleFeatureData, StringData, "E3", e3);
    return 1;
  }
  return retval;
}
// end of E3

// *** amp_drop_second_last ***
static int __amp_drop_second_last(const vector<double>& peakvoltage,
                                  vector<double>& ampdropsecondlast) {
  ampdropsecondlast.push_back(peakvoltage[1] - peakvoltage.back());
  return ampdropsecondlast.size();
}
int LibV2::amp_drop_second_last(mapStr2intVec& IntFeatureData,
                                mapStr2doubleVec& DoubleFeatureData,
                                mapStr2Str& StringData) {
  int retval;
  vector<double> peakvoltage;
  retval = getVec(DoubleFeatureData, StringData, "peak_voltage",
                        peakvoltage);
  if (retval < 3) {
    GErrorStr +=
        "At least 3 spikes needed for calculation of amp_drop_second_last.\n";
    return -1;
  }
  vector<double> ampdropsecondlast;
  retval = __amp_drop_second_last(peakvoltage, ampdropsecondlast);
  if (retval >= 0) {
    setVec(DoubleFeatureData, StringData, "amp_drop_second_last",
                 ampdropsecondlast);
  }
  return retval;
}
// end of amp_drop_second_last

// *** E4 ***
int LibV2::E4(mapStr2intVec& IntFeatureData,
              mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData) {
  int retval;
  vector<double> e4;
  retval = mean_traces_double(DoubleFeatureData, "amp_drop_second_last",
                              "APDrop", 0, e4);
  if (retval > 0) {
    setVec(DoubleFeatureData, StringData, "E4", e4);
    return 1;
  }
  return retval;
}
// end of E4

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

int LibV2::max_amp_difference(mapStr2intVec& IntFeatureData,
                              mapStr2doubleVec& DoubleFeatureData,
                              mapStr2Str& StringData) {
  int retval;
  vector<double> peakvoltage;
  retval = getVec(DoubleFeatureData, StringData, "peak_voltage",
                        peakvoltage);
  if (retval < 2) {
    GErrorStr +=
        "At least 2 spikes needed for calculation of max_amp_difference.\n";
    return -1;
  }
  vector<double> maxampdifference;
  retval = __max_amp_difference(peakvoltage, maxampdifference);
  if (retval >= 0) {
    setVec(DoubleFeatureData, StringData, "max_amp_difference",
                 maxampdifference);
  }
  return retval;
}
// end of max_amp_difference

// *** E5 ***
int LibV2::E5(mapStr2intVec& IntFeatureData,
              mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData) {
  int retval;
  vector<double> e5;
  retval = mean_traces_double(DoubleFeatureData, "max_amp_difference", "APDrop",
                              0, e5);
  if (retval > 0) {
    setVec(DoubleFeatureData, StringData, "E5", e5);
    return 1;
  }
  return retval;
}
// end of E5

// *** E8 ***
int LibV2::E8(mapStr2intVec& IntFeatureData,
              mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData) {
  int retval;
  vector<double> e8;
  retval = mean_traces_double(DoubleFeatureData, "AP_duration_half_width",
                              "APWaveForm", 0, e8);
  if (retval >= 0) {
    e8.resize(1);
    setVec(DoubleFeatureData, StringData, "E8", e8);
  }
  return retval;
}
// end of E8

// *** E9 ***
int LibV2::E9(mapStr2intVec& IntFeatureData,
              mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData) {
  int retval;
  vector<double> e9;
  retval = mean_traces_double(DoubleFeatureData, "AP_rise_time", "APWaveForm",
                              0, e9);
  if (retval >= 0) {
    e9.resize(1);
    setVec(DoubleFeatureData, StringData, "E9", e9);
  }
  return retval;
}
// end of E9

// *** E10 ***
int LibV2::E10(mapStr2intVec& IntFeatureData,
               mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData) {
  int retval;
  vector<double> e10;
  retval = mean_traces_double(DoubleFeatureData, "AP_fall_time", "APWaveForm",
                              0, e10);
  if (retval >= 0) {
    e10.resize(1);
    setVec(DoubleFeatureData, StringData, "E10", e10);
  }
  return retval;
}
// end of E10

// *** E11 ***
int LibV2::E11(mapStr2intVec& IntFeatureData,
               mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData) {
  int retval;
  vector<double> e11;
  retval = mean_traces_double(DoubleFeatureData, "AP_rise_rate", "APWaveForm",
                              0, e11);
  if (retval >= 0) {
    e11.resize(1);
    setVec(DoubleFeatureData, StringData, "E11", e11);
  }
  return retval;
}
// end of E11

// *** E12 ***
int LibV2::E12(mapStr2intVec& IntFeatureData,
               mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData) {
  int retval;
  vector<double> e12;
  retval = mean_traces_double(DoubleFeatureData, "AP_fall_rate", "APWaveForm",
                              0, e12);
  if (retval >= 0) {
    e12.resize(1);
    setVec(DoubleFeatureData, StringData, "E12", e12);
  }
  return retval;
}
// end of E12

// *** E13 ***
int LibV2::E13(mapStr2intVec& IntFeatureData,
               mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData) {
  int retval;
  vector<double> e13;
  retval =
      mean_traces_double(DoubleFeatureData, "fast_AHP", "APWaveForm", 0, e13);
  if (retval >= 0) {
    e13.resize(1);
    setVec(DoubleFeatureData, StringData, "E13", e13);
  }
  return retval;
}
// end of E13

// *** E14 ***
int LibV2::E14(mapStr2intVec& IntFeatureData,
               mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData) {
  int retval;
  vector<double> e14;
  retval = mean_traces_double(DoubleFeatureData, "peak_voltage", "APWaveForm",
                              0, e14);
  if (retval >= 0) {
    e14[0] = e14[1];
    e14.resize(1);
    setVec(DoubleFeatureData, StringData, "E14", e14);
  }
  return retval;
}
// end of E14

// *** E15 ***
int LibV2::E15(mapStr2intVec& IntFeatureData,
               mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData) {
  int retval;
  vector<double> e15;
  retval = mean_traces_double(DoubleFeatureData, "AP_duration", "APWaveForm", 0,
                              e15);
  if (retval >= 0) {
    e15[0] = e15[1];
    e15.resize(1);
    setVec(DoubleFeatureData, StringData, "E15", e15);
  }
  return retval;
}
// end of E15

// *** E16 ***
int LibV2::E16(mapStr2intVec& IntFeatureData,
               mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData) {
  int retval;
  vector<double> e16;
  retval = mean_traces_double(DoubleFeatureData, "AP_duration_half_width",
                              "APWaveForm", 0, e16);
  if (retval >= 0) {
    e16[0] = e16[1];
    e16.resize(1);
    setVec(DoubleFeatureData, StringData, "E16", e16);
  }
  return retval;
}
// end of E16

// *** E17 ***
int LibV2::E17(mapStr2intVec& IntFeatureData,
               mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData) {
  int retval;
  vector<double> e17;
  retval = mean_traces_double(DoubleFeatureData, "AP_rise_time", "APWaveForm",
                              0, e17);
  if (retval >= 0) {
    e17[0] = e17[1];
    e17.resize(1);
    setVec(DoubleFeatureData, StringData, "E17", e17);
  }
  return retval;
}
// end of E17

// *** E18 ***
int LibV2::E18(mapStr2intVec& IntFeatureData,
               mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData) {
  int retval;
  vector<double> e18;
  retval = mean_traces_double(DoubleFeatureData, "AP_fall_time", "APWaveForm",
                              0, e18);
  if (retval >= 0) {
    e18[0] = e18[1];
    e18.resize(1);
    setVec(DoubleFeatureData, StringData, "E18", e18);
  }
  return retval;
}
// end of E18

// *** E19 ***
int LibV2::E19(mapStr2intVec& IntFeatureData,
               mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData) {
  int retval;
  vector<double> e19;
  retval = mean_traces_double(DoubleFeatureData, "AP_rise_rate", "APWaveForm",
                              0, e19);
  if (retval >= 0) {
    e19[0] = e19[1];
    e19.resize(1);
    setVec(DoubleFeatureData, StringData, "E19", e19);
  }
  return retval;
}
// end of E19

// *** E20 ***
int LibV2::E20(mapStr2intVec& IntFeatureData,
               mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData) {
  int retval;
  vector<double> e20;
  retval = mean_traces_double(DoubleFeatureData, "AP_fall_rate", "APWaveForm",
                              0, e20);
  if (retval >= 0) {
    e20[0] = e20[1];
    e20.resize(1);
    setVec(DoubleFeatureData, StringData, "E20", e20);
  }
  return retval;
}
// end of E20

// *** E21 ***
int LibV2::E21(mapStr2intVec& IntFeatureData,
               mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData) {
  int retval;
  vector<double> e21;
  retval =
      mean_traces_double(DoubleFeatureData, "fast_AHP", "APWaveForm", 0, e21);
  if (retval >= 0) {
    e21[0] = e21[1];
    e21.resize(1);
    setVec(DoubleFeatureData, StringData, "E21", e21);
  }
  return retval;
}
// end of E21

// *** E22 ***
int LibV2::E22(mapStr2intVec& IntFeatureData,
               mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData) {
  int retval;
  vector<double> e22;
  retval = mean_traces_double(DoubleFeatureData, "AP_amplitude_change",
                              "APWaveForm", 0, e22);
  if (retval >= 0) {
    setVec(DoubleFeatureData, StringData, "E22", e22);
  }
  return retval;
}
// end of E22

// *** E23 ***
int LibV2::E23(mapStr2intVec& IntFeatureData,
               mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData) {
  int retval;
  vector<double> e23;
  retval = mean_traces_double(DoubleFeatureData, "AP_duration_change",
                              "APWaveForm", 0, e23);
  if (retval >= 0) {
    setVec(DoubleFeatureData, StringData, "E23", e23);
  }
  return retval;
}
// end of E23

// *** E24 ***
int LibV2::E24(mapStr2intVec& IntFeatureData,
               mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData) {
  int retval;
  vector<double> e24;
  retval = mean_traces_double(
      DoubleFeatureData, "AP_duration_half_width_change", "APWaveForm", 0, e24);
  if (retval >= 0) {
    setVec(DoubleFeatureData, StringData, "E24", e24);
  }
  return retval;
}
// end of E24

// *** E25 ***
int LibV2::E25(mapStr2intVec& IntFeatureData,
               mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData) {
  int retval;
  vector<double> e25;
  retval = mean_traces_double(DoubleFeatureData, "AP_rise_rate_change",
                              "APWaveForm", 0, e25);
  if (retval >= 0) {
    setVec(DoubleFeatureData, StringData, "E25", e25);
  }
  return retval;
}
// end of E25

// *** E26 ***
int LibV2::E26(mapStr2intVec& IntFeatureData,
               mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData) {
  int retval;
  vector<double> e26;
  retval = mean_traces_double(DoubleFeatureData, "AP_fall_rate_change",
                              "APWaveForm", 0, e26);
  if (retval >= 0) {
    setVec(DoubleFeatureData, StringData, "E26", e26);
  }
  return retval;
}
// end of E26

// *** E27 ***
int LibV2::E27(mapStr2intVec& IntFeatureData,
               mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData) {
  int retval;
  vector<double> e27;
  retval = mean_traces_double(DoubleFeatureData, "fast_AHP_change",
                              "APWaveForm", 0, e27);
  if (retval >= 0) {
    setVec(DoubleFeatureData, StringData, "E27", e27);
  }
  return retval;
}
// end of E27

// steady state of the voltage response during hyperpolarizing stimulus,
// elementary feature for E29
// *** steady_state_hyper
static int __steady_state_hyper(const vector<double>& v,
                                const vector<double>& t, double stimend,
                                vector<double>& steady_state_hyper) {
  // Find the iterator pointing to the first time value greater than or equal to stimend
  auto it_stimend = find_if(t.begin(), t.end(), 
                            [stimend](double t_val) { return t_val >= stimend; });

  // Calculate the index, ensuring you account for the offset of -5
  int i_end = distance(t.begin(), it_stimend) - 5;


  const int offset = 30;
  if (i_end < 0 || i_end < offset) {
    return -1;
  }

  size_t i_begin = static_cast<size_t>(i_end - offset);

  double sum = 0.;

  for (size_t i = i_begin; i < static_cast<size_t>(i_end); i++) {
    sum += v[i];
  }

  double mean = sum / (i_end - i_begin);
  steady_state_hyper.push_back(mean);
  return 1;
}

int LibV2::steady_state_hyper(mapStr2intVec& IntFeatureData,
                              mapStr2doubleVec& DoubleFeatureData,
                              mapStr2Str& StringData) {
  int retval;
  vector<double> v;
  retval = getVec(DoubleFeatureData, StringData, "V", v);
  if (retval < 0) return -1;
  vector<double> t;
  retval = getVec(DoubleFeatureData, StringData, "T", t);
  if (retval < 0) return -1;
  vector<double> stimend;
  retval = getVec(DoubleFeatureData, StringData, "stim_end", stimend);
  if (retval < 0) return -1;
  vector<double> steady_state_hyper;
  retval = __steady_state_hyper(v, t, stimend[0], steady_state_hyper);
  if (retval >= 0) {
    setVec(DoubleFeatureData, StringData, "steady_state_hyper",
                 steady_state_hyper);
  }
  return retval;
}

// *** E40 ***
int LibV2::E40(mapStr2intVec& IntFeatureData,
               mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData) {
  int retval;
  vector<double> e40;
  retval = mean_traces_double(DoubleFeatureData, "time_to_first_spike",
                              "IDrest", 0, e40);
  if (retval >= 0) {
    setVec(DoubleFeatureData, StringData, "E40", e40);
  }
  return retval;
}
// end of E40

//
// end of feature definition
