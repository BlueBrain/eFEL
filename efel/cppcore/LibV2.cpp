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
    // Retrieve all required features at once
    const auto& features = getFeatures(DoubleFeatureData, {"V", "T", "stim_end"});

    vector<double> steady_state_hyper;
    int retval = __steady_state_hyper(
        features.at("V"),
        features.at("T"),
        features.at("stim_end").front(),
        steady_state_hyper
    );

    if (retval >= 0) {
        setVec(DoubleFeatureData, StringData, "steady_state_hyper", steady_state_hyper);
    }
    return retval;
}


//
// end of feature definition
