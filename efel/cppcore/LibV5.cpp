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

#include "LibV5.h"

#include <math.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <deque>
#include <functional>
#include <iterator>

// slope of loglog of ISI curve
static int __ISI_log_slope(const vector<double>& isiValues,
                           vector<double>& slope, bool skip, double spikeSkipf,
                           size_t maxnSpike, bool semilog) {
  std::deque<double> skippedISIValues;

  vector<double> log_isivalues;
  vector<double> x;

  for (size_t i = 0; i < isiValues.size(); i++) {
    skippedISIValues.push_back(isiValues[i]);
  }

  if (skip) {
    // Remove n spikes given by spike_skipf or max_spike_skip
    size_t isisToRemove =
        (size_t)((isiValues.size() + 1) * spikeSkipf + .5);

    isisToRemove = std::min(maxnSpike, isisToRemove);

    // Remove spikeToRemove spike from SpikeTime list
    for (size_t i = 0; i < isisToRemove; i++) {
      skippedISIValues.pop_front();
    }
  }

  for (size_t i = 0; i < skippedISIValues.size(); i++) {
    log_isivalues.push_back(log(skippedISIValues[i]));
    if (semilog) {
      x.push_back((double)i + 1);
    } else {
      x.push_back(log((double)i + 1));
    }
  }

  if (x.size() == 0 || log_isivalues.size() == 0) return -1;

  linear_fit_result fit;
  fit = slope_straight_line_fit(x, log_isivalues);

  if (fit.slope == 0. || is_nan(fit.slope)) return -1;

  slope.push_back(fit.slope);

  return slope.size();
}

int LibV5::ISI_log_slope(mapStr2intVec& IntFeatureData,
                         mapStr2doubleVec& DoubleFeatureData,
                         mapStr2Str& StringData) {
  int size;
  if (CheckInMap(DoubleFeatureData, StringData, "ISI_log_slope", size)) {
    return size;
  }
  vector<double> isivalues;
  vector<double> slope;
  if (getVec(DoubleFeatureData, StringData, "ISI_values", isivalues) <=
      0) {
    return -1;
  }
  bool semilog = false;
  int retval = __ISI_log_slope(isivalues, slope, false, 0, 0, semilog);
  if (retval >= 0) {
    setVec(DoubleFeatureData, StringData, "ISI_log_slope", slope);
    return slope.size();
  } else {
    return retval;
  }
}

int LibV5::ISI_semilog_slope(mapStr2intVec& IntFeatureData,
                             mapStr2doubleVec& DoubleFeatureData,
                             mapStr2Str& StringData) {
  int size;
  if (CheckInMap(DoubleFeatureData, StringData, "ISI_semilog_slope",
                       size)) {
    return size;
  }
  vector<double> isivalues;
  vector<double> slope;
  if (getVec(DoubleFeatureData, StringData, "ISI_values", isivalues) <=
      0) {
    return -1;
  }
  bool semilog = true;
  int retval = __ISI_log_slope(isivalues, slope, false, 0, 0, semilog);
  if (retval >= 0) {
    setVec(DoubleFeatureData, StringData, "ISI_semilog_slope", slope);
    return slope.size();
  } else {
    return retval;
  }
}

int LibV5::ISI_log_slope_skip(mapStr2intVec& IntFeatureData,
                              mapStr2doubleVec& DoubleFeatureData,
                              mapStr2Str& StringData) {
  int size, retVal;
  vector<int> maxnSpike;
  vector<double> spikeSkipf;

  if (CheckInMap(DoubleFeatureData, StringData, "ISI_log_slope_skip",
                       size)) {
    return size;
  }
  vector<double> isivalues;
  vector<double> slope;
  if (getVec(DoubleFeatureData, StringData, "ISI_values", isivalues) <=
      0) {
    return -1;
  }
  retVal = getDoubleParam(DoubleFeatureData, "spike_skipf", spikeSkipf);
  {
    if (retVal <= 0) return -1;
  };
  // spikeSkipf is a fraction hence value should lie between >=0 and <1. [0 1)
  if ((spikeSkipf[0] < 0) || (spikeSkipf[0] >= 1)) {
    GErrorStr += "\nspike_skipf should lie between [0 1).\n";
    return -1;
  }
  retVal = getIntParam(IntFeatureData, "max_spike_skip", maxnSpike);
  {
    if (retVal <= 0) return -1;
  };

  bool semilog = false;
  retVal = __ISI_log_slope(isivalues, slope, true, spikeSkipf[0], maxnSpike[0],
                           semilog);
  if (retVal >= 0) {
    setVec(DoubleFeatureData, StringData, "ISI_log_slope_skip", slope);
    return slope.size();
  } else {
    return retVal;
  }
}

// time from stimulus start to second threshold crossing
int LibV5::time_to_second_spike(mapStr2intVec& IntFeatureData,
                                mapStr2doubleVec& DoubleFeatureData,
                                mapStr2Str& StringData) {
  int retVal, nSize;
  retVal = CheckInMap(DoubleFeatureData, StringData,
                            "time_to_second_spike", nSize);
  if (retVal) return nSize;

  vector<double> second_spike;
  vector<double> peaktime;
  vector<double> stimstart;
  retVal = getVec(DoubleFeatureData, StringData, "peak_time", peaktime);
  if (retVal < 2) {
    GErrorStr += "\n Two spikes required for time_to_second_spike.\n";
    return -1;
  }
  retVal = getVec(DoubleFeatureData, StringData, "stim_start", stimstart);
  if (retVal <= 0) return -1;

  second_spike.push_back(peaktime[1] - stimstart[0]);
  setVec(DoubleFeatureData, StringData, "time_to_second_spike",
               second_spike);
  return second_spike.size();
}

// time from stimulus start to last spike
int LibV5::time_to_last_spike(mapStr2intVec& IntFeatureData,
                              mapStr2doubleVec& DoubleFeatureData,
                              mapStr2Str& StringData) {
  int retVal, nSize;

  retVal = CheckInMap(DoubleFeatureData, StringData, "time_to_last_spike",
                            nSize);
  if (retVal) return nSize;

  vector<double> last_spike;
  vector<double> peaktime;
  vector<double> stimstart;
  retVal = getVec(DoubleFeatureData, StringData, "peak_time", peaktime);
  if (retVal < 0) {
    GErrorStr += "\n Error in peak_time calculation in time_to_last_spike.\n";
    return -1;
  } else if (retVal == 0) {
    last_spike.push_back(0.0);
    setVec(DoubleFeatureData, StringData, "time_to_last_spike",
                 last_spike);
  } else {
    retVal =
        getVec(DoubleFeatureData, StringData, "stim_start", stimstart);
    if (retVal <= 0) return -1;
    last_spike.push_back(peaktime[peaktime.size() - 1] - stimstart[0]);
    setVec(DoubleFeatureData, StringData, "time_to_last_spike",
                 last_spike);
  }
  return 1;
}

// 1.0 over first ISI (in Hz); returns 0 when no ISI
int LibV5::inv_first_ISI(mapStr2intVec& IntFeatureData,
                         mapStr2doubleVec& DoubleFeatureData,
                         mapStr2Str& StringData) {
  int retVal, nSize;
  retVal =
      CheckInMap(DoubleFeatureData, StringData, "inv_first_ISI", nSize);
  if (retVal) return nSize;

  vector<double> all_isi_values_vec;
  double inv_first_ISI;
  vector<double> inv_first_ISI_vec;
  retVal = getVec(DoubleFeatureData, StringData, "all_ISI_values",
                        all_isi_values_vec);
  if (retVal < 1) {
    inv_first_ISI = 0.0;
  } else {
    inv_first_ISI = 1000.0 / all_isi_values_vec[0];
  }
  inv_first_ISI_vec.push_back(inv_first_ISI);
  setVec(DoubleFeatureData, StringData, "inv_first_ISI",
               inv_first_ISI_vec);
  return 1;
}

// 1.0 over second ISI (in Hz); returns 0 when no ISI
int LibV5::inv_second_ISI(mapStr2intVec& IntFeatureData,
                          mapStr2doubleVec& DoubleFeatureData,
                          mapStr2Str& StringData) {
  int retVal, nSize;
  retVal =
      CheckInMap(DoubleFeatureData, StringData, "inv_second_ISI", nSize);
  if (retVal) return nSize;

  vector<double> all_isi_values_vec;
  double inv_second_ISI;
  vector<double> inv_second_ISI_vec;
  retVal = getVec(DoubleFeatureData, StringData, "all_ISI_values",
                        all_isi_values_vec);
  if (retVal < 2) {
    inv_second_ISI = 0.0;
  } else {
    inv_second_ISI = 1000.0 / all_isi_values_vec[1];
  }
  inv_second_ISI_vec.push_back(inv_second_ISI);
  setVec(DoubleFeatureData, StringData, "inv_second_ISI",
               inv_second_ISI_vec);
  return 1;
}

// 1.0 over third ISI (in Hz); returns 0 when no ISI
int LibV5::inv_third_ISI(mapStr2intVec& IntFeatureData,
                         mapStr2doubleVec& DoubleFeatureData,
                         mapStr2Str& StringData) {
  int retVal, nSize;
  retVal =
      CheckInMap(DoubleFeatureData, StringData, "inv_third_ISI", nSize);
  if (retVal) return nSize;

  vector<double> all_isi_values_vec;
  double inv_third_ISI;
  vector<double> inv_third_ISI_vec;
  retVal = getVec(DoubleFeatureData, StringData, "all_ISI_values",
                        all_isi_values_vec);
  if (retVal < 3) {
    inv_third_ISI = 0.0;
  } else {
    inv_third_ISI = 1000.0 / all_isi_values_vec[2];
  }
  inv_third_ISI_vec.push_back(inv_third_ISI);
  setVec(DoubleFeatureData, StringData, "inv_third_ISI",
               inv_third_ISI_vec);
  return 1;
}

// 1.0 over fourth ISI (in Hz); returns 0 when no ISI
int LibV5::inv_fourth_ISI(mapStr2intVec& IntFeatureData,
                          mapStr2doubleVec& DoubleFeatureData,
                          mapStr2Str& StringData) {
  int retVal, nSize;
  retVal =
      CheckInMap(DoubleFeatureData, StringData, "inv_fourth_ISI", nSize);
  if (retVal) return nSize;

  vector<double> all_isi_values_vec;
  double inv_fourth_ISI;
  vector<double> inv_fourth_ISI_vec;
  retVal = getVec(DoubleFeatureData, StringData, "all_ISI_values",
                        all_isi_values_vec);
  if (retVal < 4) {
    inv_fourth_ISI = 0.0;
  } else {
    inv_fourth_ISI = 1000.0 / all_isi_values_vec[3];
  }
  inv_fourth_ISI_vec.push_back(inv_fourth_ISI);
  setVec(DoubleFeatureData, StringData, "inv_fourth_ISI",
               inv_fourth_ISI_vec);
  return 1;
}

// 1.0 over fifth ISI (in Hz); returns 0 when no ISI
int LibV5::inv_fifth_ISI(mapStr2intVec& IntFeatureData,
                         mapStr2doubleVec& DoubleFeatureData,
                         mapStr2Str& StringData) {
  int retVal, nSize;
  retVal =
      CheckInMap(DoubleFeatureData, StringData, "inv_fifth_ISI", nSize);
  if (retVal) return nSize;

  vector<double> all_isi_values_vec;
  double inv_fifth_ISI;
  vector<double> inv_fifth_ISI_vec;
  retVal = getVec(DoubleFeatureData, StringData, "all_ISI_values",
                        all_isi_values_vec);
  if (retVal < 5) {
    inv_fifth_ISI = 0.0;
  } else {
    inv_fifth_ISI = 1000.0 / all_isi_values_vec[4];
  }
  inv_fifth_ISI_vec.push_back(inv_fifth_ISI);
  setVec(DoubleFeatureData, StringData, "inv_fifth_ISI",
               inv_fifth_ISI_vec);
  return 1;
}

// 1.0 over last ISI (in Hz); returns 0 when no ISI
int LibV5::inv_last_ISI(mapStr2intVec& IntFeatureData,
                        mapStr2doubleVec& DoubleFeatureData,
                        mapStr2Str& StringData) {
  int retVal, nSize;
  retVal =
      CheckInMap(DoubleFeatureData, StringData, "inv_last_ISI", nSize);
  if (retVal) return nSize;
  vector<double> all_isi_values_vec;
  double inv_last_ISI;
  vector<double> inv_last_ISI_vec;
  retVal = getVec(DoubleFeatureData, StringData, "all_ISI_values",
                        all_isi_values_vec);
  if (retVal < 1) {
    inv_last_ISI = 0.0;
  } else {
    inv_last_ISI = 1000.0 / all_isi_values_vec[all_isi_values_vec.size() - 1];
  }
  inv_last_ISI_vec.push_back(inv_last_ISI);
  setVec(DoubleFeatureData, StringData, "inv_last_ISI", inv_last_ISI_vec);
  return 1;
}

// 1.0 over time to first spike (in Hz); returns 0 when no spike
int LibV5::inv_time_to_first_spike(mapStr2intVec& IntFeatureData,
                                   mapStr2doubleVec& DoubleFeatureData,
                                   mapStr2Str& StringData) {
  int retVal, nSize;
  retVal = CheckInMap(DoubleFeatureData, StringData,
                            "inv_time_to_first_spike", nSize);
  if (retVal) return nSize;

  vector<double> time_to_first_spike_vec;
  double inv_time_to_first_spike;
  vector<double> inv_time_to_first_spike_vec;
  retVal = getVec(DoubleFeatureData, StringData, "time_to_first_spike",
                        time_to_first_spike_vec);
  if (retVal < 1) {
    inv_time_to_first_spike = 0.0;
  } else {
    inv_time_to_first_spike = 1000.0 / time_to_first_spike_vec[0];
  }
  inv_time_to_first_spike_vec.push_back(inv_time_to_first_spike);
  setVec(DoubleFeatureData, StringData, "inv_time_to_first_spike",
               inv_time_to_first_spike_vec);
  return 1;
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
    end_index =
        distance(t.begin(),
                 find_if(t.begin(), t.end(),
                         std::bind2nd(std::greater_equal<double>(), stim_end)));
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
int LibV5::min_AHP_indices(mapStr2intVec& IntFeatureData,
                           mapStr2doubleVec& DoubleFeatureData,
                           mapStr2Str& StringData) {
  int retVal, nSize;

  retVal = CheckInMap(IntFeatureData, StringData, "min_AHP_indices", nSize);
  if (retVal) return nSize;

  double stim_start, stim_end;
  vector<int> min_ahp_indices, strict_stiminterval_vec, peak_indices;
  vector<double> v, t, stim_start_vec, stim_end_vec, min_ahp_values;
  bool strict_stiminterval;

  // Get voltage
  retVal = getVec(DoubleFeatureData, StringData, "V", v);
  if (retVal <= 0) return -1;

  // Get time
  retVal = getVec(DoubleFeatureData, StringData, "T", t);
  if (retVal <= 0) return -1;

  // Get peak_indices
  retVal = getVec(IntFeatureData, StringData, "peak_indices", peak_indices);
  if (retVal < 1) {
    GErrorStr +=
        "\n At least one spike required for calculation of "
        "min_AHP_indices.\n";
    return -1;
  }

  // Get strict_stiminterval
  retVal = getIntParam(IntFeatureData, "strict_stiminterval",
                       strict_stiminterval_vec);
  if (retVal <= 0) {
    strict_stiminterval = false;
  } else {
    strict_stiminterval = bool(strict_stiminterval_vec[0]);
  }

  // Get stim_start
  retVal =
      getVec(DoubleFeatureData, StringData, "stim_start", stim_start_vec);
  if (retVal <= 0) {
    return -1;
  } else {
    stim_start = stim_start_vec[0];
  }

  /// Get stim_end
  retVal =
      getVec(DoubleFeatureData, StringData, "stim_end", stim_end_vec);
  if (retVal <= 0) {
    return -1;
  } else {
    stim_end = stim_end_vec[0];
  }

  retVal =
      __min_AHP_indices(t, v, peak_indices, stim_start, stim_end,
                        strict_stiminterval, min_ahp_indices, min_ahp_values);

  if (retVal > 0) {
    setVec(IntFeatureData, StringData, "min_AHP_indices", min_ahp_indices);
    setVec(DoubleFeatureData, StringData, "min_AHP_values",
                 min_ahp_values);
  }
  return retVal;
}

int LibV5::min_AHP_values(mapStr2intVec& IntFeatureData,
                          mapStr2doubleVec& DoubleFeatureData,
                          mapStr2Str& StringData) {
  int retVal, nSize;
  retVal =
      CheckInMap(DoubleFeatureData, StringData, "min_AHP_values", nSize);
  if (retVal) return nSize;
  return -1;
}

// Difference with LibV1 is that this function doesn't return -1 if there are no
// min_AHP_values
int LibV5::AHP_depth_abs(mapStr2intVec& IntFeatureData,
                         mapStr2doubleVec& DoubleFeatureData,
                         mapStr2Str& StringData) {
  int retVal, nSize;
  retVal =
      CheckInMap(DoubleFeatureData, StringData, "AHP_depth_abs", nSize);
  if (retVal) return nSize;

  vector<double> vAHP;
  retVal = getVec(DoubleFeatureData, StringData, "min_AHP_values", vAHP);
  if (retVal < 0) return -1;
  setVec(DoubleFeatureData, StringData, "AHP_depth_abs", vAHP);
  return vAHP.size();
}

// spike half width
// for spike amplitude = v_peak - v_AHP
static int __spike_width1(const vector<double>& t, const vector<double>& v,
                          const vector<int>& peak_indices,
                          const vector<int>& min_ahp_indices, double stim_start,
                          vector<double>& spike_width1) {
  int start_index =
      distance(t.begin(),
               find_if(t.begin(), t.end(),
                       std::bind2nd(std::greater_equal<double>(), stim_start)));
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
                           std::bind2nd(std::greater_equal<double>(), v_half)));
    v_dev = v_half - v[rise_index];
    delta_v = v[rise_index] - v[rise_index - 1];
    delta_t = t[rise_index] - t[rise_index - 1];
    t_dev_rise = delta_t * v_dev / delta_v;
    int fall_index = distance(
        v.begin(), find_if(v.begin() + peak_indices[i - 1],
                           v.begin() + min_ahp_indices_plus[i],
                           std::bind2nd(std::less_equal<double>(), v_half)));
    v_dev = v_half - v[fall_index];
    delta_v = v[fall_index] - v[fall_index - 1];
    delta_t = t[fall_index] - t[fall_index - 1];
    t_dev_fall = delta_t * v_dev / delta_v;
    spike_width1.push_back(t[fall_index] + t_dev_rise - t[rise_index] +
                           t_dev_fall);
  }
  return spike_width1.size();
}

int LibV5::spike_width1(mapStr2intVec& IntFeatureData,
                        mapStr2doubleVec& DoubleFeatureData,
                        mapStr2Str& StringData) {
  int retVal, nSize;
  retVal = CheckInMap(DoubleFeatureData, StringData, "spike_half_width",
                            nSize);
  if (retVal) return nSize;

  vector<int> PeakIndex, minAHPIndex;
  vector<double> V, t, dv1, dv2, spike_width1;
  vector<double> stim_start;
  retVal = getVec(DoubleFeatureData, StringData, "V", V);
  if (retVal < 0) return -1;
  retVal = getVec(DoubleFeatureData, StringData, "T", t);
  if (retVal < 0) return -1;
  retVal =
      getVec(DoubleFeatureData, StringData, "stim_start", stim_start);
  if (retVal < 0) return -1;
  retVal =
      getVec(IntFeatureData, StringData, "min_AHP_indices", minAHPIndex);
  if (retVal < 0) return -1;
  retVal = getVec(IntFeatureData, StringData, "peak_indices", PeakIndex);
  if (retVal < 0) return -1;

  // if(PeakIndex.size()<1) {GErrorStr = GErrorStr + "\nError: One spike is
  // needed for spikewidth calculation.\n"; return -1;}
  if (PeakIndex.size() == 0 || minAHPIndex.size() == 0) {
    setVec(DoubleFeatureData, StringData, "spike_half_width",
                 spike_width1);
    return 0;
  }
  // Take derivative of voltage from 1st AHPmin to the peak of the spike
  // Using Central difference derivative vec1[i] = ((vec[i+1]+vec[i-1])/2)/dx
  retVal =
      __spike_width1(t, V, PeakIndex, minAHPIndex, stim_start[0], spike_width1);
  if (retVal >= 0) {
    setVec(DoubleFeatureData, StringData, "spike_half_width",
                 spike_width1);
  }
  return retVal;
}

//
// *** AP begin indices ***
//
static int __AP_begin_indices(const vector<double>& t, const vector<double>& v,
                              double stimstart, double stimend,
                              const vector<int>& ahpi, vector<int>& apbi,
                              double dTh, int derivative_window) {
  const double derivativethreshold = dTh;
  vector<double> dvdt(v.size());
  vector<double> dv;
  vector<double> dt;
  getCentralDifferenceDerivative(1., v, dv);
  getCentralDifferenceDerivative(1., t, dt);
  transform(dv.begin(), dv.end(), dt.begin(), dvdt.begin(),
            std::divides<double>());

  /*for (size_t i = 0; i < dvdt.size(); i++) {
      printf("%d %f %f\n", i, dvdt[i]);
  }
  printf("\n");*/

  // restrict to time interval where stimulus is applied
  vector<int> minima;
  int stimbeginindex =
      distance(t.begin(),
               find_if(t.begin(), t.end(),
                       std::bind2nd(std::greater_equal<double>(), stimstart)));
  minima.push_back(stimbeginindex);
  for (size_t i = 0; i < ahpi.size(); i++) {
    if (ahpi[i] > stimbeginindex) {
      minima.push_back(ahpi[i]);
    }
    // if(t[ahpi[i]] > stimend) {
    //    break;
    //}
  }
  // if the AHP_indices are already restricted make sure that we do not miss
  // the last spike
  // if(t[minima.back()] < stimend) {
  int endindex = distance(t.begin(), t.end());
  minima.push_back(endindex);
  //}
  // printf("Found %d minima\n", minima.size());
  for (size_t i = 0; i < minima.size() - 1; i++) {
    // assure that the width of the slope is bigger than 4
    int newbegin = minima[i];
    int begin = minima[i];
    int width = derivative_window;
    bool skip = false;

    // Detect where the derivate crosses derivativethreshold, and make sure
    // this happens in a window of 'width' sampling point
    do {
      begin = distance(
          dvdt.begin(),
          find_if(
              dvdt.begin() + newbegin, dvdt.begin() + minima[i + 1],
              std::bind2nd(std::greater_equal<double>(), derivativethreshold)));
      // printf("%d %d\n", newbegin, minima[i+1]);
      if (begin == minima[i + 1]) {
        // printf("Skipping %d %d\n", newbegin, minima[i+1]);
        // could not find a spike in between these minima
        skip = true;
        break;
      }
      newbegin = begin + 1;
    } while (find_if(dvdt.begin() + begin, dvdt.begin() + begin + width,
                     std::bind2nd(std::less<double>(), derivativethreshold)) !=
             dvdt.begin() + begin + width);
    if (skip) {
      continue;
    }
    apbi.push_back(begin);
  }
  return apbi.size();
}

int LibV5::AP_begin_indices(mapStr2intVec& IntFeatureData,
                            mapStr2doubleVec& DoubleFeatureData,
                            mapStr2Str& StringData) {
  int retVal;
  int nSize;

  // Check if calculated already
  retVal = CheckInMap(IntFeatureData, StringData, "AP_begin_indices", nSize);
  if (retVal) {
    return nSize;
  }

  // Get input parameters
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

  // Get DerivativeThreshold
  vector<double> dTh;
  retVal = getDoubleParam(DoubleFeatureData, "DerivativeThreshold", dTh);
  if (retVal <= 0) {
    // derivative at peak start according to eCode specification 10mV/ms
    // according to Shaul 12mV/ms
    dTh.push_back(12.0);
  }

  // Get DerivativeWindow
  vector<int> derivative_window;
  retVal = getIntParam(IntFeatureData, "DerivativeWindow", derivative_window);
  if (retVal <= 0) {
    GErrorStr += "\nDerivativeWindow not set\n";
    return -1;
  }
  
  // Calculate feature
  retVal =
      __AP_begin_indices(t, v, stimstart[0], stimend[0], ahpi, apbi, dTh[0], derivative_window[0]);

  // Save feature value
  if (retVal >= 0) {
    setVec(IntFeatureData, StringData, "AP_begin_indices", apbi);
  }
  return retVal;
}

static int __AP_end_indices(const vector<double>& t, const vector<double>& v,
                            const vector<int>& pi, vector<int>& apei,
                            double derivativethreshold) {

  vector<double> dvdt(v.size());
  vector<double> dv;
  vector<double> dt;
  getCentralDifferenceDerivative(1., v, dv);
  getCentralDifferenceDerivative(1., t, dt);
  transform(dv.begin(), dv.end(), dt.begin(), dvdt.begin(),
            std::divides<double>());

  apei.resize(pi.size());
  vector<int> picopy(pi.begin(), pi.end());
  picopy.push_back(v.size() - 1);

  for (size_t i = 0; i < apei.size(); i++) {
    // assure that the width of the slope is bigger than 4
    apei[i] = std::distance(
        dvdt.begin(),
        std::find_if(dvdt.begin() + picopy[i] + 1, dvdt.begin() + picopy[i + 1],
                std::bind2nd(std::greater_equal<double>(), derivativethreshold)));
  }
  return apei.size();
}


int LibV5::AP_end_indices(mapStr2intVec& IntFeatureData,
                          mapStr2doubleVec& DoubleFeatureData,
                          mapStr2Str& StringData) {

  int retVal;
  int nSize;
  retVal = CheckInMap(IntFeatureData, StringData, "AP_end_indices", nSize);
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

  // Get DerivativeThreshold
  vector<double> dTh;
  retVal = getDoubleParam(DoubleFeatureData, "DownDerivativeThreshold", dTh);
  if (retVal <= 0) {
    // derivative at peak end
    dTh.push_back(-12.0);
  }

  vector<int> apei;
  retVal = __AP_end_indices(t, v, pi, apei, dTh[0]);
  if (retVal >= 0) {
    setVec(IntFeatureData, StringData, "AP_end_indices", apei);
  }
  return retVal;
}


static int __irregularity_index(vector<double>& isiValues,
                                vector<double>& irregularity_index) {
  double ISISub, iRI;
  iRI = ISISub = 0;
  if (isiValues.size() == 0) return -1;

  for (size_t i = 1; i < isiValues.size(); i++) {
    ISISub = std::abs(isiValues[i] - isiValues[i - 1]);
    iRI = iRI + (ISISub);
  }
  iRI = iRI / isiValues.size();
  irregularity_index.clear();
  irregularity_index.push_back(iRI);
  return 1;
}

int LibV5::irregularity_index(mapStr2intVec& IntFeatureData,
                              mapStr2doubleVec& DoubleFeatureData,
                              mapStr2Str& StringData) {
  int retVal, nSize;
  retVal = CheckInMap(DoubleFeatureData, StringData, "irregularity_index",
                            nSize);
  if (retVal) return nSize;

  vector<double> isiValues, irregularity_index;
  retVal = getVec(DoubleFeatureData, StringData, "ISI_values", isiValues);
  if (retVal < 0) return -1;

  retVal = __irregularity_index(isiValues, irregularity_index);
  if (retVal >= 0)
    setVec(DoubleFeatureData, StringData, "irregularity_index",
                 irregularity_index);
  return retVal;
}

static int __number_initial_spikes(vector<double>& peak_times, double stimstart,
                                   double stimend, double initial_perc,
                                   vector<int>& number_initial_spikes) {
  double initialLength = (stimend - stimstart) * initial_perc;

  int startIndex =
      distance(peak_times.begin(),
               find_if(peak_times.begin(), peak_times.end(),
                       std::bind2nd(std::greater_equal<double>(), stimstart)));
  int endIndex = distance(peak_times.begin(),
                          find_if(peak_times.begin(), peak_times.end(),
                                  std::bind2nd(std::greater_equal<double>(),
                                               stimstart + initialLength)));

  number_initial_spikes.push_back(endIndex - startIndex);

  return 1;
}

// Number of spikes in the initial_perc interval
int LibV5::number_initial_spikes(mapStr2intVec& IntFeatureData,
                                 mapStr2doubleVec& DoubleFeatureData,
                                 mapStr2Str& StringData) {
  int retVal, nSize;
  retVal = CheckInMap(DoubleFeatureData, StringData,
                            "number_initial_spikes", nSize);
  if (retVal) return nSize;

  vector<double> peak_times, initial_perc;
  vector<int> number_initial_spikes;
  retVal = getVec(DoubleFeatureData, StringData, "peak_time", peak_times);
  if (retVal < 0) return -1;

  retVal = getDoubleParam(DoubleFeatureData, "initial_perc", initial_perc);
  if (retVal <= 0) return -1;
  if ((initial_perc[0] < 0) || (initial_perc[0] >= 1)) {
    GErrorStr += "\ninitial_perc should lie between [0 1).\n";
    return -1;
  }

  vector<double> stimstart;
  retVal = getVec(DoubleFeatureData, StringData, "stim_start", stimstart);
  if (retVal < 0) return -1;

  vector<double> stimend;
  retVal = getVec(DoubleFeatureData, StringData, "stim_end", stimend);
  if (retVal < 0) return -1;

  retVal = __number_initial_spikes(peak_times, stimstart[0], stimend[0],
                                   initial_perc[0], number_initial_spikes);
  if (retVal >= 0)
    setVec(IntFeatureData, StringData, "number_initial_spikes",
              number_initial_spikes);
  return retVal;
}

// Amplitude of the first spike
int LibV5::AP1_amp(mapStr2intVec& IntFeatureData,
                   mapStr2doubleVec& DoubleFeatureData,
                   mapStr2Str& StringData) {
  int retVal, nSize;
  retVal = CheckInMap(DoubleFeatureData, StringData, "AP1_amp", nSize);
  if (retVal) return nSize;

  vector<double> AP_amplitudes, AP1_amp;
  retVal = getVec(DoubleFeatureData, StringData, "AP_amplitude",
                        AP_amplitudes);
  if (retVal < 1) {
    setVec(DoubleFeatureData, StringData, "AP1_amp", AP1_amp);
    return 0;
  } else {
    AP1_amp.push_back(AP_amplitudes[0]);
  }

  setVec(DoubleFeatureData, StringData, "AP1_amp", AP1_amp);
  return 1;
}

// Amplitude of the first spike
int LibV5::APlast_amp(mapStr2intVec& IntFeatureData,
                      mapStr2doubleVec& DoubleFeatureData,
                      mapStr2Str& StringData) {
  int retVal, nSize;
  retVal = CheckInMap(DoubleFeatureData, StringData, "APlast_amp", nSize);
  if (retVal) return nSize;
  vector<double> AP_amplitudes, APlast_amp;
  int AP_amplitude_size;

  AP_amplitude_size = getVec(DoubleFeatureData, StringData,
                                   "AP_amplitude", AP_amplitudes);
  if (AP_amplitude_size < 1) {
    setVec(DoubleFeatureData, StringData, "APlast_amp", APlast_amp);
    return 0;
  } else {
    APlast_amp.push_back(AP_amplitudes[AP_amplitude_size - 1]);
  }

  setVec(DoubleFeatureData, StringData, "APlast_amp", APlast_amp);
  return 1;
}

// Peak voltage of the first spike
int LibV5::AP1_peak(mapStr2intVec& IntFeatureData,
                    mapStr2doubleVec& DoubleFeatureData,
                    mapStr2Str& StringData) {
  int retVal, nSize;
  retVal = CheckInMap(DoubleFeatureData, StringData, "AP1_peak", nSize);
  if (retVal) return nSize;

  vector<double> peak_voltage, AP1_peak;
  retVal =
      getVec(DoubleFeatureData, StringData, "peak_voltage", peak_voltage);
  if (retVal < 1) {
    setVec(DoubleFeatureData, StringData, "AP1_peak", AP1_peak);
    return 0;
  } else {
    AP1_peak.push_back(peak_voltage[0]);
  }

  setVec(DoubleFeatureData, StringData, "AP1_peak", AP1_peak);
  return 1;
}

// Amplitude of the second spike
int LibV5::AP2_amp(mapStr2intVec& IntFeatureData,
                   mapStr2doubleVec& DoubleFeatureData,
                   mapStr2Str& StringData) {
  int retVal, nSize;
  retVal = CheckInMap(DoubleFeatureData, StringData, "AP2_amp", nSize);
  if (retVal) return nSize;

  vector<double> AP_amplitudes, AP2_amp;
  retVal = getVec(DoubleFeatureData, StringData, "AP_amplitude",
                        AP_amplitudes);
  if (retVal < 2) {
    setVec(DoubleFeatureData, StringData, "AP2_amp", AP2_amp);
    return 0;
  } else {
    AP2_amp.push_back(AP_amplitudes[1]);
  }

  setVec(DoubleFeatureData, StringData, "AP2_amp", AP2_amp);

  return 1;
}

// Peak voltage of the second spike
int LibV5::AP2_peak(mapStr2intVec& IntFeatureData,
                    mapStr2doubleVec& DoubleFeatureData,
                    mapStr2Str& StringData) {
  int retVal, nSize;
  retVal = CheckInMap(DoubleFeatureData, StringData, "AP2_peak", nSize);
  if (retVal) return nSize;

  vector<double> peak_voltage, AP2_peak;
  retVal =
      getVec(DoubleFeatureData, StringData, "peak_voltage", peak_voltage);
  if (retVal < 2) {
    setVec(DoubleFeatureData, StringData, "AP2_peak", AP2_peak);
    return 0;
  } else {
    AP2_peak.push_back(peak_voltage[1]);
  }

  setVec(DoubleFeatureData, StringData, "AP2_peak", AP2_peak);
  return 1;
}

// Difference amplitude of the second to first spike
int LibV5::AP2_AP1_diff(mapStr2intVec& IntFeatureData,
                        mapStr2doubleVec& DoubleFeatureData,
                        mapStr2Str& StringData) {
  int retVal, nSize;
  retVal =
      CheckInMap(DoubleFeatureData, StringData, "AP2_AP1_diff", nSize);
  if (retVal) return nSize;

  vector<double> AP_amplitudes, AP2_AP1_diff;
  retVal = getVec(DoubleFeatureData, StringData, "AP_amplitude",
                        AP_amplitudes);
  if (retVal < 2) {
    setVec(DoubleFeatureData, StringData, "AP2_AP1_diff", AP2_AP1_diff);
    return 0;
  } else {
    AP2_AP1_diff.push_back(AP_amplitudes[1] - AP_amplitudes[0]);
  }

  setVec(DoubleFeatureData, StringData, "AP2_AP1_diff", AP2_AP1_diff);

  return 1;
}

// Difference peak_amp of the second to first spike
int LibV5::AP2_AP1_peak_diff(mapStr2intVec& IntFeatureData,
                             mapStr2doubleVec& DoubleFeatureData,
                             mapStr2Str& StringData) {
  int retVal, nSize;
  retVal = CheckInMap(DoubleFeatureData, StringData, "AP2_AP1_peak_diff",
                            nSize);
  if (retVal) return nSize;
  vector<double> peak_voltage, AP2_AP1_peak_diff;
  retVal =
      getVec(DoubleFeatureData, StringData, "peak_voltage", peak_voltage);
  if (retVal < 2) {
    setVec(DoubleFeatureData, StringData, "AP2_AP1_peak_diff",
                 AP2_AP1_peak_diff);
    return 0;
  } else {
    AP2_AP1_peak_diff.push_back(peak_voltage[1] - peak_voltage[0]);
  }

  setVec(DoubleFeatureData, StringData, "AP2_AP1_peak_diff",
               AP2_AP1_peak_diff);

  return 1;
}

// Width of the first spike
int LibV5::AP1_width(mapStr2intVec& IntFeatureData,
                     mapStr2doubleVec& DoubleFeatureData,
                     mapStr2Str& StringData) {
  int retVal, nSize;
  retVal = CheckInMap(DoubleFeatureData, StringData, "AP1_width", nSize);
  if (retVal) return nSize;
  vector<double> spike_half_width, AP1_width;
  retVal = getVec(DoubleFeatureData, StringData, "spike_half_width",
                        spike_half_width);
  if (retVal < 1) {
    setVec(DoubleFeatureData, StringData, "AP1_width", AP1_width);
    return 0;
  } else {
    AP1_width.push_back(spike_half_width[0]);
  }

  setVec(DoubleFeatureData, StringData, "AP1_width", AP1_width);

  return 1;
}

// Width of the second spike
int LibV5::AP2_width(mapStr2intVec& IntFeatureData,
                     mapStr2doubleVec& DoubleFeatureData,
                     mapStr2Str& StringData) {
  int retVal, nSize;
  retVal = CheckInMap(DoubleFeatureData, StringData, "AP2_width", nSize);
  if (retVal) return nSize;

  vector<double> spike_half_width, AP2_width;
  retVal = getVec(DoubleFeatureData, StringData, "spike_half_width",
                        spike_half_width);
  if (retVal < 2) {
    setVec(DoubleFeatureData, StringData, "AP2_width", AP2_width);
    return 0;
  } else {
    AP2_width.push_back(spike_half_width[1]);
  }

  setVec(DoubleFeatureData, StringData, "AP2_width", AP2_width);

  return 1;
}

// Width of the last spike
int LibV5::APlast_width(mapStr2intVec& IntFeatureData,
                      mapStr2doubleVec& DoubleFeatureData,
                      mapStr2Str& StringData) {
  int retVal, nSize;
  retVal = CheckInMap(DoubleFeatureData, StringData, "APlast_width", nSize);
  if (retVal) return nSize;
  vector<double> spike_half_width, APlast_width;

  int spike_half_width_size = 
      getVec(DoubleFeatureData, StringData, "spike_half_width",
                        spike_half_width);
  
  if (spike_half_width_size < 1) {
    GErrorStr +=
        "\nError: At least one spike is needed for APlast_width.\n";  
    return -1;
  } else {
    APlast_width.push_back(spike_half_width[spike_half_width_size - 1]);
  }

  setVec(DoubleFeatureData, StringData, "APlast_width", APlast_width);
  return 1;
}


static int __AHP_time_from_peak(const vector<double>& t,
                                const vector<int>& peakIndices,
                                const vector<int>& minAHPIndices,
                                vector<double>& ahpTimeFromPeak) {
  for (size_t i = 0; i < peakIndices.size() && i < minAHPIndices.size();
       i++) {
    ahpTimeFromPeak.push_back(t[minAHPIndices[i]] - t[peakIndices[i]]);
  }
  return ahpTimeFromPeak.size();
}

int LibV5::AHP_time_from_peak(mapStr2intVec& IntFeatureData,
                              mapStr2doubleVec& DoubleFeatureData,
                              mapStr2Str& StringData) {
  int retval;
  int nsize;
  retval = CheckInMap(DoubleFeatureData, StringData, "AHP_time_from_peak",
                            nsize);
  if (retval) {
    return nsize;
  }

  vector<double> T;
  retval = getVec(DoubleFeatureData, StringData, "T", T);
  if (retval < 0) return -1;

  vector<int> peakIndices;
  retval = getVec(IntFeatureData, StringData, "peak_indices", peakIndices);
  if (retval < 0) return -1;

  vector<int> minAHPIndices;
  retval =
      getVec(IntFeatureData, StringData, "min_AHP_indices", minAHPIndices);
  if (retval < 0) return -1;

  vector<double> ahpTimeFromPeak;
  retval = __AHP_time_from_peak(T, peakIndices, minAHPIndices, ahpTimeFromPeak);
  if (retval >= 0) {
    setVec(DoubleFeatureData, StringData, "AHP_time_from_peak",
                 ahpTimeFromPeak);
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

int LibV5::AHP_depth_from_peak(mapStr2intVec& IntFeatureData,
                               mapStr2doubleVec& DoubleFeatureData,
                               mapStr2Str& StringData) {
  int retval;
  int nsize;
  retval = CheckInMap(DoubleFeatureData, StringData,
                            "AHP_depth_from_peak", nsize);

  if (retval) {
    return nsize;
  }

  vector<double> V;
  retval = getVec(DoubleFeatureData, StringData, "V", V);
  if (retval < 0) return -1;

  vector<int> peakIndices;
  retval = getVec(IntFeatureData, StringData, "peak_indices", peakIndices);
  if (retval < 0) return -1;

  vector<int> minAHPIndices;
  retval =
      getVec(IntFeatureData, StringData, "min_AHP_indices", minAHPIndices);
  if (retval < 0) return -1;

  vector<double> ahpDepthFromPeak;
  retval =
      __AHP_depth_from_peak(V, peakIndices, minAHPIndices, ahpDepthFromPeak);
  if (retval >= 0) {
    setVec(DoubleFeatureData, StringData, "AHP_depth_from_peak",
                 ahpDepthFromPeak);
  }
  return retval;
}

// AHP_depth_from_peak of first spike
int LibV5::AHP1_depth_from_peak(mapStr2intVec& IntFeatureData,
                                mapStr2doubleVec& DoubleFeatureData,
                                mapStr2Str& StringData) {
  int retVal, nSize;
  retVal = CheckInMap(DoubleFeatureData, StringData,
                            "AHP1_depth_from_peak", nSize);
  if (retVal) return nSize;

  vector<double> ahpDepthFromPeak, ahp1DepthFromPeak;
  retVal = getVec(DoubleFeatureData, StringData, "AHP_depth_from_peak",
                        ahpDepthFromPeak);
  if (retVal < 1) {
    setVec(DoubleFeatureData, StringData, "AHP1_depth_from_peak",
                 ahp1DepthFromPeak);
    return 0;
  } else {
    ahp1DepthFromPeak.push_back(ahpDepthFromPeak[0]);
  }

  setVec(DoubleFeatureData, StringData, "AHP1_depth_from_peak",
               ahp1DepthFromPeak);

  return 1;
}

// AHP_depth_from_peak of second spike
int LibV5::AHP2_depth_from_peak(mapStr2intVec& IntFeatureData,
                                mapStr2doubleVec& DoubleFeatureData,
                                mapStr2Str& StringData) {
  int retVal, nSize;
  retVal = CheckInMap(DoubleFeatureData, StringData,
                            "AHP2_depth_from_peak", nSize);
  if (retVal) return nSize;

  vector<double> ahpDepthFromPeak, ahp2DepthFromPeak;
  retVal = getVec(DoubleFeatureData, StringData, "AHP_depth_from_peak",
                        ahpDepthFromPeak);
  if (retVal < 2) {
    setVec(DoubleFeatureData, StringData, "AHP2_depth_from_peak",
                 ahp2DepthFromPeak);
    return 0;
  } else {
    ahp2DepthFromPeak.push_back(ahpDepthFromPeak[1]);
  }

  setVec(DoubleFeatureData, StringData, "AHP2_depth_from_peak",
               ahp2DepthFromPeak);

  return 1;
}

// spike width at spike start
static int __AP_begin_width(const vector<double>& t, const vector<double>& v,
                            const vector<int>& AP_begin_indices,
                            const vector<int>& min_ahp_indices,
                            vector<double>& AP_begin_width) {
  // int start_index = distance(t.begin(), find_if(t.begin(), t.end(),
  // std::bind2nd(std::greater_equal<double>(), stim_start)));
  /// vector<int> min_ahp_indices_plus(min_ahp_indices.size() + 1, start_index);
  // copy(min_ahp_indices.begin(), min_ahp_indices.end(),
  // min_ahp_indices_plus.begin());
  if (AP_begin_indices.size() < min_ahp_indices.size()) return -1;
  for (size_t i = 0; i < min_ahp_indices.size(); i++) {
    double v_start = v[AP_begin_indices[i]];
    // interpolate this one time step where the voltage is close to v_start in
    // the falling edge
    int rise_index = AP_begin_indices[i];
    int fall_index = distance(
        v.begin(),
        find_if(v.begin() + rise_index + 1, v.begin() + min_ahp_indices[i],
                std::bind2nd(std::less_equal<double>(), v_start)));
    // v_dev = v_start - v[fall_index];
    // delta_v = v[fall_index] - v[fall_index - 1];
    // delta_t = t[fall_index] - t[fall_index - 1];
    // t_dev_fall = delta_t * v_dev / delta_v;
    // printf("%f %f\n",delta_v, t_dev_fall);
    AP_begin_width.push_back(t[fall_index] - t[rise_index] /*+ t_dev_fall*/);
  }
  return AP_begin_width.size();
}

int LibV5::AP_begin_width(mapStr2intVec& IntFeatureData,
                          mapStr2doubleVec& DoubleFeatureData,
                          mapStr2Str& StringData) {
  int retVal, nSize;
  retVal =
      CheckInMap(DoubleFeatureData, StringData, "AP_begin_width", nSize);
  if (retVal) return nSize;

  vector<int> AP_begin_indices, minAHPIndex;
  vector<double> V, t, dv1, dv2, AP_begin_width;
  vector<double> stim_start;
  retVal = getVec(DoubleFeatureData, StringData, "V", V);
  if (retVal < 0) return -1;
  retVal = getVec(DoubleFeatureData, StringData, "T", t);
  if (retVal < 0) return -1;
  retVal =
      getVec(IntFeatureData, StringData, "min_AHP_indices", minAHPIndex);
  if (retVal < 0) return -1;
  retVal = getVec(IntFeatureData, StringData, "AP_begin_indices",
                     AP_begin_indices);
  if (retVal < 0) return -1;
  if (AP_begin_indices.size() < 1) {
    GErrorStr +=
        "\nError: At least one spike is needed for spikewidth calculation.\n";
    return -1;
  }
  // Take derivative of voltage from 1st AHPmin to the peak of the spike
  // Using Central difference derivative vec1[i] = ((vec[i+1]+vec[i-1])/2)/dx
  retVal =
      __AP_begin_width(t, V, AP_begin_indices, minAHPIndex, AP_begin_width);
  if (retVal >= 0) {
    setVec(DoubleFeatureData, StringData, "AP_begin_width",
                 AP_begin_width);
  }
  return retVal;
}

static int __AP_begin_time(const vector<double>& t, const vector<double>& v,
                           const vector<int>& AP_begin_indices,
                           vector<double>& AP_begin_time) {
  for (size_t i = 0; i < AP_begin_indices.size(); i++) {
    AP_begin_time.push_back(t[AP_begin_indices[i]]);
  }
  return AP_begin_time.size();
}

int LibV5::AP_begin_time(mapStr2intVec& IntFeatureData,
                         mapStr2doubleVec& DoubleFeatureData,
                         mapStr2Str& StringData) {
  int retVal, nSize;
  retVal =
      CheckInMap(DoubleFeatureData, StringData, "AP_begin_time", nSize);
  if (retVal) return nSize;

  vector<int> AP_begin_indices;
  vector<double> V, t, AP_begin_time;
  retVal = getVec(DoubleFeatureData, StringData, "V", V);
  if (retVal < 0) return -1;
  retVal = getVec(DoubleFeatureData, StringData, "T", t);
  if (retVal < 0) return -1;
  retVal = getVec(IntFeatureData, StringData, "AP_begin_indices",
                     AP_begin_indices);
  if (retVal < 0) return -1;

  retVal = __AP_begin_time(t, V, AP_begin_indices, AP_begin_time);
  if (retVal >= 0) {
    setVec(DoubleFeatureData, StringData, "AP_begin_time", AP_begin_time);
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

int LibV5::AP_begin_voltage(mapStr2intVec& IntFeatureData,
                            mapStr2doubleVec& DoubleFeatureData,
                            mapStr2Str& StringData) {
  int retVal, nSize;
  retVal = CheckInMap(DoubleFeatureData, StringData, "AP_begin_voltage",
                            nSize);
  if (retVal) return nSize;

  vector<int> AP_begin_indices;
  vector<double> V, t, AP_begin_voltage;
  retVal = getVec(DoubleFeatureData, StringData, "V", V);
  if (retVal < 0) return -1;
  retVal = getVec(DoubleFeatureData, StringData, "T", t);
  if (retVal < 0) return -1;
  retVal = getVec(IntFeatureData, StringData, "AP_begin_indices",
                     AP_begin_indices);
  if (retVal < 0) return -1;

  retVal = __AP_begin_voltage(t, V, AP_begin_indices, AP_begin_voltage);
  if (retVal >= 0) {
    setVec(DoubleFeatureData, StringData, "AP_begin_voltage",
                 AP_begin_voltage);
  }
  return retVal;
}

int LibV5::AP1_begin_voltage(mapStr2intVec& IntFeatureData,
                             mapStr2doubleVec& DoubleFeatureData,
                             mapStr2Str& StringData) {
  int retVal, nSize;
  retVal = CheckInMap(DoubleFeatureData, StringData, "AP1_begin_voltage",
                            nSize);
  if (retVal) return nSize;

  vector<double> AP_begin_voltage, AP1_begin_voltage;
  retVal = getVec(DoubleFeatureData, StringData, "AP_begin_voltage",
                        AP_begin_voltage);
  if (retVal < 1) {
    setVec(DoubleFeatureData, StringData, "AP1_begin_voltage",
                 AP1_begin_voltage);
    return 0;
  } else {
    AP1_begin_voltage.push_back(AP_begin_voltage[0]);
    setVec(DoubleFeatureData, StringData, "AP1_begin_voltage",
                 AP1_begin_voltage);
    return 1;
  }
}

int LibV5::AP2_begin_voltage(mapStr2intVec& IntFeatureData,
                             mapStr2doubleVec& DoubleFeatureData,
                             mapStr2Str& StringData) {
  int retVal, nSize;
  retVal = CheckInMap(DoubleFeatureData, StringData, "AP2_begin_voltage",
                            nSize);
  if (retVal) return nSize;

  vector<double> AP_begin_voltage, AP2_begin_voltage;
  retVal = getVec(DoubleFeatureData, StringData, "AP_begin_voltage",
                        AP_begin_voltage);
  if (retVal < 2) {
    setVec(DoubleFeatureData, StringData, "AP2_begin_voltage",
                 AP2_begin_voltage);
    return 0;
  } else {
    AP2_begin_voltage.push_back(AP_begin_voltage[1]);
    setVec(DoubleFeatureData, StringData, "AP2_begin_voltage",
                 AP2_begin_voltage);
    return 1;
  }
}

int LibV5::AP1_begin_width(mapStr2intVec& IntFeatureData,
                           mapStr2doubleVec& DoubleFeatureData,
                           mapStr2Str& StringData) {
  int retVal, nSize;
  retVal =
      CheckInMap(DoubleFeatureData, StringData, "AP1_begin_width", nSize);
  if (retVal) return nSize;

  vector<double> AP_begin_width, AP1_begin_width;
  retVal = getVec(DoubleFeatureData, StringData, "AP_begin_width",
                        AP_begin_width);
  if (retVal < 1) {
    setVec(DoubleFeatureData, StringData, "AP1_begin_width",
                 AP1_begin_width);
    return 0;
  } else {
    AP1_begin_width.push_back(AP_begin_width[0]);
    setVec(DoubleFeatureData, StringData, "AP1_begin_width",
                 AP1_begin_width);
    return 1;
  }
}

int LibV5::AP2_begin_width(mapStr2intVec& IntFeatureData,
                           mapStr2doubleVec& DoubleFeatureData,
                           mapStr2Str& StringData) {
  int retVal, nSize;
  retVal =
      CheckInMap(DoubleFeatureData, StringData, "AP2_begin_width", nSize);
  if (retVal) return nSize;

  vector<double> AP_begin_width, AP2_begin_width;
  retVal = getVec(DoubleFeatureData, StringData, "AP_begin_width",
                        AP_begin_width);
  if (retVal < 2) {
    setVec(DoubleFeatureData, StringData, "AP2_begin_width",
                 AP2_begin_width);
    return 0;
  } else {
    AP2_begin_width.push_back(AP_begin_width[1]);
    setVec(DoubleFeatureData, StringData, "AP2_begin_width",
                 AP2_begin_width);
    return 1;
  }
}

// Difference amplitude of the second to first spike
int LibV5::AP2_AP1_begin_width_diff(mapStr2intVec& IntFeatureData,
                                    mapStr2doubleVec& DoubleFeatureData,
                                    mapStr2Str& StringData) {
  int retVal, nSize;
  retVal = CheckInMap(DoubleFeatureData, StringData,
                            "AP2_AP1_begin_width_diff", nSize);
  if (retVal) return nSize;

  vector<double> AP_begin_widths, AP2_AP1_begin_width_diff;
  retVal = getVec(DoubleFeatureData, StringData, "AP_begin_width",
                        AP_begin_widths);
  if (retVal < 2) {
    setVec(DoubleFeatureData, StringData, "AP2_AP1_begin_width_diff",
                 AP2_AP1_begin_width_diff);
    return 0;
  } else {
    AP2_AP1_begin_width_diff.push_back(AP_begin_widths[1] - AP_begin_widths[0]);
  }

  setVec(DoubleFeatureData, StringData, "AP2_AP1_begin_width_diff",
               AP2_AP1_begin_width_diff);

  return 1;
}

static int __voltage_deflection_begin(const vector<double>& v,
                                      const vector<double>& t, double stimStart,
                                      double stimEnd, vector<double>& vd) {
  double deflection_range_percentage = 0.10;
  double range_begin =
      stimStart + (stimEnd - stimStart) * (deflection_range_percentage / 2);
  double range_stop =
      range_begin + (stimEnd - stimStart) * (deflection_range_percentage);
  double base = 0.;
  int base_size = 0;
  for (size_t i = 0; i < t.size(); i++) {
    if (t[i] < stimStart) {
      base += v[i];
      base_size++;
    } else {
      break;
    }
  }
  base /= base_size;
  double volt = 0;
  int volt_size = 0;
  for (size_t i = 0; i < t.size(); i++) {
    if (t[i] > range_stop) {
      break;
    }
    if (t[i] > range_begin) {
      volt += v[i];
      volt_size++;
    }
  }
  volt /= volt_size;

  vd.push_back(volt - base);
  return 1;
}

int LibV5::voltage_deflection_begin(mapStr2intVec& IntFeatureData,
                                    mapStr2doubleVec& DoubleFeatureData,
                                    mapStr2Str& StringData) {
  int retVal;
  int nSize;
  retVal = CheckInMap(DoubleFeatureData, StringData,
                            "voltage_deflection_begin", nSize);
  if (retVal) return nSize;

  vector<double> v;
  vector<double> t;
  vector<double> stimStart;
  vector<double> stimEnd;
  retVal = getVec(DoubleFeatureData, StringData, "V", v);
  if (retVal < 0) return -1;
  retVal = getVec(DoubleFeatureData, StringData, "T", t);
  if (retVal < 0) return -1;
  retVal = getVec(DoubleFeatureData, StringData, "stim_start", stimStart);
  if (retVal < 0) return -1;
  retVal = getVec(DoubleFeatureData, StringData, "stim_end", stimEnd);
  if (retVal < 0) return -1;

  vector<double> vd;
  retVal = __voltage_deflection_begin(v, t, stimStart[0], stimEnd[0], vd);
  if (retVal >= 0) {
    setVec(DoubleFeatureData, StringData, "voltage_deflection_begin", vd);
  }
  return retVal;
}

// Check if a cell is transiently stuck (i.e. not firing any spikes) at the end
// of
// retval will be -1 if the cell get's stuck, retval will be 1 if the cell
// doesn't get stuck
int LibV5::is_not_stuck(mapStr2intVec& IntFeatureData,
                        mapStr2doubleVec& DoubleFeatureData,
                        mapStr2Str& StringData) {
  int retval;
  int size;
  retval = CheckInMap(IntFeatureData, StringData, "is_not_stuck", size);
  if (retval) return size;

  vector<double> peak_time;
  vector<double> stim_start;
  vector<double> stim_end;
  retval = getVec(DoubleFeatureData, StringData, "peak_time", peak_time);
  if (retval < 0) return -1;
  retval =
      getVec(DoubleFeatureData, StringData, "stim_start", stim_start);
  if (retval < 0) return -1;
  retval = getVec(DoubleFeatureData, StringData, "stim_end", stim_end);
  if (retval < 0) return -1;

  bool stuck = true;
  for (size_t i = 0; i < peak_time.size(); i++) {
    if (peak_time[i] > stim_end[0] * 0.5 && peak_time[i] < stim_end[0]) {
      stuck = false;
      break;
    }
  }
  vector<int> tc;
  if (!stuck) {
    tc.push_back(1);
    setVec(IntFeatureData, StringData, "is_not_stuck", tc);
    return tc.size();
  } else {
    return -1;
  }
}

// The mean voltage after the stimulus in (stim_end + 25%*end_period, stim_end +
// 75%*end_period)
int LibV5::voltage_after_stim(mapStr2intVec& IntFeatureData,
                              mapStr2doubleVec& DoubleFeatureData,
                              mapStr2Str& StringData) {
  int retVal, nSize;
  retVal = CheckInMap(DoubleFeatureData, StringData, "voltage_after_stim",
                            nSize);
  if (retVal) return nSize;

  vector<double> v, t, stimEnd, vRest;
  double startTime, endTime;
  retVal = getVec(DoubleFeatureData, StringData, "V", v);
  if (retVal < 0) return -1;
  retVal = getVec(DoubleFeatureData, StringData, "T", t);
  if (retVal < 0) return -1;
  retVal = getVec(DoubleFeatureData, StringData, "stim_end", stimEnd);
  if (retVal < 0) return -1;

  startTime = stimEnd[0] + (t[t.size() - 1] - stimEnd[0]) * .25;
  endTime = stimEnd[0] + (t[t.size() - 1] - stimEnd[0]) * .75;
  int nCount = 0;
  double vSum = 0;
  // calculte the mean of voltage between startTime and endTime
  for (size_t i = 0; i < t.size(); i++) {
    if (t[i] >= startTime) {
      vSum = vSum + v[i];
      nCount++;
    }
    if (t[i] > endTime) break;
  }
  if (nCount == 0) return -1;
  vRest.push_back(vSum / nCount);
  setVec(DoubleFeatureData, StringData, "voltage_after_stim", vRest);
  return 1;
}

int LibV5::mean_AP_amplitude(mapStr2intVec& IntFeatureData,
                             mapStr2doubleVec& DoubleFeatureData,
                             mapStr2Str& StringData) {
  int retVal, nSize;
  retVal = CheckInMap(DoubleFeatureData, StringData, "mean_AP_amplitude",
                            nSize);
  if (retVal > 0) return nSize;

  vector<double> AP_amplitude;
  retVal =
      getVec(DoubleFeatureData, StringData, "AP_amplitude", AP_amplitude);

  if (retVal < 0) {
    GErrorStr += "Error calculating AP_amplitude for mean_AP_amplitude";
    return -1;
  } else if (retVal == 0) {
    GErrorStr += "No spikes found when calculating mean_AP_amplitude";
    return -1;
  } else if (AP_amplitude.size() == 0) {
    GErrorStr += "No spikes found when calculating mean_AP_amplitude";
    return -1;
  }

  vector<double> mean_AP_amplitude;
  double mean_amp = 0.0;

  for (size_t i = 0; i < AP_amplitude.size(); i++) {
    mean_amp += AP_amplitude[i];
  }

  mean_amp /= AP_amplitude.size();
  mean_AP_amplitude.push_back(mean_amp);

  setVec(DoubleFeatureData, StringData, "mean_AP_amplitude",
               mean_AP_amplitude);

  return mean_AP_amplitude.size();
}

// *** BPAPHeightLoc1 ***
int LibV5::BPAPHeightLoc1(mapStr2intVec& IntFeatureData,
                          mapStr2doubleVec& DoubleFeatureData,
                          mapStr2Str& StringData) {
  int retval;
  int nsize;
  retval =
      CheckInMap(DoubleFeatureData, StringData, "BPAPHeightLoc1", nsize);
  if (retval) return nsize;

  vector<double> peakvoltage;
  retval = getDoubleParam(DoubleFeatureData, "peak_voltage;location_dend1",
                          peakvoltage);
  // one spike required
  if (retval <= 0) return -1;

  // voltage base
  vector<double> vb_dend;
  retval =
      getDoubleParam(DoubleFeatureData, "voltage_base;location_dend1", vb_dend);
  if (retval <= 0) return -1;

  vector<double> v_dend;
  retval = getDoubleParam(DoubleFeatureData, "V;location_dend1", v_dend);
  if (retval <= 0) return -1;
  vector<double> bpapheight;

  // bpapheight.push_back(*max_element(v_dend.begin(), v_dend.end()) -
  // vb_dend[0]);
  for (size_t i = 0; i < peakvoltage.size(); i++) {
    bpapheight.push_back(peakvoltage[i] - vb_dend[0]);
    // printf("peak voltage: %f, voltage base: %f, height: %f", peakvoltage[i],
    // vb_dend[0], peakvoltage[0] - vb_dend[0]);
  }

  setVec(DoubleFeatureData, StringData, "BPAPHeightLoc1", bpapheight);
  return bpapheight.size();
}
// end of BPAPHeightLoc1

// *** BPAPAmplitudeLoc1 ***
int LibV5::BPAPAmplitudeLoc1(mapStr2intVec& IntFeatureData,
                             mapStr2doubleVec& DoubleFeatureData,
                             mapStr2Str& StringData) {
  int retval;
  int nsize;
  retval = CheckInMap(DoubleFeatureData, StringData, "BPAPAmplitudeLoc1",
                            nsize);
  if (retval) return nsize;

  vector<double> peakvoltage;
  retval = getDoubleParam(DoubleFeatureData, "peak_voltage;location_dend1",
                          peakvoltage);
  // one spike required
  if (retval <= 0) return -1;

  // voltage base
  vector<double> ap_begin_voltage_dend;
  retval = getDoubleParam(DoubleFeatureData, "AP_begin_voltage;location_dend1",
                          ap_begin_voltage_dend);
  if (retval <= 0) return -1;

  if (peakvoltage.size() > ap_begin_voltage_dend.size()) {
    GErrorStr += "More peakvoltage entries than AP begin voltages";
    return -1;
  }

  vector<double> bpapamplitude;

  for (size_t i = 0; i < peakvoltage.size(); i++) {
    bpapamplitude.push_back(peakvoltage[i] - ap_begin_voltage_dend[i]);
  }

  setVec(DoubleFeatureData, StringData, "BPAPAmplitudeLoc1",
               bpapamplitude);
  return bpapamplitude.size();
}
// end of BPAPAmplitudeLoc1

// *** BPAPAmplitudeLoc2 ***
int LibV5::BPAPAmplitudeLoc2(mapStr2intVec& IntFeatureData,
                             mapStr2doubleVec& DoubleFeatureData,
                             mapStr2Str& StringData) {
  int retval;
  int nsize;
  retval = CheckInMap(DoubleFeatureData, StringData, "BPAPAmplitudeLoc2",
                            nsize);
  if (retval) return nsize;

  vector<double> peakvoltage;
  retval = getDoubleParam(DoubleFeatureData, "peak_voltage;location_dend2",
                          peakvoltage);
  // one spike required
  if (retval <= 0) return -1;

  // voltage base
  vector<double> ap_begin_voltage_dend;
  retval = getDoubleParam(DoubleFeatureData, "AP_begin_voltage;location_dend2",
                          ap_begin_voltage_dend);
  if (retval <= 0) return -1;

  if (peakvoltage.size() > ap_begin_voltage_dend.size()) {
    GErrorStr += "More peakvoltage entries than AP begin voltages";
    return -1;
  }

  vector<double> bpapamplitude;

  for (size_t i = 0; i < peakvoltage.size(); i++) {
    bpapamplitude.push_back(peakvoltage[i] - ap_begin_voltage_dend[i]);
  }

  setVec(DoubleFeatureData, StringData, "BPAPAmplitudeLoc2",
               bpapamplitude);
  return bpapamplitude.size();
}
// end of BPAPAmplitudeLoc2

// *** BPAPHeightLoc2 ***
int LibV5::BPAPHeightLoc2(mapStr2intVec& IntFeatureData,
                          mapStr2doubleVec& DoubleFeatureData,
                          mapStr2Str& StringData) {
  int retval;
  int nsize;
  retval =
      CheckInMap(DoubleFeatureData, StringData, "BPAPHeightLoc2", nsize);
  if (retval) return nsize;

  vector<double> peakvoltage;
  retval = getDoubleParam(DoubleFeatureData, "peak_voltage;location_dend2",
                          peakvoltage);
  // one spike required
  if (retval <= 0) return -1;

  // voltage base
  vector<double> vb_dend;
  retval =
      getDoubleParam(DoubleFeatureData, "voltage_base;location_dend2", vb_dend);
  if (retval <= 0) return -1;

  vector<double> v_dend;
  retval = getDoubleParam(DoubleFeatureData, "V;location_dend2", v_dend);
  if (retval <= 0) return -1;
  vector<double> bpapheight;

  // bpapheight.push_back(*max_element(v_dend.begin(), v_dend.end()) -
  // vb_dend[0]);
  for (size_t i = 0; i < peakvoltage.size(); i++) {
    bpapheight.push_back(peakvoltage[i] - vb_dend[0]);
    // printf("peak voltage: %f, voltage base: %f, height: %f", peakvoltage[i],
    // vb_dend[0], peakvoltage[0] - vb_dend[0]);
  }
  setVec(DoubleFeatureData, StringData, "BPAPHeightLoc2", bpapheight);
  return bpapheight.size();
}
// end of BPAPHeightLoc2

// *** check_AISInitiation ***
int LibV5::check_AISInitiation(mapStr2intVec& IntFeatureData,
                               mapStr2doubleVec& DoubleFeatureData,
                               mapStr2Str& StringData) {
  int retval;
  int nsize;
  // printf("Started check_AISInitation\n");

  retval = CheckInMap(DoubleFeatureData, StringData,
                            "check_AISInitiation", nsize);
  if (retval) return nsize;

  // printf("Starting AP_begin_time\n");

  vector<double> apBeginSoma;
  retval = getDoubleParam(DoubleFeatureData, "AP_begin_time", apBeginSoma);
  if (retval <= 0) {
    printf("Error calculating AP_begin_time\n");
    return -1;
  }

  vector<double> apBeginAIS;
  retval = getDoubleParam(DoubleFeatureData, "AP_begin_time;location_AIS",
                          apBeginAIS);
  if (retval <= 0) {
    printf("Error calculating AP_begin_time\n");
    return -1;
  }

  // printf("Calculated AP_begin_time\n");

  // printf("%d, %d\n", apBeginSoma.size(), apBeginAIS.size()); fflush(stdout);

  // Not the same amount of spike in soma and AIS
  if (apBeginSoma.size() != apBeginAIS.size()) {
    GErrorStr += "\nNot the same amount of spikes in soma and AIS\n";
    return -1;
  }

  // Testing if no spike in the soma start earlier than in the dendrites
  for (size_t i = 0; i < apBeginSoma.size(); i++) {
    /// printf("%f, %f\n", apBeginSoma[i], apBeginAIS[i]); fflush(stdout);
    if (apBeginSoma[i] < apBeginAIS[i]) {
      GErrorStr =
          GErrorStr +
          "\nThere is a spike that initiates in the soma before the axon.\n";
      return -1;
    }
  }
  vector<double> returnvalues;
  returnvalues.push_back(1);
  setVec(DoubleFeatureData, StringData, "check_AISInitiation",
               returnvalues);
  return returnvalues.size();
}
// end of check_AISInitation
//

static int __AP_phaseslope(const vector<double>& v, const vector<double>& t,
                           double stimStart, double stimEnd,
                           vector<double>& ap_phaseslopes, vector<int> apbi,
                           double range) {
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
    range_min_index = apbegin_index - int(range);
    range_max_index = apbegin_index + int(range);
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
int LibV5::AP_phaseslope(mapStr2intVec& IntFeatureData,
                         mapStr2doubleVec& DoubleFeatureData,
                         mapStr2Str& StringData) {
  int retVal;
  int nSize;
  retVal =
      CheckInMap(DoubleFeatureData, StringData, "AP_phaseslope", nSize);
  if (retVal) return nSize;

  vector<double> v;
  retVal = getVec(DoubleFeatureData, StringData, "V", v);
  if (retVal < 0) return -1;

  vector<double> t;
  retVal = getVec(DoubleFeatureData, StringData, "T", t);
  if (retVal < 0) return -1;

  vector<double> stimStart;
  retVal = getVec(DoubleFeatureData, StringData, "stim_start", stimStart);
  if (retVal < 0) return -1;

  vector<double> stimEnd;
  retVal = getVec(DoubleFeatureData, StringData, "stim_end", stimEnd);
  if (retVal < 0) return -1;

  vector<double> range_param;
  retVal = getVec(DoubleFeatureData, StringData, "AP_phaseslope_range",
                        range_param);
  if (retVal < 0) return -1;

  vector<int> apbi;
  retVal = getVec(IntFeatureData, StringData, "AP_begin_indices", apbi);
  if (retVal < 0) return -1;

  vector<double> ap_phaseslopes;
  retVal = __AP_phaseslope(v, t, stimStart[0], stimEnd[0], ap_phaseslopes, apbi,
                           range_param[0]);
  if (retVal >= 0) {
    setVec(DoubleFeatureData, StringData, "AP_phaseslope",
                 ap_phaseslopes);
  }
  return retVal;
}

/// Calculate the slope of the V, dVdt plot at the beginning of every spike
/// (at the point where the derivative crosses the DerivativeThreshold)
int LibV5::AP_phaseslope_AIS(mapStr2intVec& IntFeatureData,
                             mapStr2doubleVec& DoubleFeatureData,
                             mapStr2Str& StringData) {
  int retVal;
  int nSize;
  retVal = CheckInMap(DoubleFeatureData, StringData, "AP_phaseslope_AIS",
                            nSize);
  if (retVal) return nSize;

  vector<double> ap_phaseslopes;
  retVal = getVec(DoubleFeatureData, StringData,
                        "AP_phaseslope;location_AIS", ap_phaseslopes);
  if (retVal < 0) return -1;
  setVec(DoubleFeatureData, StringData, "AP_phaseslope_AIS",
               ap_phaseslopes);
  return retVal;
}

int LibV5::BAC_width(mapStr2intVec& IntFeatureData,
                     mapStr2doubleVec& DoubleFeatureData,
                     mapStr2Str& StringData) {
  int retVal;
  int nSize;
  retVal = CheckInMap(DoubleFeatureData, StringData, "BAC_width", nSize);
  if (retVal) return nSize;

  vector<double> ap_width;
  retVal = getVec(DoubleFeatureData, StringData, "AP_width;location_epsp",
                        ap_width);
  if (retVal < 0) {
    GErrorStr += "\n AP_width calculation failed in BAC_width.\n";
    return -1;
  }
  if (retVal > 1) {
    GErrorStr +=
        "\n More than one spike found a location_epsp for BAC_width.\n";
    return -1;
  }
  setVec(DoubleFeatureData, StringData, "BAC_width", ap_width);

  return retVal;
}

int LibV5::BAC_maximum_voltage(mapStr2intVec& IntFeatureData,
                               mapStr2doubleVec& DoubleFeatureData,
                               mapStr2Str& StringData) {
  int retVal;
  int nSize;
  retVal = CheckInMap(DoubleFeatureData, StringData,
                            "BAC_maximum_voltage", nSize);
  if (retVal) return nSize;

  vector<double> ap_phaseslopes;
  retVal = getVec(DoubleFeatureData, StringData,
                        "maximum_voltage;location_epsp", ap_phaseslopes);
  if (retVal != 1) return -1;

  setVec(DoubleFeatureData, StringData, "BAC_maximum_voltage",
               ap_phaseslopes);
  return retVal;
}

int LibV5::all_ISI_values(mapStr2intVec& IntFeatureData,
                          mapStr2doubleVec& DoubleFeatureData,
                          mapStr2Str& StringData) {
  int retVal, nSize;
  retVal =
      CheckInMap(DoubleFeatureData, StringData, "all_ISI_values", nSize);
  if (retVal) return nSize;

  vector<double> VecISI, pvTime;
  retVal = getVec(DoubleFeatureData, StringData, "peak_time", pvTime);
  if (retVal < 2) {
    GErrorStr += "\n Two spikes required for calculation of all_ISI_values.\n";
    return -1;
  }
  for (size_t i = 1; i < pvTime.size(); i++) {
    VecISI.push_back(pvTime[i] - pvTime[i - 1]);
  }
  setVec(DoubleFeatureData, StringData, "all_ISI_values", VecISI);
  return VecISI.size();
}

// spike amplitude: peak_voltage - voltage_base
int LibV5::AP_amplitude_from_voltagebase(mapStr2intVec& IntFeatureData,
                                         mapStr2doubleVec& DoubleFeatureData,
                                         mapStr2Str& StringData) {
  int retVal, nSize;
  retVal = CheckInMap(DoubleFeatureData, StringData,
                            "AP_amplitude_from_voltagebase", nSize);
  if (retVal > 0) return nSize;

  vector<double> peakvoltage;
  vector<double> voltage_base_vec;
  double voltage_base;
  retVal = getVec(DoubleFeatureData, StringData, "voltage_base",
                        voltage_base_vec);
  if (retVal <= 0) {
    GErrorStr +=
        "Error calculating voltage_base for AP_amplitude_from_voltagebase";
    return -1;
  } else {
    voltage_base = voltage_base_vec[0];
  }
  retVal =
      getVec(DoubleFeatureData, StringData, "peak_voltage", peakvoltage);
  if (retVal <= 0) {
    GErrorStr +=
        "Error calculating peak_voltage for AP_amplitude_from_voltagebase";
    return -1;
  }

  vector<double> apamplitude;
  apamplitude.resize(peakvoltage.size());
  for (size_t i = 0; i < apamplitude.size(); i++) {
    apamplitude[i] = peakvoltage[i] - voltage_base;
  }
  setVec(DoubleFeatureData, StringData, "AP_amplitude_from_voltagebase",
               apamplitude);
  return apamplitude.size();
}

// min_voltage_between_spikes: minimal voltage between consecutive spikes
int LibV5::min_voltage_between_spikes(mapStr2intVec& IntFeatureData,
                                      mapStr2doubleVec& DoubleFeatureData,
                                      mapStr2Str& StringData) {
  int retVal, nSize;
  retVal = CheckInMap(DoubleFeatureData, StringData,
                            "min_voltage_between_spikes", nSize);
  if (retVal > 0) return nSize;

  vector<int> peak_indices;
  vector<double> v;
  vector<double> min_voltage_between_spikes;
  retVal = getVec(IntFeatureData, StringData, "peak_indices", peak_indices);
  if (retVal < 0) {
    GErrorStr +=
        "Error calculating peak_indices for min_voltage_between_spikes";
    return -1;
  } else if (retVal < 2) {
    setVec(DoubleFeatureData, StringData, "min_voltage_between_spikes",
                 min_voltage_between_spikes);
    return 0;
  }

  retVal = getVec(DoubleFeatureData, StringData, "V", v);
  if (retVal <= 0) {
    GErrorStr += "Error getting V for min_voltage_between_spikes";
    return -1;
  }

  for (size_t i = 0; i < peak_indices.size() - 1; i++) {
    min_voltage_between_spikes.push_back(*min_element(
        v.begin() + peak_indices[i], v.begin() + peak_indices[i + 1]));
  }
  setVec(DoubleFeatureData, StringData, "min_voltage_between_spikes",
               min_voltage_between_spikes);
  return min_voltage_between_spikes.size();
}

// return (possibly interpolate) voltage trace
int LibV5::voltage(mapStr2intVec& IntFeatureData,
                   mapStr2doubleVec& DoubleFeatureData,
                   mapStr2Str& StringData) {
  int retVal, nSize;
  retVal = CheckInMap(DoubleFeatureData, StringData, "voltage", nSize);
  if (retVal > 0) return nSize;

  vector<double> v;
  retVal = getVec(DoubleFeatureData, StringData, "V", v);
  if (retVal < 0) {
    GErrorStr += "Error getting V for voltage";
    return -1;
  }

  setVec(DoubleFeatureData, StringData, "voltage", v);
  return v.size();
}

// return (possibly interpolate) current trace
int LibV5::current(mapStr2intVec& IntFeatureData,
                   mapStr2doubleVec& DoubleFeatureData,
                   mapStr2Str& StringData) {
  int retVal, nSize;
  retVal = CheckInMap(DoubleFeatureData, StringData, "current", nSize);
  if (retVal > 0) return nSize;

  vector<double> i;
  retVal = getVec(DoubleFeatureData, StringData, "I", i);
  if (retVal < 0) {
    GErrorStr += "Error getting I for current";
    return -1;
  }

  setVec(DoubleFeatureData, StringData, "current", i);
  return i.size();
}

// return (possibly interpolate) time trace
int LibV5::time(mapStr2intVec& IntFeatureData,
                   mapStr2doubleVec& DoubleFeatureData,
                   mapStr2Str& StringData) {
  int retVal, nSize;
  retVal = CheckInMap(DoubleFeatureData, StringData, "time", nSize);
  if (retVal > 0) return nSize;

  vector<double> t;
  retVal = getVec(DoubleFeatureData, StringData, "T", t);
  if (retVal < 0) {
    GErrorStr += "Error getting T for voltage";
    return -1;
  }

  setVec(DoubleFeatureData, StringData, "time", t);
  return t.size();
}

// *** The average voltage during the last 90% of the stimulus duration. ***
int LibV5::steady_state_voltage_stimend(mapStr2intVec& IntFeatureData,
                                        mapStr2doubleVec& DoubleFeatureData,
                                        mapStr2Str& StringData) {
  int retVal;
  int nSize;
  vector<double> t, v, stimEnd, stimStart, ssv;

  retVal = CheckInMap(DoubleFeatureData, StringData,
                            "steady_state_voltage_stimend", nSize);
  if (retVal) {
    return nSize;
  }

  retVal = getVec(DoubleFeatureData, StringData, "V", v);
  if (retVal < 0) return -1;
  retVal = getVec(DoubleFeatureData, StringData, "T", t);
  if (retVal < 0) return -1;
  retVal = getVec(DoubleFeatureData, StringData, "stim_end", stimEnd);
  if (retVal < 0) return -1;
  retVal = getVec(DoubleFeatureData, StringData, "stim_start", stimStart);
  if (retVal < 0) return -1;

  double start_time = stimEnd[0] - 0.1 * (stimEnd[0] - stimStart[0]);
  size_t start_index =
      distance(t.begin(),
               find_if(t.begin(), t.end(),
                       std::bind2nd(std::greater_equal<double>(), start_time)));
  size_t stop_index =
      distance(t.begin(),
               find_if(t.begin(), t.end(),
                       std::bind2nd(std::greater_equal<double>(), stimEnd[0])));

  size_t mean_size = 0;
  double mean = 0.0;
  for (size_t i = start_index; i < stop_index; i++) {
    mean += v[i];
    mean_size++;
  }

  // Check for division by zero
  if (mean_size == 0) {
    return -1;
  } else {
    mean /= mean_size;
    ssv.push_back(mean);

    setVec(DoubleFeatureData, StringData, "steady_state_voltage_stimend",
                 ssv);

    return 1;
  }
}

int LibV5::voltage_base(mapStr2intVec& IntFeatureData,
                        mapStr2doubleVec& DoubleFeatureData,
                        mapStr2Str& StringData) {
  int retVal, nSize;
  retVal =
      CheckInMap(DoubleFeatureData, StringData, "voltage_base", nSize);
  if (retVal) return nSize;

  vector<double> v, t, stimStart, vRest, vb_start_perc_vec, vb_end_perc_vec;
  double startTime, endTime, vb_start_perc, vb_end_perc;
  retVal = getVec(DoubleFeatureData, StringData, "V", v);
  if (retVal < 0) return -1;
  retVal = getVec(DoubleFeatureData, StringData, "T", t);
  if (retVal < 0) return -1;
  retVal = getVec(DoubleFeatureData, StringData, "stim_start", stimStart);
  if (retVal < 0) return -1;
  retVal = getVec(DoubleFeatureData, StringData,
                        "voltage_base_start_perc", vb_start_perc_vec);
  if (retVal == 1) {
    vb_start_perc = vb_start_perc_vec[0];
  } else {
    vb_start_perc = 0.9;
  }
  retVal = getVec(DoubleFeatureData, StringData, "voltage_base_end_perc",
                        vb_end_perc_vec);
  if (retVal == 1) {
    vb_end_perc = vb_end_perc_vec[0];
  } else {
    vb_end_perc = 1.0;
  }

  startTime = stimStart[0] * vb_start_perc;
  endTime = stimStart[0] * vb_end_perc;

  if (startTime >= endTime) {
    GErrorStr += "\nvoltage_base: startTime >= endTime\n";
    return -1;
  }


  vector<double> precisionThreshold;
  retVal = getDoubleParam(DoubleFeatureData, "precision_threshold",
                          precisionThreshold);
  if (retVal < 0) return -1;

  std::pair<size_t, size_t> time_index = get_time_index(t, startTime, endTime,
                                                        precisionThreshold[0]);

  vector<double> subVector(v.begin()+time_index.first,
                           v.begin()+time_index.second);

  double vBase;
  std::string computation_mode;

  retVal = getStrParam(StringData, "voltage_base_mode", computation_mode);
  if (retVal < 0) return -1;


  try{
    if (computation_mode == "mean")
      vBase = vec_mean(subVector);
    else if (computation_mode == "median")
      vBase = vec_median(subVector);
    else
      throw std::invalid_argument(
        "Undefined computational mode. Only mean and median are enabled");
  }
  catch(std::exception &e) {
    GErrorStr +=
    "\nvoltage_base error:" + std::string(e.what()) + "\n";
    return -1;
    }

  vRest.push_back(vBase);
  setVec(DoubleFeatureData, StringData, "voltage_base", vRest);
  return 1;
}

int LibV5::current_base(mapStr2intVec& IntFeatureData,
                        mapStr2doubleVec& DoubleFeatureData,
                        mapStr2Str& StringData) {

  int retVal, nSize;
  retVal =
      CheckInMap(DoubleFeatureData, StringData, "current_base", nSize);
  if (retVal) return nSize;                      

  vector<double> i, t, stimStart, iRest, cb_start_perc_vec, cb_end_perc_vec;
  double startTime, endTime, cb_start_perc, cb_end_perc;
  retVal = getVec(DoubleFeatureData, StringData, "I", i);
  if (retVal < 0) return -1;
  retVal = getVec(DoubleFeatureData, StringData, "T", t);
  if (retVal < 0) return -1;
  retVal = getVec(DoubleFeatureData, StringData, "stim_start", stimStart);
  if (retVal < 0) return -1;
  retVal = getVec(DoubleFeatureData, StringData,
                        "current_base_start_perc", cb_start_perc_vec);
  if (retVal == 1) {
    cb_start_perc = cb_start_perc_vec[0];
  } else {
    cb_start_perc = 0.9;
  }
  retVal = getVec(DoubleFeatureData, StringData, "current_base_end_perc",
                        cb_end_perc_vec);
  if (retVal == 1) {
    cb_end_perc = cb_end_perc_vec[0];
  } else {
    cb_end_perc = 1.0;
  }

  startTime = stimStart[0] * cb_start_perc;
  endTime = stimStart[0] * cb_end_perc;

  if (startTime >= endTime) {
    GErrorStr += "\ncurrent_base: startTime >= endTime\n";
    return -1;
  }


  vector<double> precisionThreshold;
  retVal = getDoubleParam(DoubleFeatureData, "precision_threshold",
                          precisionThreshold);
  if (retVal < 0) return -1;


  std::pair<size_t, size_t> time_index = get_time_index(t, startTime, endTime,
                                                        precisionThreshold[0]);


  vector<double> subVector(i.begin()+time_index.first,
                           i.begin()+time_index.second);

  double iBase;
  std::string computation_mode;

  retVal = getStrParam(StringData, "current_base_mode", computation_mode);
  if (retVal < 0) return -1;


  try{
    if (computation_mode == "mean")
      iBase = vec_mean(subVector);
    else if (computation_mode == "median")
      iBase = vec_median(subVector);
    else
      throw std::invalid_argument(
        "Undefined computational mode. Only mean and median are enabled");
  }
  catch(std::exception &e) {
    GErrorStr +=
    "\ncurrent_base error:" + std::string(e.what()) + "\n";
    return -1;
    }

  iRest.push_back(iBase);
  setVec(DoubleFeatureData, StringData, "current_base", iRest);
  return 1;

}

size_t get_index(const vector<double>& times, double t) {
  return distance(times.begin(),
                  find_if(times.begin(), times.end(),
                          std::bind2nd(std::greater_equal<double>(), t)));
}

double __decay_time_constant_after_stim(const vector<double>& times,
                                        const vector<double>& voltage,
                                        const double decay_start_after_stim,
                                        const double decay_end_after_stim,
                                        const double stimStart,
                                        const double stimEnd) {
  const size_t stimStartIdx = get_index(times, stimStart);
  const size_t decayStartIdx =
      get_index(times, stimEnd + decay_start_after_stim);

  const size_t decayEndIdx = get_index(times, stimEnd + decay_end_after_stim);

  const double reference = voltage[stimStartIdx];

  vector<double> decayValues(decayEndIdx - decayStartIdx);
  vector<double> decayTimes(decayEndIdx - decayStartIdx);

  for (size_t i = 0; i != decayValues.size(); ++i) {
    const double u0 = std::abs(voltage[decayStartIdx + i] - reference);

    decayValues[i] = log(u0);
    decayTimes[i] = times[decayStartIdx + i];
  }

  if (decayTimes.size() < 1 || decayValues.size() < 1) {
      GErrorStr +=
        "\ndecay_time_constant_after_stim: no data points to calculate this feature\n";
      return -1;
  }
  else {
      linear_fit_result fit;
      fit = slope_straight_line_fit(decayTimes, decayValues);

      const double tau = -1.0 / fit.slope;
      return std::abs(tau);
  }
}

// *** Decay time constant measured during decay after the stimulus***
int LibV5::decay_time_constant_after_stim(mapStr2intVec& IntFeatureData,
                                          mapStr2doubleVec& DoubleFeatureData,
                                          mapStr2Str& StringData) {
  int retVal;
  int nSize;

  retVal = CheckInMap(DoubleFeatureData, StringData,
                            "decay_time_constant_after_stim", nSize);
  if (retVal) {
    return nSize;
  }

  vector<double> voltages;
  retVal = getVec(DoubleFeatureData, StringData, "V", voltages);
  if (retVal < 0) return -1;

  vector<double> times;
  retVal = getVec(DoubleFeatureData, StringData, "T", times);
  if (retVal < 0) return -1;

  vector<double> vect;

  retVal = getVec(DoubleFeatureData, StringData, "stim_end", vect);
  if (retVal != 1) return -1;
  const double stimEnd = vect[0];

  retVal = getVec(DoubleFeatureData, StringData, "stim_start", vect);
  if (retVal != 1) return -1;
  const double stimStart = vect[0];

  double decay_start_after_stim, decay_end_after_stim;
  retVal = getVec(DoubleFeatureData, StringData, "decay_start_after_stim",
                        vect);
  if (retVal == 1) {
    decay_start_after_stim = vect[0];
  } else {
    decay_start_after_stim = 1.0;
  }

  retVal =
      getVec(DoubleFeatureData, StringData, "decay_end_after_stim", vect);
  if (retVal == 1) {
    decay_end_after_stim = vect[0];
  } else {
    decay_end_after_stim = 10.0;
  }

  if (decay_start_after_stim >= decay_end_after_stim) {
    GErrorStr +=
        "Error decay_start_after_stim small larger than decay_end_after_stim";
    return -1;
  }

  const double val = __decay_time_constant_after_stim(
      times, voltages, decay_start_after_stim, decay_end_after_stim, stimStart,
      stimEnd);

  vector<double> dtcas;
  dtcas.push_back(val);
  setVec(DoubleFeatureData, StringData, "decay_time_constant_after_stim",
               dtcas);

  return 1;
}

/// *** Voltage deflection between voltage_base and steady_state_voltage_stimend

int LibV5::voltage_deflection_vb_ssse(mapStr2intVec& IntFeatureData,
                                      mapStr2doubleVec& DoubleFeatureData,
                                      mapStr2Str& StringData) {
  int retVal;
  int nSize;
  retVal = CheckInMap(DoubleFeatureData, StringData,
                            "voltage_deflection_vb_ssse", nSize);
  if (retVal) return nSize;

  vector<double> voltage_base;
  retVal =
      getVec(DoubleFeatureData, StringData, "voltage_base", voltage_base);
  if (retVal <= 0) return -1;

  vector<double> steady_state_voltage_stimend;
  retVal = getVec(DoubleFeatureData, StringData,
                        "steady_state_voltage_stimend",
                        steady_state_voltage_stimend);
  if (retVal <= 0) return -1;

  vector<double> voltage_deflection_vb_ssse;

  voltage_deflection_vb_ssse.push_back(steady_state_voltage_stimend[0] -
                                       voltage_base[0]);

  setVec(DoubleFeatureData, StringData, "voltage_deflection_vb_ssse",
               voltage_deflection_vb_ssse);
  retVal = 1;

  return retVal;
}

// *** ohmic input resistance based on voltage_deflection_vb_ssse***

int LibV5::ohmic_input_resistance_vb_ssse(mapStr2intVec& IntFeatureData,
                                          mapStr2doubleVec& DoubleFeatureData,
                                          mapStr2Str& StringData) {
  int retVal;
  int nSize;
  retVal = CheckInMap(DoubleFeatureData, StringData,
                            "ohmic_input_resistance_vb_ssse", nSize);
  if (retVal) return nSize;

  vector<double> voltage_deflection_vb_ssse;
  retVal =
      getVec(DoubleFeatureData, StringData, "voltage_deflection_vb_ssse",
                   voltage_deflection_vb_ssse);
  if (retVal <= 0) return -1;
  vector<double> stimulus_current;
  retVal = getDoubleParam(DoubleFeatureData, "stimulus_current", 
                          stimulus_current);

  if (retVal <= 0) return -1;
  vector<double> ohmic_input_resistance_vb_ssse;

  ohmic_input_resistance_vb_ssse.push_back(voltage_deflection_vb_ssse[0] /
                                           stimulus_current[0]);
  setVec(DoubleFeatureData, StringData, "ohmic_input_resistance_vb_ssse",
               ohmic_input_resistance_vb_ssse);
  retVal = 1;

  return retVal;
}

// *** Diff between maximum voltage during stimulus and voltage_base ***
int LibV5::maximum_voltage_from_voltagebase(mapStr2intVec& IntFeatureData,
                                            mapStr2doubleVec& DoubleFeatureData,
                                            mapStr2Str& StringData) {
  int retVal;
  int nSize;
  retVal = CheckInMap(DoubleFeatureData, StringData,
                            "maximum_voltage_from_voltagebase", nSize);
  if (retVal) return nSize;

  vector<double> maximum_voltage;
  retVal = getVec(DoubleFeatureData, StringData, "maximum_voltage",
                        maximum_voltage);
  if (retVal <= 0) return -1;

  vector<double> voltage_base;
  retVal =
      getVec(DoubleFeatureData, StringData, "voltage_base", voltage_base);
  if (retVal <= 0) return -1;

  vector<double> maximum_voltage_from_voltagebase;

  maximum_voltage_from_voltagebase.push_back(maximum_voltage[0] -
                                             voltage_base[0]);
  setVec(DoubleFeatureData, StringData,
               "maximum_voltage_from_voltagebase",
               maximum_voltage_from_voltagebase);
  retVal = 1;

  return retVal;
}

struct InInterval {
  InInterval(double lower, double upper) : lower(lower), upper(upper) {}
  bool operator()(double value) {
    return ((value >= lower) && (value <= upper));
  }

 private:
  double lower, upper;
};

// *** Spikecount_stimint ***
int LibV5::Spikecount_stimint(mapStr2intVec& IntFeatureData,
                              mapStr2doubleVec& DoubleFeatureData,
                              mapStr2Str& StringData) {
  int retval;
  int nsize;
  int spikecount_stimint_value;
  retval =
      CheckInMap(IntFeatureData, StringData, "Spikecount_stimint", nsize);
  if (retval) {
    return nsize;
  }

  // Get stimulus start
  vector<double> stimstart;
  retval = getVec(DoubleFeatureData, StringData, "stim_start", stimstart);
  if (retval <= 0) {
    GErrorStr += "\nSpikecount_stimint: stim_start not found\n";
    return -1;
  }

  // Get stimulus end
  vector<double> stimend;
  retval = getVec(DoubleFeatureData, StringData, "stim_end", stimend);
  if (retval <= 0) {
    GErrorStr += "\nSpikecount_stimint: stim_start not found\n";
    return -1;
  }

  // Get the times of the peaks
  vector<double> peaktimes;
  retval = getVec(DoubleFeatureData, StringData, "peak_time", peaktimes);

  if (retval < 0) {
    GErrorStr += "\nSpikecount_stimint: peak_time failed\n";
    return -1;
  } else {
    // Get the number of peaks between stim start and end
    spikecount_stimint_value = count_if(peaktimes.begin(), peaktimes.end(),
                                        InInterval(stimstart[0], stimend[0]));

    vector<int> spikecount_stimint(1, spikecount_stimint_value);
    setVec(IntFeatureData, StringData, "Spikecount_stimint",
              spikecount_stimint);
    return 1;
  }
}
// end of Spikecount_stimint

static int __peak_indices(double threshold, vector<double>& V,
                          vector<double>& t, vector<int>& PeakIndex,
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
  if (dnVec.size() == 0) {
    GErrorStr +=
        "\nVoltage never goes below or above threshold in spike detection.\n";
    return 0;
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
        EFEL_ASSERT(itmp >= 0, "peak_time is negative");
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

int LibV5::peak_indices(mapStr2intVec& IntFeatureData,
                        mapStr2doubleVec& DoubleFeatureData,
                        mapStr2Str& StringData) {
  int retVal, nSize;

  retVal = CheckInMap(IntFeatureData, StringData, "peak_indices", nSize);
  if (retVal) return nSize;

  vector<int> PeakIndex, strict_stiminterval_vec;
  vector<double> v, t, threshold, stim_start_vec, stim_end_vec;
  bool strict_stiminterval = false;
  double stim_start = 0.0, stim_end = 0.0;

  retVal = getVec(DoubleFeatureData, StringData, "V", v);
  if (retVal <= 0) {
    return -1;
  }

  retVal = getVec(DoubleFeatureData, StringData, "T", t);
  if (retVal <= 0) {
    return -1;
  }

  retVal = getDoubleParam(DoubleFeatureData, "Threshold", threshold);
  if (retVal <= 0) {
    return -1;
  }

  retVal = getIntParam(IntFeatureData, "strict_stiminterval",
                       strict_stiminterval_vec);
  if (retVal <= 0) {
    strict_stiminterval = false;
  } else {
    strict_stiminterval = bool(strict_stiminterval_vec[0]);
  }

  retVal =
      getVec(DoubleFeatureData, StringData, "stim_start", stim_start_vec);
  if (retVal <= 0) {
    return -1;
  } else {
    stim_start = stim_start_vec[0];
  }

  retVal =
      getVec(DoubleFeatureData, StringData, "stim_end", stim_end_vec);
  if (retVal <= 0) {
    return -1;
  } else {
    stim_end = stim_end_vec[0];
  }

  int retval = __peak_indices(threshold[0], v, t, PeakIndex,
                              strict_stiminterval, stim_start, stim_end);

  if (retval >= 0) {
    setVec(IntFeatureData, StringData, "peak_indices", PeakIndex);
  }

  return retval;
}
int LibV5::sag_amplitude(mapStr2intVec& IntFeatureData,
                                          mapStr2doubleVec& DoubleFeatureData,
                                          mapStr2Str& StringData) {
  int retVal;
  int nSize;
  retVal = CheckInMap(DoubleFeatureData, StringData,
                            "sag_amplitude", nSize);
  if (retVal) return nSize;

  // Get steady_state_voltage_stimend
  vector<double> steady_state_voltage_stimend;
  retVal =
      getVec(DoubleFeatureData, StringData, 
                   "steady_state_voltage_stimend", 
                   steady_state_voltage_stimend);
  if (retVal <= 0) return -1;

  // Get voltage_deflection_stim_ssse  
  double voltage_deflection_stim_ssse = 0.0;
  vector<double> voltage_deflection_vb_ssse_vec;
  retVal =
      getVec(DoubleFeatureData, StringData, 
                   "voltage_deflection_vb_ssse", 
                   voltage_deflection_vb_ssse_vec);
  if (retVal <= 0) {
      return -1;
  } else {
      voltage_deflection_stim_ssse = voltage_deflection_vb_ssse_vec[0];
  }
  
  // Get minimum_voltage
  vector<double> minimum_voltage;
  retVal =
      getVec(DoubleFeatureData, StringData, 
                   "minimum_voltage", 
                   minimum_voltage);
  if (retVal <= 0) return -1;
 
  // Calculate sag_amplitude
  vector<double> sag_amplitude;
  if (voltage_deflection_stim_ssse <= 0) {
      sag_amplitude.push_back(steady_state_voltage_stimend[0] 
              - minimum_voltage[0]);
  } else {
      //In case of positive voltage deflection, return an error
      GErrorStr += "\nsag_amplitude: voltage_deflection is positive\n";
      return -1;
  }
  setVec(DoubleFeatureData, StringData, "sag_amplitude",
               sag_amplitude);
  retVal = 1;
  return retVal;
}

int LibV5::sag_ratio1(mapStr2intVec& IntFeatureData,
                                          mapStr2doubleVec& DoubleFeatureData,
                                          mapStr2Str& StringData) {
  int retVal;
  int nSize;
  retVal = CheckInMap(DoubleFeatureData, StringData,
                            "sag_ratio1", nSize);
  if (retVal) return nSize;

  // Get sag_amplitude
  vector<double> sag_amplitude;
  retVal =
      getVec(DoubleFeatureData, StringData, 
                   "sag_amplitude", 
                   sag_amplitude);
  if (retVal <= 0) return -1;

  // Get voltage_base
  vector<double> voltage_base;
  retVal =
      getVec(DoubleFeatureData, StringData, 
                   "voltage_base", 
                   voltage_base);
  if (retVal <= 0) {return -1;}
  
  // Get minimum_voltage
  vector<double> minimum_voltage;
  retVal =
      getVec(DoubleFeatureData, StringData, 
                   "minimum_voltage", 
                   minimum_voltage);
  if (retVal <= 0) return -1;
 
  // Calculate sag_ratio1
  vector<double> sag_ratio1;
  if (minimum_voltage[0] == voltage_base[0]) {
      GErrorStr += "\nsag_ratio1: voltage_base equals minimum_voltage\n";
      //In case of possible division by zero return error
      return -1;
  } else {
      sag_ratio1.push_back(sag_amplitude[0] / (voltage_base[0] - minimum_voltage[0]));
  }
  setVec(DoubleFeatureData, StringData, "sag_ratio1",
               sag_ratio1);
  retVal = 1;
  return retVal;
}

int LibV5::sag_ratio2(mapStr2intVec& IntFeatureData,
                                          mapStr2doubleVec& DoubleFeatureData,
                                          mapStr2Str& StringData) {
  int retVal;
  int nSize;
  retVal = CheckInMap(DoubleFeatureData, StringData,
                            "sag_ratio2", nSize);
  if (retVal) return nSize;

  // Get voltage_base
  vector<double> voltage_base;
  retVal =
      getVec(DoubleFeatureData, StringData, 
                   "voltage_base", 
                   voltage_base);
  if (retVal <= 0) {return -1;}
  
  // Get minimum_voltage
  vector<double> minimum_voltage;
  retVal =
      getVec(DoubleFeatureData, StringData, 
                   "minimum_voltage", 
                   minimum_voltage);
  if (retVal <= 0) return -1;
 
  // Get steady_state_voltage_stimend
  vector<double> steady_state_voltage_stimend;
  retVal =
      getVec(DoubleFeatureData, StringData, 
                   "steady_state_voltage_stimend", 
                   steady_state_voltage_stimend);
  if (retVal <= 0) return -1;
  
  // Calculate sag_ratio2
  vector<double> sag_ratio2;
  if (minimum_voltage[0] == voltage_base[0]) {
      GErrorStr += "\nsag_ratio2: voltage_base equals minimum_voltage\n";
      //In case of possible division by zero return error
      return -1;
  } else {
      sag_ratio2.push_back((voltage_base[0] - steady_state_voltage_stimend[0]) / (voltage_base[0] - minimum_voltage[0]));
  }
  setVec(DoubleFeatureData, StringData, "sag_ratio2",
               sag_ratio2);
  retVal = 1;
  return retVal;
}

