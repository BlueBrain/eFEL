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

#include "EfelExceptions.h"

using std::distance;
using std::find_if;

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
    size_t isisToRemove = (size_t)((isiValues.size() + 1) * spikeSkipf + .5);

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
  const auto& isivalues = getFeature(DoubleFeatureData, {"ISI_values"});
  vector<double> slope;
  bool semilog = false;
  int retval = __ISI_log_slope(isivalues, slope, false, 0.0, 0, semilog);
  if (retval >= 0) {
    setVec(DoubleFeatureData, StringData, "ISI_log_slope", slope);
    return static_cast<int>(slope.size());
  } else {
    return retval;
  }
}

int LibV5::ISI_semilog_slope(mapStr2intVec& IntFeatureData,
                             mapStr2doubleVec& DoubleFeatureData,
                             mapStr2Str& StringData) {
  const auto& isivalues = getFeature(DoubleFeatureData, {"ISI_values"});
  vector<double> slope;
  bool semilog = true;
  int retval = __ISI_log_slope(isivalues, slope, false, 0.0, 0, semilog);
  if (retval >= 0) {
    setVec(DoubleFeatureData, StringData, "ISI_semilog_slope", slope);
    return static_cast<int>(slope.size());
  } else {
    return retval;
  }
}

int LibV5::ISI_log_slope_skip(mapStr2intVec& IntFeatureData,
                              mapStr2doubleVec& DoubleFeatureData,
                              mapStr2Str& StringData) {
  const auto& isivalues = getFeature(DoubleFeatureData, {"ISI_values"});
  const auto spikeSkipf =
      getFeature(DoubleFeatureData, {"spike_skipf"}).front();
  const auto maxnSpike = getFeature(IntFeatureData, {"max_spike_skip"}).front();

  // Check the validity of spikeSkipf value
  if (spikeSkipf < 0 || spikeSkipf >= 1)
    throw FeatureComputationError("spike_skipf should lie between [0 1).");

  vector<double> slope;
  bool semilog = false;
  int retVal =
      __ISI_log_slope(isivalues, slope, true, spikeSkipf, maxnSpike, semilog);
  if (retVal >= 0) {
    setVec(DoubleFeatureData, StringData, "ISI_log_slope_skip", slope);
    return static_cast<int>(slope.size());
  } else {
    return retVal;
  }
}

// time from stimulus start to second threshold crossing
int LibV5::time_to_second_spike(mapStr2intVec& IntFeatureData,
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

// time from stimulus start to last spike
int LibV5::time_to_last_spike(mapStr2intVec& IntFeatureData,
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

double calculateInvISI(const std::vector<double>& all_isi_values_vec,
                       size_t index) {
  if (index < all_isi_values_vec.size()) {
    return 1000.0 / all_isi_values_vec[index];
  }
  throw FeatureComputationError("inverse ISI index out of range");
}

int LibV5::inv_ISI_generic(mapStr2intVec& IntFeatureData,
                           mapStr2doubleVec& DoubleFeatureData,
                           mapStr2Str& StringData, size_t index) {
  const auto& all_isi_values_vec =
      getFeature(DoubleFeatureData, {"all_ISI_values"});
  double inv_ISI = calculateInvISI(all_isi_values_vec, index);
  std::string featureName;

  switch (index) {
    case 0:
      featureName = "inv_first_ISI";
      break;
    case 1:
      featureName = "inv_second_ISI";
      break;
    case 2:
      featureName = "inv_third_ISI";
      break;
    case 3:
      featureName = "inv_fourth_ISI";
      break;
    case 4:
      featureName = "inv_fifth_ISI";
      break;
    default:
      if (index == all_isi_values_vec.size() - 1) {
        featureName = "inv_last_ISI";
      } else {
        featureName = "inv_" + std::to_string(index + 1) + "th_ISI";
      }
      break;
  }

  vector<double> inv_ISI_vec = {inv_ISI};
  setVec(DoubleFeatureData, StringData, featureName, inv_ISI_vec);
  return 1;
}

// 1.0 over last ISI (in Hz); returns 0 when no ISI
int LibV5::inv_last_ISI(mapStr2intVec& IntFeatureData,
                        mapStr2doubleVec& DoubleFeatureData,
                        mapStr2Str& StringData) {
  const auto& all_isi_values_vec =
      getFeature(DoubleFeatureData, "all_ISI_values");

  double inv_last_ISI = calculateInvISI(
      all_isi_values_vec, all_isi_values_vec.size() - 1);  // Last ISI
  vector<double> inv_last_ISI_vec = {inv_last_ISI};
  setVec(DoubleFeatureData, StringData, "inv_last_ISI", inv_last_ISI_vec);
  return 1;
}

// 1.0 over time to first spike (in Hz); returns 0 when no spike
int LibV5::inv_time_to_first_spike(mapStr2intVec& IntFeatureData,
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
int LibV5::min_AHP_indices(mapStr2intVec& IntFeatureData,
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

int LibV5::min_AHP_values(mapStr2intVec& IntFeatureData,
                          mapStr2doubleVec& DoubleFeatureData,
                          mapStr2Str& StringData) {
  return -1;
}

// Difference with LibV1 is that this function doesn't return -1 if there are no
// min_AHP_values
int LibV5::AHP_depth_abs(mapStr2intVec& IntFeatureData,
                         mapStr2doubleVec& DoubleFeatureData,
                         mapStr2Str& StringData) {
  const auto& vAHP = getFeature(DoubleFeatureData, "min_AHP_values");
  setVec(DoubleFeatureData, StringData, "AHP_depth_abs", vAHP);
  return vAHP.size();
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

int LibV5::spike_width1(mapStr2intVec& IntFeatureData,
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

int LibV5::AP_begin_indices(mapStr2intVec& IntFeatureData,
                            mapStr2doubleVec& DoubleFeatureData,
                            mapStr2Str& StringData) {
  // Retrieve all required double and int features at once
  const auto& doubleFeatures =
      getFeatures(DoubleFeatureData, {"T", "V", "stim_start", "stim_end"});
  const auto& intFeatures =
      getFeatures(IntFeatureData, {"peak_indices", "min_AHP_indices"});

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
      intFeatures.at("peak_indices"), intFeatures.at("min_AHP_indices"), apbi,
      dTh[0], derivative_window[0]);

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
  int max_slope;
  getCentralDifferenceDerivative(1., v, dv);
  getCentralDifferenceDerivative(1., t, dt);
  transform(dv.begin(), dv.end(), dt.begin(), dvdt.begin(),
            std::divides<double>());

  apei.resize(pi.size());
  vector<int> picopy(pi.begin(), pi.end());
  picopy.push_back(v.size() - 1);

  for (size_t i = 0; i < apei.size(); i++) {
    max_slope =
        distance(dvdt.begin(), std::min_element(dvdt.begin() + picopy[i] + 1,
                                                dvdt.begin() + picopy[i + 1]));
    // assure that the width of the slope is bigger than 4
    apei[i] = distance(dvdt.begin(), find_if(dvdt.begin() + max_slope,
                                             dvdt.begin() + picopy[i + 1],
                                             [derivativethreshold](double x) {
                                               return x >= derivativethreshold;
                                             }));
  }
  return apei.size();
}

int LibV5::AP_end_indices(mapStr2intVec& IntFeatureData,
                          mapStr2doubleVec& DoubleFeatureData,
                          mapStr2Str& StringData) {
  const auto& T = getFeature(DoubleFeatureData, "T");
  const auto& V = getFeature(DoubleFeatureData, "V");
  const auto& peak_indices = getFeature(IntFeatureData, "peak_indices");

  vector<double> dTh;
  int retVal = getParam(DoubleFeatureData, "DownDerivativeThreshold", dTh);
  double downDerivativeThreshold = (retVal <= 0) ? -12.0 : dTh[0];

  vector<int> AP_end_indices;
  retVal = __AP_end_indices(T, V, peak_indices, AP_end_indices,
                            downDerivativeThreshold);
  if (retVal >= 0) {
    setVec(IntFeatureData, StringData, "AP_end_indices", AP_end_indices);
  }
  return retVal;
}

static int __irregularity_index(const vector<double>& isiValues,
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
  const vector<double>& isiValues = getFeature(DoubleFeatureData, "ISI_values");
  vector<double> irregularity_index;
  int retVal = __irregularity_index(isiValues, irregularity_index);
  if (retVal >= 0) {
    setVec(DoubleFeatureData, StringData, "irregularity_index",
           irregularity_index);
  }
  return retVal;
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
int LibV5::number_initial_spikes(mapStr2intVec& IntFeatureData,
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

// Amplitude of the first spike
int LibV5::AP1_amp(mapStr2intVec& IntFeatureData,
                   mapStr2doubleVec& DoubleFeatureData,
                   mapStr2Str& StringData) {
  const vector<double>& AP_amplitudes =
      getFeature(DoubleFeatureData, "AP_amplitude");
  vector<double> AP1_amp;
  AP1_amp.push_back(AP_amplitudes[0]);
  setVec(DoubleFeatureData, StringData, "AP1_amp", AP1_amp);
  return 1;
}

// Amplitude of the first spike
int LibV5::APlast_amp(mapStr2intVec& IntFeatureData,
                      mapStr2doubleVec& DoubleFeatureData,
                      mapStr2Str& StringData) {
  const vector<double>& AP_amplitudes =
      getFeature(DoubleFeatureData, "AP_amplitude");
  vector<double> APlast_amp;
  APlast_amp.push_back(AP_amplitudes[AP_amplitudes.size() - 1]);
  setVec(DoubleFeatureData, StringData, "APlast_amp", APlast_amp);
  return 1;
}

// Peak voltage of the first spike
int LibV5::AP1_peak(mapStr2intVec& IntFeatureData,
                    mapStr2doubleVec& DoubleFeatureData,
                    mapStr2Str& StringData) {
  const vector<double>& peak_voltage =
      getFeature(DoubleFeatureData, "peak_voltage");
  vector<double> AP1_peak;
  AP1_peak.push_back(peak_voltage[0]);
  setVec(DoubleFeatureData, StringData, "AP1_peak", AP1_peak);
  return 1;
}

// Amplitude of the second spike
int LibV5::AP2_amp(mapStr2intVec& IntFeatureData,
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

// Peak voltage of the second spike
int LibV5::AP2_peak(mapStr2intVec& IntFeatureData,
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
int LibV5::AP2_AP1_diff(mapStr2intVec& IntFeatureData,
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
int LibV5::AP2_AP1_peak_diff(mapStr2intVec& IntFeatureData,
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

// Width of the first spike
int LibV5::AP1_width(mapStr2intVec& IntFeatureData,
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
int LibV5::AP2_width(mapStr2intVec& IntFeatureData,
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
int LibV5::APlast_width(mapStr2intVec& IntFeatureData,
                        mapStr2doubleVec& DoubleFeatureData,
                        mapStr2Str& StringData) {
  const vector<double>& spike_half_width =
      getFeature(DoubleFeatureData, "spike_half_width");
  vector<double> APlast_width;
  APlast_width.push_back(spike_half_width[spike_half_width.size() - 1]);
  setVec(DoubleFeatureData, StringData, "APlast_width", APlast_width);
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

int LibV5::AHP_time_from_peak(mapStr2intVec& IntFeatureData,
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
int LibV5::AHP1_depth_from_peak(mapStr2intVec& IntFeatureData,
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
int LibV5::AHP2_depth_from_peak(mapStr2intVec& IntFeatureData,
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

int LibV5::AP_begin_width(mapStr2intVec& IntFeatureData,
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

int LibV5::AP1_begin_voltage(mapStr2intVec& IntFeatureData,
                             mapStr2doubleVec& DoubleFeatureData,
                             mapStr2Str& StringData) {
  const vector<double>& AP_begin_voltage =
      getFeature(DoubleFeatureData, "AP_begin_voltage");
  vector<double> AP1_begin_voltage;
  AP1_begin_voltage.push_back(AP_begin_voltage[0]);
  setVec(DoubleFeatureData, StringData, "AP1_begin_voltage", AP1_begin_voltage);
  return 1;
}

int LibV5::AP2_begin_voltage(mapStr2intVec& IntFeatureData,
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

int LibV5::AP1_begin_width(mapStr2intVec& IntFeatureData,
                           mapStr2doubleVec& DoubleFeatureData,
                           mapStr2Str& StringData) {
  const vector<double>& AP_begin_width =
      getFeature(DoubleFeatureData, "AP_begin_width");
  vector<double> AP1_begin_width;
  AP1_begin_width.push_back(AP_begin_width[0]);
  setVec(DoubleFeatureData, StringData, "AP1_begin_width", AP1_begin_width);
  return 1;
}

int LibV5::AP2_begin_width(mapStr2intVec& IntFeatureData,
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
int LibV5::AP2_AP1_begin_width_diff(mapStr2intVec& IntFeatureData,
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
  const vector<double>& v = getFeature(DoubleFeatureData, "V");
  const vector<double>& t = getFeature(DoubleFeatureData, "T");
  const vector<double>& stimStart = getFeature(DoubleFeatureData, "stim_start");
  const vector<double>& stimEnd = getFeature(DoubleFeatureData, "stim_end");
  vector<double> vd;
  int retVal = __voltage_deflection_begin(v, t, stimStart[0], stimEnd[0], vd);
  if (retVal > 0) {
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

// The mean voltage after the stimulus in (stim_end + 25%*end_period, stim_end +
// 75%*end_period)
int LibV5::voltage_after_stim(mapStr2intVec& IntFeatureData,
                              mapStr2doubleVec& DoubleFeatureData,
                              mapStr2Str& StringData) {
  const vector<double>& v = getFeature(DoubleFeatureData, "V");
  const vector<double>& t = getFeature(DoubleFeatureData, "T");
  const vector<double>& stimEnd = getFeature(DoubleFeatureData, "stim_end");
  double startTime = stimEnd[0] + (t.back() - stimEnd[0]) * .25;
  double endTime = stimEnd[0] + (t.back() - stimEnd[0]) * .75;
  int nCount = 0;
  double vSum = 0;

  for (size_t i = 0; i < t.size(); i++) {
    if (t[i] >= startTime) {
      vSum += v[i];
      nCount++;
    }
    if (t[i] > endTime) break;
  }

  if (nCount == 0) return -1;

  vector<double> vRest = {vSum / nCount};
  setVec(DoubleFeatureData, StringData, "voltage_after_stim", vRest);

  return 1;
}

int LibV5::mean_AP_amplitude(mapStr2intVec& IntFeatureData,
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
  const auto& doubleFeatures =
      getFeatures(DoubleFeatureData,
                  {"V", "T", "stim_start", "stim_end", "AP_phaseslope_range"});
  const auto& intFeatures = getFeatures(IntFeatureData, {"AP_begin_indices"});
  vector<double> ap_phaseslopes;
  int retVal = __AP_phaseslope(doubleFeatures.at("V"), doubleFeatures.at("T"),
                               doubleFeatures.at("stim_start")[0],
                               doubleFeatures.at("stim_end")[0], ap_phaseslopes,
                               intFeatures.at("AP_begin_indices"),
                               doubleFeatures.at("AP_phaseslope_range")[0]);

  if (retVal > 0) {
    setVec(DoubleFeatureData, StringData, "AP_phaseslope", ap_phaseslopes);
  }
  return retVal;
}

int LibV5::all_ISI_values(mapStr2intVec& IntFeatureData,
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

// spike amplitude: peak_voltage - voltage_base
int LibV5::AP_amplitude_from_voltagebase(mapStr2intVec& IntFeatureData,
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

// min_voltage_between_spikes: minimal voltage between consecutive spikes
int LibV5::min_voltage_between_spikes(mapStr2intVec& IntFeatureData,
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

// return (possibly interpolate) voltage trace
int LibV5::voltage(mapStr2intVec& IntFeatureData,
                   mapStr2doubleVec& DoubleFeatureData,
                   mapStr2Str& StringData) {
  const vector<double>& v = getFeature(DoubleFeatureData, "V");
  setVec(DoubleFeatureData, StringData, "voltage", v);
  return v.size();
}

// return (possibly interpolate) current trace
int LibV5::current(mapStr2intVec& IntFeatureData,
                   mapStr2doubleVec& DoubleFeatureData,
                   mapStr2Str& StringData) {
  const vector<double>& i = getFeature(DoubleFeatureData, "I");
  setVec(DoubleFeatureData, StringData, "current", i);
  return i.size();
}

// return (possibly interpolate) time trace
int LibV5::time(mapStr2intVec& IntFeatureData,
                mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData) {
  const vector<double>& t = getFeature(DoubleFeatureData, "T");
  setVec(DoubleFeatureData, StringData, "time", t);
  return t.size();
}

// *** The average voltage during the last 90% of the stimulus duration. ***
int LibV5::steady_state_voltage_stimend(mapStr2intVec& IntFeatureData,
                                        mapStr2doubleVec& DoubleFeatureData,
                                        mapStr2Str& StringData) {
  const auto& doubleFeatures =
      getFeatures(DoubleFeatureData, {"V", "T", "stim_end", "stim_start"});

  const vector<double>& voltages = doubleFeatures.at("V");
  const vector<double>& times = doubleFeatures.at("T");
  const double stimStart = doubleFeatures.at("stim_start")[0];
  const double stimEnd = doubleFeatures.at("stim_end")[0];

  double start_time = stimEnd - 0.1 * (stimEnd - stimStart);
  auto start_it = find_if(times.begin(), times.end(),
                          [start_time](double x) { return x >= start_time; });
  auto stop_it = find_if(times.begin(), times.end(),
                         [stimEnd](double x) { return x >= stimEnd; });

  size_t start_index = distance(times.begin(), start_it);
  size_t stop_index = distance(times.begin(), stop_it);

  double mean = accumulate(voltages.begin() + start_index,
                           voltages.begin() + stop_index, 0.0);
  size_t mean_size = stop_index - start_index;

  vector<double> ssv;
  if (mean_size > 0) {
    mean /= mean_size;
    ssv.push_back(mean);
    setVec(DoubleFeatureData, StringData, "steady_state_voltage_stimend", ssv);
    return 1;
  }
  return -1;
}

int LibV5::voltage_base(mapStr2intVec& IntFeatureData,
                        mapStr2doubleVec& DoubleFeatureData,
                        mapStr2Str& StringData) {
  const vector<double>& v = getFeature(DoubleFeatureData, "V");
  const vector<double>& t = getFeature(DoubleFeatureData, "T");
  const vector<double>& stimStart = getFeature(DoubleFeatureData, "stim_start");

  // Retrieve percentage values or use defaults.
  double vb_start_perc = 0.9;  // Default value
  double vb_end_perc = 1.0;    // Default value
  try {
    auto vb_start_perc_vec =
        getFeature(DoubleFeatureData, "voltage_base_start_perc");
    if (vb_start_perc_vec.size() == 1) vb_start_perc = vb_start_perc_vec[0];
  } catch (const EmptyFeatureError&) {
    // If there's an error, use the default value.
  }

  try {
    auto vb_end_perc_vec =
        getFeature(DoubleFeatureData, "voltage_base_end_perc");
    if (vb_end_perc_vec.size() == 1) vb_end_perc = vb_end_perc_vec[0];
  } catch (const EmptyFeatureError&) {
    // If there's an error, use the default value.
  }

  // Calculate start and end times based on stimStart and percentages.
  double startTime = stimStart[0] * vb_start_perc;
  double endTime = stimStart[0] * vb_end_perc;

  // Validate start and end times.
  if (startTime >= endTime)
    throw FeatureComputationError("voltage_base: startTime >= endTime");

  const auto& precisionThreshold =
      getFeature(DoubleFeatureData, "precision_threshold");

  // Find index range for the time vector within the specified start and end
  // times.
  std::pair<size_t, size_t> time_index =
      get_time_index(t, startTime, endTime, precisionThreshold[0]);

  // Extract sub-vector of voltages based on calculated indices.
  vector<double> subVector(v.begin() + time_index.first,
                           v.begin() + time_index.second);

  // Determine computation mode and calculate voltage base.
  std::string computation_mode;

  int retVal = getStrParam(StringData, "voltage_base_mode", computation_mode);
  if (retVal < 0) return -1;
  double vBase;

  // Perform computation based on the mode.
  if (computation_mode == "mean")
    vBase = vec_mean(subVector);
  else if (computation_mode == "median")
    vBase = vec_median(subVector);
  else
    throw FeatureComputationError("Undefined computational mode. Only mean and median are enabled.");

  vector<double> vRest = {vBase};
  setVec(DoubleFeatureData, StringData, "voltage_base", vRest);
  return 1;
}

int LibV5::current_base(mapStr2intVec& IntFeatureData,
                        mapStr2doubleVec& DoubleFeatureData,
                        mapStr2Str& StringData) {
  const auto& doubleFeatures =
      getFeatures(DoubleFeatureData, {"I", "T", "stim_start"});
  double cb_start_perc = 0.9;  // Default value
  double cb_end_perc = 1.0;    // Default value

  try {
    const double cb_start_perc =
        getFeature(DoubleFeatureData, "current_base_start_perc")[0];
  } catch (const std::runtime_error&) {
  }  // Use default value if not found or an error occurs

  try {
    const double cb_end_perc =
        getFeature(DoubleFeatureData, "current_base_end_perc")[0];
  } catch (const std::runtime_error&) {
  }  // Use default value if not found or an error occurs

  double startTime = doubleFeatures.at("stim_start")[0] * cb_start_perc;
  double endTime = doubleFeatures.at("stim_start")[0] * cb_end_perc;

  if (startTime >= endTime)
    throw FeatureComputationError("current_base: startTime >= endTime");

  vector<double> precisionThreshold;
  int retVal =
      getParam(DoubleFeatureData, "precision_threshold", precisionThreshold);
  if (retVal < 0) return -1;

  std::pair<size_t, size_t> time_index = get_time_index(
      doubleFeatures.at("T"), startTime, endTime, precisionThreshold[0]);

  vector<double> subVector(doubleFeatures.at("I").begin() + time_index.first,
                           doubleFeatures.at("I").begin() + time_index.second);

  double iBase;
  std::string computation_mode;
  retVal = getStrParam(StringData, "current_base_mode", computation_mode);
  if (retVal < 0) return -1;
  if (computation_mode == "mean")
    iBase = vec_mean(subVector);
  else if (computation_mode == "median")
    iBase = vec_median(subVector);
  else
    throw FeatureComputationError("Undefined computational mode. Only mean and median are enabled.");

  vector<double> iRest{iBase};
  setVec(DoubleFeatureData, StringData, "current_base", iRest);
  return 1;
}

size_t get_index(const vector<double>& times, double t) {
  return distance(times.begin(), find_if(times.begin(), times.end(),
                                         [t](double x) { return x >= t; }));
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
    throw FeatureComputationError("No data points to calculate decay_time_constant_after_stim");
  } 
  linear_fit_result fit;
  fit = slope_straight_line_fit(decayTimes, decayValues);

  const double tau = -1.0 / fit.slope;
  return std::abs(tau);
}

// *** Decay time constant measured during decay after the stimulus***
int LibV5::decay_time_constant_after_stim(mapStr2intVec& IntFeatureData,
                                          mapStr2doubleVec& DoubleFeatureData,
                                          mapStr2Str& StringData) {
  const auto& doubleFeatures =
      getFeatures(DoubleFeatureData, {"V", "T", "stim_start", "stim_end"});

  double decay_start_after_stim, decay_end_after_stim;

  try {
    const auto& decayStartFeatures =
        getFeatures(DoubleFeatureData, {"decay_start_after_stim"});
    decay_start_after_stim = decayStartFeatures.at("decay_start_after_stim")[0];
  } catch (const std::runtime_error&) {
    decay_start_after_stim = 1.0;  // Default value if not found
  }

  try {
    const auto& decayEndFeatures =
        getFeatures(DoubleFeatureData, {"decay_end_after_stim"});
    decay_end_after_stim = decayEndFeatures.at("decay_end_after_stim")[0];
  } catch (const std::runtime_error&) {
    decay_end_after_stim = 10.0;  // Default value if not found
  }
  // Validate decay times
  if (decay_start_after_stim >= decay_end_after_stim)
    throw FeatureComputationError("Error decay_start_after_stim small larger than decay_end_after_stim");

  // Perform calculation
  const double val = __decay_time_constant_after_stim(
      doubleFeatures.at("T"), doubleFeatures.at("V"), decay_start_after_stim,
      decay_end_after_stim, doubleFeatures.at("stim_start")[0],
      doubleFeatures.at("stim_end")[0]);

  // Store the result
  vector<double> dtcas{val};
  setVec(DoubleFeatureData, StringData, "decay_time_constant_after_stim",
         dtcas);

  return 1;
}

// Calculate the time constants after each step for a stimuli containing several
// steps, as for example SpikeRec protocols
int LibV5::multiple_decay_time_constant_after_stim(
    mapStr2intVec& IntFeatureData, mapStr2doubleVec& DoubleFeatureData,
    mapStr2Str& StringData) {
  const auto& doubleFeatures = getFeatures(
      DoubleFeatureData, {"V", "T", "multi_stim_start", "multi_stim_end"});
  vector<double> stimsEnd, stimsStart;

  stimsEnd = doubleFeatures.at("multi_stim_end");
  stimsStart = doubleFeatures.at("multi_stim_start");

  // Attempt to get decay parameters, using defaults if not found or if not
  // exactly one element
  double decay_start_after_stim = 1.0;
  double decay_end_after_stim = 10.0;
  try {
    const double decay_start_after_stim =
        getFeature(DoubleFeatureData, "decay_start_after_stim")[0];
  } catch (const std::runtime_error&) {
  }  // Use default value
  try {
    const double decay_end_after_stim =
        getFeature(DoubleFeatureData, "decay_end_after_stim")[0];
  } catch (const std::runtime_error&) {
  }  // Use default value
  vector<double> dtcas;
  for (size_t i = 0; i < stimsStart.size(); i++) {
    double ret_dtcas = __decay_time_constant_after_stim(
        doubleFeatures.at("T"), doubleFeatures.at("V"), decay_start_after_stim,
        decay_end_after_stim, stimsStart[i], stimsEnd[i]);
    dtcas.push_back(ret_dtcas);
  }
  setVec(DoubleFeatureData, StringData,
         "multiple_decay_time_constant_after_stim", dtcas);
  return 1;
}

// compute time constant for the decay from the sag to the steady_state_voltage
// noisy data is expected, so no golden section search is used
// because with noisy data, x>0 often gives a worse logarithmic fit
static int __sag_time_constant(const vector<double>& times,
                               const vector<double>& voltage,
                               const double minimum_voltage,
                               const double steady_state_v,
                               const double sag_amplitude,
                               const double stimStart, const double stimEnd,
                               vector<double>& sagtc) {
  // minimal required length of each decay (indices)
  size_t min_length = 10;

  // get start index
  const size_t decayStartIdx = distance(
      voltage.begin(),
      find_if(voltage.begin(), voltage.end(),
              [minimum_voltage](double v) { return v <= minimum_voltage; }));

  // voltage at which 90% of the sag amplitude has decayed
  double steady_state_90 = steady_state_v - sag_amplitude * 0.1;
  // get end index
  const size_t decayEndIdx = distance(
      voltage.begin(),
      find_if(voltage.begin() + decayStartIdx, voltage.end(),
              [steady_state_90](double v) { return v >= steady_state_90; }));

  // voltage reference by which the voltage (i the decay interval)
  // is going to be substracted
  // there should be no '0' in (decay_v - v_reference),
  // so no problem when the log is taken
  double v_reference = voltage[decayEndIdx];

  // decay interval
  vector<double> VInterval(&voltage[decayStartIdx], &voltage[decayEndIdx]);
  vector<double> TInterval(&times[decayStartIdx], &times[decayEndIdx]);

  // compute time constant
  vector<double> decayValues(decayEndIdx - decayStartIdx);
  for (size_t i = 0; i < VInterval.size(); ++i) {
    const double u0 = std::abs(VInterval[i] - v_reference);
    decayValues[i] = log(u0);
  }
  if (decayValues.size() < min_length) {
    throw FeatureComputationError("Not enough data points to compute time constant.");
  }
  linear_fit_result fit;
  fit = slope_straight_line_fit(TInterval, decayValues);

  // append tau
  sagtc.push_back(std::abs(1.0 / fit.slope));

  return 1;
}

// *** Decay time constant measured from minimum voltage to steady-state
// voltage***
int LibV5::sag_time_constant(mapStr2intVec& IntFeatureData,
                             mapStr2doubleVec& DoubleFeatureData,
                             mapStr2Str& StringData) {
  const auto& doubleFeatures = getFeatures(
      DoubleFeatureData, {"V", "T", "stim_end", "stim_start", "minimum_voltage",
                          "steady_state_voltage_stimend", "sag_amplitude"});
  vector<double> sagtc;
  int retVal = __sag_time_constant(
      doubleFeatures.at("T"), doubleFeatures.at("V"),
      doubleFeatures.at("minimum_voltage")[0],
      doubleFeatures.at("steady_state_voltage_stimend")[0],
      doubleFeatures.at("sag_amplitude")[0], doubleFeatures.at("stim_start")[0],
      doubleFeatures.at("stim_end")[0], sagtc);

  if (retVal > 0) {
    setVec(DoubleFeatureData, StringData, "sag_time_constant", sagtc);
  }
  return retVal;
}

/// *** Voltage deflection between voltage_base and steady_state_voltage_stimend
int LibV5::voltage_deflection_vb_ssse(mapStr2intVec& IntFeatureData,
                                      mapStr2doubleVec& DoubleFeatureData,
                                      mapStr2Str& StringData) {
  const auto& doubleFeatures = getFeatures(
      DoubleFeatureData, {"voltage_base", "steady_state_voltage_stimend"});

  vector<double> voltage_deflection_vb_ssse;
  voltage_deflection_vb_ssse.push_back(
      doubleFeatures.at("steady_state_voltage_stimend")[0] -
      doubleFeatures.at("voltage_base")[0]);
  setVec(DoubleFeatureData, StringData, "voltage_deflection_vb_ssse",
         voltage_deflection_vb_ssse);
  return 1;
}

// *** ohmic input resistance based on voltage_deflection_vb_ssse***

int LibV5::ohmic_input_resistance_vb_ssse(mapStr2intVec& IntFeatureData,
                                          mapStr2doubleVec& DoubleFeatureData,
                                          mapStr2Str& StringData) {
  const auto& doubleFeatures = getFeatures(
      DoubleFeatureData, {"voltage_deflection_vb_ssse", "stimulus_current"});
  const double stimulus_current = doubleFeatures.at("stimulus_current")[0];
  if (stimulus_current == 0)
    throw FeatureComputationError("Stimulus current is zero which will result in division by zero.");
  vector<double> ohmic_input_resistance_vb_ssse;
  ohmic_input_resistance_vb_ssse.push_back(
      doubleFeatures.at("voltage_deflection_vb_ssse")[0] / stimulus_current);
  setVec(DoubleFeatureData, StringData, "ohmic_input_resistance_vb_ssse",
         ohmic_input_resistance_vb_ssse);

  return 1;
}

// *** Diff between maximum voltage during stimulus and voltage_base ***
int LibV5::maximum_voltage_from_voltagebase(mapStr2intVec& IntFeatureData,
                                            mapStr2doubleVec& DoubleFeatureData,
                                            mapStr2Str& StringData) {
  const auto& doubleFeatures =
      getFeatures(DoubleFeatureData, {"maximum_voltage", "voltage_base"});

  vector<double> maximum_voltage_from_voltagebase;
  maximum_voltage_from_voltagebase.push_back(
      doubleFeatures.at("maximum_voltage")[0] -
      doubleFeatures.at("voltage_base")[0]);
  setVec(DoubleFeatureData, StringData, "maximum_voltage_from_voltagebase",
         maximum_voltage_from_voltagebase);
  return 1;
}

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

int LibV5::sag_amplitude(mapStr2intVec& IntFeatureData,
                         mapStr2doubleVec& DoubleFeatureData,
                         mapStr2Str& StringData) {
  const auto& doubleFeatures = getFeatures(
      DoubleFeatureData, {"steady_state_voltage_stimend",
                          "voltage_deflection_vb_ssse", "minimum_voltage"});

  vector<double> sag_amplitude;
  if (doubleFeatures.at("voltage_deflection_vb_ssse")[0] <= 0) {
    sag_amplitude.push_back(
        doubleFeatures.at("steady_state_voltage_stimend")[0] -
        doubleFeatures.at("minimum_voltage")[0]);
  } else
      throw FeatureComputationError("sag_amplitude: voltage_deflection is positive");

  if (!sag_amplitude.empty()) {
    setVec(DoubleFeatureData, StringData, "sag_amplitude", sag_amplitude);
  }
  return sag_amplitude.empty() ? -1 : 1;
}

int LibV5::sag_ratio1(mapStr2intVec& IntFeatureData,
                      mapStr2doubleVec& DoubleFeatureData,
                      mapStr2Str& StringData) {
  // Retrieve all required double features at once
  const auto& doubleFeatures = getFeatures(
      DoubleFeatureData, {"sag_amplitude", "voltage_base", "minimum_voltage"});

  vector<double> sag_ratio1;
  if (doubleFeatures.at("minimum_voltage")[0] == doubleFeatures.at("voltage_base")[0])
    throw FeatureComputationError("voltage_base equals minimum_voltage");

  sag_ratio1.push_back(doubleFeatures.at("sag_amplitude")[0] /
                        (doubleFeatures.at("voltage_base")[0] -
                        doubleFeatures.at("minimum_voltage")[0]));

  if (!sag_ratio1.empty()) {
    setVec(DoubleFeatureData, StringData, "sag_ratio1", sag_ratio1);
  }
  return sag_ratio1.empty() ? -1 : 1;
}

int LibV5::sag_ratio2(mapStr2intVec& IntFeatureData,
                      mapStr2doubleVec& DoubleFeatureData,
                      mapStr2Str& StringData) {
  // Retrieve all required double features at once
  const auto& doubleFeatures = getFeatures(
      DoubleFeatureData,
      {"voltage_base", "minimum_voltage", "steady_state_voltage_stimend"});

  vector<double> sag_ratio2;
  if (doubleFeatures.at("minimum_voltage")[0] == doubleFeatures.at("voltage_base")[0])
    throw FeatureComputationError("voltage_base equals minimum_voltage");

  sag_ratio2.push_back(
      (doubleFeatures.at("voltage_base")[0] -
        doubleFeatures.at("steady_state_voltage_stimend")[0]) /
      (doubleFeatures.at("voltage_base")[0] -
        doubleFeatures.at("minimum_voltage")[0]));

  if (!sag_ratio2.empty()) {
    setVec(DoubleFeatureData, StringData, "sag_ratio2", sag_ratio2);
  }
  return sag_ratio2.empty() ? -1 : 1;
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

int LibV5::AP_peak_upstroke(mapStr2intVec& IntFeatureData,
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

int LibV5::AP_peak_downstroke(mapStr2intVec& IntFeatureData,
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
int LibV5::min_between_peaks_indices(mapStr2intVec& IntFeatureData,
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

int LibV5::min_between_peaks_values(mapStr2intVec& IntFeatureData,
                                    mapStr2doubleVec& DoubleFeatureData,
                                    mapStr2Str& StringData) {
  return 1;
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

int LibV5::AP_width_between_threshold(mapStr2intVec& IntFeatureData,
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

int LibV5::burst_begin_indices(mapStr2intVec& IntFeatureData,
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

int LibV5::burst_end_indices(mapStr2intVec& IntFeatureData,
                             mapStr2doubleVec& DoubleFeatureData,
                             mapStr2Str& StringData) {
  return 1;
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

int LibV5::strict_burst_mean_freq(mapStr2intVec& IntFeatureData,
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

int LibV5::strict_interburst_voltage(mapStr2intVec& IntFeatureData,
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
int LibV5::ADP_peak_indices(mapStr2intVec& IntFeatureData,
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
int LibV5::ADP_peak_values(mapStr2intVec& IntFeatureData,
                           mapStr2doubleVec& DoubleFeatureData,
                           mapStr2Str& StringData) {
  return 1;
}

// strict_stiminterval should be True when using this feature
int LibV5::ADP_peak_amplitude(mapStr2intVec& IntFeatureData,
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

int LibV5::interburst_min_indices(mapStr2intVec& IntFeatureData,
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

int LibV5::interburst_min_values(mapStr2intVec& IntFeatureData,
                                 mapStr2doubleVec& DoubleFeatureData,
                                 mapStr2Str& StringData) {
  return 1;
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
    if (burst_end_indices[i] + 1 < peak_indices.size()) {
      postburst_min_index =
          min_element(v.begin() + peak_indices[burst_end_indices[i]],
                      v.begin() + peak_indices[burst_end_indices[i] + 1]) -
          v.begin();
    } else if (peak_indices[burst_end_indices[i]] < stim_end_index) {
      postburst_min_index =
          min_element(v.begin() + peak_indices[burst_end_indices[i]],
                      v.begin() + stim_end_index) -
          v.begin();
    } else {
      postburst_min_index =
          min_element(v.begin() + peak_indices[burst_end_indices[i]],
                      v.begin() + end_index) -
          v.begin();
    }

    postburst_min_indices.push_back(postburst_min_index);
    postburst_min_values.push_back(v[postburst_min_index]);
  }

  return postburst_min_indices.size();
}

int LibV5::postburst_min_indices(mapStr2intVec& IntFeatureData,
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

int LibV5::postburst_min_values(mapStr2intVec& IntFeatureData,
                                mapStr2doubleVec& DoubleFeatureData,
                                mapStr2Str& StringData) {
  return 1;
}

int LibV5::time_to_interburst_min(mapStr2intVec& IntFeatureData,
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
