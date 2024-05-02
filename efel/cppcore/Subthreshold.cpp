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

#include "Subthreshold.h"

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


// *** The average voltage during the last 90% of the stimulus duration. ***
int Subthreshold::steady_state_voltage_stimend(mapStr2intVec& IntFeatureData,
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

// steady state of the voltage response during hyperpolarizing stimulus,
// elementary feature for E29
// *** steady_state_hyper
static int __steady_state_hyper(const vector<double>& v,
                                const vector<double>& t, double stimend,
                                vector<double>& steady_state_hyper) {
  // Find the iterator pointing to the first time value greater than or equal to
  // stimend
  auto it_stimend = find_if(
      t.begin(), t.end(), [stimend](double t_val) { return t_val >= stimend; });

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

int Subthreshold::steady_state_hyper(mapStr2intVec& IntFeatureData,
                              mapStr2doubleVec& DoubleFeatureData,
                              mapStr2Str& StringData) {
  // Retrieve all required features at once
  const auto& features = getFeatures(DoubleFeatureData, {"V", "T", "stim_end"});

  vector<double> steady_state_hyper;
  int retval =
      __steady_state_hyper(features.at("V"), features.at("T"),
                           features.at("stim_end").front(), steady_state_hyper);

  if (retval > 0) {
    setVec(DoubleFeatureData, StringData, "steady_state_hyper",
           steady_state_hyper);
  }
  return retval;
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

int Subthreshold::steady_state_voltage(mapStr2intVec& IntFeatureData,
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

int Subthreshold::voltage_base(mapStr2intVec& IntFeatureData,
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

int Subthreshold::current_base(mapStr2intVec& IntFeatureData,
                        mapStr2doubleVec& DoubleFeatureData,
                        mapStr2Str& StringData) {
  const auto& doubleFeatures =
      getFeatures(DoubleFeatureData, {"I", "T", "stim_start"});
  double cb_start_perc = 0.9;  // Default value
  double cb_end_perc = 1.0;    // Default value

  try {
     cb_start_perc = getFeature(DoubleFeatureData, "current_base_start_perc")[0];
  } catch (const std::runtime_error&) {
  }  // Use default value if not found or an error occurs

  try {
    cb_end_perc = getFeature(DoubleFeatureData, "current_base_end_perc")[0];
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

int Subthreshold::time_constant(mapStr2intVec& IntFeatureData,
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
int Subthreshold::decay_time_constant_after_stim(mapStr2intVec& IntFeatureData,
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
int Subthreshold::multiple_decay_time_constant_after_stim(
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
    decay_start_after_stim = getFeature(DoubleFeatureData, "decay_start_after_stim")[0];
  } catch (const std::runtime_error&) {
  }  // Use default value
  try {
    decay_end_after_stim = getFeature(DoubleFeatureData, "decay_end_after_stim")[0];
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
int Subthreshold::sag_time_constant(mapStr2intVec& IntFeatureData,
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

int Subthreshold::sag_amplitude(mapStr2intVec& IntFeatureData,
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

int Subthreshold::sag_ratio1(mapStr2intVec& IntFeatureData,
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

int Subthreshold::sag_ratio2(mapStr2intVec& IntFeatureData,
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

// *** ohmic input resistance ***

static int __ohmic_input_resistance(double voltage_deflection,
                                    double stimulus_current,
                                    vector<double>& oir) {
  if (stimulus_current == 0)
    throw FeatureComputationError("Stimulus current is zero which will result in division by zero.");
  oir.push_back(voltage_deflection / stimulus_current);
  return 1;
}

int Subthreshold::ohmic_input_resistance(mapStr2intVec& IntFeatureData,
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

// *** ohmic input resistance based on voltage_deflection_vb_ssse***

int Subthreshold::ohmic_input_resistance_vb_ssse(mapStr2intVec& IntFeatureData,
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

/// *** Voltage deflection between voltage_base and steady_state_voltage_stimend
int Subthreshold::voltage_deflection_vb_ssse(mapStr2intVec& IntFeatureData,
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

int Subthreshold::voltage_deflection(mapStr2intVec& IntFeatureData,
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

int Subthreshold::voltage_deflection_begin(mapStr2intVec& IntFeatureData,
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

// The mean voltage after the stimulus in (stim_end + 25%*end_period, stim_end +
// 75%*end_period)
int Subthreshold::voltage_after_stim(mapStr2intVec& IntFeatureData,
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

int Subthreshold::maximum_voltage(mapStr2intVec& IntFeatureData,
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
int Subthreshold::minimum_voltage(mapStr2intVec& IntFeatureData,
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

// *** Diff between maximum voltage during stimulus and voltage_base ***
int Subthreshold::maximum_voltage_from_voltagebase(mapStr2intVec& IntFeatureData,
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
