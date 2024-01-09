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

#include "cfeature.h"

#include <time.h>

#include <cstdlib>
#include <iostream>

#include "Global.h"

using std::endl;

cFeature::cFeature(const string& strDepFile, const string& outdir)
    : logger(outdir) {
  FillFptrTable();
  mapFptrLib["LibV1"] = &FptrTableV1;
  mapFptrLib["LibV2"] = &FptrTableV2;
  mapFptrLib["LibV3"] = &FptrTableV3;
  mapFptrLib["LibV5"] = &FptrTableV5;

  fillfeaturetypes();

  cTree DepTree(strDepFile.c_str());
  // get Error
  if (DepTree.ErrorStr.length() != 0) {
    GErrorStr = DepTree.ErrorStr;
  }
  int retVal = DepTree.setFeaturePointers(mapFptrLib, &FptrTable, &fptrlookup);
  if (retVal < 0) {
    GErrorStr = DepTree.ErrorStr;
  }

  // log output
  time_t rawtime;
  time(&rawtime);
  logger << "\n" << ctime(&rawtime) << "Initializing new session." << endl;
  logger << "Using dependency file: " << strDepFile << endl;
}

template <typename T>
const vector<T> cFeature::getMapData(const string& strName,
                                     const map<string, vector<T>>& mapData) {
  auto mapItr = mapData.find(strName);
  if (mapItr == mapData.end()) {
    return vector<T>{};  // Return an empty vector without modifying the map
  }
  return mapItr->second;
}

const vector<int> cFeature::getmapIntData(string strName) {
  return getMapData(strName, mapIntData);
}

const vector<double> cFeature::getmapDoubleData(string strName) {
  return getMapData(strName, mapDoubleData);
}

void cFeature::fillfeaturetypes() {
  // initialize feature type information
  featuretypes["peak_indices"] = "int";
  featuretypes["ISI_values"] = "double";
  featuretypes["peak_voltage"] = "double";
  featuretypes["burst_ISI_indices"] = "int";
  featuretypes["mean_frequency"] = "double";
  featuretypes["peak_time"] = "double";
  featuretypes["time_to_first_spike"] = "double";
  featuretypes["min_AHP_indices"] = "int";
  featuretypes["min_AHP_values"] = "double";
  featuretypes["adaptation_index"] = "double";
  featuretypes["spike_width2"] = "double";
  featuretypes["spike_half_width"] = "double";
  featuretypes["burst_mean_freq"] = "double";
  featuretypes["interburst_voltage"] = "double";
  featuretypes["voltage_base"] = "double";
  featuretypes["AHP_depth_abs"] = "double";
  featuretypes["time_constant"] = "double";
  featuretypes["voltage_deflection"] = "double";
  featuretypes["voltage_deflection_vb_ssse"] = "double";
  featuretypes["ohmic_input_resistance"] = "double";
  featuretypes["ohmic_input_resistance_vb_ssse"] = "double";
  featuretypes["maximum_voltage"] = "double";
  featuretypes["minimum_voltage"] = "double";
  featuretypes["steady_state_voltage"] = "double";
  featuretypes["AP_begin_indices"] = "int";
  featuretypes["AP_rise_indices"] = "int";
  featuretypes["AP_end_indices"] = "int";
  featuretypes["AP_fall_indices"] = "int";
  featuretypes["AP_duration"] = "double";
  featuretypes["AP_duration_half_width"] = "double";
  featuretypes["AP_rise_time"] = "double";
  featuretypes["AP_fall_time"] = "double";
  featuretypes["AP_rise_rate"] = "double";
  featuretypes["AP_fall_rate"] = "double";
  featuretypes["fast_AHP"] = "double";
  featuretypes["AP_amplitude_change"] = "double";
  featuretypes["AP_amplitude_diff"] = "double";
  featuretypes["AHP_depth_diff"] = "double";
  featuretypes["AP_duration_change"] = "double";
  featuretypes["AP_rise_rate_change"] = "double";
  featuretypes["AP_fall_rate_change"] = "double";
  featuretypes["fast_AHP_change"] = "double";
  featuretypes["AP_duration_half_width_change"] = "double";
  featuretypes["AP_height"] = "double";
  featuretypes["AP_amplitude"] = "double";
  featuretypes["single_burst_ratio"] = "double";
  featuretypes["steady_state_hyper"] = "double";
  featuretypes["AP_width"] = "double";
  featuretypes["doublet_ISI"] = "double";
  featuretypes["adaptation_index2"] = "double";
  featuretypes["ISI_CV"] = "double";
  featuretypes["AHP_depth_abs_slow"] = "double";
  featuretypes["AHP_slow_time"] = "double";
  featuretypes["AHP_depth"] = "double";
  featuretypes["AHP_depth_slow"] = "double";
  featuretypes["amp_drop_first_second"] = "double";
  featuretypes["amp_drop_first_last"] = "double";
  featuretypes["amp_drop_second_last"] = "double";
  featuretypes["max_amp_difference"] = "double";
  featuretypes["depolarized_base"] = "double";

  // LibV5
  featuretypes["ISI_log_slope"] = "double";
  featuretypes["ISI_semilog_slope"] = "double";
  featuretypes["ISI_log_slope_skip"] = "double";
  featuretypes["time_to_second_spike"] = "double";
  featuretypes["time_to_last_spike"] = "double";
  featuretypes["inv_first_ISI"] = "double";
  featuretypes["inv_second_ISI"] = "double";
  featuretypes["inv_third_ISI"] = "double";
  featuretypes["inv_fourth_ISI"] = "double";
  featuretypes["inv_fifth_ISI"] = "double";
  featuretypes["inv_last_ISI"] = "double";
  featuretypes["inv_time_to_first_spike"] = "double";
  featuretypes["irregularity_index"] = "double";
  featuretypes["number_initial_spikes"] = "int";
  featuretypes["AP1_amp"] = "double";
  featuretypes["APlast_amp"] = "double";
  featuretypes["AP1_peak"] = "double";
  featuretypes["AP2_amp"] = "double";
  featuretypes["AP2_peak"] = "double";
  featuretypes["AP2_AP1_diff"] = "double";
  featuretypes["AP2_AP1_peak_diff"] = "double";
  featuretypes["AP1_width"] = "double";
  featuretypes["AP2_width"] = "double";
  featuretypes["APlast_width"] = "double";

  featuretypes["AHP_depth_from_peak"] = "double";
  featuretypes["AHP_time_from_peak"] = "double";
  featuretypes["AHP1_depth_from_peak"] = "double";
  featuretypes["AHP2_depth_from_peak"] = "double";

  featuretypes["AP_begin_width"] = "double";
  featuretypes["AP_begin_time"] = "double";
  featuretypes["AP_begin_voltage"] = "double";

  featuretypes["AP1_begin_width"] = "double";
  featuretypes["AP1_begin_voltage"] = "double";

  featuretypes["AP2_begin_width"] = "double";
  featuretypes["AP2_begin_voltage"] = "double";

  featuretypes["voltage_deflection_begin"] = "double";
  featuretypes["mean_AP_amplitude"] = "double";

  featuretypes["is_not_stuck"] = "int";

  featuretypes["voltage_after_stim"] = "double";
  featuretypes["AP2_AP1_begin_width_diff"] = "double";

  featuretypes["AP_phaseslope"] = "double";
  featuretypes["all_ISI_values"] = "double";
  featuretypes["AP_amplitude_from_voltagebase"] = "double";
  featuretypes["min_voltage_between_spikes"] = "double";
  featuretypes["voltage"] = "double";
  featuretypes["current"] = "double";
  featuretypes["time"] = "double";
  featuretypes["steady_state_voltage_stimend"] = "double";
  featuretypes["voltage_base"] = "double";
  featuretypes["current_base"] = "double";
  featuretypes["decay_time_constant_after_stim"] = "double";
  featuretypes["multiple_decay_time_constant_after_stim"] = "double";
  featuretypes["sag_time_constant"] = "double";
  featuretypes["maximum_voltage_from_voltagebase"] = "double";
  featuretypes["sag_amplitude"] = "double";
  featuretypes["sag_ratio1"] = "double";
  featuretypes["sag_ratio2"] = "double";
  featuretypes["AP_peak_upstroke"] = "double";
  featuretypes["AP_peak_downstroke"] = "double";
  featuretypes["min_between_peaks_indices"] = "int";
  featuretypes["min_between_peaks_values"] = "double";
  featuretypes["AP_width_between_threshold"] = "double";

  featuretypes["burst_begin_indices"] = "int";
  featuretypes["burst_end_indices"] = "int";
  featuretypes["strict_burst_mean_freq"] = "double";
  featuretypes["strict_interburst_voltage"] = "double";

  featuretypes["ADP_peak_indices"] = "int";
  featuretypes["ADP_peak_values"] = "double";
  featuretypes["ADP_peak_amplitude"] = "double";

  featuretypes["interburst_min_indices"] = "int";
  featuretypes["interburst_min_values"] = "double";
  featuretypes["postburst_min_indices"] = "int";
  featuretypes["postburst_min_values"] = "double";
  featuretypes["time_to_interburst_min"] = "double";

  // end of feature types
}

void cFeature::get_feature_names(vector<string>& feature_names) {
  feature_names.clear();
  feature_names.reserve(featuretypes.size());
  for (map<string, string>::iterator it = featuretypes.begin();
       it != featuretypes.end(); ++it) {
    feature_names.push_back(it->first);
  }
}

int cFeature::setFeatureInt(string strName, vector<int>& v) {
  logger << "Set " << strName << ":" << v << endl;
  // printf ("Setting int feature [%s] = %d\n", strName.c_str(),v[0]);
  mapIntData[strName] = v;
  return 1;
}

int cFeature::calc_features(const std::string& name) {
  // stimulus extension
  auto lookup_it = fptrlookup.find(name);
  if (lookup_it == fptrlookup.end()) {
    throw std::runtime_error(
        "Feature dependency file entry or pointer table entry for '" + name +
        "' is missing.");
  }

  bool last_failed = false;

  for (const feature_function& function : lookup_it->second) {
    setFeatureString("params", "");

    if (function(mapIntData, mapDoubleData, mapStrData) < 0) {
      last_failed = true;
    } else {
      last_failed = false;
    }
  }

  return last_failed ? -1 : 0;  // -1 if the last attempt failed, 0 otherwise
}

template <typename T>
int cFeature::getFeature(string strName, vector<T>& vec) {
  const map<string, vector<T>>* dataMap;
  if constexpr (std::is_same_v<T, int>)
    dataMap = &mapIntData;
  else  // cppcore sends only int or double, no other types
    dataMap = &mapDoubleData;

  // 1) Check if the feature is in the map.
  vec = getMapData<T>(strName, *dataMap);
  if (!vec.empty()) {
    logger << "Reusing computed value of " << strName << "." << endl;
    return vec.size();
  } else {
    // 2) If it's not in the map, compute.
    logger << "Going to calculate feature " << strName << " ..." << endl;
    int retVal = 0;
    try {
      retVal = calc_features(strName);
    } catch (const FeatureComputationError& e) {
      GErrorStr += e.what();
      return -1;
    } catch (const EmptyFeatureError& e) {
      GErrorStr += e.what();
      return -1;
    }
    if (retVal < 0) {
      logger << "Failed to calculate feature " << strName << ": " << GErrorStr
             << endl;
      return -1;
    }
    vec = getMapData<T>(strName, *dataMap);
    if (vec.empty()) GErrorStr += "Feature [" + strName + "] data is missing\n";

    logger << "Calculated feature " << strName << ":" << vec << endl;
    return vec.size();
  }
}

int cFeature::setFeatureString(const string& key, const string& value) {
  logger << "Set " << key << ": " << value << endl;
  mapStrData[key] = value;
  return 1;
}

int cFeature::getFeatureString(const string& key, string& value) {
  map<string, string>::const_iterator pstrstr(mapStrData.find(key));
  if (pstrstr != mapStrData.end()) {
    value = pstrstr->second;
    return 1;
  } else {
    GErrorStr += "String parameter [" + key + "] not in map.\n";
    return -1;
  }
}

int cFeature::setFeatureDouble(string strName, vector<double>& v) {
  if (mapDoubleData.find(strName) != mapDoubleData.end()) {
    if (strName == "V") {
      logger << "Feature \"V\" set. New trace, clearing maps." << endl;
      mapDoubleData.clear();
      mapIntData.clear();
      mapStrData.clear();
    }
  }
  mapDoubleData[strName] = v;

  // log data output
  logger << "Set " << strName << ":" << v << endl;

  return 1;
}

string cFeature::featuretype(string featurename) {
  if (featurename == "__test_efel_assertion__")  // for testing only
    throw EfelAssertionError("Test efel assertion is successfully triggered.");
  string type = featuretypes[featurename];
  if (type != "int" && type != "double")
    throw std::runtime_error("Unknown feature name: " + featurename);
  return type;
}

string cFeature::getGError() {
  string error(GErrorStr);
  GErrorStr.clear();
  return error;
}

template int cFeature::getFeature(string strName, vector<int>& vec);
template int cFeature::getFeature(string strName, vector<double>& vec);
template const vector<int> cFeature::getMapData(const string& strName, const map<string, vector<int>>& mapData);
template const vector<double> cFeature::getMapData(const string& strName, const map<string, vector<double>>& mapData);
