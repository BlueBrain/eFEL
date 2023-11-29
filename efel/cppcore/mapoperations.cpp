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

#include "mapoperations.h"

#include "EfelExceptions.h"

extern string GErrorStr;

template <typename T>
std::map<std::string, std::vector<T>> getFeatures(
    const std::map<std::string, std::vector<T>>& allFeatures,
    const std::vector<std::string>& requestedFeatures) {
  std::map<std::string, std::vector<T>> selectedFeatures;
  for (const auto& featureKey : requestedFeatures) {
    auto it = allFeatures.find(featureKey);
    if (it == allFeatures.end()) {
      throw FeatureComputationError("Feature " + featureKey + " not found");
    } else if (it->second.empty()) {
      throw EmptyFeatureError("Feature " + featureKey + " is empty");
    } else {
      selectedFeatures.insert(*it);
    }
  }
  return selectedFeatures;
}

template <typename T>
std::vector<T> getFeature(
    const std::map<std::string, std::vector<T>>& allFeatures,
    const std::string& requestedFeature) {
  // Use getFeatures to retrieve a map with the single requested feature
  auto selectedFeatures = getFeatures(allFeatures, {requestedFeature});

  // Since we requested only one feature, we can directly access the value
  // The exception handling in getFeatures ensures we have a valid, non-empty
  // vector
  return selectedFeatures.at(requestedFeature);
}

/*
 * get(Int|Double|Str)Param provides access to the Int, Double, Str map
 *
 */
template <typename T>
int getParam(std::map<std::string, std::vector<T>>& featureData,
             const std::string& param, std::vector<T>& vec) {
  auto itr = featureData.find(param);
  if (itr == featureData.end()) {
    GErrorStr += "Parameter [" + param +
                 "] is missing in the map. "
                 "In the python interface, this can be set using the "
                 "appropriate setting function\n";
    return -1;
  }
  vec = itr->second;
  return static_cast<int>(vec.size());
}

int getStrParam(mapStr2Str& StringData, const string& param, string& value) {
  mapStr2Str::const_iterator map_it(StringData.find(param));
  if (map_it == StringData.end()) {
    GErrorStr += "Parameter [" + param + "] is missing in string map\n";
    return -1;
  }
  value = map_it->second;
  return 1;
}

template <class T>
void setVec(std::map<std::string, std::vector<T>>& FeatureData,
            mapStr2Str& StringData, string key, const vector<T>& value) {
  string params;
  getStrParam(StringData, "params", params);
  key += params;
  FeatureData[key] = value;
}

/*
 * get(Int|Double)Vec provide access to the Int, Double map
 * the difference to get(Int|Double)Param is that the "params" entry is applied
 * automatically.
 * For non-elementary features the "params" entry is not empty.
 * For handling multiple traces the trace parameters are passed on to elementary
 * features in that way.
 */

template std::map<std::string, std::vector<double>> getFeatures(
    const std::map<std::string, std::vector<double>>& allFeatures,
    const std::vector<std::string>& requestedFeatures);
template std::map<std::string, std::vector<int>> getFeatures(
    const std::map<std::string, std::vector<int>>& allFeatures,
    const std::vector<std::string>& requestedFeatures);
template int getParam(std::map<std::string, std::vector<double>>& featureData,
                      const std::string& param, std::vector<double>& vec);
template int getParam(std::map<std::string, std::vector<int>>& featureData,
                      const std::string& param, std::vector<int>& vec);
template void setVec(std::map<std::string, std::vector<double>>& FeatureData,
                     mapStr2Str& StringData, string key,
                     const vector<double>& value);
template void setVec(std::map<std::string, std::vector<int>>& FeatureData,
                     mapStr2Str& StringData, string key,
                     const vector<int>& value);
template std::vector<double> getFeature(
    const std::map<std::string, std::vector<double>>& allFeatures,
    const std::string& requestedFeature);
template std::vector<int> getFeature(
    const std::map<std::string, std::vector<int>>& allFeatures,
    const std::string& requestedFeature);
