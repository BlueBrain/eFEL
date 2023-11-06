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
#include <math.h>

extern string GErrorStr;

/*
 * get(Int|Double|Str)Param provides access to the Int, Double, Str map
 *
 */
template<typename T>
int getParam(std::map<std::string, std::vector<T>>& featureData,
               const std::string& param, std::vector<T>& vec) {
    auto itr = featureData.find(param);
    if (itr == featureData.end()) {
        GErrorStr += "Parameter [" + param + "] is missing in the map. "
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
void setVec(std::map<std::string, std::vector<T> >& FeatureData, mapStr2Str& StringData,
               string key, const vector<T>& value){
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

template <class T>
int getVec(std::map<std::string, std::vector<T> >& FeatureData, mapStr2Str& StringData,
                 string strFeature, vector<T>& v){
  string params;
  getStrParam(StringData, "params", params);
  strFeature += params;

  typename std::map<std::string, std::vector<T> >::iterator
    mapstr2VecItr(FeatureData.find(strFeature));

  if (mapstr2VecItr == FeatureData.end()) {
    GErrorStr += "\nFeature [" + strFeature + "] is missing\n";
    return -1;
  }
  v = mapstr2VecItr->second;

  return (v.size());
}

template int getParam(std::map<std::string, std::vector<double>>& featureData,
               const std::string& param, std::vector<double>& vec);
template int getParam(std::map<std::string, std::vector<int>>& featureData,
                const std::string& param, std::vector<int>& vec);
template void setVec(std::map<std::string, std::vector<double> >& FeatureData, mapStr2Str& StringData,
               string key, const vector<double>& value);
template void setVec(std::map<std::string, std::vector<int> >& FeatureData, mapStr2Str& StringData,
               string key, const vector<int>& value);
template int getVec(std::map<std::string, std::vector<double> >& FeatureData, mapStr2Str& StringData,
                 string strFeature, vector<double>& v);
template int getVec(std::map<std::string, std::vector<int> >& FeatureData, mapStr2Str& StringData,
                 string strFeature, vector<int>& v);
