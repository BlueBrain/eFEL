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
int getIntParam(mapStr2intVec& IntFeatureData, const string& param,
                vector<int>& vec) {
  mapStr2intVec::iterator mapstr2IntItr(IntFeatureData.find(param));
  if (mapstr2IntItr == IntFeatureData.end()) {
    GErrorStr += "Parameter [" + param + "] is missing in int map\n";
    return -1;
  }
  vec = mapstr2IntItr->second;
  return (vec.size());
}

int getDoubleParam(mapStr2doubleVec& DoubleFeatureData, const string& param,
                   vector<double>& vec) {
  mapStr2doubleVec::iterator mapstr2DoubleItr;
  mapstr2DoubleItr = DoubleFeatureData.find(param);
  if (mapstr2DoubleItr == DoubleFeatureData.end()) {
    GErrorStr += "Parameter [" + param + "] is missing in double map\n";
    return -1;
  }
  vec = mapstr2DoubleItr->second;
  return (vec.size());
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

void setIntVec(mapStr2intVec& IntFeatureData, mapStr2Str& StringData,
               string key, const vector<int>& value) {
  string params;
  getStrParam(StringData, "params", params);
  key += params;
  IntFeatureData[key] = value;
}

void setDoubleVec(mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData,
                  string key, const vector<double>& value) {
  string params;
  getStrParam(StringData, "params", params);
  key += params;
  DoubleFeatureData[key] = value;
}

/*
 * get(Int|Double)Vec provide access to the Int, Double map
 * the difference to get(Int|Double)Param is that the "params" entry is applied
 * automatically.
 * For non-elementary features the "params" entry is not empty.
 * For handling multiple traces the trace parameters are passed on to elementary
 * features in that way.
 */
int getIntVec(mapStr2intVec& IntFeatureData, mapStr2Str& StringData,
              string strFeature, vector<int>& v) {
  string params;
  getStrParam(StringData, "params", params);
  strFeature += params;
  mapStr2intVec::iterator mapstr2IntItr(IntFeatureData.find(strFeature));
  if (mapstr2IntItr == IntFeatureData.end()) {
    GErrorStr += "\nFeature [" + strFeature + "] is missing\n";
    return -1;
  }
  v = mapstr2IntItr->second;

  return (v.size());
}

int getDoubleVec(mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData,
                 string strFeature, vector<double>& v) {
  string params;
  getStrParam(StringData, "params", params);
  strFeature += params;
  mapStr2doubleVec::iterator mapstr2DoubleItr(
      DoubleFeatureData.find(strFeature));
  if (mapstr2DoubleItr == DoubleFeatureData.end()) {
    GErrorStr += "\nFeature [" + strFeature + "] is missing\n";
    return -1;
  }
  v = mapstr2DoubleItr->second;

  return (v.size());
}

int CheckInIntmap(mapStr2intVec& IntFeatureData, mapStr2Str& StringData,
                  string strFeature, int& nSize) {
  string params;
  getStrParam(StringData, "params", params);
  strFeature += params;
  mapStr2intVec::const_iterator mapstr2IntItr(IntFeatureData.find(strFeature));
  if (mapstr2IntItr != IntFeatureData.end()) {
    nSize = mapstr2IntItr->second.size();
    return 1;
  }
  nSize = -1;
  return 0;
}

int CheckInDoublemap(mapStr2doubleVec& DoubleFeatureData,
                     mapStr2Str& StringData, string strFeature, int& nSize) {
  string params;
  getStrParam(StringData, "params", params);
  strFeature += params;
  mapStr2doubleVec::const_iterator mapstr2DoubleItr(
      DoubleFeatureData.find(strFeature));
  if (mapstr2DoubleItr != DoubleFeatureData.end()) {
    nSize = mapstr2DoubleItr->second.size();
    return 1;
  }
  nSize = -1;
  return 0;
}

/*
 *  Take a wildcard string as an argument:
 *  wildcards seperated by ';' e.g. "APWaveForm;soma"
 *  Find traces containing every required wildcard in the name
 *  e.g. "V;APWaveForm200;soma", "V;APWaveForm240;soma"
 *  Finally return a vector of all parameter strings
 *  e.g. (";APWaveForm200;soma", ";APWaveForm240;soma", ...)
 */
void getTraces(mapStr2doubleVec& mapDoubleData, const string& wildcards,
               vector<string>& params) {
  mapStr2doubleVec::const_iterator map_it;
  string featurename;
  params.clear();
  for (map_it = mapDoubleData.begin(); map_it != mapDoubleData.end();
       ++map_it) {
    featurename = map_it->first;
    // find traces
    if (featurename.find("V;") != string::npos) {
      bool match = true;
      int nextpos;
      int oldpos = 1;
      do {
        string param;
        nextpos = wildcards.find(";", oldpos + 1);
        if (nextpos == -1) {
          nextpos = wildcards.size();
        }
        param = wildcards.substr(oldpos, nextpos - oldpos - 1);
        if (featurename.find(param) == string::npos) {
          match = false;
          break;
        }
        oldpos = nextpos;
      } while (nextpos != (int)wildcards.size());
      if (match) {
        params.push_back(featurename.substr(1));
      }
    }
  }
}

// mean over all traces obtained with the same stimulus
int mean_traces_double(mapStr2doubleVec& DoubleFeatureData,
                       const string& feature, const string& stimulus_name,
                       int i_elem, vector<double>& mean) {
  double sum = 0.;
  vector<string> stim_params;
  getTraces(DoubleFeatureData, stimulus_name, stim_params);
  if (stim_params.size() > 0) {
    for (unsigned i = 0; i < stim_params.size(); i++) {
      vector<double> elem_feature;
      getDoubleParam(DoubleFeatureData, feature + stim_params[i], elem_feature);
      if (i_elem > (int)elem_feature.size() - 1 || elem_feature.size() == 0) {
        GErrorStr +=
            "mean_traces_double: feature vector of the elementary feature does "
            "not contain that many elements.\n";
      }
      if (i_elem == -1) {
        sum += elem_feature.back();
      } else {
        sum += elem_feature[i_elem];
      }
    }
    mean.push_back(sum / stim_params.size());
    return stim_params.size();
  } else {
    return -1;
  }
}

// standard deviation over all traces obtained with the same stimulus
int std_traces_double(mapStr2doubleVec& DoubleFeatureData,
                      const string& feature, const string& stimulus_name,
                      double mean, int i_elem, vector<double>& std) {
  double sum = 0.;
  double v;
  vector<string> stim_params;
  getTraces(DoubleFeatureData, stimulus_name, stim_params);
  if (stim_params.size() > 0) {
    for (unsigned i = 0; i < stim_params.size(); i++) {
      vector<double> elem_feature;
      getDoubleParam(DoubleFeatureData, feature + stim_params[i], elem_feature);
      if (i_elem > (int)elem_feature.size() - 1 ||
          (int)elem_feature.size() == 0) {
        GErrorStr +=
            "std_traces_double: feature vector of the elementary feature does "
            "not contain that many elements.\n";
      }
      if (i_elem == -1) {
        v = elem_feature.back();
      } else {
        v = elem_feature[i_elem];
      }
      double deviation = v - mean;
      sum += deviation * deviation;
    }
    std.push_back(sqrt(sum / (double)(stim_params.size() - 1)));
    return stim_params.size();
  } else {
    return -1;
  }
}
