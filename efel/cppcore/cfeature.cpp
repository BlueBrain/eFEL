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
#include "Global.h"

#include <cstdlib>
#include <iostream>

#include <time.h>

using std::cout;
using std::endl;



cFeature::cFeature(const string& strDepFile, const string& outdir)
  : logger(outdir)
{
  FillFptrTable();
  mapFptrLib["LibV1"] = &FptrTableV1;
  mapFptrLib["LibV2"] = &FptrTableV2;
  mapFptrLib["LibV3"] = &FptrTableV3;
  mapFptrLib["LibV4"] = &FptrTableV4;
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

int cFeature::setVersion(string strDepFile) {
  FptrTable.clear();
  /*
  map<string, vector< pair< fptr, string > > >::iterator mapVecItr;
  vector< pair< fptr, string > > *vecFptr;
  for(mapVecItr = pFeature->fptrlookup.begin(); mapVecItr !=
  pFeature->fptrlookup.end(); mapVecItr++){
      vecFptr   = &(mapVecItr->second);
      vecFptr->clear();
  }
  */
  fptrlookup.clear();
  cTree DepTree(strDepFile.c_str());
  DepTree.setFeaturePointers(mapFptrLib, &FptrTable, &fptrlookup);
  return 1;
}

void cFeature::fillfeaturetypes() {
  // initialize feature type information
  featuretypes["peak_indices"] = "int";
  featuretypes["trace_check"] = "int";
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
  featuretypes["E6"] = "double";
  featuretypes["E7"] = "double";
  featuretypes["AP_height"] = "double";
  featuretypes["AP_amplitude"] = "double";
  featuretypes["single_burst_ratio"] = "double";
  featuretypes["BPAPatt2"] = "double";
  featuretypes["BPAPatt3"] = "double";
  featuretypes["E39"] = "double";
  featuretypes["E39_cod"] = "double";
  featuretypes["E2"] = "double";
  featuretypes["E3"] = "double";
  featuretypes["E4"] = "double";
  featuretypes["E5"] = "double";
  featuretypes["E8"] = "double";
  featuretypes["E9"] = "double";
  featuretypes["E10"] = "double";
  featuretypes["E11"] = "double";
  featuretypes["E12"] = "double";
  featuretypes["E13"] = "double";
  featuretypes["E14"] = "double";
  featuretypes["E15"] = "double";
  featuretypes["E16"] = "double";
  featuretypes["E17"] = "double";
  featuretypes["E18"] = "double";
  featuretypes["E19"] = "double";
  featuretypes["E20"] = "double";
  featuretypes["E21"] = "double";
  featuretypes["E22"] = "double";
  featuretypes["E23"] = "double";
  featuretypes["E24"] = "double";
  featuretypes["E25"] = "double";
  featuretypes["E26"] = "double";
  featuretypes["E27"] = "double";
  featuretypes["E40"] = "double";
  featuretypes["steady_state_hyper"] = "double";
  featuretypes["AP_width"] = "double";
  featuretypes["doublet_ISI"] = "double";
  featuretypes["adaptation_index2"] = "double";
  featuretypes["ISI_CV"] = "double";
  featuretypes["AHP_depth_abs_slow"] = "double";
  featuretypes["AHP_slow_time"] = "double";
  featuretypes["Spikecount"] = "int";
  featuretypes["AHP_depth"] = "double";
  featuretypes["amp_drop_first_second"] = "double";
  featuretypes["amp_drop_first_last"] = "double";
  featuretypes["amp_drop_second_last"] = "double";
  featuretypes["max_amp_difference"] = "double";
  featuretypes["depolarized_base"] = "double";
  featuretypes["burst_number"] = "int";

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

  featuretypes["BPAPHeightLoc1"] = "double";
  featuretypes["BPAPAmplitudeLoc1"] = "double";
  featuretypes["BPAPAmplitudeLoc2"] = "double";
  featuretypes["BPAPHeightLoc2"] = "double";
  featuretypes["check_AISInitiation"] = "double";
  featuretypes["AP_phaseslope"] = "double";
  featuretypes["AP_phaseslope_AIS"] = "double";
  featuretypes["BAC_width"] = "double";
  featuretypes["BAC_maximum_voltage"] = "double";
  featuretypes["all_ISI_values"] = "double";
  featuretypes["AP_amplitude_from_voltagebase"] = "double";
  featuretypes["min_voltage_between_spikes"] = "double";
  featuretypes["voltage"] = "double";
  featuretypes["steady_state_voltage_stimend"] = "double";
  featuretypes["voltage_base"] = "double";
  featuretypes["decay_time_constant_after_stim"] = "double";
  featuretypes["maximum_voltage_from_voltagebase"] = "double";

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

vector<int>& cFeature::getmapIntData(string strName) {
  map<string, vector<int> >::iterator mapstr2IntItr;
  mapstr2IntItr = mapIntData.find(strName);
  if (mapstr2IntItr == mapIntData.end()) {
    GErrorStr += "Feature [" + strName + "] data is missing\n";
  }
  return mapstr2IntItr->second;
}
vector<double>& cFeature::getmapDoubleData(string strName) {
  map<string, vector<double> >::iterator mapstr2DoubleItr;
  mapstr2DoubleItr = mapDoubleData.find(strName);
  if (mapstr2DoubleItr == mapDoubleData.end()) {
    GErrorStr += "Feature [" + strName + "] data is missing\n";
  }
  return mapstr2DoubleItr->second;
}

/*
int cFeature::getmapfptrVec(string strName, vector<fptr> &vFptr){
    map<string, vector< fptr > >::iterator mapFptrItr;
    mapFptrItr= FptrLookup.find(strName);
    if(mapFptrItr == FptrLookup.end()) {GErrorStr = GErrorStr + string("Feature
[") + strName + "] dependency is missing\n"; return -1;}
    vFptr = mapFptrItr->second;
    return 1;
}
*/

int cFeature::printMapMember(FILE* fp) {
  map<string, vector<int> >::iterator mapstr2IntItr;
  fprintf(fin, "\n\n\n IntData.....");
  for (mapstr2IntItr = mapIntData.begin(); mapstr2IntItr != mapIntData.end();
       mapstr2IntItr++)
    fprintf(fin, "\n\t%s", mapstr2IntItr->first.c_str());
  fprintf(fin, "\n\n DoubleData..........");
  map<string, vector<double> >::iterator mapstr2DoubleItr;
  for (mapstr2DoubleItr = mapDoubleData.begin();
       mapstr2DoubleItr != mapDoubleData.end(); mapstr2DoubleItr++)
    fprintf(fin, "\n\t%s", mapstr2DoubleItr->first.c_str());
  return 1;
}

int cFeature::setFeatureInt(string strName, vector<int>& v) {
  logger << "Set " << strName << ":" << v << endl;
  // printf ("Setting int feature [%s] = %d\n", strName.c_str(),v[0]);
  mapIntData[strName] = v;
  return 1;
}

/*
 *  Take a wildcard string as an argument:
 *  wildcards seperated by ';' e.g. "APWaveForm;soma"
 *  Find traces containing every required wildcard in the name
 *  e.g. "V;APWaveForm200;soma", "V;APWaveForm240;soma"
 *  Finally return a vector of all parameter strings
 *  e.g. (";APWaveForm200;soma", ";APWaveForm240;soma", ...)
 */
void cFeature::getTraces(const string& wildcards, vector<string>& params) {
  map<string, vector<double> >::const_iterator map_it;
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

int cFeature::calc_features(const string& name) {
  // stimulus extension
  map<string, vector<featureStringPair> >::const_iterator lookup_it(
      fptrlookup.find(name));
  if (lookup_it == fptrlookup.end()) {
    fprintf(stderr,
            "\nFeature [ %s ] dependency file entry or pointer table entry is "
            "missing. Exiting\n",
            name.c_str());
    fflush(stderr);
    exit(1);
  }

  bool last_failed = false;

  for (vector<featureStringPair>::const_iterator pfptrstring =
           lookup_it->second.begin();
       pfptrstring != lookup_it->second.end(); ++pfptrstring) {
    // set parameters, for now only the wildcard 'stimulusname'
    //
    feature_function function = pfptrstring->first;
    string wildcard = pfptrstring->second;
    if (wildcard.empty()) {
      // make sure that
      //  - the feature is called only once
      //  - the feature operates on "V" and "T" if it operates on traces at all
      setFeatureString("params", "");
      if (function(mapIntData, mapDoubleData, mapStrData) < 0) {
        // GErrorStr += "\nFeature [" + name + "] called twice, or doesn't
        // operate on V and T.";
        // return -1;i
        last_failed = true;
      } else {
        last_failed = false;
      }
    } else {
      // make sure that
      //  -the feature is called once for every trace according to the wildcard
      //  -the feature operates on each trace
      vector<string> params;
      // TODO
      // read stimulus configuration file and parse additional parameters
      // such as number of required traces
      getTraces(wildcard, params);
      if (params.empty()) {
        GErrorStr += "\nMissing trace with wildcards " + wildcard;
        return -1;
      }
      for (unsigned i = 0; i < params.size(); i++) {
        // setting the "params" entry here makes sure that the required features
        // require specific traces also
        setFeatureString("params", params[i]);
        if (function(mapIntData, mapDoubleData, mapStrData) < 0) {
          last_failed = true;
        } else {
          last_failed = false;
        }
      }
    }
  }
  if (last_failed) {
    return -1;
  } else {
    // success
    return 0;
  }
}

int cFeature::getFeatureInt(string strName, vector<int>& vec) {
  logger << "Going to calculate feature " << strName << " ..." << endl;
  if (calc_features(strName) < 0) {
    logger << "Failed to calculate feature " << strName << ": " << GErrorStr
              << endl;
    return -1;
  }
  vec = getmapIntData(strName);

  logger << "Calculated feature " << strName << ":" << vec << endl;

  return vec.size();
}

int cFeature::getFeatureDouble(string strName, vector<double>& vec) {
  logger << "Going to calculate feature " << strName << " ..." << endl;
  if (calc_features(strName) < 0) {
    logger << "Failed to calculate feature " << strName << ": " << GErrorStr
              << endl;
    return -1;
  }
  vec = getmapDoubleData(strName);

  logger << "Calculated feature " << strName << ":" << vec << endl;

  return vec.size();
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
/*
 * not needed
int cFeature::ifExist(string strName, map<string, vector< int > > &mapIntData){
    map<string, vector< int > >::iterator mapItrInt;
    int n = mapIntData.size();
    string str;
    for(mapItrInt = mapIntData.begin(); mapItrInt != mapIntData.end();
mapItrInt++){
        str = mapItrInt->first;
        if( str == strName) return 1;
    }
    return 0;
}

int cFeature::ifExist(string strName, map<string, vector< double > >
&mapDoubleData){
    map<string, vector< double > >::iterator mapItrDouble;
    int n = mapDoubleData.size();
    string str;
    for(mapItrDouble = mapDoubleData.begin(); mapItrDouble !=
mapDoubleData.end(); mapItrDouble++){
        str = mapItrDouble->first;
        if( str == strName) return 1;
    }
    return 0;
}
*/

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

int cFeature::printFeature(const char* strFileName) {
  FILE* fp = fopen(strFileName, "w");
  if (fp) {
    map<string, vector<int> >::iterator mapItrInt;
    int n = mapIntData.size();
    fprintf(fp, "\n mapIntData.. Total element = [%d]", n);
    for (mapItrInt = mapIntData.begin(); mapItrInt != mapIntData.end();
         mapItrInt++) {
      string str = mapItrInt->first;
      vector<int>* v = &(mapItrInt->second);
      fprintf(fp, "\n ParameterName = [%s] size = [%d]\n\t", str.c_str(),
              (int)v->size());
      for (unsigned j = 0; j < v->size(); j++) {
        fprintf(fp, "[%d]", v->at(j));
      }
    }

    map<string, vector<double> >::iterator mapItrDouble;
    n = mapDoubleData.size();
    fprintf(fp, "\n mapDoubleData.. Total element = [%d]", n);
    for (mapItrDouble = mapDoubleData.begin();
         mapItrDouble != mapDoubleData.end(); mapItrDouble++) {
      string str = mapItrDouble->first;
      vector<double>* v = &(mapItrDouble->second);
      fprintf(fp, "\n ParameterName = [%s] size = [%d]\n\t", str.c_str(),
              (int)v->size());
      for (unsigned j = 0; j < v->size(); j++) {
        fprintf(fp, "[%f]", v->at(j));
      }
    }
    fclose(fp);
  }
  return 1;
}

double cFeature::getDistance(string strName, double mean, double std, 
        bool trace_check) {

  vector<double> feature_vec;
  vector<int> feature_veci;
  string featureType;
  int retVal, intFlag;
  double dError = 0;

  // Check if a the trace doesn't contain any spikes outside of the stimulus
  // interval
  if (trace_check) {
      retVal = getFeatureInt("trace_check", feature_veci);
      if (retVal < 0) {
          return 250.0;
      }
  }

  // check datatype of feature
  featureType = featuretype(strName);
  if (featureType.empty()) {
    printf("Error : Feature [%s] not found. Exiting..\n", strName.c_str());
    exit(1);
  }

  if (featureType == "int") {
    retVal = getFeatureInt(strName, feature_veci);
    intFlag = 1;
  } else {
    retVal = getFeatureDouble(strName, feature_vec);
    intFlag = 0;
  }
  // printf("\n Calculating distance for [%s] values [", strName.c_str());
  if (retVal <= 0) {
    // printf ("\n Error in feature calculation... [%s]\n",GErrorStr.c_str() );
    return 250;
  } else {
    if (intFlag) {
      for (unsigned i = 0; i < feature_veci.size(); i++) {
        // printf("%d\t", feature_veci[i]);
        dError = dError + fabs(feature_veci[i] - mean);
      }
      dError = dError / std / feature_veci.size();
      if (dError != dError) {
        // printf("Warning: Error distance calculation generated NaN, returning
        // 250\n");
        return 250;
      }
      // printf("] TotalError = %f\n", dError);
      return dError;
    } else {
      for (unsigned i = 0; i < feature_vec.size(); i++) {
        // printf("%f\t", feature_vec[i]);
        dError = dError + fabs(feature_vec[i] - mean);
      }
      dError = dError / std / feature_vec.size();
      if (dError != dError) {
        printf(
            "Warning: Error distance calculation generated NaN, returning "
            "250\n");
        return 250;
      }
      // printf("] TotalError = %f\n", dError);
      return dError;
    }
  }
  return 250.0;
}

string cFeature::featuretype(string featurename) {
  int npos = featurename.find(";");
  if (npos >= 0) {
    featurename = featurename.substr(0, npos);
  }
  string type(featuretypes[featurename]);
  if (type.empty()) {
    GErrorStr += featurename + "missing in featuretypes map.\n";
  }
  return type;
}

string cFeature::getGError() {
  string error(GErrorStr);
  GErrorStr.clear();
  return error;
}
