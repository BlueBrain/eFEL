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
      
#ifndef CFEATURE_H_
#define CFEATURE_H_

#include "types.h"

#include <map>
#include <string>
#include <vector>
#include <math.h>

#include "FillFptrTable.h"
#include "DependencyTree.h"
#include "eFELLogger.h"

using std::string;
using std::vector;


class cFeature {
  mapStr2intVec mapIntData;
  mapStr2doubleVec mapDoubleData;
  mapStr2Str mapStrData;
  std::map<string, string> featuretypes;
  FILE* fin;
  void fillfeaturetypes();

 public:
  std::map<string, vector<featureStringPair > > fptrlookup;
  eFELLogger logger;

  cFeature(const string& depFile, const string& outdir);
  template <typename T>
  const vector<T> getMapData(const string& strName, const map<string, vector<T>>& mapData);
  const vector<int> getmapIntData(string strName);
  const vector<double> getmapDoubleData(string strName);
  int getmapfptrVec(string strName, vector<feature_function>& vFptr);
  int calc_features(const string& name);
  template <typename T>
  int getFeature(string strName, vector<T>& vec);
  int setFeatureInt(string strName, vector<int>& intVec);
  int setFeatureDouble(string strName, vector<double>& DoubleVec);
  int setFeatureString(const string& key, const string& value);
  int getFeatureString(const string& key, string& value);
  void getTraces(const string& wildcard, vector<string>& traces);
  int printFeature(const char* strFileName);
  int printMapMember(FILE* fp);

  string featuretype(string featurename);
  string getGError();
  void get_feature_names(vector<string>& feature_names);
  int setVersion(string strDepFile);
  double getDistance(string strName, double mean, double std, 
          bool trace_check=true, double error_dist=250);

};

#endif
