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

#include <fstream>

#include <map>
#include <string>
#include <vector>
#include <math.h>
//#include "Global.h"
#include "FillFptrTable.h"
#include "DependencyTree.h"

using namespace std;

class cFeature {
  map<string, vector<int> > mapIntData;
  map<string, vector<double> > mapDoubleData;
  map<string, string> mapStrData;
  map<string, string> featuretypes;
  FILE* fin;
  void fillfeaturetypes();

 public:
  map<string, vector<pair<fptr, string> > > fptrlookup;
  vector<int>& getmapIntData(string strName);
  vector<double>& getmapDoubleData(string strName);
  fstream logfile;
  bool logging;

  cFeature(const string& depFile, const string& outdir);
  ~cFeature();
  int getmapfptrVec(string strName, vector<fptr>& vFptr);
  int calc_features(const string& name);
  int setFeatureInt(string strName, vector<int>& intVec);
  int getFeatureInt(string strName, vector<int>& vec);
  int setFeatureDouble(string strName, vector<double>& DoubleVec);
  int getFeatureDouble(string strName, vector<double>& vec);
  int setFeatureString(const string& key, const string& value);
  int getFeatureString(const string& key, string& value);
  void getTraces(const string& wildcard, vector<string>& traces);
  int printFeature(const char* strFileName);
  int printMapMember(FILE* fp);

  string featuretype(string featurename);
  string getGError();
  void get_feature_names(vector<string>& feature_names);
  int setVersion(string strDepFile);
  double getDistance(string strName, double mean, double std);

  // calculation of GA errors
  template<typename T>
  double calc_error_bio(const vector<T>& v, double bio_mean, double bio_sd)
  {
    if (v.size() != 0) {
      double error = 0.;
      for (size_t i = 0; i < v.size(); i++) {
        error += fabs(v[i] - bio_mean);
      }
      return error / bio_sd / v.size();
    } else {
      return 250.;
    }
  }

  template<typename T>
  void appendtolog(const vector<T>& v, size_t max = 10)
  {
    for (size_t i = 0; i < v.size() && i < max; i++) {
      logfile << " " << v[i];
    }
    if (v.size() > max) {
      logfile << " ...";
    }
  }
};

#endif
