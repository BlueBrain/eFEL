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


#ifndef MAPOPERATIONS_H
#define MAPOPERATIONS_H

#include "types.h"

#include <string>
#include <vector>

using std::string;
using std::vector;

extern string GErrorStr;

int getIntParam(mapStr2intVec& IntFeatureData, const string& param,
                vector<int>& vec);
int getDoubleParam(mapStr2doubleVec& DoubleFeatureData, const string& param,
                   vector<double>& vec);
int getStrParam(mapStr2Str& StringData, const string& param, string& value);
template <class T>
void setVec(std::map<std::string, std::vector<T> >& FeatureData, mapStr2Str& StringData,
               string key, const vector<T>& value);
template <class T>
int getVec(std::map<std::string, std::vector<T> >& FeatureData, mapStr2Str& StringData,
                 string strFeature, vector<T>& v);
// eCode feature convenience function
int mean_traces_double(mapStr2doubleVec& DoubleFeatureData,
                       const string& feature, const string& stimulus_name,
                       int i_elem, vector<double>& mean);
int std_traces_double(mapStr2doubleVec& DoubleFeatureData,
                      const string& feature, const string& stimulus_name,
                      double mean, int i_elem, vector<double>& std);
void getTraces(mapStr2doubleVec& DoubleFeatureData, const string& wildcard,
               vector<string>& traces);

#endif
