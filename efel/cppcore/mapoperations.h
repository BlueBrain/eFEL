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

#include <stdexcept>
#include <string>
#include <vector>

using std::string;
using std::vector;

extern string GErrorStr;

template <typename T>
std::map<std::string, std::vector<T>> getFeatures(
    const std::map<std::string, std::vector<T>>& allFeatures,
    const std::vector<std::string>& requestedFeatures);
template<typename T>
int getParam(std::map<std::string, std::vector<T>>& featureData,
               const std::string& param, std::vector<T>& vec);
int getStrParam(mapStr2Str& StringData, const string& param, string& value);
template <class T>
void setVec(std::map<std::string, std::vector<T> >& FeatureData, mapStr2Str& StringData,
               string key, const vector<T>& value);
template <class T>
int getVec(std::map<std::string, std::vector<T> >& FeatureData, mapStr2Str& StringData,
                 string strFeature, vector<T>& v);
#endif
