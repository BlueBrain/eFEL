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

#ifndef __LIBV4
#define __LIBV4
#include <iterator>
#include <algorithm>
#include <map>
#include <string>
#include <vector>
#include <list>
#include "mapoperations.h"
#include "Utils.h"

using namespace std;

typedef map<string, vector<int> > mapStr2intVec;
typedef map<string, vector<double> > mapStr2doubleVec;
typedef map<string, string> mapStr2Str;
namespace LibV4 {
int peak_indices(mapStr2intVec& IntFeatureData,
                 mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData);
int __peak_indices(const vector<double>& v, double min_spike_height,
                   double threshold, vector<int>& peakindices);
}
#endif
