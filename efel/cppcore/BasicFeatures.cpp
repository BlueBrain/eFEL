/* Copyright (c) 2015-2024, EPFL/Blue Brain Project                                   
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

#include "BasicFeatures.h"


int BasicFeatures::interpolate(mapStr2intVec& IntFeatureData,
                       mapStr2doubleVec& DoubleFeatureData,
                       mapStr2Str& StringData) {
  vector<double> V, T, VIntrpol, TIntrpol, InterpStepVec;
  T = getFeature(DoubleFeatureData, "T");
  // interp_step is a stimulus independent parameter
  int retVal = getParam(DoubleFeatureData, "interp_step", InterpStepVec);
  double InterpStep = (retVal <= 0) ? 0.1 : InterpStepVec[0];

  try  // interpolate V if it's available
  {
    V = getFeature(DoubleFeatureData, "V");
    LinearInterpolation(InterpStep, T, V, TIntrpol, VIntrpol);
    setVec(DoubleFeatureData, StringData, "V", VIntrpol);
    setVec(DoubleFeatureData, StringData, "T", TIntrpol);
  } catch (...) {
    return -1;  // interpolation failed
  }

  // also interpolate current if present
  vector<double> I, IIntrpol, TIntrpolI;
  try {
    I = getFeature(DoubleFeatureData, "I");
    LinearInterpolation(InterpStep, T, I, TIntrpolI, IIntrpol);
    setVec(DoubleFeatureData, StringData, "I", IIntrpol);
    setVec(DoubleFeatureData, StringData, "T", TIntrpol);
  } catch (...) {
  }  // pass, it is optional
  return 1;
}

// return (possibly interpolate) voltage trace
int BasicFeatures::voltage(mapStr2intVec& IntFeatureData,
                   mapStr2doubleVec& DoubleFeatureData,
                   mapStr2Str& StringData) {
  const vector<double>& v = getFeature(DoubleFeatureData, "V");
  setVec(DoubleFeatureData, StringData, "voltage", v);
  return v.size();
}

// return (possibly interpolate) current trace
int BasicFeatures::current(mapStr2intVec& IntFeatureData,
                   mapStr2doubleVec& DoubleFeatureData,
                   mapStr2Str& StringData) {
  const vector<double>& i = getFeature(DoubleFeatureData, "I");
  setVec(DoubleFeatureData, StringData, "current", i);
  return i.size();
}

// return (possibly interpolate) time trace
int BasicFeatures::time(mapStr2intVec& IntFeatureData,
                mapStr2doubleVec& DoubleFeatureData, mapStr2Str& StringData) {
  const vector<double>& t = getFeature(DoubleFeatureData, "T");
  setVec(DoubleFeatureData, StringData, "time", t);
  return t.size();
}
