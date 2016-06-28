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

#include "LibV4.h"
#include <algorithm>
#include <functional>
#include <iterator>
#include <math.h>

// new algorithm to find the spikes in the voltage trace
// based on minimum spike height
//
// get all local minima
// assume spikes between local minima
// check if the assumed spike is bigger than min_spike_height
static int __peak_indices(const vector<double>& v, double min_spike_height,
                          double threshold, vector<int>& peakindices) {
  vector<double> dv;
  vector<int> minimum_indices;
  getCentralDifferenceDerivative(1., v, dv);
  int falling_index = 0;
  int minimum_index = 0;
  do {
    minimum_indices.push_back(minimum_index);
    falling_index =
        distance(dv.begin(), find_if(dv.begin() + minimum_index, dv.end(),
                                     std::bind2nd(std::less<double>(), 0.)));
    minimum_index =
        distance(dv.begin(), find_if(dv.begin() + falling_index, dv.end(),
                                     std::bind2nd(std::greater_equal<double>(), 0.)));
  } while (dv.begin() + minimum_index != dv.end());
  minimum_indices.push_back(dv.size() - 1);

  for (unsigned i = 0; i < minimum_indices.size() - 1; i++) {
    int maximum_index =
        distance(v.begin(), max_element(v.begin() + minimum_indices[i],
                                        v.begin() + minimum_indices[i + 1]));
    double spike_height1 = v[maximum_index] - v[minimum_indices[i + 1]];
    double spike_height2 = v[maximum_index] - v[minimum_indices[i]];
    if ((v[maximum_index] > threshold && (spike_height1 > min_spike_height ||
                                          spike_height2 > min_spike_height)) ||
        (spike_height1 > min_spike_height &&
         spike_height2 > min_spike_height)) {
      // if(spike_height1 > min_spike_height && spike_height2 >
      // min_spike_height) {
      peakindices.push_back(maximum_index);
    }
  }
  return peakindices.size();
}

int LibV4::peak_indices(mapStr2intVec& IntFeatureData,
                        mapStr2doubleVec& DoubleFeatureData,
                        mapStr2Str& StringData) {
  int size;
  if (CheckInIntmap(IntFeatureData, StringData, "peak_indices", size)) {
    return size;
  }

  vector<int> peakindices;
  vector<double> v;
  vector<double> min_spike_height;
  vector<double> threshold;
  if (getDoubleVec(DoubleFeatureData, StringData, "V", v) <= 0) {
    return -1;
  }
  if (getDoubleParam(DoubleFeatureData, "min_spike_height", min_spike_height) <=
      0) {
    return -1;
  }
  if (getDoubleParam(DoubleFeatureData, "Threshold", threshold) <= 0) {
    return -1;
  }

  int retval = __peak_indices(v, min_spike_height[0], threshold[0], peakindices);
  if (retval >= 0) {
    setIntVec(IntFeatureData, StringData, "peak_indices", peakindices);
    return peakindices.size();
  }

  return retval;
}
