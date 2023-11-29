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

#include "LibV3.h"

#include <math.h>

#include <algorithm>
#include <functional>
#include <list>

using std::find_if;
using std::list;
using std::min_element;

static int __depolarized_base(const vector<double>& t, const vector<double>& v,
                              double stimstart, double stimend,
                              const vector<int>& apbi,
                              const vector<int>& apendi,
                              vector<double>& dep_base) {
  int i, n, k, startIndex, endIndex, nPt;
  double baseValue;
  // to make sure it access minimum index of both length
  if (apendi.size() < apbi.size())
    n = apendi.size();
  else
    n = apbi.size();

  if (apendi.size() == apbi.size()) n = apendi.size() - 1;

  if (n > 2) {
    dep_base.clear();
    for (i = 0; i < n; i++) {
      nPt = 0;
      baseValue = 0;
      startIndex = apendi[i];
      endIndex = apbi[i + 1];
      for (k = startIndex; k < endIndex; k++) {
        baseValue += v[k];
        nPt++;
      }
      baseValue = baseValue / nPt;
      dep_base.push_back(baseValue);
    }
    return dep_base.size();
  }
  return -1;
}

int LibV3::depolarized_base(mapStr2intVec& IntFeatureData,
                            mapStr2doubleVec& DoubleFeatureData,
                            mapStr2Str& StringData) {
  // Retrieve all required double and int features at once
  const auto& doubleFeatures =
      getFeatures(DoubleFeatureData, {"T", "V", "stim_start", "stim_end"});
  const auto& intFeatures =
      getFeatures(IntFeatureData, {"AP_end_indices", "AP_begin_indices"});

  vector<double> dep_base;
  int retVal = __depolarized_base(
      doubleFeatures.at("T"), doubleFeatures.at("V"),
      doubleFeatures.at("stim_start").front(),
      doubleFeatures.at("stim_end").front(), intFeatures.at("AP_begin_indices"),
      intFeatures.at("AP_end_indices"), dep_base);

  if (retVal > 0) {
    setVec(DoubleFeatureData, StringData, "depolarized_base", dep_base);
  }
  return retVal;
}
