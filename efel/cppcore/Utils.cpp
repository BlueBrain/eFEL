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

#include "Utils.h"
#include <algorithm>
#include <cstdio>
#include <iostream>
#include <iostream>
#include <iterator>
#include <math.h>
#include <assert.h>

int LinearInterpolation(double Stepdx,
                        const vector<double>& X,
                        const vector<double>& Y,
                        vector<double>& InterpX,
                        vector<double>& InterpY) {
  EFEL_ASSERT(X.size() == Y.size(), "X & Y have to have the same point count");
  EFEL_ASSERT(2 < X.size(), "Need at least 2 points in X");
  assert(Stepdx != 0);
  
  size_t nCount = X.size() - 1;
  size_t i = 1;

  double input = X[0];
  double dif1, dif2;

  InterpY.push_back(Y[0]);
  InterpX.push_back(X[0]);

  while(input < X[nCount]){
    input += Stepdx;

    while (X[i] < input && i < nCount) {
      i++;
    }

    dif1 = X[i] - X[i - 1];
    dif2 = input - X[i - 1];
    assert(dif1 != 0); //!=0 per definition

    InterpY.push_back(Y[i - 1] + ((Y[i] - Y[i - 1]) * dif2 / dif1));
    InterpX.push_back(input);
  }
  return 1;
}

int getCentralDifferenceDerivative(double dx, const vector<double>& v,
                                   vector<double>& dv) {
  unsigned n = v.size();
  dv.clear();
  // because formula is ((vec[i+1]+vec[i-1])/2)/dx hence it should iterate
  // through 1 to length-1
  dv.push_back((v[1] - v[0]) / dx);
  for (unsigned i = 1; i < n - 1; i++) {
    dv.push_back(((v[i + 1] - v[i - 1]) / 2) / dx);
  }
  dv.push_back((v[n - 1] - v[n - 2]) / dx);
  return 1;
}

void getfivepointstencilderivative(const vector<double>& v,
                                   vector<double>& dv) {
  dv.clear();
  dv.resize(v.size());
  dv[0] = v[1] - v[0];
  dv[1] = (v[2] - v[0]) / 2.;
  for (unsigned i = 2; i < v.size() - 2; i++) {
    dv[i] = -v[i + 2] + 8 * v[i + 1] - 8 * v[i - 1] + v[i - 2];
    dv[i] /= 12.;
  }
  dv[v.size() - 2] = (v[v.size() - 1] - v[v.size() - 3]) / 2.;
  dv[v.size() - 1] = v[v.size() - 1] - v[v.size() - 2];
}

// fit a straight line to the points (x[i], y[i]) and return the slope y'(x)
linear_fit_result
slope_straight_line_fit(const vector<double>& x,
                        const vector<double>& y
                        ) {
  EFEL_ASSERT(x.size() == y.size(), "X & Y have to have the same point count");
  EFEL_ASSERT(1 <= x.size(), "Need at least 1 points in X");

  double sum_x = 0.;
  double sum_y = 0.;
  double sum_x2 = 0.;
  double sum_xy = 0.;

  linear_fit_result result;

  for (unsigned i = 0; i < x.size(); i++) {
    sum_x += x[i];
    sum_y += y[i];
    sum_x2 += x[i] * x[i];
    sum_xy += x[i] * y[i];
  }

  double delta = x.size() * sum_x2 - sum_x * sum_x;
  result.slope = (x.size() * sum_xy - sum_x * sum_y) / delta;

  // calculate sum of squared residuals
  double yintercept = (sum_y - result.slope * sum_x) / x.size();
  double residuals = 0.;
  for (unsigned i = 0; i < x.size(); i++) {
    double res = y[i] - yintercept - result.slope * x[i];
    residuals += res * res;
  }
  result.average_rss = residuals / x.size();

  // calculate the coefficient of determination R^2
  double y_av = sum_y / x.size();
  double sstot = 0.;
  for (unsigned i = 0; i < x.size(); i++) {
    double dev = y[i] - y_av;
    sstot += dev * dev;
  }
  result.r_square = 1. - residuals / sstot;

  return result;
}
