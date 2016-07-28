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

int LinearInterpolation(double Stepdx, const vector<double>& X,
                        const vector<double>& Y, vector<double>& InterpX,
                        vector<double>& InterpY) {
  unsigned nCount = Y.size();
  int nPts = ceil((X[nCount - 1] - X[0]) / Stepdx) + 1;  // Because time is in
                                                    // millisecond and needs to
                                                    // be interpolated at 0.1 ms
                                                    // interval
  double input = X[0];
  unsigned int i = 1;
  double dif1, dif2;
  InterpY.push_back(Y[0]);
  InterpX.push_back(X[0]);
  for (int j = 1; j < nPts; j++) {
    input = input + Stepdx;
    while ((X[i] < input) && (i < nCount)) i++;
    dif1 = X[i] - X[i - 1];  //!=0 per definition
    dif2 = input - X[i - 1];
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

// fit a straight line to the points (x[i],y[i]) and return the slope y'(x)
//
// slope[0] = slope
// slope[1] = average residual sum squares
// slope[2] = coefficient of determination R^2
void slope_straight_line_fit(const vector<double>& x, const vector<double>& y,
                             vector<double>& slope) {
  slope.resize(3);
  if (x.size() != y.size()) {
    printf("Unequal vectors in straight line fit\n");
    slope[0] = 1.;
    slope[1] = 1000.;
    return;
  }
  double sum_x = 0.;
  double sum_y = 0.;
  double sum_x2 = 0.;
  double sum_xy = 0.;
  for (unsigned i = 0; i < x.size(); i++) {
    sum_x += x[i];
    sum_y += y[i];
    sum_x2 += x[i] * x[i];
    sum_xy += x[i] * y[i];
  }
  double delta = x.size() * sum_x2 - sum_x * sum_x;
  slope[0] = (x.size() * sum_xy - sum_x * sum_y) / delta;
  //
  // calculate sum of squared residuals
  double yintercept = (sum_y - slope[0] * sum_x) / x.size();
  double residuals = 0.;
  for (unsigned i = 0; i < x.size(); i++) {
    double res = y[i] - yintercept - slope[0] * x[i];
    residuals += res * res;
  }
  slope[1] = residuals / x.size();
  // calculate the coefficient of determination R^2
  double y_av = sum_y / x.size();
  double sstot = 0.;
  for (unsigned i = 0; i < x.size(); i++) {
    double dev = y[i] - y_av;
    sstot += dev * dev;
  }
  slope[2] = 1. - residuals / sstot;
}
