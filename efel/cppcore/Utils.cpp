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
  EFEL_ASSERT(Stepdx > 0, "Interpolation step needs to be strictly positive");
 
  double dx, dy, dydx;
  size_t InterpX_size;
  double x = X[0];
  double start = X[0];
  double stop = X[X.size() - 1] + Stepdx;
  
  // Inspired by the way numpy.arange works
  // Do not remove the 'ceil' in favor of < stop in for loop
  InterpX_size = ceil((stop - start)/Stepdx);

  for (size_t i = 0; i < InterpX_size; i++) {
      InterpX.push_back(x);
      x += Stepdx;
  }

  // Create the y values
  size_t j = 0;
  for (size_t i = 0; i < InterpX.size(); i++) {
    x = InterpX[i];

    EFEL_ASSERT((j+1) < X.size(), "Interpolation accessing point outside of X");
    
    while ( X[j+1] < x ) {
        j++;
        if (j+1 >= X.size()) {
            j = X.size() - 1;
            break;
        }
        EFEL_ASSERT((j+1) < X.size(), 
                "Interpolation accessing point outside of X");
    }
    


    if (j == X.size() - 1) {
        // Last point
        InterpY.push_back(Y[j]);
        break;
    } 
    else {
        EFEL_ASSERT((j+1) < X.size(), 
                "Interpolation accessing point outside of X");
        
        dx = X[j+1] - X[j];
        dy = Y[j+1] - Y[j];

        EFEL_ASSERT(dx != 0,  "Interpolation using dx == 0"); //!=0 per definition

        dydx = dy/dx;

        InterpY.push_back(Y[j] + dydx * (x - X[j]));
    }
  }

  return 1;
}

int getCentralDifferenceDerivative(double dx, const vector<double>& v,
                                   vector<double>& dv) {
  size_t n = v.size();
  dv.clear();
  // because formula is ((vec[i+1]+vec[i-1])/2)/dx hence it should iterate
  // through 1 to length-1
  dv.push_back((v[1] - v[0]) / dx);
  for (size_t i = 1; i < n - 1; i++) {
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
  for (size_t i = 2; i < v.size() - 2; i++) {
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

  for (size_t i = 0; i < x.size(); i++) {
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
  for (size_t i = 0; i < x.size(); i++) {
    double res = y[i] - yintercept - result.slope * x[i];
    residuals += res * res;
  }
  result.average_rss = residuals / x.size();

  // calculate the coefficient of determination R^2
  double y_av = sum_y / x.size();
  double sstot = 0.;
  for (size_t i = 0; i < x.size(); i++) {
    double dev = y[i] - y_av;
    sstot += dev * dev;
  }
  result.r_square = 1. - residuals / sstot;

  return result;
}

std::pair<size_t, size_t> get_time_index(std::vector<double> &t, double startTime,
                                     double endTime, double precisionThreshold) {
  /*
   * Returns the start and end index of the time array 
   * Uses a threshold to tolerate the precision loss
   * */                                    

  size_t startIndex = 0;
  size_t endIndex = t.size();

  for (size_t i = 0; i < t.size(); i++){
    if (t[i] >= startTime)
    {
      startIndex = i;
      break;
    }
  }

  for (size_t i = t.size() - 1; i > 0; i--){ // backward iterator using indices
    if (t[i] - endTime < precisionThreshold)
    {
      endIndex = i+1;
      break;
    }
  }

  return std::pair<size_t, size_t>(startIndex, endIndex);
                                      
}

template<class T>
double vec_mean(const vector<T> &v) {
  /*
   * Computes the mean of input vector v
   * Does not modify the reference v
   * Returns the mean value
   * */

    double sum = accumulate( v.begin(), v.end(), 0.0);
    size_t v_size = v.size();
    double mean = sum / v_size;
    return mean;
}

template<class T>
double vec_median(vector<T> v) {
    /*
     * param v: the input vector
     * param v is called by value since std::sort modifies the container
     * Returns the median value of the vector
     * */

    // sort using the default operator<
    std::sort(v.begin(), v.end());

    size_t n = v.size();
    if (n % 2 != 0) // odd
        return (double)v[n/2];
    return (double)(v[(n-1)/2] + v[n/2])/2.0; // even

}


template double vec_mean(const vector<double> &v);
template double vec_median(vector<double> v);
