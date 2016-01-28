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

#ifndef __UTILS
#define __UTILS

#include <vector>

using std::vector;

int LinearInterpolation(double dt, const vector<double>& X,
                        const vector<double>& Y, vector<double>& InterpX,
                        vector<double>& InterpY);
int getCentralDifferenceDerivative(double dx, const vector<double>& v,
                                   vector<double>& dv);
void getfivepointstencilderivative(const vector<double>& v, vector<double>& dv);
void slope_straight_line_fit(const vector<double>& x, const vector<double>& y,
                             vector<double>& slope);

template <class ForwardIterator>
ForwardIterator first_min_element(ForwardIterator first, ForwardIterator last) {
  ForwardIterator lowest = first;
  ForwardIterator loop = first;
  ForwardIterator lowestMin = first;
  int minCounter = 0;
  if (first == last) return last;
  while (++loop != last) {
    if (*loop < *lowest) {
      lowest = loop;
    }
    if (*loop >= *lowestMin) {
      minCounter++;
    }
    if (*loop < *lowestMin) {
      minCounter = 0;
      lowestMin = loop;
    }
    if (minCounter == 2) return lowestMin;
  }
  return lowest;
}

#endif
