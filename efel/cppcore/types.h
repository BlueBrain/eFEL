/* Copyright (c) 2016, EPFL/Blue Brain Project
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

#ifndef EFELTYPES_H
#define EFELTYPES_H

#include <map>
#include <string>
#include <utility>
#include <vector>


typedef std::map<std::string, std::vector<int> > mapStr2intVec;
typedef std::map<std::string, std::vector<double> > mapStr2doubleVec;
typedef std::map<std::string, std::string> mapStr2Str;

typedef int (*feature_function)(mapStr2intVec &,
                                mapStr2doubleVec &,
                                mapStr2Str &);

/* Name of the feature (ie: 'interpolate') -> function pointer for feature
 */
typedef std::map<std::string, feature_function> feature2function;

typedef std::pair<feature_function, std::string> featureStringPair;

#endif
