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

#ifndef GLOBAL_H
#define GLOBAL_H

#include <map>
#include <string>
#include <vector>

using namespace std;
typedef int (*fptr)(map<string, vector<int> > &, map<string, vector<double> > &,
                    map<string, string> &);

map<string, fptr> FptrTableV1;
map<string, fptr> FptrTableV2;
map<string, fptr> FptrTableV3;
map<string, fptr> FptrTableV4;
map<string, fptr> FptrTableV5;
map<string, fptr> FptrTable;
map<string, map<string, fptr> *> mapFptrLib;
string GErrorStr;

#endif
