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

#ifndef __DEPENDENCYTREE_H
#define __DEPENDENCYTREE_H
#include <iostream>
#include <stdio.h>
using namespace std;
#include <map>
#include <string>
#include <vector>
#include <list>

typedef int (*fptr)(map<string, vector<int> > &, map<string, vector<double> > &,
                    map<string, string> &);
extern map<string, fptr> FptrTableV1;

class cTree {
  vector<string> strDependencyFile;
  vector<string> vecFeature;

 public:
  string ErrorStr;
  list<string> FinalList;
  list<string> ChildList;
  cTree() { printf("\n This is constructor of cTree"); }
  cTree(const char *strFileName);
  int printTree();
  int getDependencyList(string str);
  int setFeaturePointers(map<string, map<string, fptr> *> &mapFptrLib,
                         map<string, fptr> *FptrTable,
                         map<string, vector<fptr> > *FptrLookup);
  int setFeaturePointers(map<string, map<string, fptr> *> &mapFptrLib,
                         map<string, fptr> *FptrTable,
                         map<string, vector<pair<fptr, string> > > *FptrLookup);
  int getChilds(string strLine, list<string> &childs);
  int getDependency(string strLine, string parent_stim);
  int deblank(string &str);
  int AddUniqueItem(string strFeature, list<string> &lstFinal);
  int getAllParents(vector<string> &vecFeature);
};

#endif
