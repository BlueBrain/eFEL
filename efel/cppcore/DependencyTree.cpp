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

#include "DependencyTree.h"

#include <cstdio>
#include <fstream>

cTree::cTree() {
  printf("\n This is constructor of cTree"); 
}

cTree::cTree(const char *strFileName) {
  char strLine[3000];
  std::ifstream infile(strFileName);
  if (infile.is_open()) {
    infile.getline(strLine, 3000, '\n');
    strDependencyFile.push_back(strLine);
    while (!infile.eof()) {
      // printf("%s\n", strLine );
      infile.getline(strLine, 3000, '\n');
      strDependencyFile.push_back(strLine);
    }
    infile.close();
  } else {
    ErrorStr = ErrorStr + string("\nCould not open the file ") + strFileName;
  }
  getAllParents(vecFeature);
}
int cTree::getDependencyList(string) {
  for (unsigned i = 0; i < strDependencyFile.size(); i++) {
  }
  return 1;
}

int cTree::printTree() {
  printf("\n This is inside printTree ");
  return 1;
}

int cTree::deblank(string &str) {
  size_t startpos = str.find_first_not_of(" \t");
  size_t endpos = str.find_last_not_of(" \t");
  if ((string::npos == startpos) || (string::npos == endpos))
    str = "";
  else
    str = str.substr(startpos, endpos - startpos + 1);
  return 1;
}

/**
 *
 * setFeaturePointers
 *
 * mapFptrLib :
 * FptrTable :
 * FptrLookup | vector of pairs:
 *                  first | string: feature name
 *                  second | vector of featureStringPair
 *
 */
int cTree::setFeaturePointers(map<string, feature2function *> &mapFptrLib,
                              feature2function *FptrTable,
                              map<string, vector<featureStringPair > > *FptrLookup) 
{
  list<string>::iterator lstItr;
  map<string, feature2function *>::iterator mapLibItr;
  feature2function *fptrTbl;
  feature2function::iterator mapFeatureItr;

  string strLibFeature, strLib, strFeature;
  string wildcards;

  vector<featureStringPair> vecfptr;

  int nPos;
  int wcpos;

  if (vecFeature.size() == 0) return -1;

  // vecFeature is a list with all the feature names in the first column
  // of the dependency file
  for (unsigned i = 0; i < vecFeature.size(); i++) {

    FinalList.clear();
    strLibFeature = vecFeature[i];
    // fill FinalList with all the dependencies of feature vecFeature[i]
    getDependency(strLibFeature, "");
    vecfptr.clear();
    for (lstItr = FinalList.begin(); lstItr != FinalList.end(); lstItr++) {
      // Here strLibFeature is the feature name of the dependent feature
      strLibFeature = *lstItr;
      nPos = strLibFeature.find(":");
      wcpos = strLibFeature.find(";");
      if (nPos <= 0) {
        ErrorStr += "\nLibrary version is missing in [" + strLibFeature +
                    "] expected format Lib:Feature";
        return -1;
      }
      // Library name of the feature
      strLib = strLibFeature.substr(0, nPos);

      // Put feature name in strFeature
      if (wcpos < 0) {
        strFeature = strLibFeature.substr(nPos + 1);
        wildcards = "";
      } else {
        strFeature = strLibFeature.substr(nPos + 1, wcpos - nPos - 1);
        wildcards = strLibFeature.substr(wcpos);
      }

      // Find the feature function pointer map for the library
      mapLibItr = mapFptrLib.find(strLib);
      if (mapLibItr == mapFptrLib.end()) {
        ErrorStr = ErrorStr + string("\nLibrary [") + strLib + "] is missing\n";
        return (-1);
      }

      // Find the function pointer of the feature in the library
      fptrTbl = (mapLibItr->second);
      mapFeatureItr = (fptrTbl)->find(strFeature);
      if (mapFeatureItr == (fptrTbl)->end()) {
        ErrorStr = ErrorStr + string("\nFeature [") + strFeature +
                   string("] is missing from Library [") + strLib + "]";
        return -1;
      }

      // Add the feature function pointer and wildcards to the list of dependent
      // features
      vecfptr.push_back(featureStringPair(mapFeatureItr->second, wildcards));
      FptrTable->insert(std::pair<string, feature_function>(
              strFeature, mapFeatureItr->second));
    }
    // Add the vecfptr from above to a map with as key the base featurei
    FptrLookup->insert(
        std::pair<string, vector<featureStringPair> >(strFeature, vecfptr));
  }

  return 1;
}

/**
 *
 * Fill lstFeature (vecFeature) with all the features in the first column
 * of the dependency file
 *
 */
int cTree::getAllParents(vector<string> &lstFeature) {
  string strLine, FeatureName;
  int nPos;
  for (unsigned i = 0; i < strDependencyFile.size(); i++) {
    strLine = strDependencyFile[i];
    deblank(strLine);
    nPos = int(strLine.find_first_of('#'));
    FeatureName = strLine.substr(0, nPos - 1);
    deblank(FeatureName);
    if (FeatureName.size() != 0) lstFeature.push_back(FeatureName);
  }
  return 1;
}
int cTree::getChilds(string str, list<string> &childs) {
  string strLine, FeatureName;
  int nPos;
  deblank(str);
  for (unsigned i = 0; i < strDependencyFile.size(); i++) {
    strLine = strDependencyFile[i];
    deblank(strLine);
    nPos = int(strLine.find_first_of('#'));
    FeatureName = strLine.substr(0, nPos - 1);
    deblank(FeatureName);
    if (FeatureName != str)
      continue;
    else {
      int nPos = int(strLine.find_first_of('#'));
      string Token;
      while (nPos >= 0) {
        strLine = strLine.substr(nPos + 1);
        nPos = int(strLine.find_first_of('#'));
        if (nPos > 0)
          Token = strLine.substr(0, nPos - 1);
        else
          Token = strLine;
        deblank(Token);
        childs.push_back(Token);
      }
    }
  }
  return 1;
}

/*
 *
 *  Fill FinalList with a list of all the feature matching the wildcards
 *
 */
int cTree::getDependency(string strLine, string wildcards) {
  list<string> tmpChild;

  // parse wildcards out of "LibVx:feature_name;wildcards_name"
  int wcpos = strLine.find(";");
  if (wcpos >= 0) {
    wildcards = strLine.substr(wcpos);
    strLine = strLine.substr(0, wcpos);
  }
  unsigned childCount = 0;
  getChilds(strLine, tmpChild);
  if (tmpChild.size() != 0) {
    childCount = tmpChild.size();
    for (unsigned i = 0; i < childCount; i++) {
      string str = tmpChild.front();
      tmpChild.pop_front();
      getDependency(str, wildcards);
    }
  }
  AddUniqueItem(strLine + wildcards, FinalList);
  return 0;
}

int cTree::AddUniqueItem(string strFeature, list<string> &lstFinal) {
  list<string>::iterator lstItr;
  int FoundFlag = 0;
  for (lstItr = lstFinal.begin(); lstItr != lstFinal.end(); lstItr++) {
    if (strFeature == *lstItr) {
      FoundFlag = 1;
      break;
    }
  }
  if (!FoundFlag) lstFinal.push_back(strFeature);
  return 1;
}
