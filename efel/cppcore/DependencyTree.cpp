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

#include <algorithm>  //remove
#include <cctype>     //isspace
#include <fstream>

static void removeAllWhiteSpace(string &str) {
  str.erase(std::remove_if(str.begin(), str.end(), isspace), str.end());
}

cTree::cTree(const char *strFileName) {
  std::ifstream input(strFileName);

  if (!input) {
    ErrorStr += "\nCould not open the file " + string(strFileName);
    return;
  }

  for (string line; std::getline(input, line);) {
    removeAllWhiteSpace(line);
    if (!line.empty()) {
      strDependencyFile.push_back(std::move(line));
    }
  }

  getAllParents(vecFeature);
}

/**
 *
 * setFeaturePointers
 *
 * mapFptrLib :
 * FptrTable :
 * FptrLookup | vector of pairs:
 *                  first | string: feature name
 *                  second | vector of feature_function to represent list of
 * dependent features
 *
 */
int cTree::setFeaturePointers(
    map<string, feature2function *> &mapFptrLib, feature2function *FptrTable,
    map<string, vector<feature_function> > *FptrLookup) {
  list<string>::iterator lstItr;
  map<string, feature2function *>::iterator mapLibItr;
  feature2function *fptrTbl;
  feature2function::iterator mapFeatureItr;

  string strLibFeature, strLib, strFeature;

  vector<feature_function> vecfptr;

  if (vecFeature.size() == 0) return -1;

  // vecFeature is a list with all the feature names in the first column
  // of the dependency file
  for (unsigned i = 0; i < vecFeature.size(); i++) {
    FinalList.clear();
    strLibFeature = vecFeature[i];
    // fill FinalList with all the dependencies of feature vecFeature[i]
    getDependency(strLibFeature);
    vecfptr.clear();
    for (lstItr = FinalList.begin(); lstItr != FinalList.end(); lstItr++) {
      // Here strLibFeature is the feature name of the dependent feature
      strLibFeature = *lstItr;
      size_t nPos = strLibFeature.find(":");
      if (nPos == string::npos || nPos == 0) {
        ErrorStr += "\nLibrary version is missing in [" + strLibFeature +
                    "] expected format Lib:Feature";
        return -1;
      }
      // Library name of the feature
      strLib = strLibFeature.substr(0, nPos);

      // Put feature name in strFeature
      strFeature = strLibFeature.substr(nPos + 1);

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

      vecfptr.push_back(mapFeatureItr->second);
      FptrTable->insert(std::pair<string, feature_function>(
          strFeature, mapFeatureItr->second));
    }
    // Add the vecfptr from above to a map with as key the base featurei
    FptrLookup->insert(
        std::pair<string, vector<feature_function> >(strFeature, vecfptr));
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
  for (unsigned i = 0; i < strDependencyFile.size(); i++) {
    const string &strLine = strDependencyFile[i];
    size_t nPos = strLine.find_first_of('#');
    string FeatureName = strLine.substr(0, nPos);
    if (!FeatureName.empty()) {
      lstFeature.push_back(FeatureName);
    }
  }
  return 1;
}

int cTree::getChilds(string str, list<string> &childs) {
  string strLine, FeatureName;
  for (unsigned i = 0; i < strDependencyFile.size(); i++) {
    strLine = strDependencyFile[i];
    size_t nPos = strLine.find_first_of('#');
    if (strLine.substr(0, nPos) != str)
      continue;
    else {
      size_t nPos = strLine.find_first_of('#');
      string Token;
      while (nPos != string::npos) {
        strLine = strLine.substr(nPos + 1);
        nPos = strLine.find_first_of('#');
        Token = strLine.substr(0, nPos);
        childs.push_back(Token);
      }
    }
  }
  return 1;
}

int cTree::getDependency(const string &strLine) {
  std::list<string> tmpChild;

  getChilds(strLine, tmpChild);
  for (const auto &childFeature : tmpChild) {
    getDependency(
        childFeature);  // Recursively get dependencies of the child feature.
  }
  AddUniqueItem(strLine);  // Add the feature itself to the FinalList.
  return 0;
}

void cTree::AddUniqueItem(const string &strFeature) {
  auto it = std::find(FinalList.begin(), FinalList.end(), strFeature);
  if (it == FinalList.end()) {
    FinalList.push_back(strFeature);
  }
}
