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

#include <algorithm> //remove
#include <cctype> //isspace
#include <fstream>

static void removeAllWhiteSpace(string &str) {
  str.erase(std::remove_if(str.begin(), str.end(), isspace), str.end());
}

cTree::cTree(const char *strFileName) {
  std::string line;

  std::ifstream input(strFileName);
  if (input.is_open()) {

    std::getline(input, line);
    removeAllWhiteSpace(line);
    if (!line.empty())
      strDependencyFile.push_back(line);
    while (!(input.fail() || input.eof())) {
      std::getline(input, line);
      removeAllWhiteSpace(line);
      if (!line.empty())
        strDependencyFile.push_back(line);
    }
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
      size_t nPos = strLibFeature.find(":");
      if (nPos == string::npos || nPos == 0) {
        ErrorStr += "\nLibrary version is missing in [" + strLibFeature +
                    "] expected format Lib:Feature";
        return -1;
      }
      // Library name of the feature
      strLib = strLibFeature.substr(0, nPos);

      // Put feature name in strFeature
      size_t wcpos = strLibFeature.find(";");
      if (wcpos == string::npos) {
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
  for (unsigned i = 0; i < strDependencyFile.size(); i++) {
    const string& strLine = strDependencyFile[i];
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

/*
 *
 *  Fill FinalList with a list of all the feature matching the wildcards
 *
 */
int cTree::getDependency(string strLine, string wildcards) {
  list<string> tmpChild;

  // parse wildcards out of "LibVx:feature_name;wildcards_name"
  size_t wcpos = strLine.find(";");
  if (wcpos != string::npos) {
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
  bool FoundFlag = false;
  for (lstItr = lstFinal.begin(); lstItr != lstFinal.end(); lstItr++) {
    if (strFeature == *lstItr) {
      FoundFlag = true;
      break;
    }
  }
  if (!FoundFlag) lstFinal.push_back(strFeature);
  return 1;
}
