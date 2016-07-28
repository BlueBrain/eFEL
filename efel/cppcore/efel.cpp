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

#include <time.h>
#include "efel.h"
#include "cfeature.h"

cFeature *pFeature = NULL;

int Initialize(const char *strDepFile, const char *outdir) {
  if (pFeature != NULL) {
    delete pFeature;
  }

  pFeature = new cFeature(string(strDepFile), string(outdir));
  if (pFeature == NULL) {
    return -1;
  } else {
    return 1;
  }
}

int setVersion(const char *strDepFile) {
  return pFeature->setVersion(string(strDepFile));
}

int FeaturePrint(const char *strFile) {
  pFeature->printFeature(strFile);
  return 1;
}
int setFeatureInt(const char *strName, int *A, unsigned nValue) {
  vector<int> v(nValue);
  for (unsigned i = 0; i < nValue; i++) {
    v[i] = A[i];
  }
  pFeature->setFeatureInt(string(strName), v);
  return 1;
}

int setFeatureDouble(const char *strName, double *A, unsigned nValue) {
  // printf("\nInside featureLibrary.. Before setdouble [%s = %f add = %d nVal=
  // %d]\n", strName, A[0], pFeature, nValue);
  vector<double> v(nValue);
  for (unsigned i = 0; i < nValue; i++) {
    v[i] = A[i];
  }
  // mapDoubleData.insert(pair<string, vector< double > > (string(strName), v));
  pFeature->setFeatureDouble(string(strName), v);
  // printf("\nInside featureLibrary.. After setdouble [%s = %f]\n", strName,
  // A[0]);
  return 1;
}

int getFeatureInt(const char *strName, int **A) {
  vector<int> vec;
  if (pFeature->getFeatureInt(string(strName), vec) < 0) {
    return -1;
  }
  *A = new int[vec.size()];
  for (unsigned i = 0; i < vec.size(); i++) {
    (*A)[i] = vec[i];
  }
  return vec.size();
}

int getFeatureDouble(const char *strName, double **A) {
  vector<double> vec;
  // printf("\nInside featureLibrary.. Before getdouble [%s ]\n", strName);
  if (pFeature->getFeatureDouble(string(strName), vec) < 0) {
    return -1;
  }
  *A = new double[vec.size()];
  for (unsigned i = 0; i < vec.size(); i++) {
    (*A)[i] = vec[i];
  }
  // printf("\nInside featureLibrary.. After getdouble [%s= %f ]\n", strName,
  // (*A)[0]);
  return vec.size();
}

int getFeatureString(const char *key, char *&value) {
  string strvalue;
  pFeature->getFeatureString(key, strvalue);
  value = new char[strvalue.length() + 1];
  copy(strvalue.begin(), strvalue.end(), value);
  value[strvalue.length()] = '\0';
  return 1;
}

int setFeatureString(const char *key, const char *value) {
  pFeature->setFeatureString(key, value);
  return 1;
}

char *getgError() {
  string error = GErrorStr + pFeature->getGError();
  GErrorStr.clear();
  return (char *)error.c_str();
}

double getDistance(const char *strName, double mean, double std, 
        bool trace_check) {
  double value;
  value = pFeature->getDistance(string(strName), mean, std, trace_check);
  return value;
}

int printFptr() {
  printf("\n size of fptrlookup %d", (int)pFeature->fptrlookup.size());
  return 1;
}
