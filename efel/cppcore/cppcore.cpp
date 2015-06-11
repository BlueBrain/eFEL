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

/*
 * =============================================================================
 *
 *       Filename:  CppCore.cpp
 *
 *    Description:  Python wrapper file, to be compiled as shared library
 *
 *        Version:  1.0
 *        Created:  05/06/2015
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Werner Van Geit,
 *        Company:  BBP
 *
 * =============================================================================
 */

#include <cfeature.h>
#include <efel.h>
#include <Python.h>

using namespace std;

extern cFeature* pFeature;


static PyObject* CppCoreInitialize(PyObject * self, PyObject * args) {

    char * depfilename, * outfilename;
    if (!PyArg_ParseTuple(args, "ss", &depfilename, &outfilename)) {
        return NULL;
    }

    Initialize(depfilename, outfilename);
    return Py_BuildValue("");
}


static vector<int> PyList_to_vectorint(PyObject *input) {
    vector<int> result_vector;
    int list_size;
    int index;
      
    list_size = PyList_Size(input);
    for (index = 0; index < list_size; index++) {
        result_vector.push_back(PyInt_AsLong(PyList_GetItem(input, index)));
    } 
    return result_vector;
}

void PyList_from_vectorint(vector<int> input, PyObject *output) {
    int vector_size;
    int index;

    vector_size = input.size();
    for (index = 0; index < vector_size; index++) {
        PyList_Append(output, Py_BuildValue("i", input[index]));
    } 
}

static vector<double> PyList_to_vectordouble(PyObject *input) {
    vector<double> result_vector;
    int list_size;
    int index;
      
    list_size = PyList_Size(input);
    for (index = 0; index < list_size; index++) {
        result_vector.push_back(PyFloat_AsDouble(PyList_GetItem(input, index)));
    } 
    return result_vector;
}

void PyList_from_vectordouble(vector<double> input, PyObject *output) {
    int vector_size;
    int index;

    vector_size = input.size();
    for (index = 0; index < vector_size; index++) {
        PyList_Append(output, Py_BuildValue("f", input[index]));
    } 
}

void PyList_from_vectorstring(vector<string> input, PyObject *output) {
    int vector_size;
    int index;

    vector_size = input.size();
    for (index = 0; index < vector_size; index++) {
        PyList_Append(output, Py_BuildValue("s", input[index].c_str()));
    } 
}

static PyObject* setfeatureint(PyObject * self, PyObject * args) {
  char * feature_name;
  PyObject * py_values;
  vector<int> values; 
  int return_value;
  if (!PyArg_ParseTuple(args, "sO!", &feature_name, &PyList_Type, &py_values)) {
          return NULL;
  }
  
  values = PyList_to_vectorint(py_values);
  return_value = pFeature->setFeatureInt(string(feature_name), values);
  
  return Py_BuildValue("i", return_value); 
}

static PyObject* getfeatureint(PyObject * self, PyObject * args) {
  char * feature_name;
  PyObject * py_values;
  vector<int> values; 
  int return_value;
  if (!PyArg_ParseTuple(args, "sO!", &feature_name, &PyList_Type, &py_values)) {
          return NULL;
  }
 
 
  return_value = pFeature->getFeatureInt(string(feature_name), values);
  PyList_from_vectorint(values, py_values);
  
  return Py_BuildValue("i", return_value); 
}

static PyObject* setfeaturedouble(PyObject * self, PyObject * args) {
  char * feature_name;
  PyObject * py_values;
  vector<double> values; 
  int return_value;
  if (!PyArg_ParseTuple(args, "sO!", &feature_name, &PyList_Type, &py_values)) {
          return NULL;
  }
  
  values = PyList_to_vectordouble(py_values);
  return_value = pFeature->setFeatureDouble(string(feature_name), values);
  
  return Py_BuildValue("f", return_value); 
}

static PyObject* getfeaturedouble(PyObject * self, PyObject * args) {
  char * feature_name;
  PyObject * py_values;
  vector<double> values; 
  int return_value;
  if (!PyArg_ParseTuple(args, "sO!", &feature_name, &PyList_Type, &py_values)) {
          return NULL;
  }
 
  return_value = pFeature->getFeatureDouble(string(feature_name), values);
  PyList_from_vectordouble(values, py_values);

  return Py_BuildValue("i", return_value); 
}

static PyObject* getFeatureNames(PyObject * self, PyObject * args) {
  vector<string> feature_names;
  PyObject * py_feature_names;
  if (!PyArg_ParseTuple(args, "O!", &PyList_Type, &py_feature_names)) {
          return NULL;
  }
 
  pFeature->get_feature_names(feature_names);
  PyList_from_vectorstring(feature_names, py_feature_names);

  return Py_BuildValue("");
}

static PyObject* getDistance(PyObject * self, PyObject * args) {
  char * feature_name;  
  double mean, std, distance;
  
  if (!PyArg_ParseTuple(args, "sdd", &feature_name, &mean, &std)) {
          return NULL;
  }
 
  distance = pFeature->getDistance(feature_name, mean, std);

  return Py_BuildValue("d", distance);
}

static PyObject* featuretype(PyObject * self, PyObject * args) {
  char * feature_name;
  string feature_type;
  
  if (!PyArg_ParseTuple(args, "s", &feature_name)) {
          return NULL;
  }
  
  feature_type = pFeature->featuretype(string(feature_name));

  return Py_BuildValue("s", feature_type.c_str());
}

static PyObject* getgerrorstr(PyObject * self, PyObject * args) {
    return Py_BuildValue("s", pFeature->getGError().c_str());
}


static PyMethodDef CppCoreMethods[] = {                                      
            {"Initialize",  CppCoreInitialize, METH_VARARGS,                      
                               "Initialise CppCore."},                               
            {"setFeatureInt",  setfeatureint, METH_VARARGS,                      
                               "Set a integer feature."},                               
            {"getFeatureInt",  getfeatureint, METH_VARARGS,                      
                               "Set a integer feature."},                               
            {"setFeatureDouble",  setfeaturedouble, METH_VARARGS,                      
                               "Set a double feature."},                               
            {"getFeatureDouble",  getfeaturedouble, METH_VARARGS,                      
                               "Get a double feature."},                               
            {"featuretype",  featuretype, METH_VARARGS,                      
                               "Get the type of a feature"},                               
            {"getgError",  getgerrorstr, METH_VARARGS,                      
                               "Get CppCore error string"},                               
            {"getFeatureNames",  getFeatureNames, METH_VARARGS,                      
                               "Get the names of all the available features"},                               
            {"getDistance",  getDistance, METH_VARARGS,                      
                    "Get the distance between a feature and experimental data"},                               
            {NULL, NULL, 0, NULL}        /* Sentinel */                              
};                                                                               
                                                                                 
PyMODINIT_FUNC initcppcore(void) {                                           
        (void) Py_InitModule("cppcore", CppCoreMethods);                     
}           
