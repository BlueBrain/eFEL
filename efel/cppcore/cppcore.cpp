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

#include <Python.h>

#include <cstddef>
#include <cfeature.h>
#include <efel.h>

extern cFeature* pFeature;

static PyObject* CppCoreInitialize(PyObject* self, PyObject* args) {

  char* depfilename, *outfilename;
  if (!PyArg_ParseTuple(args, "ss", &depfilename, &outfilename)) {
    return NULL;
  }

  Initialize(depfilename, outfilename);
  return Py_BuildValue("");
}

static vector<int> PyList_to_vectorint(PyObject* input) {
  vector<int> result_vector;
  int list_size;
  int index;

  list_size = PyList_Size(input);
  for (index = 0; index < list_size; index++) {
    result_vector.push_back(PyLong_AsLong(PyList_GetItem(input, index)));
  }
  return result_vector;
}

static void PyList_from_vectorint(vector<int> input, PyObject* output) {
  size_t vector_size = input.size();

  for (size_t index = 0; index < vector_size; index++) {
    PyObject *obj = Py_BuildValue("i", input[index]);
    PyList_Append(output, obj);
    Py_DECREF(obj);
  }
}

static vector<double> PyList_to_vectordouble(PyObject* input) {
  vector<double> result_vector;
  int list_size;
  int index;

  list_size = PyList_Size(input);
  for (index = 0; index < list_size; index++) {
    result_vector.push_back(PyFloat_AsDouble(PyList_GetItem(input, index)));
  }
  return result_vector;
}

static void PyList_from_vectordouble(vector<double> input, PyObject* output) {
  size_t vector_size = input.size();

  for (size_t index = 0; index < vector_size; index++) {
    PyObject *obj = Py_BuildValue("f", input[index]);
    PyList_Append(output, obj);
    Py_DECREF(obj);
  }
}

static void PyList_from_vectorstring(vector<string> input, PyObject* output) {
  size_t vector_size = input.size();

  for (size_t index = 0; index < vector_size; index++) {
    PyObject *obj = Py_BuildValue("s", input[index].c_str());
    PyList_Append(output, obj);
    Py_DECREF(obj);
  }
}

static PyObject*
_getfeature(PyObject* self, PyObject* args, const string &input_type) {
  char* feature_name;
  PyObject* py_values;

  int return_value;
  if (!PyArg_ParseTuple(args, "sO!", &feature_name, &PyList_Type, &py_values)) {
    PyErr_SetString(PyExc_TypeError, "Unexpected argument type provided.");
    return NULL;
  }

  try {
    string feature_type = pFeature->featuretype(string(feature_name));

    if (!input_type.empty() && feature_type != input_type){  // when types do not match
      PyErr_SetString(PyExc_TypeError, "Feature type does not match");
      return NULL;
    }

    if (feature_type == "int") {
      vector<int> values;
      return_value = pFeature->getFeature<int>(string(feature_name), values);
      PyList_from_vectorint(values, py_values);
      return Py_BuildValue("i", return_value);
    } else { // double
      vector<double> values;
      return_value = pFeature->getFeature<double>(string(feature_name), values);
      PyList_from_vectordouble(values, py_values);
      return Py_BuildValue("i", return_value);
    }
  }
  catch(EfelAssertionError& e) {  // more specialised exception
    PyErr_SetString(PyExc_AssertionError, e.what());
    return NULL;
  }
  catch(const std::runtime_error& e) {
    PyErr_SetString(PyExc_RuntimeError, e.what());
    return NULL;
  }
  catch(...) {
    PyErr_SetString(PyExc_RuntimeError, "Unknown error");
    return NULL;
  }
}


static PyObject*
_getmapdata(PyObject* self, PyObject* args, const string &type) {
  char* data_name;
  PyObject* py_values = PyList_New(0);

  if (!PyArg_ParseTuple(args, "s", &data_name)) {
    return NULL;
  }

  if (type == "int") {
    vector<int> return_values;
    return_values = pFeature->getmapIntData(string(data_name));
    PyList_from_vectorint(return_values, py_values);
  } else if (type == "double") {
    vector<double> return_values;
    return_values = pFeature->getmapDoubleData(string(data_name));
    PyList_from_vectordouble(return_values, py_values);
  } else {
    PyErr_SetString(PyExc_TypeError, "Unknown data name");
    return NULL;
  }

  return py_values;
}

static PyObject* getfeature(PyObject* self, PyObject* args) {
  const string empty("");
  return _getfeature(self, args, empty);
}

static PyObject* setfeatureint(PyObject* self, PyObject* args) {
  char* feature_name;
  PyObject* py_values;
  vector<int> values;
  int return_value;
  if (!PyArg_ParseTuple(args, "sO!", &feature_name, &PyList_Type, &py_values)) {
    return NULL;
  }

  values = PyList_to_vectorint(py_values);
  return_value = pFeature->setFeatureInt(string(feature_name), values);

  return Py_BuildValue("i", return_value);
}

static PyObject* getfeatureint(PyObject* self, PyObject* args) {
  const string type("int");
  return _getfeature(self, args, type);
}
static PyObject* getmapintdata(PyObject* self, PyObject* args) {
  const string type("int");
  return _getmapdata(self, args, type);
}
static PyObject* getmapdoubledata(PyObject* self, PyObject* args) {
  const string type("double");
  return _getmapdata(self, args, type);
}

static PyObject* setfeaturedouble(PyObject* self, PyObject* args) {
  char* feature_name;
  PyObject* py_values;
  vector<double> values;
  int return_value;
  if (!PyArg_ParseTuple(args, "sO!", &feature_name, &PyList_Type, &py_values)) {
    return NULL;
  }

  values = PyList_to_vectordouble(py_values);
  return_value = pFeature->setFeatureDouble(string(feature_name), values);

  return Py_BuildValue("f", return_value);
}

static PyObject* getfeaturedouble(PyObject* self, PyObject* args) {
  const string type ("double");
  return _getfeature(self, args, type);
}

static PyObject* setfeaturestring(PyObject* self, PyObject* args){
  char* feature_name;
  char* py_value;

  int return_value;
  if (!PyArg_ParseTuple(args, "ss", &feature_name, &py_value)) {
    return NULL;
  }

  return_value = pFeature->setFeatureString(string(feature_name), py_value);

  return Py_BuildValue("i", return_value);

}

static PyObject* getFeatureNames(PyObject* self, PyObject* args) {
  vector<string> feature_names;
  PyObject* py_feature_names;
  if (!PyArg_ParseTuple(args, "O!", &PyList_Type, &py_feature_names)) {
    return NULL;
  }

  pFeature->get_feature_names(feature_names);
  PyList_from_vectorstring(feature_names, py_feature_names);

  return Py_BuildValue("");
}

static PyObject* getDistance_wrapper(PyObject* self,
                                     PyObject* args,
                                     PyObject* kwds) {
  char* feature_name;
  double mean, std, distance, error_dist=250;
  int trace_check = 1;

  const char *kwlist[] = {"feature_name", "mean", "std", "trace_check", 
      "error_dist", NULL};

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "sdd|id",
                                   const_cast<char**>(kwlist),
                                   &feature_name, &mean, &std,
                                   &trace_check, &error_dist)) {
    return NULL;
  }

  distance = pFeature->getDistance(feature_name, mean, std,
                                   trace_check, error_dist);

  return Py_BuildValue("d", distance);
}

static PyObject* featuretype(PyObject* self, PyObject* args) {
  char* feature_name;
  string feature_type;

  if (!PyArg_ParseTuple(args, "s", &feature_name)) {
    return NULL;
  }

  feature_type = pFeature->featuretype(string(feature_name));

  return Py_BuildValue("s", feature_type.c_str());
}

static PyObject* getgerrorstr(PyObject* self, PyObject* args) {
  return Py_BuildValue("s", pFeature->getGError().c_str());
}

static PyMethodDef CppCoreMethods[] = {
    {"Initialize", CppCoreInitialize, METH_VARARGS,
      "Initialise CppCore."},

    {"getFeature", getfeature, METH_VARARGS,
      "Get a values associated with a feature. Takes a list() to be filled."},
    {"getFeatureInt", getfeatureint, METH_VARARGS,
      "Get a integer feature."},
    {"getFeatureDouble", getfeaturedouble, METH_VARARGS,
      "Get a double feature."},
    {"getMapIntData", getmapintdata, METH_VARARGS,
      "Get a int data."},
    {"getMapDoubleData", getmapdoubledata, METH_VARARGS,
      "Get a double data."},

    {"setFeatureInt", setfeatureint, METH_VARARGS,
      "Set a integer feature."},
    {"setFeatureDouble", setfeaturedouble, METH_VARARGS,
      "Set a double feature."},
    {"setFeatureString", setfeaturestring, METH_VARARGS,
      "Set a string feature."},

    {"featuretype", featuretype, METH_VARARGS,
      "Get the type of a feature"},
    {"getgError", getgerrorstr, METH_VARARGS,
      "Get CppCore error string"},
    {"getFeatureNames", getFeatureNames, METH_VARARGS,
      "Get the names of all the available features"},

    {"getDistance", (PyCFunction)getDistance_wrapper, METH_VARARGS|METH_KEYWORDS,
      "Get the distance between a feature and experimental data"},
    {NULL, NULL, 0, NULL} /* Sentinel */
};

struct module_state {
  PyObject* error;
};

#define GETSTATE(m) ((struct module_state*)PyModule_GetState(m))
static int cppcore_traverse(PyObject* m, visitproc visit, void* arg) {
  Py_VISIT(GETSTATE(m)->error);
  return 0;
}

static int cppcore_clear(PyObject* m) {
  Py_CLEAR(GETSTATE(m)->error);
  return 0;
}

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,       "cppcore",      NULL,
    sizeof(struct module_state), CppCoreMethods, NULL,
    cppcore_traverse,            cppcore_clear,  NULL};

extern "C" PyObject* PyInit_cppcore(void) {
  PyObject* module = PyModule_Create(&moduledef);
  return module;
}
