[![Build Status](https://travis-ci.org/BlueBrain/eFEL.png?branch=master)](https://travis-ci.org/BlueBrain/eFEL)

Introduction
============

The Electrophys Feature Extract Library (eFEL) allows neuroscientists
to automatically extract features from time series data recorded from neurons 
(both in vitro and in silico). 
Examples are the action potential width and amplitude in voltage traces recorded
during whole-cell patch clamp experiments.
The user of the library provides a set of traces and selects the features to
be calculated. The library will then extract the requested features and return
the values to the user.

The core of the library is written in C++, and a Python wrapper is included.
At the moment we provide a way to automatically compile and install the library
as a Python module. Soon instructions will be added on how to link C++ code 
directly with the eFEL.

Requirements
============

* [Python 2.7+](https://www.python.org/download/releases/2.7/) (not compatible with Python 3)
* [Pip](https://pip.pypa.io) (installed by default in newer version of Python)
* [Numpy](http://www.numpy.org) (will be installed automatically by pip)

Installation
============

The easiest way to install eFEL is to use [pip](https://pip.pypa.io)

```bash
pip install git+git://github.com/BlueBrain/eFEL
```
In case you don't have administrator access this command might fail with a 
permission error. In that case you could install eFEL in your home directory

```bash
pip install --user git+git://github.com/BlueBrain/eFEL
```

Or you could use a [python virtual environment](https://virtualenv.pypa.io)

```bash
virtualenv pythonenv
. ./pythonenv/bin/activate
pip install git+git://github.com/BlueBrain/eFEL
```

Quick Start
===========

First you need to import the module

```python
import efel
```

To get a list with all the available feature names

``python
efel.getFeatureNames()
```
