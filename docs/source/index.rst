.. eFEL documentation master file, created by
   sphinx-quickstart on Mon May 11 14:40:15 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Electrophys Feature Extraction Library
======================================
The Electrophys Feature Extract Library (eFEL) allows neuroscientists to 
automatically extract eFeatures from time series data recorded from neurons 
(both in vitro and in silico). Examples are the action potential width and 
amplitude in voltage traces recorded during whole-cell patch clamp experiments. 
The user of the library provides a set of traces and selects the eFeatures to be 
calculated. The library will then extract the requested eFeatures and return the 
values to the user.

The core of the library is written in C++, and a Python wrapper is included. 
At the moment we provide a way to automatically compile and install the library 
as a Python module. Soon instructions will be added on how to link C++ code 
directly with the eFEL.

The source code of the eFEL is located on github: 
`BlueBrain/eFEL <https://github.com/BlueBrain/eFEL>`_

.. toctree::
   :maxdepth: 3

   installation
   examples
   eFeatures
   api

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

