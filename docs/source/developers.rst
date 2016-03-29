=================
Developer's Guide
=================

.. contents::

Requirements
============
As a developer you will need some extra requirements

* To get the latest source code: `Git <https://git-scm.com/>`_
* To run the tests: `Nose <https://nose.readthedocs.org/en/latest/>`_
* To build the documentation: `Sphinx <http://sphinx-doc.org/>`_, and pdflatex
  (e.g. from `Mactex <https://tug.org/mactex/>`_)

Forking and cloning the git repository
======================================
To make changes to the eFEL, one first needs to fork the eFEL::

    https://help.github.com/articles/fork-a-repo/

Then one creates a local clone of the git repository on your computer::

    git clone https://github.com/yourgithubusername/eFEL.git

After changes are made, they should be pushed back to your github account.
Then a pull request can be created::

    https://help.github.com/articles/using-pull-requests/

Makefile
========
To simplify certain tasks for developers, a Makefile is provided in the root of
the eFEL project. This Makefile has the following targets

* **install**: installs the eFEL using pip from the working directory
* **test**: run the installation and all the Nose tests
* **doc**: build the sphinx and latex documentation
* **clean**: clean up the build directories
* **pypi**: run test target and upload to pypi
* **push**: clean the build, update the version from the git hash, install eFEL,
  run the tests, build the doc, and push the documentation and source to github

Adding a new eFeature
=====================
Adding a new eFeature requires several steps.

Picking a name
--------------
Try to be specific in the name of the eFeature, because in the future you or
somebody else might want to develop an eFeature with slightly different
behavior. Don't be afraid to use long names, e.g. 'min_voltage_between_spikes'
is perfectly ok.

Creating a branch
-----------------
Create a git branch with the name of the new eFeature::

    git checkout -b your_efeaturename

Implementation
--------------
All the eFeatures in the eFEL are coded in C++. Thanks to an
`eFeatures dependency settings file <https://github.com/BlueBrain/eFEL/blob/
master/efel/DependencyV5.txt>`_,
several implementation of the same eFeature name can coexist. E.g.
`this <https://github.com/BlueBrain/eFEL/blob/master/efel/cppcore/LibV5.cpp>`_
is the file with the implementations of all 'V5' features.
You can implement the new eFeature by extending one of the current LibV* files,
or by creating your own.
You might want to consider starting the implementation by writing a test for
the eFeature (see below for instruction on how to do that).

Updating relevant files
-----------------------
Apart from the implementation in the LibV*.cpp file, other files have to be
changed to accomodate the new eFeature

* efel/cppcore/LibV5.h: Declare your feature
* efel/DependencyV5.txt: Add your eFeature and its dependencies to this file
* efel/cppcore/FillFptrTable.cpp: Add a reference to the eFeature in the
  relevant table
* efel/cppcore/cfeature.cpp: Add the type of the eFeature
* AUTHORS.txt: If your name isn't there yet, add yourself to the authors list

You can confirm everything compiles correctly by executing::

    make test

Adding a test
-------------
Most eFeatures are fairly easy to implement in Python, so it is advised to first
write a Python implementation of your eFeature, and to add it Nose tests.
Then, while you are implementing the code in C++ you can easily compare the
results to the Nose test.

The Nose tests of the individual eFeatures are
`here <https://github.com/BlueBrain/eFEL/blob/master/efel/tests/
test_basic.py>`_
.Just add your own test by defining a new function 'test_yourfeature()'.

Some test data is available
`here <https://github.com/BlueBrain/eFEL/tree/master/efel/tests/
testdata/basic>`_
, but you can of course add your own traces.

The easiest way to run the tests is by executing::

    make test

Add documentation
-----------------
Add the documentation of the new eFeature to this file:

https://github.com/BlueBrain/eFEL/blob/master/docs/source/eFeatures.rst

Please provide some pseudo-Python code for the eFeature.

The documentation can be built by::

    make doc

It can be viewed by opening::

    docs/build/html/index.html

To build the documentation, pdflatex has to be present on the system. On a Mac
this can be installed using `Mactex <https://tug.org/mactex/>`_. On Ubuntu one
can use::

    sudo apt-get install texlive-latex-base texlive-latex-extra xzdec
    tlmgr install helvetic

Pull request
-------------
When all the above steps were succesfull, you can push the
new eFeature branch to your github repository::

    git commit -a
    git push origin your_efeaturename

Finally create a pull request:

https://help.github.com/articles/using-pull-requests/
