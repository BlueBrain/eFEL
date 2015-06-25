Installation                                                                     
============                                                                     
            
Requirements
------------
* `Python 2.6+ <https://www.python.org/download/releases/2.7/>`_ (not compatible yet with Python 3)
* `Pip <https://pip.pypa.io>`_ (installed by default in newer versions of Python)  
* `Numpy <http://www.numpy.org>`_ (will be installed automatically by pip)         
* The instruction below are written assuming you have access to a command shell on Linux / UNIX / MacOSX / Cygwin                                                
                                                                                 
Installation using pip
----------------------                                                                                

The easiest way to install eFEL is to use `pip <https://pip.pypa.io>`_::

    pip install git+git://github.com/BlueBrain/eFEL                                  
                                                                              
In case you don't have administrator access this command might fail with a       
permission error. In that case you could install eFEL in your home directory::

    pip install --user git+git://github.com/BlueBrain/eFEL                           
                                                                                 
Or you could use a `python virtual environment <https://virtualenv.pypa.io>`_::

    virtualenv pythonenv                                                             
    . ./pythonenv/bin/activate                                                       
    pip install git+git://github.com/BlueBrain/eFEL

Installing the C++ standalone library
-------------------------------------

If your system doesn't have it, install `CMake <http://www.cmake.org/>`_.

Make a new build directory::

    mkdir build_cmake

Configure the build, replace YOURINSTALLDIR with the directory in which you want
to install the efel library (e.g. /usr/local)::

    cd build_cmake
    cmake .. -DCMAKE_INSTALL_PREFIX=YOURINSTALLDIR

Run the compilation and installation::

    make install

This will have installed a static and shared library as::
    
    YOURINSTALLDIR/lib/libefel.a
    YOURINSTALLDIR/lib/libefel.so
