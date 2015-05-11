Requirements                                                                     
============                                                                     
                                                                                 
* [Python 2.7+](https://www.python.org/download/releases/2.7/) (not compatible yet with Python 3)
* [Pip](https://pip.pypa.io) (installed by default in newer versions of Python)  
* [Numpy](http://www.numpy.org) (will be installed automatically by pip)         
* The instruction below are written assuming you have access to a command shell  
on Linux / UNIX / MacOSX / Cygwin                                                
                                                                                 
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
