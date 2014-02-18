Installation
========================================================================

.. index:: Software requirements

Software requirements
------------------------------------------------------------------------

The diffpy.srxplanar requires Python 2.6 or 2.7 and the following software:

* setuptools	>=0.61(https://pypi.python.org/pypi/setuptools)
* numpy       >=1.60(http://www.numpy.org/)
* scipy		>=1.10(www.scipy.org/)
* FabIO		>=0.80(http://sourceforge.net/projects/fable/files/fabio/)

If your python version < 2.7 (these two packages are included in 2.7 but not in 2.6)
* ordereddict	(https://pypi.python.org/pypi/ordereddict)
* argparse	(https://pypi.python.org/pypi/argparse)

On Ubuntu Linux the part of required software can be easily installed using
the system package manager:

    sudo aptitude install \
        python-setuptools python-numpy python-scipy
        
For Mac OS X machine with the MacPorts package manager one could do

    sudo port install \
        python27 py27-setuptools py27-numpy py27-scipy

When installing with MacPorts, make sure the MacPorts bin directory is the
first in the system PATH and that python27 is selected as the default
Python version in MacPorts:

    sudo port select --set python python27
    
For other Linux distributions use their respective package manager; note
the packages may have slightly different names. diffpy.srxplanar should work
on other Unix-like operating systems as well.  Please, search the
web for instructions how to install external dependencies on your particular
system.

For other packages, please go to the webpage list above to download and install. 

.. index:: SrXplanar installation

SrXplanar installation
------------------------------------------------------------------------

To install the diffpy.srxplanar package:

    python setup.py install

By default the files are installed in the system directories, which are
usually only writeable by the root.  See the usage info 
"./setup.py install --help" for options to install as a normal user under
different location.  Note that installation to non-standard directories may
require adjustments to the PATH and PYTHONPATH environment variables.

DEVELOPMENT and CONTRIBUTION

diffpy.srxplanar is an open-source software developed at the Columbia University
The diffpy.srxplanar sources are hosted at

    https://github.com/diffpy/diffpy.srxplanar

Feel free to fork the project and contribute.  To install diffpy.srxplanar
in a development mode, where the sources are directly used by Python
rather than copied to a system directory, use

    python setup.py develop --user

.. rubric:: Footnotes

.. [#fpy3] Python 3 is not supported yet.
