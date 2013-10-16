diffpy.SrXplanar

diffpy.SrXplanar package provides 2D diffraction image integration using
non splitting pixel algorithm. And it can estimate and propagate statistic
uncertainty of raw counts and integrated intensity. 

To learn more about diffpy.Structure library, see the examples directory
included in this distribution or the API documentation at

    http://docs.danse.us/diffraction/diffpy.srxplanar


REQUIREMENTS

The diffpy.Structure requires Python 2.6 and the following external software:

    setuptools	>=0.61(https://pypi.python.org/pypi/setuptools)
    numpy       >=1.60(http://www.numpy.org/)
    scipy		>=1.10(www.scipy.org/)
    FabIO		>=0.80(http://sourceforge.net/projects/fable/files/fabio/)

If your python version < 2.7 (these two packages are included in 2.7 but not in 2.6)
	ordereddict	(https://pypi.python.org/pypi/ordereddict)
	argparse	(https://pypi.python.org/pypi/argparse)

We recommand to install EPD (Enthought Python Distribution)from 
http://www.enthought.com/ which include most of them (except FabIO) by default.

INSTALLATION

To install the diffpy.SrXplanar package:

    python setup.py install

By default the files are installed in the system directories, which are
usually only writeable by the root.  See the usage info 
"./setup.py install --help" for options to install as a normal user under
different location.  Note that installation to non-standard directories may
require adjustments to the PATH and PYTHONPATH environment variables.

The Python setuptools library provides "easy_install" script, which can
be used to update diffpy.Structure installation or even to perform a new
install without explicit need to download and unzip the code:

    easy_install -U diffpy.srxplanar

This checks the package repository at http://www.diffpy.org/packages/
for any newer releases of diffpy.Structure and if present, it updates the
installation.  The easy_install can be also used to get in sync with the
latest development sources in the subversion repository:

    easy_install -U \
        svn://svn@danse.us/diffraction/diffraction/diffpy.srxplanar/trunk


CONTACTS

For more information on diffpy.Structure please visit the project web-page:

    http://www.diffpy.org/

or email Prof. Simon Billinge at sb2896@columbia.edu

