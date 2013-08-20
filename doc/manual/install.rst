Installation
========================================================================


Software requirements
------------------------------------------------------------------------

SrXplanar has been written in Python programming language, therefore
to use the software, you must have Python 2.6 or
Python 2.7 installed. [#fpy3]_
In addition, the following third-party Python libraries are
also required:

* setuptools - >=0.61, tools for installing Python packages
* numpy - >=1.60, library for scientific computing with Python
* scipy - >=1.10, Scientific Library for Python
* FabIO - >=0.80, 'Image IO for fable <http://sourceforge.net/projects/fable/files/fabio/>'_

If your python version < 2.7 (these two packages are included in 2.7 but not in 2.6)
* ordereddict - https://pypi.python.org/pypi/ordereddict
* argparse - https://pypi.python.org/pypi/argparse

Standard Python releases can be obtained from
http://www.python.org/download/.
The third-party libraries can be found at the
`Python Package Index <http://pypi.python.org/pypi>`_
or using any Internet search engine.

Another, more convenient option is to obtain science-oriented Python
distributions, such as `PythonXY <https://code.google.com/p/pythonxy/>`_
or `Enthought Canopy <http://www.enthought.com/>`_.  These distributions
already include all the necessary libraries, so the required Python
software can be all installed in one step.

On Linux operating systems the third-party libraries are usually
included in a system software package repository.  For example on an
Ubuntu Linux computer the software dependencies can be all installed
with a single shell command ::

  sudo apt-get install \
    python-distribute python-numpy python-scipy

This may be, of course, just as well accomplished using the GUI
driven Synaptic package manager.  Other Linux
distributions may use different software management tools,
but the names of the necessary packages should be very similar
to those above.

On Windows operating system, it may be necessary to add the
``C:\Python27`` directory and the scripts directory
``C:\Python27\Scripts`` to the system :envvar:`!PATH`.
Some Python distributions may already do that as a part of their
installation process.  The easiest way to check is to start the
:program:`Command Prompt`, type there ``python`` and see if this
starts the Python interpreter.


PDFgetX3 installation
------------------------------------------------------------------------

FIXME: update the URL

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

    easy_install -U diffpy.SrXplanar

This checks the package repository at http://www.diffpy.org/packages/
for any newer releases of diffpy.Structure and if present, it updates the
installation.  The easy_install can be also used to get in sync with the
latest development sources in the subversion repository:

    easy_install -U \
        svn://svn@danse.us/diffraction/diffraction/diffpy.srxplanar/trunk

.. rubric:: Footnotes

.. [#fpy3] Python 3 is not supported yet.
