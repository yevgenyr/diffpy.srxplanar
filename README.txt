INSTALLATION

To install the diffpy.srxplanar package:

    python setup.py install

By default the files are installed in the system directories, which are
usually only writeable by the root.  See the usage info "./setup.py install --help" for options to install as a normal user under a different location.  Note that installation to non-standard directories you may require adjustments to the PATH and PYTHONPATH environment variables.

DEPENDENCIES


diffpy.srxplanar requires the third-party Python libraries, 
numpy >1.60,
scipy >1.10, 
FabIO(http://sourceforge.net/projects/fable/files/fabio/)

if python version < 2.7 (these two packages are included in 2.7 but not in 2.6)
ordereddict  
argparse

We recommand to install EPD (Enthought Python Distribution)from http://www.enthought.com/ which include most of them by default.

------------------------------------------------------------------------------

For more information email Prof. Simon Billinge at sb2896@columbia.edu
