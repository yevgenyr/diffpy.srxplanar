#!/usr/bin/env python
##############################################################################
#
# diffpy.srxplanar  by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2010 Trustees of the Columbia University
#                   in the City of New York.  All rights reserved.
#
# File coded by:    Xiaohao Yang
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSENOTICE.txt for license information.
#
##############################################################################

# package version
from diffpy.srxplanar.version import __version__

# some convenience imports
from diffpy.srxplanar.srxplanar import SrXplanar
from diffpy.srxplanar.srxplanarconfig import SrXplanarConfig

# unit tests
def test():
    '''
    Execute all unit tests for the diffpy.pdfgetx package.
    Return a unittest TestResult object.
    '''
    from diffpy.srxplanar.tests import test
    return test()


# End of file
