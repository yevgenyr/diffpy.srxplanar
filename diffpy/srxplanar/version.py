#!/usr/bin/env python
##############################################################################
#
# diffpy.srxplanar  by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2012 Trustees of the Columbia University
#                   in the City of New York.  All rights reserved.
#
# File coded by:    Xiaohao Yang
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
##############################################################################

"""Definition of __version__ and __date__ for the srsig2d package.
"""

# obtain version information
from pkg_resources import get_distribution
_pkgname = __name__.rsplit('.', 1)[0]
__version__ = get_distribution(_pkgname).version

# we assume that tag_date was used and __version__ ends in YYYYMMDD
__date__ = __version__[-8:-4] + '-' + \
           __version__[-4:-2] + '-' + __version__[-2:]

# End of file
