#!/usr/bin/env python
##############################################################################
#
# diffpy.pdflive    by DANSE Diffraction group
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

import numpy as np

def _configPropertyRad(nm):
    '''helper function of options delegation, rad 2 degree'''
    rv = property(fget = lambda self: np.radians(getattr(self, nm)), 
                  fset = lambda self, val: setattr(self, nm, np.degrees(val)), 
                  fdel = lambda self: delattr(self, nm))
    return rv

def _configPropertyR(name):
    '''Create a property that forwards self.name to self.config.name.
    '''
    rv = property(fget = lambda self: getattr(self.config, name),
            doc='attribute forwarded to self.config, read-only')
    return rv

def _configPropertyRW(name):
    '''Create a property that forwards self.name to self.config.name.
    '''
    rv = property(fget = lambda self: getattr(self.config, nm), 
                  fset = lambda self, value: setattr(self.config, nm, value),
                  fdel = lambda self: delattr(self, nm),
                  doc='attribute forwarded to self.config, read/write')
    return rv
