#!/usr/bin/env python
##############################################################################
#
# diffpy.confutils  by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2013 Trustees of the Columbia University
#                   in the City of New York.  All rights reserved.
#
# File coded by:    Xiaohao Yang
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSENOTICE.txt for license information.
#
##############################################################################

"""Unit tests for diffpy.confutils
"""


def test():
    '''Execute all unit tests for the diffpy.confutils package.
    Return a unittest TestResult object.
    '''
    import unittest
    modulenames = '''
        diffpy.pdfgetx.confutils.test.testconfig
    '''.split()
    suite = unittest.TestSuite()
    loader = unittest.defaultTestLoader
    for mname in modulenames:
        exec ('import %s as mobj' % mname)
        suite.addTests(loader.loadTestsFromModule(mobj))
    runner = unittest.TextTestRunner()
    result = runner.run(suite)
    return result


# End of file
