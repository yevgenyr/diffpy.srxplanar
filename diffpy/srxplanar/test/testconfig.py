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


"""Unit tests for diffpy.confutils
"""

import unittest
import numpy as np
import ConfigParser
import re, os, sys, fnmatch
from functools import partial
import argparse

from diffpy.confutils.config import ConfigBase, initConfigClass
from diffpy.confutils.tools import _configPropertyRad, _configPropertyR, _configPropertyRW

##############################################################################

initConfigClass(ConfigBase)

class TestConfig(unittest.TestCase):
    testconfigdata = '''
[Experiment]
tifdirectory = ./
integrationspace = qspace
wavelength = 0.222
rotationd = 0.3223

[Beamline]
includepattern = *.tif, *.png
excludepattern = *.non.tif, *.test.tif
fliphorizontal = True

[Others]
regulartmatrixenable = True
maskedges = 2, 3, 4, 5, 50
'''
    
    testargsdata = '--tifdirectory ./ --integrationspace qspace --wavelength 0.222 --rotationd 0.3223 \
        --includepattern *.tif *.png --excludepattern  *.non.tif *.test.tif --fliphorizontal True \
        --regulartmatrixenable True --maskedges 2 3 4 5 50'.split()
    
    lines = testconfigdata.splitlines()
    lines[0] = '# title #'
    headershort = "\n".join(lines[:-1]+['']) + "\n"
    headerfull = "\n".join(lines+['']) + "\n"
    
    def setUp(self):
        self.config = ConfigBase()
        return
    
    def tearDown(self):
        self.config = None
        filelist = os.listdir('.')
        filelist = fnmatch.filter(filelist, '*.cfg')
        for f in filelist:
            os.remove(f)
        return
    
    def assertConfig(self, config=None, mode='full'):
        config = self.config if config==None else config
        self.assertEqual(config.tifdirectory, './')
        self.assertEqual(config.integrationspace, 'qspace')
        self.assertEqual(config.wavelength, 0.222)
        self.assertEqual(config.rotationd, 0.3223)
        self.assertEqual(config.includepattern, ['*.tif', '*.png'])
        self.assertEqual(config.excludepattern, ['*.non.tif', '*.test.tif'])
        self.assertEqual(config.fliphorizontal, True)
        self.assertEqual(config.regulartmatrixenable, True)
        if mode=='full':
            self.assertEqual(config.maskedges, [2, 3, 4, 5, 50])
        return
    
    def test_createConfig(self):
        '''test create config file
        '''
        self.config.updateConfig(args=self.testargsdata)
        self.assertConfig()
        
        self.config.writeConfig('testshort.cfg', mode='short')
        self.config.writeConfig('testfull.cfg', mode='full')
        self.assertTrue(os.path.exists('testshort.cfg'))
        self.assertTrue(os.path.exists('testfull.cfg'))
        
        config1 = ConfigBase(filename='testshort.cfg')
        self.assertConfig(config1, mode='short')
        config2 = ConfigBase(filename='testfull.cfg')
        self.assertConfig(config2, mode='full')
        return
    
    def test_configfile(self):
        '''test config file
        '''
        f = open('test.cfg', 'w')
        f.write(self.testconfigdata)
        f.close()
        
        self.config.updateConfig(filename='test.cfg')
        self.assertConfig()
        return
    
    def test_configfile_args(self):
        '''test read config from args
        '''
        f = open('test.cfg', 'w')
        f.write(self.testconfigdata)
        f.close()
        
        self.config.updateConfig(args=['-c', 'test.cfg'])
        self.assertConfig()
        return

    def test_args(self):
        '''test args
        '''
        self.config.updateConfig(args=self.testargsdata)
        self.assertConfig()
        return

    def test_kwargs(self):
        '''test kwargs
        '''
        self.config.updateConfig(
            tifdirectory = './',
            integrationspace = 'qspace',
            wavelength = 0.222,
            rotationd = 0.3223,
            includepattern = ['*.tif', '*.png'],
            excludepattern = ['*.non.tif', '*.test.tif'],
            fliphorizontal = True,
            regulartmatrixenable = True,
            maskedges = [2, 3, 4, 5, 50])
        self.assertConfig()
        return
    
    def test_headers(self):
        '''test header
        '''
        self.config.updateConfig(args=self.testargsdata)
        self.assertConfig()
        headershort = self.config.getHeader('title', 'short')
        self.assertEqual(headershort, self.headershort)
        headerfull = self.config.getHeader('title', 'full')
        self.assertEqual(headerfull, self.headerfull)
        return
    

# End of class TestRoutines

if __name__ == '__main__':
    unittest.main()

# End of file
