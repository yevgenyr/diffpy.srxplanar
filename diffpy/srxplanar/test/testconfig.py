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


"""Unit tests for diffpy.srxplanar
"""

import unittest
import numpy as np
import ConfigParser
import re, os, sys, fnmatch
from functools import partial
import argparse

from diffpy.srxplanar.srxplanarconfig import SrXplanarConfig

##############################################################################

class TestConfig(unittest.TestCase):
    testconfigdatafull = '''
[Control]
summation = True

[Experiment]
fit2dconfig = 
tifdirectory = ./
savedirectory = ./
backgroundfile = bkg.tif
addmask = edgemask, deadpixel
integrationspace = qspace
wavelength = 0.222
xbeamcenter = 2000.0
ybeamcenter = 2000.0
distance = 100.0
rotationd = 2.0
tiltd = 2.0
tthstepd = 0.01
qstep = 0.01

[Beamline]
includepattern = *.tif
excludepattern = *.dark.tif, *.test.tif
fliphorizontal = True
flipvertical = False
xdimension = 3045
ydimension = 3045
xpixelsize = 0.1
ypixelsize = 0.1

[Others]
uncertaintyenable = True
sacorrectionenable = False
polcorrectionenable = False
polcorrectf = 0.50
selfcorrenable = False
gsasoutput = esd
filenameplus = test
maskedges = 20, 20, 20, 20, 200
'''

    testconfigdatashort = '''
[Control]
summation = True

[Experiment]
fit2dconfig = 
tifdirectory = ./
savedirectory = ./
backgroundfile = bkg.tif
addmask = edgemask, deadpixel
integrationspace = qspace
wavelength = 0.222
xbeamcenter = 2000.0
ybeamcenter = 2000.0
distance = 100.0
rotationd = 2.0
tiltd = 2.0
tthstepd = 0.01
qstep = 0.01

[Beamline]
xdimension = 3045
ydimension = 3045
xpixelsize = 0.1
ypixelsize = 0.1

[Others]
uncertaintyenable = True
gsasoutput = esd
filenameplus = test
'''
    testargsdata = '--summation True --tifdirectory ./ --savedirectory ./ --backgroundfile bkg.tif --addmask edgemask deadpixel \
        --integrationspace qspace --wavelength 0.222 --xbeamcenter 2000.0 --ybeamcenter 2000.0 \
        --distance 100.0 --rotationd 2.0 --tiltd 2.0 --tthstepd 0.01 --qstep 0.01 --includepattern *.tif\
        --excludepattern *.dark.tif *.test.tif --fliphorizontal True --flipvertical False \
        --xdimension 3045 --ydimension 3045 --xpixelsize 0.1 --ypixelsize 0.1 --uncertaintyenable True\
        --sacorrectionenable False --polcorrectionenable False --polcorrectf 0.5 --selfcorrenable False\
        --gsasoutput esd --filenameplus test --maskedges 20 20 20 20 200'.split()
    
    lines = testconfigdatashort.splitlines()
    lines[0] = '# SrXplanar configration #'
    headershort = "\n".join(lines+['']) + "\n"
    lines = testconfigdatafull.splitlines()
    lines[0] = '# SrXplanar configration #'
    headerfull = "\n".join(lines+['']) + "\n"
    
    def setUp(self):
        self.config = SrXplanarConfig()
        return
    
    def tearDown(self):
        self.config = None
        filelist = os.listdir('.')
        for extname in ('*.cfg', '*.tif', '*.chi', '*.fq', '*.gr'):
            filelist1 = fnmatch.filter(filelist, extname)
            for f in filelist1:
                os.remove(f)
        return
    
    def assertConfig(self, config=None, mode='full'):
        config = self.config if config==None else config
        
        self.assertEqual(config.summation, True)
        self.assertEqual(config.tifdirectory, './')
        self.assertEqual(config.savedirectory, './')
        self.assertEqual(config.backgroundfile, 'bkg.tif')
        self.assertEqual(config.addmask, ['edgemask', 'deadpixel'])
        self.assertEqual(config.integrationspace, 'qspace')
        self.assertEqual(config.wavelength, 0.222)
        self.assertEqual(config.xbeamcenter, 2000.0)
        self.assertEqual(config.ybeamcenter, 2000.0)
        self.assertEqual(config.distance, 100.0)
        self.assertEqual(config.rotationd, 2.0)
        self.assertEqual(config.tiltd, 2.0)
        self.assertEqual(config.tthstepd, 0.01)
        self.assertEqual(config.qstep, 0.01)
        self.assertEqual(config.xdimension, 3045)
        self.assertEqual(config.ydimension, 3045)
        self.assertEqual(config.xpixelsize, 0.1)
        self.assertEqual(config.ypixelsize, 0.1)
        self.assertEqual(config.uncertaintyenable, True)
        self.assertEqual(config.gsasoutput, 'esd')
        self.assertEqual(config.filenameplus, 'test')
        if mode=='full':
            self.assertEqual(config.includepattern, ['*.tif'])
            self.assertEqual(config.excludepattern, ['*.dark.tif', '*.test.tif'])
            self.assertEqual(config.fliphorizontal, True)
            self.assertEqual(config.flipvertical, False)
            self.assertEqual(config.sacorrectionenable, False)
            self.assertEqual(config.polcorrectionenable, False)
            self.assertEqual(config.polcorrectf, 0.50)
            self.assertEqual(config.selfcorrenable, False)
            self.assertEqual(config.maskedges, [20, 20, 20, 20, 200])
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
        
        config1 = SrXplanarConfig(filename='testshort.cfg')
        self.assertConfig(config1, mode='short')
        config2 = SrXplanarConfig(filename='testfull.cfg')
        self.assertConfig(config2, mode='full')
        return
    
    def test_configfile(self):
        '''test config file
        '''
        #full
        f = open('test.cfg', 'w')
        f.write(self.testconfigdatafull)
        f.close()
        self.config.updateConfig(filename='test.cfg')
        self.assertConfig()
        #short
        f = open('test.cfg', 'w')
        f.write(self.testconfigdatashort)
        f.close()
        config = SrXplanarConfig(filename='test.cfg')
        self.assertConfig(config, mode='short')
        return
    
    def test_configfile_args(self):
        '''test read config from args
        '''
        f = open('test.cfg', 'w')
        f.write(self.testconfigdatafull)
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
            summation = True,
            tifdirectory = './',
            savedirectory = './',
            backgroundfile = 'bkg.tif',
            addmask = ['edgemask', 'deadpixel'],
            integrationspace = 'qspace',
            wavelength = 0.222,
            xbeamcenter = 2000.0,
            ybeamcenter = 2000.0,
            distance = 100.0,
            rotationd = 2.0,
            tiltd = 2.0,
            tthstepd = 0.01,
            qstep = 0.01,
            includepattern = ['*.tif'],
            excludepattern = ['*.dark.tif', '*.test.tif'],
            fliphorizontal = True,
            flipvertical = False,
            xdimension = 3045,
            ydimension = 3045,
            xpixelsize = 0.1,
            ypixelsize = 0.1,
            uncertaintyenable = True,
            sacorrectionenable = False,
            polcorrectionenable = False,
            polcorrectf = 0.50,
            selfcorrenable = False,
            gsasoutput = 'esd',
            filenameplus = 'test',
            maskedges = [20, 20, 20, 20, 200])
        self.assertConfig()
        return
    
    def tttest_headers(self):
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
