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
import re, os, sys, fnmatch, itertools
from functools import partial
import argparse

from diffpy.srxplanar.srxplanarconfig import SrXplanarConfig
from diffpy.srxplanar.srxplanar import SrXplanar

##############################################################################

class TestConfig(unittest.TestCase):
    testargsdata = '--summation True --tifdirectory ./ --savedirectory ./ --backgroundfile bkg.tif --addmask edgemask deadpixel \
        --integrationspace qspace --wavelength 0.222 --xbeamcenter 2000.0 --ybeamcenter 2000.0 \
        --distance 100.0 --rotationd 2.0 --tiltd 2.0 --tthstepd 0.01 --qstep 0.01 --includepattern *.tif\
        --excludepattern *.dark.tif *.test.tif --fliphorizontal True --flipvertical False \
        --xdimension 3045 --ydimension 3045 --xpixelsize 0.1 --ypixelsize 0.1 --uncertaintyenable True\
        --sacorrectionenable False --polcorrectionenable False --polcorrectf 0.5 --selfcorrenable False\
        --gsasoutput esd --filenameplus test --maskedges 20 20 20 20 200'.split()
    

    def setUp(self):
        #self.srxplanar = SrXplanar()
        self.picnames = ['pic_001', 'pic_002', 'pic_003', 'Pic_001', 'Pic_002', 'Pic_003']
        self.picnames.sort()
        self.extnames = ['.tif', '.raw.tif', '.dark.tif', '.tif.metadata', '.png']
        os.makedirs('raw')
        for pic, ext in itertools.product(self.picnames, self.extnames):
            f = open(pic+ext, 'w')
            f.close()
            f = open('raw/'+pic+ext, 'w')
            f.close()
        return
    
    def tearDown(self):
        #self.srxplanar = None
        filelist = os.listdir('.')
        for extname in ('*.cfg', '*.tif', '*.chi', '*.fq', '*.gr','*.metadata', '*.png'):
            filelist1 = fnmatch.filter(filelist, extname)
            for f in filelist1:
                os.remove(f)
                os.remove('raw/'+f)
        os.rmdir('raw')
        return
    
    def test_filter(self):
        self.srxplanar= SrXplanar(args = 
            '--tifdirectory ./ --includepattern *.tif\
            --excludepattern *.dark.tif *.raw.tif'.split())
        
        filelist = self.srxplanar.loadimage.genFileList()
        filelist.sort()
        self.assertEqual(filelist, [pp+'.tif' for pp in self.picnames])
        
        self.srxplanar.updateConfig(args = 
            '--tifdirectory ./raw --includepattern *.tif\
            --excludepattern *.dark.tif *.raw.tif'.split())
        filelist = self.srxplanar.loadimage.genFileList()
        filelist.sort()
        self.assertEqual(filelist, [pp+'.tif' for pp in self.picnames])
        
        self.srxplanar.updateConfig(args = 
            'pic*'.split())
        filelist = self.srxplanar.loadimage.genFileList()
        filelist.sort()
        self.assertEqual(filelist, ['pic_001.tif', 'pic_002.tif', 'pic_003.tif'])       
        return

# End of class TestRoutines

if __name__ == '__main__':
    unittest.main()

# End of file
