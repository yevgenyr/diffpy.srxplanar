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
'''SrSig2D main modular
'''

import numpy as np
import scipy.sparse as ssp
import os, sys
#import time

from diffpy.srxplanar.srxplanarconfig import SrXplanarConfig
from diffpy.srxplanar.calculate import Calculate
from diffpy.srxplanar.loadimage import LoadImage
from diffpy.srxplanar.mask import Mask
from diffpy.srxplanar.saveresults import SaveResults

class SrXplanar(object):
    '''
    '''
    def __init__(self, srxplanarconfig=None):
        if srxplanarconfig!=None:
            if type(srxplanarconfig) == str:
                self.config = SrXplanarConfig(srxplanarconfig)
            else:
                self.config = srxplanarconfig
        else:
            self.config = SrXplanarConfig()

        #init modulars
        self.loadimage = LoadImage(self.config)
        self.calculate = Calculate(self.config)
        self.mask = Mask(self.config, self.calculate)
        self.saveresults = SaveResults(self.config)
        #init variables
        self.prepareCalculation()
        return
    
    def updateConfig(self, configfile=None):
        '''update config, rerun all prepareCalculation() for each modulars,
        usually used after you changed some paramters.  
        :configfile str: you can specify a configfile, program will read the config file and update
        '''
        if type(configfile)==str:
            self.config.loadFromFile(configfile)
        #update instances
        self.config.updateConfig()
        self.calculate.prepareCalculation()
        self.loadimage.prepareCalculation()
        self.mask.prepareCalculation()
        self.saveresults.prepareCalculation()
        return
    
    def prepareCalculation(self):
        '''prepare data used in calculation
        '''
        masknormal = self.mask.normalMask()
        self.correction = self.calculate.genCorrectionMatrix()
        self.calculate.genIntegrationInds(masknormal)
        return
    
    def integrate(self, image, filename=None, savefile=True):
        if type(image)==str:
            self.pic = self.loadimage.loadImage(image)
            self.filename = image
        else:
            self.pic = image
            self.filename = filename if filename!=None else 'output'
        
        if self.config.uncertaintyenable:
            self.pic = self.pic * self.correction
            picvar = self.calculate.varianceLocal(self.pic)
            self.chi = self.calculate.intensity(self.pic, picvar)
        else:
            self.pic = self.pic * self.correction
            self.chi = self.calculate.intensity(self.pic)
        if savefile:
            self.saveresults.save(self.chi, self.filename)
        return self.chi
    
    def integrateAll(self):
        '''integrate all image in self.tifdirectory
        '''
        filelist = self.loadimage.genFileList()
        filelistfull = map(lambda name: os.path.normpath(self.config.tifdirectory+'/'+name), filelist)
        for file1 in filelistfull:
            self.integrate(file1)
        return
    
def main1():
    configfile = sys.argv[-1] if len(sys.argv)>1 else ''
    if os.path.exists(configfile):
        xplanar = SrXplanar(configfile)
    xplanar.integrateAll()
    return

def main():
    '''read config and integrate all images
    '''
    configfile = None
    # start with a default configuration file if it exists
    if os.path.isfile('srxplanarconfig.cfg'):
        configfile = 'srxplanarconfig.cfg'
    # use user's file no matter what she provided
    if len(sys.argv) > 1:
        configfile = sys.argv[1]
    # check for -h, --help options; this should be replaced with optparse
    helprequested = set(['-h', '--help']).intersection(sys.argv[1:])
    # here configfile is None if the default does not exist
    # and user did not give any argument
    if configfile is None or helprequested:
        print 'usage: %s [srxplanarconfig.cfg]'
        print 'Please provide configuration file.'
        sys.exit()
    # configfile is set to something here
    xplanar = SrXplanar(configfile)
    xplanar.integrateAll()
    return

if __name__=='__main__':
    sys.exit(main())
