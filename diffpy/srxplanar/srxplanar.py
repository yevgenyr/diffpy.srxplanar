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

from srxplanarconfig import SrXPlanarConfig
from calculate import Calculate
from loadimage import LoadImage
from mask import Mask
from saveresults import SaveResults

class SrXPlanar(object):
    '''
    '''
    def __init__(self, srxplanarconfig=None):
        if srxplanarconfig!=None:
            if type(srxplanarconfig) == str:
                self.config = SrXPlanarConfig(srxplanarconfig)
            else:
                self.config = srxplanarconfig
        else:
            self.config = SrXPlanarConfig()

        #init modulars
        self.mask = Mask(self.config)
        self.loadimage = LoadImage(self.config)
        self.calculate = Calculate(self.config)
        self.saveresults = SaveResults(self.config)
        #init variables
        self.hasimage = False
        self.prepareCalculation()
        return
    
    @staticmethod
    def configProperty(nm):
        '''helper function that make delegate parameters from srsig2dconfig.
        '''
        rv = property(fget = lambda self: getattr(self.config, nm),
                      fset = lambda self, value: setattr(self.config, nm, value))
        return rv
    
    def updateConfig(self, configfile=None):
        '''update config, rerun all prepareCalculation() for each modulars,
        usually used after you changed some paramters.  
        :configfile str: you can specify a configfile, program will read the config file and update
        '''
        if type(configfile)==str:
            self.config.loadFromFile(configfile)
        #update instances
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
            picvar = self.calculate.varianceLocal(self.pic)
            self.pic = self.pic * self.correction
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
    
def main1(argv=sys.argv):
    #print argv
    #config = SrXPlanarConfig('test.cfg')
    #config.updateConfig()
    xplanar = SrXPlanar('test_s_sr.cfg')
    xplanar.updateConfig()
    #sig2d.newImageFile('KFe2As2-00838.tif')
    xplanar.integrate('CeO2.tif')
    return

def main():
    '''read config and integrate all images
    '''
    configfile = sys.argv[-1] if len(sys.argv)>1 else ''
    if os.path.exists(configfile):
        xplanar = SrXPlanar(configfile)
    elif os.path.exists('srxplanarconfig.cfg'):
        xplanar = SrXPlanar('srxplanarconfig.cfg')
    else:
        print 'please provide config file'
    xplanar.integrateAll()
    return

if __name__=='__main__':
    sys.exit(main())
