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
    '''main modular for srxplanar
    '''
    
    def __init__(self, srxplanarconfig=None, configfile=None, args=None, **kwargs):
        '''init srxplanar form a SrXplanarConfig instance, or config file, or args passed from cmd
        or kwargs. If both SrXplanarConfig instance and other configfile/args/kwargs is specified, 
        it will first init from config instance then update using configfile/args/kwargs
        
        param srxplanarconfig: SrXplanarConfig, init srxplanar from a config instance
        param configfile: string, name of config file
        param args: list of str, usually be sys.argv
        param **kwargs: you can use like 'xbeamcenter=1024' or a dict to update the value of xbeamcenter
        '''
        if srxplanarconfig!=None:
            self.config = srxplanarconfig
            self.config.updateConfig(filename=configfile, args=args, **kwargs)
        else:
            self.config = SrXplanarConfig(filename=configfile, args=args, **kwargs)
        #init modulars
        if not self.config.nocalculation:
            self.mask = Mask(self.config)
            self.loadimage = LoadImage(self.config)
            self.calculate = Calculate(self.config)
            self.saveresults = SaveResults(self.config)                
            self.prepareCalculation(False)
        return
    
    def updateConfig(self, filename=None, args=None, **kwargs):
        '''update config using configfile/args/kwargs, then rerun all prepareCalculation()
        
        param configfile: string, name of config file
        param args: list of str, usually be sys.argv
        param **kwargs: you can use like 'xbeamcenter=1024' or a dict to update the value of xbeamcenter
        '''
        self.config.updateConfig(filename=filename, args=args, **kwargs)
        #update instances
        self.mask.prepareCalculation()
        self.loadimage.prepareCalculation()
        self.calculate.prepareCalculation()
        self.saveresults.prepareCalculation()
        self.prepareCalculation(False)
        return
        
    def prepareCalculation(self, reloadimage=True):
        '''prepare data used in calculation
        
        param reloadimage: boolean, if True, then recalculate data related to image data such as dynamic mask
        '''
        self.staticmask = self.mask.staticMask()
        self.correction = self.calculate.genCorrectionMatrix()
        if reloadimage:
            dynamicmask = self.mask.dynamicMask(self.pic)
            if dynamicmask != None:
                mask = np.logical_or(self.staticmask, dynamicmask)
                self.calculate.genIntegrationInds(mask)
            else:
                self.calculate.genIntegrationInds(self.staticmask)
        else:
            self.calculate.genIntegrationInds(self.staticmask)
        return
    
    def _picChanged(self):
        '''update all pic related data (such as dynamic mask) when a new image is read
        '''
        dynamicmask = self.mask.dynamicMask(self.pic)
        if dynamicmask != None:
            mask = np.logical_or(self.staticmask, dynamicmask)
            self.calculate.genIntegrationInds(mask)
        return
    
    
    def _getSaveFileName(self, imagename=None, filename=None):
        '''get the save file name, the priority order is self.output> filename> imagename > 'output'(default name)
        
        param imagename: string, filename/path of image file (drop this term if it is an image array)
        param filename: string, 
        
        return: string, string, name of file to be saved 
        '''
        rv = 'output'
        if self.config.output!=None and self.config.output!='':
            rv = self.config.output
        elif filename!=None:
            rv = filename
        elif imagename!=None and type(imagename)==str:
            rv = imagename
        return rv
    
       
    def integrate(self, image, savename=None, savefile=True, flip=True):
        ''' integrate 2d image to 1d diffraction pattern, then save to disk
        
        param image: str or 2d array, 
            if str, then read image file using it as file name
            if 2d array, integrate this 2d array
        param savename: str, name of file to save
        param savefile: boolean, if True, save file to disk, if False, do not save file to disk
        param flip: boolean, if True and 'image' is a 2d array, flip this array and integrate it
            if False and 'image' is a 2d array, directly integrate it. 
        
        return: dict, rv['chi'] is a 2d array of integrated intensity, shape is (2, len of intensity) 
            or (3, len of intensity) in [tth or q, intensity, (uncertainty)]. rv['filename'] is the 
            name of file to save to disk
        '''
        rv = {}
        if type(image)==str:
            self.pic = self.loadimage.loadImage(image)
        elif flip:
            self.pic = self.loadimage.flipImage(image)
        else:
            self.pic = image
        rv['filename'] = self._getSaveFileName(imagename= image, filename=savename)
        self._picChanged()
        #calculate
        if self.config.uncertaintyenable:
            self.pic = self.pic * self.correction
            picvar = self.calculate.varianceLocal(self.pic)
            rv['chi'] = self.chi = self.calculate.intensity(self.pic, picvar)
        else:
            self.pic = self.pic * self.correction
            rv['chi'] = self.chi = self.calculate.intensity(self.pic)
        if savefile:
            self.saveresults.save(rv)
        return rv
    
    def createMask(self, filename= None, pic=None, addmask=None):
        '''create and save a mask according to addmask, pic, 1 stands for masked pixel in saved file
        
        param filename: name of mask file to save, 'mask.npy' if it is None
        param pic: 2d image array, may used in generating dynamic mask, Be careful if this one is flipped or not
        param addmask: list of str, control how to generate mask, see Mask module for detail
        
        return: 2d array, 0 stands for masked pixel here
        '''
        filename = self.config.createmask if filename==None else filename
        filename = 'mask.npy' if filename =='' else filename
        addmask = self.config.addmask if addmask==None else addmask
        if not hasattr(self, 'mask'):
            self.mask = Mask(self.config)
        if not hasattr(self, 'loadimage'):
            self.loadimage = LoadImage(self.config)
        if pic==None:
            filelist = self.loadimage.genFileList()
            if hasattr(self, 'pic'):
                if self.pic!=None:
                    pic =self.pic
                else:
                    pic = self.loadimage.loadImage(filelist[0]) if len(filelist)>0 else None
            else:
                pic = self.loadimage.loadImage(filelist[0]) if len(filelist)>0 else None
        rv = self.mask.saveMask(filename, pic, addmask)
        return rv
    
    def process(self):
        '''process the images according to filenames/includepattern/excludepattern/summation
        by default, it will scan current/tifdirectory and integrate all files match 
        includepattern/excludepattern and/or filenames.
        
        Usually this one is called from cmd line rather then script.
        '''
        if not self.config.nocalculation:
            filelist = self.loadimage.genFileList()
            if (self.config.summation)and(len(filelist)>1):
                image = np.zeros((self.config.ydimension, self.config.xdimension))
                for imagefile in filelist:
                    rv += self.loadimage.loadImage(imagefile)
                self.integrate(rv, imagefile)
            else:
                for imagefile in filelist:
                    self.integrate(imagefile)
        #mask creating
        elif self.config.createmask!='':
            self.createMask()
        return


def main():
    '''read config and integrate images
    '''
    srxplanar = SrXplanar(args=sys.argv[1:])
    srxplanar.process()
    return

if __name__=='__main__':
    sys.exit(main())
