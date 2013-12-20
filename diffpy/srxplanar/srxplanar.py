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
        '''
        init srxplanar form a SrXplanarConfig instance, or config file, or args passed from cmd
        or kwargs. If both SrXplanarConfig instance and other configfile/args/kwargs is specified, 
        it will first init from config instance then update using configfile/args/kwargs
        
        :param srxplanarconfig: SrXplanarConfig, init srxplanar from a config instance
        :param configfile: string, name of config file
        :param args: list of str, usually be sys.argv
        :param kwargs: you can use like 'xbeamcenter=1024' or a dict to update the value of xbeamcenter
        '''
        if srxplanarconfig!=None:
            self.config = srxplanarconfig
            self.config.updateConfig(filename=configfile, args=args, **kwargs)
        else:
            self.config = SrXplanarConfig(filename=configfile, args=args, **kwargs)
        #init modulars
        self.mask = Mask(self.config)
        self.loadimage = LoadImage(self.config)
        self.calculate = Calculate(self.config)
        self.saveresults = SaveResults(self.config)
        return
    
    def updateConfig(self, filename=None, args=None, **kwargs):
        '''
        update config using configfile/args/kwargs, then rerun all prepareCalculation()
        
        :param configfile: string, name of config file
        :param args: list of str, usually be sys.argv
        :param kwargs: you can use like 'xbeamcenter=1024' or a dict to update the value of xbeamcenter
        
        :return: None
        '''
        self.config.updateConfig(filename=filename, args=args, **kwargs)
        #update instances
        self.calculate.prepareCalculation()
        self.saveresults.prepareCalculation()
        return
        
    def prepareCalculation(self, pic=None):
        '''
        prepare data used in calculation
        
        :param reloadimage: boolean, if True, then recalculate data related to image data such as dynamic mask
        
        :return: None
        '''
        self.staticmask = self.mask.staticMask()
        self.correction = self.calculate.genCorrectionMatrix()
        # if a pic is provided, then generate one-time dynamicmask 
        if pic != None:
            image = self._getPic(pic)
            image *= self.correction
            dymask = self.mask.dynamicMask(image, addmask = ['dead', 'bright'])
            dymask = np.logical_or(self.staticmask, dymask)
            self.calculate.genIntegrationInds(dymask)
            chi = self.calculate.intensity(image)
            index = np.rint(self.calculate.tthorqmatrix / self.config.tthorqstep).astype(int)
            index[index>=len(chi[1])] = len(chi[1]-1)
            avgimage = chi[1][index.ravel()].reshape(index.shape)
            mask = np.logical_or(image<avgimage*0.7, image>avgimage*1.5)
            self.staticmask = np.logical_or(np.logical_or(self.staticmask, mask), dymask)
        
        self.calculate.genIntegrationInds(self.staticmask)
        return
    
    def _picChanged(self):
        '''
        update all pic related data (such as dynamic mask) when a new image is read
        
        :return: None
        '''
        dynamicmask = self.mask.dynamicMask(self.pic)
        if dynamicmask != None:
            mask = np.logical_or(self.staticmask, dynamicmask)
            self.calculate.genIntegrationInds(mask)
        return
    
    def _getSaveFileName(self, imagename=None, filename=None):
        '''
        get the save file name, the priority order is self.output> filename> imagename > 'output'(default name)
        
        :param imagename: string, filename/path of image file (drop this term if it is an image array)
        :param filename: string, 
        
        :return: string, name of file to be saved 
        '''
        rv = 'output'
        if self.config.output!=None and self.config.output!='':
            rv = self.config.output
        elif filename!=None:
            rv = filename
        elif imagename!=None and isinstance(imagename, (str, unicode)):
            rv = imagename
        return rv
    
    def _getPic(self, image, flip=True):
        if isinstance(image, list):
            rv = np.zeros((self.config.ydimension, self.config.xdimension))
            for imagefile in image:
                rv += self._getPic(imagefile)
        elif isinstance(image, (str, unicode)):
            rv = self.loadimage.loadImage(image)
        elif flip:
            rv = self.loadimage.flipImage(image)
        else:
            rv = image
        return rv
    
    def integrate(self, image, savename=None, savefile=True, flip=True):
        '''
        integrate 2d image to 1d diffraction pattern, then save to disk
        
        :param image: str or 2d array, 
            if str, then read image file using it as file name.
            if 2d array, integrate this 2d array.
        :param savename: str, name of file to save
        :param savefile: boolean, if True, save file to disk, if False, do not save file to disk
        :param flip: boolean, if True and 'image' is a 2d array, flip this array and integrate it
            if False and 'image' is a 2d array, directly integrate it. 
        
        :return: dict, rv['chi'] is a 2d array of integrated intensity, shape is (2, len of intensity) 
            or (3, len of intensity) in [tth or q, intensity, (uncertainty)]. rv['filename'] is the 
            name of file to save to disk
        '''
        rv = {}
        self.pic = self._getPic(image, flip)
        rv['filename'] = self._getSaveFileName(imagename= image, filename=savename)
        self._picChanged()
        #calculate
        rv['chi'] = self.chi = self.calculate.intensity(self.pic)
        #save
        if savefile:
            rv['filename'] = self.saveresults.save(rv)
        return rv
    
    def integrateFilelist(self, filelist, summation=None, filename=None):
        '''
        process all file in filelist, integrate them separately or together
        
        :param filelist: list of string, file list (full path)
        :param summation: bool or None, sum all files together or not, if None,
            use self.config.summation
        :param filename: file name of output file 
        '''
        summation =  self.config.summation if summation == None else summation
        if (summation)and(len(filelist)>1):
            image = self._getPic(filelist)
            filename = os.path.splitext(filelist[-1])[0]+'_sum.chi' if filename == None else filename
            rv = [self.integrate(image, savename = filename, flip=False)]
        else:
            i = 0
            rv = []
            for imagefile in filelist:
                if filename==None:
                    rvv = self.integrate(imagefile)
                else:
                    rvv = self.integrate(imagefile, savename = filename+'%03d'%i)
                rv.append(rvv)
        return rv
    
    def process(self):
        '''
        process the images according to filenames/includepattern/excludepattern/summation
        by default, it will scan current/tifdirectory and integrate all files match 
        includepattern/excludepattern and/or filenames.
        
        Usually this one is called from cmd line rather then script.
        
        :return: None
        '''
        if not self.config.nocalculation:
            filelist = self.loadimage.genFileList()
            if len(filelist)>0:
                self.prepareCalculation(pic = filelist[0])
                self.integrateFilelist(filelist)
            else:
                print 'No input files or configurations'
                self.config.args.print_help()
        #mask creating
        elif self.config.createmask!='':
            self.createMask()
        #if no config is passed to srxplanar
        else:
            print 'No input files or configurations'
            self.config.args.print_help()
        return
    
    def createMask(self, filename= None, pic=None, addmask=None):
        '''
        create and save a mask according to addmask, pic, 1 stands for masked pixel in saved file
        
        :param filename: name of mask file to save, 'mask.npy' if it is None
        :param pic: 2d image array, may used in generating dynamic mask, Be careful if this one is flipped or not
        :param addmask: list of str, control how to generate mask, see Mask module for detail
        
        :return: 2d array, 0 stands for masked pixel here
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
    


def main():
    '''
    read config and integrate images
    '''
    srxplanar = SrXplanar(args=sys.argv[1:])
    srxplanar.process()
    return

if __name__=='__main__':
    sys.exit(main())
