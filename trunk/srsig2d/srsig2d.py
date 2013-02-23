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
#import scipy as sp
import scipy.sparse as ssp
#import fabio, fabio.openimage
import os, sys
#import time

from srsig2dconfig import SrSig2DConfig
from calculate import Calculate
from loadimage import LoadImage
from mask import Mask
from saveresults import SaveResults
from tmatrix import Tmatrix

class SrSig2D(object):
    '''
    '''
    def __init__(self, srsig2dconfig=None):
        if srsig2dconfig!=None:
            if type(srsig2dconfig) == str:
                self.config = SrSig2DConfig(srsig2dconfig)
            else:
                self.config = srsig2dconfig
        else:
            self.config = SrSig2DConfig()
        
        self.configlist = self.config.configlist
        #for optionname in self.configlist:
        #    setattr(self.__class__, optionname, self.configProperty(optionname))
        
        #init modulars
        self.tm = Tmatrix(self.config)
        self.mask = Mask(self.config)
        self.loadimage = LoadImage(self.config)
        self.calculate = Calculate(self.config)
        self.saveresults = SaveResults(self.config)
        #init variables
        self.hasimage = False
        self.pic = np.zeros((self.config.ydimension, self.config.xdimension))
        #init Tmatrix
        if srsig2dconfig!=None:
            self.calTmatrix(False)
        return
    
    @staticmethod
    def configProperty(nm):
        '''helper function that make delegate parameters from srsig2dconfig.
        '''
        rv = property(fget = lambda self: getattr(self.config, nm),
                      fset = lambda self, value: setattr(self.config, nm, value))
        return rv
    
    ##############################################
    #Tmatrix
    ##############################################
    ''' tm:             Tmatrix Instance, class that do tmatrix calculation
        tmatrixorg:     csr matrix, store the original tmatrix (the normal mask has already been applied on this tmatrix
        tmatrixvarorg:  csr matrix, store the original tmatrixvar (the normal mask has already been applied on this tmatrixvar
        tmatrix:        csr matrix, store the actual tmatrix used in 2D integration. Every time a new image comes, 
                        it will update selfcorrmask and apply on this tmatrix.
        tmatrixvar:     csr matrix, store the tmatrix used in error calculation.
        tmatrixselfcorr:csr matrix, store the tmatrix used in self correction mask calculation
        number:         1d ndarray, helper array that store the sum of each row of tmatrix
        number:         1d ndarray, helper array that store the sum of each row of tmatrixvar
        xrdxgrid:       1d ndarray, store the x grid of xrd integration either in 2theta space (degree) or qspace(A^-1), only used in output
    '''

    def _genTmatrix(self):
        '''get tmatrix from tmatrixorg (and/or self corr mask) and calculate number array'''
        if (self.config.maskselfcorrenable)and self.hasimage:
            result = np.dot(self.tmatrixorg, self.maskselfcorrdiag).tocsr()
        else:
            result = self.tmatrixorg
        self.number = self.calculate.number(result)
        self.tmatrix = result
    def _genTmatrixvar(self):
        '''get tmatrixvar from tmatrixvarorg (and/or self corr mask) and calculate numbervar array'''
        if self.config.xrduncertaintyenable:
            if (self.config.maskselfcorrenable)and self.hasimage:
                result = np.dot(self.tmatrixvarorg, self.maskselfcorrdiag).tocsr()
            else:
                result = self.tmatrixvarorg
            self.numbervar = self.calculate.number(result)
            self.tmatrixvar = result
    def _genTmatrixorg(self):
        '''generate tmatrixorg and apply normal mask'''
        self.tmatrixorg = self.tm.genTmatrix()
        if self.config.maskenable:
            self.tmatrixorg = np.dot(self.tmatrixorg, self.masknormaldiag)
    def _genTmatrixvarorg(self):
        '''generate tmatrixvarorg and apply normal mask'''
        self.tmatrixvarorg = self.tm.genTmatrixvar()
        if self.config.maskenable:
            self.tmatrixvarorg = np.dot(self.tmatrixvarorg, self.masknormaldiag)
    def _genTmatrixselfcorr(self):
        '''generate tmatrixselfcorr and apply normal mask'''
        self.tmatrixselfcorr = self.tm.genTmatrixselfcorr()
        if self.config.maskenable:
            self.tmatrixselfcorr = np.dot(self.tmatrixselfcorr, self.masknormaldiag)
    def _genXrdxgrid(self):
        '''update xrdxgrid'''
        config = self.config
        if config.integrationspace =='twotheta':
            rv = np.arange(0.0, config.xrdtthmaxd, config.xrdtthstepd)
        elif config.integrationspace == 'qspace':
            rv = np.arange(0.0, config.xrdqmax, config.xrdqstep)
        self.xrdxgrid = rv
        return   
 
    ##############################################
    #2d image
    ##############################################
    ''' store the image in different formats used in calculation
        pic:            2d ndarray, store the latest raw image data array, this pic array is already cropped and flipped.
        picflat:        1d ndarray, flattened image data array, picflat[y*xdimension+x] = pic[y,x]
        picvector:      csc matrix, vertical vector of picflat, shape is (len(picflat),1)
        picvectordiag:  csr matrix, a matrix with diagonal is picflat, shape is (len(picflat), len(picflat))  
    '''
    def _pic_changed(self):
        '''update all pic related data (include self corr mask) when a new image is read
        '''
        self.picflat = self.pic.ravel()
        self.picvector = ssp.csc_matrix (self.picflat.reshape(len(self.picflat),1))
        self.picvectordiag = ssp.spdiags(self.picflat, [0], len(self.picflat), len(self.picflat))
        if self.config.maskselfcorrenable:
            self.maskselfcorr = self.mask.selfcorrMask(self.picvector, self.picflat, self.tmatrixselfcorr)
            self.maskselfcorrdiag = self.mask.selfcorrMaskDiag()
        self.hasimage = True
        #generate tmatrix in the case selfcorr mask changed
        self._genTmatrix()
        self._genTmatrixvar()
        return

    ##############################################
    #mask
    ##############################################
    ''' perpare the mask related data
        mask:           Mask Instance, do mask calculation
        masknormal:     2d boolean ndarray, store the normal mask (predefined mask), 0 stands for masked pixel
        masknormaldiag: csr matrix, matrix with flattened normal mask array stored in diagonal
        maskselfcorr:   2d boolean ndarray, store the self corr mask (mask calculated according to every frame)
                        , 0 stands for masked pixel
        maskselfcorrdiag: csr matrix, matrix with flattened selfcorr mask array stored in diagonal
    '''
    def _genMaskNormal(self):
        '''generate normal mask(predefined mask)'''
        self.masknormal = self.mask.normalMask()
        self.masknormaldiag = self.mask.normalMaskDiag()
        return
    
    ##############################################
    #job
    ##############################################
        
    def calTmatrix(self, reloadimage=True):
        '''calculate base tmatrixs used in calculation. Some tmatrix will automaticly generated
        during the process.''' 
    
        if self.config.maskenable:
            self._genMaskNormal()
        self._genTmatrixorg()
        self._genXrdxgrid()
        if self.config.xrduncertaintyenable:
            self._genTmatrixvarorg()
        if self.config.maskselfcorrenable:
            self._genTmatrixselfcorr()
        if reloadimage and self.hasimage:
            self._pic_changed()
        return
    
    def calChi(self):
        '''calculate xrd intensity
        read the image from self.pic and return chi array and/or xrd VC matirx array
        :return: list, return self.chi(xrd intensity) and self.xrdvc (VC matrix)
                self.chi: self.chi[0], xgrid of xrd intensity, could be twotheta grid or q-grid;
                self.chi[1], xrd intensity; self.chi[2], uncertainty of intensity(if calculated)
                self.xrdvc: scipy.csr_matrix sparse matrix
        '''
        #chi calculation
        chi = self.calculate.intensity(self.picvector, self.tmatrix, self.number)
        if self.config.xrduncertaintyenable:
            vcmatrix = self.calculate.vcmatrixLocal(self.picflat, self.tmatrixvar, 
                            self.tmatrix, self.number)
            variance = np.array(vcmatrix.diagonal())
            std = np.sqrt(variance)
            self.xrdvc = ssp.csr_matrix(vcmatrix)
            self.chi = np.vstack([self.xrdxgrid, chi, std])
            rv = [self.chi, self.xrdvc]
        else:
            self.chi = np.vstack([self.xrdxgrid, chi])
            rv = [self.chi]
        return rv
    
    def newImageFile(self, imagefile):
        '''process new image file, integrate image and save the results
        :imagefile str: filename(path) of new image file
        '''
        self.pic = self.loadimage.loadImage(imagefile)
        self.filename = os.path.splitext(os.path.split(imagefile)[1])[0] + self.config.filenameplus
        self._pic_changed()
        rv = self.calChi()
        self.saveContent()
        return rv 
    
    def newImage(self, imagearray, name):
        '''process new image array, integrate and save
        :imagearray 2darray: image array to be integrated
        :name str: filename used in saving results
        '''
        self.pic = imagearray
        self.filename = os.path.splitext(name)[0] + self.config.filenameplus
        self._pic_changed()
        rv = self.calChi()
        self.saveContent()
        return rv
    
    def saveContent(self, filename=None):
        '''save the content into files
        :filename str: file name, if None, read the file name from self.filename 
        '''
        filename = self.filename if filename==None else filename
        filename = filename+'_'+self.config.integrationspace+'_'+self.config.method
        filepath = os.path.normpath(self.config.savedirectory + '/' + filename)
        self.saveresults.saveChi(self.chi, filepath+'.chi')
        if hasattr(self, 'xrdvc'):
            self.saveresults.saveXrdvc(self.xrdvc, filepath+'.mtx')
        return
    
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
        self.tm.prepareCalculation()
        self.calTmatrix(False)
        return
    
    def integrateAll(self):
        '''integrate all image in self.tifdirectory
        '''
        filelist = self.loadimage.genFileList()
        filelistfull = map(lambda name: os.path.normpath(self.config.tifdirectory+'/'+name), filelist)
        for file1 in filelistfull:
            self.newImageFile(file1)
        return
    
def main1(argv=sys.argv):
    #print argv
    #config = SrSig2DConfig('test.cfg')
    #config.xrduncertaintyenable = True
    #config.updateConfig()
    sig2d = SrSig2D('test.cfg')
    #sig2d.newImageFile('KFe2As2-00838.tif')
    sig2d.integrateAll()
    return

def main():
    '''read config and integrate all images
    '''
    configfile = sys.argv[-1] if len(sys.argv)>1 else ''
    if os.path.exists(configfile):
        sig2d = SrSig2D(configfile)
    elif os.path.exists('srsig2dconfig.cfg'):
        sig2d = SrSig2D('srsig2dconfig.cfg')
    try:
        sig2d.integrateAll()
    except:
        print 'please provide config file'
    return

if __name__=='__main__':
    sys.exit(main())
