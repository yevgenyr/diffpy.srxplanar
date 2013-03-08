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
import scipy.sparse as ssp
import fabio.openimage
import scipy.ndimage.filters as snf
import scipy.ndimage.morphology as snm
import os
from diffpy.srxplanar.srxplanarconfig import _configPropertyR

class Mask(object):
    '''module to provide mask support. it provide following functions:
    creat mask from fit2d mask file/tif file
    creat mask by cake cutting or box cutting
    creat mask by filting out the pixels with too high or too low intensity (selfcorr)
    for all masks, 1 stands for unmasked pixel, 0 stands for masked pixel
    
    call *Mask() to return a 2d ndarray 
    call *MaskDiag() to return a scipy.sparse matrix with main diagonal equal to the 
    
    normalMask() & normalMaskDiag() to return normal mask
    selfcorrMask() & selfcorrMaskDiag() to return selfcorr mask
    '''
    
    xdimension = _configPropertyR('xdimension')
    ydimension = _configPropertyR('ydimension')
    fliphorizontal = _configPropertyR('fliphorizontal')
    flipvertical = _configPropertyR('flipvertical')
    maskfit2d = _configPropertyR('maskfit2d')
    maskedges = _configPropertyR('maskedges')
    
    def __init__(self, p, cal):
        self.config = p
        self.calculate = cal
        self.prepareCalculation()
        return
    
    def prepareCalculation(self):
        self.tthorqmatrix = self.calculate.tthorqmatrix
        self.tthorqstep = self.calculate.tthorqstep
        return
    
    def configProperty(self, nm):
        '''helper function of property delegation
        '''
        rv = property(fget = lambda self: getattr(self.config, nm))
        return rv
    
    def normalMask(self):
        """create a mask file which indicate the dead pixel
        return a 2d ndarray with boolean (1 stands for unmasked pixel, 0 stands for masked pixel)
        
        fit2d:read the fit2d mask file
        """
        rv = np.zeros((self.ydimension, self.xdimension))
        #right here, '1' stands for masked pixel, in the actual mask array, '1' stands for unmasked pixel 
        if self.maskfit2d != None:
            if os.path.exists(self.maskfit2d):
                immask = fabio.openimage.openimage(self.maskfit2d)
                rv = self.flipImage(immask.data)
        if np.sum(self.maskedges)!=0:
            rv = rv + self.edgeMask(self.maskedges)
        #rv = self.edgeMaskHalf(rv)
        #rv[0:100, 0:100] = 1
        self.mask = (rv == 0)
        #self.mask = np.zeros_like(rv)
        #self.mask[::2,::2] = 1
        #self.mask[1::2,1::2] = 1
        return self.mask
    
    def selfMask(self, pic):
        '''return: 0 for masked pixel
        '''
        avgpic = np.average(pic)
        ks = np.ones((5,5))
        ks1 = np.ones((7,7))
        picb = snf.percentile_filter(pic, 5, 3) < avgpic/10
        picb = snm.binary_dilation(picb, structure=ks)
        picb = snm.binary_erosion(picb, structure=ks1)
        picb = np.logical_not(picb)
        return picb

    def edgeMaskHalf(self, curmask):
        '''input: 1 for masked pixel
        return: 1 for masked pixel
        '''
        curmask1 = np.zeros_like(curmask, dtype=bool) + curmask
        curmask1[int(self.ydimension * 0.2):int(self.ydimension*0.8),
                 int(self.xdimension * 0.2):int(self.xdimension*0.8)] = 0
        curmask1[0,:] = 1
        curmask1[-1,:] = 1
        curmask1[:,0] = 1
        curmask1[:,-1] = 1
        ks = np.ones((3,3))
        curmask2 = snm.binary_dilation(curmask1, ks)
        
        #curmask2 = np.logical_and(curmask2, np.logical_not(curmask1))
        tth = (self.tthorqmatrix / self.tthorqstep) % 1.0
        tth = tth > 0.5
        mask = np.logical_and(tth, curmask2)
        mask = np.logical_or(curmask, mask)
        
        '''from pylab import *
        imshow(mask, interpolation='none')
        show()'''
        
        return mask
        

    def edgeMask(self, edges=None):
        '''number in edges stands for the number of masked pixels
        left, right, top, bottom
        return: 1 for masked pixel
        '''
        edges = self.maskedges if edges==None else edges
        rv = np.zeros((self.ydimension, self.xdimension))
        if edges[0]!=0:
            rv[:,:edges[0]] = 1
        if edges[1]!=0:
            rv[:,-edges[1]:] = 1
        if edges[2]!=0:
            rv[-edges[2]:,:] = 1
        if edges[3]!=0:
            rv[:edges[3]:,:] = 1
        return rv
    
    def flipImage(self, pic):
        '''flip image if configured in config 
        '''
        if self.fliphorizontal:
            pic = pic[:,::-1]
        if self.flipvertical:
            pic = pic[::-1,:]
        return pic
    