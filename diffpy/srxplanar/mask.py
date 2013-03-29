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
    maskedges = _configPropertyR('maskedges')
    wavelength = _configPropertyR('wavelength')
    addmask = _configPropertyR('addmask')
    
    def __init__(self, p):
        self.config = p
        self.prepareCalculation()
        return
    
    def prepareCalculation(self):
        self.xydimension = self.xdimension * self.ydimension
        return

    def normalMask(self, addmask = None):
        """create a mask file which indicate the dead pixel
        return a 2d ndarray with boolean (1 stands for masked pixel)
        
        param addmask: list of string
        """
        addmask = self.addmask if addmask==None else addmask
        
        #right here, '1' stands for masked pixel, in the actual mask array, '1' stands for unmasked pixel
        rv = np.zeros((self.ydimension, self.xdimension))
        #fit2d mask
        maskfit2ds = filter(lambda msk: msk.endswith('.msk'), addmask)
        if len(maskfit2ds)>0:
            for maskfit2d in maskfit2ds:
                if os.path.exists(maskfit2d):
                    immask = fabio.openimage.openimage(maskfit2d)
                    rv += immask.data
        #.npy mask
        npymasks = filter(lambda msk: msk.endswith('.npy'), addmask)
        if len(npymasks)>0:
            for npymask in npymasks:
                if os.path.exists(npymask):
                    rv += np.load(npymask)
        #.tif mask
        tifmasks = filter(lambda msk: msk.endswith('.tif'), addmask)
        if len(tifmasks)>0:
            for tifmask in tifmasks:
                if os.path.exists(tifmask):
                    immask = fabio.openimage.openimage(tifmask)
                    rv += self.flipImage(immask.data)
        #edge mask 
        edgemask = filter(lambda msk: msk.startswith('edge'), addmask)
        if len(edgemask)>0:
            if np.sum(self.maskedges)!=0:
                rv += self.edgeMask(self.maskedges)
        
        self.mask = (rv > 0)
        return self.mask
    
    def dynamicMask(self, pic, addmask = None):
        '''dyamic mask generated according to pic itself.
        
        return: 2d array of boolean, 1 stands for masked pixel
        '''
        addmask = self.addmask if addmask==None else addmask
        #right here, '1' stands for masked pixel, in the actual mask array, '1' stands for unmasked pixel
        rv = np.zeros((self.ydimension, self.xdimension))
        flag = False
        #deadpixel mask
        dpmask = filter(lambda msk: msk.startswith('deadpixel'), addmask)
        if len(dpmask)>0:
            rv += self.deadPixelMask(pic)
            flag = True
        #spot mask
        spmask = filter(lambda msk: msk.startswith('spot'), addmask)
        if len(spmask)>0:
            rv += self.spotMask(pic)
            flag = True
        #return None if none mask applied
        if flag:
            self.dynamicmask = (rv>0)
        else:
            self.dynamicmask = None
        return self.dynamicmask
    
    def deadPixelMask(self, pic):
        '''mask the dead pixel
        return: 1 for masked pixel
        '''
        avgpic = np.average(pic)
        ks = np.ones((5,5))
        ks1 = np.ones((7,7))
        picb = snf.percentile_filter(pic, 5, 3) < avgpic/10
        picb = snm.binary_dilation(picb, structure=ks)
        picb = snm.binary_erosion(picb, structure=ks1)
        return picb
    
    def spotMask(self, pic, size=5, r = 1.2):
        '''mask the spot in image
        return: 1 for masked pixel
        '''
        rank = snf.rank_filter(pic, -size, size)
        ind = snm.binary_dilation(pic>rank*r, np.ones((3,3)))
        return ind
    
    def edgeMask(self, edges=None):
        '''number in edges stands for the number of masked pixels
        left, right, top, bottom, corner
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
        
        ra = edges[4]
        ball = np.zeros((ra*2, ra*2))
        radi = (np.arange(ra*2)-ra).reshape((1, ra*2))**2 + \
                (np.arange(ra*2)-ra).reshape((ra*2, 1)) ** 2
        radi = np.sqrt(radi)
        ind = radi > ra
        rv[edges[3]:edges[3]+ra, edges[0]:edges[0]+ra] = ind[:ra,:ra]
        rv[edges[3]:edges[3]+ra, -edges[1]-ra:-edges[1]] = ind[:ra,-ra:]
        rv[-edges[2]-ra:-edges[2], edges[0]:edges[0]+ra] = ind[-ra:, :ra]
        rv[-edges[2]-ra:-edges[2], -edges[1]-ra:-edges[1]] = ind[-ra:,-ra:]
        return rv
    
    def undersample(self, undersamplerate):
        '''a special mask used for undesampling image. It will create a mask that
        discard (total number*(1-undersamplerate)) pixels
        undersamplerate: 0~1
        '''
        ind = np.random.permutation(len(self.picflat))[:len(self.picflat)*undersamplerate]
        n = self.xdimension * self.ydimension
        mask = np.zeros_like(n, dtype = bool)
        mask[ind] = True
        mask1 = ssp.spdiags((mask), [0], n, n).tocsr()
        return mask1
        
    def flipImage(self, pic):
        '''flip image if configured in config 
        '''
        if self.fliphorizontal:
            pic = pic[:,::-1]
        if self.flipvertical:
            pic = pic[::-1,:]
        return pic
    
    def saveMask(self, filename, pic=None, addmask=None):
        '''generate a mask according to the addmask and pic. save it to .npy. in .npy, 1 stands for masked pixel
        the mask has same ort as the pic, which means if the pic is fliped, the mask is fliped
        '''
        if not hasattr(self, 'mask'):
            self.normalMask(addmask)
        if (not hasattr(self, 'dynamicmask')) and (pic!=None):
            self.dynamicMask(pic, addmask=addmask)
        tmask = np.logical_or(self.mask, self.dynamicmask) if pic!=None else self.mask
        np.save(filename, tmask)
        return tmask