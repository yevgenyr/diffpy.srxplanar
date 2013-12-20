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

import numpy as np
import scipy.sparse as ssp
import fabio.openimage
import scipy.ndimage.filters as snf
import scipy.ndimage.morphology as snm
import os
from diffpy.srxplanar.srxplanarconfig import _configPropertyR

class Mask(object):
    '''
    provide methods for mask generation, including:
    
    static mask: fit2d (.msk) mask, tif mask, npy mask, masking edge pixels, 
    dymanic mask: masking dead pixels, bright pixels
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
        return
    
    def staticMask(self, addmask = None):
        '''
        create a static mask according existing mask file. This mask remain unchanged for different images
        
        :param addmask: list of string, file name of mask and/or 'edge', 
            mask file supported: .msk, .npy, .tif file, ATTN: mask array should be already flipped, 
            and 1 (or larger) stands for masked pixels, 0(<0) stands for unmasked pixels
            if 'edge' is specified here. it will create a mask that mask the pixel near the edge of detector,
            require self.maskedges 
        
        :return: 2d array of boolean, 1 stands for masked pixel
        '''
        addmask = self.addmask if addmask==None else addmask
        
        rv = np.zeros((self.ydimension, self.xdimension))
        #fit2d mask
        maskfit2ds = filter(lambda msk: msk.endswith('.msk'), addmask)
        if len(maskfit2ds)>0:
            for maskfit2d in maskfit2ds:
                if os.path.exists(maskfit2d):
                    immask = fabio.openimage.openimage(maskfit2d)
                    #rv += self.flipImage(immask.data)
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
        
        self.staticmask = (rv > 0)
        return self.staticmask
    
    def dynamicMask(self, pic, addmask = None):
        '''
        create a dynamic mask according to image array. This mask changes for different images
        
        :param pic: 2d array, image array to be processed
        :param addmask: list of string, ['deadpixel', 'brightpixel']  
            deadpixel: pixels with much lower intensity compare to adjacent pixels will be masked
            brightpixel: pixels with much higher intensity compare to adjacent pixels will be masked 
             
        :return: 2d array of boolean, 1 stands for masked pixel
        '''
        
        addmask = self.addmask if addmask==None else addmask
        rv = np.zeros((self.ydimension, self.xdimension))
        flag = False
        #deadpixel mask
        dpmask = filter(lambda msk: msk.startswith('dead'), addmask)
        if len(dpmask)>0:
            rv += self.deadPixelMask(pic)
            flag = True
        #bright pixel mask
        bpmask = filter(lambda msk: msk.startswith('bright'), addmask)
        if len(bpmask)>0:
            rv += self.brightPixelMask(pic)
            flag = True
        #return None if none mask applied
        if flag:
            self.dynamicmask = (rv>0)
        else:
            self.dynamicmask = None
        return self.dynamicmask
    
    def deadPixelMask(self, pic):
        '''
        pixels with much lower intensity compare to adjacent pixels will be masked
        
        :param pic: 2d array, image array to be processed
        
        :return: 2d array of boolean, 1 stands for masked pixel
        '''
        avgpic = np.average(pic)
        ks = np.ones((5,5))
        ks1 = np.ones((7,7))
        picb = snf.percentile_filter(pic, 5, 3) < avgpic/10
        picb = snm.binary_dilation(picb, structure=ks)
        picb = snm.binary_erosion(picb, structure=ks1)
        return picb
    
    def brightPixelMask(self, pic, size=5, r = 1.2):
        '''
        pixels with much higher intensity compare to adjacent pixels will be masked,
        this mask is used when there are some bright spots/pixels whose intensity is higher 
        than its neighbors but not too high. Only use this on a very good powder averaged 
        data. Otherwise it may mask wrong pixels. 
        
        This mask has similar functions as 'selfcorr' function. However, this mask will only 
        consider pixels' local neighbors pixels and tend to mask more pixels. While 'selfcorr' 
        function compare one pixel to other pixels in same bin.
        
        :param pic: 2d array, image array to be processed
        :param size: int, size of local area to test if a pixel is a bright pixel
        :param r: float, a threshold for masked pixels   
        
        :return: 2d array of boolean, 1 stands for masked pixel
        '''
        rank = snf.rank_filter(pic, -size, size)
        ind = snm.binary_dilation(pic>rank*r, np.ones((3,3)))
        return ind
    
    def edgeMask(self, edges=None):
        '''
        mask the pixels near edge and around corner
        
        :param edges: list of int (length of 5), first 4 are numbers of pixels masked at each edge
            in (left, right, top, bottom), last one is the radius of round cut at the corner
        
        :return: 2d array of boolean, 1 stands for masked pixel
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
        '''
        a special mask used for undesampling image. It will create a mask that
        discard (total number*(1-undersamplerate)) pixels
        :param undersamplerate: float, 0~1, ratio of pixels to keep
        
        :return: 2d array of boolean, 1 stands for masked pixel
        '''
        mask = np.random.rand(self.ydimension, self.xdimension) < undersamplerate
        return mask
        
    def flipImage(self, pic):
        '''
        flip image if configured in config
        
        :param pic: 2d array, image array
        
        :return: 2d array, flipped image array
        '''
        if self.fliphorizontal:
            pic = pic[:,::-1]
        if self.flipvertical:
            pic = pic[::-1,:]
        return pic
    
    def saveMask(self, filename, pic=None, addmask=None):
        '''
        generate a mask according to the addmask and pic. save it to .npy. 1 stands for masked pixel
        the mask has same order as the pic, which means if the pic is flipped, the mask is fliped
        (when pic is loaded though loadimage, it is flipped)
        
        :param filename: str, filename of mask file to be save
        :param pic: 2d array, image array
        :param addmask: list of str, control which mask to generate
        
        :return: 2d array of boolean, 1 stands for masked pixel
        '''
        if not hasattr(self, 'mask'):
            self.normalMask(addmask)
        if (not hasattr(self, 'dynamicmask')) and (pic!=None):
            self.dynamicMask(pic, addmask=addmask)
        tmask = self.mask
        if hasattr(self, 'dynamicmask'):
            if self.dynamicmask!=None:
                tmask = np.logical_or(self.mask, self.dynamicmask) if pic!=None else self.mask
        np.save(filename, tmask)
        return tmask