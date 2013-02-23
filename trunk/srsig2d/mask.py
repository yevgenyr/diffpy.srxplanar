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
import os

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
    def __init__(self, p):
        self.config = p
        self.configlist = ['xdimension',
                           'ydimension',
                           'wavelength',
                           'xpixelsize',
                           'ypixelsize',
                           'xbeamcenter',
                           'ybeamcenter',
                           'distance',
                           'rotation',
                           'tilt',
                           'fliphorizontal',
                           'flipvertical',
                           'maskfit2denable',
                           'maskfit2d',
                           'masktiffenable',
                           'masktiff',
                           'maskcakeenable',
                           'maskcake',
                           'maskboxenable',
                           'maskbox',
                           'maskselfcorrmethod',
                           'maskselfcorrcounthi',
                           'maskselfcorrcountlow',
                           'maskselfcorrpercentagehi',
                           'maskselfcorrpercentagelow',
                           'maskselfcorraccuenable']
        
        for optionname in self.configlist:
            if hasattr(self.config, 'xrd'+optionname):
                setattr(self.__class__, optionname, self.configProperty('xrd'+optionname))
            elif hasattr(self.config, optionname):
                setattr(self.__class__, optionname, self.configProperty(optionname))
        self.prepareCalculation()
        return
    
    def prepareCalculation(self):
        self.xydimension = self.xdimension * self.ydimension
        self.tthmatrix = self.genTTHMatrix()
        
        xr = (np.arange(self.xdimension, dtype=float) - self.xbeamcenter) * self.xpixelsize
        yr = (np.arange(self.ydimension, dtype=float) - self.ybeamcenter) * self.ypixelsize
        self.azimuthmatrix = np.arctan2(yr.reshape(self.ydimension,1),xr.reshape(1,self.xdimension))
        self.mask = np.zeros((self.ydimension, self.xdimension), dtype = bool)
        self.selfcorrmask = np.zeros((self.ydimension, self.xdimension), dtype = bool)
        if self.maskselfcorraccuenable:
            self.selfcorrmaskhistory = np.zeros((self.ydimension, self.xdimension), dtype = int)
        data = np.ones(self.xydimension)
        self.onevector = ssp.csc_matrix(data.reshape(self.xydimension,1), dtype = bool)
        return
    
    def configProperty(self, nm):
        '''helper function of property delegation
        '''
        rv = property(fget = lambda self: getattr(self.config, nm))
        return rv
    
    def normalMaskDiag(self):
        '''return a scipy.sparse matrix with main diagonal equal to the flattened mask array
        '''
        #rv = ssp.spdiags(self.mask.ravel(), [0], self.xydimension, self.xydimension).tocsr()
        rv = ssp.spdiags(self.mask.ravel(), [0], self.xydimension, self.xydimension)
        return rv
    
    def normalMask(self):
        """create a mask file which indicate the dead pixel
        return a 2d ndarray with boolean (1 stands for unmasked pixel, 0 stands for masked pixel)
        
        fit2d:
            read the fit2d mask file
        tiff:
            read a tiff image as a mask, 1 in tiff file stands for masked pixels
        cake:
            cake cut of image file, 
            [radius, start azimuth value of cut area, end azimuth value of cut area]
            note: the azimuth range of image is between (-pi,pi)
        box:
            box cut of image file,could be a series of boxes that desired to cutout
            [[i][0]:d[i][1],d[i][2]:d[i][3]] ->  [x1,x2,y1,y2] of ith box
            pic[x1:x2,y1:y2] will be cut 
        """
        re = np.zeros((self.xdimension, self.ydimension), dtype = bool)
        #right here, '1' stands for masked pixel, in the actual mask array, '1' stands for unmasked pixel 
        if (self.maskfit2d != '')and(os.path.exists(self.maskfit2d)):
            immask = fabio.openimage.openimage(self.maskfit2d)
            fit2dim = self.cropAndFlipImage(immask.data)
            re += fit2dim
        if self.masktiffenable:
            immask = fabio.openimage.openimage(self.masktiff)
            tiffim = self.cropAndFlipImage(immask.data)
            re += tiffim
        if self.maskcakeenable:
            cake = self.tthmatrix > self.maskcake[0]
            azimuth1 = np.radians(self.maskcake[1])
            azimuth2 = np.radians(self.maskcake[2])
            azimuthmin = np.min((azimuth1, azimuth2))
            azimuthmax = np.max((azimuth1, azimuth2))
            cake += np.logical_and((self.azimuthmatrix > azimuthmin), (self.azimuthmatrix < azimuthmax))
            re += cake
        if self.maskboxenable:
            d = self.maskbox
            box = np.zeros_like(re)
            for i in range(len(d)/4):
                x = i * 4
                box[d[x]:d[x+1],d[x+2]:d[x+3]] = 1
            re += box
        self.mask = (re == 0)
        return self.mask
    
    def selfcorrMask(self, picvector, picflat, tmatrix):
        '''return 2d ndarray selfcorr mask using method determined by self.selfcorrmethod
        supported method:'Fast', 'Filter by raw counts', 'Filter by percentage', 
        'Filter by raw counts/percentage'.
        
        'Fast':     First calculate the average raw counts with all pixels in one ring, then 
                    filter out the pixels with raw counts larger or lower than intensity*(1+-filter)
        'Raw counts': First calculate the median number of raw counts in one ring,
                    then filter out pixels with raw counts larger or lower than intensity*(1+-filter)
        'Percentage': First sort the raw counts in one ring, then keep the pixels with 
                    raw counts fall in range (low%~high%) (for example 5%~95%).
        'Raw counts/Percentage': both raw counts and percentage are used, the overlap of 
                    masked pixels calculated with two method will be masked.
        
        picvector:  ssp.csr_matrix, a flattened pic vector with shape(xdimension*ydimension,1)
        picflat:    ndarray, a 1D flattened pic array
        tmatrix:    tmatrix for selfcorr mask calculation. Could be generated with Tmatrix.genTmatrixselfcorr()
                    it is a tmatrix calculated using non-splitting method.
        
        return    2d ndarray of boolean
        '''
        if self.maskselfcorrmethod=='Fast':
            re = self.selfcorrNormal(picvector, picflat, tmatrix)
        else:
            re = self.selfcorrMMM(picvector, picflat, tmatrix)
        # accumulate mode
        if self.maskselfcorraccuenable:
            self.selfcorrmaskhistory = self.selfcorrmaskhistory + re
            historymax = np.max(self.selfcorrmaskhistory)
            re = self.selfcorrmaskhistory > (historymax / 2.0)
        self.selfcorrmask = re
        return self.selfcorrmask
    
    def selfcorrMaskDiag(self):
        '''return a scipy.sparse matrix with main diagonal equal to the flattened mask array
        '''
        return ssp.spdiags(self.selfcorrmask.ravel(), [0], self.xydimension, self.xydimension)
    
    def selfcorrMMM(self, picvector, picflat, tmatrix):
        '''calculate the selfcorr mask of three modes other than 'Fast'
        return    2d ndarray of boolean
        '''
        picvectordiag = ssp.spdiags(picflat, [0], len(picflat), len(picflat)).tocsr()
        tmatrixpic = np.dot(tmatrix, picvectordiag).tocsr()
        maxint = np.zeros(tmatrix.shape[0])
        minint = np.zeros(tmatrix.shape[0])
        for i in range(tmatrix.shape[0]):
            ss = np.sort(tmatrixpic.data[tmatrixpic.indptr[i]:tmatrixpic.indptr[i+1]])
            if len(ss)>0:
                if self.maskselfcorrmethod == 'Raw counts/Percentage':
                    mid = float(ss[len(ss)/2])
                    xxhi = int(len(ss) * (1 - self.maskselfcorrPercentagehi))
                    xxlow = int(len(ss) * self.maskselfcorrPercentagelow)
                    midhi = mid * (1 + self.maskselfcorrCounthi)
                    midlow = mid * (1 - self.maskselfcorrCountlow)
                    if xxhi>0:
                        maxint[i] = min(ss[-xxhi], midhi)
                    else:
                        maxint[i] = midhi
                    if xxlow>0:
                        minint[i] = max(ss[xxlow], midlow)
                    else:
                        minint[i] = midlow
                elif self.maskselfcorrmethod == 'Percentage':
                    xxhi = max(int(len(ss) * (1 - self.maskselfcorrPercentagehi)),1)
                    xxlow = max(int(len(ss) * self.maskselfcorrPercentagelow), 0)
                    maxint[i] = ss[-xxhi]
                    minint[i] = ss[xxlow]
                elif self.maskselfcorrmethod == 'Raw counts':
                    mid = float(ss[len(ss)/2])
                    maxint[i] = mid * (1 + self.maskselfcorrcounthi)
                    minint[i] = mid * (1 - self.maskselfcorrcountlow)
        maxflat = np.dot(tmatrix.transpose().tocsr(),
                         ssp.csc_matrix(np.array([maxint]).transpose()))
        maxflat = np.array(maxflat.todense().transpose())[0]
        minflat = np.dot(tmatrix.transpose().tocsr(),
                         ssp.csc_matrix(np.array([minint]).transpose()))
        minflat = np.array(minflat.todense().transpose())[0]
        mask1 = np.logical_and((picflat < maxflat), (picflat > minflat))
        mask1 = mask1.reshape((self.ydimension, self.xdimension))
        return mask1
    
    def selfcorrNormal(self, picvector, picflat, tmatrix):
        '''calculate selfcorr mask of 'Fast' mode
        ''' 
        #intensity 1st run
        number = np.dot(tmatrix, self.onevector)
        intensity = np.dot(tmatrix, picvector)
        intensity = intensity / number
        #auto intensity filting
        intensityflat = np.dot(tmatrix.transpose().tocsr(),intensity.tocsc())
        intensityflat = np.array(intensityflat.todense().transpose())[0]
        mask1 = np.logical_and((picflat < (intensityflat * (1+self.maskselfcorrcounthi))),
                              (picflat > (intensityflat * (1-self.maskselfcorrcountlow))))
        #mask1 = ssp.spdiags(mask1, [0], self.xydimension, self.xydimension).tocsr()
        mask1 = mask1.reshape((self.ydimension, self.xdimension))
        return mask1
    
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
    
    def cropImage(self, pic):
        '''crop image if xrange/yrange is specified in config
        '''
        if self.xdimensionorg != self.xdimension:
            pic = pic[:, self.xrangelow:self.xrangehi]
        if self.ydimensionorg != self.ydimension:
            pic = pic[self.yrangelow:self.yrangehi, :]
        return pic
    
    def cropAndFlipImage(self, pic):
        '''flip image first then crop it.
        '''
        pic = self.flipImage(pic)
        pic = self.cropImage(pic)
        return pic

    def genTTHMatrix(self):
        ''' Calculate the two theta matrix which stores the two theta value of
        each pixel in the 2D-detector.
        
        Return:    2d ndarray of float, shape(ydimension,xdimension), store tth value of each pixel
        '''
        
        xr = (np.arange(self.xdimension, dtype=float) - self.xbeamcenter) * self.xpixelsize
        yr = (np.arange(self.ydimension, dtype=float) - self.ybeamcenter) * self.ypixelsize
        
        sinr = np.sin(self.rotation)
        cosr = np.cos(self.rotation)
        sint = np.sin(self.tilt)
        cost = np.cos(self.tilt)
        sourcexr = self.distance * sint * cosr
        sourceyr = self.distance * sint * sinr
        sourcezr = self.distance * cost
        
        dmatrix = np.zeros((self.ydimension, self.xdimension), dtype=float)
        dmatrix += ((xr - sourcexr) ** 2).reshape(1, self.xdimension)
        dmatrix += ((yr - sourceyr) ** 2).reshape(self.ydimension, 1)
        dmatrix += sourcezr ** 2
        dmatrix = np.sqrt(dmatrix)
        
        tthmatrix1 = np.zeros((self.ydimension,self.xdimension), dtype=float)
        tthmatrix1 += (-(xr - sourcexr) * sint * cosr).reshape(1, self.xdimension)
        tthmatrix1 += (-(yr - sourceyr) * sint * sinr).reshape(self.ydimension, 1)
        tthmatrix1 += sourcezr * cost
        tthmatrix = np.arccos(tthmatrix1 / dmatrix)
        return tthmatrix
