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
import scipy.ndimage.filters as snf
import scipy.ndimage.morphology as snm
from diffpy.srxplanar.srxplanarconfig import _configPropertyR


class Calculate(object):
    '''provide methods for integration, variance calculation and distance/Q matrix calculation etc.
    '''
    # define configuration properties that are forwarded to self.config
    xdimension = _configPropertyR('xdimension')
    ydimension = _configPropertyR('ydimension')
    xpixelsize = _configPropertyR('xpixelsize')
    ypixelsize = _configPropertyR('ypixelsize')
    xbeamcenter = _configPropertyR('xbeamcenter')
    ybeamcenter = _configPropertyR('ybeamcenter')
    rotation = _configPropertyR('rotation')
    tilt = _configPropertyR('tilt')
    distance = _configPropertyR('distance')
    wavelength = _configPropertyR('wavelength')
    integrationspace = _configPropertyR('integrationspace')
    qmax = _configPropertyR('qmax')
    qstep = _configPropertyR('qstep')
    tthmax = _configPropertyR('tthmax')
    tthstep = _configPropertyR('tthstep')
    tthmaxd = _configPropertyR('tthmaxd')
    tthstepd = _configPropertyR('tthstepd')
    tthorqstep = _configPropertyR('tthorqstep')
    tthorqmax = _configPropertyR('tthorqmax')
    uncertaintyenable = _configPropertyR('uncertaintyenable')
    sacorrectionenable = _configPropertyR('sacorrectionenable')
    polcorrectionenable = _configPropertyR('polcorrectionenable')
    polcorrectf = _configPropertyR('polcorrectf')
    selfcorrenable = _configPropertyR('selfcorrenable')


    def __init__(self, p):
        #create parameter proxy, so that parameters can be accessed by self.parametername in read-only mode
        self.config = p
        self.prepareCalculation()
        return
    
    def prepareCalculation(self):
        '''prepare data for calculation
        '''
        self.xydimension = self.xdimension * self.ydimension
        self.xr = (np.arange(self.xdimension, dtype=float) - self.xbeamcenter + 0.5) * self.xpixelsize
        self.yr = (np.arange(self.ydimension, dtype=float) - self.ybeamcenter + 0.5) * self.ypixelsize
        self.dmatrix = self.genDistanceMatrix()
        self.azimuthmatrix = np.arctan2(self.yr.reshape(self.ydimension,1), 
                                        self.xr.reshape(1,self.xdimension))
        #set tth or q grid
        if self.integrationspace == 'twotheta':
            self.tthorq = self.tth = np.arange(0.0, self.tthmax, self.tthstep)
            self.tthorqoutput = np.arange(0.0, self.tthmaxd, self.tthstepd)
            self.tthorqmatrix = self.genTTHMatrix()
        elif self.integrationspace == 'qspace':
            self.tthorq = self.q = np.arange(0.0, self.qmax, self.qstep)
            self.tthorqoutput = self.tthorq
            self.tthorqmatrix = self.genQMatrix()
        return
    
    def genIntegrationInds(self, mask=None):
        '''generate index used in integration
        picflat[ind[indlow[i]]:ind[indhigh[i]]] belong to one tth or q bin
        
        param mask: mask 2D array, same dimension as image, 1 for masked pixel
        
        return: self.ind, self.indlow, self.indhigh: see usage before
        '''
        if mask == None:
            mask = np.zeros((self.ydimension, self.xdimension), dtype=boolean)
        tthorqmatrix = self.tthorqmatrix
        tthorqmatrix[mask] = 1000.0
        tthorqmatrix = np.rint(tthorqmatrix / self.tthorqstep).astype(int)
        tthorqflat = tthorqmatrix.ravel()
        
        self.ind = np.argsort(tthorqflat)
        self.sortedtthorqindflat = tthorqflat[self.ind]
        sind = np.nonzero(np.diff(self.sortedtthorqindflat))[0] + 1
        self.indlow = np.concatenate([[0], sind])
        self.indhigh = np.concatenate([sind, [self.xdimension*self.ydimension]])
        return self.ind, self.indlow, self.indhigh
    
    def intensity(self, pic, picvar=None):
        '''2D to 1D image integration, intensity of pixels are binned and then take average,
        if self.selfcorrenable is True, then pixels whose intensity are too high/low will be dropped.
        
        self.ind, self.indlow, self.indhigh should be calculated first
        uncertainty will be calculated if self.unceratintyenable=True
        
        param pic: 2D array, array of raw counts, raw counts should be corrected
        param picvar: 2D array, variance of raw counts, if None but required in calculation, use zero matrix instead.
        
        retrun: 2d array, [tthorq, intensity, unceratinty] or [tthorq, intensity]
        '''
        import matplotlib.pyplot as plt 
        
        ind = self.ind
        indlow = self.indlow
        indhigh = self.indhigh
        if self.uncertaintyenable:
            picflat = pic.ravel()[ind]
            picvarflat = picvar.ravel()[ind] if picvar!=None else np.zeros_like(picflat)
            intensity = np.zeros_like(self.tthorqoutput)
            std = np.zeros_like(self.tthorqoutput)
            leng = len(intensity)
            for i in xrange(len(indlow)):
                dataind = self.sortedtthorqindflat[indlow[i]]
                if dataind<leng:
                    data = picflat[indlow[i]:indhigh[i]]
                    datavar = picvarflat[indlow[i]:indhigh[i]]
                    azimuth = azimuthflat[indlow[i]:indhigh[i]]
                    if self.selfcorrenable:
                        medianint = np.median(data)
                        ind1 = np.logical_and(medianint*0.2<data, data<medianint*5)
                        data = data[ind1]
                        datavar = datavar[ind1]
                    intensity[dataind] = np.mean(data)
                    std[dataind] = np.sqrt(np.mean(datavar)/len(datavar))
            rv = np.vstack([self.tthorqoutput, intensity, std])
            rv[np.isnan(rv)] = 0
        else:
            picflat = pic.ravel()[ind]
            intensity = np.zeros_like(self.tthorqoutput)
            std = np.zeros_like(self.tthorqoutput)
            leng = len(intensity)
            for i in xrange(len(indlow)):
                dataind = self.sortedtthorqindflat[indlow[i]]
                if dataind<leng:
                    data = picflat[indlow[i]:indhigh[i]]
                    if self.selfcorrenable:
                        medianint = np.median(data)
                        ind1 = np.logical_and(medianint*0.2<data, data<medianint*5)
                        data = data[ind1]
                    intensity[dataind] = np.mean(data)
            rv = np.vstack([self.tthorqoutput, intensity])
            rv[np.isnan(rv)] = 0    
        return rv
    
    def varianceLocal(self, pic):
        '''calculate the variance of raw counts of each pixel are calculated according to their 
        loacl variance.
        
        param picflat: 1d array, flattend image array
        
        return: 2d array, variance of each pixel
        '''
        
        picavg = snf.uniform_filter(pic, 5, mode='wrap')
        pics2 = (pic - picavg)**2
        pvar = snf.uniform_filter(pics2, 5, mode = 'wrap')       
        
        gain = pvar/pic
        gain[np.isnan(gain)] = 0
        gain[np.isinf(gain)] = 0
        gainmedian = np.median(gain)
        var = pic * gainmedian
        return var
    
    def genDistanceMatrix(self):
        '''Calculate the distance matrix
        
        return: 2d array, distance between source and each pixel'''
        sinr = np.sin(self.rotation)
        cosr = np.cos(self.rotation)
        sint = np.sin(self.tilt)
        cost = np.cos(self.tilt)
        sourcexr = self.distance * sint * cosr
        sourceyr = -self.distance * sint * sinr
        sourcezr = self.distance * cost
        
        dmatrix = np.zeros((self.ydimension, self.xdimension), dtype=float)
        dmatrix += ((self.xr - sourcexr) ** 2).reshape(1, self.xdimension)
        dmatrix += ((self.yr - sourceyr) ** 2).reshape(self.ydimension, 1)
        dmatrix += sourcezr ** 2
        self.dmatrix = np.sqrt(dmatrix)
        return self.dmatrix

    def genTTHMatrix(self):
        '''Calculate the diffraction angle matrix 
        
        return: 2d array, two theta angle of each pixel's center'''
    
        sinr = np.sin(self.rotation)
        cosr = np.cos(self.rotation)
        sint = np.sin(self.tilt)
        cost = np.cos(self.tilt)
        sourcexr = self.distance * sint * cosr
        sourceyr = -self.distance * sint * sinr
        sourcezr = self.distance * cost
        
        tthmatrix1 = np.zeros((self.ydimension,self.xdimension), dtype=float)
        tthmatrix1 += ((-self.xr + sourcexr) * sourcexr).reshape(1, self.xdimension)
        tthmatrix1 += ((-self.yr + sourceyr) * sourceyr).reshape(self.ydimension, 1)
        tthmatrix1 += sourcezr * sourcezr
        tthmatrix = np.arccos(tthmatrix1 / self.dmatrix / self.distance)
        return tthmatrix
    
    def genQMatrix(self):
        '''Calculate the q matrix 
        
        return: 2d array, q value of each pixel's center'''
        sinr = np.sin(self.rotation)
        cosr = np.cos(self.rotation)
        sint = np.sin(self.tilt)
        cost = np.cos(self.tilt)
        sourcexr = self.distance * sint * cosr
        sourceyr = -self.distance * sint * sinr
        sourcezr = self.distance * cost
        
        tthmatrix1 = np.zeros((self.ydimension,self.xdimension), dtype=float)
        tthmatrix1 += ((-self.xr + sourcexr) * sourcexr).reshape(1, self.xdimension)
        tthmatrix1 += ((-self.yr + sourceyr) * sourceyr).reshape(self.ydimension, 1)
        tthmatrix1 += sourcezr * sourcezr
        tthmatrix = np.arccos(tthmatrix1 / self.dmatrix / self.distance)
        Q = 4 * np.pi * np.sin(tthmatrix / 2.0) / self.wavelength
        return Q
    
    def genCorrectionMatrix(self):
        '''apply solid angle correction or polarization correction (powder) on tmatrix
        return:     2darray, correction to apply on the image
        '''
        rv = self._solidAngleCorrection() * self._polarizationCorrection()
        return rv

    def _solidAngleCorrection(self):
        '''generate correction matrix of soild angle correction for 2D flat detector. 
        return:    2darray, float, correction matrix'''
        if self.sacorrectionenable:
            sourcezr = self.distance * np.cos(self.tilt)
            correction = (self.dmatrix / sourcezr) 
        else:
            correction = np.ones((self.ydimension, self.xdimension))
        return correction
    
    def _polarizationCorrection(self):
        '''generate correction matrix of polarization correction for powder diffraction for 2D flat detector.
        require the polarization factor in configuration. 

        return:    2darray, float, correction matrix'''
        if self.polcorrectionenable:
            tthmatrix = self.tthorqmatrix if self.integrationspace == 'twotheta' else self.genTTHMatrix() 
            azimuthmatrix = self.azimuthmatrix
            p = 0.5 * (1 + (np.cos(tthmatrix))**2)
            p1 = 0.5 * self.polcorrectf * np.cos(2*azimuthmatrix) *  (np.sin(tthmatrix))**2
            p = 1.0 / (p-p1)
        else:
            p = np.ones((self.ydimension, self.xdimension))
        return p
