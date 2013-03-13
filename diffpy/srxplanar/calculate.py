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

'''module that provides the various integration method'''

import numpy as np
import scipy.sparse as ssp
import scipy.ndimage.filters as snf
import scipy.ndimage.morphology as snm
from diffpy.srxplanar.srxplanarconfig import _configPropertyR


class Calculate(object):
    '''provide calculate methods that do integration, transformations and
    calculation
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
    selfcorrenable = _configPropertyR('selfcorrenable')
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
        self.xydimension = self.xdimension * self.ydimension
        self.xr = (np.arange(self.xdimension, dtype=float) - self.xbeamcenter + 0.5) * self.xpixelsize
        self.yr = (np.arange(self.ydimension, dtype=float) - self.ybeamcenter + 0.5) * self.ypixelsize
        self.dmatrix = self.genDistanceMatrix()
        
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
        '''
        if mask == None:
            mask = np.zeros((self.ydimension, self.xdimension), dtype=boolean)
        else:
            mask = np.logical_not(mask)
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
        '''pic should be corrected
        '''
        ind = self.ind
        indlow = self.indlow
        indhigh = self.indhigh
        if self.uncertaintyenable:
            picflat = pic.ravel()[ind]
            picvarflat = picvar.ravel()[ind]
            intensity = np.zeros_like(self.tthorqoutput)
            std = np.zeros_like(self.tthorqoutput)
            leng = len(intensity)
            for i in xrange(len(indlow)):
                dataind = self.sortedtthorqindflat[indlow[i]]
                if dataind<leng:
                    data = picflat[indlow[i]:indhigh[i]]
                    datavar = picvarflat[indlow[i]:indhigh[i]]
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
        loacl variance.  (This one is fast, the result is similar to mode1 with spotty=True)
        picflat:        1d array, flattend pic array
        
        return:         2d array of variance of each pixel
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
        '''Calculate the distance matrix which stores the distance between the pixels and source.
        Return:    ndarray, shape(ydimension, xdimension), float'''
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
        """ Calculate the two theta matrix which stores the two theta value of each pixel.
        Return:    ndarray, shape(ydimension, xdimension), float"""
    
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
        """ Calculate the q matrix which stores the q value of each pixel.
        Return:    ndarray, shape(ydimension, xdimension), float"""
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
        '''apply solid angle correction or polarization correction on tmatrix
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
            tthmatrix = self.genTTHMatrix()
            a = (1.0 - self.polcorrectf) / (1.0 + self.polcorrectf)
            p = 1.0 / ((1.0 + a * (np.cos(tthmatrix)) ** 2) / (1.0 + a))
        else:
            p = np.ones((self.ydimension, self.xdimension))
        return p
