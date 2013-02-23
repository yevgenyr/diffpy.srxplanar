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

class Calculate(object):
    '''provide calculate methods that do integration, transformations and calculation
    '''
    def __init__(self, p):
        '''set the runtime parameter'''
        self.config = p
        self.configlist = ['xdimension', 
                   'ydimension', 
                   'xpixelsize',
                   'ypixelsize',
                   'xbeamcenter',
                   'ybeamcenter',
                   'integrationspace', 
                   'spotty',
                   'qmax', 
                   'qstep', 
                   'tthmax', 
                   'tthstep',
                   ]
        # FIXME PJ:  Here you are changing the class definition as a part of
        # instance initialization.  This is hard to understand and hard to maintain
        # in the future.  The class properties should be defined just once
        # in the class code.  The __init__ method should leave the class alone
        # (unless really justified which is not here).
        for optionname in self.configlist:
            if hasattr(self.config, 'xrd'+optionname):
                setattr(self.__class__, optionname, self.configProperty('xrd'+optionname))
            elif hasattr(self.config, optionname):
                setattr(self.__class__, optionname, self.configProperty(optionname))
        self.prepareCalculation()
        return
    
    def configProperty(self, nm):
        '''helper function of property delegation
        '''
        rv = property(fget = lambda self: getattr(self.config, nm))
        return rv
    
    def prepareCalculation(self):
        self.xydimension = self.xdimension * self.ydimension
        self.data = np.ones(self.xydimension)
        self.onevector = ssp.csc_matrix(self.data.reshape(len(self.data),1))        
        if self.integrationspace == 'twotheta':
            self.tthorq = self.tth = np.arange(0.0, self.tthmax, self.tthstep)
        elif self.integrationspace == 'qspace':
            self.tthorq = self.q = np.arange(0.0, self.qmax, self.qstep)
        xr = (np.arange(self.xdimension, dtype=float) - self.xbeamcenter) * self.xpixelsize
        yr = (np.arange(self.ydimension, dtype=float) - self.ybeamcenter) * self.ypixelsize
        azimuthmatrix = np.arctan2(yr.reshape(self.ydimension,1), 
                               xr.reshape(1,self.xdimension))
        self.azimuthflat = azimuthmatrix.ravel()
        return
    

    # FIXME PJ: functions that do not use self should be defined as
    # normal standalone functions.  Another option is to define them as
    # classmethods if you feel they really belong to the class namespace.
    def intensity(self, picvector, tmatrix, number):
        '''calculate the integrated intensity with transformation matrix tmatrix.
        picvector:         a vector of flattened image shape=(n,1)
        tmatrix, number:   transfromation matrix and its corresponding row summation array
        
        return:            integrated intensity array (1d numpy array)'''
        intensity = np.dot(tmatrix, picvector)
        intensity = np.array(intensity.todense().transpose())[0]
        ind = np.nonzero(number)
        intensity[ind] = intensity[ind] / number[ind]
        return intensity
    
    def varianceAll(self, picvector, tmatrix, number, intensityarray):
        '''calculate the variance of intensity, according to the all pixels in one ring
        picvector:         a vector of flattened image shape=(n,1)
        tmatrix, number:   transfromation matrix and its corresponding row summation array
        intensityarray:    1d array, integrated intensity 
        
        return:
        var:               variace of intensity
        varflat:           variance of each pixel in a flattened array'''
        
        intensity = ssp.csc_matrix(intensityarray.reshape(len(intensityarray),1))
        intensityflat = np.dot(tmatrix.transpose().tocsr(),intensity)
        var = picvector - intensityflat
        var.data = var.data ** 2
        var = np.dot(tmatrix, var)
        var = np.array(var.todense().transpose())[0]
        ind = np.nonzero(var)
        var[ind] = var[ind] / number[ind]
        varvector = ssp.csc_matrix(var.reshape(len(var),1))
        varflat = np.dot(tmatrix.transpose().tocsr(),varvector)
        varflat = np.array(varflat.todense().transpose())[0]
        return var, varflat
    
    def standardDeviation(self, variance, number):
        '''calculate the standard deviation of intensity, take std by sqrt of the variance, 
        this method only work with normal mode
        variance, number:     variance of intensity and row summation array of 
                              corresponding transformation matrix
    
        return:               1d ndarray, standard deviation of intensity'''
        
        ind = (number - 1) > 0
        std = np.zeros_like(variance)
        std[ind] = np.sqrt(variance[ind] / (number[ind] - 1))
        ind = np.logical_not(ind)
        std[ind] = np.sqrt(variance[ind])
        tempvar = (np.isfinite(std) == False)
        std[tempvar] = 0.0
        return std
    
    def vcmatrixLocal(self, picflat, tmatrixvar, tmatrix, number):
        '''calculate the VC matrix of integrated intensity. the variance of raw counts of each 
        pixel are calculated according to their local variance.
        picflat:        1d array, flattend pic array
        tmatrixvar:     ssp.csr_matrix, tmatrix used for error calculation
        tamtrix:        ssp.csr_matrix, tmatrix used for intensity calculation
        number:         row summation array of tmatrix
        
        return:         ssp.csc_matrix, vcmatrix of xrd intensity
        '''
        varflat = self.varianceLocal(picflat, tmatrixvar)    
        vardiag = ssp.spdiags(varflat, [0], len(varflat), len(varflat)).tocsc()
        ind = np.nonzero(number)
        number1 = np.zeros_like(number)
        number1[ind] = (1 / number[ind])
        numberdiag = ssp.spdiags(number1, [0], len(number), len(number)).tocsr()
        tmatrix1 = np.dot(numberdiag, tmatrix.tocsc()).tocsr()
        vcmatrix = np.dot(tmatrix1, vardiag)
        vcmatrix = np.dot(vcmatrix, tmatrix1.transpose()).tocsc()
        return vcmatrix
    
    def varianceLocal_mode1(self, picflat, tmatrix):
        '''calculate the variance of raw counts of each pixel are calculated according to their 
        loacl variance. (This one is slow)
        picflat:        1d array, flattend pic array
        tamtrix:        ssp.csr_matrix, tmatrix used for intensity calculation
        *spotty:        If self.spotty is True, the raw counts variance of spotty pixels will be calculated 
                        from the average detector gain determined by rest pixels
        
        return:         flattened array of variance of each pixel
        '''
        picvectordiag = ssp.spdiags(picflat, [0], len(picflat), len(picflat)).tocsr()
        n = self.xydimension
        tmatrixdotint = np.dot(tmatrix, picvectordiag).tocsr()
        data = np.ones(n)
        col = np.arange(n)
        azimuthind = np.argsort(self.azimuthflat)
        sortmatrix = ssp.csr_matrix((data, (azimuthind, col)), shape = (n,n))
        transformed = np.dot(tmatrixdotint, sortmatrix).tocsr()
        transformed.sort_indices()

        tn = 20
        for i in range(tmatrixdotint.shape[0]):
            a = transformed[i]
            if a.data.shape[0]>0:
                p = a.data
                pavg = snf.uniform_filter1d(a.data, tn, mode='wrap')
                ps2 = (p - pavg)**2
                pvar = snf.uniform_filter1d(ps2, tn, mode = 'wrap')                
                transformed.data[transformed.indptr[i]:transformed.indptr[i+1]] = pvar

        result = np.dot(transformed.tocsr(), sortmatrix.transpose().tocsc())
        onediag = ssp.csr_matrix(np.ones(result.shape[0]))
        varflat = np.array(np.dot(onediag, result).todense()[0])[0]
        
        if self.spotty:
            gainflat = varflat/picflat
            gainflat[np.isnan(gainflat)] = 0
            gainflat[np.isinf(gainflat)] = 0
            gainavg = np.average(gainflat)
            ind = gainflat > (gainavg) * 5.0
            #gainavg = np.average(gainflat[np.logical_not(ind)])
            #ind = gainflat > (gainavg)
            gain = np.average(gainflat[np.logical_not(ind)])
            varflat[ind] = picflat[ind] * gain
        #if hasattr(self, 'savevar'):
        #    np.save('std.npy', np.sqrt(varflat.reshape(self.ydimension, self.xdimension)))
        return varflat
    
    def varianceLocal(self, picflat, tmatrix):
        '''calculate the variance of raw counts of each pixel are calculated according to their 
        loacl variance.  (This one is fast, the result is similar to mode1 with spotty=True)
        picflat:        1d array, flattend pic array
        tamtrix:        ssp.csr_matrix, tmatrix used for intensity calculation, not used in this
                        method.
        
        return:         flattened array of variance of each pixel
        '''
        pic = picflat.reshape((self.ydimension, self.xdimension))
        k = np.ones((5,5))
        k[0,0] = k[4,4] = k[0,4] = k[4,0] = 0
        picavg = snf.uniform_filter1d(pic, 5, mode='wrap')
        pics2 = (pic - picavg)**2
        pvar = snf.uniform_filter1d(pics2, 5, mode = 'wrap')
        varflat = pvar.ravel()       
        
        gainflat = varflat/picflat
        gainflat[np.isnan(gainflat)] = 0
        gainflat[np.isinf(gainflat)] = 0
        gainmedian = np.median(gainflat)
        ind = gainflat > (gainmedian) * 100.0
        gain = np.average(gainflat[np.logical_not(ind)])
        varflat[ind] = picflat[ind] * gain
        #if hasattr(self, 'savevar'):
        #    np.save('std.npy', np.sqrt(varflat.reshape(self.ydimension, self.xdimension)))
        return varflat
    
    def sqToTheta(self, azimuthmaprange, tm, anglen):
        '''calculate theta value of q sequence
        '''
        q = tm.q[azimuthmaprange[0]:azimuthmaprange[1]]
        azimuth = tm.azimuth
        tth = 2 * np.arcsin(q * self.xraywavelength / 4 /np.pi)
        theta1 = (np.pi - tth) / 2
        theta2 = np.roll(theta1, -anglen)
        theta1 = theta1.reshape(len(theta1),1)
        theta2 = theta2.reshape(len(theta2),1)
        temp = np.sin(theta1) * np.sin(theta2) * np.cos(azimuth) + np.cos(theta1) * np.cos(theta2)
        temp[temp>1] = 1
        temp[temp<-1] = -1
        theta = np.arccos(temp)
        #np.save('theta.npy', theta)
        return theta
        
    def number(self, tmatrix):
        '''calculate the row summation array used in integration
        '''
        number = np.dot(tmatrix, self.onevector)
        number = np.array(number.todense().transpose())[0]
        return number
