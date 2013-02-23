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

class Tmatrix(object):
    '''generate transformation matrix used in 2d image integration'''
    def __init__(self, p):
        '''pass configs from pdfliveconfig
        '''
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
                           'tthstep',
                           'tthmax',
                           'azimuthstep',
                           'qmax',
                           'qstep',
                           'tthvarstep',
                           'qvarstep',
                           'tthorqmax',
                           'tthorqstep',
                           'xrduncertaintyenable',
                           'sacorrectionenable',
                           'polcorrectionenable',
                           'polcorrectf',
                           'integrationspace',
                           'method',
                           'regulartmatrixenable',
                           'azimuthmapenable',
                           'azimuthmapmode',
                           'azimuthmaprangelow',
                           'azimuthmaprangehi',
                           'azimuthmaprangelown',
                           'azimuthmaprangehin']
        
        for optionname in self.configlist:
            if hasattr(self.config, 'xrd'+optionname):
                setattr(self.__class__, optionname, self.configProperty('xrd'+optionname))
            elif hasattr(self.config, optionname):
                setattr(self.__class__, optionname, self.configProperty(optionname))
        self.prepareCalculation()
        return

    def configProperty(self, nm):
        '''helper function of options delegation'''
        rv = property(fget = lambda self: getattr(self.config, nm))
        return rv
    
    def updateConfig(self,p):
        '''call perparecalculation to update the metadata when the configruation changes
        '''
        self.prepareCalculation()
        return
    
    def prepareCalculation(self):
        """prepare the data for integration
        xr,yr:         relative position of each pixel
        dmatrix:       distance matrix
        tthorq:        tth(q) value of tth(q) bin sequence
        tthmatrix:     tth/tthstep value of each pixel
        qmatrix:       q/qstep of each pixel 
        tthorqflat:    flattened tth/qmatrix
        tthorqvar:     tth(q) value of tth(q) bin sequence used in error calculation
        tthvarmatrix:  tthmatrix used in error calculation
        qvarmatrix:    qmatrix used in error calculation
        tthorqvarflat: flattened tthvarmatrix(qvarmatrix)
        azimuth:       azimuth value of azimuth bin sequence
        azimuthmatrix: store the azimuth value of each pixel
        azimuthflat:   flattened azimuthmatrix 
        """
        self.xr = (np.arange(self.xdimension, dtype=float) - self.xbeamcenter) * self.xpixelsize
        self.yr = (np.arange(self.ydimension, dtype=float) - self.ybeamcenter) * self.ypixelsize
        self.dmatrix = self.genDistanceMatrix(self.xr, self.yr)
        
        if self.integrationspace == 'twotheta':
            self.tthorq = self.tth = np.arange(0.0, self.tthmax, self.tthstep)
            tthmatrix = self.genTTHMatrix(self.xr, self.yr, self.dmatrix)
            self.tthmatrix = (tthmatrix / self.tthstep).astype(int)
            self.tthorqflat = self.tthflat = self.tthmatrix.ravel()
            if self.xrduncertaintyenable:
                self.tthorqvar = self.tthvar = np.arange(0.0, self.tthmax, self.tthvarstep)
                self.tthvarmatrix = (tthmatrix / self.tthvarstep).astype(int)
                self.tthorqvarflat = self.tthvarflat = self.tthvarmatrix.ravel()
            
        elif self.integrationspace == 'qspace':
            self.tthorq = self.q = np.arange(0.0, self.qmax, self.qstep)
            qmatrix = self.genQmatrix(self.xr, self.yr, self.dmatrix, self.wavelength)
            self.qmatrix = (qmatrix / self.qstep).astype(int)
            self.tthorqflat = self.qflat = self.qmatrix.ravel()
            if self.xrduncertaintyenable:
                self.tthorqvar = self.qvar = np.arange(0.0, self.qmax, self.qvarstep)
                self.qvarmatrix = (qmatrix / self.qvarstep).astype(int)
                self.tthorqvarflat = self.qvarflat = self.qvarmatrix.ravel()
            
        self.azimuth = np.arange(0.0, np.radians(360), self.azimuthstep)
        self.azimuthmatrix = np.arctan2(self.yr.reshape(self.ydimension,1),
                                    self.xr.reshape(1,self.xdimension))
        self.azimuthflat = self.azimuthmatrix.ravel()
        return None
    
    def clean(self):
        '''release memory
        '''
        for data in ['xr', 'yr', 'dmatrix', 'tthorq', 'tthmatrix',
                     'tthorqflat', 'tthorqvar', 'tthvarmatrix', 'tthorqvarflat',
                     'azimuth', 'azimuthmatrix', 'azimuthflat',
                     'xrflat', 'yrflat', 'col0', 'yswapind', 'xswapind', 'swapind',
                     'tth1flat', 'tth4flat', 'q1flat', 'q4flat']:
            if hasattr(self, data):
                delattr(self, data)
        return
    
    def genTmatrix(self):
        '''General method to calculate Tmatrix used in intensity calculation 
        generate the tramsformation matrix according to the integration space and method, then apply 
        basic correction (solid angle and polarization) if enabled in config
        
        return:    ssp.csr_matrix, transformation matrix
        '''
        if self.method == 'normal':
            re = self.genTmatrixNormal(self.tthorq, self.tthorqflat)
        elif self.method == 'split':
            re = self.genTmatrixSplit()
        re = self.tmatrix = self.applyCorrection(re)
        return re
    
    def genTmatrixvar(self):
        '''calculate Tmatrix used in error calculation. The method is generally same as genTmatrix, 
        but use tthvarstep/qvarstep instead of tthstep/qstep.
        
        return:    ssp.csr_matrix, transformation matrix
        '''
        if not self.xrduncertaintyenable:
            print 'Mode change!'
            self.prepareCalculation()
        re = self.genTmatrixNormal(self.tthorqvar, self.tthorqvarflat)
        #re = self.tmatrixvar = self.applyCorrection(re)
        self.tmatrixvar = re
        return re
    
    def genTmatrixselfcorr(self):
        '''calculate the tmatrix used in self correction calculation. This tmatrix use non-splitting method.
        
        return:    ssp.csr_matrix, transformation matrix
        '''
        re = self.genTmatrixNormal(self.tthorq, self.tthorqflat)
        #re = self.tmatrixselfcorr = self.applyCorrection(re)
        self.tmatrixselfcorr = re
        return re
    
    def genDistanceMatrix(self, xr, yr):
        '''Calculate the distance matrix which stores the distance between the pixels and source.
        xr, yr:    ndarray, float, relative position array of input 2D pixel array
        
        Return:    ndarray, shape(ydimension, xdimension), float'''
        sinr = np.sin(self.rotation)
        cosr = np.cos(self.rotation)
        sint = np.sin(self.tilt)
        cost = np.cos(self.tilt)
        sourcexr = self.distance * sint * cosr
        sourceyr = self.distance * sint * sinr
        sourcezr = self.distance * cost
        ydimension = len(yr)
        xdimension = len(xr)
        
        dmatrix = np.zeros((ydimension, xdimension), dtype=float)
        dmatrix += ((xr - sourcexr) ** 2).reshape(1, xdimension)
        dmatrix += ((yr - sourceyr) ** 2).reshape(ydimension, 1)
        dmatrix += sourcezr ** 2
        dmatrix = np.sqrt(dmatrix)
        return dmatrix

    def genTTHMatrix(self, xr, yr, dmatrix):
        """ Calculate the two theta matrix which stores the two theta value of each pixel.
        xr, yr:    ndarray, float, relative position array of input 2D pixel array
        dmatrix:   2d ndarray, float, distance matrix of pixels
        
        Return:    ndarray, shape(ydimension, xdimension), float"""
    
        sinr = np.sin(self.rotation)
        cosr = np.cos(self.rotation)
        sint = np.sin(self.tilt)
        cost = np.cos(self.tilt)
        sourcexr = self.distance * sint * cosr
        sourceyr = self.distance * sint * sinr
        sourcezr = self.distance * cost
        ydimension = len(yr)
        xdimension = len(xr)
        
        tthmatrix1 = np.zeros((ydimension,xdimension), dtype=float)
        tthmatrix1 += (-(xr - sourcexr) * sint * cosr).reshape(1, xdimension)
        tthmatrix1 += (-(yr - sourceyr) * sint * sinr).reshape(ydimension, 1)
        tthmatrix1 += sourcezr * cost
        tthmatrix = np.arccos(tthmatrix1 / dmatrix)
        return tthmatrix
    
    def genQmatrix(self, xr, yr, dmatrix, xraywavelength):
        """ Calculate the q matrix which stores the q value of each pixel.
        xr, yr:    ndarray, float, relative position array of input 2D pixel array
        dmatrix:   2d ndarray, float, distance matrix of pixels
        xraywavelength: float, wavelength of xray
        
        Return:    ndarray, shape(ydimension, xdimension), float"""
        sinr = np.sin(self.rotation)
        cosr = np.cos(self.rotation)
        sint = np.sin(self.tilt)
        cost = np.cos(self.tilt)
        sourcexr = self.distance * sint * cosr
        sourceyr = self.distance * sint * sinr
        sourcezr = self.distance * cost
        ydimension = len(yr)
        xdimension = len(xr)
        
        tthmatrix1 = np.zeros((ydimension,xdimension), dtype=float)
        tthmatrix1 += (-(xr - sourcexr) * sint * cosr).reshape(1, xdimension)
        tthmatrix1 += (-(yr - sourceyr) * sint * sinr).reshape(ydimension, 1)
        tthmatrix1 += sourcezr * cost
        tthmatrix = np.arccos(tthmatrix1 / dmatrix)
        Q = 4 * np.pi * np.sin(tthmatrix / 2.0) / xraywavelength
        return Q

    def applyCorrection(self, tmatrix):
        '''apply solid angle correction or polarization correction on tmatrix
        tmatrix:    ssp.csr_matrix, float, tmatrix that to be corrected
        
        return:     ssp.csr_matrix, float, corrected tmatrix
        '''
        if self.sacorrectionenable:           
            corr = self._solidAngleCorrection()
            tmatrix = np.dot(tmatrix, corr).tocsr()
        if self.polcorrectionenable:
            corr = self._polarizationCorrection()
            tmatrix = np.dot(tmatrix, corr).tocsr()
        return tmatrix

    def _solidAngleCorrection(self):
        """generate correction matrix of soild angle correction for 2D flat detector. 
        Tmatrix(corrected) = Tmatrix dot CorrectionMatrix

        return:    ssp.csc_matrix, shape(xydimension ,xydimension), float, correction matrix"""
        sourcezr = self.distance * np.cos(self.tilt)
        correction = (self.dmatrix / sourcezr) ** 3
        lens = correction.shape[0] * correction.shape[1] 
        correction = ssp.spdiags(correction.ravel(), [0], lens, lens).tocsc()
        return correction
    
    def _polarizationCorrection(self):
        """generate correction matrix of polarization correction for powder diffraction for 2D flat detector.
        require the polarization factor in configuration. 
        Tmatrix(corrected) = Tmatrix dot CorrectionMatrix

        return:    ssp.csc_matrix, shape(xydimension ,xydimension), float, correction matrix"""
        tthmatrix = self.genTTHMatrix(self.xr, self.yr, self.dmatrix)
        a = (1.0 - self.polcorrectf) / (1.0 + self.polcorrectf)
        p = 1.0 / ((1.0 + a * (np.cos(tthmatrix)) ** 2) / (1.0 + a))
        lens = p.shape[0] * p.shape[1] 
        correction = ssp.spdiags(p.ravel(), [0], lens, lens).tocsc()
        return correction    
    
    def genTmatrixNormal(self, tthorq, tthorqflat):
        '''calculate the tramsformation matrix of normal method
        tthorq:     tth(q) value of tth(q) bin sequence
        tthorqflat: flattened tth/qmatrix which store the tth(q)/tthstep(qstep)
                    
        returm:     ssp.csr_matrix, float, all entries equal to 1 (non-splitting method)
        '''
        data = np.zeros(self.xdimension * self.ydimension) + 1
        col = np.arange((self.xdimension * self.ydimension),dtype = int)
        row = tthorqflat
        tmatrix = ssp.csr_matrix((data, (row,col)), shape = (len(tthorq), (self.xdimension * self.ydimension)))
        return tmatrix
    
    def genTmatrixSplit(self):
        '''calculate the tramsformation matrix of pixel splitting method
        require integration space in configration 
                    
        return:     ssp.csr_matrix, float, entries value is the coverage percentage of bin and pixel (between 0 and 1)

        note:            
        xrsp, yrsp:     x,y coordinate of pixels' corner
        dmatrixsp:      distance of pixels' corner and diffraction postion
        tthmatrixsp:    twotheta matrix of pixels' corner
        ld, lu, rd, ru: twotheta value of left down, left up, right down, right up corner
        azimuthind:     azimuth index of each pixel
        r1, r2, r3, r4: twotheta value, the first, second, third, fourth corners of each pixel arranged by their twotheta value
        x1, x2, x3, x4: index of twotheta value, the first, second, third, fourth corners
        opind:          operation index of different types of coverage
        sincos2:        1/(sin azimuth * cos azimuth)
        sincos:         max(1/ sin azimuth , 1/ cos azimuth)
        q***:           same as *** but in q-space
        '''
        xrsp = (np.arange(self.xdimension + 1, dtype=float) -
                self.xbeamcenter - 0.5) * self.xpixelsize
        yrsp = (np.arange(self.ydimension + 1, dtype=float) -
                self.ybeamcenter - 0.5) * self.ypixelsize
        dmatrixsp = self.genDistanceMatrix(xrsp, yrsp)
        tthmatrixsp = self.genTTHMatrix(xrsp, yrsp, dmatrixsp)
        ld, lu, rd, ru = self._genldlurdru(tthmatrixsp)
        azimuthind = self._genAzimuthInd18(self.azimuthmatrix)
        r1, r2, r3, r4 = self._genr1r2r3r4(azimuthind, ld, lu, rd, ru)
        sincos2, sincos = self._genSinCosMatrix(self.azimuthmatrix, azimuthind)

        if self.integrationspace=='twotheta':
            x1, x2, x3, x4 = self._genx1x2x3x4(r1, r2, r3, r4, self.tthstep)
            opind = self._genOperationInd(x1, x2, x3, x4)  
            tmatrixsp = self._genTmatrixspOP(opind, x1, x2, x3, x4, r1, r2, r3, r4, sincos2, sincos)
            
        if self.integrationspace=='qspace':
            qmatrixsp = self.genQmatrix(xrsp, yrsp, dmatrixsp, self.wavelength)
            qld, qlu, qrd, qru = self._genldlurdru(qmatrixsp)
            qazimuthind = self._genAzimuthInd18(self.azimuthmatrix)
            qr1, qr2, qr3, qr4 = self._genr1r2r3r4(qazimuthind, qld, qlu, qrd, qru)
            x1, x2, x3, x4 = self._genx1x2x3x4(qr1, qr2, qr3, qr4, self.qstep)
            opind = self._genOperationInd(x1, x2, x3, x4)
            tmatrixsp = self._genTmatrixspOP(opind, x1, x2, x3, x4, r1, r2, r3, r4, sincos2, sincos)
        if self.regulartmatrixenable:
            tmatrixsp = self._regularTmatrix(tmatrixsp)
        return tmatrixsp
    
    def _regularTmatrix(self, tmatrix):
        '''Normalize tmatrix so that no intensity losing in integration
        '''
        temp = tmatrix.tocsc()
        ones = ssp.csr_matrix(np.ones(temp.shape[0]))
        sumtm = np.dot(ones, tmatrix).transpose().tocsc()
        sumtm.data[sumtm.data == 0] = 1.0
        sumtm.data = 1 / sumtm.data
        sumtmdiag = ssp.spdiags(sumtm.data, [0], len(sumtm.data), len(sumtm.data))
        rv = np.dot(tmatrix, sumtmdiag)
        return rv.tocsr()
    
    def _genldlurdru(self, matrixsp):
        '''calculate twotheta or q value of left down, left up, right down, right up corner
        matrixsp: 2d ndarray, twotheta or q matrix of pixels' corner

        return:
        ld,lu,rd,ru: 2d ndarray, twotheta or q matrix of each pixels' 4 corners'''
        
        xdimension = matrixsp.shape[1]
        ydimension = matrixsp.shape[0]
        ld = matrixsp[0:ydimension - 1, 0:xdimension - 1]
        lu = matrixsp[1:ydimension, 0:xdimension - 1]
        rd = matrixsp[0:ydimension - 1, 1:xdimension]
        ru = matrixsp[1:ydimension, 1:xdimension]
        return ld,lu,rd,ru
    
    def _genAzimuthInd18(self, azimuthmatrix):
        '''calcualte azimuth index of all pixels
        azimuthmatrix: 2d ndarray, azimuth value of each pixel
        
        return:
        ind[0~7]: 2d boolean array
        ind[x]: True for x*pi/4 < azimuth < (x+1)*pi/4
        (for pi<azimuth<2pi -> -pi<azimuth<0) 
        '''
        ind = [[]] * 8
        ind00 = azimuthmatrix >= 0
        ind34 = azimuthmatrix >= (np.pi / 4 * 3)
        ind14 = azimuthmatrix >= (np.pi / 4)
        ind12 = azimuthmatrix >= (np.pi / 2)
        ind[0] = np.logical_and(ind00, np.logical_not(ind14))
        ind[1] = np.logical_and(ind14, np.logical_not(ind12))
        ind[2] = np.logical_and(ind12, np.logical_not(ind34))
        ind[3] = ind34
        ind00 = azimuthmatrix < 0
        ind34 = azimuthmatrix <= -(np.pi / 4 * 3)
        ind14 = azimuthmatrix <= -(np.pi / 4)
        ind12 = azimuthmatrix <= -(np.pi / 2)
        ind[4] = ind34
        ind[5] = np.logical_and(ind12, np.logical_not(ind34))
        ind[6] = np.logical_and(ind14, np.logical_not(ind12))
        ind[7] = np.logical_and(ind00, np.logical_not(ind14))
        return ind
    
    def _genr1r2r3r4(self, ind, ld, lu, rd, ru):
        '''calculate the 1st, 2nd, 3rd, 4th corner of each pixel according to the corners' twotheta or q value
        ind: azimuth ind generated by genazimuthind18
        ld,lu,rd,ru: twotheta or q value of each pixels' corner
        
        return:
        r1, r2, r3, r4: 2d array, twotheta or q value of 1st, 2nd, 3rd, 4th corner of each pixel
        '''
        ind1 = np.logical_or(ind[0], ind[1])
        ind2 = np.logical_or(ind[2], ind[3])
        ind3 = np.logical_or(ind[4], ind[5])
        ind4 = np.logical_or(ind[6], ind[7])
        r1 = ind1 * ld + ind2 * rd + ind3 * ru + ind4 * lu
        r4 = ind1 * ru + ind2 * lu + ind3 * ld + ind4 * rd
        ind1 = np.logical_or(ind[0], ind[5])
        ind2 = np.logical_or(ind[1], ind[4])
        ind3 = np.logical_or(ind[2], ind[7])
        ind4 = np.logical_or(ind[3], ind[6])
        r2 = ind1 * lu + ind2 * rd + ind3 * ld + ind4 * ru
        r3 = ind1 * rd + ind2 * lu + ind3 * ru + ind4 * ld
        
        tempind = np.nonzero((r2-r1)<0)
        temp = r1[tempind]
        r1[tempind] = r2[tempind]
        r2[tempind] = temp
        tempind = np.nonzero((r3-r2)<0)
        temp = r2[tempind]
        r2[tempind] = r3[tempind]
        r3[tempind] = temp
        tempind = np.nonzero((r4-r3)<0)
        temp = r3[tempind]
        r3[tempind] = r4[tempind]
        r4[tempind] = temp
        return r1,r2,r3,r4

    def _genx1x2x3x4(self, r1, r2, r3, r4, step):
        '''calculate the twotheta or q index of 1st, 2nd, 3rd, 4th corner of each pixel
        r1, r2, r3, r4: 2d array of twotheta or q value of 1st~4th corner
        step: twotheta or q step
        
        return:
        x1, x2, x3, x4: 2d array, twotheta or q index of 1st~4th corner of each pixel
        '''
        x1 = np.floor(r1 / step).astype(int)
        x2 = np.floor(r2 / step).astype(int)
        x3 = np.floor(r3 / step).astype(int)
        x4 = np.floor(r4 / step).astype(int)
        return x1,x2,x3,x4
    
    def _genOperationInd(self, x1, x2, x3, x4):
        '''generate the operation index of each pixel according to the relationship of x1,x2,x3,x4
        x1, x2, x3, x4: two theta or q index
        
        return:
        ind[0~7]: operation index
        [0]: x1=x2=x3=x4
        [1]: x1=x2=x3<x4
        [2]: x1=x2<x3=x4
        [3]: x1=x2<x3<X4
        [4]: x1<x2=x3=x4
        [5]: x1<x2=x3<x4
        [6]: x1<x2<x3=x4
        [7]: x1<x2<x3<x4
        '''
        ind12 = (x1 == x2)
        ind23 = (x2 == x3)
        ind34 = (x3 == x4)
        nind12 = np.logical_not(ind12)
        nind23 = np.logical_not(ind23)
        nind34 = np.logical_not(ind34)
        ind = [[]] * 8
        ind[0] = np.logical_and(np.logical_and(ind12, ind23), ind34)
        ind[1] = np.logical_and(np.logical_and(ind12, ind23), nind34)
        ind[2] = np.logical_and(np.logical_and(ind12, nind23), ind34)
        ind[3] = np.logical_and(np.logical_and(ind12, nind23), nind34)
        ind[4] = np.logical_and(np.logical_and(nind12, ind23), ind34)
        ind[5] = np.logical_and(np.logical_and(nind12, ind23), nind34)
        ind[6] = np.logical_and(np.logical_and(nind12, nind23), ind34)
        ind[7] = np.logical_and(np.logical_and(nind12, nind23), nind34)
        return ind
    
    def _genSinCosMatrix(self, azimuthmatrix, azimuthind):
        '''calculate the auxiliary sincos matrix
        azimuthmatrix: azimuth value of each pixel
        ind: azimuth ind generated by genazimuthind18
        
        return:
        sincos2: 1 / (2 * sin azimuth * cos azimuth)
        sincos: 1 / max(sin azimuth, cos azimuth)
        '''
        sincos2 = np.sin(azimuthmatrix) * np.cos(azimuthmatrix) * 2
        sincos2 = np.abs(1 / sincos2)
        sin1 = np.abs(np.sin(azimuthmatrix))
        cos1 = np.abs(np.cos(azimuthmatrix))
        sinind = sin1 > cos1
        cosind = np.logical_not(sinind)
        sincos = np.zeros_like(sincos2, dtype = float)
        sincos[sinind] = sin1[sinind]
        sincos[cosind] = cos1[cosind]
        sincos = 1 / sincos
        return sincos2, sincos
        
    def _genTmatrixspOP(self, opind, x11, x22, x33, x44, r11, r22, r33, r44, sincos22, sincos11):
        '''calculate the transformation matrix, modified consequence
        return:     ssp.csr_matrix, float, transformation matrix
        '''
        data0 = np.zeros(self.xdimension * self.ydimension) + 1
        col0 = np.arange((self.xdimension * self.ydimension),dtype = int)
        x1 = x11.ravel()
        x2 = x22.ravel()
        x3 = x33.ravel()
        x4 = x44.ravel()
        r1 = r11.ravel()
        r2 = r22.ravel()
        r3 = r33.ravel()
        r4 = r44.ravel()
        sincos2 = sincos22.ravel()
        sincos = sincos11.ravel()
        if self.integrationspace=='twotheta':
            r = self.tth
            maxi = int(np.ceil(self.xpixelsize / (self.distance * self.tthstep)))
        elif self.integrationspace=='qspace':
            r = 2 * np.arcsin(self.q * self.wavelength / 4 / np.pi)
            maxi = int(np.ceil(self.xpixelsize * 2 * np.pi/ (self.distance * self.wavelength * self.qstep)))
        tmatrixshape = (len(r), self.xdimension * self.ydimension)
        #cos matrix
        sourcezr = self.distance * np.cos(self.tilt)
        cosflat = (1 / (sourcezr/self.dmatrix)).ravel()
        dispix = (self.dmatrix / self.xpixelsize).ravel()
        corr = cosflat * dispix
        corr2 = corr ** 2
        
        #1 x1=x2=x3=x4
        ind = opind[0].ravel()
        data = data0[ind]
        col = col0[ind]
        row = x1[ind]
        tmatrix = ssp.csr_matrix((data, (row, col)), shape = tmatrixshape)
        
        #2 x1=x2=x3<x4
        ind = opind[1].ravel()
        col = col0[ind]
        temp = (r4[ind] - r[x4][ind]) ** 2 * sincos2[ind] * corr2[ind] 
        tmatrix = tmatrix + ssp.csr_matrix(((1-temp), (x1[ind], col)), shape = tmatrixshape)
        tmatrix = tmatrix + ssp.csr_matrix((temp, (x4[ind], col)), shape = tmatrixshape)
        
        #5 x1<x2=x3=x4
        ind = opind[4].ravel()
        col = col0[ind]
        temp = (r[x1+1][ind] - r1[ind]) ** 2 * sincos2[ind] * corr2[ind] 
        tmatrix = tmatrix + ssp.csr_matrix((temp, (x1[ind], col)), shape = tmatrixshape)
        tmatrix = tmatrix + ssp.csr_matrix(((1-temp), (x4[ind], col)), shape = tmatrixshape)
        
        #x1
        #3 x1=x2<x3=x4    #4 x1=x2<x3<x4
        ind = np.logical_or(opind[2], opind[3]).ravel()
        col = col0[ind]
        temp = (r2[ind] - r1[ind]) ** 2 * sincos2[ind] * corr2[ind] + (r[x2+1][ind] - r2[ind]) * sincos[ind] * corr[ind]
        tmatrix = tmatrix + ssp.csr_matrix((temp, (x1[ind], col)), shape = tmatrixshape)
        #6 x1<x2=x3<x4    #7 x1<x2<x3=x4    #8 x1<x2<x3<x4
        ind = np.logical_or(np.logical_or(opind[5], opind[6]), opind[7]).ravel()
        col = col0[ind]
        temp = (r[x1+1][ind] - r1[ind]) ** 2 * sincos2[ind] * corr2[ind] 
        tmatrix = tmatrix + ssp.csr_matrix((temp, (x1[ind], col)), shape = tmatrixshape)
        
        #x2
        #6 x1<x2=x3<x4
        ind = opind[5].ravel()
        col = col0[ind]
        temp = ((r2[ind] - r1[ind]) ** 2 - (r[x2][ind] - r1[ind]) ** 2 + (r4[ind] - r3[ind]) ** 2 - (r4[ind] - r[x3+1][ind]) ** 2) \
               * sincos2[ind] * corr2[ind] + (r3[ind] - r2[ind]) * sincos[ind] * corr[ind]
        tmatrix = tmatrix + ssp.csr_matrix((temp, (x2[ind], col)), shape = tmatrixshape)
        #7 x1<x2<x3=x4   #8 x1<x2<x3<x4
        ind = np.logical_or(opind[6], opind[7]).ravel()
        col = col0[ind]
        temp = ((r2[ind] - r1[ind]) ** 2 - (r[x2][ind] - r1[ind]) ** 2) * sincos2[ind] * corr2[ind] +\
               (r[x2+1][ind] - r2[ind]) * sincos[ind] * corr[ind]
        tmatrix = tmatrix + ssp.csr_matrix((temp, (x2[ind], col)), shape = tmatrixshape)
        
        #x3
        #4 x1=x2<x3<x4    #8 x1<x2<x3<x4
        ind = np.logical_or(opind[3], opind[7]).ravel()
        col = col0[ind]
        temp = (r3[ind] - r[x3][ind]) * sincos[ind] * corr[ind] + \
               ((r4[ind] - r3[ind]) ** 2 - (r4[ind] - r[x3+1][ind]) ** 2) * sincos2[ind] * corr2[ind] 
        tmatrix = tmatrix + ssp.csr_matrix((temp, (x3[ind], col)), shape = tmatrixshape)
        
        #x4
        #3 x1=x2<x3=x4    #7 x1<x2<x3=x4
        ind = np.logical_or(opind[2], opind[6]).ravel()
        col = col0[ind]
        temp = (r4[ind] - r3[ind]) ** 2 * sincos2[ind] * corr2[ind] + (r3[ind] - r[x3][ind]) * sincos[ind] * corr[ind]
        tmatrix = tmatrix + ssp.csr_matrix((temp, (x4[ind], col)), shape = tmatrixshape)
        #4 x1=x2<x3<x4    #6 x1<x2=x3<x4    #8 x1<x2<x3<x4
        ind = np.logical_or(np.logical_or(opind[3],opind[5]), opind[7]).ravel()
        col = col0[ind]
        temp = (r4[ind] - r[x4][ind]) ** 2 * sincos2[ind] * corr2[ind]
        tmatrix = tmatrix + ssp.csr_matrix((temp, (x4[ind], col)), shape = tmatrixshape)
        
        #x1+1~x2
        #6 x1<x2=x3<x4    #7 x1<x2<x3=x4    #8 x1<x2<x3<x4
        ind = np.logical_or(np.logical_or(opind[5], opind[6]), opind[7]).ravel()
        col = col0[ind]
        for i in range(maxi):
            tempx = x1+1+i
            ind1 = tempx < x2
            ind1 = np.logical_and(ind1, ind)
            col = col0[ind1]
            temp = ((r[tempx+1][ind1] - r1[ind1]) ** 2 - (r[tempx][ind1] - r1[ind1]) ** 2) * sincos2[ind1] * corr2[ind1]
            tmatrix = tmatrix + ssp.csr_matrix((temp, (tempx[ind1], col)), shape = tmatrixshape)
        
        #x2+1~x3
        #3 x1=x2<x3=x4    #4 x1=x2<x3<x4    #7 x1<x2<x3=x4    #8 x1<x2<x3<x4
        ind = np.logical_or(np.logical_or(np.logical_or(opind[2], opind[3]), opind[6]), opind[7]).ravel()
        col = col0[ind]
        for i in range(maxi):
            tempx = x2+1+i
            ind1 = tempx < x3
            ind1 = np.logical_and(ind1, ind)
            col = col0[ind1]
            temp = (r[tempx+1][ind1] - r[tempx][ind1]) * sincos[ind1] * corr[ind1]
            tmatrix = tmatrix + ssp.csr_matrix((temp, (tempx[ind1], col)), shape = tmatrixshape)
        
        #x3+1~x4
        #6 x1<x2=x3<x4    #8 x1<x2<x3<x4
        ind = np.logical_or(opind[5], opind[7]).ravel()
        col = col0[ind]
        for i in range(maxi):
            tempx = x3+1+i
            ind1 = tempx < x4
            ind1 = np.logical_and(ind1, ind)
            col = col0[ind1]
            temp = ((r4[ind1] - r[tempx][ind1]) ** 2 - (r4[ind1] - r[tempx+1][ind1]) ** 2) * sincos2[ind1] * corr2[ind1]
            tmatrix = tmatrix + ssp.csr_matrix((temp, (tempx[ind1], col)), shape = tmatrixshape)
        return tmatrix
        
    def _transEllipseParm(self, p):
        """calcuate the standard form of ellipse from given parameter
        input:
        p:[xc,yc,x0,y0,z0,tth]
        ((xc-x0)**2+(yc-y0)**2+z0**2)((x-x0)**2+(y-y0)**2+z0**2)*Cos(tth)**2 =
        ((xc-x0)(x-x0)+(yc-y0)(y-y0)+z0**2)**2
        return:
        normalized a,b,c,d,e,f:
        a*x**2 + b*y**2 + c*x*y + d*x + e*y + f = 0
        """
        xc, yc, x0, y0, z0, tth=p
        cc = (np.cos(tth)) ** 2
        cc1 = ((x0 - xc) ** 2 + (y0 - yc) ** 2 + z0 ** 2)
        cc2 = (x0 * (x0 - xc) + y0 * (y0 - yc) + z0 ** 2)
        a = -(x0 - xc) ** 2 + cc * cc1
        b = -(y0 - yc) ** 2 + cc * cc1
        c = - 2 * (yc - y0) * (xc - x0)
        d = 2 * (x0 - xc) * cc2 - 2 * x0 * cc * cc1
        e = 2 * (y0 - yc) * cc2 - 2 * y0 * cc * cc1
        f = - cc2 ** 2 + cc * cc1 * (x0 ** 2 + y0 ** 2 + z0 ** 2)
        re = [a, b, c, d, e, f]
        return re
    
