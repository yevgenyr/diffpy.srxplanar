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
import scipy.io
import os
from diffpy.srxplanar.srxplanarconfig import _configPropertyR

class SaveResults(object):
    '''module to save results into files
    '''
    integrationspace = _configPropertyR('integrationspace')
    savedirectory = _configPropertyR('savedirectory')
    gsasoutput = _configPropertyR('gsasoutput')
    filenameplus = _configPropertyR('filenameplus')
    
    def __init__(self, p):
        self.config = p
        self.prepareCalculation()
        return
    
    def prepareCalculation(self):
        if not os.path.exists(self.savedirectory):
                os.makedirs(self.savedirectory)
        return

    def getFilePathWithoutExt(self, filename):
        '''get the normalized full path of filename to be saved
        param filename: string, could be full path/ file name only/ with/without ext, only the base part of filename is used.
        '''
        filebase = os.path.splitext(os.path.split(filename)[1])[0]
        if self.filenameplus!='' and self.filenameplus!=None:
            filenamep = '_'.join([filebase, self.filenameplus, self.integrationspace])
        else:
            filenamep = '_'.join([filebase, self.integrationspace])
        filepathwithoutext = self.savedirectory + '/' + filenamep
        return filepathwithoutext
        
    def save(self, rv):
        '''save diffraction intensity
        '''
        self.saveChi(rv['chi'], rv['filename'])
        if self.gsasoutput:
            if self.gsasoutput in set(['std', 'esd', 'fxye']):
                self.saveGSAS(rv['chi'], rv['filename'])
        return
    
    def saveChi(self, xrd, filename):
        '''save xrd diffraction pattern in .chi
        '''
        filepath = self.getFilePathWithoutExt(filename) + '.chi'
        f = open(filepath, 'wb')
        f.write(self.config.getHeader(mode='short'))
        np.savetxt(f, xrd.transpose(), fmt='%g')
        f.close()
        return
    
    def saveGSAS(self, xrd, filename):
        '''save xrd diffraction pattern in .gsas
        '''
        filepath = self.getFilePathWithoutExt(filename) + '.gsas'
        f = open(filepath,'wb')
        f.write(self.config.getHeader(mode='short'))
        if xrd.shape[0]==3:
            s = writeGSASStr(os.path.splitext(path)[0], self.gsasoutput, xrd[0], xrd[1], xrd[2])
        elif xrd.shape[0]==2:
            s = writeGSASStr(os.path.splitext(path)[0], self.gsasoutput, xrd[0], xrd[1])
        f.write(s)
        f.close()
        return
  
def writeGSASStr(name, mode, tth, iobs, esd=None):
    """Return string of integrated intensities in GSAS format.
    mode:    string, gsas file type, could be 'std', 'esd', 'fxye' (gsas format)
    tth:     ndarray, two theta angle
    iobs:    ndarray, Xrd intensity
    esd:     ndarray, optional error value of intensity
    
    return:  string, a string to be saved to file
    """
    maxintensity = 999999
    logscale = numpy.floor(numpy.log10(maxintensity / numpy.max(iobs)))
    logscale = min(logscale, 0)
    scale = 10 ** int(logscale)
    lines = []
    ltitle = 'Angular Profile'
    ltitle += ': %s' % name
    ltitle += ' scale=%g' % scale
    if len(ltitle) > 80:    ltitle = ltitle[:80]
    lines.append("%-80s" % ltitle)
    ibank = 1
    nchan = len(iobs)
    # two-theta0 and dtwo-theta in centidegrees
    tth0_cdg = tth[0] * 100
    dtth_cdg = (tth[-1] - tth[0]) / (len(tth) - 1) * 100
    if esd==None: mode='std'
    if mode=='std':
        nrec = int(numpy.ceil(nchan/10.0))
        lbank = "BANK %5i %8i %8i CONST %9.5f %9.5f %9.5f %9.5f STD" % \
                (ibank, nchan, nrec, tth0_cdg, dtth_cdg, 0, 0)
        lines.append("%-80s" % lbank)
        lrecs = [ "%2i%6.0f" % (1, ii * scale) for ii in iobs ]
        for i in range(0, len(lrecs), 10):
            lines.append("".join(lrecs[i:i+10]))
    if mode=='esd':
        nrec = int(numpy.ceil(nchan/5.0))
        lbank = "BANK %5i %8i %8i CONST %9.5f %9.5f %9.5f %9.5f ESD" % \
                (ibank, nchan, nrec, tth0_cdg, dtth_cdg, 0, 0)
        lines.append("%-80s" % lbank)
        lrecs = [ "%8.0f%8.0f" % (ii, ee * scale) for ii, ee in zip(iobs, esd) ]
        for i in range(0, len(lrecs), 5):
            lines.append("".join(lrecs[i:i+5]))
    if mode=='fxye':
        nrec = nchan
        lbank = "BANK %5i %8i %8i CONST %9.5f %9.5f %9.5f %9.5f FXYE" % \
                (ibank, nchan, nrec, tth0_cdg, dtth_cdg, 0, 0)
        lines.append("%-80s" % lbank)
        lrecs = [ "%22.10f%22.10f%24.10f" % (xx * scale, yy * scale, ee * scale) for xx, yy, ee in zip(tth, iobs, esd) ]
        for i in range(len(lrecs)):
            lines.append("%-80s" % lrecs[i])
    lines[-1] = "%-80s" % lines[-1]
    rv = "\r\n".join(lines) + "\r\n"
    return rv