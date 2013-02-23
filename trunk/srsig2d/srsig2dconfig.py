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
import ConfigParser
import re
import os

class SrSig2DConfig(ConfigParser.ConfigParser):
    '''Class used for storing the configuration value. It bases on the configparser class provided by python core'''
    def __init__(self, filename=None):
        ConfigParser.ConfigParser.__init__(self)
        rootpath = os.getcwd()
        #files
        self.add_section('BaseSection')
        self.addBase('fit2dconfig', '')
        self.addBase('tifdirectory', rootpath)
        self.addBase('includepattern', '*.tif')
        self.addBase('excludepattern', ['*.dark.tif'])
        #image process
        self.addBase('fliphorizontal', False)
        self.addBase('flipvertical', True)
        #experiment conditions
        self.addBase('backgroundfile', '')
        #self.addBase('backgroundenable', False)
        #experiment condition
        self.addBase('wavelength', 0.1000)
           
        #geometry parameters        
        self.addBase('xdimension', 2048)
        self.addBase('ydimension', 2048)
        self.addBase('xpixelsize', 200e-3)
        self.addBase('ypixelsize', 200e-3)
        self.addBase('xbeamcenter', 1024.0)
        self.addBase('ybeamcenter', 1024.0)
        self.addBase('distance', 200.0)
        self.addBase('rotationd', 0.0)
        self.addBase('tiltd', 0.0)
        self.addBase('xrdtthstepd', 0.04)
        self.addBase('xrdtthmaxd', 40.0)
        self.addBase('xrdqmax', 40.0)
        self.addBase('xrdqstep', 0.04)
        self.addBase('xrdtthvarstepd', 0.01)
        self.addBase('xrdqvarstep', 0.01)
        self.addBase('xrdazimuthstepd', 0.5)
        self.addBase('xrduncertaintyenable', False)
        self.addBase('xrdspotty', False)
        
        for name in ['rotation', 'tilt', 'xrdtthstep', 'xrdtthmax', 'xrdtthvarstep', 'xrdazimuthstep']:
            setattr(self.__class__, name, self.configPropertyR(name+'d'))

        #correction
        self.addBase('sacorrectionenable', False)
        self.addBase('polcorrectionenable', False)
        self.addBase('polcorrectf', 0.95)
        self.addBase('maskselfcorrenable', False)
        self.addBase('maskselfcorrmethod', "Filter by raw counts") #"Filter by raw counts", "Filter by percentage", "Both"
        self.addBase('maskselfcorrcounthi', 0.5)
        self.addBase('maskselfcorrcountlow', 0.5)
        self.addBase('maskselfcorrpercentagehi', 0.95)
        self.addBase('maskselfcorrpercentagelow', 0.05)
        self.addBase('maskselfcorraccuenable', False)
        
        #integration method)
        self.addBase('integrationspace', 'twotheta') #'twotheta','qspace'
        self.addBase('method', 'normal') #'normal','split'
        self.addBase('regulartmatrixenable', False)
        
        #mask
        self.addBase('maskenable', True)
        self.addBase('maskfit2denable', False)
        self.addBase('maskfit2d', '')
        self.addBase('masktiffenable', False)
        self.addBase('masktiff', '')
        self.addBase('maskcakeenable', False)
        self.addBase('maskcake', [])
        self.addBase('maskboxenable', False)
        self.addBase('maskbox', [])
        
        #save
        self.addBase('savedirectory', rootpath+'/re')
        self.addBase('gsasoutput', 'None') #'None', 'std', 'esd', 'fxye'
        self.addBase('filenameplus', '')
        
        self.configlist = self.options('BaseSection')
        self.configlist.extend(['rotation', 'tilt', 'xrdtthstep', 'xrdtthmax', 'xrdtthvarstep', 'xrdazimuthstep'])
        self.strlistoption = set(['excludepattern'])
        
        if filename!=None:
            self.loadFromFile(filename)
        if (self.fit2dconfig!='')and os.path.exists(self.fit2dconfig):
            self.loadFromFit2D(self.fit2dconfig)
        return
    
    def configPropertyR(self, nm):
        '''helper function of options delegation, rad 2 degree'''
        rv = property(fget = lambda self: np.radians(getattr(self, nm)), 
                      fset = lambda self, val: setattr(self, nm, np.radians(val)), 
                      fdel = lambda self: delattr(self, nm))
        return rv    
    
    def addBase(self, optionsname, optionsvalue):
        '''add options to base section with automatically determined type
        '''
        if type(optionsvalue) in set([int, float, str, bool]):
            self.set('BaseSection', optionsname, str(optionsvalue))
            setattr(self, optionsname, optionsvalue)
        elif type(optionsvalue) == list:
            self.set('BaseSection', optionsname, str(optionsvalue)[1:-1])
            setattr(self, optionsname, optionsvalue)
        return
    
    def loadFromFile(self, filename):
        '''load config from file
        '''
        self.read(filename)
        self._updateConfigtoSelf()
        return
    
    def _updateConfigtoSelf(self):
        '''copy the values in configParser container to self.'options'
        '''
        for sectionname in self.sections():
            for optionname in self.options(sectionname):
                if type(getattr(self, optionname)) == str:
                    setattr(self, optionname, self.get(sectionname, optionname))
                elif type(getattr(self, optionname)) == int:
                    setattr(self, optionname, self.getint(sectionname, optionname))
                elif type(getattr(self, optionname)) == float:
                    setattr(self, optionname, self.getfloat(sectionname, optionname))
                elif type(getattr(self, optionname)) == bool:
                    setattr(self, optionname, self.getboolean(sectionname, optionname))
                elif type(getattr(self, optionname)) == list:
                    if optionname=='maskcake':      
                        setattr(self, optionname, self._getFloatList(sectionname, optionname))
                    if optionname=='maskbox':
                        setattr(self, optionname, self._getIntList(sectionname, optionname))
                    if optionname in self.strlistoption:      
                        setattr(self, optionname, self._getStrList(sectionname, optionname))
        return
    
    def updateConfig(self):
        '''update the options value, 
        then write the values in the self.'options' to configParser container
        '''
        self.checkMax()
        self.checkStep()
        if self.integrationspace == 'twotheta':
            self.xrdtthorqmax = self.xrdtthmax
            self.xrdtthorqstep = self.xrdtthstep
            self.qmaxfromxrd = 4 * np.pi * np.sin(self.xrdtthmax/2) / self.wavelength
        elif self.integrationspace == 'qspace':
            self.xrdtthorqmax = self.xrdqmax
            self.xrdtthorqstep = self.xrdqstep
            self.qmaxfromxrd = self.xrdqmax
        
        for sectionname in self.sections():
            for optionname in self.options(sectionname):
                if type(getattr(self, optionname)) in set([str, int, bool, float]):
                    self.set(sectionname, optionname, str(getattr(self, optionname)))
                elif type(getattr(self, optionname)) == list:
                    liststring = ''
                    for listitem in getattr(self, optionname):
                        liststring = liststring + str(listitem) + ', '
                    self.set(sectionname, optionname, liststring[:-2])
        return
    
    def checkMax(self):
        '''check tthmax and qmax, change them if they are smaller than actual tthmax or qmax 
        '''
        xr = (np.array([0, self.xdimension]) - self.xbeamcenter) * self.xpixelsize
        yr = (np.array([0, self.ydimension]) - self.ybeamcenter) * self.ypixelsize
        sinr = np.sin(self.rotation)
        cosr = np.cos(self.rotation)
        sint = np.sin(self.tilt)
        cost = np.cos(self.tilt)
        sourcexr = self.distance * sint * cosr
        sourceyr = self.distance * sint * sinr
        sourcezr = self.distance * cost
        
        dmatrix = ((xr - sourcexr) ** 2).reshape(1, 2) + \
                  ((yr - sourceyr) ** 2).reshape(2, 1) + sourcezr ** 2
        dmatrix = np.sqrt(dmatrix)
        tthmatrix1 = (-(xr - sourcexr) * sint * cosr).reshape(1, 2) + \
                     (-(yr - sourceyr) * sint * sinr).reshape(2, 1) + sourcezr * cost
        tthmatrix = np.arccos(tthmatrix1 / dmatrix)
        qmatrix = 4 * np.pi * np.sin(tthmatrix / 2.0) / self.wavelength
        
        if np.max(tthmatrix)>self.xrdtthmax:
            self.xrdtthmax = np.max(tthmatrix) + 0.1
        if np.max(qmatrix)>self.xrdqmax:
            self.xrdqmax = np.max(qmatrix) + 1.0
        return
    
    def checkStep(self):
        tthstep = self.xpixelsize / self.distance 
        qstep = 4 * np.pi * np.sin(tthstep / 2.0) / self.wavelength
        if np.abs(tthstep - self.xrdtthstep)/tthstep > 0.05:
            self.xrdtthstep = tthstep
        if np.abs(qstep - self.xrdqstep)/qstep > 0.05:
            self.xrdqstep = qstep
        return
        
    def writeConfig(self, filename):
        '''write config to file
        '''
        self.updateConfig()
        conffile = open(filename, 'w')
        self.write(conffile)
        conffile.close()
        return
    
    def loadFromFit2D(self, filename):
        '''load parameters from fit2d calibration information. copy/paste the fit2d calibration 
        results to a txt file. this function will load xbeamcenter, ybeamceter... from the file
        '''
        def findFloat(line):
            temp = re.findall('[-+]?\d*\.\d+|[-+]?\d+', line)
            return map(float, temp)
            
        f = open(filename, 'r')
        lines = f.readlines()
        for line in lines:
            if re.search('Refined Beam centre.*pixels', line):
                self.xbeamcenter, self.ybeamcenter = findFloat(line)
            elif re.search('Refined sample to detector distance', line):
                self.distance = findFloat(line)[0]
            elif re.search('Refined wavelength', line):
                self.wavelength = findFloat(line)[0]
            elif re.search('Refined tilt plane rotation angle', line):
                self.rotationd = findFloat(line)[0]
            elif re.search('Refined tilt angle', line):
                self.tiltd = findFloat(line)[0]
            elif re.search('Refined wavelength', line):
                self.wavelength = findFloat(line)[0]
        f.close()
        self.updateConfig()
        return
            
    
    def _getIntList(self, sectionname, optionname):
        '''helpler function that get a list of int from one string
        '''
        temp = re.split('\s*,\s*', self.get(sectionname, optionname))
        if len(temp)>1:
            intlist = map(int, temp)
        else:
            intlist = []
        return intlist
    
    def _getFloatList(self, sectionname, optionname):
        '''helpler function that get a list of float from one string
        '''
        temp = re.split('\s*,\s*', self.get(sectionname, optionname))
        if len(temp)>1:
            floatlist = map(float, temp)
        else:
            floatlist = []
        return floatlist
    
    def _getStrList(self, sectionname, optionname):
        '''helpler function that get a list of string from one string(separated by ',')
        '''
        temp = re.split('\s*,\s*', self.get(sectionname, optionname))
        if len(temp)>1:
            strlist = temp
        else:
            strlist = []
        return strlist
    
if __name__=='__main__':
    a = SrSig2DConfig()
    a.loadFromFile('test.cfg')
    a.updateConfig()
    a.writeConfig('test.cfg')
