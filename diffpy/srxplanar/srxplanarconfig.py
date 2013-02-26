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

class SrXPlanarConfig(ConfigParser.ConfigParser):
    '''Class used for storing the configuration value. It bases on the configparser class provided by python core'''
    def __init__(self, filenames=None):
        ConfigParser.ConfigParser.__init__(self)
        rootpath = os.getcwd()
        
        self.add_section('Experiment')
        self._addExp('beamline', 'X17A')
        self._addExp('fit2dconfig', '')
        self._addExp('tifdirectory', rootpath)
        self._addExp('savedirectory', rootpath+'/re')
        self._addExp('backgroundfile', '')
        self._addExp('maskfit2d', '')
        self._addExp('integrationspace', 'twotheta') #'twotheta','qspace'
        self._addExp('wavelength', 0.1000)   
        self._addExp('xbeamcenter', 1024.0)
        self._addExp('ybeamcenter', 1024.0)
        self._addExp('distance', 200.0)
        self._addExp('rotationd', 0.0)
        self._addExp('tiltd', 0.0)
        self._addExp('tthstepd', 0.04)
        self._addExp('tthmaxd', 40.0)
        self._addExp('qmax', 40.0)
        self._addExp('qstep', 0.04)
        for name in ['rotation', 'tilt', 'tthstep', 'tthmax']:
            setattr(self.__class__, name, self.configPropertyR(name+'d'))
        
        self.add_section('Beamline')
        self._addBeamline('includepattern', '*.tif')
        self._addBeamline('excludepattern', ['*.dark.tif'])
        self._addBeamline('fliphorizontal', False)
        self._addBeamline('flipvertical', True)
        self._addBeamline('xdimension', 2048)
        self._addBeamline('ydimension', 2048)
        self._addBeamline('xpixelsize', 200e-3)
        self._addBeamline('ypixelsize', 200e-3)
        
        self.add_section('Others')
        self._addOthers('uncertaintyenable', False)
        self._addOthers('sacorrectionenable', True)
        self._addOthers('polcorrectionenable', False)
        self._addOthers('polcorrectf', 0.95)
        self._addOthers('selfcorrenable', True)
        self._addOthers('gsasoutput', 'None') #'None', 'std', 'esd', 'fxye'
        self._addOthers('filenameplus', '')        
        
        self.configlistexperiment = self.options('Experiment')
        self.configlistbeamline = self.options('Beamline')
        self.configlistothers = self.options('Others')
        self.configlistexperiment.extend(['rotation', 'tilt', 'tthstep', 'tthmax'])
        
        self.strlistoption = set(['excludepattern'])
        self.intlistoption = set([])
        self.floatlistoption = set([])
        
        #init from config file(s)
        self.loadConfig(filenames)
        return
    
    def configPropertyR(self, nm):
        '''helper function of options delegation, rad 2 degree'''
        rv = property(fget = lambda self: np.radians(getattr(self, nm)), 
                      fset = lambda self, val: setattr(self, nm, np.degrees(val)), 
                      fdel = lambda self: delattr(self, nm))
        return rv    
    
    def _addOpt(self, sectionname, optionsname, optionsvalue):
        '''add options to section with automatically determined type
        '''
        if type(optionsvalue) in set([int, float, str, bool]):
            self.set(sectionname, optionsname, str(optionsvalue))
            setattr(self, optionsname, optionsvalue)
        elif type(optionsvalue) == list:
            self.set(sectionname, optionsname, str(optionsvalue)[1:-1])
            setattr(self, optionsname, optionsvalue)
        return
    
    def _addExp(self, optionsname, optionsvalue):
        '''add options to Experiment section with automatically determined type
        '''
        self._addOpt('Experiment', optionsname, optionsvalue)
        return
    
    def _addBeamline(self, optionsname, optionsvalue):
        '''add options to Beamline section with automatically determined type
        '''
        self._addOpt('Beamline', optionsname, optionsvalue)
        return
    
    def _addOthers(self, optionsname, optionsvalue):
        '''add options to Others section with automatically determined type
        '''
        self._addOpt('Others', optionsname, optionsvalue)
        return
    
    def loadConfig(self, filenames):
        '''load config from config file or list of config file
        :param filenames: str or list of str, config file name
        '''
        configfilelist = filenames if type(filenames)==list else [filenames]
        for filename in configfilelist:
            self.loadFromFile(filename)
        #init from fit2d config if specified
        self.loadFromFit2D(self.fit2dconfig)
        #load default beamline config if specifed (only X17A is support by now)
        if self.beamline.lower in set(['x17a']):
            self.loadFromFile('x17a.cfg')
        else:
            self.loadFromFile(self.beamline)
        return
    
    def loadFromFile(self, filename):
        '''load config from file
        '''
        if filename!=None:
            if os.path.exists(filename):
                self.read(filename)
                self._updateConfigtoSelf()
                self.updateConfig()
        return
    
    def loadFromFit2D(self, filename):
        '''load parameters from fit2d calibration information. copy/paste the fit2d calibration 
        results to a txt file. this function will load xbeamcenter, ybeamceter... from the file
        '''
        def findFloat(line):
            temp = re.findall('[-+]?\d*\.\d+|[-+]?\d+', line)
            return map(float, temp)
        if filename != None:
            if os.path.exists(filename):
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
                    if optionname in self.strlistoption:      
                        setattr(self, optionname, self._getStrList(sectionname, optionname))
                    elif optionname in self.intlistoption:      
                        setattr(self, optionname, self._getIntList(sectionname, optionname))
                    elif optionname in self.floatlistoption:      
                        setattr(self, optionname, self._getFloatList(sectionname, optionname))
        return
    
    def updateConfig(self):
        '''update the options value, 
        then write the values in the self.'options' to configParser container
        '''
        self._checkMax()
        self._checkStep()
        if self.integrationspace == 'twotheta':
            self.tthorqmax = self.tthmax
            self.tthorqstep = self.tthstep
        elif self.integrationspace == 'qspace':
            self.tthorqmax = self.qmax
            self.tthorqstep = self.qstep
        
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
            
    def writeConfig(self, filename):
        '''write config to file
        '''
        self.updateConfig()
        conffile = open(filename, 'w')
        self.write(conffile)
        conffile.close()
        return

    def _checkMax(self):
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
        
        if np.max(tthmatrix)>self.tthmax:
            self.tthmax = np.max(tthmatrix) + 0.05
        if np.max(qmatrix)>self.qmax:
            self.qmax = np.max(qmatrix) + 0.1
        return
    
    def _checkStep(self):
        tthstep = self.xpixelsize / self.distance 
        qstep = 4 * np.pi * np.sin(tthstep / 2.0) / self.wavelength
        if np.abs(tthstep - self.tthstep)/tthstep > 0.05:
            self.xrdtthstep = tthstep
        if np.abs(qstep - self.qstep)/qstep > 0.05:
            self.xrdqstep = qstep
        return
    
    def _getIntList(self, sectionname, optionname):
        '''helpler function that get a list of int from one string
        '''
        temp = re.split('\s*,\s*', self.get(sectionname, optionname))
        if len(temp)>0:
            intlist = map(int, temp)
        else:
            intlist = []
        return intlist
    
    def _getFloatList(self, sectionname, optionname):
        '''helpler function that get a list of float from one string
        '''
        temp = re.split('\s*,\s*', self.get(sectionname, optionname))
        if len(temp)>0:
            floatlist = map(float, temp)
        else:
            floatlist = []
        return floatlist
    
    def _getStrList(self, sectionname, optionname):
        '''helpler function that get a list of string from one string(separated by ',')
        '''
        temp = re.split('\s*,\s*', self.get(sectionname, optionname))
        if len(temp)>0:
            strlist = temp
        else:
            strlist = []
        return strlist
    
if __name__=='__main__':
    a = SrXPlanarConfig()
    a.loadFromFile('test.cfg')
    a.updateConfig()
    a.writeConfig('test.cfg')
