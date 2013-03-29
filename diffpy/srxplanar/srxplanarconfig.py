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
import re, os, sys
from functools import partial
import argparse
from diffpy.confutils.config import ConfigBase, _configPropertyRad, _configPropertyR, _configPropertyRW

class SrXplanarConfig(ConfigBase):
    '''Class used for storing the configuration value. It bases on the configparser class provided by python core'''
    '''Class used for storing the configuration value and process cmd command
    config: ConfigParser.ConfigParser
    args: argparse
    '''
    
    # optdata contains these keys:
    # full(f), short(s), help(h), type(t), action(a), nargs(n), default(d), choices(c), required(r), dest, const
    
    '''optional args: 
        'config', set False to hide this option in self.config
        'args', set False to hide this option in self.args
        'header', set False to hide this option in header output
        
        add options will process in the following order:
        _optdatalistpre, _optdatalist(defined in ConfigBase), _optdatalistext
        this is also the order that options present in self.args and self.config
    '''
    _optdatalistpre = [
        ['filenames',{'sec':'Control', 'config':'n', 'header':'n',
            'f':'filename',
            'h':'filename or list of filenames or filename pattern or list of filename pattern',
            'n':'*',
            'd':[],}],
        ['output',{'sec':'Experiment', 'config':'n', 'header':'n',
            's':'o',
            'h':'base filename of output file',
            'd':'',}],
        ]
    
    _optdatalist = [
        ['configfile',{'sec':'Control', 'config':'n', 'header':'n',
            's':'c',
            'h':'name of input config file',
            'd':'',}],
        ['createconfig',{'sec':'Control', 'config':'n', 'header':'n',
            'h':'create a config file according to default or current values',
            'd':'',}],
        ['createconfigfull',{'sec':'Control', 'config':'n', 'header':'n',
            'h':'create a full configurable config file',
            'd':'',}],
        ]
    
    _optdatalistext = [
        # control group
        ['summation',{'sec':'Control', 'header':'n',
            's':'s',
            'h':'sum all the image and then integrate',
            'n':'?',
            'co':True,
            'd':False,}],
        #Expeiment gropu
        ['fit2dconfig',{'sec':'Experiment',
            's':'fit2d',
            'h':'fit2d calibration file name. It contains the calibration results copy from fit2d cmd windows',
            'd':'',}],
        ['tifdirectory',{'sec':'Experiment', 'header':'n',
            's':'tifdir',
            'h':'directory of raw tif files',
            'd':'currentdir',}],
        ['savedirectory',{'sec':'Experiment', 'header':'n',
            's':'savedir',
            'h':'directory of save files',
            'd':'currentdir',}],
        ['backgroundfile',{'sec':'Experiment',
            'h':'background file name, should be a image file',
            'd':'',}],
        ['addmask',{'sec':'Experiment',
            'h':'options to control masks, add "deadpixel" to mask the deadpixel, add "spot" \
                to mask the spot, or names of mask file, support .npy and .tif and .msk(fit2dmsk) \
                add "selfcorr" to mask the pixel whose intensity is too high/low compare to pixels in the same bin',
            'n':'*',
            'd':[],}],
        ['createmask',{'sec':'Control', 'config':'n','header':'n',
            'h':'create a mask file according to current image file and value of addmask, maskedges',
            'd':'',}],
        ['integrationspace',{'sec':'Experiment',
            'h':'integration space, could be twotheta or qspace',
            'd':'twotheta',
            'c':['twotheta','qspace'],}],
        ['wavelength',{'sec':'Experiment',
            'h':'wavelength of x-ray, in A',
            'd':0.1000,}],
        ['xbeamcenter',{'sec':'Experiment',
            's':'xc',
            'h':'beamcenter in x axis, in pixel',
            'd':1024.0,}],
        ['ybeamcenter',{'sec':'Experiment',
            's':'yc',
            'h':'beamcenter in y axis, in pixel',
            'd':1024.0,}],
        ['distance',{'sec':'Experiment',
            's':'dis',
            'h':'distance between detector and sample, in mm',
            'd':200.0,}],
        ['rotationd',{'sec':'Experiment',
            's':'rot',
            'h':'rotation angle of tilt plane, in degree',
            'd':0.0,}],
        ['tiltd',{'sec':'Experiment',
            's':'tilt',
            'h':'tilt angle of tilt plane, in degree',
            'd':0.0,}],
        ['tthstepd',{'sec':'Experiment',
            's':'ts',
            'h':'integration step in twotheta space, in degree',
            'd':0.02,}],
        ['qstep',{'sec':'Experiment',
            's':'qs',
            'h':'integration step in q space, in A^-1',
            'd':0.02,}],
        #Beamline group
        ['includepattern',{'sec':'Beamline','header':'n','config':'f',
            's':'ipattern',
            'h':'file name pattern for included files',
            'n':'*',
            'd':['*.tif'],}],
        ['excludepattern',{'sec':'Beamline','header':'n','config':'f',
            's':'epattern',
            'h':'file name pattern for excluded files',
            'n':'*',
            'd':['*.dark.tif', '*.raw.tif'],}],
        ['fliphorizontal',{'sec':'Beamline','header':'n','config':'f',
            'h':'filp the image horizontally',
            'n':'?',
            'co':True,
            'd':False,}],
        ['flipvertical',{'sec':'Beamline','header':'n','config':'f',
            'h':'filp the image vertically',
            'n':'?',
            'co':True,
            'd':True,}],
        ['xdimension',{'sec':'Beamline',
            's':'xd',
            'h':'detector dimension in x axis, in pixel',
            'd':2048,}],
        ['ydimension',{'sec':'Beamline',
            's':'yd',
            'h':'detector dimension in y axis, in pixel',
            'd':2048,}],
        ['xpixelsize',{'sec':'Beamline',
            's':'xp',
            'h':'detector pixel size in x axis, in mm',
            'd':0.2,}],
        ['ypixelsize',{'sec':'Beamline',
            's':'yp',
            'h':'detector pixel size in y axis, in mm',
            'd':0.2,}],
        #Others Group
        ['uncertaintyenable',{'sec':'Others',
            's':'error',
            'h':'propagate uncertainty',
            'n':'?',
            'co':True,
            'd':False,}],
        ['sacorrectionenable',{'sec':'Others','config':'f',
            's':'sacorr',
            'h':'enable solid angle correction',
            'n':'?',
            'co':True,
            'd':True,}],
        ['polcorrectionenable',{'sec':'Others','config':'f',
            's':'polarcorr',
            'h':'enable polarization correction',
            'n':'?',
            'co':True,
            'd':True,}],        
        ['polcorrectf',{'sec':'Others','config':'f',
            's':'polarf',
            'h':'polarization correction factor',
            'd':0.99,}],
        ['selfcorrenable',{'sec':'Others','config':'f',
            's':'selfcorr',
            'h':'mask the pixels whose intensity is too high/low compare to other pixels in the same bin',
            'n':'?',
            'co':True,
            'd':True,}],
        ['gsasoutput',{'sec':'Others','header':'n',
            'h':'select if want to output gsas format file',
            'c':['None', 'std', 'esd', 'fxye'],
            'd':'None',}],
        ['filenameplus',{'sec':'Others','header':'n',
            'h':'str that added behind normal file name',
            'd':'',}],
        ['maskedges',{'sec':'Others','config':'f',
            'h':'mask the edge pixels, first four means the number of pixels masked in each edge \
                (left, right, top, bottom), the last one is the radius of a region masked around the corner',
            'n':5,
            'd':[10,10,10,10,100],}],
        ['nocalculation',{'sec':'Others','config':'f','header':'n',
            'h':'set True to disable all calculation, will automaticly set True if createconfig or createmask',
            'n':'?',
            'co':True,
            'd':False,}],
        ]
    
    #options that consist of list of value
    _listoptdata = {
        'strlist':['excludepattern'],
        'floatlist':[],
        'intlist':['maskedges'],
        'strlist':[],
        }
    #configlist, store the options name for each sections
    #_configlist = {} 
    
    #default config file path and name, overload it for your config class
    _defaultconfigpath = ['srxplanar.cfg', 'SrXplanar.cfg']
    
    def _additionalInit(self):
        '''method called in init process, overload it
        this method will be called after all options in self._optdata are processced and before reading config from file/args/kwargs
        '''
        for name in ['rotation', 'tilt', 'tthstep', 'tthmax']:
            setattr(self.__class__, name, _configPropertyRad(name+'d'))
        self._configlist['Experiment'].extend(['rotation', 'tilt', 'tthstep', 'tthmax'])
        return
    
    def _additionalUpdataSelf(self, **kwargs):
        '''additional process called in self._updateSelf, this method is called before 
        self._copySelftoConfig()
        '''
        #test if tthmax or qmax has changed
        self._checkMax()
        #self._checkStep()
        if self.integrationspace == 'twotheta':
            self.tthorqmax = self.tthmax
            self.tthorqstep = self.tthstep
        elif self.integrationspace == 'qspace':
            self.tthorqmax = self.qmax
            self.tthorqstep = self.qstep
        return
    
    def _additionalPostProcessing(self, **kwargs):
        self._loadFromFit2D(self.fit2dconfig)
        if self.createconfig!='' or self.createconfigfull!='':
            self.nocalculation = True
        if self.createmask!='':
            self.nocalculation = True
        return
    
    def _loadFromFit2D(self, filename):
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
                self._updateSelf()
        return
    
    def _checkMax(self):
        '''check and change tthmax and qmax to actual tthmax and qmax of image 
        '''
        xr = (np.array([0, self.xdimension+1]) - self.xbeamcenter) * self.xpixelsize
        yr = (np.array([0, self.ydimension+1]) - self.ybeamcenter) * self.ypixelsize
        sinr = np.sin(self.rotation)
        cosr = np.cos(self.rotation)
        sint = np.sin(self.tilt)
        cost = np.cos(self.tilt)
        sourcexr = self.distance * sint * cosr
        sourceyr = -self.distance * sint * sinr
        sourcezr = self.distance * cost
        
        dmatrix = ((xr - sourcexr) ** 2).reshape(1, 2) + \
                  ((yr - sourceyr) ** 2).reshape(2, 1) + sourcezr ** 2
        dmatrix = np.sqrt(dmatrix)
        tthmatrix1 = ((-xr + sourcexr) * sourcexr).reshape(1, 2) + \
                     ((-yr + sourceyr) * sourceyr).reshape(2, 1) + sourcezr * sourcezr
        tthmatrix = np.arccos(tthmatrix1 / dmatrix / self.distance)
        qmatrix = 4 * np.pi * np.sin(tthmatrix / 2.0) / self.wavelength
        
        self.tthmaxd = np.degrees(np.max(tthmatrix)) + 1.0
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

if __name__=='__main__':
    a = SrXPlanarConfig()
    a.loadFromFile('test.cfg')
    a.updateConfig()
    a.writeConfig('test.cfg')
