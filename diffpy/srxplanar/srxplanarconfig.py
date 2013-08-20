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
import ConfigParser
import re, os, sys
from functools import partial
import argparse

from diffpy.confutils.config import ConfigBase, initConfigClass
from diffpy.confutils.tools import _configPropertyRad, _configPropertyR, _configPropertyRW

class SrXplanarConfig(ConfigBase):
    '''
    config class, based on ConfigBase class in diffpy.confutils
    '''
    
    # Text to display before the argument help
    _description = \
    '''
    SrXplanar -- integrate 2D powder diffraction image to 1D with unceratinty propagation
    ''' 
    # Text to display after the argument help
    _epilog = \
    '''
Examples:

    srxplanar KFe2As2-00838.tif -c test.cfg
    --integration using config file test.cfg

    srxplanar KFe2As2-00838.tif -fit2d fit2d.txt
    --integration using calibration information in fit2d.txt
    
    srxplanar *.tif -c test.cfg -s
    --integration all .tif image and sum them into one
    
    srxplanar KFe2As2-00838.tif -fit2d fit2d.txt --integrationspace twotheta
    -- integrate using calibration information in fit2d.txt and integrate into two theta space
    
    srxplanar --createconfig config.cfg
    --create default (short) config file using all default value
    
    srxplanar --createconfigfull configfull.cfg -fit2d fit2d.txt
    --create a complete config file using calibration information in fit2d.txt 
    '''
    
    # optdata contains these keys:
    # full(f), short(s), help(h), type(t), action(a), nargs(n), default(d), choices(c), required(r), dest, const
    
    '''optional args: 
        'args': default is 'a'
            if 'a', this option will be available in self.args
            if 'n', this option will not be available in self.args
        'config': default is 'a'
            if 'f', this option will present in self.config and be written to config file only in full mode
            if 'a', this option will present in self.config and be written to config file both in full and short mode
            if 'n', this option will not present in self.config
        'header', default is 'a'
            if 'f', this option will be written to header only in full mode
            if 'a', this option will be written to header both in full and short mode
            if 'n', this option will not be written to header
        so in short mode, all options with 'a' will be written, in full mode, all options with 'a' or 'f' will be written
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
            'd':['edgemask'],}],
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
            'd':['*.tif', '*.tif.bz2'],}],
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
            'd':True,}],
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
            'd':[20,20,20,20,100],}],
        ['nocalculation',{'sec':'Others','config':'n','header':'n',
            'h':'set True to disable all calculation, will automaticly set True if createconfig or createmask',
            'n':'?',
            'co':True,
            'd':False,}],
        ]
    
    #default config file path and name
    _defaultconfigpath = ['srxplanar.cfg', 'SrXplanar.cfg']
    
    #default first list for header
    _defaultheaderline = 'SrXplanar configration'
    
    def _additionalInit(self):
        '''
        method called in init process
        this method will be called after all options in self._optdata are processed, i.e. all options are created. 
        and before reading config from file/args/kwargs
        
        add degree/rad delegation for rotation, tilt, tthstep, tthmax
        '''
        for name in ['rotation', 'tilt', 'tthstep', 'tthmax']:
            setattr(self.__class__, name, _configPropertyRad(name+'d'))
        self._configlist['Experiment'].extend(['rotation', 'tilt', 'tthstep', 'tthmax'])
        return
    
    def _additionalUpdataSelf(self, **kwargs):
        '''
        additional process called in self._updateSelf, this method is called before 
        self._copySelftoConfig(), i.e. before copy options value to self.config (config file)
        
        check the tthmaxd and qmax, and set tthorqmax, tthorqstep according to integration space
        
        :param kwargs: optional kwargs
        '''
        self.tthmaxd, self.qmax = checkMax(self)
        if self.integrationspace == 'twotheta':
            self.tthorqmax = self.tthmax
            self.tthorqstep = self.tthstep
        elif self.integrationspace == 'qspace':
            self.tthorqmax = self.qmax
            self.tthorqstep = self.qstep
        return
    
    def _additionalPostProcessing(self, nofit2d=False, **kwargs):
        '''
        post processing after parse args or kwargs, this method is called after 
        in self._postPocessing and before creating config file action  
        
        load fit2d config if specified in config, and set nocalculatio flag when create 
        config or create mask
        
        :param nofit2d: boolean, if True, it will skip loading fit2d calibration, this is useful
            when you reload/update some parameters but don't want to reload fit2d calibration 
            results.
        :param kwargs: optional kwargs
        '''
        if not nofit2d:
            self._loadFromFit2D(self.fit2dconfig)
        if self.createconfig!='' or self.createconfigfull!='':
            self.nocalculation = True
        if self.createmask!='':
            self.nocalculation = True
        return
    
    def _loadFromFit2D(self, filename):
        '''
        load parameters from fit2d calibration information. copy/paste the fit2d calibration 
        results to a txt file. this function will load xbeamcenter, ybeamceter... from the file
    
        :param filename: str, file name (with full path if not in current dir) of fit2d file,
            or a string containing the calibraiton parameters copy from fit2d.
        '''
        rv = parseFit2D(filename)
        for optname in rv.keys():
            setattr(self, optname, rv[optname])
        self._updateSelf()
        return

def checkMax(config):
    '''
    calculate the max twotheta angle (and q) of a detector with current geometry
    
    :param config: SrXplanarConfig, config instance stores the geometry parameters
    
    :return: [tthmaxd, qmax], max twotheta angle(in degree) and max q value of current
        detector.
    '''
    xdimension = getattr(config, 'xdimension')
    ydimension = getattr(config, 'ydimension')
    xbeamcenter = getattr(config, 'xbeamcenter')
    ybeamcenter = getattr(config, 'ybeamcenter')
    xpixelsize = getattr(config, 'xpixelsize')
    ypixelsize = getattr(config, 'ypixelsize')
    rotation = getattr(config, 'rotation')
    tilt = getattr(config, 'tilt')
    distance = getattr(config, 'distance')
    wavelength = getattr(config, 'wavelength')
    
    
    xr = (np.array([0, xdimension+1]) - xbeamcenter) * xpixelsize
    yr = (np.array([0, ydimension+1]) - ybeamcenter) * ypixelsize
    sinr = np.sin(rotation)
    cosr = np.cos(rotation)
    sint = np.sin(tilt)
    cost = np.cos(tilt)
    sourcexr = distance * sint * cosr
    sourceyr = -distance * sint * sinr
    sourcezr = distance * cost
    
    dmatrix = ((xr - sourcexr) ** 2).reshape(1, 2) + \
              ((yr - sourceyr) ** 2).reshape(2, 1) + sourcezr ** 2
    dmatrix = np.sqrt(dmatrix)
    tthmatrix1 = ((-xr + sourcexr) * sourcexr).reshape(1, 2) + \
                 ((-yr + sourceyr) * sourceyr).reshape(2, 1) + sourcezr * sourcezr
    tthmatrix = np.arccos(tthmatrix1 / dmatrix / distance)
    qmatrix = 4 * np.pi * np.sin(tthmatrix / 2.0) / wavelength
    
    tthmaxd = np.degrees(np.max(tthmatrix)) + 0.5
    qmax = np.max(qmatrix) + 0.1
    return tthmaxd, qmax

def parseFit2D(filename):
    '''
    load parameters from fit2d calibration information. copy/paste the fit2d calibration 
    results to a txt file. this function will load xbeamcenter, ybeamceter... from the file
    
    :param filename: str, file name (with full path if not in current dir) of fit2d file,
        or a string containing the calibraiton parameters copy from fit2d.
        
    :return: dict, including 'xbeamcenter', 'ybeamcenter', 'wavelength', 'rotationd'(angle of ratation), 
        'tiltd'(angle of tilt rotation)
    '''
    rv = {}
    def findFloat(line):
        temp = re.findall('[-+]?\d*\.\d+|[-+]?\d+', line)
        return map(float, temp)
    if filename != None:
        if os.path.exists(filename):
            f = open(filename, 'r')
            lines = f.readlines()
            f.close()
        else:
            lines = filename.split() 
        for line in lines:
            if re.search('Refined Beam centre.*pixels', line):
                rv['xbeamcenter'], rv['ybeamcenter'] = findFloat(line)
            elif re.search('Refined sample to detector distance', line):
                rv['distance'] = findFloat(line)[0]
            elif re.search('Refined wavelength', line):
                rv['wavelength'] = findFloat(line)[0]
            elif re.search('Refined tilt plane rotation angle', line):
                rv['rotationd'] = findFloat(line)[0]
            elif re.search('Refined tilt angle', line):
                rv['tiltd'] = findFloat(line)[0]
    return rv

initConfigClass(SrXplanarConfig)

if __name__=='__main__':
    a = SrXplanarConfig()
    a.updateConfig()
    a.writeConfig('test.cfg')
