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
# See LICENSE.txt for license information.
#
##############################################################################

import numpy as np
import ConfigParser
import re, os, sys
from functools import partial
import argparse

from diffpy.confutils.config import ConfigBase
from diffpy.confutils.tools import _configPropertyRad, _configPropertyR, _configPropertyRW

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

_optdatalist = [
        # control group
        ['filenames', {'sec':'Control', 'config':'n', 'header':'n',
            'f':'filename',
            'h':'filename or list of filenames or filename pattern or list of filename pattern',
            'n':'*',
            'd':[], }],
        ['output', {'sec':'Experiment', 'config':'n', 'header':'n',
            's':'o',
            'h':'basename of output file',
            'd':'', }],
        ['summation', {'sec':'Control', 'config':'n', 'header':'n',
            's':'s',
            'h':'sum all the image and then integrate',
            'n':'?',
            'co':True,
            'd':False, }],
        # Expeiment gropu
        ['fit2dconfig', {'sec':'Experiment', 'config':'n', 'header':'n',
            's':'fit2d',
            'h':'fit2d calibration file name. It contains the calibration results copy from fit2d cmd windows',
            'd':'', }],
        ['opendirectory', {'sec':'Control', 'header':'n',
            's':'opendir',
            'h':'directory of input 2D image files',
            'd':'currentdir',
            'tt':'directory'}],
        ['savedirectory', {'sec':'Control', 'header':'n',
            's':'savedir',
            'h':'directory of output files',
            'd':'currentdir',
            'tt':'directory'}],
        ['addmask', {'sec':'Experiment',
            'h':
'''list of masks to apply on the image. Specify filename of mask file (support .msk, .npy, and .tif).
add "edgemask" to mask the pixels at edges, add "deadpixel" to mask dead pixels, add "brightpixel"\
to mask bright pixels.''',
            'n':'*',
            # 'd':['edgemask'], }],
            'd':[], }],
        ['createmask', {'sec':'Control', 'config':'n', 'header':'n',
            'h':'create a mask file according to current image file and value of addmask',
            'd':'', }],
        ['integrationspace', {'sec':'Experiment',
            'h':'the x-grid of integrated 1D diffraction data',
            'd':'twotheta',
            'c':['qspace', 'twotheta'], }],
        ['wavelength', {'sec':'Experiment',
            'h':'wavelength of x-ray, in Angstrom',
            'd':0.1000, }],
        ['xbeamcenter', {'sec':'Experiment',
            's':'xc',
            'h':'beamcenter in x axis, in pixel',
            'd':1024.0, }],
        ['ybeamcenter', {'sec':'Experiment',
            's':'yc',
            'h':'beamcenter in y axis, in pixel',
            'd':1024.0, }],
        ['distance', {'sec':'Experiment',
            's':'dis',
            'h':'distance between detector and sample, in mm',
            'd':200.0, }],
        ['rotationd', {'sec':'Experiment',
            's':'rot',
            'h':'rotation angle of tilt plane, in degree',
            'd':0.0, }],
        ['tiltd', {'sec':'Experiment',
            's':'tilt',
            'h':'tilt angle of tilt plane, in degree',
            'd':0.0, }],
        ['tthstepd', {'sec':'Experiment',
            's':'ts',
            'h':'integration step in twotheta space, in degree',
            'd':0.02, }],
        ['qstep', {'sec':'Experiment',
            's':'qs',
            'h':'integration step in q space, in Angstrom^-1',
            'd':0.02, }],
        # Beamline group
        ['includepattern', {'sec':'Beamline', 'header':'n', 'config':'f',
            's':'ipattern',
            'h':'list of string, file name patterns for included files',
            'n':'*',
            'd':['*.tif', '*.tif.bz2'], }],
        ['excludepattern', {'sec':'Beamline', 'header':'n', 'config':'f',
            's':'epattern',
            'h':'list of string, file name patterns for excluded files',
            'n':'*',
            'd':['*.dark.tif', '*.raw.tif'], }],
        ['fliphorizontal', {'sec':'Beamline',
            'h':'filp the image horizontally',
            'n':'?',
            'co':True,
            'd':False, }],
        ['flipvertical', {'sec':'Beamline',
            'h':'filp the image vertically',
            'n':'?',
            'co':True,
            'd':True, }],
        ['xdimension', {'sec':'Beamline',
            's':'xd',
            'h':'detector dimension in x axis, in pixel',
            'd':2048, }],
        ['ydimension', {'sec':'Beamline',
            's':'yd',
            'h':'detector dimension in y axis, in pixel',
            'd':2048, }],
        ['xpixelsize', {'sec':'Beamline',
            's':'xp',
            'h':'detector pixel size in x axis, in mm',
            'd':0.2, }],
        ['ypixelsize', {'sec':'Beamline',
            's':'yp',
            'h':'detector pixel size in y axis, in mm',
            'd':0.2, }],
        # Others Group
        ['uncertaintyenable', {'sec':'Others',
            's':'error',
            'h':'enable uncertainty propagation',
            'n':'?',
            'co':True,
            'd':True, }],
        ['sacorrectionenable', {'sec':'Others',
            's':'sacorr',
            'h':'enable solid angle correction',
            'n':'?',
            'co':True,
            'd':True, }],
        ['polcorrectionenable', {'sec':'Others',
            's':'polarcorr',
            'h':'enable polarization correction',
            'n':'?',
            'co':True,
            'd':True, }],
        ['polcorrectf', {'sec':'Others',
            's':'polarf',
            'h':'polarization correction factor',
            'd':0.99, }],
        ['gsasoutput', {'sec':'Others', 'header':'n',
            'h':'select if want to output gsas format file',
            'c':['None', 'std', 'esd', 'fxye'],
            'd':'None', }],
        ['filenameplus', {'sec':'Others', 'header':'n',
            'h':'string appended to the output filename',
            'd':'', }],
        ['maskedges', {'sec':'Others',
            'h':'mask the edge pixels, first four means the number of pixels masked in each edge \
(left, right, top, bottom), the last one is the radius of a region masked around the corner',
            'n':5,
            'tt':'array',
            't':'intlist',
            'd':[20, 20, 20, 20, 100], }],
        ['cropedges', {'sec':'Others',
            'h':'crop the edge pixels, first four means the number of pixels masked in each edge \
(left, right, top, bottom)',
            'n':4,
            'tt':'array',
            't':'intlist',
            'd':[20, 20, 20, 20], }],
        ['nocalculation', {'sec':'Others', 'config':'n', 'header':'n',
            'h':'set True to disable all calculation, will automaticly set True if createconfig or createmask',
            'n':'?',
            'co':True,
            'd':False, }],
        ]

_defaultdata = {'configfile': ['srxplanar.cfg', 'SrXplanar.cfg'],
                'headertitle': 'SrXplanar configration'
                }

class SrXplanarConfig(ConfigBase):
    '''
    config class, based on ConfigBase class in diffpy.confutils
    '''

    # Text to display before the argument help
    _description = _description

    # Text to display after the argument help
    _epilog = _epilog

    _optdatalist = _optdatalist

    _defaultdata = _defaultdata

    def _preInit(self, **kwargs):
        '''
        method called in init process, overload it!
        
        this method will be called before reading config from file/args/kwargs
        
        add degree/rad delegation for rotation, tilt, tthstep, tthmax
        '''

        for name in ['rotation', 'tilt', 'tthstep', 'tthmax']:
            setattr(self.__class__, name, _configPropertyRad(name + 'd'))
        # cls._configlist['Experiment'].extend(['rotation', 'tilt', 'tthstep', 'tthmax'])
        return

    def _preUpdateSelf(self, **kwargs):
        '''
        additional process called in self._updateSelf, this method is called
        before self._copySelftoConfig(), i.e. before copy options value to
        self.config (config file)
        
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

    # def _postUpdateConfig(self, nofit2d=False, **kwargs):
    def _postUpdateConfig(self, **kwargs):
        '''
        post processing after parse args or kwargs, this method is called after 
        in self._postPocessing and before creating config file action  
        
        load fit2d config if specified in config, and set nocalculatio flag when create 
        config or create mask
        
        :param kwargs: optional kwargs
        '''

        self._loadFromFit2D(self.fit2dconfig)
        if (self.createconfig != '')and(self.createconfig != None):
            self.nocalculation = True
        if (self.createconfigfull != '')and(self.createconfigfull != None):
            self.nocalculation = True
        if self.createmask != '':
            self.nocalculation = True
        return

    def _loadFromFit2D(self, filename):
        '''
        load parameters from fit2d calibration information. copy/paste the fit2d calibration 
        results to a txt file. this function will load xbeamcenter, ybeamceter... from the file
        
        :param filename: file name of fit2d result file
        '''
        rv = parseFit2D(filename)
        if len(rv.values()) > 0:
            for optname in rv.keys():
                setattr(self, optname, rv[optname])
            self.fit2dconfig = ''
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


    xr = (np.array([0, xdimension + 1]) - xbeamcenter) * xpixelsize
    yr = (np.array([0, ydimension + 1]) - ybeamcenter) * ypixelsize
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

SrXplanarConfig.initConfigClass()

if __name__ == '__main__':
    a = SrXplanarConfig()
    a.updateConfig()
    a.writeConfig('test.cfg')
