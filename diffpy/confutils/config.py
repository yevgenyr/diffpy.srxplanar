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

'''for python 2.6, argparse is required
'''


import numpy as np
import ConfigParser
import re, os, sys
from functools import partial
import argparse

def _configPropertyRad(nm):
    '''helper function of options delegation, rad 2 degree'''
    rv = property(fget = lambda self: np.radians(getattr(self, nm)), 
                  fset = lambda self, val: setattr(self, nm, np.degrees(val)), 
                  fdel = lambda self: delattr(self, nm))
    return rv


class ConfigBase(object):
    '''Class used for storing the configuration value and process cmd command
    config: ConfigParser.ConfigParser
    args: argparse
    '''
    config = ConfigParser.ConfigParser()
    args = argparse.ArgumentParser(description='SrSig2d Configuration')
    
    # optdata contains these keys:
    # full(f), short(s), help(h), type(t), action(a), nargs(n), default(d), choices(c), required(r), dest, const
    _optdatanamedic = {'h':'help',
                       't':'type',
                       'a':'action',
                       'n':'nargs',
                       'd':'default',
                       'c':'choices',
                       'r':'required',
                       'dest':'dest',
                       'const':'const'}
    #optional args: 'ctrl':[True, True], Set False to hide thie options in [self.config, self.args]
    
    #examples, overload it
    _optdatalist = [
        ['fit2dconfig',{'sec':'Experiment',
            's':'fit2d',
            'h':'fit2d calibration file name. It contains the calibration results copy from fit2d cmd windows',
            'd':'',}],
        ['integrationspace',{'sec':'Experiment',
            'h':'integration space, could be twotheta or qspace',
            'd':'twotheta',
            'c':['twotheta','qspace'],}],
        ['wavelength',{'sec':'Experiment',
            'h':'wavelength of x-ray, in A',
            'd':0.1000,}],
        ['rotationd',{'sec':'Experiment',
            's':'rot',
            'h':'rotation angle of tilt plane, in degree',
            'd':0.0,}],
        ['includepattern',{'sec':'Beamline',
            's':'ipattern',
            'h':'file name pattern for included files',
            'd':'*.tif',}],
        ['excludepattern',{'sec':'Beamline',
            's':'epattern',
            'h':'file name pattern for excluded files',
            'n':'*',
            'd':['*.dark.tif', '*.raw.tif'],}],
        ['fliphorizontal',{'sec':'Beamline',
            'h':'filp the image horizontally',
            'd':False,}],
        ['regulartmatrixenable',{'sec':'Others',
            'h':'normalize tmatrix in splitting method',
            'd':False,}],
        ['maskedges',{'sec':'Others',
            'h':'mask the edge pixels, first four means the number of pixels masked in each edge \
                (left, right, top, bottom), the last one is the radius of a region masked around the corner',
            'n':5,
            'd':[1,1,1,1,50],}],
        ]
    
    _optdata = dict(_optdatalist)
    #options that consist of list of value
    #examples, overload it
    _listoptdata = {
        'strlist':['excludepattern'],
        'floatlist':[],
        'intlist':['maskedges'],
        'strlist':[],
        }
    #configlist, store the options name for each sections
    _configlist = {} 
    
    #example, overload it
    def __init__(self, filename=None, args=None, **kwargs):
        self._detectAddSections()
        for opt in self._optdatalist:
            key = opt[0]
            self._addOpt(key)
        #degree to rad
        
        #updata config
        self.updateConfig(filename, args, **kwargs)
        return
    
    #example, overload it
    def _additionalInit(self):
        '''method called in init process, overload it
        this method will be called after all options in self._optdata are processced and before reading config from file/args/kwargs
        '''
        for name in ['rotation']:
            setattr(self.__class__, name, _configPropertyRad(name+'d'))
        self._configlist['Experiment'].extend(['rotation'])
        return
    
    #example, overload it
    def _updateSelf(self, optnames=None):
        '''update the options value, then copy the values in the self.'options' to self.config
        param optnames: str or list of str, name of options whose value has been changed, if None, update all options
        '''
        #so some check right here
        pass
        #copy value to self.config
        self._copySelftoConfig(optnames)
        return
    
    def _getType(self, optname):
        '''detect the type of option
        param optname: str, name of option
        
        return: string of type 
        '''
        optdata = self._optdata[optname]
        if optdata.has_key('t'):
            opttype = optdata['t']
        else:
            value = optdata['d']
            if type(value)==str:
                opttype = 'str'
            elif type(value)==float:
                opttype = 'float'
            elif type(value)==int:
                opttype = 'int'
            elif type(value)==bool:
                opttype = 'bool'
            elif type(value)==list:
                if type(value[0])==str:
                    opttype = 'strlist'
                elif type(value[0])==float:
                    opttype = 'floatlist'
                elif type(value[0])==int:
                    opttype = 'intlist'
                elif type(value[0])==bool:
                    opttype = 'boollist'
        return opttype
    
    def _detectAddSections(self):
        '''detect sections present in self._optdata and add them to self.config
        also add it to self._configlist
        '''
        seclist = [self._optdata[key]['sec'] for key in self._optdata.keys()]
        for sec in set(seclist):
            self.config.add_section(sec)
            self._configlist[sec] = []
        return
    
    def _addOpt(self, optname):
        '''add options to self.config and self.args and self.'optname'
        
        param optname: string, name of option
        '''
        optdata = self._optdata[optname]
        opttype = self._getType(optname)
        if opttype.endswith('list'):
            self.config.set(optdata['sec'], optname, str(optdata['d'])[1:-1])
        else:
            self.config.set(optdata['sec'], optname, optdata['d'])
        #add to self.configlist
        self._configlist[optdata['sec']].append(optname)
        #add to self.'optname'
        setattr(self, optname, optdata['d'])
        #transform optdata to a dict that can pass to add_argument method
        pargs = dict()
        for key in optdata.keys():
            if self._optdatanamedic.has_key(key):
                pargs[self._optdatanamedic[key]] = optdata[key]
        #replace currentdir in default to os.getcwd()
        if pargs['default'] =='currentdir':
            pargs['default'] == os.getcwd()
        #add args
        if optdata.has_key('s'):
            self.args.add_argument('--'+optname, '-'+optdata['s'], **pargs)
        else:
            self.args.add_argument('--'+optname, **pargs)
        return
    
    
    def _copyConfigtoSelf(self, optnames=None):
        '''copy the values in self.config to self.'options'
        
        param optname: str or list of str, names of options whose value copied from self.config to self.'optname'. Set None to update all
        '''
        if not hasattr(self, 'mdict'):
            self.mdict = {
                'str': self.config.get,
                'int': self.config.getint,
                'float': self.config.getfloat,
                'bool': self.config.getboolean,
                'strlist': partial(self._getListValue, vtype='strlist'),
                'intlist': partial(self._getListValue, vtype='intlist'),
                'floatlist': partial(self._getListValue, vtype='floatlist'),
                'boollist': partial(self._getListValue, vtype='boollist'),
                }
        if optnames!=None:
            optnames = optnames if type(optnames)==list else [optnames]
            for optname in optnames:
                if self._optdata.has_key(optname):
                    secname = self._optdata[optname]['sec']
                    opttype = self._getType(optname)
                    setattr(self, optname, self.mdict[opttype](secname, optname))
        else:
            for secname in self.config.sections():
                for optname in self.config.options(secname):
                    if self._optdata.has_key(optname):
                        opttype = self._getType(optname)
                        setattr(self, optname, self.mdict[opttype](secname, optname))
        return
    
    def _copySelftoConfig(self, optnames=None):
        '''copy the value in self.'options' to self.config
        
        param optname: str or list of str, names of options whose value copied from self.'optname' to self.config. Set None to update all
        '''
        if optnames!=None:
            optnames = optnames if type(optnames)==list else [optnames]
            for optname in optnames:
                if self._optdata.has_key(optname):
                    secname = self._optdata[optname]['sec']
                    opttype = self._getType(optname)
                    if opttype.endswith('list'):
                        liststr = map(str, getattr(self, optname))
                        liststring = ', '.join(liststr)
                        self.config.set(secname, optname, liststring)
                    else:
                        self.config.set(secname, optname, str(getattr(self, optname)))
        else:
            for secname in self.config.sections():
                for optname in self.config.options(secname):
                    if self._optdata.has_key(optname):
                        opttype = self._getType(optname)
                        if opttype.endswith('list'):
                            liststr = map(str, getattr(self, optname))
                            liststring = ', '.join(liststr)
                            self.config.set(secname, optname, liststring)
                        else:
                            self.config.set(secname, optname, str(getattr(self, optname)))
        return
    
    def _getListValue(self, sectionname, optionname, vtype):
        '''helper function that get a list of value for string read from self.config
        
        param sectionname: string, name of section in self.config
        param optionname: string, name of option in self.config
        vtype: string, type of list, can be strlist, intlist, floatlist, boollist
        
        return: list of value in desired format, or [] if options value is empty
        ''' 
        temp = re.split('\s*,\s*', self.config.get(sectionname, optionname))
        if len(temp)>0:
            if vtype.startswith('str'):
                rv = temp
            elif  vtype.startswith('int'):
                rv = map(int, temp)
            elif  vtype.startswith('float'):
                rv = map(float, temp)
            elif  vtype.startswith('bool'):
                rv = map(lambda s: (s.lower()=='true')or(s.lower()=='yes'), temp)
        else:
            rv = []
        return rv
    
    def parseArgs(self, pargs):
        '''parse args and update the value in self.'optname', this will call the self.args to parse args,
        
        param pargs: list of string, pargs to parse, usually be sys.argv
        '''
        changedargs = filter(lambda args: args.startswith('-'), pargs)
        changedargs = map(lambda arg: arg.replace('-', ''), changedargs)
        if len(changedargs)>0:
            obj = self.args.parse_args(pargs)
            for optname in changedargs:
                if self._optdata.has_key(optname):
                    setattr(self, optname, getattr(obj, optname))
            #update self
            self._updateSelf(changedargs)
        return
    
    def parseKwargs(self, **kwargs):
        '''update self.'optname' values according to the kwargs
        '''
        if kwargs:
            changedargs = []
            for optname, optvalue in kwargs.iteritems():
                if self._optdata.has_key(optname):
                    setattr(self, optname, optvalue)
                    changedargs.append(optname)
            #update self
            self._updateSelf(changedargs)
        return

    def parseConfigFile(self, filename):
        '''read a config file and update the self.'optname'
        '''
        if filename!=None:
            if os.path.exists(filename):
                self.config.read(filename)
                self._copyConfigtoSelf()
                self._updateSelf()
                self._loadFromFit2D(self.fit2dconfig)
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
    
    def updateConfig(self, filename=None, args=None, **kwargs):
        '''update config according to config file, args(from sys.argv) or **kwargs
        
        param filename: string, name of config file
        param args: list of str, usually be sys.argv
        param **kwargs: you can use like 'xbeamcenter=1024' to update the value of xbeamcenter
        '''
        if filename!=None:
            self.parseConfigFile(filename)
        if args!=None:
            self.parseArgs(args)
        if kwargs!={}:
            self.parseKwargs(**kwargs)
        if (filename==None)and(args==None)and(kwargs==None):
            self._updateSelf()
        return
    
    def writeConfig(self, filename):
        '''write config to file
        
        param filename: string, name of file
        '''
        self._updateSelf()
        conffile = open(filename, 'w')
        self.config.write(conffile)
        conffile.close()
        return

# Helper Routines ------------------------------------------------------

def _configPropertyR(name):
    '''Create a property that forwards self.name to self.config.name.
    '''
    rv = property(fget = lambda self: getattr(self.config, name),
            doc='attribute forwarded to self.config, read-only')
    return rv

def _configPropertyRW(name):
    '''Create a property that forwards self.name to self.config.name.
    '''
    rv = property(fget = lambda self: getattr(self.config, nm), 
                  fset = lambda self, value: setattr(self.config, nm, value),
                  fdel = lambda self: delattr(self, nm),
                  doc='attribute forwarded to self.config, read/write')
    return rv
