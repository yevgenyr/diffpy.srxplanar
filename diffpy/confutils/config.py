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

'''package for organizing program configurations. It can read/write configurations file, 
parse arguments from command lines, and also parse arguments passed from method/function calling
inside python.

Note: for python 2.6, argparse and orderedDict is required, install them with easy_install
'''

import numpy as np
import ConfigParser
import re, os, sys
from functools import partial
import argparse
if sys.version_info < (2,7):
    from ordereddict import OrderedDict
else:
    from collections import OrderedDict

from diffpy.confutils.tools import _configPropertyRad, _configPropertyR, _configPropertyRW

class ConfigBase(object):
    '''_optdatalistpre, _optdatalist, _optdatalistext are metadata used to initialize the options, see below for examples
    
    options presents in --help (in cmd), config file, and headers have the same order in there three list, so arrange them 
    in right order here. 
    
    optional args to control if the options presents in args, config file or file header: 
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
    
    '''optdata contains these keys:
    full(f, positional), short(s), help(h), type(t), action(a), nargs(n), default(d), choices(c), required(r), dest, const
    used in argparse'''
    _optdatanamedic = {'h':'help',
                       't':'type',
                       'a':'action',
                       'n':'nargs',
                       'd':'default',
                       'c':'choices',
                       'r':'required',
                       'de':'dest',
                       'co':'const'}
    
    #examples, overload it
    _optdatalistpre = []
    
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
    #examples, overload it
    _optdatalistext = [
        ['tifdirectory',{'sec':'Experiment', 'header':'n',
            's':'tifdir',
            'h':'directory of raw tif files',
            'd':'currentdir',}],
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
            'n':'?',
            'co':True,
            'd':False,}],
        ['regulartmatrixenable',{'sec':'Others',
            'h':'normalize tmatrix in splitting method',
            'n':'?',
            'co':True,
            'd':False,}],
        ['maskedges',{'sec':'Others',
            'h':'mask the edge pixels, first four means the number of pixels masked in each edge \
                (left, right, top, bottom), the last one is the radius of a region masked around the corner',
            'n':5,
            'd':[1,1,1,1,50],}],
        ]
    
    #default config file path and name, overload it for your config class
    _defaultconfigpath = ['config.cfg']
    
    def __init__(self, filename=None, args=None, **kwargs):
        '''
        param filename: str, file name of config file
        param filename: list of str, args passed from cmd
        param kwargs: optional kwargs
        '''
        self._additionalInit()
        #updata config, first detect if a default config should be load
        filename = self._findDefaultConfigFile(filename, args, **kwargs)
        self.updateConfig(filename, args, **kwargs)
        return
    
    #example, overload it
    def _additionalInit(self):
        '''method called in init process
        this method will be called after all options in self._optdata are processed, i.e. all options are created. 
        and before reading config from file/args/kwargs
        '''
        for name in ['rotation']:
            setattr(self.__class__, name, _configPropertyRad(name+'d'))
        self._configlist['Experiment'].extend(['rotation'])
        return
    
    #example, overload it
    def _additionalPostProcessing(self, **kwargs):
        '''post processing after parse args or kwargs, this method is called after 
        in self._postPocessing and before creating config file action  
        '''
        #Overload it!
        return
    
    #example, overload it
    def _additionalUpdataSelf(self, **kwargs):
        '''additional process called in self._updateSelf, this method is called before 
        self._copySelftoConfig(), i.e. before copy options value to self.config (config file)
        '''
        return
    
    ##########################################################################
    
    def _updateSelf(self, optnames=None, **kwargs):
        '''update the options value, then copy the values in the self.'options' to self.config
        param optnames: str or list of str, name of options whose value has been changed, if None, update all options
        '''
        #so some check right here
        self._additionalUpdataSelf(**kwargs)
        #copy value to self.config
        self._copySelftoConfig(optnames)
        return
    
    
    def _findDefaultConfigFile(self, filename=None, args=None, **kwargs):
        '''find default config file, if any config is specified in filename/args/kwargs
        then return the filename of config. 
        
        param filename: str, file name of config file
        param filename: list of str, args passed from cmd
        param kwargs: optional kwargs
        
        return: name of config file if found, otherwise None 
        '''
        rv = None
        flag = False
        if (filename!=None):
            flag == True
            rv = filename
        if (args!=None):
            if ('--configfile' in args) or ('-c' in args):
                flag = True
        if kwargs.has_key('configfile'):
                flag = True
        if not flag:
            for dconf in self._defaultconfigpath:
                if (os.path.exists(dconf))and(rv==None):
                    rv = dconf
        return rv
            
    
    def _getType(self, optname):
        '''detect the type of option
        param optname: str, name of option
        
        return: string of type 
        '''
        opttype =self._getTypeC(optname)
        return opttype
    
    @classmethod
    def _getTypeC(cls, optname):
        '''detect the type of option
        param optname: str, name of option
        
        return: string of type 
        '''
        optdata = cls._optdata[optname]
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
                if len(value)==0:
                    opttype = 'strlist'
                elif type(value[0])==str:
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
        self._detectAddSectionsC(self)
        return
    
    @classmethod
    def _detectAddSectionsC(cls):
        '''detect sections present in self._optdata and add them to self.config
        also add it to self._configlist
        param cls: class, class to configurate
        '''
        #seclist = [self._optdata[key]['sec'] for key in self._optdata.keys()]
        seclist = [cls._optdata[opt[0]]['sec'] for opt in cls._optdatalist]
        secdict = OrderedDict.fromkeys(seclist)
        #for sec in set(seclist):
        for sec in secdict.keys():
            cls.config.add_section(sec)
            cls._configlist[sec] = []
        return
    
    def _addOpt(self, optname):
        '''add options to self.config and self.args and self.'optname',
        this will read metadata in self._optdatalist
        
        param optname: string, name of option
        '''
        self._addOptC(self, optname)
        return
    
    @classmethod
    def _addOptC(cls, optname):
        '''add options to self.config and self.args and self.'optname',
        this will read metadata in self._optdatalist
        
        param optname: string, name of option
        '''
        optdata = cls._optdata[optname]
        opttype = cls._getTypeC(optname)
        
        #replace currentdir in default to os.getcwd()
        if optdata['d'] =='currentdir':
            optdata['d'] = os.getcwd()
        
        #add to cls.'optname'
        setattr(cls, optname, optdata['d'])
        
        #add to cls.config
        secname =  optdata['sec'] if optdata.has_key('sec') else 'Others'
        cls._configlist[secname].append(optname)
        if optdata.get('config', 'a')!='n':
            strvalue = ', '.join(map(str, optdata['d'])) if type(optdata['d'])==list else str(optdata['d'])
            cls.config.set(secname, optname, strvalue)
        #add to cls.args
        if optdata.get('args', 'a')!='n':
            #transform optdata to a dict that can pass to add_argument method
            pargs = dict()
            for key in optdata.keys():
                if cls._optdatanamedic.has_key(key):
                    pargs[cls._optdatanamedic[key]] = optdata[key]
            pargs['default'] = argparse.SUPPRESS
            #add args
            if optdata.has_key('f'):
                cls.args.add_argument(optname, **pargs)
            elif optdata.has_key('s'):
                cls.args.add_argument('--'+optname, '-'+optdata['s'], **pargs)
            else:
                cls.args.add_argument('--'+optname, **pargs)
        return
    
    def _copyConfigtoSelf(self, optnames=None):
        '''copy the values in self.config to self.'options'
        
        param optname: str or list of str, names of options whose value copied from self.config to self.'optname'. 
        Set None to update all
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
        
        param optname: str or list of str, names of options whose value copied from self.'optname' to self.config. 
        Set None to update all
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
        obj = self.args.parse_args(pargs)
        changedargs = obj.__dict__.keys()
        for optname in changedargs:
            if self._optdata.has_key(optname):
                setattr(self, optname, getattr(obj, optname))
        #update self
        if len(changedargs)>0:
            self._updateSelf(changedargs)
        return obj
    
    def parseKwargs(self, **kwargs):
        '''update self.'optname' values according to the kwargs
        
        param kwargs: keywords=value
        '''
        if kwargs!={}:
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
        
        param filename: str, file name of config file
        '''
        if filename!=None:
            if os.path.exists(filename):
                self.config.read(filename)
                self._copyConfigtoSelf()
                self._updateSelf()
        return
    
    def updateConfig(self, filename=None, args=None, **kwargs):
        '''update config according to config file, args(from sys.argv) or **kwargs
        
        param filename: string, name of config file
        param args: list of str, usually be sys.argv
        param **kwargs: you can use like 'xbeamcenter=1024' or a dict to update the value of xbeamcenter
        '''
        if filename!=None:
            rv = self.parseConfigFile(filename)
        if args!=None:
            rv = self.parseArgs(args)
        if kwargs!={}:
            rv = self.parseKwargs(**kwargs)
        
        #read a configfile if specified in args or kwargs
        if (self.configfile!='')and(self.configfile!=None):
            self.parseConfigFile(filename=self.configfile)
            self.configfile = ''
            rv = self.parseArgs(args) if args!=None else None
            rv = self.parseKwargs(**kwargs) if kwargs!={} else None
            
        if (filename==None)and(args==None)and(kwargs=={}):
            rv = self._updateSelf()
        self._postProcessing(**kwargs)
        return rv
    
    def writeConfig(self, filename, mode='short'):
        '''write config to file
        
        param filename: string, name of file
        param mode: string, 'short' or 'full' ('s' or 'f'). 
            'short', options with 'config'=='f' will not be written into config file
            'full', all available options in self.config will be written to config file
        '''
        self._updateSelf()
        #func decide if wirte the option to config according to mode
        #options not present in self._optdata will not be written to config
        if mode.startswith('s'):
            mcond = lambda optname: self._optdata.get(optname, {'config':'n'}).get('config', 'a')=='a'
        else:
            mcond = lambda optname: self._optdata.get(optname, {'config':'n'}).get('config', 'a')!='n'
        
        lines = []        
        for section in self.config._sections:
            tlines = []
            for (key, value) in self.config._sections[section].items():
                if (key != "__name__") and mcond(key):
                    tlines.append("%s = %s" %(key, str(value).replace('\n', '\n\t')))
            if len(tlines)>0:
                lines.append("[%s]" % section)
                lines.extend(tlines)
                lines.append('')
        rv = "\n".join(lines) + "\n"
        fp = open(filename, 'w')
        fp.write(rv)
        fp.close()
        return
    
    def getHeader(self, title=None, mode='short'):
        '''get a header of configurations values, by default, all options are written to a header, 
        It can be disabled by set 'header' to False in self._optdata
        
        param title: str, title of header
        param mode: string, 'short' or 'full' ('s' or 'f'). 
            'short', options with 'config'=='f' will not be written into config file
            'full', all available options in self.config will be written to config file
        
        return: string, lines that can be directly writen to a text file
        '''
        
        lines = []
        title = '#Configuration information#' if title == None else '# %s #' % title
        lines.append(title)
        #func decide if wirte the option to header according to mode
        #options not present in self._optdata will not be written to header
        if mode.startswith('s'):
            mcond = lambda optname: self._optdata.get(optname, {'config':'n'}).get('config', 'a')=='a'
        else:
            mcond = lambda optname: self._optdata.get(optname, {'config':'n'}).get('config', 'a')!='n'
                
        lines = []
        for secname in self._configlist.keys():
            tlines = []
            for optname in self._configlist[secname]:
                if mcond(optname):
                    value = getattr(self, optname)
                    strvalue = ', '.join(map(str, value)) if type(value)==list else str(value)
                    tlines.append("%s = %s" % (optname, strvalue))
            if len(tlines)>0:
                lines.append("[%s]" % secname)
                lines.extend(tlines)
                lines.append('')
        lines.append('# data #')         
        rv = "\n".join(lines) + "\n"
        return rv
    
    def _postProcessing(self, **kwargs):
        '''post processing after parse args or kwargs or config file
        this method is called in self.updateConfig, after all file/args/kwargs
        some options are reset to their default to pervent process twice
        
        param kwargs: keywords=value
        '''
        
        self._additionalPostProcessing(**kwargs)
        
        if (self.createconfig!='')and(self.createconfig!=None):
            self.writeConfig(self.createconfig, 'short')
            self.createconfig = ''
        if (self.createconfigfull!='')and(self.createconfigfull!=None):
            self.writeConfig(self.createconfigfull, 'full')
            self.createconfigfull = ''
        return
    
#VERY IMPORTANT!!!
def initConfigClass(config):
    '''init config class and add options to class
    this funtion should be executed to generate options according to metadata defined in class
    '''
    config.config = ConfigParser.ConfigParser(dict_type = OrderedDict)
    config.args = argparse.ArgumentParser(description='Configuration')
    config._configlist = OrderedDict({})
    
    config._optdatalist = config._optdatalistpre + config._optdatalist + config._optdatalistext
    config._optdata = dict(config._optdatalist)
    config._detectAddSectionsC()
    for opt in config._optdatalist:
        key = opt[0]
        config._addOptC(key)
    return

#VERY IMPORTANT!!!
#add options to class
#initConfigClass(ConfigBase)
