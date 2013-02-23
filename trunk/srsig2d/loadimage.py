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
import fabio, fabio.openimage
import time
import os
import fnmatch
import threading
import sys
import zlib
import hashlib
import scipy.io

class LoadImage(threading.Thread):
    def __init__(self, p):
        threading.Thread.__init__(self)
        self.config = p
        self.configlist = ['xdimension',
                           'ydimension',
                           'tifdirectory',
                           'includepattern',
                           'excludepattern',
                           'fliphorizontal',
                           'flipvertical',
                           'backgroundfile',
                           #'backgroundenable',
                           'method',
                           'integrationspace',
                           'savedirectory',
                           'filenameplus']
        for optionname in self.configlist:
            setattr(self.__class__, optionname, self.configProperty(optionname))
        self.prepareCalculation()
        return
    
    def configProperty(self, nm):
        '''helper function of options delegation'''
        rv = property(fget = lambda self: getattr(self.config, nm))
        return rv
    
    def prepareCalculation(self):
        if (self.backgroundfile != '') and (os.path.exists(self.backgroundfile)):
            temp = fabio.openimage.openimage(self.backgroundfile)
            self.backgroundpic = self.cropAndFlipImage(temp.data)
            self.backgroundenable = True
        else:
            self.backgroundenable = False
        return
    
        
    def flipImage(self, pic):
        '''flip image if configured in config 
        '''
        if self.fliphorizontal:
            pic = pic[:,::-1]
        if self.flipvertical:
            pic = pic[::-1,:]
        return pic
    
    def cropAndFlipImage(self, pic):
        '''flip image first then crop it.
        '''
        pic = self.flipImage(pic)
        return pic
    
    def loadImage(self, filename):
        '''load image, then deduct the background.
        filename:       str, image file name
        
        return:         2d ndarray, 2d image data
        '''
        if os.path.exists(filename):
            filenamefull = filename
        else:
            filenamefull = os.path.normpath(self.tifdirectory+'/'+filename)
        image = fabio.openimage.openimage(filenamefull)
        image = self.cropAndFlipImage(image.data)
        image[image<0] = 0
        if self.backgroundenable:
            image = image - self.backgroundpic
        return image
    
    def genFileList(self, opendir=None, includepattern=None, excludepattern=None):
        '''generate the list of file in opendir according to include/exclude pattern
        opendir:        string, the directory of files
        includepattern: string, wildcard of files that will be loaded into PDFLive
        excludepattern: list of string, a list of wildcard of files that will be blocked
        
        return:         list of string, a list of filenames
        '''
        opendir = self.tifdirectory if opendir == None else opendir
        includepattern = self.includepattern if includepattern == None else includepattern
        excludepattern = self.excludepattern if excludepattern == None else excludepattern
        fileset = self.genFileSet(opendir, includepattern, excludepattern)
        return sorted(list(fileset))
    
    def genFileSet(self, opendir=None, includepattern=None, excludepattern=None):
        '''generate the set of file in opendir according to include/exclude pattern
        opendir:        string, the directory of files
        includepattern: string, wildcard of files that will be loaded into PDFLive
        excludepattern: list of string, a list of wildcard of files that will be blocked
        
        return:         set of string, a set of filenames
        '''
        opendir = self.tifdirectory if opendir == None else opendir
        includepattern = self.includepattern if includepattern == None else includepattern
        excludepattern = self.excludepattern if excludepattern == None else excludepattern
        filelist = fnmatch.filter(os.listdir(opendir), includepattern)
        fileset = set(filelist)
        if len(excludepattern) > 0:
            for excludep in excludepattern:
                excludefilelist = fnmatch.filter(filelist, excludep)
                fileset = fileset - set(excludefilelist)
        return fileset
    def checkCRC32(self, filename):
        '''calculate the crc32 value of file'''
        try:
            fd = open(filename, 'rb')
        except:
            return 'Read error'
        eachLine = fd.readline()
        prev = 0
        while eachLine:
            prev = zlib.crc32(eachLine, prev)
            eachLine = fd.readline()
        fd.close()
        return prev
    
    def checkMD5(self, filename, blocksize=65536):
        '''calculate the MD5 value of file'''
        try:
            fd = open(filename, 'rb')
        except:
            return 'Read error'
        buf = fd.read(blocksize)
        md5 = hashlib.md5()
        while len(buf) > 0:
            md5.update(buf)
            buf = fd.read(blocksize)
        fd.close()
        return md5.hexdigest()
