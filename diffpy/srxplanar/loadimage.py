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
import os
import fnmatch
import sys

class LoadImage(object):
    def __init__(self, p):
        self.config = p
        self.configlist = ['xdimension',
                           'ydimension',
                           'tifdirectory',
                           'includepattern',
                           'excludepattern',
                           'fliphorizontal',
                           'flipvertical',
                           'backgroundfile',
                           ]
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
        image = self.flipImage(image.data)
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
