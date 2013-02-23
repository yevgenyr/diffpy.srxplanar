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
'''create a configuration template or create configuration file with interactively 
'''

class configtemplate(object):
    def __init__(self):
        self.conftemplate='''
[BaseSection]
#########################################
fit2dconfig = 
tifdirectory = ./
savedirectory = ./re
backgroundfile = 
maskfit2d =
#########################################
integrationspace = qspace
method = normal
xrduncertaintyenable = True
#########################################
xdimension = 2048
ydimension = 2048
xpixelsize = 0.2
ypixelsize = 0.2
#########################################
wavelength = 0.1078
xbeamcenter = 1298.930
ybeamcenter = 1010.190
distance = 369.579
rotationd = 49.715
tiltd = -0.480
#########################################
xrdtthstepd = 0.03
xrdtthmaxd = 50.0
xrdqmax = 50.0
xrdqstep = 0.03
xrdtthvarstepd = 0.01
xrdqvarstep = 0.01
#########################################

includepattern = *.tif
excludepattern = *.dark.tif
fliphorizontal = False
flipvertical = True

xrdazimuthstepd = 0.5
xrdspotty = False
sacorrectionenable = True
polcorrectionenable = False
polcorrectf = 0.95
maskenable = True
maskselfcorrenable = True
maskselfcorrmethod = Fast
maskselfcorrcounthi = 0.5
maskselfcorrcountlow = 0.5
maskselfcorrpercentagehi = 0.95
maskselfcorrpercentagelow = 0.05
maskselfcorraccuenable = False
regulartmatrixenable = False

maskfit2denable = False
masktiffenable = False
masktiff = 
maskcakeenable = False
maskcake = 
maskboxenable = False
maskbox = 

gsasoutput = None
filenameplus = 
'''
        self.optionlist=[]
        return
    
    def addoption(self, name, defaultvalue, datatype, description):
        self.optionlist.append([name, defaultvalue, datatype, description])
        return
    
    def createnewtemplate(self, filename):
        fp = open(filename, 'w')
        fp.write(self.conftemplate)
        fp.close()
        return
    
temp = configtemplate()
temp.createnewtemplate('newconf.cfg')
temp.createnewtemplatei('newconfi.cfg')
