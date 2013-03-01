#!/usr/bin/env python
'''class for calibration geometry parameters, will add calibration functino in the future.
Now it only read the calibration results obtained from Fit2D.
'''
import re
from diffpy.srxplanar.srxplanarconfig import _configPropertyRW

class Calibration(object):
    
    # define configuration properties that are forwarded to self.config
    fit2dconfig = _configPropertyRW('fit2dconfig')
    xbeamcenter = _configPropertyRW('xbeamcenter')
    ybeamcenter = _configPropertyRW('ybeamcenter')
    rotationd = _configPropertyRW('rotationd')
    tiltd = _configPropertyRW('tiltd')
    distance = _configPropertyRW('distance')
    wavelength = _configPropertyRW('wavelength')
    
    def __init__(self, p):
        self.config = p
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
                self.config.updateConfig()
        return