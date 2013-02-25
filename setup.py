#!/usr/bin/env python

# Installation script for diffpy.pdfgetx

from setuptools import setup, find_packages

# define distribution
dist = setup(
        name = 'SrXPlanar',
        version = '0.1',
        namespace_packages = ['diffpy'],
        packages = find_packages(),
        #entry_points = {
        #    # define console_scripts here, see setuptools docs for details.
        #    'console_scripts' : [],
        #},
        author = 'Simon J.L. Billinge',
        author_email = 'sb2896@columbia.edu',
        description = '2D diffraction image integration and error propagation',
        license = 'BSD',
        url = '',
        keywords = '',
)

# End of file
