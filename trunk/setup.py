#!/usr/bin/env python

# Installation script for diffpy.pdfgetx

from setuptools import setup, find_packages

# define distribution
dist = setup(
        name = 'SrSig2D',
        version = '0.1',
        packages = ['srsig2d'],
        entry_points = {
            # define console_scripts here, see setuptools docs for details.
            'console_scripts' : [
                'srsig2d = srsig2d.srsig2d:main',
            ],
        },
        author = 'Simon J.L. Billinge',
        author_email = 'sb2896@columbia.edu',
        description = '2D diffraction image integration and error propagation',
        license = 'BSD',
        url = '',
        keywords = '',
)

# End of file
