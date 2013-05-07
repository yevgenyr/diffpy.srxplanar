#!/usr/bin/env python

# Installation script for diffpy.srxplanar

from setuptools import setup, find_packages

# define distribution
dist = setup(
        name = 'diffpy.srxplanar',
        version = '0.2',
        namespace_packages = ['diffpy'],
        packages = find_packages(),
        entry_points = {
            # define console_scripts here, see setuptools docs for details.
            'console_scripts' : ['srxplanar = diffpy.srxplanar.srxplanar:main'
                                 ],
        },
        author = 'Simon J.L. Billinge',
        author_email = 'sb2896@columbia.edu',
        description = '2D diffraction image integration and uncertainty propagation',
        license = 'BSD',
        url = '',
        keywords = '',
)

# End of file
