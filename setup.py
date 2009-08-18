#! /usr/bin/env python
# Last Change: Wed Nov 05 11:00 AM 2008 J

# Copyright (C) 2008 Cournapeau David <cournape@gmail.com>

descr   = """\
The spinwaves module does spinwave stuff.
"""

import os
import sys

DESCRIPTION         = 'Spinwave Calculator'
LONG_DESCRIPTION    = descr
MAINTAINER          = ''
MAINTAINER_EMAIL    = ''
LICENSE             = 'BSD'
#URL                 = ''
#DOWNLOAD_URL        = URL
VERSION             = '0.2'

classifiers = [
              'Development Status :: '+VERSION,
              'Environment :: Console',
              'Intended Audience :: Developers',
              'Intended Audience :: Science/Research',
              'License :: OSI Approved :: BSD License',
              'Topic :: Scientific/Engineering',
               'Operating System :: Microsoft :: Windows',
               'Operating System :: POSIX',
               'Operating System :: Unix',
               'Operating System :: MacOS',
              ]

import setuptools
from numpy.distutils.core import setup
from os.path import join

def configuration(parent_package='', top_path=None, package_name='spinwaves'):
    if os.path.exists('MANIFEST'): os.remove('MANIFEST')

    from numpy.distutils.misc_util import Configuration
    config = Configuration(package_name,
                           parent_package,
                           top_path,
                           #namespace_packages=['scikits'],
                           version = VERSION,
                           maintainer  = MAINTAINER,
                           maintainer_email = MAINTAINER_EMAIL,
                           description = DESCRIPTION,
                           license = LICENSE,
                           #url = URL,
                           #download_url = DOWNLOAD_URL,
                           long_description = LONG_DESCRIPTION)

    sources = [join('lib',f) for f in ['main1.c','dSFMT.c']]
    config.add_extension(join('MonteCarlo','_monteCarlo'),
                         sources=sources,
                         define_macros=[('DSFMT_MEXP',19937)])

    scripts = [join('bin',f) for f in ['spinwaves_ui.py']]
    config.add_scripts(scripts)
    return config

def main():
    setup(
        install_requires = ['numpy'],
        #namespace_packages = ['scikits'],
        packages = setuptools.find_packages(),
        include_package_data = True,
        #test_suite="nose.collector",
        zip_safe = False, # the package can run out of an .egg file
        classifiers = classifiers,
        platforms = ["Windows", "Linux", "Solaris", "Mac OS-X", "Unix"],
        configuration = configuration,
    )

if __name__ == "__main__":
    main()

