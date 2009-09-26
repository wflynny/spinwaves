"""This is a setup script to install spinwaves with setup tools.
Author: Tom 9/25/09"""
from os.path import join

#Bootstrap setuptools for systems that do not have it installed
from ez_setup import use_setuptools
use_setuptools()

#For compiling C extensions
#from distutils.core import Extension

from setuptools import setup, find_packages, Extension
setup(
	name = "Spinwaves",
	version = "0.1",
	url = 'http://spinwaves.googlecode.com',
	packages = find_packages(exclude = ['sympy_WORKING']),

	#Add C Extensions
	ext_modules=[Extension(join('MonteCarlo','_monteCarlo'), [join('lib',f) for f in ['main1.c','dSFMT.c']])],
	

	entry_points = {'gui_scripts':['spinwaves = spinwaves.vtkModel.wxGUI.GUI_Main:main']},

	install_requires = ['numpy', 'wxPython', 'sympy', 'matplotlib'],

	zip_safe = False
)

#Extension('spinwaves.MonteCarlo._monteCarlo.so', ['main1.c', 'dSFMT.c'], include_dirs=['lib'])
