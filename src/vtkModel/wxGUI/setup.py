from distutils.core import setup
import glob
import py2exe

#This is used by py2exe to create a windows executable
#In DOS prompt:
#python setup.py py2exe

setup(windows=["wxvtktest.py"])
