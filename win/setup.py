from distutils.core import setup
import glob
import py2exe
import matplotlib


#This is used by py2exe to create a windows executable
#In DOS prompt:
#python setup.py py2exe



python_dir = "C:\Python25"
#Matplotlib code taken from: http://www.py2exe.org/index.cgi/MatPlotLib

# We need to exclude matplotlib backends not being used by this executable.  You may find
# that you need different excludes to create a working executable with your chosen backend.
# We also need to include include various numerix libraries that the other functions call.

opts = {
    'py2exe': { "includes" : ["sip", "matplotlib.backends",  "matplotlib.backends.backend_qt4agg",                                "matplotlib.figure","pylab", "numpy", "matplotlib.numerix.fft",
                               "matplotlib.numerix.linear_algebra", "matplotlib.numerix.random_array",
                               "matplotlib.backends.backend_tkagg"],
                'excludes': ['_gtkagg', '_tkagg', '_agg2', '_cairo', '_cocoaagg',
                             '_fltkagg', '_gtk', '_gtkcairo', ],
                'dll_excludes': ['libgdk-win32-2.0-0.dll',
                                 'libgobject-2.0-0.dll']
              }
       }
 
# Save matplotlib-data to mpl-data ( It is located in the matplotlib\mpl-data
# folder and the compiled programs will look for it in \mpl-data
# note: using matplotlib.get_mpldata_info
data_files = [(r'mpl-data', glob.glob(python_dir + r'\Lib\site-packages\matplotlib\mpl-data\*.*')),
                    # Because matplotlibrc does not have an extension, glob does not find it (at least I think that's why)
                    # So add it manually here:
                  (r'mpl-data', [python_dir + r'\Lib\site-packages\matplotlib\mpl-data\matplotlibrc']),
                  (r'mpl-data\images',glob.glob(python_dir + r'\Lib\site-packages\matplotlib\mpl-data\images\*.*')),
                  (r'mpl-data\fonts',glob.glob(python_dir + r'\Lib\site-packages\matplotlib\mpl-data\fonts\*.*'))]
#for inno setup
data_files.append("screen.ico")
data_files.append("..\spinwaves\MonteCarlo\_monteCarlo.pyd")

setup(windows=[{"script" : "Spinwaves.py", "icon_resources": [(0x0004, "screen.ico")]}], 
      options=opts,   data_files=data_files)
