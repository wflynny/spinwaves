#!/usr/bin/env python

import sys
from multiprocessing import freeze_support
import spinwaves
print spinwaves.__path__
sys.path.append(spinwaves.__path__[0])
import vtkModel

from spinwaves.vtkModel.wxGUI.GUI_Main import main

import warnings
warnings.simplefilter('ignore', UserWarning)

if __name__ == "__main__":
    freeze_support()
    main()
