#!/usr/bin/env python

import sys

import spinwaves
print spinwaves.__path__
sys.path.append(spinwaves.__path__[0])
import vtkModel
from spinwaves.vtkModel.wxGUI.GUI_Main import main

if __name__ == "__main__":
    main()
