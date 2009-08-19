#!/usr/bin/env python

import sys
import os
#Add the main folder (up one level) to the path
spinwaves_path = os.path.split(os.path.dirname(os.path.abspath(__file__)))[0]
sys.path.append(spinwaves_path)
import spinwaves
#print spinwaves.__path__
#sys.path.append(spinwaves.__path__[0])
from spinwaves.vtkModel.wxGUI.GUI_Main import main

if __name__ == "__main__":
    main()
