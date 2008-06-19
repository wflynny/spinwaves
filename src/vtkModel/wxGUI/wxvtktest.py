#!/usr/bin/env python

"""This is a test of wxVTKRenderWindow"""

import wx
from picker import Picker
from wxVTKRenderWindowInteractor import *
from vtkModel.SpaceGroups import *
from vtkModel.VTKDrawer import *
from vtkModel.MagneticCellClass import *
from vtkModel.CellClass import *
import random

class Frame(wx.Frame):
    def __init__(self, parent, id):
        wx.Frame.__init__(self, parent, id, 'Magnetic Cell', size= (900,500))
        
        self.initVTKWindow()
        
        
        #Add Menus
        menuBar = wx.MenuBar()
        
        #Add File Menu
        fileMenu = wx.Menu()
        fileMenu.Append(wx.NewId(), "&Open")
        fileMenu.Append(wx.NewId(), "&Save")
        fileMenu.Append(wx.NewId(), "&Quit")
        menuBar.Append(fileMenu, "&File")
        
        self.SetMenuBar(menuBar)

    
    
    def initVTKWindow(self):
        #Code from wxRenderWindowInterActor Sample
       
        self.window = wxVTKRenderWindowInteractor(self, -1, self.afterRender)
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.window, 1, wx.EXPAND)
        self.SetSizer(sizer)
        self.Layout()
        # It would be more correct (API-wise) to call widget.Initialize() and
        # widget.Start() here, but Initialize() calls RenderWindow.Render().
        # That Render() call will get through before we can setup the
        # RenderWindow() to render via the wxWidgets-created context; this
        # causes flashing on some platforms and downright breaks things on
        # other platforms.  Instead, we call widget.Enable().  This means
        # that the RWI::Initialized ivar is not set, but in THIS SPECIFIC CASE,
        # that doesn't matter.
 #       self.window.Enable(1)
    
        self.window.AddObserver("ExitEvent", lambda o,e,f=self: f.Close())
    
        #My Code
        Space_Group = sg50
        unitcell = Cell(Space_Group)
        atomPos = [.25, .25, .5]
    
        #Create the unit cell
        randGen = random.Random()
        unitcell.generateAtoms(atomPos, "atom1" , .05, randGen.uniform(0,1), randGen.uniform(0,1), randGen.uniform(0,1))
        
        #Create the Magnetic Cell
        self.MagCell = MagneticCell(unitcell, 2,2,3, Space_Group)
        AllAtoms = self.MagCell.getAllAtoms()
        for i in range(0, len(AllAtoms)):
            print i, AllAtoms[i]
        self.MagCell.addBond(AllAtoms[0], AllAtoms[1])
         
        # a renderer for the data
        ren1 = vtkRenderer()
        ren1.SetBackground(1,1,1)
        
        #Add the renderer to the window
        self.window.GetRenderWindow().AddRenderer(ren1)
            
        #Create vtkDrawer
        self.drawer = vtkDrawer(ren1)
        
        #Add my picker
        Picker(self.drawer, self.window._Iren, ren1)
        
        #Draw the Magnetic Cell
        self.drawer.drawMagneticCell(self.MagCell)
 
        
    def afterRender(self):
        """Does everythinng that must be done after hte first render,
        such as creating the axes and atom labels"""
        
        self.drawer.addAxes()
        self.drawer.labelAtoms(self.MagCell)



class App(wx.App):
    
    def OnInit(self, redirect = False, filename = None):
        wx.App.__init__(self, redirect, filename)
        self.frame = Frame(None, -1)
        self.frame.Show()
        self.SetTopWindow(self.frame)
        return True
    
if __name__ == '__main__':

    app = App(False)
    app.MainLoop()