#!/usr/bin/env python

"""This is a test of wxVTKRenderWindow"""

import wx
import wx.grid
from picker import Picker
from wxVTKRenderWindowInteractor import *
from vtkModel.SpaceGroups import *
from vtkModel.VTKDrawer import *
from vtkModel.MagneticCellClass import *
from vtkModel.CellClass import *
#import random


class atomPanel(wx.Panel):
    def __init__(self, parent, id):
        wx.Panel.__init__(self, parent, id, size = (200,200))
        self.SetSizer()

class atomListGrid(wx.grid.Grid):
    def __init__(self, parent, id):
        wx.grid.Grid.__init__(self, parent, id)
        self.CreateGrid(3,5)
        self.SetColLabelValue(0, "Name")
        self.SetColLabelValue(1, "Atomic Number")
        self.SetColLabelValue(2, "x")
        self.SetColLabelValue(3, "y")
        self.SetColLabelValue(4, "z")
        self.AutoSize()

class vtkPanel(wx.Panel):
    def __init__(self, parent, id):
        wx.Panel.__init__(self, parent)
        self.initVTKWindow()

    
    def initVTKWindow(self):
        #Code from wxRenderWindowInterActor Sample
       
        self.window = wxVTKRenderWindowInteractor(self, -1)
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
        self.window.Render()
        
#        self.initializeVTKData()
#        self.draw()
        self.window.Bind(wx.EVT_KEY_DOWN, self.OnKeyEvent)
 
    
    
    def OnKeyEvent(self, event):
        event.Skip()
        
        #Delete if key pressed is Del key  Del = 127
        if event.GetKeyCode() == 127:
            self.OnDelete(event)
        
    
    
    def OnDelete(self,event):
        event.Skip()
        selectedObj = self.drawer.getObjFromActor(self.picker.getPicked())
        if isinstance(selectedObj, Bond):
            self.MagCell.deleteBond(selectedObj)
            #handle atoms next
        
        self.draw()
    
#    def initializeVTKData(self):
                #My Code
#        Space_Group = sg141
#        unitcell = Cell(Space_Group)
#        atomPos = [0, 0, 0]
    
        #Create the unit cell
#        unitcell.generateAtoms(atomPos, "atom1")

        
        #Create the Magnetic Cell
#        self.MagCell = MagneticCell(unitcell, 1,1,1, Space_Group)
        
#        AllAtoms = self.MagCell.getAllAtoms()
#        for i in range(0, len(AllAtoms)):
#            print i, AllAtoms[i]
#        self.MagCell.addBond(AllAtoms[0], AllAtoms[1])
        
        
        
    def draw(self):
         
        # a renderer for the data
        ren1 = vtkRenderer()
        ren1.SetBackground(1,1,1)
        
        #Add the renderer to the window
        self.window.GetRenderWindow().AddRenderer(ren1)
            
        #Create vtkDrawer
        self.drawer = vtkDrawer(ren1)
        
        #Add my picker
        self.picker = Picker(self.drawer, self.window._Iren, ren1)

        #Draw the Magnetic Cell
        self.drawer.drawMagneticCell(self.MagCell)
        
        self.window.setUpRender()    

        self.drawer.addAxes()
        self.drawer.labelAtoms(self.MagCell)
        self.window.Render()
        
    
    def openCif(self, filename):
        self.MagCell = magneticCellFromCif(filename)
        self.draw()
    
    
    def getStatusText(self):
        return self.picker.getPicked()




class Frame(wx.Frame):
    def __init__(self, parent, id):
        wx.Frame.__init__(self, parent, id, 'Magnetic Cell', size= (900,600))

        self.vtkPanel = vtkPanel(self, -1)
        
        #Add Menus
        self.AddMenus()
        
        #Add Tool Bar
#        self.AddToolBar()
        
        #Add Status Bar                     
#        self.AddStatusBar()

   
   
#    def AddToolBar(self):
#        toolbar = self.CreateToolBar()
#        toolbar.AddSimpleTool()
        
#    def AddStatusBar(self):
#        self.statusBar = self.CreateStatusBar()
#        self.Bind(wx.EVT_LEFT_UP, self.OnClick)
#        
#    
#    def OnClick(self, event):
#        self.statusBar.SetStatusText("This")
#        #self.statusBar.SetStatusText(self.vtkPanel.getStatusText())
#        event.Skip()
        
    def AddMenus(self):
                #Add Menus
        menuBar = wx.MenuBar()
        
        #Add File Menu
        fileMenu = wx.Menu()
        newMenuItem = fileMenu.Append(wx.NewId(), "&New Magnetic Cell")
        openMenuItem = fileMenu.Append(wx.NewId(), "&Open")
        saveMenuItem = fileMenu.Append(wx.NewId(), "&Save")
        quitMenuItem = fileMenu.Append(wx.NewId(), "&Quit")
        menuBar.Append(fileMenu, "&File")
        
        #Add Model Menu
        modelMenu = wx.Menu()
        addCellMenuItem = modelMenu.Append(wx.NewId(), "Add Atom")
        addBondMenuItem = modelMenu.Append(wx.NewId(), "Add Bond")
        deleteMenuItem = modelMenu.Append(wx.NewId(), "Delete")
        menuBar.Append(modelMenu, "Model")
        self.SetMenuBar(menuBar)
        
        #Bind Events
        self.Bind(wx.EVT_MENU, self.OnCloseMe, quitMenuItem)
        self.Bind(wx.EVT_MENU, self.OnSave, saveMenuItem)
        self.Bind(wx.EVT_MENU, self.OnOpenFile, openMenuItem)
        self.Bind(wx.EVT_MENU, self.vtkPanel.OnDelete, deleteMenuItem)
#        self.Bind(wx.EVT_MENU, self.OnNew, newMenuItem)
    
 #   def OnNew(self, event):
        #Add drawing panel
#        self.vtkPanel.draw()
#        self.GetEventHandler().ProcessEvent(wx.SizeEvent())
    
    def OnCloseMe(self, event):
        self.Close(True)
    
    def OnSave(self, event):
        saveDialog = wx.FileDialog(self, "Save File", style = wx.SAVE)
        if saveDialog.ShowModal() == wx.ID_OK:
            print saveDialog.GetPath()
        saveDialog.Destroy()
        
    def OnOpenFile(self, event):
        saveDialog = wx.FileDialog(self, "Open File", style = wx.OPEN, wildcard = "*.cif")
        if saveDialog.ShowModal() == wx.ID_OK:
            self.vtkPanel.openCif(saveDialog.GetPath())
        saveDialog.Destroy()
        



class App(wx.App):
    def __init__(self, redirect = False, filename = None):
        wx.App.__init__(self, redirect, filename)
    
    def OnInit(self):
        self.frame = Frame(None, -1)
        self.frame.Show()
        self.SetTopWindow(self.frame)
        frame1 = wx.Frame(self.frame, -1)
        atomListGrid(frame1, -1)
        frame1.Show()
        return True
    


if __name__ == '__main__':
    app = App()
    app.MainLoop()