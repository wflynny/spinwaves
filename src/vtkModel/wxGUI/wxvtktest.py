#!/usr/bin/env python

"""This is a test of wxVTKRenderWindow"""

import wx
from picker import Picker
from wxVTKRenderWindowInteractor import *
from vtkModel.AtomGeneratorSample import *

class Frame(wx.Frame):
    def __init__(self, parent, id):
        wx.Frame.__init__(self, parent, id, 'Magnetic Cell', size= (900,500))
        
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
    
        
        #My Code
        Space_Group = sg65
        unitcell = Cell(Space_Group)
        atomPos = [.25, .25, .5]
    
        #Create the unit cell
        randGen = random.Random()
        unitcell.generateAtoms(atomPos, "atom1" , .05, randGen.uniform(0,1), randGen.uniform(0,1), randGen.uniform(0,1))
        
        #Create the Magnetic Cell
        MagCell = MagneticCell(unitcell, 1,2,3, Space_Group)
        AllAtoms = MagCell.getAllAtoms()
        MagCell.addInterCellularBond(AllAtoms[0], AllAtoms[6])
         
        # a renderer for the data
        ren1 = vtkRenderer()
        ren1.SetBackground(1,1,1)
    
        # a render window to display the contents
        renWin = self.window.GetRenderWindow()
        renWin.AddRenderer(ren1)
        
        #Add my picker
        Picker(MagCell, self.window._Iren, ren1)
            
        #Draw the Magnetic Cell
        MagCell.drawCell(ren1)
        
        #Add Axes
        axes = vtkAxes()
        axes.SetOrigin(0,0,0)
        axesMapper = vtkPolyDataMapper()
        axesMapper.SetInputConnection(axes.GetOutputPort())
        axesActor = vtkActor()
        axesActor.SetMapper(axesMapper)
        ren1.AddActor(axesActor)
        xLabel = vtkVectorText()
        yLabel = vtkVectorText()
        zLabel = vtkVectorText()
        xLabel.SetText("x")
        yLabel.SetText("y")
        zLabel.SetText("z")
        xLabelMapper = vtkPolyDataMapper()
        yLabelMapper = vtkPolyDataMapper()
        zLabelMapper = vtkPolyDataMapper()
        xLabelMapper.SetInputConnection(xLabel.GetOutputPort())
        yLabelMapper.SetInputConnection(yLabel.GetOutputPort())
        zLabelMapper.SetInputConnection(zLabel.GetOutputPort())
        xLabelActor = vtkFollower()
        yLabelActor = vtkFollower()
        zLabelActor = vtkFollower()
        xLabelActor.SetMapper(xLabelMapper)
        yLabelActor.SetMapper(yLabelMapper)
        zLabelActor.SetMapper(zLabelMapper)
        xLabelActor.SetScale(0.1,0.1,0.1)
        yLabelActor.SetScale(0.1,0.1,0.1)
        zLabelActor.SetScale(0.1,0.1,0.1)
        xLabelActor.AddPosition(1,0,0)
        yLabelActor.AddPosition(0,1,0)
        zLabelActor.AddPosition(0,0,1)
        xLabelActor.GetProperty().SetColor(0,0,0)
        yLabelActor.GetProperty().SetColor(0,0,0)
        zLabelActor.GetProperty().SetColor(0,0,0)
        ren1.AddActor(xLabelActor)
        ren1.AddActor(yLabelActor)
        ren1.AddActor(zLabelActor)
        

        
 #       xLabelActor.SetCamera(ren1.GetActiveCamera())
#        yLabelActor.SetCamera(ren1.GetActiveCamera())
#        zLabelActor.SetCamera(ren1.GetActiveCamera())

#        self.Bind(wx.EVT_BUTTON, self.ButtonOnVTK, self.window)
#        self.Bind(wx.EVT_CLOSE, self.OnCloseWindow)
    
 #   def ButtonOnVTK(self, event):
 #       iren.InvokeEvent("LeftButtonPressEvent")
 #       print "Button Pressed"
        
 #   def OnCloseWindow(self, event):
 #       self.Destroy()

class App(wx.App):
    
    def OnInit(self, redirect = False, filename = None):
        wx.App.__init__(self, redirect, filename)
        self.frame = Frame(None, -1)
        self.frame.Show()
        self.SetTopWindow(self.frame)
        return True
    
if __name__ == '__main__':

    app = App()
    app.MainLoop()