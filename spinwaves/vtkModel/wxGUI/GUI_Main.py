#!/usr/bin/env python

import time
import gc
from multiprocessing import Process, Pipe
import copy
import wx
import wx.grid
from wx.py.dispatcher import connect, send
from scipy.optimize.optimize import fmin, fmin_cg
from sympy import pi, cos
import numpy

from picker import Picker
from wxVTKRenderWindowInteractor import *
from spinwaves.vtkModel.SpaceGroups import *
from spinwaves.vtkModel.VTKDrawer import *
from spinwaves.vtkModel.MagneticCellClass import *
from spinwaves.vtkModel.CellClass import *
#import random
import spinwaves.vtkModel.SpaceGroups
from Session import Session
from spinwaves.cross_section.general_case2 import run_cross_section, run_eval_cross_section
import spinwaves.MonteCarlo.CSim as CSim
import spinwaves.spinwavecalc.spinwavepanel as spinwavepanel
import spinwaves.spinwavecalc.spinwave_calc_file as spinwave_calc_file

import spinwaves.cross_section.util.printing as printing

from spinwaves.vtkModel.BondClass import JParam
from spinwaves.vtkModel.Parameter_Manager import Fitter
import spinwaves.cross_section.util.printing as printing
from spinwaves.utilities.fitting import ShowFittingFrame
from spinwaves.utilities.Processes import ProcessManager
from spinwaves.vtkModel.BondClass import bondPanel
#gc.enable()
#Atom and cell info window

class atomPanel(wx.Panel):
    """This is the panel in which the user enters information about the space group,
    dimensions of the unit cell and the atom information."""
    def __init__(self, parent, id, session):
        wx.Panel.__init__(self, parent, id)
        
        self.session = session
        
        #Add Space Group
        spaceGroupLabel = wx.StaticText(self, -1, "Space Group:")
        self.spaceGroupSpinner = wx.SpinCtrl(self, -1, "")
        self.spaceGroupSpinner.SetRange(1,230)
        self.spaceGroupSpinner.SetValue(1)
        
        #Add Atom List
        self.atomList = atomListGrid(self, -1, session = self.session)
         
        #Add a button on upper right to generate new image
        self.genButton = wx.Button(self, -1, "Generate")
        
        #Add a, b, c, Alpha, gamma, beta, Na, Nb, Nc
        aLabel = wx.StaticText(self, -1, "a:")
        self.aText = wx.TextCtrl(self, -1, size = (60, -1), style = wx.TE_RICH2)
        aSizer = wx.BoxSizer(wx.VERTICAL)
        aSizer.Add(aLabel, 0)
        aSizer.Add(self.aText, 0)
        
        bLabel = wx.StaticText(self, -1, "b:")
        self.bText = wx.TextCtrl(self, -1, size = (60, -1), style = wx.TE_RICH2)
        bSizer = wx.BoxSizer(wx.VERTICAL)
        bSizer.Add(bLabel, 0)
        bSizer.Add(self.bText, 0)
        
        cLabel = wx.StaticText(self, -1, "c:")
        self.cText = wx.TextCtrl(self, -1, size = (60, -1), style = wx.TE_RICH2)
        cSizer = wx.BoxSizer(wx.VERTICAL)
        cSizer.Add(cLabel, 0)
        cSizer.Add(self.cText, 0)
        
        alphaLabel = wx.StaticText(self, -1, "alpha:")
        self.alphaText = wx.TextCtrl(self, -1, size = (60, -1), style = wx.TE_RICH2)
        alphaSizer = wx.BoxSizer(wx.VERTICAL)
        alphaSizer.Add(alphaLabel, 0)
        alphaSizer.Add(self.alphaText, 0)
        
        betaLabel = wx.StaticText(self, -1, "beta:")
        self.betaText = wx.TextCtrl(self, -1, size = (60, -1), style = wx.TE_RICH2)
        betaSizer = wx.BoxSizer(wx.VERTICAL)
        betaSizer.Add(betaLabel, 0)
        betaSizer.Add(self.betaText, 0)
        
        gammaLabel = wx.StaticText(self, -1, "gamma:")
        self.gammaText = wx.TextCtrl(self, -1, size = (60, -1), style = wx.TE_RICH2)
        gammaSizer = wx.BoxSizer(wx.VERTICAL)
        gammaSizer.Add(gammaLabel, 0)
        gammaSizer.Add(self.gammaText, 0)
        
        #Magnetic Cell
        naLabel = wx.StaticText(self, -1, "Na:")
        self.naText = wx.TextCtrl(self, -1, size = (60, -1), style = wx.TE_RICH2)
        naSizer = wx.BoxSizer(wx.VERTICAL)
        naSizer.Add(naLabel, 0)
        naSizer.Add(self.naText, 0)
        
        nbLabel = wx.StaticText(self, -1, "Nb:")
        self.nbText = wx.TextCtrl(self, -1, size = (60, -1), style = wx.TE_RICH2)
        nbSizer = wx.BoxSizer(wx.VERTICAL)
        nbSizer.Add(nbLabel, 0)
        nbSizer.Add(self.nbText, 0)
        
        ncLabel = wx.StaticText(self, -1, "Nc:")
        self.ncText = wx.TextCtrl(self, -1, size = (60, -1), style = wx.TE_RICH2)
        ncSizer = wx.BoxSizer(wx.VERTICAL)
        ncSizer.Add(ncLabel, 0)
        ncSizer.Add(self.ncText, 0)
        
        #Cutoff
        cutoffNaLabel = wx.StaticText(self, -1, "Na:")
        self.cutoffNaText = wx.TextCtrl(self, -1, size = (60, -1), style = wx.TE_RICH2)
        cutoffNaSizer = wx.BoxSizer(wx.VERTICAL)
        cutoffNaSizer.Add(cutoffNaLabel, 0)
        cutoffNaSizer.Add(self.cutoffNaText, 0)
        
        cutoffNbLabel = wx.StaticText(self, -1, "Nb:")
        self.cutoffNbText = wx.TextCtrl(self, -1, size = (60, -1), style = wx.TE_RICH2)
        cutoffNbSizer = wx.BoxSizer(wx.VERTICAL)
        cutoffNbSizer.Add(cutoffNbLabel, 0)
        cutoffNbSizer.Add(self.cutoffNbText, 0)
        
        cutoffNcLabel = wx.StaticText(self, -1, "Nc:")
        self.cutoffNcText = wx.TextCtrl(self, -1, size = (60, -1), style = wx.TE_RICH2)
        cutoffNcSizer = wx.BoxSizer(wx.VERTICAL)
        cutoffNcSizer.Add(cutoffNcLabel, 0)
        cutoffNcSizer.Add(self.cutoffNcText, 0)
        
        
        #Organize a,b,c and alpha, gamma, beta , Na, Nb, Nc into a grid
        dimSizer = wx.GridSizer(cols = 3, hgap = 15, vgap = 5)
        dimSizer.Add(aSizer)
        dimSizer.Add(bSizer)
        dimSizer.Add(cSizer)
        dimSizer.Add(alphaSizer)
        dimSizer.Add(betaSizer)
        dimSizer.Add(gammaSizer)

        
        unitCellBox = wx.StaticBox(self, -1, "Unit Cell")
        unitCellSizer = wx.StaticBoxSizer(unitCellBox, wx.VERTICAL)
        unitCellSizer.Add(dimSizer)
        
        leftTopSizer = wx.GridBagSizer(2,2)
        leftTopSizer.Add(unitCellSizer, (0,0), (1,2))
        leftTopSizer.Add(wx.StaticText(self, -1, "Atoms:"), (1,0))
        self.atomSpinner = wx.SpinCtrl(self, -1, "")
        self.atomSpinner.SetRange(1,100)
        self.atomSpinner.SetValue(1)
        self.atomSpinner.Bind(wx.EVT_TEXT, self.OnGridResize, self.atomSpinner)
        leftTopSizer.Add(self.atomSpinner, (1,1))
        
        magCellBox = wx.StaticBox(self, -1, "Magnetic Cell")
        magCellSizer = wx.StaticBoxSizer(magCellBox, wx.HORIZONTAL)
        magCellSizer.Add(naSizer)
        magCellSizer.Add(nbSizer)
        magCellSizer.Add(ncSizer)
        
        cutoffBox = wx.StaticBox(self, -1, "Cutoff")
        cutoffSizer = wx.StaticBoxSizer(cutoffBox, wx.HORIZONTAL)
        cutoffSizer.Add(cutoffNaSizer)
        cutoffSizer.Add(cutoffNbSizer)
        cutoffSizer.Add(cutoffNcSizer)
        
        
        spaceGroupSizer = wx.BoxSizer(wx.HORIZONTAL)
        spaceGroupSizer.Add(spaceGroupLabel, 0)
        spaceGroupSizer.Add(self.spaceGroupSpinner, 0)
        
        MagCutoffSizer = wx.BoxSizer(wx.VERTICAL)
        MagCutoffSizer.Add(magCellSizer)
        MagCutoffSizer.Add(cutoffSizer)
        
        topSizer = wx.FlexGridSizer(cols = 2, hgap = 5, vgap = 5)
        topSizer.Add(spaceGroupSizer)
        topSizer.Add(self.genButton, 0, wx.ALIGN_RIGHT)
        topSizer.Add(leftTopSizer)
        topSizer.Add(MagCutoffSizer)

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(topSizer, 0)
        sizer.Add(self.atomList, 1, wx.EXPAND)
        
        self.SetSizer(sizer)
        self.Fit()
        self.GetParent().Fit()
        self.GetParent().SetMinSize(self.GetParent().GetSize())
        
        #execute the function OnGenerate when the generate button is pressed
        self.Bind(wx.EVT_BUTTON, self.OnGenerate, self.genButton)
        connect(self.OnFileLoad, signal = "File Load")

    
    def OnFileLoad(self, spaceGroup, a, b, c, alpha, beta, gamma, magNa, magNb, magNc, cutNa, cutNb, cutNc):
        """This is run when a message is received from the session that a file
        was loaded.  All the values in the panel are changed to reflect the
        values that were read in."""
        self.spaceGroupSpinner.SetValue(spaceGroup)
        self.aText.SetValue(a)
        self.bText.SetValue(b)
        self.cText.SetValue(c)
        self.alphaText.SetValue(alpha)
        self.betaText.SetValue(beta)
        self.gammaText.SetValue(gamma)
        self.naText.SetValue(magNa.__str__())
        self.nbText.SetValue(magNb.__str__())
        self.ncText.SetValue(magNc.__str__())
        self.cutoffNaText.SetValue(cutNa.__str__())
        self.cutoffNbText.SetValue(cutNb.__str__())
        self.cutoffNcText.SetValue(cutNc.__str__())
        self.atomSpinner.SetValue(self.atomList.GetNumberRows())
        #self.atomList.ForceRefresh()
        self.atomList.Refresh()
        self.atomList.AutoSize()
        self.Fit()
        self.GetParent().Fit()
#        for i in range(len(atomData)):
#            self.atomList.SetCellValue(i, 0, atomData[i][0].__str__())
#            self.atomList.SetCellValue(i, 1, atomData[i][1].__str__())
#            self.atomList.SetCellValue(i, 2, atomData[i][2].__str__())
#            self.atomList.SetCellValue(i, 3, atomData[i][3].__str__())
#            self.atomList.SetCellValue(i, 4, atomData[i][4].__str__())
        
    
    def OnGenerate(self, event):
        """This is executed when the generate button is pressed.  First the data
        entered by the used is validated to make sure it is all the correct
        type.  Then session.cellChange is called to change the model."""
        failed, a, b, c, alpha, beta, gamma, magNa, magNb, magNc, cutNa, cutNb, cutNc, atomData = self.validate()
        if failed:
            return
        spaceGroup = self.spaceGroupSpinner.GetValue() #int
#        send(signal = "Cell Change", sender = "Cell Panel",
#             spaceGroup = spaceGroup,
#             a = a, b = b, c = c,
#             alpha = alpha, beta = beta, gamma = gamma,
#             magNa = magNa, magNb = magNb, magNc = magNc,
#             cutNa = cutNa, cutNb = cutNb, cutNc = cutNc,
#             atomData = atomData)
        self.session.cellChange(spaceGroupInt = spaceGroup,
             a = a, b = b, c = c,
             alpha = alpha, beta = beta, gamma = gamma,
             magNa = magNa, magNb = magNb, magNc = magNc,
             cutNa = cutNa, cutNb = cutNb, cutNc = cutNc,
             atomData = atomData)

    
    def validate(self):
        """Checks that all values are the right type

        Any field that is not of the right type will be turned pink."""
        a = self.aText.GetValue() #str - float
        b = self.bText.GetValue() #str - float
        c = self.cText.GetValue() #str - float
        alpha = self.alphaText.GetValue() #str float
        beta = self.betaText.GetValue() #str float
        gamma = self.gammaText.GetValue() #str float
        magNa = self.naText.GetValue() #int
        magNb = self.nbText.GetValue() #int
        magNc = self.ncText.GetValue() #int
        cutNa = self.cutoffNaText.GetValue()
        cutNb = self.cutoffNbText.GetValue()
        cutNc = self.cutoffNcText.GetValue()
 
         
        bgColor = "pink"
        failed = False
        #Validate a(must be a float)
        numA = None
        try:
            numA = float(a)
            self.aText.SetBackgroundColour("white")
        except:
            self.aText.SetBackgroundColour(bgColor)
            failed = True
        
        #Validate b(must be a float)
        numB = None
        try:
            numB = float(b)
            self.bText.SetBackgroundColour("white")
        except:
            self.bText.SetBackgroundColour(bgColor)
            failed = True
        
        #Validate c(must be a float)
        numC = None
        try:
            numC = float(c)
            self.cText.SetBackgroundColour("white")
        except:
            self.cText.SetBackgroundColour(bgColor)
            failed = True
        
        #Validate alpha(must be a float)
        numAlpha = None
        try:
            numAlpha = float(alpha)
            self.alphaText.SetBackgroundColour("white")
        except:
            self.alphaText.SetBackgroundColour(bgColor)
            failed = True
        
        #Validate beta(must be a float)
        numBeta = None
        try:
            numBeta = float(beta)
            self.betaText.SetBackgroundColour("white")
        except:
            self.betaText.SetBackgroundColour(bgColor)
            failed = True
         
        #Validate gamma(must be a float)
        numGamma = None
        try:
            numGamma = float(gamma)
            self.gammaText.SetBackgroundColour("white")
        except:
            self.gammaText.SetBackgroundColour(bgColor)
            failed = True
        
        
        #Validate Magnetic Cell Na(must be an int)
        numMagNa = None
        try:
            numMagNa = int(magNa)
            self.naText.SetBackgroundColour("white")
        except:
            self.naText.SetBackgroundColour(bgColor)
            failed = True
            
        #Validate Magnetic Cell Nb(must be an int)
        numMagNb = None
        try:
            numMagNb = int(magNb)
            self.nbText.SetBackgroundColour("white")
        except:
            self.nbText.SetBackgroundColour(bgColor)
            failed = True
            
        #Validate Magnetic Cell Nc(must be an int)
        numMagNc = None
        try:
            numMagNc = int(magNc)
            self.ncText.SetBackgroundColour("white")
        except:
            self.ncText.SetBackgroundColour(bgColor)
            failed = True
        
        #Validate cutoff Na(must be a int)
        numCutNa = None
        try:
            numCutNa = int(cutNa)
            self.cutoffNaText.SetBackgroundColour("white")
        except:
            self.cutoffNaText.SetBackgroundColour(bgColor)
            failed = True
        
        #Validate cutoff Nb(must be a int)
        numCutNb = None
        try:
            numCutNb = int(cutNb)
            self.cutoffNbText.SetBackgroundColour("white")
        except:
            self.cutoffNbText.SetBackgroundColour(bgColor)
            failed = True
        
        #Validate cutoff Nc(must be a int)
        numCutNc = None
        try:
            numCutNc = int(cutNc)
            self.cutoffNcText.SetBackgroundColour("white")
        except:
            self.cutoffNcText.SetBackgroundColour(bgColor)
            failed = True
        
        
        #Validate atom data in table
        data = []
        for row in range(self.atomList.GetNumberRows()):
            
            atomicNum = None
            try:
                atomicNum = int(self.atomList.GetCellValue(row, 1))
                attr = wx.grid.GridCellAttr()
                attr.SetBackgroundColour("white")
                self.atomList.SetAttr(row, 1, attr)
            except:
                attr = wx.grid.GridCellAttr()
                attr.SetBackgroundColour(bgColor)
                self.atomList.SetAttr(row, 1, attr)
                failed = True
            
            #Valence at 2
            valence = self.atomList.GetCellValue(row, 2)
            
            numXCoord = None
            try:
                numXCoord = float(self.atomList.GetCellValue(row, 3))
                attr = wx.grid.GridCellAttr()
                attr.SetBackgroundColour("white")
                self.atomList.SetAttr(row, 3, attr)
            except:
                attr = wx.grid.GridCellAttr()
                attr.SetBackgroundColour(bgColor)
                self.atomList.SetAttr(row, 3, attr)
                failed = True
        
            numYCoord = None
            try:
                numYCoord = float(self.atomList.GetCellValue(row, 4))
                attr = wx.grid.GridCellAttr()
                attr.SetBackgroundColour("white")
                self.atomList.SetAttr(row, 4, attr)
            except:
                attr = wx.grid.GridCellAttr()
                attr.SetBackgroundColour(bgColor)
                self.atomList.SetAttr(row, 4, attr)
                failed = True
        
        
            numZCoord = None
            try:
                numZCoord = float(self.atomList.GetCellValue(row, 5))
                attr = wx.grid.GridCellAttr()
                attr.SetBackgroundColour("white")
                self.atomList.SetAttr(row, 5, attr)
            except:
                attr = wx.grid.GridCellAttr()
                attr.SetBackgroundColour(bgColor)
                self.atomList.SetAttr(row, 5, attr)
                failed = True
                
            
            numDx = None
            try:
                numDx = float(self.atomList.GetCellValue(row,6))
                attr = wx.grid.GridCellAttr()
                attr.SetBackgroundColour("white")
                self.atomList.SetAttr(row, 6, attr)
            except:
                attr = wx.grid.GridCellAttr()
                attr.SetBackgroundColour(bgColor)
                self.atomList.SetAttr(row, 6, attr)
                failed = True
                
            
            numDy = None
            try:
                numDy = float(self.atomList.GetCellValue(row,7))
                attr = wx.grid.GridCellAttr()
                attr.SetBackgroundColour("white")
                self.atomList.SetAttr(row, 7, attr)
            except:
                attr = wx.grid.GridCellAttr()
                attr.SetBackgroundColour(bgColor)
                self.atomList.SetAttr(row, 7, attr)
                failed = True
                
            
            numDz = None
            try:
                numDz = float(self.atomList.GetCellValue(row,8))
                attr = wx.grid.GridCellAttr()
                attr.SetBackgroundColour("white")
                self.atomList.SetAttr(row, 8, attr)
            except:
                attr = wx.grid.GridCellAttr()
                attr.SetBackgroundColour(bgColor)
                self.atomList.SetAttr(row, 8, attr)
                failed = True
                
            spinMag = None
            try:
                spinMag = float(self.atomList.GetCellValue(row, 9))
                attr = wx.grid.GridCellAttr()
                attr.SetBackgroundColour("white")
                self.atomList.SetAttr(row, 9, attr)
            except:
                attr = wx.grid.GridCellAttr()
                attr.SetBackgroundColour(bgColor)
                self.atomList.SetAttr(row, 9, attr)
                failed = True
            
            name = self.atomList.GetCellValue(row, 0)
            
            
            self.atomList.AutoSize()  #There may be a better way to do this, but
            #this rerenders the cells to they show the color change
            
            
            data.append([name, atomicNum, numXCoord, numYCoord, numZCoord, numDx, numDy, numDz, spinMag, valence])
        
        return failed, numA, numB, numC, numAlpha, numBeta, numGamma, numMagNa, numMagNb, numMagNc, numCutNa, numCutNb, numCutNc, data
         
    
    
    def OnGridResize(self, event):
        """Resizes the grid when the spinner value is changed."""
        rows = self.atomSpinner.GetValue()
        self.atomList.SetNumberRows(rows)
#        self.atomList.GetTable().SetNumberRows(rows)
        self.Fit()
        self.GetParent().Fit()

        event.Skip()


class atomListGrid(wx.grid.Grid):
    """This is the table of atom values.  It displays values in the atom table
    stored by the session."""
    def __init__(self, parent, id, session):
        wx.grid.Grid.__init__(self, parent, id)
        self.session = session
        self.SetTable(self.session.getAtomTable())
        self.table = session.getAtomTable()
        self.AutoSize()
    def SetNumberRows(self, num):
        diff = num - self.table.GetNumberRows()
        if diff > 0:
            for i in range(diff):
                self.table.AppendRow() #Add a blank row
        elif diff < 0:
            self.table.DeleteRows(self.table.GetNumberRows()+diff, -diff)
#        self.AutoSize()
        return diff


class vtkPanel(wx.Panel):
    """This is a the main panel which displays the 3D vtk rendering."""
    def __init__(self, parent, id, session):
        wx.Panel.__init__(self, parent)
        self.session = session
        self.initVTKWindow()
        self.bindEvents()
        self.mode = None
        self.picker = None
        self.procManager = ProcessManager(self)
  
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
#        self.window.Render()
        
#        self.initializeVTKData()
#        self.draw
   
    def bindEvents(self):
        self.window.Bind(wx.EVT_KEY_DOWN, self.OnKeyEvent)
#        connect(self.OnCellChange, signal = "Cell Change")
#        connect(self.OnBondChange, signal = "Bond Change")
        connect(self.OnPick, signal = "Pick Event")
        connect(self.OnModelChange, signal = "Model Change")
        connect(self.OnSaveImage, signal = "Save Image")
        connect(self.OnAnalyticDispCalc, signal = "Analytic Dispersion Complete")
    

    def OnAnalyticDispCalc(self, answer):
        eig_frame = printing.LaTeXDisplayFrame(self, answer, 'Dispersion Eigenvalues')
        eig_frame.Show()
        
    def OnSaveImage(self, path):
        """Saves a tiff image of the current screen."""
        w2i = vtkWindowToImageFilter()
        w2i.SetInput(self.window.GetRenderWindow())
        tiffWriter = vtkTIFFWriter()
        tiffWriter.SetInput(w2i.GetOutput())
        tiffWriter.SetFileName(path)
        tiffWriter.Write()
        
    def OnModelChange(self):
        self.draw()
    
    def OnKeyEvent(self, event):
        """Handles the Del Key.  If a bond is selected and hte Del key is
        pressed, the bond and all symmetry equivalent bonds are deleted."""
        event.Skip()
        
        #Delete if key pressed is Del key  Del = 127
        if event.GetKeyCode() == 127:
            self.OnDelete(event)
         
    def OnDelete(self,event):
        """Handles deleteion of a bond and all it's symmetry equivalent bonds."""
        event.Skip()
        selectedObj = self.drawer.getObjFromActor(self.picker.getPicked())
        if isinstance(selectedObj, Bond):
            self.session.MagCell.deleteBond(selectedObj)
        
        self.draw()
    
           
    def draw(self, progDialog = None):
        """Re-renders the vtkRender window to reflect any changes."""
        #remove old renderer
#        if self.ren1 != None:
 #           self.window.GetRenderWindow().RemoveRenderer(self.ren1)
        
#        a = wx.ClientDC(self.window)
        
        #Create a progress bar
        destroy = False #So that a pre-existing progress bar could be used and Destroyed by the calling function
        if not progDialog:
            progDialog = wx.ProgressDialog("Progress", "Rendering Model...", parent = self, style = wx.PD_APP_MODAL | wx.PD_AUTO_HIDE)
            destroy = True

        
        # a renderer for the data
        ren1 = vtkRenderer()
        ren1.SetBackground(1,1,1)
        ren1.SetAllocatedRenderTime(3)
        
        progDialog.Update(5, "Rendering Model...")
        
        #Add the renderer to the window
        self.window.GetRenderWindow().AddRenderer(ren1)
            
        #Create vtkDrawer
        self.drawer = vtkDrawer(ren1)
           
        #Set up trackball mode
        interactor = vtk.vtkInteractorStyleSwitch()
        interactor.SetCurrentStyleToTrackballCamera()
        self.window._Iren.SetInteractorStyle(interactor)
        
        #Add my picker
        if self.picker:
            self.picker.removeObserver()
        self.picker = Picker(self.drawer, self.window._Iren, ren1)

        progDialog.Update(10) #This is a very rough estimation

        #Draw the Magnetic Cell
#        self.drawer.drawMagneticCell(self.session.getMagneticCell())
        self.drawer.drawCutoffCell(self.session.getCutoffCell())
        
        progDialog.Update(30)
        
        self.window.setUpRender()   
        

        self.drawer.addAxes()
        self.drawer.labelAtoms(self.session.getMagneticCell())
        
        
        progDialog.Update(100)

        if destroy:
            progDialog.Destroy()

        #Rendering does not work when the window is disabled which it seems
        #to be when the progress dialog exits

        ren1.ResetCamera()
#        self.window.Render()
#        self.window.SetRenderWhenDisabled(True)
        self.window.setUpRender()
#        self.window.SetRenderWhenDisabled(False)


    def OnPick(self, obj):
        if self.mode:
            self.mode.OnPick(obj)
#        self.mode = None

            
    def OnChooseBondMode(self, event):
        self.mode = BondMode()
        
    def OnChooseNormalMode(self, event):
        self.mode = None

class BondMode():
    """This class handles pick events in bond mode, by creating an interaction
    when two atoms are clicked."""
    def __init__(self):
        self.atom1 = None
    
    def OnPick(self, obj):
        if isinstance(obj, Atom):
            if self.atom1==None:
                self.atom1 = obj
            else:
                send(signal = "Bond Added", sender = "VTK Window", atom1 = self.atom1, atom2 = obj)
                self.atom1 = None


#def chi_sq(x,fitter,err,domain):
#    chi = 0
#    fitter.fit_list = x
#    calc_vals = fitter.GetResult()
#    print "\neigen values = ", calc_vals
#    for i in range(len(fitter.fit_list)):
#        diff = calc_vals[i] - (4 - 4*cos(domain[i][2]))
#        diff=diff.evalf()#test
#        print diff
#        num = diff**2
#        den = err[i]**2
#        print den
#        chi += num#/den
#    print "\nChi^2 = ", chi
#    return chi


class Frame(wx.Frame):
    """This is the main frame containing the vtkPanel."""
    def __init__(self, parent, id, session):
        wx.Frame.__init__(self, parent, id, 'Magnetic Cell', size= (700,700))

        self.session = session
        self.vtkPanel = vtkPanel(self, -1, session)
        self.procManager = self.vtkPanel.procManager
        
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
#        newMenuItem = fileMenu.Append(wx.NewId(), "&New Magnetic Cell")
        openMenuItem = fileMenu.Append(wx.NewId(), "&Open")
        saveMenuItem = fileMenu.Append(wx.NewId(), "&Save")
        saveImageMenuItem = fileMenu.Append(wx.NewId(), "Save Image")
        quitMenuItem = fileMenu.Append(wx.NewId(), "&Quit")
        menuBar.Append(fileMenu, "&File")
        
        #Add Monte Carlo Menu
        monteCarloMenu = wx.Menu()
        exportMenuItem = monteCarloMenu.Append(wx.NewId(), "Export for Monte Carlo")
        runSimulationMenuItem = monteCarloMenu.Append(wx.NewId(), "Launch Simulation")
        loadSpinsMenuItem = monteCarloMenu.Append(wx.NewId(), "Load Spins from file")
#        outputSnapshotsMenuItem = monteCarloMenu.Append(wx.NewId(), "Output snapshots")
        calculateSpinwavesMenuItem=monteCarloMenu.Append(wx.NewId(), "Perform Spinwave Calculation")
        crossSectionMenuItem = monteCarloMenu.Append(wx.NewId(), "Perform Cross Section Calculation")
        fittingMenuItem = monteCarloMenu.Append(wx.NewId(), "Fit Parameters")
        #This is not yet working
        #fitParametersMenuItem = monteCarloMenu.Append(wx.NewId(), "Fit Parameters")
        menuBar.Append(monteCarloMenu, "Monte Carlo")
        
        
        #Add Model Menu
        modelMenu = wx.Menu()
#        addCellMenuItem = modelMenu.Append(wx.NewId(), "Add Atom")
#        addBondMenuItem = modelMenu.Append(wx.NewId(), "Add Bond")
        deleteMenuItem = modelMenu.Append(wx.NewId(), "Delete")
        menuBar.Append(modelMenu, "Model")
        self.SetMenuBar(menuBar)
        
        #Add Mode Menu
        modeMenu = wx.Menu()
        normalModeMenuItem = modeMenu.Append(wx.NewId(), "Normal Selection")
        bondModeMenuItem = modeMenu.Append(wx.NewId(), "Bond Creation")
        menuBar.Append(modeMenu, "Mode")
        
        #Bind Events
        self.Bind(wx.EVT_MENU, self.OnCloseMe, quitMenuItem)
        self.Bind(wx.EVT_MENU, self.OnSave, saveMenuItem)
        self.Bind(wx.EVT_MENU, self.OnOpenFile, openMenuItem)
        self.Bind(wx.EVT_MENU, self.vtkPanel.OnDelete, deleteMenuItem)
        self.Bind(wx.EVT_MENU, self.vtkPanel.OnChooseNormalMode, normalModeMenuItem)
        self.Bind(wx.EVT_MENU, self.vtkPanel.OnChooseBondMode, bondModeMenuItem)
#        self.Bind(wx.EVT_MENU, self.OnNew, newMenuItem)
        self.Bind(wx.EVT_MENU, self.OnExport, exportMenuItem)
        self.Bind(wx.EVT_MENU, self.OnLoadSpins, loadSpinsMenuItem)
        self.Bind(wx.EVT_MENU, self.OnSaveImage, saveImageMenuItem)
        self.Bind(wx.EVT_MENU, self.OnLaunchSim, runSimulationMenuItem)
        self.Bind(wx.EVT_MENU, self.OnLaunchSpinWave, calculateSpinwavesMenuItem)
        self.Bind(wx.EVT_MENU, self.OnLaunchCrossSection, crossSectionMenuItem)
        self.Bind(wx.EVT_MENU, self.OnFitParameters, fittingMenuItem)
        #Doesn't work yet
#        self.Bind(wx.EVT_MENU, self.OnFitParameters, fitParametersMenuItem)
          
#        self.Bind(wx.EVT_MENU, self.createMonteCarloVideo, outputSnapshotsMenuItem)
        
    
 #   def OnNew(self, event):
        #Add drawing panel
#        self.vtkPanel.draw()
#        self.GetEventHandler().ProcessEvent(wx.SizeEvent())
    
    
    def createMonteCarloVideo(self, evt):
        def imageName(imageNum):
            """Generates the name of hte image given the image number.
            I want the images to automatically arrange in alphabetical order
            to make the video creation easier.
            
            Just using the number does not work becuase image2 comes after image10

            Therefore this will create numbers with precending 0's
            For example:

                image001
                image002
                ...
                image010
                image011"""
            imageStr = str(imageNum)
            
            val = imageNum
            totalDigits = 5 #max of 99999 pics
            for i in range(totalDigits):
                val = val / 10
                if val == 0:
                    imageStr = "0" + imageStr
            imageStr = "image" + imageStr + ".tiff"
            return imageStr
        
        
        def imageOutputFunction(spinsFile, imageNum):
            """This funstion is passed to CSim.createVideo to handle creation of
            the images."""
            imagePath = "C:\\monteCarloSnapshots\\" + imageName(imageNum)
            #Load the spins and render
            self.session.loadSpinFile(spinsFile)
            #Output to image
            send("Save Image", sender = "Main Frame", path = imagePath)
        
        CSim.createVideo(imageOutputFunction, "C:\\Spins.txt", "C:\\Export.txt")

            
    
    def OnFitParameters(self, evt):
        ShowFittingFrame(self.session, self.procManager)
#        #domain = []
#        #test case
#        domain = [(0,0,0,),
#                  (0,0,.25*pi),
#                  (0,0,.5*pi),
#                  (0,0,.75 *pi),
#                  (0,0, pi),
#                  (0,0,1.25 *pi),
#                  (0,0,1.5 *pi),
#                  (0,0,1.75 * pi)]
#        err = [1e-1]*8
#        
#        fitter = Fitter(self.session, domain)
#        #ans = fmin(chi_sq, fitter.fit_list,args=(fitter,err,domain),maxiter=30,ftol=.01,maxfun=30)
#        #lower=fitter.min_range_list
#        #upper=fitter.max_range_list
#        #lower[0]=0.8
#        #upper[0]=1.2
#        #ans=scipy.optimize.anneal(chi_sq,fitter.fit_list,args=(fitter,err,domain),
#        #                          lower=lower,upper=upper)
#        def Anneal(func, guess, min_list, max_list):
#            k = 20
#            tMax = 50
#            tMin = .01
#            tFactor = .95
#            T = tMax
#            domain = guess
#            result = func(guess)
#            while T > tMin:
#                numSwitched = 0
#                randGen = random.Random()
#                for i in range(k):
#                    #Create a new domain
#                    newVals = []
#                    for index in range(len(domain)):
#                        newVals.append(randGen.uniform(min_list[index],
#                                                       max_list[index]))
#                    newResult=func(newVals)
#                    if newResult < result:
#                        result = newResult
#                        domain = newVals
#                        numSwitched += 1
#                    else:
#                        deltE = newResult - result
#                        probChange = math.exp(-deltE/T)
#                        if randGen.random() < probChange:
#                            result = newResult
#                            domain = newVals
#                            numSwitched += 1
#                time.sleep(5)
#                T = T*tFactor
#                #print "Temp: ", T, " %flipped: ", (float(numSwitched*100)/k),"%\n"
#                #break#Temp
#            
#            return domain
#        
#        def chi_sq(x):
#            chi = 0
#            fitter.fit_list = x
#            calc_vals = fitter.GetResult()
#            print "\neigen values = ", calc_vals
#            for i in range(len(fitter.fit_list)):
#                diff = calc_vals[i] - (4 - 4*cos(domain[i][2]))
#                diff=diff.evalf()#test
#                print diff
#                num = diff**2
#                den = err[i]**2
#                print den
#                chi += num#/den
#            print "\nChi^2 = ", chi
#            return chi
#        
#        ans = Anneal(chi_sq, fitter.fit_list, fitter.min_range_list, fitter.max_range_list)
#        #ans=ans[0]
#        print ans
#        wrange=[]
#        qrange = []
#        for a in domain:
#            qrange.append(a[2])
#        wrange=(fitter.GetResult())
#        print numpy.__dict__
#        #datavals=(4 - 4*numpy.cos(numpy.array(qrange)))
#        #pylab.plot(qrange,datavals,'s')
#        pylab.plot(qrange,wrange)
#        pylab.show()
                
                
    
    def OnLaunchSim(self, evt):
        """Runs the simulation from this app."""
        CSim.ShowSimulationFrame()
        
    def OnLaunchSpinWave(self,evt):
        myparent=self
        #myparent=None
        frame1 = wx.Frame(myparent, -1, "Spinwaves")
        dlg=spinwavepanel.SpinwavePanel(procManager = self.procManager, parent=frame1, id=-1)
        #dlg=spinwavepanel.FormDialog(parent=None,id=-1)
        self.SetExtraStyle(wx.WS_EX_VALIDATE_RECURSIVELY)
        frame1.Fit()
#No longer modal
#        result=dlg.ShowModal()
#        print 'result', result
#        if result==wx.ID_OK:
#            #self.Validate()
#            dlg.Validate()
#            print "OK"
#            dlg.TransferDataFromWindow()
#            print 'data',dlg.data
#            print dlg.data['step']
#            print dlg.interactionfile
#            print dlg.spinfile
#            spinwave_calc_file.driver(dlg.spinfile,dlg.interactionfile,dlg.data,dlg.data['step'])
#        else:
#            print "Cancel"
#        dlg.Destroy()
        
        frame1.Show()
        frame1.Refresh()
        
    def OnLaunchCrossSection(self, evt):
        myparent = self
        #frame_1 = wx.Frame(myparent, -1, "Cross-section")
        #dlg = Cross_Section(parent = frame_1,id=-1)
        frame_1 = Cross_Section(self.procManager, self, -1, "")
        #self.SetExtraStyle(wx.WS_EX_VALIDATE_RECURSIVELY)
        frame_1.Show()
        frame_1.Refresh()
    
    def OnSaveImage(self, evt):
        """Saves an image of the current rendering.  Currently only .tiff
        format is used."""
        saveDialog = wx.FileDialog(self, "Save Image", style = wx.SAVE, wildcard = "*.tiff")

        #If the user clicked OK in the save dialog, use the filename they chose
        if saveDialog.ShowModal() == wx.ID_OK:
            send("Save Image", sender = "Main Frame", path = saveDialog.GetPath())
        saveDialog.Destroy()
    
    def OnExport(self, evt):
        """Exports the interactions to a file for use in the monte carlo
        simulation to find the ground state."""
        #Maximum size is set to 25 right now, becuase that is the largest size 
        #that this computer seems to be able to reasonably handle with the current algorithm
        #The minimum size has been set to 3, so that a completely surrounded
        #Interaction cell can always be selected for hte spinwave calculation
        size = wx.GetNumberFromUser("How many times would you like to translate the cutoff cell in the a,b, and c directions?", prompt = "size:", caption = "Monte Carlo Simulation Size", value = 3, min = 3, max=25, parent = self)
        if size != None:
            saveDialog = wx.FileDialog(self, "Save File", style = wx.SAVE, wildcard = "*.txt")
            if saveDialog.ShowModal() == wx.ID_OK:
                self.session.exportForMonteCarlo(saveDialog.GetPath(), size)
            saveDialog.Destroy()
        else:
            print None
            
    
    def OnCloseMe(self, event):
        self.Close(True)
    
    def OnSave(self, event):
        saveDialog = wx.FileDialog(self, "Save File", style = wx.SAVE, wildcard = "*.xml")
        if saveDialog.ShowModal() == wx.ID_OK:
            self.session.saveSessionToXML(saveDialog.GetPath())
        saveDialog.Destroy()
        
    def OnOpenFile(self, event):
        openDialog = wx.FileDialog(self, "Open File", style = wx.OPEN, wildcard = "XML Session (*.xml)|*.xml|Crystallographic Information File (*.cif)|*.cif")
        if openDialog.ShowModal() == wx.ID_OK:
            index = openDialog.GetFilterIndex()
            if index == 0: #If they chose .xml format
                self.session.openXMLSession(openDialog.GetPath())
            if index == 1: #If they chose .cif format
                self.session.openCif(openDialog.GetPath())
        openDialog.Destroy()
    
    def OnLoadSpins(self, event):
        """Loads a text (.txt) file with a list of spins at atom positions.
        The monte Carlo simulation outputs a spin file like this."""
        openDialog = wx.FileDialog(self, "Open Spin File", style = wx.OPEN, wildcard = "*.txt")
        if openDialog.ShowModal() == wx.ID_OK:
            self.session.loadSpinFile(openDialog.GetPath())
        openDialog.Destroy()



class Cross_Section(wx.Frame):
    def __init__(self, processManager, *args, **kwds):
    #def __init__(self, parent, id):
        #wx.Panel.__init__(self, parent, id)
         #begin wxGlade: Cross_Section.__init__
        kwds["style"] = wx.DEFAULT_FRAME_STYLE
        wx.Frame.__init__(self, *args, **kwds)
        self.sizer_4_staticbox = wx.StaticBox(self, -1, "Spins File")
        self.sizer_5_staticbox = wx.StaticBox(self, -1, "Interactions File")
        self.text_ctrl_1 = wx.TextCtrl(self, -1, "")
        self.interactionsFileBrowse = wx.Button(self, -1, "Browse")
        self.text_ctrl_2 = wx.TextCtrl(self, -1, "")
        self.button_1 = wx.Button(self, -1, "Browse")
        self.launchButton = wx.Button(self, -1, "Launch")
        #self.process_list = []
        self.procManager = processManager

        
        
        #self.sizer_4_staticbox = wx.StaticBox(self, -1, "Spins File")
        #self.sizer_5_staticbox = wx.StaticBox(self, -1, "Interactions File")
        #self.text_ctrl_1 = wx.TextCtrl(self, -1, "")
        #self.interactionsFileBrowse = wx.Button(self, -1, "Browse")
        #self.text_ctrl_2 = wx.TextCtrl(self, -1, "")
        #self.button_1 = wx.Button(self, -1, "Browse")
        #self.launchButton = wx.Button(self, -1, "Launch")
        #self.process_list = []
        
        ##Fit this window
        #self.Fit()
        #self.GetParent().Fit()#Fit the frame containing this panel
        #self.GetParent().SetMinSize(self.GetParent().GetSize())

        self.__set_properties()
        self.__do_layout()
        # end wxGlade

    def __set_properties(self):
        # begin wxGlade: Cross_Section.__set_properties
        self.SetTitle("Cross-Section")
        self.text_ctrl_1.SetMinSize((160, 27))
        self.text_ctrl_2.SetMinSize((160, 27))
        # end wxGlade

    def __do_layout(self):
        # begin wxGlade: Cross_Section.__do_layout
        sizer_2 = wx.BoxSizer(wx.VERTICAL)
        sizer_4 = wx.StaticBoxSizer(self.sizer_4_staticbox, wx.HORIZONTAL)
        sizer_5 = wx.StaticBoxSizer(self.sizer_5_staticbox, wx.HORIZONTAL)
        sizer_5.Add(self.text_ctrl_1, 0, 0, 0)
        sizer_5.Add(self.interactionsFileBrowse, 0, 0, 0)
        sizer_2.Add(sizer_5, 1, wx.EXPAND, 0)
        sizer_4.Add(self.text_ctrl_2, 0, 0, 0)
        sizer_4.Add(self.button_1, 0, 0, 0)
        sizer_2.Add(sizer_4, 1, wx.EXPAND, 0)
        sizer_2.Add(self.launchButton, 0, wx.ALIGN_CENTER_HORIZONTAL, 0)
        self.SetSizer(sizer_2)
        sizer_2.Fit(self)
        self.Layout()
        # end wxGlade
        self.Bind(wx.EVT_BUTTON, self.OnBrowseInteractions, self.interactionsFileBrowse)
        self.Bind(wx.EVT_BUTTON, self.OnBrowseSpins, self.button_1)
        self.Bind(wx.EVT_BUTTON, self.OnLaunch, self.launchButton)

    #def OnCloseWindow(self, evt):
        #for pro in self.process_list:
            #pro.terminate()
        #self.Destroy()
        
    def OnBrowseInteractions(self, evt):
        confBase = wx.ConfigBase.Create()
        confBase.SetStyle(wx.CONFIG_USE_LOCAL_FILE)
        defaultDir=confBase.Get().GetPath()
        wildcard="files (*.txt)|*.txt|All files (*.*)|*.*"
        dlg = wx.FileDialog(
            self, message="Choose an interaction file",
            defaultDir=defaultDir,
            defaultFile="",
            wildcard=wildcard,
            style=wx.OPEN | wx.CHANGE_DIR
            )

        # Show the dialog and retrieve the user response. If it is the OK response,
        # process the data.
        if dlg.ShowModal() == wx.ID_OK:
            # This returns a Python list of files that were selected.
            paths = dlg.GetPaths()
            self.text_ctrl_1.SetValue(paths[0])
            
    
    def OnBrowseSpins(self, evt):
        confBase = wx.ConfigBase.Create()
        confBase.SetStyle(wx.CONFIG_USE_LOCAL_FILE)
        defaultDir=confBase.Get().GetPath()
        wildcard="files (*.txt)|*.txt|All files (*.*)|*.*"
        dlg = wx.FileDialog(
            self, message="Choose a spin file",
            defaultDir=defaultDir,
            defaultFile="",
            wildcard=wildcard,
            style=wx.OPEN | wx.CHANGE_DIR
            )

        # Show the dialog and retrieve the user response. If it is the OK response,
        # process the data.
        if dlg.ShowModal() == wx.ID_OK:
            # This returns a Python list of files that were selected.
            paths = dlg.GetPaths()
            self.text_ctrl_2.SetValue(paths[0])
    
    def OnLaunch(self, evt):
        #check to see if the two files paths are openable files.
        try:
            f = open(self.text_ctrl_1.GetValue())
            f.close()
        except:
            wx.MessageBox("The interactions file cannot be opened.")
            return
        
        
        try:
            f = open(self.text_ctrl_2.GetValue())
            f.close()
        except:
            wx.MessageBox("The spins file cannot be opened.")
            return
        
        self.procManager.startAnalyticCrossSection(self.text_ctrl_1.GetValue(), self.text_ctrl_2.GetValue())
        #N_atoms_uc,csection,kaprange,qlist,tau_list,eig_list,kapvect,wtlist = run_cross_section(self.text_ctrl_1.GetValue(), self.text_ctrl_2.GetValue())
        
        #print 'create pipe'
        #left_conn, right_conn = Pipe()
        #p = Process(target = printing.create_latex, args = (right_conn, csection, "eigs"))
        #print 'process starting'
        #p.start()
        #print 'process started'
        #self.process_list.append(p)

        #printing.process_info('main line')
        ##q = Process(target = run_eval_cross_section, args = (N_atoms_uc,csection,kaprange,qlist,tau_list,eig_list,kapvect,wtlist))
        #run_eval_cross_section(N_atoms_uc,csection,kaprange,qlist,tau_list,eig_list,kapvect,wtlist)
        ##q.start()
        
        #p.join()
        #print 'displaying window'
        #eig_frame = printing.LaTeXDisplayFrame(None, p.pid, left_conn.recv(), 'Cross-section')
        #if eig_frame.PID == p.pid:
            #eig_frame.Show()
            #self.process_list.remove(p)
            #p.terminate()
        #else:
            #if p in self.process_list:
                #p.terminate()
                #raise Exception('process messed up')
        #p.terminate()
        #print 'p terminated'
        ##q.join()
        ##print 'q done'
        ##q.terminate()
        #print 'q terminated'
# end of class Cross_Section



class App(wx.App):
    def __init__(self, redirect = False, filename = None):
        wx.App.__init__(self, redirect, filename)
    
    def OnInit(self):
        #Start a new session
        session = Session()
        #Create the main frame
        self.frame = Frame(None, -1, session = session)
        self.frame.Show()
        self.SetTopWindow(self.frame)
        #Create the atom frame
        frame1 = wx.Frame(self.frame, -1, "Atoms")
        #frame1.SetMinSize((500,245))
        atomPanel(frame1, -1, session = session)
        frame1.Show()
        #Create the bond frame
        frame2 = wx.Frame(self.frame, -1, 'Bonds')
        #frame2.SetMinSize((655, 140))
        bondPanel(frame2, -1, session = session)
        #bondPanel.Fit()
        #frame2.Fit()
        #frame2.SetMinSize(frame2.GetSize())
        frame2.Show()

        return True


def main():

    app = App(False)
    app.MainLoop()

if __name__ == '__main__':
    main()

    
    


