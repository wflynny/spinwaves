#!/usr/bin/env python

"""This is the main GUI."""

import sys
import os
import time

import wx
import wx.grid
from wx.py.dispatcher import connect, send

from picker import Picker
from wxVTKRenderWindowInteractor import *
from vtkModel.SpaceGroups import *
from vtkModel.VTKDrawer import *
from vtkModel.MagneticCellClass import *
from vtkModel.CellClass import *
#import random
import vtkModel.SpaceGroups
from Session import Session

#It could not find MonteCarlo package (import MonteCarlo.CSim)
#sys.path.append(mainPath +"\\MonteCarlo")
import MonteCarlo.CSim as CSim
import spinwavecalc.spinwavepanel as spinwavepanel
import spinwavecalc.spinwave_calc_file as spinwave_calc_file

from vtkModel.BondClass import JParam

#Atom and cell info window

class atomPanel(wx.Panel):
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
        self.atomList.ForceRefresh()
        self.atomList.AutoSize()
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

        Any field that is not of the right type will be turned pink.  The
        textCtrl will turn the backcolor of the writting pink, while whole cells
        in the grid will turn pink.  Therefore, if there is no writting in a
        textCtrl field, there will be no color either unitl the user types."""
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
                atomicNum = float(self.atomList.GetCellValue(row, 1))
                attr = wx.grid.GridCellAttr()
                attr.SetBackgroundColour("white")
                self.atomList.SetAttr(row, 1, attr)
            except:
                attr = wx.grid.GridCellAttr()
                attr.SetBackgroundColour(bgColor)
                self.atomList.SetAttr(row, 1, attr)
                failed = True
            
            
            numXCoord = None
            try:
                numXCoord = float(self.atomList.GetCellValue(row, 2))
                attr = wx.grid.GridCellAttr()
                attr.SetBackgroundColour("white")
                self.atomList.SetAttr(row, 2, attr)
            except:
                attr = wx.grid.GridCellAttr()
                attr.SetBackgroundColour(bgColor)
                self.atomList.SetAttr(row, 2, attr)
                failed = True
        
            numYCoord = None
            try:
                numYCoord = float(self.atomList.GetCellValue(row, 3))
                attr = wx.grid.GridCellAttr()
                attr.SetBackgroundColour("white")
                self.atomList.SetAttr(row, 3, attr)
            except:
                attr = wx.grid.GridCellAttr()
                attr.SetBackgroundColour(bgColor)
                self.atomList.SetAttr(row, 3, attr)
                failed = True
        
        
            numZCoord = None
            try:
                numZCoord = float(self.atomList.GetCellValue(row, 4))
                attr = wx.grid.GridCellAttr()
                attr.SetBackgroundColour("white")
                self.atomList.SetAttr(row, 4, attr)
            except:
                attr = wx.grid.GridCellAttr()
                attr.SetBackgroundColour(bgColor)
                self.atomList.SetAttr(row, 4, attr)
                failed = True
                
            
            numDx = None
            try:
                numDx = float(self.atomList.GetCellValue(row,5))
                attr = wx.grid.GridCellAttr()
                attr.SetBackgroundColour("white")
                self.atomList.SetAttr(row, 5, attr)
            except:
                attr = wx.grid.GridCellAttr()
                attr.SetBackgroundColour(bgColor)
                self.atomList.SetAttr(row, 5, attr)
                failed = True
                
            
            numDy = None
            try:
                numDy = float(self.atomList.GetCellValue(row,6))
                attr = wx.grid.GridCellAttr()
                attr.SetBackgroundColour("white")
                self.atomList.SetAttr(row, 6, attr)
            except:
                attr = wx.grid.GridCellAttr()
                attr.SetBackgroundColour(bgColor)
                self.atomList.SetAttr(row, 6, attr)
                failed = True
                
            
            numDz = None
            try:
                numDz = float(self.atomList.GetCellValue(row,7))
                attr = wx.grid.GridCellAttr()
                attr.SetBackgroundColour("white")
                self.atomList.SetAttr(row, 7, attr)
            except:
                attr = wx.grid.GridCellAttr()
                attr.SetBackgroundColour(bgColor)
                self.atomList.SetAttr(row, 7, attr)
                failed = True
            
            
            name = self.atomList.GetCellValue(row, 0)
            
            
            self.atomList.AutoSize()  #There may be a better way to do this, but
            #this rerenders the cells to they show the color change
            
            
            data.append([name, atomicNum, numXCoord, numYCoord, numZCoord, numDx, numDy, numDz])
        
        return failed, numA, numB, numC, numAlpha, numBeta, numGamma, numMagNa, numMagNb, numMagNc, numCutNa, numCutNb, numCutNc, data
         
    
    
    def OnGridResize(self, event):
        """Resizes the grid when the spinner value is changed."""
        rows = self.atomSpinner.GetValue()
        self.atomList.SetNumberRows(rows)
#        self.atomList.GetTable().SetNumberRows(rows)
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


#Bond information window
class bondListGrid(wx.grid.Grid):
    """This is the table of bonds displayed in the bond panel.  It displays
    values stored in the bond table stored by the session."""
    def __init__(self, parent, id, session):
        wx.grid.Grid.__init__(self, parent, id)
        self.table = session.getBondTable()
        self.SetTable(self.table)
        self.AutoSize()
        
        #Set up last cell for clicks only
#        wx.grid.Grid.EnableEditing(self,False)
        attr=wx.grid.GridCellAttr()
        attr.SetReadOnly(True)
        self.SetColAttr(9,attr)
        self.Bind(wx.grid.EVT_GRID_CELL_LEFT_CLICK, self.OnLeftClick ,self)
        
        #to toggle view between parameter values and names
        self.Bind(wx.grid.EVT_GRID_CMD_LABEL_LEFT_CLICK, self.OnLabelClick, self)
        
    
    def OnLabelClick(self, evt):
        #True = value is displayed, False = name is displayed
        self.table.valueView = not self.table.valueView
        #For now the event will not be skipped because highlighting of rows or columns
        #is unnecessary
        self.Refresh()
        
    def SetNumberRows(self, num):
        diff = num - self.table.GetNumberRows()
        if diff > 0:
            for i in range(diff):
                self.table.AppendRow() #Add a blank row
        elif diff < 0:
            self.table.DeleteRows(self.table.GetNumberRows()+diff, -diff)
#        self.AutoSize()
        return diff
    
    def OnLeftClick(self, evt):
        """If there is a left lick in the rightmost cell(the cell that turns the
        bond on or off), an X is placed there to represent that is is on, or
        removed to show that it is off.

        If there is a click in column 8(Jij Matrix), then a dialog is opened to
        allow the user to enter a Jij Matrix."""
        evt.Skip()
        col=evt.GetCol()
        row=evt.GetRow()

        if col>=9 and row >=0: # ON/Off
            currval=self.table.GetValue(row,9)
            if currval=='':
                self.table.SetValue(row,9,'X')
            else:
                self.table.SetValue(row,9,'')
        elif col==8 and row >=0:  #Jij matrix
            dialog = jijDialog(self.table.GetActualValue(row,8))#Pass current Jij value
            result = dialog.ShowModal()
            if result == wx.ID_OK:
                #self.table.SetValue(row, 8, numpy.array(dialog.getMatrix()))
                dialog.setMatrix()
                self.Refresh()#values should have been updated

            dialog.Destroy()

        self.AutoSize()
        #if self.CanEnableCellControl():
        #    self.EnableCellEditControl()
#        wx.grid.Grid.ForceRefresh(self)

class bondPanel(wx.Panel):
    """This panel allows the user to create bonds."""
    def __init__(self, parent, id, session):
        print 1
        wx.Panel.__init__(self, parent, id)
        print "initializing bondtable"
        
        self.session = session
        
        #Create the table of bonds
        self.bondList = bondListGrid(self, -1, session)
        print "set bondlistgrid"
        
        #Create the spinner which controls the length of the bond list
        self.bondSpinner = wx.SpinCtrl(self, -1, "")
        self.bondSpinner.SetRange(1,100)
        self.bondSpinner.SetValue(1)
        self.bondSpinner.Bind(wx.EVT_TEXT, self.OnGridResize, self.bondSpinner)
        
        #Create the generate button
        self.genButton = wx.Button(self, -1, "Generate")

        topSizer = wx.BoxSizer(wx.HORIZONTAL)
        topSizer.Add(self.bondSpinner)
        topSizer.Add(self.genButton)
        
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(topSizer)
        sizer.Add(self.bondList, 1, wx.EXPAND)
        self.SetSizer(sizer)
        
        self.Bind(wx.EVT_BUTTON, self.OnGenerate, self.genButton)
        connect(self.OnBondAddition, signal = "Bond Added")
        connect(self.OnFileLoad, signal = "File Load")

        #Create the parameter editing frame
        paramFrame = wx.Frame(self.GetParent(), -1, "Parameters")
        param_panel = ParameterPanel(self.bondList.table, paramFrame, -1)
        paramFrame.Show()


    def OnFileLoad(self):
        """Executed when the session sends a message that a file was loaded."""
        self.bondSpinner.SetValue(self.bondList.GetNumberRows())
        self.bondList.AutoSize()
#    def OnFileLoad(self, bondData):
#        for i in range(len(bondData)):
#            self.bondList.SetCellValue(i, 0, bondData[i][0])
#            self.bondList.SetCellValue(i, 1, bondData[i][1])
#            self.bondList.SetCellValue(i, 2, bondData[i][2])
#            self.bondList.SetCellValue(i, 3, bondData[i][3])
#            self.bondList.SetCellValue(i, 4, bondData[i][4])
#            self.bondList.SetCellValue(i, 5, bondData[i][5])
#            self.bondList.SetCellValue(i, 6, bondData[i][6])
#            self.bondList.SetCellValue(i, 7, bondData[i][7])
#            self.bondList.SetCellValue(i, 8, bondData[i][8])
#            self.bondList.SetCellValue(i, 9, bondData[i][9])
    
    def OnBondAddition(self, atom1, atom2):
        """This method is called when another window adds a bond and calls
        send(signal = "Bond Added"..."""
        cell1 = atom1.getUnitCell()
        cell2 = atom2.getUnitCell()
        
        index1 = atom1.getIndexNumber()
        index2 = atom2.getIndexNumber()
        pos1 = cell1.getPosition()
        pos2 = cell2.getPosition()

        row = self.bondList.GetNumberRows() #Append row and add bond to last row
        self.bondList.SetCellValue(row, 0, str(index1+1))
        self.bondList.SetCellValue(row, 1, str(pos1[0]))
        self.bondList.SetCellValue(row, 2, str(pos1[1]))
        self.bondList.SetCellValue(row, 3, str(pos1[2]))
        self.bondList.SetCellValue(row, 4, str(index2+1))
        self.bondList.SetCellValue(row, 5, str(pos2[0]))
        self.bondList.SetCellValue(row, 6, str(pos2[1]))
        self.bondList.SetCellValue(row, 7, str(pos2[2]))
        self.bondList.SetCellValue(row, 9, 'X') #Turn the bond on
        
        self.bondSpinner.SetValue(self.bondList.GetNumberRows())
        
        
        self.clearEmptyRows()
#        self.OnGenerate(None)
        self.bondList.AutoSize()
        
    def clearEmptyRows(self):
        numRows = self.bondList.GetNumberRows()
        for r in range(0,numRows):
            row = numRows-1-r
            if self.RowIsEmpty(row):
                self.bondList.table.DeleteRows(row, 1)
        self.bondSpinner.SetValue(self.bondList.GetNumberRows())
        
    def RowIsEmpty(self, row):
        """Returns true if no column of the row has been filled in, other than On/Off,
        which does not matter.  Used to clear empty rows in OnBondAddition."""
        for i in range(8):
            if self.bondList.GetCellValue(row, i) != '':
                return False
        
        #if self.bondList.GetCellValue(row, 8) == numpy.array([[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]]):
        #    return False
        
        #Check if the J matrix is different from the default
        matrix = self.bondList.table.GetActualValue(row, 8)
        for i in range(3):
            for j in range(3):
                if not matrix[i][j].isDefault():
                    return False
        
        return True #(regardless if it's on or not)
            
    
    def OnGenerate(self, event):
        failed, bondData = self.validate()
        print failed
        if failed:
            return
#        send(signal = "Bond Change", sender = "Bond Panel", bondData = bondData)
        self.session.changeBonds(bondData)
      
    
    def validate(self):
        """Currently checks that all values are the right type"""
         
        bgColor = "pink"
        failed = False

         #Validate bond data in table
        bondData = []
        for row in range(self.bondList.GetNumberRows()):
            
            check = self.bondList.GetCellValue(row, 9)
            if check == 'X':
                atom1Num = None
                try:
                    atom1Num = int(self.bondList.GetCellValue(row, 0))
                    attr = wx.grid.GridCellAttr()
                    attr.SetBackgroundColour("white")
                    self.bondList.SetAttr(row, 0, attr)
                except:
                    attr = wx.grid.GridCellAttr()
                    attr.SetBackgroundColour(bgColor)
                    self.bondList.SetAttr(row, 0, attr)
                    failed = True
                
                Na1 = None
                try:
                    Na1 = int(self.bondList.GetCellValue(row, 1))
                    attr = wx.grid.GridCellAttr()
                    attr.SetBackgroundColour("white")
                    self.bondList.SetAttr(row, 1, attr)
                except:
                    attr = wx.grid.GridCellAttr()
                    attr.SetBackgroundColour(bgColor)
                    self.bondList.SetAttr(row, 1, attr)
                    failed = True
    
                Nb1 = None
                try:
                    Nb1 = int(self.bondList.GetCellValue(row, 2))
                    attr = wx.grid.GridCellAttr()
                    attr.SetBackgroundColour("white")
                    self.bondList.SetAttr(row, 2, attr)
                except:
                    attr = wx.grid.GridCellAttr()
                    attr.SetBackgroundColour(bgColor)
                    self.bondList.SetAttr(row, 2, attr)
                    failed = True
            
                Nc1 = None
                try:
                    Nc1 = int(self.bondList.GetCellValue(row, 3))
                    attr = wx.grid.GridCellAttr()
                    attr.SetBackgroundColour("white")
                    self.bondList.SetAttr(row, 3, attr)
                except:
                    attr = wx.grid.GridCellAttr()
                    attr.SetBackgroundColour(bgColor)
                    self.bondList.SetAttr(row, 3, attr)
                    failed = True
                    
                
                atom2Num = None
                try:
                    atom2Num = int(self.bondList.GetCellValue(row, 4))
                    attr = wx.grid.GridCellAttr()
                    attr.SetBackgroundColour("white")
                    self.bondList.SetAttr(row, 4, attr)
                except:
                    attr = wx.grid.GridCellAttr()
                    attr.SetBackgroundColour(bgColor)
                    self.bondList.SetAttr(row, 4, attr)
                    failed = True
                
                Na2 = None
                try:
                    Na2 = int(self.bondList.GetCellValue(row, 5))
                    attr = wx.grid.GridCellAttr()
                    attr.SetBackgroundColour("white")
                    self.bondList.SetAttr(row, 5, attr)
                except:
                    attr = wx.grid.GridCellAttr()
                    attr.SetBackgroundColour(bgColor)
                    self.bondList.SetAttr(row, 5, attr)
                    failed = True
    
                Nb2 = None
                try:
                    Nb2 = int(self.bondList.GetCellValue(row, 6))
                    attr = wx.grid.GridCellAttr()
                    attr.SetBackgroundColour("white")
                    self.bondList.SetAttr(row, 6, attr)
                except:
                    attr = wx.grid.GridCellAttr()
                    attr.SetBackgroundColour(bgColor)
                    self.bondList.SetAttr(row, 6, attr)
                    failed = True
            
                Nc2 = None
                try:
                    Nc2 = int(self.bondList.GetCellValue(row, 7))
                    attr = wx.grid.GridCellAttr()
                    attr.SetBackgroundColour("white")
                    self.bondList.SetAttr(row, 7, attr)
                except:
                    attr = wx.grid.GridCellAttr()
                    attr.SetBackgroundColour(bgColor)
                    self.bondList.SetAttr(row, 7, attr)
                    failed = True
                    

                jij = self.bondList.table.GetActualValue(row, 8)#numpy.array
                
                bondData.append([atom1Num, Na1,Nb1,Nc1, atom2Num, Na2,Nb2,Nc2, jij])
            else: #If the row is not checked, all cells should be white
                attr = wx.grid.GridCellAttr()
                attr.SetBackgroundColour("white")
                self.bondList.SetAttr(row, 0, attr)
                attr = wx.grid.GridCellAttr()
                attr.SetBackgroundColour("white")
                self.bondList.SetAttr(row, 1, attr)
                attr = wx.grid.GridCellAttr()
                attr.SetBackgroundColour("white")
                self.bondList.SetAttr(row, 2, attr)
                attr = wx.grid.GridCellAttr()
                attr.SetBackgroundColour("white")
                self.bondList.SetAttr(row, 3, attr)
                attr = wx.grid.GridCellAttr()
                attr.SetBackgroundColour("white")
                self.bondList.SetAttr(row, 4, attr)
                attr = wx.grid.GridCellAttr()
                attr.SetBackgroundColour("white")
                self.bondList.SetAttr(row, 5, attr)
                attr = wx.grid.GridCellAttr()
                attr.SetBackgroundColour("white")
                self.bondList.SetAttr(row, 6, attr)
                attr = wx.grid.GridCellAttr()
                attr.SetBackgroundColour("white")
                self.bondList.SetAttr(row, 7, attr)
                    
            
        self.bondList.AutoSize()  #There may be a better way to do this, but
        #this re-renders the cells so they show the color change
        
        return failed, bondData
        
        
    def OnGridResize(self, event):
        rows = self.bondSpinner.GetValue()
        self.bondList.SetNumberRows(rows)
        self.bondList.AutoSize()
#        self.atomList.GetTable().SetNumberRows(rows)
        event.Skip()
 

class ParamTable(wx.grid.PyGridTableBase):
    """This is the table base for the editable parameter tables in
    ParameterPanel.  It only acts as an intermediary between the bondTable,
    where the matrices of JParam objects are stored and the 3*3 parameter grid."""
    
    def __init__(self, bond_table_base, row_num):
        wx.grid.PyGridTableBase.__init__(self)
        self.bond_table = bond_table_base
        self.row_num = row_num
        
    def GetNumberRows(self):
        return 3
    
    def GetNumberCols(self):
        return 3
    
    def GetValue(self, row, col):
        """Returns a String representation of the value in the given cell."""
        return self.bond_table.GetJMatrix(self.row_num)[row][col].__str__()
        
    def SetValue(self, row, col, value):
        """Attempts to change the JParam in the bond table base."""
        self.bond_table.GetJMatrix(self.row_num)[row][col].value = float(value)


class ParameterPanel(wx.Panel):
    def __init__(self, bond_table_base, *args, **kwds):
        # begin wxGlade: ParameterPanel.__init__
        kwds["style"] = wx.TAB_TRAVERSAL
        wx.Panel.__init__(self, *args, **kwds)
        self.bond_num_col_label = wx.StaticText(self, -1, "Interaction Number")
        self.tie_col_label = wx.StaticText(self, -1, "Tie Parameters")
        self.edit_col_label = wx.StaticText(self, -1, "Edit Parameters")
        #self.Type_of_Param_RadioBox = wx.RadioBox(self, -1, "Type of Parameter", choices=["Fixed Value", "Variable"], majorDimension=1, style=wx.RA_SPECIFY_ROWS)
        self.__set_properties()
        self.__do_layout()
        #self.Bind(wx.EVT_RADIOBOX, self.OnTypeChange, self.Type_of_Param_RadioBox)
        
        # end wxGlade
        self.bond_table_base = bond_table_base
        self.sized_items = []#keep a reference to all items added to the main
        #sizer so that they can be removed later.
        self.bond_table_base.AddParamPanel(self)

    def __set_properties(self):
        # begin wxGlade: ParameterPanel.__set_properties
        self.tie_col_label.SetToolTipString("To tie parameters, click a parameter, then hold Ctrl and click other paramters to tie to.")
        self.edit_col_label.SetToolTipString("Click a parameter to edit it.")
        #self.Type_of_Param_RadioBox.SetToolTipString("Is the value known(fixed), or would you like to solve for it?")
        #self.Type_of_Param_RadioBox.SetSelection(0)
   
        # end wxGlade

    def __do_layout(self):
        # begin wxGlade: ParameterPanel.__do_layout
        self.main_grid_sizer = wx.FlexGridSizer(1, 3, 10, 10)#Initially has no matrices
        sizer_1 = wx.BoxSizer(wx.VERTICAL)
        self.main_grid_sizer.Add(self.bond_num_col_label, 0, 0, 0)
        self.main_grid_sizer.Add(self.tie_col_label, 0, wx.ALIGN_CENTER_HORIZONTAL, 0)
        sizer_1.Add(self.edit_col_label, 0, wx.ALIGN_CENTER_HORIZONTAL, 0)
        #sizer_1.Add(self.Type_of_Param_RadioBox, 0, wx.ALIGN_CENTER_HORIZONTAL|wx.ALIGN_CENTER_VERTICAL|wx.FIXED_MINSIZE, 1)
        self.main_grid_sizer.Add(sizer_1, 1, wx.EXPAND, 0)
        self.main_grid_sizer.AddGrowableCol(2)
        self.SetSizer(self.main_grid_sizer)
        self.main_grid_sizer.Fit(self)
        # end wxGlade

    def AddRow(self, row_num):
        """Adds another row (to the end) with a new bond label, tie matrix,
        and parameter editing matrix corresponding to the entry in the bond
        table at row_num."""
        bond_label = wx.StaticText(self, -1, "Bond " + str(row_num+1))
        tie_grid = wx.grid.Grid(self, -1)
        edit_grid = wx.grid.Grid(self, -1)
        self.Bind(wx.grid.EVT_GRID_CMD_CELL_LEFT_CLICK, self.OnTieCellClick,
                  tie_grid)
        self.Bind(wx.grid.EVT_GRID_CMD_CELL_LEFT_CLICK, self.OnEditCellClick,
                  edit_grid)
        
        tie_grid.CreateGrid(3, 3)
        tie_grid.SetRowLabelSize(0)
        tie_grid.SetColLabelSize(0)
        #fill in the table with values
        array = self.bond_table_base.GetJMatrix(row_num)
        for i in range(3):
            for j in range(3):
                tie_grid.SetCellValue(i,j,array[i][j].getName())
        
        #self.tie_grid.SetColLabelValue(0, "")
        tie_grid.SetToolTipString("To tie parameters, click a parameter, then\
                                   hold Ctrl and click other parameters to tie\
                                   to.")
        #edit_grid.CreateGrid(3, 3)
        edit_grid.SetTable(ParamTable(self.bond_table_base, row_num))
        edit_grid.SetRowLabelSize(0)
        edit_grid.SetColLabelSize(0)
        edit_grid.SetToolTipString("Click a parameter to edit it.")
        
        self.main_grid_sizer.Add(bond_label, 0, wx.ALIGN_CENTER_HORIZONTAL |
                                 wx.ALIGN_CENTER_VERTICAL, 0)
        self.sized_items.append(bond_label)
        self.main_grid_sizer.Add(tie_grid, 1, wx.EXPAND, 0)
        self.sized_items.append(tie_grid)
        self.main_grid_sizer.Add(edit_grid, 1, wx.EXPAND, 0)
        self.sized_items.append(edit_grid)

        self.main_grid_sizer.AddGrowableRow(row_num)
        
        self.main_grid_sizer.Fit(self)
        #self.Fit()
        self.GetParent().Fit()
        self.GetParent().SetMinSize(self.GetParent().GetSize())
        
    
    def RemoveRows(self, pos, numRows):
        """remove numRows starting from the row at pos."""
        print "remove ", pos, numRows
        for i in range(pos*3, (pos+numRows)*3):#There are 3 columns
            j = (pos+numRows)*3 -i + pos*3 -1#reverse the order
            print j, len(self.sized_items)
            print self.sized_items
            item = self.sized_items.pop(j)
            self.main_grid_sizer.Remove(item)
            item.Destroy()
        self.main_grid_sizer.Layout()
        self.main_grid_sizer.Fit(self)
        self.GetParent().Fit()
        self.GetParent().SetMinSize(self.GetParent().GetSize())


    def OnTypeChange(self, event): # wxGlade: ParameterPanel.<event_handler>
        print "Event handler `OnTypeChange' not implemented!"
        event.Skip()

    def OnTieCellClick(self, event): # wxGlade: ParameterPanel.<event_handler>
        print "Event handler `OnTieCellClick' not implemented!"
        event.Skip()

    def OnEditCellClick(self, event): # wxGlade: ParameterPanel.<event_handler>
        print "Event handler `OnEditCellClick' not implemented!"
        event.Skip()

# end of class ParameterPanel


class jijDialog(wx.Dialog):
    """This dialog is displayed when the user clicks on the Jij Cell in the bond
    grid.  It allows them to enter a Jij Matrix."""
    def __init__(self, currentVal):
        wx.Dialog.__init__(self, None, -1, 'Jij Matrix', size = (300,300))
        self.matrix = currentVal
        
        okButton = wx.Button(self, wx.ID_OK, "OK", pos = (25, 225), size = (100, 25))
        okButton.SetDefault()
        cancelButton = wx.Button(self, wx.ID_CANCEL, "Cancel",  pos = (175, 225), size = (100, 25))     
#        buttonSizer = wx.BoxSizer(wx.HORIZONTAL)
#        buttonSizer.Add(okButton, 1, wx.ALIGN_CENTER_HORIZONTAL)
#        buttonSizer.Add(cancelButton, 1, wx.ALIGN_CENTER_HORIZONTAL)
#        self.SetSizer(buttonSizer)

        #Add radio buttons to switch between fixed values and variables
        self.fixedValues = True
        
        #if the currentVal contains any variable parameters, the radio button will be set
        #to 'Variable', otherwise the default will be Fixed only
        for i in range(3):
            for j in range(3):
                if currentVal[i][j].fit:
                     self.fixedValues = False
        
        self.rb = wx.RadioBox(self, -1, "", (40,30), wx.DefaultSize,['Fixed Values Only ', 'Variable Values'],
        2, wx.RA_SPECIFY_COLS)
        if not self.fixedValues:
            self.rb.SetSelection(1)
        
        self.Bind(wx.EVT_RADIOBOX, self.EvtRadioBox, self.rb)
        self.rb.SetToolTip(wx.ToolTip("Are the values known, or will they be solved for?"))

        self.grid = wx.grid.Grid(self, -1, pos = (40,80))
        self.grid.CreateGrid(3,3)
        self.grid.SetColLabelValue(0,"     a     ")
        self.grid.SetColLabelValue(1,"     b     ")
        self.grid.SetColLabelValue(2,"     c     ")
        self.grid.SetRowLabelValue(0,"a")
        self.grid.SetRowLabelValue(1,"b")
        self.grid.SetRowLabelValue(2,"c")
        #Fill the table in with the current Jmatrix value
        self.updateTable()
        
        #For validating when 'ok' button is pressed
        self.Bind(wx.EVT_BUTTON, self.OnOk, okButton)
        
        #When the user clicks on a cell
        self.Bind(wx.grid.EVT_GRID_CELL_LEFT_CLICK, self.OnLeftClick ,self.grid)
        
    def updateTable(self):
        """Updates displayed values to those in the parameter matrix"""
        self.grid.SetCellValue(0,0,str(self.matrix[0][0]))
        self.grid.SetCellValue(0,1,str(self.matrix[0][1]))
        self.grid.SetCellValue(0,2,str(self.matrix[0][2]))
        self.grid.SetCellValue(1,0,str(self.matrix[1][0]))
        self.grid.SetCellValue(1,1,str(self.matrix[1][1]))
        self.grid.SetCellValue(1,2,str(self.matrix[1][2]))
        self.grid.SetCellValue(2,0,str(self.matrix[2][0]))
        self.grid.SetCellValue(2,1,str(self.matrix[2][1]))
        self.grid.SetCellValue(2,2,str(self.matrix[2][2]))
        self.grid.AutoSize()
        
    def OnLeftClick(self, evt):
        """When there is a left lick on a cell, if the radio button is set to fixed
        values, the event will be passed on and the user will be able to enter values.
        If the radio button is set to variable, then the click will open another dialog
        in which the user can enter parameter information."""
        if self.fixedValues:
            evt.Skip()
        else:
            #Open a paramDialog
            col=evt.GetCol()
            row=evt.GetRow()
    
            dialog = ParamDialog(self.matrix[row][col], self, -1)#Pass current Jij value
            dialog.ShowModal()
            self.updateTable()
            dialog.Destroy()

#            result = dialog.ShowModal()
#            if result == wx.ID_OK:
                #since the parameter will be passed by reference, and keeping these references
                #is important for each parameter which has a list of tied parameters, the dialog
                #will set the values in the given JParam instance, rather than returning the value
#                dialog.setParam()
#                self.updateTable()
    
            #dialog.Destroy()
        
#        self.AutoSize()

    
    
    def EvtRadioBox(self, evt):
        """called when radio buttons toggling between fixed values and variable are changed."""
        self.fixedValues = (evt.GetInt() == 0)
        if not self.fixedValues:
            #for values that are correctly entered set the JParam.value to the float value
            for i in range(3):
                for j in range(3):
                    try:
                        self.matrix[i][j].value = float(self.grid.GetCellValue(i,j))
                    except:
                        """do nothing"""    
        self.updateTable()
    
    def OnOk(self, event):
        """If the OK button is pressed, and the data is all of the right type,
        this event will be passed on, which will close the window."""
        if self.validate():
            event.Skip()
    
    def validate(self):
        """If the radio button is set to fixed values:
        Makes sure all the cell values are of type float, and turns the cell
        pink if one is not."""
        bgColor = "pink"
        failed = False
        if self.fixedValues:
            for i in range(3):
                for j in range(3):
                    try:
                        float(self.grid.GetCellValue(i,j))
                        attr = wx.grid.GridCellAttr()
                        attr.SetBackgroundColour("white")
                        self.grid.SetAttr(i,j,attr)
                    except:
                        attr = wx.grid.GridCellAttr()
                        attr.SetBackgroundColour(bgColor)
                        self.grid.SetAttr(i,j, attr)
                        failed = True
            
        
        #self.grid.AutoSize() #To refresh cells
        self.grid.Refresh()
        return not failed

    def setMatrix(self):
        if self.fixedValues:
            for i in range(3):
                for j in range(3):
                    self.matrix[i][j].value = float(self.grid.GetCellValue(i,j))
                    self.matrix[i][j].fit = False


class ParamDialog(wx.Dialog):
    """This is a dialog box in which the user can enter information about a single
    parameter."""
    def __init__(self, param, *args, **kwds):
        tiedStr = str(param.tied).replace('[', '').replace(']','')
        # begin wxGlade: ParamDialog.__init__
        wx.Dialog.__init__(self, *args, **kwds)
        self.Type_of_Param_RadioBox = wx.RadioBox(self, -1, "Type of Parameter", choices=["Fixed Value", "Variable"], majorDimension=1, style=wx.RA_SPECIFY_ROWS)
        self.label_1 = wx.StaticText(self, -1, "Value:")
        self.fixed_value_TxtCtrl = wx.TextCtrl(self, -1, str(param.value))
        self.label_2 = wx.StaticText(self, -1, "Min:")
        self.min_val_TxtCtrl = wx.TextCtrl(self, -1, param.min)
        self.label_3 = wx.StaticText(self, -1, "Max:")
        self.max_val_TxtCtrl = wx.TextCtrl(self, -1, param.max)
        self.label_4 = wx.StaticText(self, -1, "Tied to Parameters:")
        self.tiedParamsTxtCtrl = wx.TextCtrl(self, -1, tiedStr)
        self.okButton_copy = wx.Button(self, wx.ID_OK, "Ok")#Changed ID for Modal
        self.cancelButton_copy = wx.Button(self, wx.ID_CANCEL, "Cancel")#Made cancel ID for modal

        self.__set_properties()
        self.__do_layout()

        self.Bind(wx.EVT_RADIOBOX, self.OnTypeChange, self.Type_of_Param_RadioBox)
        self.Bind(wx.EVT_BUTTON, self.OnOk, self.okButton_copy)
        # end wxGlade
        self.param = param
        #Default Radio control
        self.fixedVals = not self.param.fit
        self.fixed_value_TxtCtrl.Enable(self.fixedVals)
        self.min_val_TxtCtrl.Enable(not self.fixedVals)
        self.max_val_TxtCtrl.Enable(not self.fixedVals)
        self.tiedParamsTxtCtrl.Enable(not self.fixedVals)
        if self.fixedVals:
            self.Type_of_Param_RadioBox.SetSelection(0)
        else:
            self.Type_of_Param_RadioBox.SetSelection(1)

    def __set_properties(self):
        # begin wxGlade: ParamDialog.__set_properties
        self.SetTitle("Parameter")
        self.Type_of_Param_RadioBox.SetToolTipString("Is the value known(fixed), or would you like to solve for it?")
        self.Type_of_Param_RadioBox.SetSelection(0)
        self.tiedParamsTxtCtrl.SetToolTipString("write the number of the other paramters which must be equal to this parameter.")
        # end wxGlade

    def __do_layout(self):
        # begin wxGlade: ParamDialog.__do_layout
        sizer_1 = wx.BoxSizer(wx.VERTICAL)
        sizer_8_copy = wx.BoxSizer(wx.HORIZONTAL)
        sizer_7 = wx.BoxSizer(wx.VERTICAL)
        sizer_2 = wx.BoxSizer(wx.HORIZONTAL)
        sizer_4 = wx.BoxSizer(wx.HORIZONTAL)
        sizer_6 = wx.BoxSizer(wx.VERTICAL)
        sizer_5 = wx.BoxSizer(wx.VERTICAL)
        sizer_3 = wx.BoxSizer(wx.VERTICAL)
        sizer_1.Add(self.Type_of_Param_RadioBox, 0, wx.ALIGN_CENTER_HORIZONTAL|wx.ALIGN_CENTER_VERTICAL|wx.FIXED_MINSIZE, 1)
        sizer_3.Add(self.label_1, 0, 0, 0)
        sizer_3.Add(self.fixed_value_TxtCtrl, 0, 0, 0)
        sizer_2.Add(sizer_3, 1, wx.EXPAND, 0)
        sizer_5.Add(self.label_2, 0, 0, 0)
        sizer_5.Add(self.min_val_TxtCtrl, 0, 0, 0)
        sizer_4.Add(sizer_5, 1, wx.EXPAND, 0)
        sizer_6.Add(self.label_3, 0, 0, 0)
        sizer_6.Add(self.max_val_TxtCtrl, 0, 0, 0)
        sizer_4.Add(sizer_6, 1, wx.EXPAND, 0)
        sizer_2.Add(sizer_4, 1, wx.EXPAND, 0)
        sizer_1.Add(sizer_2, 1, wx.EXPAND, 0)
        sizer_7.Add(self.label_4, 0, 0, 0)
        sizer_7.Add(self.tiedParamsTxtCtrl, 0, wx.EXPAND, 0)
        sizer_1.Add(sizer_7, 1, wx.EXPAND, 0)
        sizer_1.Add((0, 15), 0, 0, 0)
        sizer_8_copy.Add(self.okButton_copy, 0, 0, 0)
        sizer_8_copy.Add(self.cancelButton_copy, 0, 0, 0)
        sizer_1.Add(sizer_8_copy, 1, wx.ALIGN_CENTER_HORIZONTAL|wx.ALIGN_CENTER_VERTICAL, 0)
        self.SetSizer(sizer_1)
        sizer_1.Fit(self)
        self.Layout()
        # end wxGlade

    def OnTypeChange(self, event): # wxGlade: ParamDialog.<event_handler>
        print "Event handler `OnTypeChange' not implemented"
        event.Skip()
        if event.GetInt()==0:
            self.fixedVals = True
            self.min_val_TxtCtrl.Enable(False)
            self.max_val_TxtCtrl.Enable(False)
            self.tiedParamsTxtCtrl.Enable(False)
            self.fixed_value_TxtCtrl.Enable(True)
        else:
            self.fixedVals = False
            self.min_val_TxtCtrl.Enable(True)
            self.max_val_TxtCtrl.Enable(True)
            self.tiedParamsTxtCtrl.Enable(True)
            self.fixed_value_TxtCtrl.Enable(False)

    def OnOk(self, event): # wxGlade: ParamDialog.<event_handler>
        valid, fixed, val, min, max, tiedParams = self.validate()
        if valid:
            self.param.fit = not fixed
            if self.param.fit:
                self.param.min = min
                self.param.max = max
                if self.param.fit:
                    self.param.tieToMany(tiedParams)
            else:
                self.param.tieToMany([])#untie all if it is set to a fixed value
                self.param.value = val

            event.Skip()#Exit
            
            #test
            for param in self.param.manager.parameters:
                print param.tied
        
    #def setParam(self):
    #    """sets the values in the JParam object to those entered in this window"""
    #    print "setParam not yet implemented"

        
    def validate(self):
        """validate all information entered in textControls to make sure it is of the 
        correct type/format"""
        success = True
        bgColor = "pink"
        val = None
        min = '-inf'
        max = '+inf'
        tiedparams = None
        tiedParamInts =[]
        if self.fixedVals:
            try:
                valStr = self.fixed_value_TxtCtrl.GetValue()
                val = float(valStr)
                self.fixed_value_TxtCtrl.SetBackgroundColour("white")
            except:
                self.fixed_value_TxtCtrl.SetBackgroundColour(bgColor)
                success = False
        else:
            try:
                minStr = self.min_val_TxtCtrl.GetValue()
                min = float(minStr)
                self.min_val_TxtCtrl.SetBackgroundColour("white")
            except:
                if minStr.strip() != "-inf":
                    self.min_val_TxtCtrl.SetBackgroundColour(bgColor)
                    success = False
            try:
                maxStr = self.max_val_TxtCtrl.GetValue()
                max = float(maxStr)
                self.max_val_TxtCtrl.SetBackgroundColour("white")
                #check that min < max
                if min != '-inf':
                    if min >= max:
                        self.min_val_TxtCtrl.SetBackgroundColour(bgColor)
                        self.max_val_TxtCtrl.SetBackgroundColour(bgColor)
                        success = False    
            except:
                if maxStr.strip() != "+inf":
                    self.max_val_TxtCtrl.SetBackgroundColour(bgColor)
                    success = False
            #validate list of tied parameters
            tiedStr = self.tiedParamsTxtCtrl.GetValue()
            tiedStr2 = tiedStr.replace(' ', '')#remove all spaces
            tiedparams = tiedStr2.split(',')
            #now resolve each element in tiedparams
            #acceptable formats include  p1,p2,... or 1,2,...
            print tiedparams
            if tiedparams != ['']:
                for param in tiedparams:
                    try:
                        if param[0] == 'p':
                            param = param[1:]#remove the 'p'
                        paramInt = int(param)
                        tiedParamInts.append(paramInt)
                        #check if a coinciding parameter to this index number actually exits
                        if not self.param.manager.validIndex(paramInt):
                            self.tiedParamsTxtCtrl.SetBackgroundColour(bgColor)
                            self.Refresh()
                            print 'invalid index'
                            return False, self.fixedVals, val, min, max, tiedParamInts
                    except:
                        self.tiedParamsTxtCtrl.SetBackgroundColour(bgColor)
                        self.Refresh()
                        return False, self.fixedVals, val, min, max, tiedParamInts
            self.tiedParamsTxtCtrl.SetBackgroundColour("white")
        self.Refresh()
        return success, self.fixedVals, val, str(min), str(max), tiedParamInts

# end of class ParamDialog



class vtkPanel(wx.Panel):
    """This is a the main panel which displays the 3D vtk rendering."""
    def __init__(self, parent, id, session):
        wx.Panel.__init__(self, parent)
        self.session = session
        self.initVTKWindow()
        self.bindEvents()
        self.mode = None
        self.picker = None
  
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

class Frame(wx.Frame):
    """This is the main frame containing the vtkPanel."""
    def __init__(self, parent, id, session):
        wx.Frame.__init__(self, parent, id, 'Magnetic Cell', size= (700,700))

        self.session = session
        self.vtkPanel = vtkPanel(self, -1, session)
        
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
            
    
    def OnLaunchSim(self, evt):
        """Runs the simulation from this app."""
        CSim.ShowSimulationFrame()
        
    def OnLaunchSpinWave(self,evt):
        frame1 = wx.Frame(self, -1, "Spinwaves")
        dlg=spinwavepanel.FormDialog(parent=frame1,id=-1)
        #dlg=spinwavepanel.FormDialog(parent=None,id=-1)
        self.SetExtraStyle(wx.WS_EX_VALIDATE_RECURSIVELY)
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
        #Maximum size is set to 25 right now, becuase that is the largest size that this computer seems to be able to reasonably handle with the current algorithm
        size = wx.GetNumberFromUser("How many times would you like to translate the cutoff cell in the a,b, and c directions?", prompt = "size:", caption = "Monte Carlo Simulation Size", value = 2, min = 1, max=25, parent = self)
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
        frame1 = wx.Frame(self.frame, -1, "Atoms", size = (500,245))
        frame1.SetMinSize((500,245))
        atomPanel(frame1, -1, session = session)
        frame1.Show()
        #Create the bond frame
        frame2 = wx.Frame(self.frame, -1, 'Bonds', size = (655,200))
        frame2.SetMinSize((655, 140))
        bondPanel(frame2, -1, session = session)
        frame2.Show()

        return True
    

def main():
    app = App(False)
    app.MainLoop()

if __name__ == '__main__':
    main()

    
    

