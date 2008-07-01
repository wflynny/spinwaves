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
from wx.py.dispatcher import connect, send
import vtkModel.SpaceGroups
import time
from Session import Session







#________________________________________________
#My Code
###########################################################



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
        
        #For now, the generate button will create a new magnetic cell
        #which will then be drawn by the vtkPanel
        self.Bind(wx.EVT_BUTTON, self.OnGenerate, self.genButton)

    
    def OnGenerate(self, event):
        failed, a, b, c, alpha, beta, gamma, magNa, magNb, magNc, cutNa, cutNb, cutNc, atomData = self.validate()
        if failed:
            return
        spaceGroup = self.spaceGroupSpinner.GetValue() #int
        send(signal = "Cell Change", sender = "Cell Panel",
             spaceGroup = spaceGroup,
             a = a, b = b, c = c,
             alpha = alpha, beta = beta, gamma = gamma,
             magNa = magNa, magNb = magNb, magNc = magNc,
             cutNa = cutNa, cutNb = cutNb, cutNc = cutNc,
             atomData = atomData)

    
    def validate(self):
        """Currently checks that all values are the right type"""
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
            self.aText.SetStyle(0, len(a), wx.TextAttr(colBack = "white"))
        except:
            self.aText.SetStyle(0, len(a), wx.TextAttr(colBack = bgColor))
            failed = True
        
        #Validate b(must be a float)
        numB = None
        try:
            numB = float(b)
            self.bText.SetStyle(0, len(b), wx.TextAttr(colBack = "white"))
        except:
            self.bText.SetStyle(0, len(b), wx.TextAttr(colBack = bgColor))
            failed = True
        
        #Validate c(must be a float)
        numC = None
        try:
            numC = float(c)
            self.cText.SetStyle(0, len(c), wx.TextAttr(colBack = "white"))
        except:
            self.cText.SetStyle(0, len(c), wx.TextAttr(colBack = bgColor))
            failed = True
        
        #Validate alpha(must be a float)
        numAlpha = None
        try:
            numAlpha = float(alpha)
            self.alphaText.SetStyle(0, len(alpha), wx.TextAttr(colBack = "white"))
        except:
            self.alphaText.SetStyle(0, len(alpha), wx.TextAttr(colBack = bgColor))
            failed = True
        
        #Validate beta(must be a float)
        numBeta = None
        try:
            numBeta = float(beta)
            self.betaText.SetStyle(0, len(beta), wx.TextAttr(colBack = "white"))
        except:
            self.betaText.SetStyle(0, len(beta), wx.TextAttr(colBack = bgColor))
            failed = True
         
        #Validate gamma(must be a float)
        numGamma = None
        try:
            numGamma = float(gamma)
            self.gammaText.SetStyle(0, len(gamma), wx.TextAttr(colBack = "white"))
        except:
            self.gammaText.SetStyle(0, len(gamma), wx.TextAttr(colBack = bgColor))
            failed = True
        
        
        #Validate Magnetic Cell Na(must be a int)
        numMagNa = None
        try:
            numMagNa = int(magNa)
            self.naText.SetStyle(0, len(magNa), wx.TextAttr(colBack = "white"))
        except:
            self.naText.SetStyle(0, len(magNa), wx.TextAttr(colBack = bgColor))
            failed = True
            
        #Validate Magnetic Cell Nb(must be a int)
        numMagNb = None
        try:
            numMagNb = int(magNb)
            self.nbText.SetStyle(0, len(magNb), wx.TextAttr(colBack = "white"))
        except:
            self.nbText.SetStyle(0, len(magNb), wx.TextAttr(colBack = bgColor))
            failed = True
            
        #Validate Magnetic Cell Nc(must be a int)
        numMagNc = None
        try:
            numMagNc = int(magNc)
            self.ncText.SetStyle(0, len(magNc), wx.TextAttr(colBack = "white"))
        except:
            self.ncText.SetStyle(0, len(magNc), wx.TextAttr(colBack = bgColor))
            failed = True
        
        #Validate cutoff Na(must be a int)
        numCutNa = None
        try:
            numCutNa = int(cutNa)
            self.cutoffNaText.SetStyle(0, len(cutNa), wx.TextAttr(colBack = "white"))
        except:
            self.cutoffNaText.SetStyle(0, len(cutNa), wx.TextAttr(colBack = bgColor))
            failed = True
        
        #Validate cutoff Nb(must be a int)
        numCutNb = None
        try:
            numCutNb = int(cutNb)
            self.cutoffNbText.SetStyle(0, len(cutNb), wx.TextAttr(colBack = "white"))
        except:
            self.cutoffNbText.SetStyle(0, len(cutNb), wx.TextAttr(colBack = bgColor))
            failed = True
        
        #Validate cutoff Nc(must be a int)
        numCutNc = None
        try:
            numCutNc = int(cutNc)
            self.cutoffNcText.SetStyle(0, len(cutNc), wx.TextAttr(colBack = "white"))
        except:
            self.cutoffNcText.SetStyle(0, len(cutNc), wx.TextAttr(colBack = bgColor))
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
                
            name = self.atomList.GetCellValue(row, 0)
            
            
            self.atomList.AutoSize()  #There may be a better way to do this, but
            #this rerenders the cells to they show the color change
            
            
            data.append([name, atomicNum, numXCoord, numYCoord, numZCoord])
        
        return failed, numA, numB, numC, numAlpha, numBeta, numGamma, numMagNa, numMagNb, numMagNc, numCutNa, numCutNb, numCutNc, data
         
    
    
    def OnGridResize(self, event):
        rows = self.atomSpinner.GetValue()
        self.atomList.SetNumberRows(rows)
#        self.atomList.GetTable().SetNumberRows(rows)
        event.Skip()


class atomListGrid(wx.grid.Grid):
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
        evt.Skip()
#        print 'LeftClick'
        col=evt.GetCol()
        row=evt.GetRow()
#        print row, col
#        if col<=0 and row >=0:
#            currval=self.table.GetValue(row,0)
        if col>=9 and row >=0: # ON/Off
            currval=self.table.GetValue(row,9)
            if currval=='':
                self.table.SetValue(row,9,'X')
            else:
                self.table.SetValue(row,9,'')
        elif col==8 and row >=0:  #Jij matrix
            dialog = jijDialog()
            result = dialog.ShowModal()
            if result == wx.ID_OK:
                print dialog.getMatrix()
                self.table.SetValue(row, 8, numpy.array(dialog.getMatrix()))

            dialog.Destroy()

        self.AutoSize()
        #if self.CanEnableCellControl():
        #    self.EnableCellEditControl()
#        wx.grid.Grid.ForceRefresh(self)

class bondPanel(wx.Panel):
    def __init__(self, parent, id, session):
        wx.Panel.__init__(self, parent, id)
        
        self.session = session
        
        self.bondList = bondListGrid(self, -1, session)
        
        self.bondSpinner = wx.SpinCtrl(self, -1, "")
        self.bondSpinner.SetRange(1,100)
        self.bondSpinner.SetValue(1)
        self.bondSpinner.Bind(wx.EVT_TEXT, self.OnGridResize, self.bondSpinner)
        
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


    def OnBondAddition(self, atom1, atom2):
        """This method is called when another window adds a bond and calls send(signal = "Bond Added"..."""
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
        
        self.OnGenerate(None)
    
    def OnGenerate(self, event):
        failed, bondData = self.validate()
        print failed
        if failed:
            return
        send(signal = "Bond Change", sender = "Bond Panel", bondData = bondData)
      
    
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
                    
#                jij = self.bondList.GetCellValue(row, 8)
#                if not jij == '':
#                    attr = wx.grid.GridCellAttr()
#                    attr.SetBackgroundColour("white")
#                    self.bondList.SetAttr(row, 8, attr)
#                else:
#                    attr = wx.grid.GridCellAttr()
#                    attr.SetBackgroundColour(bgColor)
#                    self.bondList.SetAttr(row, 8, attr)
#                    failed = True
                jij = self.bondList.GetCellValue(row, 8)  #allow a bond to be made with no Jij Matrix
                if jij == '':
                    jij = None
            
                
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
        #this rerenders the cells to they show the color change
        
        return failed, bondData
        
        
    def OnGridResize(self, event):
        rows = self.bondSpinner.GetValue()
        self.bondList.SetNumberRows(rows)
#        self.atomList.GetTable().SetNumberRows(rows)
        event.Skip()
 


class jijDialog(wx.Dialog):
    def __init__(self):
        wx.Dialog.__init__(self, None, -1, 'Jij Matrix', size = (300,300))
        okButton = wx.Button(self, wx.ID_OK, "OK", pos = (25, 225), size = (100, 25))
        okButton.SetDefault()
        cancelButton = wx.Button(self, wx.ID_CANCEL, "Cancel",  pos = (175, 225), size = (100, 25))     
#        buttonSizer = wx.BoxSizer(wx.HORIZONTAL)
#        buttonSizer.Add(okButton, 1, wx.ALIGN_CENTER_HORIZONTAL)
#        buttonSizer.Add(cancelButton, 1, wx.ALIGN_CENTER_HORIZONTAL)
#        self.SetSizer(buttonSizer)

        self.grid = wx.grid.Grid(self, -1, pos = (40,50))
        self.grid.CreateGrid(3,3)
        self.grid.SetColLabelValue(0,"    a    ")
        self.grid.SetColLabelValue(1,"    b    ")
        self.grid.SetColLabelValue(2,"    c    ")
        self.grid.SetRowLabelValue(0,"a")
        self.grid.SetRowLabelValue(1,"b")
        self.grid.SetRowLabelValue(2,"c")
        self.grid.AutoSize()
        
        #For validating when 'ok' button is pressed
        self.Bind(wx.EVT_BUTTON, self.OnOk, okButton)
    
    
    def OnOk(self, event):
        if self.validate():
            event.Skip()
    
    def validate(self):
        bgColor = "pink"
        failed = False
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
        
        self.grid.AutoSize() #To refresh cells
        return not failed
    
    def getMatrix(self):
        jMatrix = []
        for i in range(3):#rows
            row = []
            for j in range(3):
                val = float(self.grid.GetCellValue(i,j))
                row.append(val)
            jMatrix.append(row)
        
        return jMatrix
            

       

class vtkPanel(wx.Panel):
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
        connect(self.OnCellChange, signal = "Cell Change")
        connect(self.OnBondChange, signal = "Bond Change")
        connect(self.OnPick, signal = "Pick Event")
    
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
           
    def draw(self, progDialog = None):
        #remove old renderer
#        if self.ren1 != None:
 #           self.window.GetRenderWindow().RemoveRenderer(self.ren1)
        
#        a = wx.ClientDC(self.window)
        
        destroy = False
        if not progDialog:
            progDialog = wx.ProgressDialog("Progress", "Rendering Model...", parent = self, style = wx.PD_APP_MODAL | wx.PD_AUTO_HIDE)
            destroy = True
        
        #progDialog disables the render window, but we still want it to render
#       self.window.SetRenderWhenDisabled(True)
        
        # a renderer for the data
        ren1 = vtkRenderer()
        ren1.SetBackground(1,1,1)
        ren1.SetAllocatedRenderTime(3)
        
        progDialog.Update(5, "Rendering Model...")
        
        #Add the renderer to the window
        self.window.GetRenderWindow().AddRenderer(ren1)
            
        #Create vtkDrawer
        self.drawer = vtkDrawer(ren1)
           
        #Add my picker
        if self.picker:
            self.picker.removeObserver()
        self.picker = Picker(self.drawer, self.window._Iren, ren1)

        progDialog.Update(10)

        #Draw the Magnetic Cell
        self.drawer.drawMagneticCell(self.session.getMagneticCell())
        
        progDialog.Update(30)
        
        self.window.setUpRender()   
        

        self.drawer.addAxes()
        self.drawer.labelAtoms(self.session.getMagneticCell())
        
        
        progDialog.Update(100)
#        print "reset Cam"
        if destroy:
            progDialog.Destroy()

        #Rendering does not work when the window is disabled which it seems
        #to be when the progress dialog exits

#        self.window.SetRenderWhenDisabled(False)
        ren1.ResetCamera()
#        self.window.Render()
        self.window.SetRenderWhenDisabled(True)
        self.window.setUpRender()
        self.window.SetRenderWhenDisabled(False)
#        self.window._Iren.LeftButtonPressEvent() #This trigers a render, there may be better ways to do this
#        self.window._Iren.LeftButtonReleaseEvent()
        
    
    #def getStatusText(self):
      #  return self.picker.getPicked()
    
    def OnCellChange(self,spaceGroup,a,b,c,alpha, beta, gamma, magNa, magNb, magNc, cutNa, cutNb, cutNc, atomData):
        """For now this just creates a new Magnetic Cell and draws it
        This could be transfered to the Session class but a progress bar would be
        more difficult then"""
        print "Making Magentic Cell"
        
        progDialog = wx.ProgressDialog("Progress", "Generating Magnetic Cell...", parent = self, style = wx.PD_APP_MODAL | wx.PD_AUTO_HIDE)
        
        spaceGroup = SpaceGroups.GetSpaceGroup(spaceGroup)
        
        unitcell = Cell(spaceGroup, 0,0,0, a, b, c, alpha, gamma, beta)
        
        for i in range(len(atomData)):
            unitcell.generateAtoms((float(atomData[i][2]), float(atomData[i][3]), float(atomData[i][4])), atomData[i][0])

        progDialog.Update(40)
        
        print "Unit Cell Generated"
        
        #Create a Magnetic Cell
        self.session.setMagneticCell(MagneticCell(unitcell, magNa, magNb, magNc, spaceGroup))
        print "Magcell created"
        
        progDialog.Update(100)
        
        progDialog.Destroy()
        
        #Regenerate Bonds as well
        try:
            self.OnBondChange(self.bondData)
        except:#If there are not bonds yet
            print "No Bonds yet"
            self.draw()
        print "Drawn"
        

    def OnBondChange(self, bondData):
        """Checks if each bond exists, if not, it generates them""" 
        print "Creating Bonds"

        progDialog = wx.ProgressDialog("Progress", "Creating Bonds...", parent = self, style = wx.PD_APP_MODAL | wx.PD_AUTO_HIDE)
        
        self.bondData = bondData  #added this so that bonds can be regenerated later if there is a change in the cells
        self.session.getMagneticCell().clearAllBonds()
        
        progDialog.Update(5,"Creating Bonds...")
        progFrac = 99/(len(bondData) + .01)
        progress = 0
        for i in range(len(bondData)):
            cell1 = self.session.getMagneticCell().cellAtPosition((bondData[i][1], bondData[i][2], bondData[i][3]))
            cell2 = self.session.getMagneticCell().cellAtPosition((bondData[i][5], bondData[i][6], bondData[i][7]))
            
            atom1 = cell1.atomAtIndex(bondData[i][0] - 1)
            atom2 = cell2.atomAtIndex(bondData[i][4] - 1)
            
            progDialog.Update(progress)

            if not self.session.getMagneticCell().hasBond(atom1, atom2, bondData[i][8]):
                self.session.getMagneticCell().addBond(atom1, atom2, bondData[i][8])
            
            
            progress += progFrac
            print progress
            progDialog.Update(progress)

        
        progDialog.Update(100)
        progDialog.Destroy()
        print "drawing"
        self.draw()

    def OnPick(self, obj):
        print "OnPick"
        if self.mode:
            self.mode.OnPick(obj)
#        self.mode = None

            
    def OnChooseBondMode(self, event):
        self.mode = BondMode()
        
    def OnChooseNormalMode(self, event):
        self.mode = None

class BondMode():
    def __init__(self):
        self.atom1 = None
    
    def OnPick(self, obj):
        print "Pick:", obj
        if isinstance(obj, Atom):
            if self.atom1==None:
                self.atom1 = obj
                print "atom 1 =", self.atom1 
            else:
                send(signal = "Bond Added", sender = "VTK Window", atom1 = self.atom1, atom2 = obj)
#                print "OK!!!!!!!!!!!!!"
                self.atom1 = None

class Frame(wx.Frame):
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
    
 #   def OnNew(self, event):
        #Add drawing panel
#        self.vtkPanel.draw()
#        self.GetEventHandler().ProcessEvent(wx.SizeEvent())
    
    def OnCloseMe(self, event):
        self.Close(True)
    
    def OnSave(self, event):
        saveDialog = wx.FileDialog(self, "Save File", style = wx.SAVE, wildcard = "*.xml")
        if saveDialog.ShowModal() == wx.ID_OK:
            self.session.saveSessionToXML(saveDialog.GetPath())
        saveDialog.Destroy()
        
    def OnOpenFile(self, event):
        saveDialog = wx.FileDialog(self, "Open File", style = wx.OPEN, wildcard = "*.cif")
        if saveDialog.ShowModal() == wx.ID_OK:
            self.session.openCif(saveDialog.GetPath())
        saveDialog.Destroy()
    



class App(wx.App):
    def __init__(self, redirect = False, filename = None):
        wx.App.__init__(self, redirect, filename)
    
    def OnInit(self):
        session = Session()
        self.frame = Frame(None, -1, session = session)
        self.frame.Show()
        self.SetTopWindow(self.frame)
        frame1 = wx.Frame(self.frame, -1, size = (460,245))
        atomPanel(frame1, -1, session = session)
        frame1.Show()
        frame2 = wx.Frame(self.frame, -1, 'Bonds', size = (655,200))
        bondPanel(frame2, -1, session = session)
        frame2.Show()
#        dialog = jijDialog()
#        result = dialog.ShowModal()
#        if result == wx.ID_OK:
#            print "OK"
#        else:
#            print "CANCEL"
#        dialog.Destroy()
        
        
        #Williams
#       frame2 = AtomFrame(self.frame, -1)
#        frame2.Show()
        
        return True
    


if __name__ == '__main__':
    app = App(False)
    app.MainLoop()
    
    
    


