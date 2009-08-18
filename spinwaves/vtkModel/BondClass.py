import AtomClass
#from vtk import *
import numpy
import wx
from wx.py.dispatcher import connect

class JParam():
    """This class represents one value in a 3*3 J matrix.  In the simplest case it would
    contain a single float value, but it can also be a variable tied to other variables
    for fitting purposes."""
    def __init__(self, manager, fit = False, value = 0., min = '-inf', max = '+inf', default = 0.):
        """-fit is a boolean value specifying whether this parameter is variable(True)
        or a fixed float(False)
        -manager is the instance of ParamManager that is to contain this JParam.
        This JParam will be appended to the manager list.
        -value is the float value if there is one(if fit is False)
        -min is a string of the float number representing the minimum possible value
        (only used if fit is TRUE) '-inf' will be used for negative infinite.
        -max is a string of the float number representing the maximim possible value
        (only used if fit is TRUE) '+inf' will be used for positive infinite.
        -default is the value the monte carlo ground state calculator file 
        export will default to if the value is set to variable."""
        self.manager = manager
        self.manager.addParam(self)
        
        
        self.fit = fit
        self.value = value
        self.min = min
        self.max = max
        self.default = default
        self.tied = []#only the manager should touch this
        self.group = None
        self.manager.AssignNewGroup(self)
        
    def tieTo(self, index):
        """Ties this parameter to the JParam in the manager at the given index."""
        self.manager.tie(self, index)
        
    def untie(self, index):
        """Will remove the tie between this parameter and the one given by index."""
        self.manager.untie(self, index)
    
    
    def tieToMany(self, list):
        """will tie this parameter to all the parameters in the given by the list of indices.
        Any ties that currently exist that are not in the list will be removed."""
        for p in self.tied:#untie all
            self.untie(p)
        for p in list:#tie to all in the new list
            self.tieTo(p)
    
    def isDefault(self):
        """Returns true if the values are the in this parameter are the default values, or
        if False if they have been set to something else."""
        if not self.fit:
            if self.value == 0.0:
                return True
        return False
    
    def getName(self):
        """returns the name of this parameter in the form p#, where # is the index number
        in the manager's list of parameters."""
        return "p" + str(self.manager.getIndex(self))
    
    def isTiedTo(self, index):
        """returns True if this parameter is tied to the parameter given by index."""
        try:#check if the list already contains the index
            self.tied.index(index)
            return True
        except:
            return False
    
    def GetIndex(self):
        return self.manager.getIndex(self)
    
    def __str__(self):
        """This will return the value if fit==False, or the range of values if fit==True"""
        if not self.fit:
            return str(self.value)
        return self.min + " - " + self.max



class Bond():
    """This class represents interactions, or "bonds" between atoms."""
    
    def __init__(self, Atom1, Atom2, jMatrix = None, r = 0,g = 0,b = 1):
        """Cell is the Unit Cell if this bond is in a unit Cell or None Otherwise.
        Atom1 and Atom2 are the atoms that this bond connects.
        r,g,b are the color of the actor. jMatrix is a numpy 2D array, such as
        [[1, 0, 0],
         [0, 1, 0],
         [0, 0, 1]]   for a ferromagnetic interaction.  The values in the jMatrix
         will be JParam objects, so that the values can be variable."""

        #Color
        self.r = r
        self.g = g
        self.b = b
        
        self.Atom1 = Atom1
        self.Atom2 = Atom2
        
        self.jMat = jMatrix #a numpy 2D array
        
        #New attributes added for fitting purposes

    
    def getAtom1(self):
        return self.Atom1
    
    def getAtom2(self):
        return self.Atom2
    
    def setAtom1(self, atom):
        self.Atom1 = atom
    
    def setAtom2(self, atom):
        self.Atom2 = atom

    def sameBond(self, otherBond):
        """returns true if otherBond connects the same 2 atoms and false otherwise"""
        if self.getAtom1() == otherBond.getAtom1() or self.getAtom1() == otherBond.getAtom2():
            if self.getAtom2() == otherBond.getAtom1() or self.getAtom2() == otherBond.getAtom2():
                return True
        return False
    
    def __str__(self):
        str = "Bond between " + self.Atom1.__str__() + " and " + self.Atom2.__str__()
        if self.jMat != None:
            str += ":\n" + self.jMat.__str__()
        return str
    
    def getRGBColor(self):
        return self.r, self.g, self.b
    
    def getJMatrix(self):
        return self.jMat

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
        wx.Panel.__init__(self, parent, id)
        
        self.session = session
        #Create the table of bonds
        self.bondList = bondListGrid(self, -1, session)
        
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
        
        #Fit this window
        self.Fit()
        self.GetParent().Fit()#Fit the frame containing this panel
        self.GetParent().SetMinSize(self.GetParent().GetSize())
        
        connect(self.OnParamChange, signal = "Parameter Values Changed")


    def OnParamChange(self):
        self.bondList.Refresh()
    
    def OnFileLoad(self):
        """Executed when the session sends a message that a file was loaded."""
        self.bondSpinner.SetValue(self.bondList.GetNumberRows())
        self.bondList.AutoSize()
        
        self.Fit()
        self.GetParent().Fit()#Fit the frame containing this panel
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
        self.Fit()
        self.GetParent().Fit()#Fit the frame containing this panel
        
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
        self.Fit()
        self.GetParent().Fit()#Fit the frame containing this panel
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
        try:#sometimes it is not Destroyed fast enough
            return self.bond_table.GetJMatrix(self.row_num)[row][col].__str__()
        except:
            return ''
        
    def SetValue(self, row, col, value):
        """Attempts to change the JParam in the bond table base."""
        self.bond_table.GetJMatrix(self.row_num)[row][col].value = float(value)
        
    def GetParameter(self, row, col):
        return self.bond_table.GetJMatrix(self.row_num)[row][col]


class ParameterPanel(wx.Panel):
    """This is the panel in which the user can edit parameters(values in the Jij
    martices)."""
    def __init__(self, bond_table_base, *args, **kwds):
        kwds["style"] = wx.TAB_TRAVERSAL
        wx.Panel.__init__(self, *args, **kwds)
        self.bond_num_col_label = wx.StaticText(self, -1, "Interaction Number")
        self.tie_col_label = wx.StaticText(self, -1, "Tie Parameters")
        self.edit_col_label = wx.StaticText(self, -1, "Edit Parameters")
        #self.Type_of_Param_RadioBox = wx.RadioBox(self, -1, "Type of Parameter", choices=["Fixed Value", "Variable"], majorDimension=1, style=wx.RA_SPECIFY_ROWS)
        
        self.main_grid_sizer = wx.FlexGridSizer(1, 3, 10, 10)#Initially has no matrices
        self.__set_properties()
        self.__do_layout()
        #self.Bind(wx.EVT_RADIOBOX, self.OnTypeChange, self.Type_of_Param_RadioBox)
        
        self.bond_table_base = bond_table_base
        #self.sized_items = []#keep a reference to all items added to the main
        #sizer so that they can be removed later.
        self.labels = []
        self.tie_grids = []
        self.edit_grids = []
        self.selected_parameter = None #for click tying
        self.bond_table_base.AddParamPanel(self)
        
        self.colors = []
        self.currentColor = 0
        self.colorMappings = []
        self.__populateColorList()
        
        connect(self.OnParamChange, signal = "Parameter Values Changed")
        connect(self.UpdateTables, signal = "File Load")
    
    def OnParamChange(self):
        for edit_grid in self.edit_grids:
            edit_grid.Refresh()
        
    def __populateColorList(self):
        #Create the list of colors used for visualization of the ties between parameters
        for r in range(0,256,64):
            for g in range(0,256,64):
                for b in range(0,256,64):
                    if(r+g+b > 128):#get rid of dark colors
                        self.colors.append((r,g,b))
        self.colors.reverse()#lighter colors first
        
        #remove similar colors
        i = 0
        while i < len(self.colors):
            c1 = self.colors[i]
            i2 = i+1
            while i2 < len(self.colors):
                c2 = self.colors[i2]
                if c2[0] + c2[1] + c2[2] > 192:
                    if c2[0] == c1[0] and c2[1] == c1[1] and abs(c2[2] - c1[2]) < 85:
                        if c2[0]>128 or c2[1]>128:
                            self.colors.pop(i2)
                    elif c2[0] == c1[0] and c2[2] == c1[2] and abs(c2[1] - c1[1]) < 85:
                        if c2[2]>128 or c2[0]>128:
                            self.colors.pop(i2)
                    elif c2[2] == c1[2] and c2[1] == c1[1] and abs(c2[0] - c1[0]) < 85:
                        if c2[2]>128 or c2[1]>128:
                            self.colors.pop(i2)
                i2 += 1
            i += 1
        
        self.colors.pop(12)
        #self.colors.pop()
        
        #shuffle them
        for i in range(0,len(self.colors),2):
            self.colors.append(self.colors.pop(i))
        #print self.colors
        
    def _getColor(self, group):
        for mapping in self.colorMappings:
            if mapping[0] == group:
                return self.colors[mapping[1]]
        
        if(self.currentColor+1 < len(self.colors)):
            self.currentColor +=1
        else:
            self.currentColor = 0
            dialog = wx.MessageDialog(self, "There are no more unique colors, colors must now be repeated.", "Out of Colors")
            dialog.ShowModal()
        self.colorMappings.append((group, self.currentColor))
        return self.colors[self.currentColor]

    def __set_properties(self):
        # begin wxGlade: ParameterPanel.__set_properties
        self.tie_col_label.SetToolTipString("To tie parameters, click a parameter, then hold Ctrl and click other parameters to tie to.")
        self.edit_col_label.SetToolTipString("Click a parameter to edit it.")
        #self.Type_of_Param_RadioBox.SetToolTipString("Is the value known(fixed), or would you like to solve for it?")
        #self.Type_of_Param_RadioBox.SetSelection(0)
   
        # end wxGlade

    def __do_layout(self):
        # begin wxGlade: ParameterPanel.__do_layout
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
        self.Bind(wx.grid.EVT_GRID_CMD_CELL_LEFT_DCLICK, self.OnEditCellDClick, 
                  tie_grid)
        
        tie_grid.CreateGrid(3, 3)
        tie_grid.SetRowLabelSize(0)
        tie_grid.SetColLabelSize(0)
        #fill in the table with values
        array = self.bond_table_base.GetJMatrix(row_num)
        for i in range(3):
            for j in range(3):
                tie_grid.SetCellValue(i,j,array[i][j].getName())
        
        #self.tie_grid.SetColLabelValue(0, "")
        tie_grid.SetToolTipString("To tie parameters, click a parameter, then \
hold Ctrl and click other parameters to tie to.")
        
        #edit_grid.CreateGrid(3, 3)
        edit_grid.SetTable(ParamTable(self.bond_table_base, row_num))
        edit_grid.SetRowLabelSize(0)
        edit_grid.SetColLabelSize(0)
        edit_grid.SetToolTipString("Click a parameter to edit it.")
        
        self.main_grid_sizer.Add(bond_label, 0, wx.ALIGN_CENTER_HORIZONTAL |
                                 wx.ALIGN_CENTER_VERTICAL, 0)
        self.labels.append(bond_label)
        self.main_grid_sizer.Add(tie_grid, 1, wx.EXPAND, 0)
        self.tie_grids.append(tie_grid)
        self.main_grid_sizer.Add(edit_grid, 1, wx.EXPAND, 0)
        self.edit_grids.append(edit_grid)

        self.main_grid_sizer.AddGrowableRow(row_num+1)
        
        self.main_grid_sizer.Layout()
        self.main_grid_sizer.Fit(self)
        self.GetParent().Fit()
        self.GetParent().SetMinSize(self.GetParent().GetSize())
        
    
    def RemoveRows(self, pos, numRows):
        """remove numRows starting from the row at pos."""
        #print "remove ", pos, numRows
        for i in range(pos, (pos+numRows)):#There are 3 columns
            j = (pos+numRows) -i + pos -1#reverse the order
            #print j, len(self.labels), len(self.tie_grids), len(self.edit_grids)
            label = self.labels.pop(j)
            tie_grid = self.tie_grids.pop(j)
            edit_grid = self.edit_grids.pop(j)
            self.main_grid_sizer.Remove(label)
            self.main_grid_sizer.Remove(tie_grid)
            self.main_grid_sizer.Remove(edit_grid)
            label.Destroy()
            tie_grid.Destroy()
            edit_grid.Destroy()
            #decrease the row num in the existing objects
            for index in range(j, len(self.labels)):
                self.main_grid_sizer.Remove(self.labels[index])
                self.labels[index] = wx.StaticText(self, -1, "Bond " + str(index+1))
                self.main_grid_sizer.Insert((index+1)*3, self.labels[index], 0,
                                            wx.ALIGN_CENTER_HORIZONTAL |
                                            wx.ALIGN_CENTER_VERTICAL, 0)
                self.edit_grids[index].GetTable().row_num -= 1
                array = self.bond_table_base.GetJMatrix(index)
                for row in range(3):
                    for col in range(3):
                        self.tie_grids[index].SetCellValue(row, col, array[row][col].getName())
        self.main_grid_sizer.Layout()
        self.main_grid_sizer.Fit(self)
        #Increasing the min size fits the frame nicely
        self.GetParent().SetMinSize((1,1))
        self.GetParent().Fit()
        self.GetParent().SetMinSize(self.GetParent().GetSize())


    def OnTypeChange(self, event): # wxGlade: ParameterPanel.<event_handler>
        print "Event handler `OnTypeChange' not implemented!"
        event.Skip()

    def OnTieCellClick(self, event):
        """When a cell in the parameter tying matrices is clicked and the
        Control Key is held, it will be tied to the last parameter to be
        clicked.  Otherwise, if the Control Key is not held, the parameter will
        be tied to any parameters which are clicked later while holding the
        Control Key."""
        
        col=event.GetCol()
        row=event.GetRow()
        grid = event.GetEventObject()
        #Get the parameter object from the associated edit_grid
        edit_grid = self.edit_grids[self.tie_grids.index(grid)]
        param = edit_grid.GetTable().GetParameter(row, col)
        if not wx.GetMouseState().ControlDown():
            self.selected_parameter = param
        else:
            self.selected_parameter.tieTo(param.GetIndex())
        self.UpdateTables()
        send(signal = "Parameter Values Changed", sender = "Jij Dialog")


    def OnEditCellClick(self, event):

        #If the radio button is later added, the event will be skipped if it is
        #set to fixed values
        #event.Skip()
        #Open a paramDialog
        col=event.GetCol()
        row=event.GetRow()
        grid = event.GetEventObject()

        dialog = ParamDialog(grid.GetTable().GetParameter(row, col), self, -1)#Pass current Jij value
        dialog.ShowModal()
        #grid.Refresh()
        self.UpdateTables()
        dialog.Destroy()
        
    def OnEditCellDClick(self, evt):
        """Untie the parameter when it is double clicked."""
        index = self.tie_grids.index(evt.GetEventObject())
        param = self.edit_grids[index].GetTable().GetParameter(evt.GetRow(), 
                                                               evt.GetCol())
        param.tieToMany([])
        self.UpdateTables()
        
    def UpdateTables(self):
        """Updates the grid displays to match the values of the parameters they
        represent."""
        for i in range(len(self.labels)):
            #Set the colors of the cells in tie grids
            for row in range(3):
                for col in range(3):
                    #Get the parameter object from the associated edit_grid
                    edit_grid = self.edit_grids[i]
                    param = edit_grid.GetTable().GetParameter(row, col)
                    if len(param.tied) == 0:
                        color = (255, 255, 255)#white if not tied
                    else:
                        color = self._getColor(param.group)
                    attr = wx.grid.GridCellAttr()
                    attr.SetBackgroundColour(color)
                    self.tie_grids[i].SetAttr(row,col,attr)
                    #print param.group, color
            self.tie_grids[i].Refresh()
            self.tie_grids[i].ClearSelection()
            self.edit_grids[i].Refresh()

# end of class ParameterPanel


class jijDialog(wx.Dialog):
    """This dialog is displayed when the user clicks on the Jij Cell in the bond
    grid.  It allows them to enter a Jij Matrix."""
    def __init__(self, currentVal):
        wx.Dialog.__init__(self, None, -1, 'Jij Matrix')
        self.matrix = currentVal
       
        mainSizer = wx.BoxSizer(wx.VERTICAL)   

        #Add radio buttons to switch between fixed values and variables
        self.fixedValues = True
        
        #if the currentVal contains any variable parameters, the radio button will be set
        #to 'Variable', otherwise the default will be Fixed only
        for i in range(3):
            for j in range(3):
                if currentVal[i][j].fit:
                    self.fixedValues = False
        
        self.rb = wx.RadioBox(self, -1, "", wx.DefaultPosition, wx.DefaultSize,['Fixed Values Only ', 'Variable Values'],
        2, wx.RA_SPECIFY_COLS)
        mainSizer.Add(self.rb, 0,wx.ALIGN_CENTER_HORIZONTAL | wx.BOTTOM ,10)
        if not self.fixedValues:
            self.rb.SetSelection(1)
        
        self.Bind(wx.EVT_RADIOBOX, self.EvtRadioBox, self.rb)
        self.rb.SetToolTip(wx.ToolTip("Are the values known, or will they be solved for?"))

        self.grid = wx.grid.Grid(self, -1)
        self.grid.CreateGrid(3,3)
        self.grid.SetColLabelValue(0,"     a     ")
        self.grid.SetColLabelValue(1,"     b     ")
        self.grid.SetColLabelValue(2,"     c     ")
        self.grid.SetRowLabelValue(0,"a")
        self.grid.SetRowLabelValue(1,"b")
        self.grid.SetRowLabelValue(2,"c")
        #Fill the table in with the current Jmatrix value
        self.updateTable()
        mainSizer.Add(self.grid, 0, wx.ALIGN_CENTER_HORIZONTAL, 0)
        
        okButton = wx.Button(self, wx.ID_OK, "OK")
        okButton.SetDefault()
        cancelButton = wx.Button(self, wx.ID_CANCEL, "Cancel")
        btnSizer = wx.BoxSizer(wx.HORIZONTAL)
        btnSizer.Add(okButton)
        btnSizer.Add(cancelButton)
        mainSizer.Add((1,15),0,0,0)
        mainSizer.Add(btnSizer, 0, wx.ALIGN_CENTER_HORIZONTAL, 0)
        #self.SetSizer(mainSizer)
        border = wx.BoxSizer(wx.HORIZONTAL)
        border.Add(mainSizer, 0, wx.ALL, 20)
        self.SetSizer(border)
        self.Fit()
        
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
                        self.matrix[i][j].default = self.matrix[i][j].value
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
        """If the radio button is set to fixed values, then the values of the
        parameters will be set to the values entered in the grid.  The numbers
        should be validated first before this function is called."""
        if self.fixedValues:
            for i in range(3):
                for j in range(3):
                    self.matrix[i][j].value = float(self.grid.GetCellValue(i,j))
                    self.matrix[i][j].default = self.matrix[i][j].value
                    self.matrix[i][j].fit = False
                    send(signal = "Parameter Values Changed", sender = "Jij Dialog")


class ParamDialog(wx.Dialog):
    """This is a dialog box in which the user can enter information about a single
    parameter."""
    def __init__(self, param, *args, **kwds):
        tiedStr = str(param.tied).replace('[', '').replace(']','')
        # begin wxGlade: ParamDialog.__init__
        wx.Dialog.__init__(self, *args, **kwds)
        self.Param_Type_RB = wx.RadioBox(self, -1, "Type of Parameter",
                                         choices=["Fixed Value", "Variable"],
                                         majorDimension=1,
                                         style=wx.RA_SPECIFY_ROWS)
        self.label_1 = wx.StaticText(self, -1, "Value:")
        self.fixed_value_TxtCtrl = wx.TextCtrl(self, -1, str(param.value))
        self.label_2 = wx.StaticText(self, -1, "Min:")
        self.min_val_TxtCtrl = wx.TextCtrl(self, -1, param.min)
        self.label_3 = wx.StaticText(self, -1, "Max:")
        self.max_val_TxtCtrl = wx.TextCtrl(self, -1, param.max)
        self.label_4 = wx.StaticText(self, -1, "Tied to Parameters:")
        self.tiedParamsTxtCtrl = wx.TextCtrl(self, -1, tiedStr)
        self.okButton_copy = wx.Button(self, wx.ID_OK, "Ok")
        self.cancelButton_copy = wx.Button(self, wx.ID_CANCEL, "Cancel")
        self.default_label = wx.StaticText(self, -1, "Default")
        self.defaultTxtCtrl = wx.TextCtrl(self, -1, str(param.default))

        self.__set_properties()
        self.__do_layout()

        self.Bind(wx.EVT_RADIOBOX, self.OnTypeChange, self.Param_Type_RB)
        self.Bind(wx.EVT_BUTTON, self.OnOk, self.okButton_copy)
        # end wxGlade
        self.param = param
        #Default Radio control
        self.fixedVals = not self.param.fit
        self.fixed_value_TxtCtrl.Enable(self.fixedVals)
        self.min_val_TxtCtrl.Enable(not self.fixedVals)
        self.max_val_TxtCtrl.Enable(not self.fixedVals)
        self.tiedParamsTxtCtrl.Enable(not self.fixedVals)
        self.defaultTxtCtrl.Enable(not self.fixed_value_TxtCtrl)
        if self.fixedVals:
            self.Param_Type_RB.SetSelection(0)
        else:
            self.Param_Type_RB.SetSelection(1)

    def __set_properties(self):
        # begin wxGlade: ParamDialog.__set_properties
        self.SetTitle("Parameter")
        self.Param_Type_RB.SetToolTipString("Is the value known(fixed), or would you like to solve for it?")
        self.Param_Type_RB.SetSelection(0)
        self.tiedParamsTxtCtrl.SetToolTipString("write the number of the other parameters which must be equal to this parameter.")
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
        sizer_1.Add(self.Param_Type_RB, 0, wx.ALIGN_CENTER_HORIZONTAL|wx.ALIGN_CENTER_VERTICAL|wx.FIXED_MINSIZE, 1)
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
        
        sizer_1.Add(self.default_label, 0, 0, 0)
        sizer_1.Add(self.defaultTxtCtrl, 0, 0, 0)
        
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
            self.defaultTxtCtrl.Enable(False)
        else:
            self.fixedVals = False
            self.min_val_TxtCtrl.Enable(True)
            self.max_val_TxtCtrl.Enable(True)
            self.tiedParamsTxtCtrl.Enable(True)
            self.fixed_value_TxtCtrl.Enable(False)
            self.defaultTxtCtrl.Enable(True)

    def OnOk(self, event): # wxGlade: ParamDialog.<event_handler>
        valid, fixed, val, min, max, tiedParams, default = self.validate()
        if valid:
            self.param.fit = not fixed
            if self.param.fit:
                self.param.min = min
                self.param.max = max
                self.param.default = default
                if self.param.fit:
                    self.param.tieToMany(tiedParams)
            else:
                self.param.tieToMany([])#untie all if it is set to a fixed value
                self.param.value = val
                self.param.default = val

            event.Skip()#Exit
            
            #test
            for param in self.param.manager.parameters:
                print param.tied
        
        send(signal = "Parameter Values Changed", sender = "Parameter Dialog")
        
    #def setParam(self):
    #    """sets the values in the JParam object to those entered in this window"""
    #    print "setParam not yet implemented"

        
    def validate(self):
        """validate all information entered in textControls to make sure it is of the 
        correct type/format"""
        success = True
        bgColor = "pink"
        val = None
        default = None
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
            try:
                defaultStr = self.defaultTxtCtrl.GetValue()
                default = float(defaultStr)
                self.defaultTxtCtrl.SetBackgroundColour("white")
            except:
                self.defaultTxtCtrl.SetBackgroundColour(bgColor)
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
        return success, self.fixedVals, val, str(min), str(max), tiedParamInts, default

# end of class ParamDialog
