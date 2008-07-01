import xml.dom.minidom
import xml.dom.ext
import wx.grid

class Session():
    """Stores information about a user session"""
    def __init__(self):
        self.bondTable = bondTable()
        self.atomTable = atomTable()
        self.MagCell = None
        
    def getAtomTable(self):
        return self.atomTable
    
    def getBondTable(self):
        return self.bondTable
    
    def getMagneticCell(self):
        return self.MagCell
    
    def setMagneticCell(self, magCell):
        #I am lettin ghe vtkPanel create magnetic cells rightt now
        #so that it can keep track of progress fo rthe progress bar.
        self.MagCell = magCell
    
    #Should add info to cell and bond windows
    def openCif(self, filename):
        self.MagCell = magneticCellFromCif(filename)
    
    def saveSessionToXML(self, filename):
        #Create the document
        doc = xml.dom.minidom.Document()
        
        #Create the main element
        main_element = doc.createElement('Spinwaves Session')
        doc.appendChild(main_element)
        
        #Currently this will save cell info from the last tme the 
        #generate button was pressed, but atom and bond info that is
        #currently displayed in the tables
        
        #Add Unit Cell Attributes
        unitCell = doc.createElement('unit cell')
        main_element.appendChild(unitCell)
        unitCell.setAttribute('a', str(self.MagCell.getAllUnitCells()[0].getA()))
        unitCell.setAttribute('b', str(self.MagCell.getAllUnitCells()[0].getB()))
        unitCell.setAttribute('c', str(self.MagCell.getAllUnitCells()[0].getC()))
        unitCell.setAttribute('alpha', str(self.MagCell.getAllUnitCells()[0].getAlpha()))
        unitCell.setAttribute('beta', str(self.MagCell.getAllUnitCells()[0].getBeta()))
        unitCell.setAttribute('gamma', str(self.MagCell.getAllUnitCells()[0].getGamma()))
        
        #Add Magnetic Unit Cell Attributes
        magCell = doc.createElement('magnetic cell')
        main_element.appendChild(magCell)
        unitCell.setAttribute('Na', str(self.MagCell.getNa()))
        unitCell.setAttribute('Nb', str(self.MagCell.getNb()))
        unitCell.setAttribute('Nc', str(self.MagCell.getNc()))
        
        #Cutoff Cell Not yet being used
        
        
        #Add Atoms
        atomsElement = doc.createElement('atoms')
        main_element.appendChild(atomsElement)
        for i in range(self.atomTable.GetNumberRows()):
            atomsElement = doc.createElement('atom')
            atomsElement.setAttribute('Name', str(self.atomTable.GetValue(i, 0)))
            atomsElement.setAttribute('Atomic Number', str(self.atomTable.GetValue(i, 1)))
            atomsElement.setAttribute('x', str(self.atomTable.GetValue(i, 2)))
            atomsElement.setAttribute('y', str(self.atomTable.GetValue(i, 3)))
            atomsElement.setAttribute('z', str(self.atomTable.GetValue(i, 4)))
          
        #Add Bonds
        bondsElement = doc.createElement('bonds')
        main_element.appendChild(bondsElement)
        for i in range(self.bondTable.GetNumberRows()):
            bondsElement = doc.createElement('bond')
            bondsElement.setAttribute('Atom1 Number', str(self.bondTable.GetValue(i, 0)))
            bondsElement.setAttribute('Atom1 Na', str(self.bondTable.GetValue(i, 1)))
            bondsElement.setAttribute('Atom1 Nb', str(self.bondTable.GetValue(i, 2)))
            bondsElement.setAttribute('Atom1 Nc', str(self.bondTable.GetValue(i, 3)))
            bondsElement.setAttribute('Atom2 Number', str(self.bondTable.GetValue(i, 4)))
            bondsElement.setAttribute('Atom2 Na', str(self.bondTable.GetValue(i, 5)))
            bondsElement.setAttribute('Atom2 Nb', str(self.bondTable.GetValue(i, 6)))
            bondsElement.setAttribute('Atom2 Nc', str(self.bondTable.GetValue(i, 7)))
            bondsElement.setAttribute('On', str(self.bondTable.GetValue(i, 9)))
          
            
            #Right now I allow bond creation without Jmatrix         
            matrix = self.bondTable.GetValue(i, 8)
            print matrix
            if matrix != None:
                jMatrix = doc.createElement('jMatrix')
                jMatrix.setAttribute('j11', str(matrix[0][0]))
                jMatrix.setAttribute('j12', str(matrix[0][1]))
                jMatrix.setAttribute('j13', str(matrix[0][2]))
                jMatrix.setAttribute('j21', str(matrix[1][0]))
                jMatrix.setAttribute('j22', str(matrix[1][1]))
                jMatrix.setAttribute('j23', str(matrix[1][2]))
                jMatrix.setAttribute('j31', str(matrix[2][0]))
                jMatrix.setAttribute('j32', str(matrix[2][1]))
                jMatrix.setAttribute('j33', str(matrix[2][2]))
                bondsElement.appendChild(jMatrix)
            
        #Write to screen
        xml.dom.ext.PrettyPrint(doc)
        
        #Write to the file
        xml.dom.ext.PrettyPrint(doc, open(filename, 'w'))
        
        
       
    
class atomTable(wx.grid.PyGridTableBase):
    def __init__(self):
        wx.grid.PyGridTableBase.__init__(self)
        self.colLabels = ['   Name   ', 'Atomic Number','       x       ', '       y       ','       z       ']
        self.rowLabels=['Atom 1']
        
        self.data = [
                     ['','','','','']#Row 1
                     ]
    
    def GetNumberRows(self):
        return len(self.data)
    def AppendRows(self, num):
        for i in range(num):
            self.AppendRow()
        return True
    def GetNumberCols(self):
        return len(self.colLabels)
    def GetColLabelValue(self, col):
        return self.colLabels[col]
    def GetRowLabelValue(self, row):
        return self.rowLabels[row]
    def IsEmptyCell(self, row, col):
        try:
            return not self.data[row][col]
        except IndexError:
            return True
    # Get/Set values in the table.  The Python version of these
    # methods can handle any data-type, (as long as the Editor and
    # Renderer understands the type too,) not just strings as in the
    # C++ version.
    def GetValue(self, row, col):
        try:
            return self.data[row][col]
        except IndexError:
            return ''
    def SetValue(self, row, col, value):
        try:
            self.data[row][col] = value
        except IndexError:
            # add a new row
            self.AppendRow()
            self.data[row][col]=value
        return
    
    def AppendRow(self):
            self.data.append([''] * self.GetNumberCols())
            self.rowLabels.append('Atom ' + str(self.GetNumberRows()))

            # tell the grid we've added a row
            msg = wx.grid.GridTableMessage(self,            # The table
                    wx.grid.GRIDTABLE_NOTIFY_ROWS_APPENDED, # what we did to it
                    1                                       # how many
                    )
            self.GetView().ProcessTableMessage(msg)
            return True

    def DeleteRows(self,pos=0,numRows=1):
        if numRows>=0 and numRows<=self.GetNumberRows():
            del self.data[pos:pos+numRows]
            msg = wx.grid.GridTableMessage(self,            # The table
            wx.grid.GRIDTABLE_NOTIFY_ROWS_DELETED, # what we did to it
            pos,numRows                                     # how many
            )
            #msg = wx.grid.GridTableMessage(self, 0, numRows)
            self.GetView().ProcessTableMessage(msg)
            
            self.UpdateValues()
            return True
        else:
            return False
        
    def UpdateValues( self ):
            """Update all displayed values"""
            msg =wx.grid.GridTableMessage(self, wx.grid.GRIDTABLE_REQUEST_VIEW_GET_VALUES)
            self.GetView().ProcessTableMessage(msg)



class bondTable(wx.grid.PyGridTableBase):
    def __init__(self):
        wx.grid.PyGridTableBase.__init__(self)
        self.colLabels = ['Atom1 Number', '   Na   ','   Nb   ', '   Nc   ', 'Atom2 Number', '   Na   ','   Nb   ', '   Nc   ', ' Jij Matrix ', 'On']
        self.rowLabels=['Bond 1']
        
        self.data = [
                     ['','','','','','','','','','']#Row 1
                     ]
    
    def GetNumberRows(self):
        return len(self.data)
    def AppendRows(self, num):
        for i in range(num):
            self.AppendRow()
        return True
    def GetNumberCols(self):
        return len(self.colLabels)
    def GetColLabelValue(self, col):
        return self.colLabels[col]
    def GetRowLabelValue(self, row):
        return self.rowLabels[row]
    def IsEmptyCell(self, row, col):
        try:
            return not self.data[row][col]
        except IndexError:
            return True
        except ValueError: #Jij Matrix
            return False
    # Get/Set values in the table.  The Python version of these
    # methods can handle any data-type, (as long as the Editor and
    # Renderer understands the type too,) not just strings as in the
    # C++ version.
    def GetValue(self, row, col):
        try:
            return self.data[row][col]
        except IndexError:
            return ''
    def SetValue(self, row, col, value):
        try:
            self.data[row][col] = value
        except IndexError:
            # add a new row
            self.AppendRow()
            self.data[row][col]=value
        return
    
 #   def AppendRows(self, numRows):
 #       for i in range(numRows):
 #           self.AppendRow()
 #       return True
    
    def AppendRow(self):
            row = [''] * (self.GetNumberCols())
            self.data.append(row)
            self.rowLabels.append('Bond ' + str(self.GetNumberRows()))

            # tell the grid we've added a row
            msg = wx.grid.GridTableMessage(self,            # The table
                    wx.grid.GRIDTABLE_NOTIFY_ROWS_APPENDED, # what we did to it
                    1                                       # how many
                    )
            self.GetView().ProcessTableMessage(msg)
            return True
    
    def DeleteRows(self,pos=0,numRows=1):
        if numRows>=0 and numRows<=self.GetNumberRows():
            del self.data[pos:pos+numRows]
            msg = wx.grid.GridTableMessage(self,            # The table
            wx.grid.GRIDTABLE_NOTIFY_ROWS_DELETED, # what we did to it
            pos,numRows                                     # how many
            )
            #msg = wx.grid.GridTableMessage(self, 0, numRows)
            self.GetView().ProcessTableMessage(msg)
            
            self.UpdateValues()
            return True
        else:
            return False
        
    def UpdateValues( self ):
            """Update all displayed values"""
            msg =wx.grid.GridTableMessage(self, wx.grid.GRIDTABLE_REQUEST_VIEW_GET_VALUES)
            self.GetView().ProcessTableMessage(msg)

