import xml.dom.minidom
import xml.dom.ext
import wx.grid
from wx.py.dispatcher import send
import CifFile
from vtkModel import SpaceGroups
from vtkModel.CellClass import Cell
from vtkModel.MagneticCellClass import MagneticCell
import numpy

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
    
    
    def openXMLSession(self, filename):
        doc = xml.dom.minidom.parse(filename)
        
        
        #Get Main node
        nodes = doc.getElementsByTagName('Spinwaves_Session')
        if len(nodes) == 1:
            mainNode = nodes[0]
        else:#There should only be one of these elements
            raise Exception('not valid XML Session file')
        
        #Get unit cell node
        nodes = mainNode.getElementsByTagName('unit_cell')
        if len(nodes) == 1:
            unitCellNode = nodes[0]
        else:#There should only be one of these elements
            raise Exception('not valid XML Session file')
        
        #Get magnetic cell node
        nodes = mainNode.getElementsByTagName('magnetic_cell')
        if len(nodes) == 1:
            magCellNode = nodes[0]
        else:#There should only be one of these elements
            raise Exception('not valid XML Session file')
        
        #Get atoms node
        nodes = mainNode.getElementsByTagName('atoms')
        if len(nodes) == 1:
            atomsNode = nodes[0]
        else:#There should only be one of these elements
            raise Exception('not valid XML Session file')
        
        #Get bonds node
        nodes = mainNode.getElementsByTagName('bonds')
        if len(nodes) == 1:
            bondsNode = nodes[0]
        else:#There should only be one of these elements
            raise Exception('not valid XML Session file')
        
        
        spaceGroupInt = int(unitCellNode.getAttribute('space_group'))
        a = unitCellNode.getAttribute('a')
        b = unitCellNode.getAttribute('b')
        c = unitCellNode.getAttribute('c')
        alpha = unitCellNode.getAttribute('alpha')
        beta = unitCellNode.getAttribute('beta')
        gamma = unitCellNode.getAttribute('gamma')
        Na = magCellNode.getAttribute('Na')
        Nb = magCellNode.getAttribute('Nb')
        Nc = magCellNode.getAttribute('Nc')
        
        #Get atoms
        atomNodes = atomsNode.getElementsByTagName('atom')
        atomData = []
        for atomNode in atomNodes:
            name = atomNode.getAttribute('Name')
            atomicNum = atomNode.getAttribute('Atomic_Number')
            x = atomNode.getAttribute('x')
            y = atomNode.getAttribute('y')
            z = atomNode.getAttribute('z')
            atomData.append([name, atomicNum, x,y,z])
        
        bondNodes = bondsNode.getElementsByTagName('bond')
        bondData = []
        for bondNode in bondNodes:
            print "here"
            atom1Num = bondNode.getAttribute('Atom1_Number')
            atom1Na = bondNode.getAttribute('Atom1_Na')
            atom1Nb = bondNode.getAttribute('Atom1_Nb')
            atom1Nc = bondNode.getAttribute('Atom1_Nc')
            atom2Num = bondNode.getAttribute('Atom2_Number')
            atom2Na = bondNode.getAttribute('Atom2_Na')
            atom2Nb = bondNode.getAttribute('Atom2_Nb')
            atom2Nc = bondNode.getAttribute('Atom2_Nc')
            on = bondNode.getAttribute('On')
            
            #Jmatrix 
            nodes = bondNode.getElementsByTagName('jMatrix')
            if len(nodes) > 0:   
                if len(nodes) == 1:
                    jNode = nodes[0]
                else:#There should only be one of these elements
                    raise Exception('not valid XML Session file')
                
                j11 = jNode.getAttribute('j11')
                j12 = jNode.getAttribute('j12')
                j13 = jNode.getAttribute('j13')
                j21 = jNode.getAttribute('j21')
                j22 = jNode.getAttribute('j22')
                j23 = jNode.getAttribute('j23')
                j31 = jNode.getAttribute('j31')
                j32 = jNode.getAttribute('j32')
                j33 = jNode.getAttribute('j33')
                
                jMatrix = numpy.array([[j11,j12,j12],
                                       [j21,j22,j23],
                                       [j31,j32,j33]])
            else:
                jMatrix = ''
            print "here2"  
            bondData.append([atom1Num, atom1Na, atom1Nb, atom1Nc, atom2Num, atom2Na, atom2Nb, atom2Nc, jMatrix, on])
            
        #Create the Magnetic Cell
        spaceGroup = SpaceGroups.GetSpaceGroup(spaceGroupInt)
        unitcell = Cell(spaceGroup, 0,0,0, a, b, c, alpha, gamma, beta)
        for i in range(len(atomData)):
            unitcell.generateAtoms((float(atomData[i][2]), float(atomData[i][3]), float(atomData[i][4])), atomData[i][0])
        self.MagCell = MagneticCell(unitcell, 1,1,1, spaceGroup)
        self.changeBonds(bondData) 
        
        #test
        print len(bondData)
        for i in range(len(bondData)):
            for j in bondData[i]:
                print j
        
        #Send Message to GUI
        send(signal = "File Load", sender = "Session", spaceGroup = spaceGroupInt, a = a, b = b, c = c, alpha = alpha, beta = beta, gamma = gamma, magNa = Na, magNb = Nb, magNc = Nc, cutNa = Na, cutNb = Nb, cutNc = Nc, atomData = atomData, bondData = bondData)
             
        
           
    def changeBonds(self, bondData):
        self.MagCell.clearAllBonds()
        for i in range(len(bondData)):
            if len(bondData[i]) < 10 or bondData[i][9] == 'X': #When this method is called by the load file method, it includes the rows that are off
                cell1 = self.MagCell.cellAtPosition((int(bondData[i][1]), int(bondData[i][2]), int(bondData[i][3])))
                cell2 = self.MagCell.cellAtPosition((int(bondData[i][5]), int(bondData[i][6]), int(bondData[i][7])))
                
                atom1 = cell1.atomAtIndex(int(bondData[i][0]) - 1)
                atom2 = cell2.atomAtIndex(int(bondData[i][4]) - 1)
    
                #bondData[i][8] can be '' instead of None
                if not isinstance(bondData[i][8], str):
                    if not self.MagCell.hasBond(atom1, atom2, bondData[i][8]):
                        self.MagCell.addBond(atom1, atom2, bondData[i][8])
                else:
                    if not self.MagCell.hasBond(atom1, atom2):
                        self.MagCell.addBond(atom1, atom2)
                    
        
        
    def openCif(self, filename):
        cf = CifFile.ReadCif(filename)
        
        #Assuming all data is in one outter block like NIST examples:
        data = cf[cf.keys()[0]]
        
        #Create a Crystollographic Unit Cell
        a = data['_cell_length_a']
        b = data['_cell_length_b']
        c = data['_cell_length_c']
        
        alpha = data['_cell_angle_alpha']
        gamma = data['_cell_angle_gamma']
        beta = data['_cell_angle_beta']
        
        spaceGroupInt = int(data['_symmetry_Int_Tables_number'])
        spaceGroup = SpaceGroups.GetSpaceGroup(spaceGroupInt)
        
        unitcell = Cell(spaceGroup, 0,0,0, a, b, c, alpha, gamma, beta)
        
        atomLabels = data['_atom_site_label']
    #Not Currently used        atomType = data['_atom_site_type_symbol']
        xPositions = data['_atom_site_fract_x']
        yPositions = data['_atom_site_fract_y']
        zPositions = data['_atom_site_fract_z']
        
        atoms = [] #for the cell window
        for i in range(len(atomLabels)):
            unitcell.generateAtoms((float(xPositions[i]), float(yPositions[i]), float(zPositions[i])), atomLabels[i])
            atoms.append([atomLabels[i], 0, float(xPositions[i]), float(yPositions[i]), float(zPositions[i])])
    
        #Create a Magnetic Cell
        self.MagCell = MagneticCell(unitcell, 1,1,1, spaceGroup)


        Na = 1  #Cif files only contain 1 unit cell
        Nb = 1
        Nc = 1
        
        #send signal to the cell window to show the info that has been loaded and to vtkWindow to draw it
        send(signal = "File Load", sender = "Session", spaceGroup = spaceGroupInt, a = a, b = b, c = c, alpha = alpha, beta = beta, gamma = gamma, magNa = Na, magNb = Nb, magNc = Nc, cutNa = Na, cutNb = Nb, cutNc = Nc, atomData = atoms, bondData = [])

        #set the values in the tables
        self.atomTable.SetValue(i, 0, atomData)

    
    def saveSessionToXML(self, filename):
        #Create the document
        doc = xml.dom.minidom.Document()
        
        #Create the main element
        main_element = doc.createElement('Spinwaves_Session')
        doc.appendChild(main_element)
        
        #Currently this will save cell info from the last tme the 
        #generate button was pressed, but atom and bond info that is
        #currently displayed in the tables
        
        #Add Unit Cell Attributes
        unitCell = doc.createElement('unit_cell')
        main_element.appendChild(unitCell)
        unitCell.setAttribute('space_group', str(self.MagCell.getAllUnitCells()[0].getSpaceGroup().number))
        unitCell.setAttribute('a', str(self.MagCell.getAllUnitCells()[0].getA()))
        unitCell.setAttribute('b', str(self.MagCell.getAllUnitCells()[0].getB()))
        unitCell.setAttribute('c', str(self.MagCell.getAllUnitCells()[0].getC()))
        unitCell.setAttribute('alpha', str(self.MagCell.getAllUnitCells()[0].getAlpha()))
        unitCell.setAttribute('beta', str(self.MagCell.getAllUnitCells()[0].getBeta()))
        unitCell.setAttribute('gamma', str(self.MagCell.getAllUnitCells()[0].getGamma()))
        
        #Add Magnetic Unit Cell Attributes
        magCell = doc.createElement('magnetic_cell')
        main_element.appendChild(magCell)
        magCell.setAttribute('Na', str(self.MagCell.getNa()))
        magCell.setAttribute('Nb', str(self.MagCell.getNb()))
        magCell.setAttribute('Nc', str(self.MagCell.getNc()))
        
        #Cutoff Cell Not yet being used
        
        
        #Add Atoms
        atomsElement = doc.createElement('atoms')
        main_element.appendChild(atomsElement)
        for i in range(self.atomTable.GetNumberRows()):
            atomElement = doc.createElement('atom')
            atomElement.setAttribute('Name', str(self.atomTable.GetValue(i, 0)))
            atomElement.setAttribute('Atomic_Number', str(self.atomTable.GetValue(i, 1)))
            atomElement.setAttribute('x', str(self.atomTable.GetValue(i, 2)))
            atomElement.setAttribute('y', str(self.atomTable.GetValue(i, 3)))
            atomElement.setAttribute('z', str(self.atomTable.GetValue(i, 4)))
            atomsElement.appendChild(atomElement)
            
          
        #Add Bonds
        bondsElement = doc.createElement('bonds')
        main_element.appendChild(bondsElement)
        for i in range(self.bondTable.GetNumberRows()):
            bondElement = doc.createElement('bond')
            bondElement.setAttribute('Atom1_Number', str(self.bondTable.GetValue(i, 0)))
            bondElement.setAttribute('Atom1_Na', str(self.bondTable.GetValue(i, 1)))
            bondElement.setAttribute('Atom1_Nb', str(self.bondTable.GetValue(i, 2)))
            bondElement.setAttribute('Atom1_Nc', str(self.bondTable.GetValue(i, 3)))
            bondElement.setAttribute('Atom2_Number', str(self.bondTable.GetValue(i, 4)))
            bondElement.setAttribute('Atom2_Na', str(self.bondTable.GetValue(i, 5)))
            bondElement.setAttribute('Atom2_Nb', str(self.bondTable.GetValue(i, 6)))
            bondElement.setAttribute('Atom2_Nc', str(self.bondTable.GetValue(i, 7)))
            bondElement.setAttribute('On', str(self.bondTable.GetValue(i, 9)))
            bondsElement.appendChild(bondElement)
          
            
            #Right now I allow bond creation without Jmatrix         
            matrix = self.bondTable.GetValue(i, 8)
            print matrix
            if matrix != '':
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
                bondElement.appendChild(jMatrix)
            
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

