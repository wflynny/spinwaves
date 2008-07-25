import xml.dom.minidom
import xml.dom.ext
import wx.grid
from wx.py.dispatcher import send
import CifFile
from vtkModel import SpaceGroups
from vtkModel.CellClass import Cell
from vtkModel.MagneticCellClass import MagneticCell
import numpy
import time
import datetime

class Session():
    """Stores information about a user session"""
    def __init__(self):
        self.bondTable = bondTable()
        self.atomTable = atomTable()
        self.MagCell = None
        #For now the magnetic cell class will be used for the cutoff cell since
        #they are exactly the same
#        self.cutoffCell = self.MagCell
        
    def getAtomTable(self):
        return self.atomTable
    
    def getBondTable(self):
        return self.bondTable
    
    def getMagneticCell(self):
        return self.MagCell
    
    #For now magnetic cell and cuttoff cell are exaclty the same so magnetic cell will
    #be treated as cutoff cell
    def getCutoffCell(self):
        return self.MagCell
    
#    def setMagneticCell(self, magCell):
        #I am lettin ghe vtkPanel create magnetic cells right now
        #so that it can keep track of progress fo rthe progress bar.
#        self.MagCell = magCell
    
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
        Na = int(magCellNode.getAttribute('Na'))
        Nb = int(magCellNode.getAttribute('Nb'))
        Nc = int(magCellNode.getAttribute('Nc'))
        
        #Get atoms
        atomNodes = atomsNode.getElementsByTagName('atom')
        atomData = []
        for i in range(len(atomNodes)):
            name = atomNodes[i].getAttribute('Name')
            atomicNum = atomNodes[i].getAttribute('Atomic_Number')
            x = atomNodes[i].getAttribute('x')
            y = atomNodes[i].getAttribute('y')
            z = atomNodes[i].getAttribute('z')
            atomData.append([name, atomicNum, x,y,z])
            self.atomTable.SetValue(i, 0, name)
            self.atomTable.SetValue(i, 1, atomicNum)
            self.atomTable.SetValue(i, 2, x)
            self.atomTable.SetValue(i, 3, y)
            self.atomTable.SetValue(i, 4, z)
        
        bondNodes = bondsNode.getElementsByTagName('bond')
#        bondData = []
        for i in range(len(bondNodes)):
            atom1Num = bondNodes[i].getAttribute('Atom1_Number')
            atom1Na = bondNodes[i].getAttribute('Atom1_Na')
            atom1Nb = bondNodes[i].getAttribute('Atom1_Nb')
            atom1Nc = bondNodes[i].getAttribute('Atom1_Nc')
            atom2Num = bondNodes[i].getAttribute('Atom2_Number')
            atom2Na = bondNodes[i].getAttribute('Atom2_Na')
            atom2Nb = bondNodes[i].getAttribute('Atom2_Nb')
            atom2Nc = bondNodes[i].getAttribute('Atom2_Nc')
            on = bondNodes[i].getAttribute('On')
            self.bondTable.SetValue(i, 0, atom1Num)
            self.bondTable.SetValue(i, 1, atom1Na)
            self.bondTable.SetValue(i, 2, atom1Nb)
            self.bondTable.SetValue(i, 3, atom1Nc)
            self.bondTable.SetValue(i, 4, atom2Num)
            self.bondTable.SetValue(i, 5, atom2Na)
            self.bondTable.SetValue(i, 6, atom2Nb)
            self.bondTable.SetValue(i, 7, atom2Nc)
            self.bondTable.SetValue(i, 9, on)
            
            #Jmatrix 
            nodes = bondNodes[i].getElementsByTagName('jMatrix')
            if len(nodes) > 0:   
                if len(nodes) == 1:
                    jNode = nodes[0]
                else:#There should only be one of these elements
                    raise Exception('not valid XML Session file')
                
                j11 = float(jNode.getAttribute('j11'))
                j12 = float(jNode.getAttribute('j12'))
                j13 = float(jNode.getAttribute('j13'))
                j21 = float(jNode.getAttribute('j21'))
                j22 = float(jNode.getAttribute('j22'))
                j23 = float(jNode.getAttribute('j23'))
                j31 = float(jNode.getAttribute('j31'))
                j32 = float(jNode.getAttribute('j32'))
                j33 = float(jNode.getAttribute('j33'))
                
                jMatrix = numpy.array([[j11,j12,j12],
                                       [j21,j22,j23],
                                       [j31,j32,j33]])
            else:
                jMatrix = ''
            
            self.bondTable.SetValue(i, 8, jMatrix)
#            bondData.append([atom1Num, atom1Na, atom1Nb, atom1Nc, atom2Num, atom2Na, atom2Nb, atom2Nc, jMatrix, on])
            
        #Create the Magnetic Cell
#        spaceGroup = SpaceGroups.GetSpaceGroup(spaceGroupInt)
#        unitcell = Cell(spaceGroup, 0,0,0, a, b, c, alpha, gamma, beta)
#        for i in range(len(atomData)):
#            unitcell.generateAtoms((float(atomData[i][2]), float(atomData[i][3]), float(atomData[i][4])), atomData[i][0])
#        self.MagCell = MagneticCell(unitcell, 1,1,1, spaceGroup)
#        self.changeBonds(bondData) 
        
        self.cellChange(spaceGroupInt, a, b, c, alpha, beta, gamma, Na, Nb, Nc, Na, Nb, Nc, atomData)
        
        #Send Message to GUI
        send(signal = "File Load", sender = "Session", spaceGroup = spaceGroupInt, a = a, b = b, c = c, alpha = alpha, beta = beta, gamma = gamma, magNa = Na, magNb = Nb, magNc = Nc, cutNa = Na, cutNb = Nb, cutNc = Nc)
             
                
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
                        
        send(signal = "Model Change", sender = "Session")
                    
        
    def cellChange(self,spaceGroupInt,a,b,c,alpha, beta, gamma, magNa, magNb, magNc, cutNa, cutNb, cutNc, atomData):    
        spaceGroup = SpaceGroups.GetSpaceGroup(spaceGroupInt)
        
        unitcell = Cell(spaceGroup, 0,0,0, a, b, c, alpha, gamma, beta)
        
        
        for i in range(len(atomData)):
            #Since I am adding anisotropy later, some things(files) will not have it
            anisotropy = None
            try:
                anisotropy = (atomData[i][5], atomData[i][6], atomData[i][7])
            except IndexError:
                print "Anisotropy not included!"
            unitcell.generateAtoms((float(atomData[i][2]), float(atomData[i][3]), float(atomData[i][4])), atomData[i][0], anisotropy = anisotropy)
        
        #Create a Magnetic Cell
        self.MagCell = MagneticCell(unitcell, magNa, magNb, magNc, spaceGroup)
        
        #Regenerate Bonds as well
        self.changeBonds(self.bondTable.data)
    
    
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
            self.atomTable.SetValue(i, 0, atomLabels[i])
            self.atomTable.SetValue(i, 2, xPositions[i])
            self.atomTable.SetValue(i, 3, yPositions[i])
            self.atomTable.SetValue(i, 4, zPositions[i])
    
        #Create a Magnetic Cell
        self.MagCell = MagneticCell(unitcell, 1,1,1, spaceGroup)


        Na = 1  #Cif files only contain 1 unit cell
        Nb = 1
        Nc = 1
        
        self.cellChange(spaceGroupInt, a, b, c, alpha, beta, gamma, magNa = Na, magNb = Nb, magNc = Nc, cutNa = Na, cutNb = Nb, cutNc = Nc, atomData = atoms)
        #send signal to the cell window to show the info that has been loaded and to vtkWindow to draw it
        n =  self.atomTable.GetNumberRows()
        for i in range(n):
            print self.atomTable.GetValue(i, 0)
        send(signal = "File Load", sender = "Session", spaceGroup = spaceGroupInt, a = a, b = b, c = c, alpha = alpha, beta = beta, gamma = gamma, magNa = Na, magNb = Nb, magNc = Nc, cutNa = Na, cutNb = Nb, cutNc = Nc)
        
        

        #set the values in the tables
#        self.atomTable.SetValue(i, 0, atomData)


    #Need to add anisotropy to this
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
        
       
       
    #Not used by me, but may be useful for other people, should add this to GUi  
    def export(self, filename):
        size = 10
        
        initialTime = time.clock()
        
        file = open(filename, 'w')
#        file.write("#Atoms\n#Number X Y Z\n")
#        for i in range(size):
#            for j in range(size):
#                for k in range(size):
#                    for atom in self.getCutoffCell().getAllAtoms():
#                        num = atom.getIndexNumber()
#                        pos = atom.getPosition()
#                        x = pos[0] + (i * self.getCutoffCell().getNa())
#                        y = pos[1] + (j * self.getCutoffCell().getNa())
#                        z = pos[2] + (k * self.getCutoffCell().getNa())
#                        line = str(num) + " " + str(x) + " " + str(y) + " " + str(z)
#                        file.write(line + "\n")
        #Right now this will be extremely slow (8 nested for lops)       
        class SimpleBond():
            def __init__(self, pos1, pos2, jMatrix):
                self.pos1 = pos1
                self.pos2 = pos2
                self.jMatrix = jMatrix
                
            def sameBond(self, bond2):
                if self.pos1 == bond2.pos1 or self.pos1 == bond2.pos2:
                    if self.pos2 == bond2.pos2 or self.pos2 == bond2.pos1:
                        return True
                return False
        
        class SimpleBondList():
            def __init__(self):
                self.list = []
            
            def addBond(self, bond):
                if not self.containsBond(bond):
                    self.list.append(bond)
                    
            def containsBond(self, bond):
                for eachBond in self.list:
                    if eachBond.sameBond(bond):
                        return True
                return False
        
        Na = self.getCutoffCell().getNa()
        Nb = self.getCutoffCell().getNb()
        Nc = self.getCutoffCell().getNc()
        
        
        def contains(list, element):
            for item in list:
                if (item == element).all():
                    return True
            return False
        
        def indexOf(list, item):
            for i in range(len(list)):
                if item.all() == list[i].all():
                    return i
            return -1
        
        
        matrices = []
        for bond in self.getCutoffCell().getBonds():
#           pos1 = bond.getAtom1().getPosition()
#           pos2 = bond.getAtom2().getPosition()
            jMat = bond.getJMatrix()
#            count = matrices.count(jMat)
            if not contains(matrices, jMat):
                matrices.append(jMat)
        
        
        simpleCellBonds = []
        for bond in self.getCutoffCell().getBonds():
            pos1 = bond.getAtom1().getPosition()
            pos2 = bond.getAtom2().getPosition()
            jMat = bond.getJMatrix()
            newBond = SimpleBond(pos1, pos2, indexOf(matrices,jMat))
            simpleCellBonds.append(newBond)
        
        
        simpleBonds = SimpleBondList()
        for bond in simpleCellBonds:
            pos1 = bond.pos1
            pos2 = bond.pos2
            jMatInt = bond.jMatrix
            for i in range(2):
                for j in range(2):
                    for k in range(2):
                        for a in range(Na):
                            for b in range(Nb):
                                for c in range(Nc):
                                    x1 = pos1[0] + a + (Na * i)
                                    y1 = pos1[1] + b + (Nb * j)
                                    z1 = pos1[2] + c + (Nc * k)
                                    
                                    x2 = pos2[0] + a + (Na * i)
                                    y2 = pos2[1] + b + (Nb * j)
                                    z2 = pos2[2] + c + (Nc * k)    
                                    bond = SimpleBond( (x1,y1,z1), (x2,y2,z2), jMatInt )
                                    simpleBonds.addBond(bond)                              

        #Pick out bonds that link first cutoff cell and another
        interCutoffBonds = []
        for eachBond in simpleBonds.list:
            #Check if one bond is in the first cutoff cell
            if eachBond.pos1[0] < Na or eachBond.pos2[0] < Na: #x
                if eachBond.pos1[1] < Nb or eachBond.pos2[1] < Nb: #y
                    if eachBond.pos1[2] < Nc or eachBond.pos2[2] < Nc: #z
                        #check if the second bond is not in the first cutoff cell
                        if (not eachBond.pos1[0] < Na) or (not eachBond.pos2[0] < Na): #x
                            if (not eachBond.pos1[1] < Nb) or (not eachBond.pos2[1] < Nb): #y
                                if (not eachBond.pos1[2] < Nc) or (not eachBond.pos2[2] < Nc): #z
                                    interCutoffBonds.append(eachBond)
                                    
        
        
        
#        finalBondList = []
        #Translate all bonds within the cutoff cell

            
        file.write("#J Matrices\n#Number J11 J12 J13 J21 J22 J23 J31 J32 J33\n")
        for i in range(len(matrices)):
            jMat = matrices[i]
            jStr = str(i) + " " + str(jMat[0][0]) + " " + str(jMat[0][1]) + " " + str(jMat[0][2]) + " " + str(jMat[1][0]) + " " + str(jMat[1][1]) + " " + str(jMat[1][2]) + " " + str(jMat[2][0]) + " " + str(jMat[2][1]) + " " + str(jMat[2][2])
            file.write(jStr + "\n")
        
        
        file.write("#Bonds\n#X1 Y1 Z1 X2 Y2 Z2 J\n")
        for bond in simpleCellBonds:   
            pos1 = bond.pos1
            pos2 = bond.pos2
#            jStr += " " 
#            jStr += str(jMat[0][1]) 
#            jStr += " " + str(jMat[0][2]) 
#            jStr += " " + str(jMat[1][0]) + " " + str(jMat[1][1]) + " " + str(jMat[1][2]) + " " + str(jMat[2][0]) + " " + str(jMat[2][1]) + " " + str(jMat[2][2])
            for i in range(size):
                for j in range(size):
                    for k in range(size):
                        x1 = pos1[0] + (Na * i)
                        y1 = pos1[1] + (Nb * j)
                        z1 = pos1[2] + (Nc * k)
                        
                        x2 = pos2[0] + (Na * i)
                        y2 = pos2[1] + (Nb * j)
                        z2 = pos2[2] + (Nc * k)    
#                        smplBond = SimpleBond( (x1,y1,z1), (x2,y2,z2), jMat )
 #                       finalBondList.append(smplBond)
                        pos1Str = str(x1) + " " + str(y1) + " " + str(z1)
                        pos2Str = str(x2) + " " + str(y2) + " " + str(z2)
                        jStr = str(bond.jMatrix)
#There should always be a jMatrix                        if jMat != None:
                        file.write(pos1Str + " " + pos2Str + " " + jStr + "\n")
#                    else:
#                        file.write(pos1Str + " " + pos2Str + "\n")
        
        for smplBond in interCutoffBonds:
            pos1 = smplBond.pos1
            pos2 = smplBond.pos2
            jMat = smplBond.jMatrix
            jStr = str(jMat)
#            jStr = str(jMat[0][0]) + " " + str(jMat[0][1]) + " " + str(jMat[0][2]) + " " + str(jMat[1][0]) + " " + str(jMat[1][1]) + " " + str(jMat[1][2]) + " " + str(jMat[2][0]) + " " + str(jMat[2][1]) + " " + str(jMat[2][2])
            aDisp = abs(pos1[0] - pos2[0])
            bDisp = abs(pos1[1] - pos2[1])
            cDisp = abs(pos1[2] - pos2[2])
#            if pos1[0] > pos2[0]:
#                aDisp = pos1[0]
#            else:
#                aDisp = pos2[0]
#            if pos1[1] > pos2[1]:
#                bDisp = pos1[1]
#            else:
#                bDisp = pos2[1]
#            if pos1[2] > pos2[2]:
#                cDisp = pos1[2]
#            else:
#                cDisp = pos2[2]
            for i in range(size - aDisp):
                for j in range(size - bDisp):
                    for k in range(size - cDisp):
                        x1 = pos1[0] + (Na * i)
                        y1 = pos1[1] + (Nb * j)
                        z1 = pos1[2] + (Nc * k)
                        
                        x2 = pos2[0] + (Na * i)
                        y2 = pos2[1] + (Nb * j)
                        z2 = pos2[2] + (Nc * k)    
#                        smplBond = SimpleBond( (x1,y1,z1), (x2,y2,z2), jMat )
 #                       finalBondList.append(smplBond)
                        pos1Str = str(x1) + " " + str(y1) + " " + str(z1)
                        pos2Str = str(x2) + " " + str(y2) + " " + str(z2)
#There should always be a jMatrix                        if jMat != None:
                        file.write(pos1Str + " " + pos2Str + " " + jStr + "\n")
        
        #Check for reapeats in finalBond List just for testing
#        def isRepeat(finalBondList):
#            for i in range(0, len(finalBondList)):
#                for j in range(i + 1, len(finalBondList)):
#                    if finalBondList[i].sameBond(finalBondList[j]):
#                        return True
#            return False
#        
#        if isRepeat(finalBondList):
#            print "There is a repeat!"
#        else:
#            print "NO repeats!"
        
        

        
        file.close()
        seconds = time.clock() - initialTime
        minutes = int(seconds)/60
        seconds -= (minutes*60)
        hours = minutes/60
        minutes -= hours*60
        print "Done\nTime:", hours, "hours", minutes, "minutes", seconds, "seconds" 
        
       
        
    #This one orgnizes it by atom, not bond
    
    
    def exportForMonteCarlo(self, filename):
        size = 2
        
        timer = Timer()
        
        file = open(filename, 'w')
   
        class SimpleBond():
            def __init__(self, pos1, pos2, jMatrix, anisotropy1 = None, anisotropy2 = None):
                self.pos1 = pos1
                self.anisotropy1 = anisotropy1
                self.pos2 = pos2
                self.anisotropy2 = anisotropy2
                self.jMatrix = jMatrix
                
            def sameBond(self, bond2):
                if self.pos1 == bond2.pos1 or self.pos1 == bond2.pos2:
                    if self.pos2 == bond2.pos2 or self.pos2 == bond2.pos1:
                        return True
                return False
        
        Na = self.getCutoffCell().getNa()
        Nb = self.getCutoffCell().getNb()
        Nc = self.getCutoffCell().getNc()
        
        class SimpleAtom():
            def __init__(self, pos, anisotropy):
                self.anisotropy = anisotropy
                self.pos = pos
#                self.cellPos = []
#                self.cellPosX = int(pos[0])/Na
#                self.cellPosY = int(pos[1])/Nb
#                self.cellPosZ = int(pos[2])/Nc
                self.interactions = []
                #self.interactions[position of other atom] = j number
                
            #might want to change this to position later when all atoms wont be created in same list
            def addInteraction(self, atom2, jMat):
#                self.interactions[atom2] = jMat
                self.interactions.append([atom2, jMat])
            
            def __eq__(self, other):
                if self.pos[0] == other.pos[0]:
                    if self.pos[1] == other.pos[1]:
                        if self.pos[2] == other.pos[2]:
                            return True
                return False
        
        
        class SimpleBondList():
            def __init__(self):
                self.list = []
                
            def addBond(self, bond):
                if not self.containsBond(bond):
                    self.list.append(bond)
#                else:
#                    print "Duplicate Bonds!" #should not get here
                    
            def containsBond(self, bond):
                for eachBond in self.list:
                    if eachBond.sameBond(bond):
                        return True
                return False
        
        
        
        def contains(list, element):
            for item in list:
                if (item == element).all():
                    return True
            return False
        
        def atomListContains(list, element):
            for item in list:
                if item == element:
                    return True
            return False
        
        def indexOf(list, item):
            for i in range(len(list)):
                if (item == list[i]).all():
                    return i
            return -1


        
        #Create list of matrices  
        matrices = []
        for bond in self.getCutoffCell().getBonds():
#           pos1 = bond.getAtom1().getPosition()
#           pos2 = bond.getAtom2().getPosition()
            jMat = bond.getJMatrix()
#            count = matrices.count(jMat)
            if not contains(matrices, jMat):
                matrices.append(jMat)
        
        
        #create simple bonds within cutoff cell
        simpleCellBonds = []
        for bond in self.getCutoffCell().getBonds():
            pos1 = bond.getAtom1().getPosition()
            anisotropy1 = bond.getAtom1().getAnisotropy()
            print anisotropy1
            pos2 = bond.getAtom2().getPosition()
            anisotropy2 = bond.getAtom2().getAnisotropy()
            print anisotropy2
            time.sleep(.1)
            jMat = bond.getJMatrix()
            newBond = SimpleBond(pos1, pos2, indexOf(matrices,jMat), anisotropy1, anisotropy2)
            simpleCellBonds.append(newBond)
        
        
        def PosInFirstCutoff(pos):
            return (pos[0] < Na and pos[1] < Nb and pos[2] < Nc)
        
        
        cellAtoms = []
        cellBonds = SimpleBondList()
        for bond in simpleCellBonds:
            pos1 = bond.pos1
            anisotropy1 = bond.anisotropy1
            pos2 = bond.pos2
            anisotropy2 = bond.anisotropy2
            jMatInt = bond.jMatrix
            for i in range(2):
                for j in range(2):
                    for k in range(2):
                        for a in range(Na):
                            for b in range(Nb):
                                for c in range(Nc):
                                    x1 = pos1[0] + a + (Na * i)
                                    y1 = pos1[1] + b + (Nb * j)
                                    z1 = pos1[2] + c + (Nc * k)
                                    
                                    x2 = pos2[0] + a + (Na * i)
                                    y2 = pos2[1] + b + (Nb * j)
                                    z2 = pos2[2] + c + (Nc * k)  
                                    newPos1 = (x1,y1,z1)
                                    newPos2 = (x2,y2,z2)
                                    if PosInFirstCutoff(newPos1):
                                        newAtom = SimpleAtom(newPos1, anisotropy1)
                                        bond = SimpleBond( (x1,y1,z1), (x2,y2,z2), jMatInt )
                                        if not atomListContains(cellAtoms, newAtom):
                                            cellAtoms.append(newAtom)
                                        cellBonds.addBond(bond)
                                    if PosInFirstCutoff(newPos2):
                                        newAtom = SimpleAtom(newPos2, anisotropy2)
                                        bond = SimpleBond( (x1,y1,z1), (x2,y2,z2), jMatInt )
                                        if not atomListContains(cellAtoms, newAtom):
                                            cellAtoms.append(newAtom)
                                        cellBonds.addBond(bond)
                
        

         
        #Write the matrix list to the file
        file.write("#J Matrices\n#Number J11 J12 J13 J21 J22 J23 J31 J32 J33\n")
        for i in range(len(matrices)):
            jMat = matrices[i]
            jStr = str(i) + " " + str(jMat[0][0]) + " " + str(jMat[0][1]) + " " + str(jMat[0][2]) + " " + str(jMat[1][0]) + " " + str(jMat[1][1]) + " " + str(jMat[1][2]) + " " + str(jMat[2][0]) + " " + str(jMat[2][1]) + " " + str(jMat[2][2])
            file.write(jStr + "\n")
        
        
        allAtoms = []
        numAtomsPerCell = len(cellAtoms)
        
        print "atoms in cell: ", numAtomsPerCell
        
        for i in range(size):
            for j in range(size):
                for k in range(size):
                    for index in range(len(cellAtoms)):
                        pos = cellAtoms[index].pos
                        anisotropy = cellAtoms[index].anisotropy
                        x = pos[0] + (Na * i)
                        y = pos[1] + (Nb * j)
                        z = pos[2] + (Nc * k)
                        newAtom = SimpleAtom((x,y,z), anisotropy)
#                        print (len(allAtoms)%numAtomsPerCell == index)#just a check, should always be true
                        allAtoms.append(newAtom)

 
        #Add bonds cellBonds to allAtoms (currently not most efficient way, but its a short list)
        for bond in cellBonds.list:
            #Make sure each position is not yet represented the easy but inefficient way
            #It's not actually that inefficient since it should only go through a fraction of the list
            pos1 = bond.pos1
            pos2 = bond.pos2
            pos1Index = -1
            pos2Index = -1
            for i in range(len(allAtoms)):
                currentPos = allAtoms[i].pos
                if currentPos == pos1:
                    pos1Index = i
                    break
                
            for i in range(len(allAtoms)):
                currentPos = allAtoms[i].pos  
                if currentPos == pos2:
                    pos2Index = i
                    break
            
            if pos1Index < 0 or pos2Index < 0:
                print "Atom list does not contain all atoms!" 
            else:
                allAtoms[pos1Index].addInteraction(pos2Index, bond.jMatrix)
                allAtoms[pos2Index].addInteraction(pos1Index, bond.jMatrix)


        timer.printTime()
        print"translating..."
        
        #translate bonds
        for i in range(len(allAtoms)- numAtomsPerCell):
            for interaction in allAtoms[i].interactions:
                newInteraction = interaction[0] + numAtomsPerCell
                if newInteraction < len(allAtoms):
                    allAtoms[i+numAtomsPerCell].addInteraction(newInteraction, interaction[1])
        
        
        
#        print "done translating, checking list"
#        timer.printTime()
        
        #Check for reapeats in finalBond List just for testing
#        def isRepeat(finalBondList):
#            for i in range(0, len(finalBondList)):
#                for j in range(i + 1, len(finalBondList)):
#                    if finalBondList[i].sameBond(finalBondList[j]):
#                        return True
#            return False
#        
#        if isRepeat(finalBondList):
#            print "There is a repeat!"
#        else:
#            print "NO repeats!"
#            
#        timer.printTime()
        
        
        #Check the simple atom list
        def atomBalanced(atomIndex):
            atom = allAtoms[atomIndex]
            for otherAtomIndex in range(len(allAtoms)):
                otherAtom = allAtoms[otherAtomIndex]
                if atomInteractsWithAtom(atomIndex, otherAtom):
                    if atomInteractsWithAtom(otherAtomIndex, atom):
                        return True
                    else:
                        return False
            return False
        
        def atomInteractsWithAtom(atomIndex, otherAtom):
            for interaction in otherAtom.interactions:
                if atomIndex == interaction[0]:
                    return True
            return False
        
        
#        for atomIndex in range(len(allAtoms)):
#            if not atomBalanced(atomIndex):
#                print "Not Balanced!!!"
#                break
#        else:
#            print "Balanced!"
            
        
        timer.printTime()
        print "number of atoms: ", len(allAtoms), "\n writing to disk..."
            
        
        #print out the simple atom list
        file.write("#AtomNumber AtomPosition(X Y Z) Anisotropy(X Y Z) OtherIndex Jmatrix OtherIndex Jmatrix...\n")
        for atomIndex in range(len(allAtoms)):
            atom = allAtoms[atomIndex]
            atomStr = str(atomIndex) + " " + str(atom.pos[0]) + " " + str(atom.pos[1]) + " " + str(atom.pos[2])
            atomStr += " " + str(atom.anisotropy[0]) + " " + str(atom.anisotropy[1]) + " " + str(atom.anisotropy[2])
            for interaction in atom.interactions:
                otherAtom = interaction[0]
                jMat = interaction[1]
                atomStr += " " + str(otherAtom)
                atomStr += " " + str(jMat)
            file.write(atomStr + "\n")
        
        file.close()
        print "done"
        timer.printTime()
            
        
    def loadSpinFile(self, fileName):
        def getSpin(position, lines):
            x = str(position[0])
            y = str(position[1])
            z = str(position[2])
            
            spin = None
            
            for line in lines:
                if not line.startswith('#'):
                    values = line.split()
                    if x == values[1] and y == values[2] and z == values[3]:
                        spin = (float(values[4]), float(values[5]), float(values[6]))
                        break
            
            return spin
        
        
        file = open(fileName, 'r')
        lines = file.readlines()
        file.close()
        
        for atom in self.getCutoffCell().getAllAtoms():
            newSpin = getSpin(atom.getPosition(), lines)
            if newSpin != None:
                atom.setSpin(newSpin)
                
        #for test purposes
        print "test:"
        for atom in self.getCutoffCell().getAllAtoms():
            print atom.getSpin()
        
        send(signal = "Model Change", sender = "Session")
        
        
                
            
        
        



class Timer():
    """for diagnostics"""
    def __init__(self):
        self.initialTime = time.clock()
    
    def printTime(self):
        seconds = time.clock() - self.initialTime
        minutes = int(seconds)/60
        seconds -= (minutes*60)
        hours = minutes/60
        minutes -= hours*60
        print "Time:", hours, "hours", minutes, "minutes", seconds, "seconds" 
 
 
    
class atomTable(wx.grid.PyGridTableBase):
    def __init__(self):
        wx.grid.PyGridTableBase.__init__(self)
        self.colLabels = ['   Name   ', 'Atomic Number','    x    ', '    y    ','    z    ', '  Dx  ', '  Dy  ', '  Dz  ']
        self.rowLabels=['Atom 1']
        
        self.data = [
                     ['','','','','', 0.0, 0.0, 0.0]#Row 1
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
            self.SetValue(self.GetNumberRows() - 1, 5, 0.0)
            self.SetValue(self.GetNumberRows() - 1, 6, 0.0)
            self.SetValue(self.GetNumberRows() - 1, 7, 0.0)

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

