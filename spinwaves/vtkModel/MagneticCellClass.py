from BondClass import Bond
import numpy
import CifFile
import SpaceGroups
from CellClass import Cell

class MagneticCell():
    """This class is instantiated with a crystallographic unit cell (complete with atoms)
    and the dimensions of the magnetic unit cell.  Bonds are then added using
    the addBond() method.  In short, atoms are added to and stored in the
    crystallographic cell (Cell class), but "bonds" are added to and stored in
    this class.

    Right now this class is being treated as the cutoff cell.  We don't
    do anything at the moment where we need to know the magnetic cell, so i
    haven't gotten around to changing the names yet.  It seems like it might
    suffice to have the cutoff cell(which would be identical to this except
    for the name, and then simply store the dimensions of the magnetic cell."""
    
    def __init__(self, Unit_Cell, Na, Nb, Nc, spaceGroup):
        """The magnetic cell class is created from an existing crystallographic
        unit cell.  That unit cell is then translated by the dimensions
        Na, Nb, Nc to create the magnetic cell.
        
        Unit_Cell is the initial crystallographic unit cell.
        SpaceGroup is an instance of the SpaceGroup class.
        Na, Nb, Nc are the dimensions(integers) of the magnetic unit cell in
        units of the the crystallographic unit cell axes(a,b,c).(the number
        of times the crystallographic unit celll is translated in each direction
        to create the magnetic unit cell)"""

        self.unit_cell = Unit_Cell
        
        #This is a list of all the crystallographic unit cells that make up
        #this magnetic unit cell
        self.AllUnitCells = [Unit_Cell]
        
        #bonds (interactions) are stored in this class since it is treated
        #as the cutoff cell.(The cell containing all unique bonds)
        self.Bonds = []
        
        self.space_Group = spaceGroup
        
        self.Na = Na
        self.Nb = Nb
        self.Nc = Nc
        
        #generate the magnetic Cell by translating the unit cell
        for i in range(0, Na):
            for j in range(0, Nb):
                for k in range(0, Nc):
                    if i !=0 or j != 0 or k != 0: #to not duplicate original unit cell
                        self.AllUnitCells.append(self.unit_cell.translateCell(i,j,k))
    

        #Recording bonds that are mapped back onto themselves by a symOp
        self.bondConstraints = []


    def addBond(self, Atom1, Atom2, jMatrix = None):     
        """Adds a bond and all symmettry equivalent bonds within the magnetic Cell
        
        Performs every symmetry operation and every possible translation for each 
        symmetry operation to find all possible bonds"""
        
        #Create Symmetry Bonds
        
        xyz = Atom1.getPosition()
        xyz2 = Atom2.getPosition()
        
        for symop in self.space_Group.iter_symops():
        # operate on coordinates in non-shifted spacegroup
            pos1 = symop(xyz)
            pos2 = symop(xyz2)
            
            #Recording symops that map bond back to itself (for symmetry constraints)
            if xyz[0] == pos1[0] and xyz2[0] == pos2[0]: #x component
                if xyz[1] == pos1[1] and xyz2[1] == pos2[1]: #y component
                    if xyz[2] == pos1[2] and xyz2[2] == pos2[2]: #z component
                        self.bondConstraints.append(BondConstraint(xyz, xyz2, symop))
            
            
            mask1 = numpy.logical_or(pos1 < 0.0, pos1 >= 1.0)
            translation = numpy.floor(pos1[mask1])  #translates the first atom back to cell at (0,0,0)
            pos1[mask1] -= translation
            pos2[mask1] -= translation  #Uses same translation to translate other atom
                 
            
            #translate new Bond to all cells 
            
            #Find the number of times to translate the bond in each direction so
            #that it will will cover the whole magnetic cell but not go outside
            if pos1[0]>pos2[0]:
                Xiter = self.Na - pos1[0]
            else:
                Xiter = self.Na - pos2[0]
            
            if pos1[1] > pos2[1]:
                Yiter = self.Nb - pos1[1]
            else:
                Yiter = self.Nb - pos2[1]
            
            if pos1[2] > pos2[2]:
                Ziter = self.Nc - pos1[2]
            else:
                Ziter = self.Nc - pos2[2]
            
            Xiter = int(Xiter) + 1 
            Yiter = int(Yiter) + 1
            Ziter = int(Ziter) + 1
            
            #iterate through each possible translation and check if there are atoms there that could
            #be bonded; if so, add the bond
            for i in range(0, Xiter): #translate in x direction (Na - Cell X position) times
                for j in range(0, Yiter): #translate in y direction (Nb - Cell Y position) times
                    for k in range(0, Ziter): #translate in z direction (Nc - Cell Z position) times
                        translatedAtom1 = self.atomAtPosition((i + pos1[0],j + pos1[1],k + pos1[2]))
                        translatedAtom2 = self.atomAtPosition((i + pos2[0],j + pos2[1],k + pos2[2]))
                        
                        if translatedAtom1 != None and translatedAtom2 != None:
                            newBond = Bond(translatedAtom1, translatedAtom2, jMatrix)
                        
                            #make sure bond does not already exist
                            for currentBond in self.Bonds:
                                if newBond.sameBond(currentBond):
                                    break
                            else:  #if not, add the bond to the list of unique bonds
                                self.Bonds.append(newBond)
  
    def deleteBond(self, bond):
        """Removes a bond and all symmetry equivalent bonds within the magnetic Cell
        
        Performs every symmetry operation and every possible translation for each 
        symmetry operation to find all possible bonds

        This works pretty much the same way as addBond()."""
        
        #Find Symmetry Bonds
        xyz = bond.getAtom1().getPosition()
        xyz2 = bond.getAtom2().getPosition()
        
        
        for symop in self.space_Group.iter_symops():
        # operate on coordinates in non-shifted spacegroup
            pos1 = symop(xyz)
            pos2 = symop(xyz2)
            
            mask1 = numpy.logical_or(pos1 < 0.0, pos1 >= 1.0)
            translation = numpy.floor(pos1[mask1])  #translates the first atom back to cell at (0,0,0)
            pos1[mask1] -= translation
            pos2[mask1] -= translation  #Uses same translation to translate other atom
            
            
            #translate new Bond to all cells 
            
            if pos1[0]>pos2[0]:
                Xiter = self.Na - pos1[0]
            else:
                Xiter = self.Na - pos2[0]
            
            if pos1[1] > pos2[1]:
                Yiter = self.Nb - pos1[1]
            else:
                Yiter = self.Nb - pos2[1]
            
            if pos1[2] > pos2[2]:
                Ziter = self.Nc - pos1[2]
            else:
                Ziter = self.Nc - pos2[2]
            
            Xiter = int(Xiter) + 1 
            Yiter = int(Yiter) + 1
            Ziter = int(Ziter) + 1
            
            for i in range(0, Xiter): #translate in x direction (Na - Cell X position) times
                for j in range(0, Yiter): #translate in y direction (Nb - Cell Y position) times
                    for k in range(0, Ziter): #translate in z direction (Nc - Cell Z position) times
                        translatedAtom1 = self.atomAtPosition((i + pos1[0],j + pos1[1],k + pos1[2]))
                        translatedAtom2 = self.atomAtPosition((i + pos2[0],j + pos2[1],k + pos2[2]))
                        

                        
                        if translatedAtom1 != None and translatedAtom2 != None:
                            #Find bond that connects these and delete
                            for bond in self.Bonds:
                                if bond.getAtom1() == translatedAtom1 or bond.getAtom1() == translatedAtom2:
                                    if bond.getAtom2() == translatedAtom1 or bond.getAtom2() == translatedAtom2:
                                         self.Bonds.remove(bond)
    
    def clearAllBonds(self):
        self.Bonds = []
    
    def positionsInSameCell(self, pos1, pos2):
        """Returns true if the two positions (x,y,z) are in the same crystallographic unit cell"""
        x1,y1,z1 = pos1
        x2,y2,z2 = pos2
        if int(x1) == int(x2):
            if int(y1) == int(y2):
                if int(z1) == int(z2):
                    return True
    
    def getAllUnitCells(self):
        return self.AllUnitCells
    
    def cellAtPosition(self, Position):
        """Returns the crystallographic unit cell containing the given position,
        or None if there is no cell at that position (in this magnetic cell)."""
        for cell in self.AllUnitCells:
            cellPos = cell.getPosition()
            if cellPos[0] <= Position[0] and (cellPos[0]+1) > Position[0]: #check x
                if cellPos[1] <= Position[1] and (cellPos[1]+1) > Position[1]: #check y
                    if cellPos[2] <= Position[2] and (cellPos[2]+1) > Position[2]: #check z
                        return cell
        return None  #No Cell At this Position
    
    def atomAtPosition(self, Position):
        """Returns the atom located at Position (a,b,c), or None if there is no
        atom at that position or if the position is not in this magnetic cell."""
        cell = self.cellAtPosition(Position)
        if cell != None:
            return cell.atomAtPosition(Position)
        return None #If there is no Cell there are no Atoms
            
    def getAllAtoms(self):
        """Returns a list of all the atoms in this magnetic cell"""
        atoms = []
        for cell in self.getAllUnitCells():
            atoms += cell.getAtoms()
        return atoms
    
    def getBonds(self):
        return self.Bonds
    
    def hasBond(self, atom1, atom2, jMatrix = None):
        """Returns true if this magnetic cell contains a bond linking atom1 and
        atom2 with JMatrix or false otherwise"""
        for eachBond in self.getBonds():
            if eachBond.getAtom1() == atom1 or eachBond.getAtom2() == atom1:
                if eachBond.getAtom1() == atom2 or eachBond.getAtom2() == atom2:
                    if jMatrix != None and eachBond.getJMatrix() != None:
                        if eachBond.getJMatrix().all() == jMatrix.all():
                            print "True"
                            return True
                    elif (not jMatrix) and (not eachBond.getJMatrix()):
                        print "True"
                        return True
        return False

    def getNa(self):
        return self.Na
    
    def getNb(self):
        return self.Nb
    
    def getNc(self):
        return self.Nc
    

class BondConstraint():
    def __init__(self, pos1, pos2, symop):
        self.pos1 = pos1
        self.pos2 = pos2
        self.symop = symop
  
  
  
#def magneticCellFromCif(filename):
#    cf = CifFile.ReadCif(filename)
#    
#    #Assuming all data is in one outter block like NIST examples:
#    data = cf[cf.keys()[0]]
#    
#    #Create a Crystollographic Unit Cell
#    a = data['_cell_length_a']
#    b = data['_cell_length_b']
#    c = data['_cell_length_c']
#    
#    alpha = data['_cell_angle_alpha']
#    gamma = data['_cell_angle_gamma']
#    beta = data['_cell_angle_beta']
#    
#    spaceGroup = SpaceGroups.GetSpaceGroup(int(data['_symmetry_Int_Tables_number']))
#    
#    unitcell = Cell(spaceGroup, 0,0,0, a, b, c, alpha, gamma, beta)
#    
#    atomLabels = data['_atom_site_label']
##Not Currently used        atomType = data['_atom_site_type_symbol']
#    xPositions = data['_atom_site_fract_x']
#    yPositions = data['_atom_site_fract_y']
#    zPositions = data['_atom_site_fract_z']
#    
#    for i in range(len(atomLabels)):
#        unitcell.generateAtoms((float(xPositions[i]), float(yPositions[i]), float(zPositions[i])), atomLabels[i])

    #Create a Magnetic Cell
#    return MagneticCell(unitcell, 1,1,1, spaceGroup)
    
