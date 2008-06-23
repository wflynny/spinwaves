import copy
from BondClass import *
import numpy

class MagneticCell():
    
    def __init__(self, Unit_Cell, Na, Nb, Nc, spaceGroup):
        self.unit_cell = Unit_Cell
        
        self.AllUnitCells = [Unit_Cell]
        
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
    

        self.bondConstraints = []
    
    
    #does not currently support anything but default bond colors
#    def addInterCellularBond(self, Atom1, Atom2):
     
        #Make Sure the Atoms are not in the same cell
 #       if Atom1.getUnitCell() == Atom2.getUnitCell():
 #           raise Exception("These atoms are in the same Unit Cell:" + Atom1.__str__() + ", " + Atom2.__str__())
    def addBond(self, Atom1, Atom2):     
        """Adds a bond and all symmettry equivalent bonds within the magnetic Cell
        
        Performs every symmetry operaation and every possible translation for each 
        symmetry operation to find all possible bonds"""
        
        #Create Symmetry Bonds
        
        xyz = Atom1.getPosition()
        xyz2 = Atom2.getPosition()
        
        
        for symop in self.space_Group.iter_symops():
        # operate on coordinates in non-shifted spacegroup
            pos1 = symop(xyz)
            pos2 = symop(xyz2)
            
            #Recording symops that map bond back to itself
            if xyz[0] == pos1[0] and xyz2[0] == pos2[0]: #x component
                if xyz[1] == pos1[1] and xyz2[1] == pos2[1]: #y component
                    if xyz[2] == pos1[2] and xyz2[2] == pos2[2]: #z component
                        self.bondConstraints.append(BondConstraint(xyz, xyz2, symop))
            
            
            
            mask1 = numpy.logical_or(pos1 < 0.0, pos1 >= 1.0)
            translation = numpy.floor(pos1[mask1])  #translates the first atom back to cell at (0,0,0)
            pos1[mask1] -= translation
            pos2[mask1] -= translation  #Uses same translation to translate other atom
                 

                
            
#            atom1Index = self.unit_cell.atomAtPosition(pos1).getIndexNumber()
#            atom2Index = self.atomAtPosition(pos2).getIndexNumber()
            
            
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
                            newBond = Bond(translatedAtom1, translatedAtom2)
                        
                            #make sure bond does not already exist
                            for currentBond in self.Bonds:
                                if newBond.sameBond(currentBond):
                                    break
                            else:  #if not, add the bond to the list of unique bonds
                                self.Bonds.append(newBond)

                            
    
    def deleteBond(self, bond):
        """Removes a bond and all symmetry equivalent bonds within the magnetic Cell
        
        Performs every symmetry operaation and every possible translation for each 
        symmetry operation to find all possible bonds"""
        
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
    
    def positionsInSameCell(self, pos1, pos2):
        x1,y1,z1 = pos1
        x2,y2,z2 = pos2
        if int(x1) == int(x2):
            if int(y1) == int(y2):
                if int(z1) == int(z2):
                    return True
    
    def getAllUnitCells(self):
        return self.AllUnitCells
    
    def cellAtPosition(self, Position):
        for cell in self.AllUnitCells:
            cellPos = cell.getPosition()
            if cellPos[0] <= Position[0] and (cellPos[0]+1) > Position[0]: #check x
                if cellPos[1] <= Position[1] and (cellPos[1]+1) > Position[1]: #check y
                    if cellPos[2] <= Position[2] and (cellPos[2]+1) > Position[2]: #check z
                        return cell
        return None  #No Cell At this Position
    
    def atomAtPosition(self, Position):
        cell = self.cellAtPosition(Position)
        if cell != None:
            return cell.atomAtPosition(Position)
        return None #If there is no Cell there are no Atoms
            
    def getAllAtoms(self):
        atoms = []
        for cell in self.getAllUnitCells():
            atoms += cell.getAtoms()
        return atoms
    
    def getBonds(self):
        return self.Bonds
    
class BondConstraint():
    def __init__(self, pos1, pos2, symop):
        self.pos1 = pos1
        self.pos2 = pos2
        self.symop = symop
