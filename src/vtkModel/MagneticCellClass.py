import copy
from BondClass import *

class MagneticCell():
    
    def __init__(self, Unit_Cell, Na, Nb, Nc, spaceGroup):
        self.unit_cell = Unit_Cell
        
        self.AllUnitCells = [Unit_Cell]
        
        self.IntercellularBonds = []
        
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
    
    
    
    #does not currently support anythin but default bond colors
    def addInterCellularBond(self, Atom1, Atom2):
        
        #Make Sure the Atoms are not in the same cell
        if Atom1.getUnitCell() == Atom2.getUnitCell():
            raise Exception("These atoms are in the same Unit Cell:" + Atom1.__str__() + ", " + Atom2.__str__())
        
        newBonds = [] #this bond and bonds created by symmetry
        original = Bond(None, Atom1, Atom2)
        
        #Create Symmetry Bonds
        
        xyz = original.getAtom1().getPosition()
        xyz2 = original.getAtom2().getPosition()
        count = 0
        for symop in self.space_Group.iter_symops():
        # operate on coordinates in non-shifted spacegroup
            pos1 = symop(xyz)
            pos2 = symop(xyz2)
            mask1 = numpy.logical_or(pos1 < 0.0, pos1 >= 1.0)
            translation = numpy.floor(pos1[mask1])  #translates the first atom back to cell at (0,0,0)
            pos1[mask1] -= translation
            pos2[mask1] -= translation  #Uses same translation to translate other atom
            
            #if the second atom is not in the magntic cell, do another translation so that it is
            if pos2[0] < 0 or pos2[1] < 0 or pos2[2] < 0:
                mask2 = pos2 < 0.0
                translation = numpy.floor(pos2[mask2])  #translates the first atom back to cell at (0,0,0)
                pos1[mask2] -= translation
                pos2[mask2] -= translation  #Uses same translation to translate other atom
                atomAtPos1 = self.unit_cell.atomAtPosition(pos2) #Second Position was translated to the first cell (0,0,0)
                atomAtPos2 = self.atomAtPosition(pos1)
            else:
                atomAtPos1 = self.unit_cell.atomAtPosition(pos1) #first Position was translated to the first cell (0,0,0)
                atomAtPos2 = self.atomAtPosition(pos2)
            
            print "symmop:\n", symop
            print "position1:", pos1
            if atomAtPos1 != None:
                print atomAtPos1.getPosition()
            print "position2:", pos2
            if atomAtPos2 != None:
                print atomAtPos2.getPosition()
#            translation = numpy.floor(pos1)
#            pos1 -= numpy.floor(pos1)
#            pos2 -= numpy.floor(pos1)
            
            
            #Right Now this allows Bonds to be created within a Unit Cell and stored as intercellular
            if atomAtPos2 != None and atomAtPos1 != None:  #Both Atoms could be translated outside the bounds of the Magnetic Cell
                newBond = Bond(None, atomAtPos1, atomAtPos2)
                #check if the bond already exists
                for currentBond in newBonds:
                    if newBond.sameBond(currentBond):
                        break
                else:  #if not, add the bond to the list of unique bonds
                    newBonds.append(newBond)
                    
                    #translate new Bond to all cells
                    originalAtom1 = newBond.getAtom1()
                    originalAtom2 = newBond.getAtom2()
                    origPos1 = originalAtom1.getPosition()
                    origPos2 = originalAtom2.getPosition()
                    #Using the fact that the translated Cells store the atoms in the same order:
                    for i in range(0, self.Na - newBond.getAtom2().getUnitCell().getPosition()[0]): #translate in x direction (Na - Cell X position) times
                        for j in range(0, self.Nb - newBond.getAtom2().getUnitCell().getPosition()[1]): #translate in y direction (Nb - Cell Y position) times
                            for k in range(0, self.Nc - newBond.getAtom2().getUnitCell().getPosition()[2]): #translate in z direction (Nc - Cell Z position) times
                                translatedCell1 = self.cellAtPosition((i + origPos1[0],j + origPos1[1],k + origPos1[2]))
                                translatedCell2 = self.cellAtPosition((i + origPos2[0],j + origPos2[1],k + origPos2[2]))
                    
                                translatedAtom1 = translatedCell1.atomAtIndex(originalAtom1.getIndexNumber())
                                translatedAtom2 = translatedCell2.atomAtIndex(originalAtom2.getIndexNumber())
                                self.IntercellularBonds.append(Bond(None,translatedAtom1, translatedAtom2))

                                
    
    def drawCell(self, Renderer):
        #draw all Cells
        for cell in self.AllUnitCells:
            cell.drawCell(Renderer)
        self.drawIntercellularBonds(Renderer)
    
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
    
    def drawIntercellularBonds(self, renderer):
        for bond in self.IntercellularBonds:
            renderer.AddActor(bond.getActor())
            
    def getAllAtoms(self):
        atoms = []
        for cell in self.getAllUnitCells():
            atoms += cell.getAtoms()
        return atoms
    
    def getIntercellularBonds(self):
        return self.IntercellularBonds


