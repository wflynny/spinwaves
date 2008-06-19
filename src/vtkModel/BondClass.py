import AtomClass
from vtk import *
import numpy



class Bond():
    def __init__(self, Atom1, Atom2, r = 0,g = 0,b = 1):
        """Cell is the Unit Cell if this bond is in a unit Cell or None Otherwise.
        Atom1 and Atom2 are the atoms that this bond connects.
        r,g,b are the color of the actor."""
        self.r = r#color
        self.g = g
        self.b = b
        
        self.Atom1 = Atom1
        self.Atom2 = Atom2
#        self.Cell = Cell

    
    def getAtom1(self):
        return self.Atom1
    
    def getAtom2(self):
        return self.Atom2
    
    def setAtom1(self, atom):
        self.Atom1 = atom
    
    def setAtom2(self, atom):
        self.Atom2 = atom



#Handled in MagneticCell Class now

#    def createSymmetryBonds(self, space_Group):
#        """returns list of new bonds within the same cell formed from space_group symOps
#        Can Only be used on a Bond that is contained within a cell (not for bonds between cells)
#        -Should only really be used on original unit Cell"""
#        newBonds = []
#        
#        xyz = self.Atom1.getPosition()
#        xyz2 = self.Atom2.getPosition()
#        for symop in space_Group.iter_symops():
#        # operate on coordinates in non-shifted spacegroup
#            pos1 = symop(xyz)
#            pos2 = symop(xyz2)
#            mask1 = numpy.logical_or(pos1 < 0.0, pos1 >= 1.0)
#            mask2 = numpy.logical_or(pos2 < 0.0, pos2 >= 1.0)
#            pos1[mask1] -= numpy.floor(pos1[mask1])
#            pos2[mask2] -= numpy.floor(pos2[mask2])
#            newBond = Bond(self.Cell, self.Cell.atomAtPosition(pos1), self.Cell.atomAtPosition(pos2), self.r, self.g, self.b)
#            #check if the bond already exists
#            for currentBond in newBonds:
#                if newBond.sameBond(currentBond):
#                    break
#            else:  #if not, add the bond to the list
#                newBonds.append(newBond)
#            
#        return newBonds


    def sameBond(self, otherBond):
        """returns true if otherBond connects the same atoms"""
        if self.getAtom1() == otherBond.getAtom1() or self.getAtom1() == otherBond.getAtom2():
            if self.getAtom2() == otherBond.getAtom1() or self.getAtom2() == otherBond.getAtom2():
                return True
        return False
    
    def __str__(self):
        return "Bond between " + self.Atom1.__str__() + " and " + self.Atom2.__str__()
    
    def getRGBColor(self):
        return self.r, self.g, self.b

