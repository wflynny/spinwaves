import AtomClass
from vtk import *
import numpy



class Bond():
    def __init__(self, Atom1, Atom2, jMatrix = None, r = 0,g = 0,b = 1):
        """Cell is the Unit Cell if this bond is in a unit Cell or None Otherwise.
        Atom1 and Atom2 are the atoms that this bond connects.
        r,g,b are the color of the actor."""
        self.r = r#color
        self.g = g
        self.b = b
        
        self.Atom1 = Atom1
        self.Atom2 = Atom2
        
        self.jMat = jMatrix

    
    def getAtom1(self):
        return self.Atom1
    
    def getAtom2(self):
        return self.Atom2
    
    def setAtom1(self, atom):
        self.Atom1 = atom
    
    def setAtom2(self, atom):
        self.Atom2 = atom

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

