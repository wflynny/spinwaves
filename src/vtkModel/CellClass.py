import SymmetryUtilities
import numpy
from AtomClass import *
from BondClass import *
import random

class Cell():
    def  __init__(self, Space_Group, PosX = 0, PosY = 0, PosZ = 0, a=1, b=1, c=1, alpha=90, gamma=90, beta=90):
        """PosX, PosY, PosZ are the fractional coordinates of the cell - they should all be integers"""
        self.Space_Group = Space_Group
        
        #These should be integers
        self.PosX = PosX#position world coordinates in vtk renderer
        self.PosY = PosY
        self.PosZ = PosZ
        
        #Dimensions (Currently not used)
        self.a = a
        self.b = b
        self.c = c
        self.alpha = alpha
        self.gamma = gamma
        self.beta = beta
        
        self.Atoms = []
        
        
    
    #functions
    def getSpaceGroup(self):
        return self.Space_Group
    
    def atomAtPosition(self,position):
        """Returns the atom at the position if one exists, None otherwise"""
        if self.positionIsInCell(position):
            positionList = []
            for atom in self.Atoms:
                positionList.append(atom.getPosition())
            closest = self.Atoms[SymmetryUtilities.nearestSiteIndex(positionList,position)]
            if SymmetryUtilities.equalPositions(closest.getPosition(), position):
                return closest
        return None
    
    def positionIsInCell(self, position):
        if self.PosX <= position[0] and (self.PosX+1) > position[0]: #check x
            if self.PosY <= position[1] and (self.PosY+1) > position[1]: #check y
                if self.PosZ <= position[2] and (self.PosZ+1) > position[2]: #check z
                    return True
        return False
    
    def addAtom(self, Atom):
        self.Atoms.append(Atom)
    
#    def addBond(self, Bond):
#        self.Bonds.append(Bond)
    
    def getAtoms(self):
        return self.Atoms
    
#    def getBonds(self):
#        return self.Bonds
    
    def setPosX(self, x):
        self.PosX = x
    
    def setPosY(self, y):
        self.PosY = y
    
    def setPosZ(self, z):
        self.PosZ = z
        
    def getPosition(self):
        return (self.PosX, self.PosY, self.PosZ)
                   
    def translateCell(self, a, b, c):
        new_cell = Cell(self.Space_Group,a,b,c)
        for atomn in self.Atoms:  #should preserve order of Atoms
            position = atomn.getPosition()
            color = atomn.getColor()
            new_cell.addAtom(Atom(new_cell, position[0], position[1], position[2], atomn.getDescription(), atomn.getRadius(), color[0], color[1], color[2]))

#Bonds no longer associated with cells
#        for bondn in self.Bonds:
#            newAtom1 = new_cell.atomAtIndex( self.getAtomIndex(bondn.getAtom1()) )
#            newAtom2 = new_cell.atomAtIndex( self.getAtomIndex(bondn.getAtom2()) )
#            new_cell.addBond(Bond(new_cell, newAtom1, newAtom2, bondn.getRGBColor()))
        
        return new_cell

    def __str__(self):
        return "unit cell at (" + str(self.PosX) + ", " + str(self.PosY) + ", " + str(self.PosZ) + ")"
    
    def atomAtIndex(self, i):
        return self.Atoms[i]
    
    def getAtomIndex(self, atom):
        return self.Atoms.index(atom)
    
    def generateAtoms(self, position, description, radius = .05):
        locations = SymmetryUtilities.expandPosition(self.Space_Group, numpy.array([position[0],position[1], position[2]]))[0]
        randGen = random.Random()
        r = random.uniform(0,1)
        g = random.uniform(0,1)
        b = random.uniform(0,1)
        for coord in locations:
    #        r,g,b = atom1.getColor()
            atom = Atom(self, coord[0], coord[1], coord[2], description, radius, r,g,b)
            self.addAtom(atom)

    def getA(self):
        return self.a
    
    def getB(self):
        return self.b
    
    def getC(self):
        return self.c
    
    def getAlpha(self):
        return self.alpha
    
    def getBeta(self):
        return self.beta
    
    def getGamma(self):
        return self.gamma
    

