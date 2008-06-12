import SymmetryUtilities
import numpy
from AtomClass import *
from BondClass import *

class Cell():
    def  __init__(self, PosX = 0, PosY = 0, PosZ = 0, DimX = 1, DimY = 1, DimZ = 1):
        
        #These should be integers
        self.PosX = PosX#position world coordinates in vtk renderer
        self.PosY = PosY
        self.PosZ = PosZ
        
        #Dimensions
        self.DimX = DimX
        self.DimY = DimY
        self.DimZ = DimZ
        
        self.Atoms = []
        self.Bonds = []
        


    
    #functions
    def atomAtPosition(self,position):  #untested
        positionList = []
        for atom in self.Atoms:
            positionList.append(atom.getPosition())
        return self.Atoms[SymmetryUtilities.nearestSiteIndex(positionList,position)]
    
    def addAtom(self, Atom):
        self.Atoms.append(Atom)
    
    def addBond(self, Bond):
        self.Bonds.append(Bond)
    
    def getAtoms(self):
        return self.Atoms
    
    def getBonds(self):
        return self.Bonds
    
    def setPosX(self, x):
        self.PosX = x
    
    def setPosY(self, y):
        self.PosY = y
    
    def setPosZ(self, z):
        self.PosZ = z
        
    def getPosition(self):
        return (self.PosX, self.PosY, self.PosZ)
         
    def drawCell(self, renderer):
        for atomn in self.Atoms:
            renderer.AddActor(atomn.getActor())
        for bondn in self.Bonds:
#           x,y,z = bondn.getActor().GetPosition()
#            bondn.getActor().SetPosition(x + self.PosX, y + self.PosY, z + self.PosZ)
            renderer.AddActor(bondn.getActor())
            
    def translateCell(self, a, b, c):
        new_cell = Cell(a,b,c)
        for atomn in self.Atoms:  #should preserve order of Atoms
            print "this"
            position = atomn.getPosition()
            color = atomn.getActor().GetProperty().GetColor()
            new_cell.addAtom(Atom(new_cell, position[0], position[1], position[2], atomn.getDescription(), atomn.getSource().GetRadius(), color[0], color[1], color[2]))
        
        print "here  in translate cell atoms done"
        for bondn in self.Bonds:
            newAtom1 = new_cell.atomAtIndex( self.getAtomIndex(bondn.getAtom1()) )
            newAtom2 = new_cell.atomAtIndex( self.getAtomIndex(bondn.getAtom2()) )
            new_cell.addBond(Bond(new_cell, newAtom1, newAtom2, bondn.getRGBColor()))
        
        return new_cell

    def __str__(self):
        return "unit cell at (" + str(self.PosX) + ", " + str(self.PosY) + ", " + str(self.PosZ) + ")"
        
    def test(self):
        return Cell()
    
    def atomAtIndex(self, i):
        return self.Atoms[i]
    
    def getAtomIndex(self, atom):
        return self.Atoms.index(atom)
