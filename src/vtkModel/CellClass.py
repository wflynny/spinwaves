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
        #draw very Light Box 
        #Create "Cube" Source
        box = vtkCubeSource()
        box.SetXLength(1)
        box.SetYLength(1)
        box.SetZLength(1)
        
        boxMap = vtkPolyDataMapper()
        boxMap.SetInput(box.GetOutput())
        
        #create actor
        abox = vtkActor()
        abox.SetMapper(boxMap)
        abox.GetProperty().SetColor(0,.1,.6)
        abox.GetProperty().SetOpacity(.1)
        abox.SetPosition(self.PosX + .5, self.PosY + .5, self.PosZ + .5)
        abox.PickableOff()
        
        #Add the actor to the renderer 
        renderer.AddActor(abox)
        
        #Add the Atoms and Bonds to the renderer
        for atomn in self.Atoms:
            renderer.AddActor(atomn.getActor())
        for bondn in self.Bonds:
            renderer.AddActor(bondn.getActor())
            
            
            
    def translateCell(self, a, b, c):
        new_cell = Cell(a,b,c)
        for atomn in self.Atoms:  #should preserve order of Atoms
            position = atomn.getPosition()
            color = atomn.getActor().GetProperty().GetColor()
            new_cell.addAtom(Atom(new_cell, position[0], position[1], position[2], atomn.getDescription(), atomn.getSource().GetRadius(), color[0], color[1], color[2]))

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
