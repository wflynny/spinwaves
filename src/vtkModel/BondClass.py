import AtomClass
from vtk import *
import scipy
import math
import numpy


class Bond():
    def __init__(self, Cell, Atom1, Atom2, r = 0,g = 0,b = 1):
        """Cell is the Unit Cell if this bond is in a unit Cell or None Otherwise.
        Atom1 and Atom2 are the atoms that this bond connects.
        r,g,b are the color of the actor."""
        self.r = r#color
        self.g = g
        self.b = b
        
        self.Atom1 = Atom1
        self.Atom2 = Atom2
        self.Cell = Cell
        self.cylinder_actor = self.makeCylinder(Atom1.getActor(), Atom2.getActor(), Atom1.getSource().GetRadius(), Atom2.getSource().GetRadius())
    

    defualtRadius = .025

    
    def crossProduct(self, x1,y1,z1,x2,y2,z2):
        """x,y,z represent vector coordinates"""
        i = (y1*z2) - (z1*y2)
        j = (x1*z2) - (z1*x2)
        k = (x1*y2) - (y1*x2)
        return i,j,k

    def makeCylinder (self, SphereOne, SphereTwo, rad_One, rad_Two):
        """SphereOne and SphereTwo are two spherical actors rad_One and rad_Two are the radii of the spheres
        returns a cylindrical Actor"""
        
        posOne = SphereOne.GetPosition()
        posTwo = SphereTwo.GetPosition()
        #coordinates of vector pointing from SphereTwo to SphereOne
        x = posTwo[0] - posOne[0]
        y = posTwo[1] - posOne[1]
        z = posTwo[2] - posOne[2]
        distance = ((x**2 + y**2 + z**2)**.5 -rad_One - rad_Two)
        
        #create cylinder
        cylinder = vtkCylinderSource()
        cylinder.SetRadius(self.defualtRadius)
        cylinder.SetResolution(100)
        
        #map to another map
        cylinderMap = vtkPolyDataMapper()
        cylinderMap.SetInput(cylinder.GetOutput())
        
        #create cylinder actor
        aCylinder = vtkActor()
        aCylinder.SetMapper(cylinderMap)
        aCylinder.GetProperty().SetColor(0,.1,.6)
        aCylinder.SetScale(.2,distance,.2)
        
        #to get the center point between the surfaces of the two spheres
        #create a unit vector and multiply it by the distance/2
        
        #create unit vetor
        vectLength = (x**2 + y**2 + z**2)**.5
        print vectLength
        print posOne
        print posTwo
        unitX = x/vectLength
        unitY = y/vectLength
        unitZ = z/vectLength
    
        #get center coordinates (center point between two sphere surfaces)
        centerX = (unitX * (rad_One + distance/2)) + posOne[0]
        centerY = (unitY * (rad_One + distance/2)) + posOne[1]
        centerZ = (unitZ * (rad_One + distance/2)) + posOne[2]
        aCylinder.SetPosition(centerX, centerY, centerZ)
        
        
        #Angle
        #cross product of Unit vectors cylinder direction and desired orientation
        i,j,k = self.crossProduct(0,1,0,unitX,unitY,unitZ)  #default orientation for the cylinder is along y axis
        theta = scipy.arcsin((i**2+j**2+k**2)**.5)*180/math.pi
        
        #if the angle is obtuse, theta must be corrected
        if self.dotProduct(0,1,0,unitX,unitY,unitZ) >= 0:
            aCylinder.RotateWXYZ(theta,i,j,k)
        else:
            aCylinder.RotateWXYZ(theta,-i,-j,-k)
        
        return aCylinder
    
    def getActor(self):
        return self.cylinder_actor
    
    def getAtom1(self):
        return self.Atom1
    
    def getAtom2(self):
        return self.Atom2
    
    def setAtom1(self, atom):
        self.Atom1 = atom
    
    def setAtom2(self, atom):
        self.Atom2 = atom

    def createSymmetryBonds(self, space_Group):
        """returns list of new bonds within the same cell formed from space_group symOps
        Can Only be used on a Bond that is contained within a cell (not for bonds between cells)
        -Should only really be used on original unit Cell"""
        newBonds = []
        
        xyz = self.Atom1.getPosition()
        xyz2 = self.Atom2.getPosition()
        for symop in space_Group.iter_symops():
        # operate on coordinates in non-shifted spacegroup
            pos1 = symop(xyz)
            pos2 = symop(xyz2)
            mask1 = numpy.logical_or(pos1 < 0.0, pos1 >= 1.0)
            mask2 = numpy.logical_or(pos2 < 0.0, pos2 >= 1.0)
            pos1[mask1] -= numpy.floor(pos1[mask1])
            pos2[mask2] -= numpy.floor(pos2[mask2])
            print pos1
            print pos2
            newBond = Bond(self.Cell, self.Cell.atomAtPosition(pos1), self.Cell.atomAtPosition(pos2), self.r, self.g, self.b)
            #check if the bond already exists
            for currentBond in newBonds:
                if newBond.sameBond(currentBond):
                    break
            else:  #if not, add the bond to the list
                newBonds.append(newBond)
            
        return newBonds

    def sameBond(self, otherBond):
        """returns true if otherBond connects the same atoms"""
        if self.getAtom1() == otherBond.getAtom1() or self.getAtom1() == otherBond.getAtom2():
            if self.getAtom2() == otherBond.getAtom1() or self.getAtom2() == otherBond.getAtom2():
                return True
        return False
    
    def __str__(self):
        return "Bond between " + self.Atom1.__str__() + " and " + self.Atom2.__str__()
    
    def dotProduct(self, x1,y1,z1,x2,y2,z2):
        """x,y,z represent vector coordinates"""
        return (x1*x2 + y1*y2 + z1*z2)
    
    def getRGBColor(self):
        return self.r, self.g, self.b

