import SymmetryUtilities
import numpy
from AtomClass import Atom
from BondClass import Bond
import random

class Cell():
    """This class models a single crystallographic unit cell. (The entire
    lattice could consist of many of these.)

    The class is instatiated with the space group and dimensions(although
    the dimensions are not currently used for anything).
    Then atoms can be added using the generateAtoms() method which will add
    all symmetry equivalent atoms."""
    
    def  __init__(self, Space_Group, PosX = 0, PosY = 0, PosZ = 0, a=1, b=1,
                  c=1, alpha=90, gamma=90, beta=90):
        """PosX, PosY, PosZ are the fractional coordinates of the cell - they
        should all be integers
        Space_Group is an instance of the SpaceGroup class(SpaceGroups.py)

        Bonds are not stored in this class becuase some span more than one
        crystallographic unit cell.  They are stored in the magnetic or cutoff
        cell classes."""
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
        
        #List of atom contained in this unit cell
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
        """Returns true if the given position is in this unit cell and false
        otherwise"""
        if self.PosX <= position[0] and (self.PosX+1) > position[0]: #check x
            if self.PosY <= position[1] and (self.PosY+1) > position[1]: #check y
                if self.PosZ <= position[2] and (self.PosZ+1) > position[2]: #check z
                    return True
        return False
    
    def addAtom(self, Atom):
        self.Atoms.append(Atom)
    
    def getAtoms(self):
        return self.Atoms
    
    def setPosX(self, x):
        self.PosX = x
    
    def setPosY(self, y):
        self.PosY = y
    
    def setPosZ(self, z):
        self.PosZ = z
        
    def getPosition(self):
        return (self.PosX, self.PosY, self.PosZ)
                   
    def translateCell(self, a, b, c):
        """Returns a new unit cell translated by a,b and c in those respective
        directions.

        The new cell will have the translated coordinates and a list of atoms
        with the same indeces(which are used as an identifying characteristic)
        as their cooresponding atoms in this cell."""
        new_cell = Cell(self.Space_Group,a,b,c)
        for atomn in self.Atoms:  #should preserve order of Atoms
            position = atomn.getPosition()
            color = atomn.getColor()
            new_cell.addAtom(Atom(new_cell, position[0], position[1], position[2], atomn.getDescription(), atomn.getRadius(), color[0], color[1], color[2], anisotropy = atomn.getAnisotropy()))
        
        return new_cell

    def __str__(self):
        return "unit cell at (" + str(self.PosX) + ", " + str(self.PosY) + ", " + str(self.PosZ) + ")"
    
    def atomAtIndex(self, i):
        return self.Atoms[i]
    
    def getAtomIndex(self, atom):
        return self.Atoms.index(atom)
    
    def generateAtoms(self, position, description, anisotropy = (0,0,0), radius = .05):
        """Given the information of one atom and the space group associated with
        this Cell object, this method creates all the symmetry equivalent atoms
        and adds them to the list."""
        
        locations = SymmetryUtilities.expandPosition(self.Space_Group, numpy.array([position[0],position[1], position[2]]))[0]

        #Create a random color for the atoms; In the future the color will be chosen
        #from a list based on hte atomic number along with the radius
        randGen = random.Random()
        r = random.uniform(0,1)
        g = random.uniform(0,1)
        b = random.uniform(0,1)
        for coord in locations:
            atom = Atom(self, coord[0], coord[1], coord[2], description, radius, r,g,b, anisotropy = anisotropy)
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
    

