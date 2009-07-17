class Atom():
    def __init__(self, unit_Cell, x,y,z, atomicNum, description = "", valence = "", 
                 radius=.05,r = 1,g=0,b = 0, spin = None, anisotropy = (0,0,0),
                 spinMagnitude = 1):
        """
        x,y,z are fractional coordinates in the unit cell
        unit_Cell is the unit cell containing the atom (instance of Cell class)
        radius is the radius of the sphere actor (for drawing purposes later)
        r,g,b are the red, green, blue values of the actor
        spin is a tuple(Sx, Sy, Sz)
        anisotropy is the single ion anisotropy of the atom (Dx, Dy, Dz)
        description is string describing the atom such as a name.
        """
#===============================================================================
#        print '\n\nunit_cell:', unit_Cell
#        print 'x,y,z:', x,y,z
#        print 'atomicNum:', atomicNum
#        print 'description:', description
#        print 'radius:', radius
#        print 'r,g,b:', r,g,b
#        print 'spin:' , spin
#        print 'anisotropy: ', anisotropy
#        print 'spinMag: ', spinMagnitude
#===============================================================================
        self.anisotropy = anisotropy
        self.description = description
        self.radius = radius
        self.color = (r,g,b)
        self.spin = spin
        self.spinMagnitude = spinMagnitude
        self.atomicNumber = atomicNum
        self.valence = valence
        
        if x<1 and y<1 and z<1: 
            #coordinates  (within cell) - Actor will contain world coordinates in vtk renderer
            self.a = x
            self.b = y
            self.c = z
        else:
            raise Exception("a,b,c are fractional coordinates within the unit cell.  They should be less than 1.")
        
        #Unit Cell Containing this Atom
        self.unit_Cell = unit_Cell
        

    def getAnisotropy(self):
        return self.anisotropy
    
    def getSpin(self):
        return self.spin
    
    def getSpinMagnitude(self):
        return self.spinMagnitude
    
    def setSpin(self, spin):
        self.spin = spin
    
    def getRadius(self):
        return self.radius
    
    def getColor(self):
        """returns (r,g,b) color"""
        return self.color
    
    def getPosition(self):
        """returns (a,b,c) fractional coordinates within Unit Cell + (a,b,c) Unit Cell Position  (integers)"""
        cellPos = self.unit_Cell.getPosition()
        return (self.a + cellPos[0], self.b + cellPos[1], self.c + cellPos[2])
    
    def getDescription(self):
        return self.description
    
    def getUnitCell(self):
        return self.unit_Cell
    
    def __str__(self):
        return self.getDescription().rstrip() + " at " + str(self.getPosition()) + " in " + self.unit_Cell.__str__() + "  Spin = " + str(self.spinMagnitude)
    
    def getIndexNumber(self):
        """Returns the Atom's Index Number in the Unit Cell
        The index is used as a unique identifier when translating the cell (to associate new and old) and when drawing"""
        return self.unit_Cell.getAtomIndex(self) 
        
