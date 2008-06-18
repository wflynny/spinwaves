class Atom():
    def __init__(self, unit_Cell, x,y,z, description = "", radius=.05,r = 1,g=0,b = 0):
        """
        x,y,z are fractional coordinates in the unit cell
        unit_Cell is the unit cell containing the atom (instance of Cell class)
        radius is the radius of the sphere actor
        r,g,b are the red, green, blue values of the actor
        """
        self.description = description
        self.radius = radius
        self.color = [r,g,b]
        
        if x<1 and y<1 and z<1: 
            #coordinates  (within cell) - Actor will contain world coordinates in vtk renderer
            self.x = x
            self.y = y
            self.z = z
        else:
            raise Exception("x,y,z are fractional coordinates within the unit cell.  They should be less than 1.")
        
        #Unit Cell Containing this Atom
        self.unit_Cell = unit_Cell
        

    def getRadius(self):
        return self.radius
    
    def getColor(self):
        """returns [r,g,b] color"""
        return self.color
    
    def getPosition(self):
        """returns (x,y,z) fractional coordinates within Unit Cell + (x,y,z) Unit Cell Position  (integers)"""
        cellPos = self.unit_Cell.getPosition()
        return (self.x + cellPos[0],self.y + cellPos[1],self.z + cellPos[2])
    
    def getDescription(self):
        return self.description
    
    def getUnitCell(self):
        return self.unit_Cell
    
    def __str__(self):
        return self.getDescription().rstrip() + " at " + str(self.getPosition()) + " in " + self.unit_Cell.__str__()
    
    def getIndexNumber(self):
        """Returns the Atom's Index Number in the Unit Cell"""
        return self.unit_Cell.getAtomIndex(self) 
        
