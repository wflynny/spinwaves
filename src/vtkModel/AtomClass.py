from vtk import *

class Atom():
    def __init__(self, unit_Cell, x,y,z, description = "", radius=.05,r = 1,g=0,b = 0):
        """x,y,z are fractional coordinates in the unit cell"""
        self.description = description
        
        if x<1 and y<1 and z<1: 
            #coordinates  (within cell) - Actor will contain world coordinates in vtk renderer
            self.x = x
            self.y = y
            self.z = z
        else:
            raise Exception("x,y,z are fractional coordinates within the unit cell.  They should be less than 1.")
        
        #Unit Cell Containing this Atom
        self.unit_Cell = unit_Cell
        
        #sphere geometry
        self.sphere_Source = vtkSphereSource()
        self.sphere_Source.SetRadius(radius)
        self.sphere_Source.SetThetaResolution(self.defualt_res)
        self.sphere_Source.SetPhiResolution(self.defualt_res)
        
        # map to graphics objects
        sphereMap = vtkPolyDataMapper()
        sphereMap.SetInput(self.sphere_Source.GetOutput())
        
        # actor coordinates geometry, properties, transformation
        self.sphere_Actor = vtkActor()
        self.sphere_Actor.SetMapper(sphereMap)
        self.sphere_Actor.GetProperty().SetColor(r,g,b)
        
        cellPosition = unit_Cell.getPosition()
        self.sphere_Actor.SetPosition(x + cellPosition[0],y + cellPosition[1],z + cellPosition[2])
    

    
    
    #Default Sphere Resolution Theta and Phi
    defualt_res = 100


    
    def getSource(self):
        return self.sphere_Source
    
    def getActor(self):
        return self.sphere_Actor
    
    def getPosition(self):
        """returns (x,y,z) fractional coordinates within Unit Cell + (x,y,z) Unit Cell Position  (integers)"""
        cellPos = self.unit_Cell.getPosition()
        return (self.x + cellPos[0],self.y + cellPos[1],self.z + cellPos[2])
    
    def getDescription(self):
        return self.description
    
    def __str__(self):
        return self.getDescription().rstrip() + " at " + str(self.getPosition()) + " in " + self.unit_Cell.__str__()
        
