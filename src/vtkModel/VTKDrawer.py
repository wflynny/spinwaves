from vtk import *
import scipy
import math
import numpy


class vtkDrawer():
    def __init__(self, renderer):
        self.ren1 = renderer
        self.actors = {}# asscociate actors with their model objects for picking
    
    
    #Default Sphere resolution
    defualt_Sphere_Res = 20
    defualtBondRadius = .025
    
    
    def getObjFromActor(self, actor):
        return self.actors[actor]
    
    
    def drawAtom(self, atom):   
        #sphere geometry
        sphere_Source = vtkSphereSource()
        sphere_Source.SetRadius(atom.getRadius())
        sphere_Source.SetThetaResolution(self.defualt_Sphere_Res)
        sphere_Source.SetPhiResolution(self.defualt_Sphere_Res)
            
        #map to graphics objects
        sphereMap = vtkPolyDataMapper()
        sphereMap.SetInput(sphere_Source.GetOutput())
            
        #actor coordinates geometry, properties, transformation
        sphere_Actor = vtkLODActor()
        sphere_Actor.SetMapper(sphereMap)
        sphere_Actor.GetProperty().SetColor(atom.getColor())
            
        sphere_Actor.SetPosition(atom.getPosition())
        
        self.ren1.AddActor(sphere_Actor)
        
        self.actors[sphere_Actor] = atom
        
        
        #Add arrow
        if atom.getSpin() != None:
            arrowSource = vtkArrowSource()
    #        arrowSource.SetShaftRadius(arrowSource.GetShaftRadius()/5)
    #        arrowSource.SetTipLength(arrowSource.GetTipLength()/5)
    #        arrowSource.SetTipRadius(arrowSource.GetTipRadius()/5)
    
            aTransform = vtkTransform()
            aTransform.Scale(.2,.2,.2)
            transform = vtkTransformPolyDataFilter()
            transform.SetTransform(aTransform)
            transform.SetInputConnection(arrowSource.GetOutputPort())
            
            arrowMap = vtkPolyDataMapper()
            arrowMap.SetInput(transform.GetOutput())
            
            arrowActor = vtkLODActor()
            arrowActor.SetMapper(arrowMap)
            arrowActor.SetPosition(atom.getPosition())
            arrowActor.GetProperty().SetColor((1,0,0))
            
            
            #Spin Vector
            x = atom.getSpin()[0]
            y = atom.getSpin()[1]
            z = atom.getSpin()[2]
            
            #create unit vetor
            vectLength = (x**2 + y**2 + z**2)**.5
            unitX = x/vectLength
            unitY = y/vectLength
            unitZ = z/vectLength
            
            #Angle
            #cross product of Unit vectors arrow direction and desired orientation
            i, j, k = self.crossProduct(1, 0, 0, unitX, unitY, unitZ)  #default orientation for the arrow is along x axis
            theta = scipy.arcsin((i**2+j**2+k**2)**.5)*180/math.pi
#            print i,j,k
            
            #if the angle is obtuse, theta must be corrected
            if self.dotProduct(1, 0, 0, unitX, unitY, unitZ) >= 0:
                print "acute", atom.getSpin(), theta
                arrowActor.RotateWXYZ(theta, i, j, k)
            else:
                print "obtuse", atom.getSpin(), theta
                arrowActor.RotateWXYZ(180 - theta, i, j, k)
            
            self.ren1.AddActor(arrowActor)
        
        
        
    
    
    
    
    def drawBond(self, bond):
        cylinder = self.makeCylinder(bond.getAtom1().getPosition(), bond.getAtom2().getPosition(), bond.getAtom1().getRadius(), bond.getAtom2().getRadius())
        self.actors[cylinder] = bond
        self.ren1.AddActor(cylinder)

    def makeCylinder (self, posOne, posTwo, rad_One, rad_Two):
        """returns a cylindrical Actor to represent a bond
        
        PosOne and PosTwo are the positions of the two spherical actors
        rad_One and rad_Two are the radii of the spheres"""

        #coordinates of vector pointing from SphereTwo to SphereOne
        x = posTwo[0] - posOne[0]
        y = posTwo[1] - posOne[1]
        z = posTwo[2] - posOne[2]
        distance = ((x**2 + y**2 + z**2)**.5 -rad_One - rad_Two)
        
        #create cylinder
        cylinder = vtkCylinderSource()
        cylinder.SetRadius(self.defualtBondRadius)
        cylinder.SetResolution(100)
        
        #map to another map
        cylinderMap = vtkPolyDataMapper()
        cylinderMap.SetInput(cylinder.GetOutput())
        
        #create cylinder actor
        aCylinder = vtkLODActor()
        aCylinder.SetMapper(cylinderMap)
        aCylinder.GetProperty().SetColor(0, .1, .6)
        aCylinder.SetScale(.2, distance, .2)
        
        #to get the center point between the surfaces of the two spheres
        #create a unit vector and multiply it by the distance/2
        
        #create unit vetor
        vectLength = (x**2 + y**2 + z**2)**.5
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
        i, j, k = self.crossProduct(0, 1, 0, unitX, unitY, unitZ)  #default orientation for the cylinder is along y axis
        theta = scipy.arcsin((i**2+j**2+k**2)**.5)*180/math.pi
        
        #if the angle is obtuse, theta must be corrected
        if self.dotProduct(0, 1, 0, unitX, unitY, unitZ) >= 0:
            aCylinder.RotateWXYZ(theta, i, j, k)
        else:
            aCylinder.RotateWXYZ(theta, -i, -j, -k)
        
        return aCylinder
    
    def dotProduct(self, x1, y1, z1, x2, y2, z2):
        """x,y,z represent vector coordinates - used by makeCylinder"""
        return (x1*x2 + y1*y2 + z1*z2)
    
    def crossProduct(self, x1, y1, z1, x2, y2, z2):
        """x,y,z represent vector coordinates - used by makeCylinder"""
        i = (y1*z2) - (z1*y2)
        j = -(x1*z2) + (z1*x2)
        k = (x1*y2) - (y1*x2)
        return i, j, k
    
    
   
    def drawUnitCell(self, cell):
        #draw very Light Box 
        #Create "Cube" Source
        box = vtkCubeSource()
        box.SetXLength(1)
        box.SetYLength(1)
        box.SetZLength(1)
        
        boxMap = vtkPolyDataMapper()
        boxMap.SetInput(box.GetOutput())
        
        #create actor
        abox = vtkLODActor()
        abox.SetMapper(boxMap)
        abox.GetProperty().SetColor(0, .1, .6)
        abox.GetProperty().SetOpacity(.1)
        a, b, c = cell.getPosition()
        abox.SetPosition(a + .5, b + .5, c + .5)
        abox.PickableOff()
        
        #Add the actor to the renderer 
        self.ren1.AddActor(abox)
        
        #Add the Atoms to the renderer
        for atom in cell.getAtoms():
            self.drawAtom(atom)
            
            
    
    def labelAtoms(self, magneticCell):
        """This should only be called after the image has been renderered, otherwise strange camera effects have been cuased"""
        for cell in magneticCell.getAllUnitCells():
            AtomList = cell.getAtoms()
            for index in range(0, len(AtomList)):
                atom = AtomList[index]
                label = vtkVectorText()
                label.SetText(str(index + 1))
                labelMapper = vtkPolyDataMapper()
                labelMapper.SetInputConnection(label.GetOutputPort())
                labelActor = vtkFollower()
                labelActor.SetMapper(labelMapper)
                labelActor.SetScale(0.05, 0.05, 0.05)
                x, y, z = atom.getPosition()
                x += atom.getRadius()  #display the label on the +x side of sphere
                labelActor.AddPosition(x, y, z)
                labelActor.GetProperty().SetColor(0, 0, 0)
    
                self.ren1.AddActor(labelActor)
                labelActor.SetCamera(self.ren1.GetActiveCamera())
        
    
    def drawMagneticCell(self, MagCell):
        for cell in MagCell.getAllUnitCells():
            self.drawUnitCell(cell)
            #Draw intercellular bonds
            for bond in MagCell.getBonds():
                self.drawBond(bond)

    #Replacing magnetic cell
    def drawCutoffCell(self, cutoffCell):
        for cell in cutoffCell.getAllUnitCells():
            self.drawUnitCell(cell)
            for bond in cutoffCell.getBonds():
                self.drawBond(bond)
            
    def addAxes(self):
        """needs to be run after the window is rendered"""
        
        #Add Axes
        axes = vtkAxes()
        axes.SetOrigin(0, 0, 0)
        axesMapper = vtkPolyDataMapper()
        axesMapper.SetInputConnection(axes.GetOutputPort())
        axesActor = vtkLODActor()
        axesActor.PickableOff()
        axesActor.SetMapper(axesMapper)
        self.ren1.AddActor(axesActor)
        xLabel = vtkVectorText()
        yLabel = vtkVectorText()
        zLabel = vtkVectorText()
        xLabel.SetText("x")
        yLabel.SetText("y")
        zLabel.SetText("z")
        xLabelMapper = vtkPolyDataMapper()
        yLabelMapper = vtkPolyDataMapper()
        zLabelMapper = vtkPolyDataMapper()
        xLabelMapper.SetInputConnection(xLabel.GetOutputPort())
        yLabelMapper.SetInputConnection(yLabel.GetOutputPort())
        zLabelMapper.SetInputConnection(zLabel.GetOutputPort())
        xLabelActor = vtkFollower()
        yLabelActor = vtkFollower()
        zLabelActor = vtkFollower()
        xLabelActor.SetMapper(xLabelMapper)
        yLabelActor.SetMapper(yLabelMapper)
        zLabelActor.SetMapper(zLabelMapper)
        xLabelActor.SetScale(0.1, 0.1, 0.1)
        yLabelActor.SetScale(0.1, 0.1, 0.1)
        zLabelActor.SetScale(0.1, 0.1, 0.1)
        xLabelActor.AddPosition(1, 0, 0)
        yLabelActor.AddPosition(0, 1, 0)
        zLabelActor.AddPosition(0, 0, 1)
        xLabelActor.GetProperty().SetColor(0, 0, 0)
        yLabelActor.GetProperty().SetColor(0, 0, 0)
        zLabelActor.GetProperty().SetColor(0, 0, 0)
        self.ren1.AddActor(xLabelActor)
        self.ren1.AddActor(yLabelActor)
        self.ren1.AddActor(zLabelActor)
        
        #These must be called after ren1 has a camera, otherwise it makes one that is not satisfactory
        #Without them though, 'x','y','z' labels will not be visible from some angles
        xLabelActor.SetCamera(self.ren1.GetActiveCamera())
        yLabelActor.SetCamera(self.ren1.GetActiveCamera())
        zLabelActor.SetCamera(self.ren1.GetActiveCamera())
