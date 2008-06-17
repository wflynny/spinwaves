from vtk import *
from AtomClass import *
from BondClass import *
from SymmetryUtilities import *
from SpaceGroups import *
import random
from CellClass import *
from MagneticCellClass import *


def generateAtoms(Space_Group, unitcell, position, description, radius, r,g,b):
    locations = expandPosition(Space_Group, numpy.array([position[0],position[1], position[2]]))[0]
    for coord in locations:
        print coord[0], coord[1], coord[2]
#        r,g,b = atom1.getActor().GetProperty().GetColor()
        atom = Atom(unitcell, coord[0], coord[1], coord [2], description, radius, r,g,b)
        unitcell.addAtom(atom)


def pick(obj, Event):
    global SelectedActor
    Mouse_Position = iren.GetEventPosition()
    picker.PickProp(Mouse_Position[0],Mouse_Position[1], ren1)
    if(SelectedActor == picker.GetActor()):
        return
    if(SelectedActor != None):
        SelectedActor.GetProperty().SetAmbient(0)
    SelectedActor = picker.GetActor()
    if SelectedActor != None:
        SelectedActor.GetProperty().SetAmbient(1)
        renWin.Render()
        #find the Atom at this position and print its description
        for celln in MagCell.getAllUnitCells():   
            for atom in celln.getAtoms():
                if atom.getActor() == SelectedActor:
                    print atom
                    break
            else: # the atom is not found so check bonds
                for bond in celln.getBonds():
                    if bond.getActor() == SelectedActor:
                        print bond
                        break
        else: #the unit cell bond is not found so check intercellular bonds
            for bond1 in MagCell.getIntercellularBonds():
                if bond1.getActor() == SelectedActor:
                    print bond1
                    break



def menu():
    #generate Unit Cell
    randGen = random.Random()
    choice = int(raw_input("""Unit Cell:\nEnter Option Number:
    1) Add Atom
    2) Add Bond
    3) Generate Magenetic Cell"""))
    
    
    #Gather Information
    while choice < 3 :
        if choice == 1 :
            generateAtoms([float(raw_input("x coordinate:")), float(raw_input("y coordinate:")), float(raw_input("z coordinate:"))], str(raw_input("description:")) , .05, randGen.uniform(0,1), randGen.uniform(0,1), randGen.uniform(0,1))
        
        if choice == 2 :
            atoms = unitcell.getAtoms()
            i = 0
            for singleatom in atoms:
                print i, ") ", singleatom
                i +=1
            bond = Bond(unitcell, atoms[int(raw_input("first atom number:\n"))], atoms[int(raw_input("second atom number\n"))], randGen.uniform(0,1), randGen.uniform(0,1), randGen.uniform(0,1))
            ren1.AddActor(bond.getActor())
            unitcell.addBond(bond)
            for eachBond in bond.createSymmetryBonds(Space_Group):
                unitcell.addBond(eachBond)
        
        choice = int(raw_input("""Unit Cell:\nEnter Option Number:
    1) Add Atom
    2) Add Bond
    3) Generate Magnetic Cell"""))
    
    na = int(raw_input("N(a):"))
    nb = int(raw_input("N(b):"))
    nc = int(raw_input("N(c):"))
    
    #create the Magnetic Cell and add it to the Renderer
    MagCell = MagneticCell(unitcell, na, nb, nc, Space_Group)
    
    
    #Add Bonds between Cells
    choice = int(raw_input("""Magnetic Cell:
     1) Add Bond Between Cells
     2) Draw"""))
    
    AllAtoms = MagCell.getAllAtoms()
    #Gather Information
    while choice < 2 :
        if choice == 1 :
            for i in range(0,len(AllAtoms)):
                print i, ") " + AllAtoms[i].__str__()
            MagCell.addInterCellularBond(AllAtoms[int(raw_input("First atom number:\n"))], AllAtoms[int(raw_input("Second atom number:\n"))])
            
        
        choice = int(raw_input("""Magnetic Cell:
     1) Add Bond Between Cells
     2) Draw"""))
    
    MagCell.drawCell(ren1)
    for bond in MagCell.getIntercellularBonds():
        print bond 
    
    return MagCell



if __name__=='__main__':
    
    
    unitcell = Cell()
    Space_Group = sg65
    atomPos = [.25, .25, .5]

    # a renderer for the data
    ren1 = vtkRenderer()
    ren1.SetBackground(1,1,1)

    # a render window to display the contents
    renWin = vtkRenderWindow()
    renWin.AddRenderer(ren1)
    renWin.SetSize(700,700)

    # an interactor to allow control of the objects
    iren = vtkRenderWindowInteractor()
    iren.SetRenderWindow(renWin)
    
    
#    MagCell = menu()
    randGen = random.Random()
    generateAtoms(Space_Group, unitcell, atomPos, "atom1" , .05, randGen.uniform(0,1), randGen.uniform(0,1), randGen.uniform(0,1))
#    bond = Bond(unitcell, atoms[0], atoms[5], randGen.uniform(0,1), randGen.uniform(0,1), randGen.uniform(0,1))
#    ren1.AddActor(bond.getActor())
 #   unitcell.addBond(bond)
 #   for eachBond in bond.createSymmetryBonds(Space_Group):
  #      unitcell.addBond(eachBond)
    
    
    MagCell = MagneticCell(unitcell, 1, 2, 3, Space_Group)
    AllAtoms = MagCell.getAllAtoms()
    MagCell.addInterCellularBond(AllAtoms[0], AllAtoms[6])
    MagCell.drawCell(ren1)
    
    
    interactor = vtkInteractorStyleSwitch()
    interactor.SetCurrentStyleToTrackballCamera()
    iren.SetInteractorStyle(interactor)
    picker = vtkPropPicker()
    iren.SetPicker(picker)
    
    #Add Axes
    axes = vtkAxes()
    axes.SetOrigin(0,0,0)
    axesMapper = vtkPolyDataMapper()
    axesMapper.SetInputConnection(axes.GetOutputPort())
    axesActor = vtkActor()
    axesActor.SetMapper(axesMapper)
    ren1.AddActor(axesActor)
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
    xLabelActor.SetScale(0.1,0.1,0.1)
    yLabelActor.SetScale(0.1,0.1,0.1)
    zLabelActor.SetScale(0.1,0.1,0.1)
    xLabelActor.AddPosition(1,0,0)
    yLabelActor.AddPosition(0,1,0)
    zLabelActor.AddPosition(0,0,1)
    xLabelActor.GetProperty().SetColor(0,0,0)
    yLabelActor.GetProperty().SetColor(0,0,0)
    zLabelActor.GetProperty().SetColor(0,0,0)
    ren1.AddActor(xLabelActor)
    ren1.AddActor(yLabelActor)
    ren1.AddActor(zLabelActor)
    
    
    SelectedActor = None  #used by the picker
    iren.AddObserver("LeftButtonPressEvent", pick)
    
    renWin.Render()    
    xLabelActor.SetCamera(ren1.GetActiveCamera())
    yLabelActor.SetCamera(ren1.GetActiveCamera())
    zLabelActor.SetCamera(ren1.GetActiveCamera())
    iren.Start()