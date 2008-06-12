from vtk import *
from AtomClass import *
from BondClass import *
from SymmetryUtilities import *
from SpaceGroups import *
import random
from CellClass import *
from MagneticCellClass import *


def generateAtoms(position, description, radius, r,g,b):
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
                    for bond in MagCell.getIntercellularBonds():
                        if bond.getActor() == SelectedActor:
                            print bond
                            break






if __name__=='__main__':
    
    
    unitcell = Cell()
    Space_Group = sg65

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
     
    
    interactor = vtkInteractorStyleSwitch()
    interactor.SetCurrentStyleToTrackballCamera()
    iren.SetInteractorStyle(interactor)
    picker = vtkPropPicker()
    iren.SetPicker(picker)
    
    
    SelectedActor = None  #used by the picker
    iren.AddObserver("LeftButtonPressEvent", pick)
    
    
    renWin.Render()
    iren.Start()