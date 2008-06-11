from vtk import *
from AtomClass import *
from BondClass import *




atoms = []
bonds = []

# a renderer for the data
ren1 = vtkRenderer()
ren1.SetBackground(1,1,1)

# a render window to display the contents
renWin = vtkRenderWindow()
renWin.AddRenderer(ren1)
renWin.SetSize(400,400)

# an interactor to allow control of the objects
iren = vtkRenderWindowInteractor()
iren.SetRenderWindow(renWin)






choice = int(raw_input("""Enter Option Number:
1) Add Atom
2) Add Bond
3) Quit"""))
print choice <3
while(choice < 3):
    if(choice == 1):
        atom = Atom(float(raw_input("x coordinate:")), float(raw_input("y coordinate:")), float(raw_input("z coordinate:")))
        ren1.AddActor(atom.getActor())
        print atom.getActor()
        atoms.append(atom)
        
    if(choice == 2):
        bond = Bond(atoms[int(raw_input("first atom number:" + atoms.__str__()))], atoms[int(raw_input("second atom number" + atoms.__str__()))])
        ren1.AddActor(bond.getActor())
        bonds.append(bond)
    
    # trigger the rendering and start the interaction
#    renWin.Render()
    
    
    choice = int(raw_input("""Enter Option Number:
1) Add Atom
2) Add Bond
3) Quit"""))




renWin.Render()
iren.Start()




#iren.Initialize()
#don't really know what this does or if it is necessary
#add my observer for picking
#picker = iren.GetPicker()
#def pickEvent(obj, event):
#    if(picker.GetActor()==aSphere):
#        print "aSphere"
#    if(picker.GetActor()==aSphere_Two):
 #       print "aSphere_Two"
#    print picker.GetActor()

#picker.AddObserver("PickEvent", pickEvent)
#iren.Enable()


