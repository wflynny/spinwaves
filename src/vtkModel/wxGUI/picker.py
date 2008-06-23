import vtk

class Picker():
    def __init__(self, vtkDrawer, iren, renderer):
        self.iren = iren
        self.ren1 = renderer
        self.drawer = vtkDrawer
        #Set up trackball mode (not really picking)
        interactor = vtk.vtkInteractorStyleSwitch()
        interactor.SetCurrentStyleToTrackballCamera()
        iren.SetInteractorStyle(interactor)
        
        #set the picker so props can be picked
        self.picker = vtk.vtkPropPicker()
        iren.SetPicker(self.picker)
        
        #Add my own pick function
        iren.AddObserver("LeftButtonPressEvent", self.pick)
        
        self.SelectedActor = None  #used by the pick()
    
    def pick(self, obj, event):
        Mouse_Position = self.iren.GetEventPosition()
        self.picker.PickProp(Mouse_Position[0],Mouse_Position[1], self.ren1)
        if(self.SelectedActor == self.picker.GetActor()): #the actor is already picked
            return
        if(self.SelectedActor != None):
            self.SelectedActor.GetProperty().SetAmbient(0) #unhighlight old picked actor
        self.SelectedActor = self.picker.GetActor()
        if self.SelectedActor != None:
            self.SelectedActor.GetProperty().SetAmbient(1)#make the selected actor stand out
            self.iren.GetRenderWindow().Render()
            #find the Atom at this position and print its description
            print self.drawer.getObjFromActor(self.SelectedActor)
    
    def getPicked(self):
        return self.SelectedActor