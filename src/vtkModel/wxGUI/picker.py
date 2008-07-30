import vtk
from wx.py.dispatcher import send

class Picker():
    
    """This class is in charge of listening to pick events in the render window in
    the and handling them.  The actor is highlighted by turning the ambient lighting
    on.  Then the object (Atom or Bond; others are not supported) is found using
    the vtkDrawer.  The object's __str__ is printed"""

    def __init__(self, vtkDrawer, iren, renderer):
        self.iren = iren
        self.ren1 = renderer
        self.drawer = vtkDrawer
        
        #moved this to draw method of main panel
        #Set up trackball mode (not really picking)
#        interactor = vtk.vtkInteractorStyleSwitch()
#        interactor.SetCurrentStyleToTrackballCamera()
#        iren.SetInteractorStyle(interactor)
        
        #set the picker so props can be picked
        self.picker = vtk.vtkPropPicker()
        iren.SetPicker(self.picker)
        
        #Add my own pick function
        self.observerNum = iren.AddObserver("LeftButtonPressEvent", self.pick)
        
        self.SelectedActor = None   #used by the pick() so that only one item is
                                    #picked at a time and it is not repicked
    
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
            modelObj = self.drawer.getObjFromActor(self.SelectedActor)
#            print "sending signal pick..."
            send(signal = "Pick Event", sender = "Picker", obj = modelObj)
            actors = self.ren1.GetActors()
#            for i in range(actors.GetNumberOfItems()):
#                print actors.GetNextActor()
            print modelObj
#            print self.SelectedActor
    
    def removeObserver(self):
        """removes this picker from the render window
        If this is not called and another picker is added, both will be active"""
        
        self.iren.RemoveObserver(self.observerNum)
    
    def getPicked(self):
        return self.SelectedActor