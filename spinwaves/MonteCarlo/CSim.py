import os
import ctypes
from ctypes import c_float, c_int, c_long
import sys
import time

import numpy as N
import wx
#import pylab

from simple import readFile, Timer, simpleAtom

#dllpath=r'C:\mytripleaxisproject\trunk\eclipse\src\spinwaves\C code'
#dllpath=r'C:\tom winter\code\C code'#NIST
#linux_so_path = r'/home/tom/Desktop/workCode/C code'#home laptop
#dllpath=r'C:\Users\Tom\Documents\spinwaves\C code\monteCarlo.dll'

if sys.platform in ('darwin'):
    ext = '.dylib'
elif sys.platform in ('win32','cygwin'):
    ext = '.pyd'
else:
    ext = '.so'
dllpath=os.path.join(os.path.dirname(__file__),'_monteCarlo'+ext)
monteCarloDll = ctypes.cdll[dllpath]



def createVideo(spinsToImageFunction, outFilePath, inFilePath):
    """This method is used to create snapshots of the monte carlo simulation.

    This is almost the same as simulate, but the outer loops
    are controlled in python so that snapshots can be taken"""
    
    timer = Timer()
    
    #C code and dll should be put into a folder in close to this code

    if sys.platform=='win32':
        print 'win32'
        monteCarloDll = N.ctypeslib.load_library('monteCarlo', 'C:/spinwaves')
#    elif sys.platform=='mac':
#        monteCarloDll = N.ctypeslib.load_library('libpolarization2.so', '.')
#    else:
#        monteCarloDll = N.ctypeslib.load_library('libpolarization2.so', '.') #linux
        
    atoms, jMatrices = readFile(inFilePath)
    

    successCode = c_int(0)
    atomListPointer = monteCarloDll.new_atom_list(c_int(len(atoms)), ctypes.byref(successCode))
#    print "Success Code: ", successCode
    if successCode.value != 1:
        raise Exception("Problem Allocating memory, simulation size may be too large")

    #to keep pointers
    matListList = []
    nbr_ListList = []

    #Create atom list in C
    s1Class = c_float * 3
    for i in range(len(atoms)):
#    for i in range(5):
        atom = atoms[i]
        numInteractions = len(atom.interactions)
        matListClass = c_int * numInteractions
        nbrListClass = c_int * numInteractions
        
        matList = matListClass()
        neighbors = nbrListClass()
        for j in range(numInteractions):
            neighbors[j] = c_int(atom.interactions[j][0])
            matList[j] = c_int(atom.interactions[j][1])
#            print "interaction: " , atom.interactions[j][0], atom.interactions[j][1],  " = " , neighbors[j], matList[j]
        
        s1 = s1Class()
        s1[0] = c_float(1)
        s1[1] = c_float(0)
        s1[2] = c_float(0)
        
        anisotropyClass = c_float * 3
        anisotropy = anisotropyClass()
        anisotropy[0] = c_float(atom.anisotropy[0])
        anisotropy[1] = c_float(atom.anisotropy[1])
        anisotropy[2] = c_float(atom.anisotropy[2])
        
        #to keep pointers
        matListList.append(matList)
        nbr_ListList.append(neighbors)

        monteCarloDll.set_atom(atomListPointer, c_int(i), anisotropy, matList, neighbors, c_int(numInteractions), s1)
    
    print "atoms added"
    timer.printTime()

    #Create JMatrix List in C
    matPointer = monteCarloDll.new_jMatrix_list(c_int(len(jMatrices)), ctypes.byref(successCode))
    if successCode.value != 1:
        raise Exception("Problem Allocating memory, simulation size may be too large")


    for i in range(len(jMatrices)):
        j = jMatrices[i]
        j11 = c_float(j[0][0])
        j12 = c_float(j[0][1])
        j13 = c_float(j[0][2])
        j21 = c_float(j[1][0])
        j22 = c_float(j[1][1])
        j23 = c_float(j[1][2])
        j31 = c_float(j[2][0])
        j32 = c_float(j[2][1])
        j33 = c_float(j[2][2])
        monteCarloDll.add_matrix(matPointer, c_int(i),j11, j12, j13, j21, j22, j23, j31, j32, j33)

    
    #Until this point, code is almost the same as simulate(), now flipspins() must be called directly from python
    totalMagX = 0
    totalMagY = 0
    totalMagZ = 0
    magnetizations = []#to make graph after
    temperatures = []
    T = 10#tMax
    imageNum = 0
    spins = None
    while T > .005:#tMin
        for i in range(10):
            for j in range(10):
                monteCarloDll.flipSpins(atomListPointer, c_int(len(atoms)), matPointer, c_float(T), ctypes.byref(c_int(0)))#last parameter not used\
            #output spins to file
            spins = []
            for i in range(len(atoms)):
                spin = s1Class()
                monteCarloDll.getSpin(atomListPointer, c_int(i), ctypes.byref(spin))
                spins.append(spin)

            print "writing spins to file..."
            timer.printTime()
            
            #output the spins to a file
            outFile = open(outFilePath, 'w')
            outFile.write("#Atom_Number Position_X Position_Y Position_Z Spin_X Spin_Y Spin_Z\n")
            for i in range(len(atoms)):
                atom = atoms[i]
                Posx = str(atom.pos[0])
                Posy = str(atom.pos[1])
                Posz = str(atom.pos[2])
                spin = spins[i]
                Spinx = str(spin[0])
                Spiny = str(spin[1])
                Spinz = str(spin[2])
                atomStr = str(i) + " " + Posx + " " + Posy + " " + Posz + " " + Spinx + " " + Spiny + " " + Spinz + "\n"
                outFile.write(atomStr)
            
            outFile.close()
            #draw the spins
            spinsToImageFunction(outFilePath, imageNum)
            imageNum += 1
        
        for spin in spins:
            totalMagX += spin[0]
            totalMagY += spin[1]
            totalMagZ += spin[2]
            
        magnetizations.append( (totalMagX**2 + totalMagY**2 + totalMagZ**2)**(.5) )
        print "magnetization:", magnetizations[len(magnetizations)-1]  
        temperatures.append(T)
        totalMagX = 0
        totalMagY = 0
        totalMagZ = 0
        T = T*.9#tFactor
                                


    monteCarloDll.del_jMat(matPointer)
    monteCarloDll.del_atom(atomListPointer)
    
    
    #drawing magnetization graphs
    #find max value
    max = 0
    min = 0
    for num in magnetizations:
        if num > max:
            max = num
        elif num < min:
            min = num
    
    #draw magnetization vs. temp
#    pylab.plot(temperatures, magnetizations, 'ro')
#    pylab.axis([temperatures[len(temperatures)-1], temperatures[0], min, max])
#    pylab.show()
    
    #draw magnetization vs. iteration number
#    pylab.plot(range(len(magnetizations)), magnetizations, 'ro')
#    pylab.axis([len(temperatures)-1, 0, min, max])
    #    savefig('secondfig.png')
#    pylab.show()
    
    
    timer.printTime()
    print "done"


def simulate(k, tMax, tMin, tFactor, inFilePath, outFilePath):
    """Runs the monte carlo simulation written in C.

    k is the number of steps per temperature.
    tFactor is the factor which is multiplied by the temperature to decrease it
    at each temperature step.
    i.e. .9 would be a good factor."""
    timer = Timer()
    

    atoms, jMatrices = readFile(inFilePath)

    

    #Create the atom list in C
    successCode = c_int(0)
    atomListPointer = monteCarloDll.new_atom_list(c_int(len(atoms)), ctypes.byref(successCode))
#    print "Success Code: ", successCode
    if successCode.value != 1:
        raise Exception("Problem Allocating memory, simulation size may be too large")

    #to keep pointers
    matListList = []
    nbr_ListList = []

    #Populate atom list in C
    s1Class = c_float * 3
    for i in range(len(atoms)):
        atom = atoms[i]
        numInteractions = len(atom.interactions)
        matListClass = c_int * numInteractions
        nbrListClass = c_int * numInteractions
        
        matList = matListClass()
        neighbors = nbrListClass()
        for j in range(numInteractions):
            neighbors[j] = c_int(atom.interactions[j][0])
            matList[j] = c_int(atom.interactions[j][1])
   
        s1 = s1Class()
        s1[0] = c_float(1)
        s1[1] = c_float(0)
        s1[2] = c_float(0)
        
        anisotropyClass = c_float * 3
        anisotropy = anisotropyClass()
        anisotropy[0] = c_float(atom.anisotropy[0])
        anisotropy[1] = c_float(atom.anisotropy[1])
        anisotropy[2] = c_float(atom.anisotropy[2])
        
        #to keep pointers
        matListList.append(matList)
        nbr_ListList.append(neighbors)

        monteCarloDll.set_atom(atomListPointer, c_int(i), anisotropy, matList, neighbors, c_int(numInteractions), s1)
    
    print "atoms added"
    timer.printTime()

    #Create JMatrix List in C
    matPointer = monteCarloDll.new_jMatrix_list(c_int(len(jMatrices)), ctypes.byref(successCode))
#    print "Matrix Success Code: ", successCode
#   time.sleep(5)
    if successCode.value != 1:
        raise Exception("Problem Allocating memory, simulation size may be too large")

    for i in range(len(jMatrices)):
        j = jMatrices[i]
        j11 = c_float(j[0][0])
        j12 = c_float(j[0][1])
        j13 = c_float(j[0][2])
        j21 = c_float(j[1][0])
        j22 = c_float(j[1][1])
        j23 = c_float(j[1][2])
        j31 = c_float(j[2][0])
        j32 = c_float(j[2][1])
        j33 = c_float(j[2][2])
        monteCarloDll.add_matrix(matPointer, c_int(i),j11, j12, j13, j21, j22, j23, j31, j32, j33)
    

    #Run the simulation
    monteCarloDll.simulate(atomListPointer, c_int(len(atoms)), matPointer, c_int(k), c_float(tMax), c_float(tMin), c_float(tFactor))

    
    #Get the spins from the C code
    spins = []
    for i in range(len(atoms)):
        spin = s1Class()
        monteCarloDll.getSpin(atomListPointer, c_int(i), ctypes.byref(spin))
        spins.append(spin)
    
   
    #Free the memory used by the C code
    monteCarloDll.del_jMat(matPointer)
    monteCarloDll.del_atom(atomListPointer)
    
    
    print "writing spins to file..."
    timer.printTime()
    
    #output the spins to a file
    outFile = open(outFilePath, 'w')
    outFile.write("#Atom_Number Position_X Position_Y Position_Z Spin_X Spin_Y Spin_Z\n")
    for i in range(len(atoms)):
        atom = atoms[i]
        Posx = str(atom.pos[0])
        Posy = str(atom.pos[1])
        Posz = str(atom.pos[2])
        spin = spins[i]
        Spinx = str(spin[0])
        Spiny = str(spin[1])
        Spinz = str(spin[2])
        atomStr = str(i) + " " + Posx + " " + Posy + " " + Posz + " " + Spinx + " " + Spiny + " " + Spinz + "\n"
        outFile.write(atomStr)
    
    outFile.close()
    
    
    
    timer.printTime()
    print "done"



class MonteCarloPanel(wx.Panel):
    """This is a simple GUI for the monte Carlo simulation."""
    def __init__(self, parent, id, defaultSteps, defaultTMax, defaultTMin, defaultTFactor, defaultInPath, defaultOutPath):
        wx.Panel.__init__(self, parent, id)
        
        #Create the Flex Grid Sizer
        topSizer = wx.FlexGridSizer(4,2,2,2)
        
        #Add Temperacure Max text
        tMaxLabel = wx.StaticText(self, -1, "Max Temperature:")
        self.tMaxText = wx.TextCtrl(self, -1, value = str(defaultTMax), size = (60, -1), style = wx.TE_RICH2)
        topSizer.Add(tMaxLabel)
        topSizer.Add(self.tMaxText)

        #Add Temperacure Min text
        tMinLabel = wx.StaticText(self, -1, "Min Temperature:")
        self.tMinText = wx.TextCtrl(self, -1, value = str(defaultTMin), size = (60, -1), style = wx.TE_RICH2)
        topSizer.Add(tMinLabel)
        topSizer.Add(self.tMinText)
        
        #Add Temperature Factor Text
        tFactorLabel = wx.StaticText(self, -1, "Temperature Factor:")
        self.tFactorText = wx.TextCtrl(self, -1, value = str(defaultTFactor), size = (60, -1), style = wx.TE_RICH2)
        topSizer.Add(tFactorLabel, 0)
        topSizer.Add(self.tFactorText, 0)
        
        #Add Steps per Temperature Level
        stepsLabel = wx.StaticText(self, -1, "Steps Per Temperacture Level:")
        self.stepsText = wx.TextCtrl(self, -1, value = str(defaultSteps), size = (60, -1), style = wx.TE_RICH2)
        topSizer.Add(stepsLabel, 0)
        topSizer.Add(self.stepsText, 0)
        
        
        bottomSizer = wx.FlexGridSizer(2,3,2,2)
        
        #Add Bond File Path
        inPathLabel = wx.StaticText(self, -1, "Bond File Path:")
        self.inPathText = wx.TextCtrl(self, -1, value= defaultInPath, size = (120, -1), style = wx.TE_RICH2)
        bottomSizer.Add(inPathLabel, 0)
        bottomSizer.Add(self.inPathText, 0)
        self.inPathBrowseButton = wx.Button(self, -1, "Browse")
        bottomSizer.Add(self.inPathBrowseButton)
        
        #Add Spin Output File Path
        outPathLabel = wx.StaticText(self, -1, "Spin Output Path:")
        self.outPathText = wx.TextCtrl(self, -1, value = defaultOutPath, size = (120, -1), style = wx.TE_RICH2)
        bottomSizer.Add(outPathLabel, 0)
        bottomSizer.Add(self.outPathText, 0)
        self.outPathBrowseButton = wx.Button(self, -1, "Browse")
        bottomSizer.Add(self.outPathBrowseButton)
        
        
        self.runButton = wx.Button(self, -1, "Run Simulation")
        
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(topSizer, 0)
        sizer.Add((1,20), 0, wx.EXPAND)#Spacer between filenames and simulation parameters
        sizer.Add(bottomSizer, 0)
        sizer.Add((1,15), 0, wx.EXPAND)#Spacer before run button
        sizer.Add(self.runButton,0, wx.ALIGN_CENTER)
        
        self.SetSizer(sizer)
        
        #Add event handlers for the buttons
        self.Bind(wx.EVT_BUTTON, self.OnInFileBrowse, self.inPathBrowseButton)
        self.Bind(wx.EVT_BUTTON, self.OnOutFileBrowse, self.outPathBrowseButton)
        self.Bind(wx.EVT_BUTTON, self.OnRun, self.runButton)
        
    def OnInFileBrowse(self, event):
        openDialog = wx.FileDialog(self, "Open File", style = wx.OPEN, wildcard = "Bond Export File(*.txt)|*.txt")
        if openDialog.ShowModal() == wx.ID_OK:
            self.inPathText.SetValue(openDialog.GetPath())
        openDialog.Destroy()
        event.Skip()
        
    def OnOutFileBrowse(self, event):
        saveDialog = wx.FileDialog(self, "Save File", style = wx.SAVE, wildcard = "Spin File(*.txt)|*.txt")
        if saveDialog.ShowModal() == wx.ID_OK:
            self.outPathText.SetValue(saveDialog.GetPath())
        saveDialog.Destroy()
        event.Skip()
        
    def OnRun(self, event):
        event.Skip()
        failed, tMax, tMin, tFactor, steps = self.validate()
        if not failed:
#            print steps, tMax, tMin, tFactor
#            time.sleep(10)
            try:
                simulate(steps, tMax, tMin, tFactor, self.inPathText.GetValue(), self.outPathText.GetValue())
                wx.MessageDialog(self, "Simulation Done!", caption = "Complete", style = wx.OK).ShowModal()
            except IOError:
                wx.MessageDialog(None, "Bad File Name!", "File Name Error", wx.OK).ShowModal()
        
    def validate(self):
        """Validates the info entered for tMin, tMax, tFactor, and k to make sure
        they are of hte right types."""
        tMax = self.tMaxText.GetValue()
        tMin = self.tMinText.GetValue()
        tFactor = self.tFactorText.GetValue()
        steps = self.stepsText.GetValue()
        #For now, I will not validate the file path
        #inPath = self.inPathtext.GetValue()#Can check if it is readable
        
        bgColor = "pink"
        failed = False
        #Validate tMax(must be a float)
        numTmax = None
        try:
            numTmax = float(tMax)
            self.tMaxText.SetStyle(0, len(tMax), wx.TextAttr(colBack = "white"))
        except:
            self.tMaxText.SetStyle(0, len(tMax), wx.TextAttr(colBack = bgColor))
            failed = True
        
        #Validate tMin(must be a float)
        numTmin = None
        try:
            numTmin = float(tMin)
            self.tMinText.SetStyle(0, len(tMin), wx.TextAttr(colBack = "white"))
        except:
            self.tMinText.SetStyle(0, len(tMin), wx.TextAttr(colBack = bgColor))
            failed = True
        
        #Validate tFactor(must be a float)
        numTfactor = None
        try:
            numTfactor = float(tFactor)
            self.tFactorText.SetStyle(0, len(tFactor), wx.TextAttr(colBack = "white"))
        except:
            self.tFactorText.SetStyle(0, len(tFactor), wx.TextAttr(colBack = bgColor))
            failed = True
            
        #Validate steps(must be an int)
        numSteps = None
        try:
            numSteps = int(steps)
            self.stepsText.SetStyle(0, len(steps), wx.TextAttr(colBack = "white"))
        except:
            self.stepsText.SetStyle(0, len(steps), wx.TextAttr(colBack = bgColor))
            failed = True
            
        return failed, numTmax, numTmin, numTfactor, numSteps 
 
         
def ShowSimulationFrame():
    """Creates and displays a simple frame containing the MonteCarloPanel."""
    #Defaults
    k = 100
    tMax = 10
    tMin = .01
    tFactor = .90
    inFilePath = "C:\montecarlo.txt"
    outFilePath = "C:\Spins.txt"
        
    
    frame = wx.Frame(None, -1, title = "Monte Carlo Simulation", size = (300,250))
    MonteCarloPanel(frame, -1, k, tMax, tMin, tFactor, inFilePath, outFilePath)
    frame.Show()
    return frame
            


 


class App(wx.App):
    def __init__(self, redirect = False, filename = None):
        wx.App.__init__(self, redirect, filename)
    
    def OnInit(self):
        self.SetTopWindow(ShowSimulationFrame())
        return True


if __name__ == '__main__':       
    app = App(False)
    app.MainLoop()
    
