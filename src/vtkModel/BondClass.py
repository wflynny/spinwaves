import AtomClass
#from vtk import *
import numpy

class JParam():
    """This class represents one value in a 3*3 J matrix.  In the simplest case it would
    contain a single float value, but it can also be a variable tied to other variables
    for fitting purposes."""
    def __init__(self, manager, fit = False, value = 0., min = '-inf', max = '+inf'):
        """-fit is a boolean value specifying whether this parameter is variable(True)
        or a fixed float(False)
        -manager is the instance of ParamManager that is to contain this JParam
        -value is the float value if there is one(if fit is False)
        -min is a string of the float number representing the minimum possible value
        (only used if fit is TRUE) '-inf' will be used for negative infinite.
        -max is a string of the float number representing the maximim possible value
        (only used if fit is TRUE) '+inf' will be used for positive infinite."""
        self.manager = manager
        self.manager.addParam(self)
        
        
        self.fit = fit
        self.value = value
        self.min = min
        self.max = max
        self.tied = []#only the manager should touch this
        
    def tieTo(self, index):
        """Ties this parameter to the JParam in the manager at the given index."""
        self.manager.tie(self, index)
        
    def untie(self, index):
        """Will remove the tie between this parameter and the one given by index."""
        self.manager.untie(self, index)
    
    
    def tieToMany(self, list):
        """will tie this parameter to all the parameters in the given by the list of indices.
        Any ties that currently exist that are not in the list will be removed."""
        for p in self.tied:#untie all
            self.untie(p)
        for p in list:#tie to all in the new list
            self.tieTo(p)
    
    def isDefault(self):
        """Returns true if the values are the in this parameter are the defualt values, or
        if False if they have been set to something else."""
        if not self.fit:
            if self.value == 0.0:
                return True
        return False
    
    def getName(self):
        return "p" + str(self.manager.getIndex(self))
    
    def isTiedTo(self, index):
        """returns True if this parameter is tied to the parameter given by index."""
        try:#check if the list already contains the index
            self.tied.index(index)
            return True
        except:
            return False
    
    def __str__(self):
        """This will return the value if fit==False, or the range of values if fit==True"""
        if not self.fit:
            return str(self.value)
        return self.min + " - " + self.max



class Bond():
    """This class represents interactions, or "bonds" between atoms."""
    
    def __init__(self, Atom1, Atom2, jMatrix = None, r = 0,g = 0,b = 1):
        """Cell is the Unit Cell if this bond is in a unit Cell or None Otherwise.
        Atom1 and Atom2 are the atoms that this bond connects.
        r,g,b are the color of the actor. jMatrix is a numpy 2D array, such as
        [[1, 0, 0],
         [0, 1, 0],
         [0, 0, 1]]   for a ferromagnetic interaction.  The values in the jMatrix
         will be JParam objects, so that the values can be variable."""

        #Color
        self.r = r
        self.g = g
        self.b = b
        
        self.Atom1 = Atom1
        self.Atom2 = Atom2
        
        self.jMat = jMatrix #a numpy 2D array
        print "\nnew Bond JMat =\n",jMatrix
        
        #New attributes added for fitting purposes

    
    def getAtom1(self):
        return self.Atom1
    
    def getAtom2(self):
        return self.Atom2
    
    def setAtom1(self, atom):
        self.Atom1 = atom
    
    def setAtom2(self, atom):
        self.Atom2 = atom

    def sameBond(self, otherBond):
        """returns true if otherBond connects the same 2 atoms and false otherwise"""
        if self.getAtom1() == otherBond.getAtom1() or self.getAtom1() == otherBond.getAtom2():
            if self.getAtom2() == otherBond.getAtom1() or self.getAtom2() == otherBond.getAtom2():
                return True
        return False
    
    def __str__(self):
        str = "Bond between " + self.Atom1.__str__() + " and " + self.Atom2.__str__()
        if self.jMat != None:
            str += ":\n" + self.jMat.__str__()
        return str
    
    def getRGBColor(self):
        return self.r, self.g, self.b
    
    def getJMatrix(self):
        return self.jMat

