import sys
import sympy
import numpy
import matplotlib
matplotlib.use('WXAgg')
import pylab
from sympy import pi
from numpy import PINF, NINF
from spinwaves.MonteCarlo.simple import simpleAtom
from spinwaves.MonteCarlo.CSim import Sim_Aux
from spinwaves.spinwavecalc.readfiles import findmat, atom as SpinwaveAtom
from spinwaves.spinwavecalc.spinwave_calc_file import calculate_dispersion, calc_eigs

class ParamManager():
    def __init__(self):
        self.parameters = []
        #each group of tied parameters will be given a group number
        self.groups = []#list of tuples (num(=index), bool(active/inactive))
    
    def AssignNewGroup(self, param):
        for group in self.groups:
            if not group[1]:
                param.group = group[0]
                group[1] = True
                return
        num = len(self.groups)
        self.groups.append( [num,True] )
        param.group = num
        #print "assigned ", param.group
    
    def addParam(self, param):
        """Appends a new parameter to this list."""
        self.parameters.append(param)
        
    def validIndex(self, index):
        """Checks if the integer index is a valid index(identifier) of a parameter."""
        return (index < len(self.parameters))
    
    def tie(self, paramObj, param2_index):
        """Adds param2_index to the list of tied indices on paramObj, and vice versa.
        Also changes the min, max, value, and fit(bool) values in the parameter given
        by param2_index to equal those of paramObj."""
        index1 = self.parameters.index(paramObj)
        param1 = paramObj
        index2 = param2_index
        param2 = self.parameters[index2]
        if index1 != index2:
            #Add second index to first param
            if not param1.isTiedTo(index2):
                param1.tied.append(index2)
            #Add first index to second param
            if not param2.isTiedTo(index1):
                param2.tied.append(index1)
            
        #Change values on second parameter to equal those of paramObj
        param2.fit = param1.fit
        param2.value = param1.value
        param2.min = param1.min
        param2.max = param1.max
        param2.group = param1.group
        param2.default = param1.default
        #param1.fit = param2.fit
        #param1.value = param2.value
        #param1.min = param2.min
        #param1.max = param2.max
        #param1.group = param2.group
        
        self.__updateGroups()
        #tie to any parameters tied to the second parameter
        for eachIndex in param2.tied:
            if not param1.isTiedTo(eachIndex):
                self.tie(param1, eachIndex)
        
        for eachIndex in param1.tied:
            if not param2.isTiedTo(eachIndex):
                self.tie(param2, eachIndex)
                
        #to look nicer for the user
        param1.tied.sort()
        param2.tied.sort()
    
    def __updateGroups(self):
        groups = []
        for param in self.parameters:
            for group in groups:
                if group == param.group:
                    break
            else:#Creating a list of unique groups that are represented
                groups.append(param.group)
        
        for i in range(len(self.groups)):
            for rep_group in groups:
                if rep_group == i:
                    break
            else:
                self.groups[i][1] = False
                
        #shuffle groups down to eliminate empty group numbers
        #for i in range(len(self.groups)):
        #    if not self.groups[i][1]:
        #        for param in self.parameters:
        #            if param.group > i:
        #                param.group -= 1
        #    self.groups[i][1] = True
        
    
    def untie(self, paramObj, index):
        """will untie the given Jparam object with the parameter given by index as well as
        all parameters tied to the parameter at index."""
        self.AssignNewGroup(paramObj)
        if paramObj.isTiedTo(index):
            paramObj.tied.remove(index)
            self.parameters[index].tied.remove(self.parameters.index(paramObj))
        
            for p in self.parameters[index].tied:
                self.untie(paramObj, p)
        self.__updateGroups()
    
    def removeParam(self, param):
        """removes the given JParam object from the list of parameters and corrects all
        indices and tied lists."""
        
        index = self.parameters.index(param)
        self.parameters.pop(index)
        for parameter in self.parameters:
            for tiedIndex in parameter.tied:
                if tiedIndex>=index:
                    tiedIndex-=1
        self.__updateGroups()
                    
    def getIndex(self, paramObj):
        """Returns the index of the object in the parameter list, or an error will occur if the
        parameter object is not in this manager."""
        return self.parameters.index(paramObj)
    
    def GetGroupedList(self):
        """returns a list of lists, where each of the lists contained in the
        main list contains parameters all of the same group."""
        list = []
        for param in self.parameters:
            for group in list:
                if group[0].group == param.group:
                    group.append(param)
                    break
            else:
                list.append([param])
        return list
    


class Fitter():
    def __init__(self, session, spinwave_domain = [], size=3, k = 100,
                 tFactor = .95):
        """The optimizer will read min_range_list, max_range_list, and change
        fit_list values to be between the minimum and maximum values at the
        same index.  Then GetResult() is called to return the list of
        eigenvalues.
        -spinwave_domain is a list of tuples (kx,ky,kz)"""
        self.spinwave_domain = spinwave_domain
        self._k = k
        self.tMin = .01
        self.tMax = 10
        self._tFactor = tFactor
        self._interactionCellDimensions = (session.MagCell.Na, session.MagCell.Nb, session.MagCell.Nc)
        #matrices is a list of 2D numpy arrays of JParam objects
        #allAtoms is a list of simpleAtoms as defined in the method Export_Aux
        self._matrices, allAtoms = session.Export_Aux(size)
        #convert to the form used by the monteCarlo simulation
        self._simpleAtoms = []
        for atom in allAtoms:
            new_atom = simpleAtom(atom.pos, atom.anisotropy, atom.spinMag)
            new_atom.interactions = atom.interactions
            new_atom.interactions.extend(atom.interCellInteractions)
            self._simpleAtoms.append(new_atom)
        
        groupList = session.parameter_manager.GetGroupedList()
        self._paramGroupList = []
        self.fit_list = []
        self.min_range_list = []
        self.max_range_list = []
        print "\ngrouplist: ", groupList
        for group in groupList:
            param = group[0]
            print "\nfit = ", param.fit
            if param.fit:
                self._paramGroupList.append(group)
                self.fit_list.append(param.default)
                try:
                    self.min_range_list.append(float(param.min))
                except ValueError:
                    self.min_range_list.append(NINF)#should be -inf 
                try:
                    self.max_range_list.append(float(param.max))
                except ValueError:
                    self.max_range_list.append(PINF)#should be +inf
        
    
    def GetResult(self):
        #Propagate any value changes through all tied parameters
        for i in range(len(self.fit_list)):
            group = self._paramGroupList[i]
            for param in group:
                #This will change the values in the matrices
                param.value = self.fit_list[i]
        monteCarloMats = []
        #Used later for spinwave calculation
        jnums = []
        jmats = []
        for i in range(len(self._matrices)):
            j11 = self._matrices[i][0][0].value
            j12 = self._matrices[i][0][1].value
            j13 = self._matrices[i][0][2].value
            j21 = self._matrices[i][1][0].value
            j22 = self._matrices[i][1][1].value
            j23 = self._matrices[i][1][2].value
            j31 = self._matrices[i][2][0].value
            j32 = self._matrices[i][2][1].value
            j33 = self._matrices[i][2][2].value
            monteCarloMats.append(numpy.array([[j11,j12,j13],
                                               [j21,j22,j23],
                                               [j31,j32,j33]]))
            #Used later for spinwave calculation
            jij=sympy.matrices.Matrix([[j11,j12,j13],
                                       [j21,j22,j23],
                                       [j31,j32,j33]])
            jnums.append(i)
            jmats.append(jij)
        
        
        
        #Run the monteCarlo simulation
        #temporary (bad) method of picking tMin/tMax
        jAvg = 0
        for mat in monteCarloMats:
            jSum = 0
            for i in range(3):
                for j in range(3):
                    val = abs(mat[i][j])
                    jSum += val
            jAvg += jSum/9
        
        self._tMax = jAvg*12
        self._tMin = jAvg/1000000000
        self._tMax = 10
        self._tMin = 1E-12
        
        
        spins = Sim_Aux(self._k, self._tMax, self._tMin, self._tFactor,
                        self._simpleAtoms, monteCarloMats)
        
        
        
        def inUnitCell(atom):
            """Returns true of the atom is in the unit cell at Na, Nb, Nc, where
            those are the dimensions of the interaction cell.  That cell is being
            used to ensure that it is completely surrounded and no interactions
            are left out.  Handles simpleAtom type"""
            if atom.pos[0] >= self._interactionCellDimensions[0] and atom.pos[0] < (self._interactionCellDimensions[0] + 1):
                if atom.pos[1] >= self._interactionCellDimensions[1] and atom.pos[1] < (self._interactionCellDimensions[1] + 1):
                    if atom.pos[2] >= self._interactionCellDimensions[2] and atom.pos[2] < (self._interactionCellDimensions[2] + 1):
                        return True
            return False
            
        
        def inInteractionCell(atoms, atom):
            """Handles simpleAtom type."""
            #First check if the atom is in the first crystallographic cell
#            if atom.pos[0] < 1.0 and atom.pos[1] < 1.0 and atom.pos[2] < 1.0:
#                return True
#Now the crystallographic unit cell being used is at Na, Nb, Nc, not 0,0,0
            if inUnitCell(atom):
                return True
                        
            #If not, check if it bonds to an atom that is
            for i in range(len(atom.interactions)):
                if inUnitCell(atoms[atom.interactions[i][0]]):
                    return True
            return False
        
        #Do the spinwave calculation
        spinwave_atoms = []
        first_cell_atoms = 0
        for i in range(len(self._simpleAtoms)):
            atom = self._simpleAtoms[i]
            if inInteractionCell(self._simpleAtoms, atom):
                Dx = atom.anisotropy[0]
                Dy = atom.anisotropy[1]
                Dz = atom.anisotropy[2]
                spinwave_atom = SpinwaveAtom(pos=atom.pos, Dx=Dx,Dy=Dy,Dz=Dz,
                                             orig_Index = i)
                rmat = findmat(spins[i])#indices match up
                spinwave_atom.spinRmatrix = rmat
                
                for interaction in atom.interactions:
                    spinwave_atom.neighbors.append(interaction[0])
                    spinwave_atom.interactions.append(interaction[1])
                spinwave_atoms.append(spinwave_atom)
                #x,y,z = spinwave_atom.pos
                #if x<1 and y<1 and z<1:
                #    first_cell_atoms +=1
        
        #Find atoms in desired cell and organize list so they come first
        tmp = []
        for i in range(len(spinwave_atoms)):
            if inUnitCell(spinwave_atoms[i]):
                first_cell_atoms += 1
                tmp.append(spinwave_atoms.pop(i))
        for a in spinwave_atoms:
            tmp.append(a)
        spinwave_atoms = tmp       
        
        #change interaction indices to match indices in new list
        for a in spinwave_atoms:
            neighborList = a.neighbors
            i = 0
            while i < len(neighborList):
                neighbor_index = neighborList[i]
                for atom_index in range(len(spinwave_atoms)):
                    atom = spinwave_atoms[atom_index]
                    if atom.origIndex == neighbor_index:
                        a.neighbors[i] = atom_index
                        i +=1
                        break
                else:#This neighbor index is not in the interaction cell
                    neighborList.pop(i)
                    a.interactions.pop(i)
        
        N_atoms=len(spinwave_atoms)
        Hsave=calculate_dispersion(spinwave_atoms,first_cell_atoms,N_atoms,jmats
                                   ,showEigs=True)
        qrange = []
        wrange = []
        q = 0
        for point in self.spinwave_domain:
            #print "point:", point
            wrange.append(calc_eigs(Hsave,point[0], point[1], point[2]))
            qrange.append(q)
            q += 1
     
        wrange=numpy.real(wrange)
        wrange=numpy.array(wrange)
        wrange=numpy.real(wrange.T)
        
        #for wrange1 in wrange:
        #    pylab.plot(qrange,wrange1,'s')
        #pylab.show()
        
        return wrange[0]#for now just one set

        
            