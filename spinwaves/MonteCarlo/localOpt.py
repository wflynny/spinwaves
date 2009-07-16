'''
Created on Jul 16, 2009

@author: tsarvey
'''
import copy
from ctypes import c_float, c_int
from numpy import cos, sin, arctan, arccos, pi
from spinwaves.utilities.mpfit import mpfit
from CSim import passAtoms, passMatrices, loadLib
from simple import readFile


#Populate C atom list
def optimizeSpins(interactionFile, inFile, outFile):
    """This will use the local optimizer, mpfit, to perfect the values of ground
    state spins.  'inFile is a spin file with the same format that the 
    montecarlo simulation creates.  'outFile' is the path of the file to output,
    with locally optimized spins.  It is the same format as inFile."""
    
    atoms, jMatrices = readFile(interactionFile)
    atomListPointer, nbr_pointers, mat_pointers = passAtoms(atoms)
    matPointer = passMatrices(jMatrices)
    c_lib = loadLib()
    
    def hamiltonian(p, fjac=None, y=None, err=None):
        """p is the parameters and will be formatted as:
        [theta1, phi1, theta2, phi2, ...] where theta and phi specify the
        angular coordinates of the spin."""
        #copy new spin values to C spinwave list
        for i in range(len(atoms)):
            theta = p[2*i]
            phi = p[(2*i) + 1]
            r = atoms[i].spinMag
            
            Sx = r*cos(theta)*sin(phi)
            Sy = r*sin(theta)*sin(phi)
            Sz = r*cos(phi)
            
            c_lib.setSpin(atomListPointer, c_int(i), c_float(Sx), c_float(Sy), c_float(Sz))
        
        status = 0
        result = c_lib.totalE(atomListPointer, len(atoms), matPointer)
        print float(result)
        return [status, [result]]
    
    #populate initial p list
    p0 = []
    for atom in atoms:
        x = atom.s[0]
        y = atom.s[1]
        z = atom.s[2]
        r = atom.spinMag
        theta = arctan(y/x)#-pi/2 to pi/2
        #range will now be from -pi to pi
        if(x < 0 and y < 0):
            theta -= pi
        elif(x<0):
            theta += pi
        phi = arccos(z/r)#0 to pi
        p0.append(theta)
        p0.append(phi)
        
    parbase={'value':0., 'limited':[0,0], 'limits':[0.,0.]}
    parinfo=[]
    for i in range(len(p0)):
        parinfo.append(copy.deepcopy(parbase))
    for i in range(len(p0)):
        parinfo[i]['value']=p0[i]
        if(i/2 == 0):#even index means theta
            parinfo[i]['limits'][0]= -pi
            parinfo[i]['limits'][1] = pi
            parinfo[i]['limited'][0] = 1
            parinfo[i]['limited'][1] = 1
        else:#phi
            parinfo[i]['limits'][0]= 0
            parinfo[i]['limits'][1] = pi
            parinfo[i]['limited'][0] = 1
            parinfo[i]['limited'][1] = 1
        
    
    m = mpfit(hamiltonian, p0, parinfo = parinfo)
    return m.params

if __name__ == '__main__':
    print optimizeSpins('C:\\export.txt', 'C:\\spins.txt', 'C:\\spins2.txt')