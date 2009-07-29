'''
Created on Jul 16, 2009

@author: tsarvey
'''
import copy
#from ctypes import c_float, c_int
from numpy import cos, sin, arctan, arccos, pi
from spinwaves.utilities.mpfit import mpfit
#from CSim import passAtoms, passMatrices, loadLib
#from simple import readFile
import spinwaves.spinwavecalc.readfiles as rf
import numpy as np
import sympy as sp
from scipy.sparse import bsr_matrix
from scipy.optimize import fmin_l_bfgs_b


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

def gen_Ham(atom_list):
    """ Creates a scipy bsr sparse array """
    N_atoms = len(atom_list)
    jij_values = []
    jij_columns = []
    jij_rowIndex = []   
    
    # Counts total number of interactions: needed for row indexing
    num_inters = 0
    # Scan through atom_list
    for i in range(N_atoms):
        print 'atom %i'%(i,)
        ints_list = atom_list[i].interactions
        nbrs_list = atom_list[i].neighbors
        nbrs_ints = [ (nbrs_list[x],ints_list[x]) for x in range(len(nbrs_list)) ]
        nbrs_ints.sort()

        # Now we have a sorted list of (nbr,intr) tuples from lowest neighbor to highest neighbor
        # Scan through interactions
        for j in range(len(nbrs_ints)):
            nbr = nbrs_ints[j][0]
            intr = nbrs_ints[j][1]
            
            print 'atom',i, 'nbr',nbr, 'inter',intr

            #Get an interaction matrix
            curr_mat = jmats[intr].tolist()
            curr_mat = np.array(curr_mat, dtype=np.float64)

            # Values   = current matrix
            # Columns  = the current neighbor
            # RowIndex = total number of interactions 
            jij_values.append(curr_mat)
            jij_columns.append(nbr)
            if j == 0:
                jij_rowIndex.append(num_inters)
            
            # Increase number of total interactions
            num_inters = num_inters + 1
    # Add dummy index to rowIndex
    jij_rowIndex.append(len(jij_values))

    # Convert to numpy arrays
    jij_values = np.array(jij_values)
    jij_columns = np.array(jij_columns)
    jij_rowIndex = np.array(jij_rowIndex)
    
    # Create Sparse Array
    jij = bsr_matrix( (jij_values,jij_columns,jij_rowIndex), shape=(3*N_atoms,3*N_atoms) ).todense()
    sij = gen_spinVector(atom_list)
    sijT = sij.reshape((1,9))
    anis = gen_anisotropy(atom_list)
    
    res1 = sijT * jij
    Hij = np.dot(res1,sij).flat[0]

    print 'Hij', Hij
    print 'anis', anis
    
    Ham = - Hij - anis

    return Ham

def gen_spinVector(atom_list):
    N_atoms = len(atom_list)
    vect_list = []
    for i in range(N_atoms):
        spinvect = atom_list[i].spin
        vect_list.append(spinvect[0])
        vect_list.append(spinvect[1])
        vect_list.append(spinvect[2])
    return np.array(vect_list)

def gen_anisotropy(atom_list):
    N_atoms = len(atom_list)
    anis_vect = []
    for i in range(N_atoms):
        anis_vect.append(atom_list[i].Dx)
        anis_vect.append(atom_list[i].Dy)
        anis_vect.append(atom_list[i].Dz)
    anis_vect = np.array(anis_vect)
    spin_vect = gen_spinVector(atom_list)
    return np.dot(anis_vect, spin_vect**2)

def local_optimizer(interfile, spinfile, outfile):
    atom_list, jnums, jmats,N_atoms_uc=rf.readFiles(interfile,spinfile)
    N_atoms = len(atom_list)
    
    def hamiltonian(p, fjac=None, y=None, err=None):
        N_atoms = len(atom_list)
        fake_atom_list = []

        for i in range(N_atoms):
            S = atom_list[i].spinMagnitude
            theta = p[2*i]
            phi   = p[2*i + 1]
            
            Sx = S*cos(theta)*cos(phi)
            Sy = S*cos(theta)*sin(phi)
            Sz = S*sin(theta)
            
            nbrs  = atom_list[i].neighbors
            intrs = atom_list[i].interactions
            new_atom = rf.atom( spinMag = S, spin = np.array([Sx,Sy,Sz]), neighbors = nbrs, interactions = intrs )
            fake_atom_list.append( new_atom )

        status = 0

        result = gen_Ham(fake_atom_list)
        print float(result)
        return [status, [result]]        
    # Ham
    # [S cos(t)cos(p), S cos(t)sin(p), S sin(t)] Jij [S cos(t)cos(p), S cos(t)sin(p), S sin(t)] - [Dx,Dy,Dz].[S^2 cos^2(t)cos^2(p), S^2 cos^2(t)sin^2(p), S^2 sin^2(t)]
    # derivative t
    # [-S sin(t)cos(p), -S sin(t)sin(p), S cos(t)] Jij [-S sin(t)cos(p), -S sin(t)sin(p), S cos(t)] - [Dx,Dy,Dz].[-2 S^2 cos(t)sin(t)cos^2(p), -2 S^2 cos(t)sin(t)sin^2(p), 2 S^2 cos(t)sin(t)]
    # derivative p
    # [-S cos(t)sin(p), S cos(t)cos(p), 0] Jij [-S cos(t)sin(p), S cos(t)cos(p), 0] - [Dx,Dy,Dz].[-2 S^2 cos^2(t)sin(p)cos(p), 2 S^2 cos^2(t)cos()sin(p), 0]
    
    def dEdphi():
        N_atoms = len(atom_list)
        fake_atom_list = []
        
        for i in range(N_atoms):
            S = atom_list[i].spinMagnitude
            tehta = p[2*i]
            phi   = p[2*i + 1]
            
            Sx = -S*cos(theta)*sin(phi)
            Sy = S*cos(theta)*cos(phi)
            Sz = 0
            
            nbrs  = atom_list[i].neighbors
            intrs = atom_list[i].interactions
            new_atom = rf.atom( spinMag = S, spin = np.array([Sx,Sy,Sz]), neighbors = nbrs, interactions = intrs )
            fake_atom_list.append( new_atom )
        status = 0

        result = gen_Ham(fake_atom_list)
        print float(result)
        return [status, [result]] 

    def dEdtheta():
        N_atoms = len(atom_list)
        fake_atom_list = []
        
        for i in range(N_atoms):
            S = atom_list[i].spinMagnitude
            tehta = p[2*i]
            phi   = p[2*i + 1]
            
            Sx = -S*sin(theta)*cos(phi)
            Sy = -S*sin(theta)*sin(phi)
            Sz = S*cos(theta)
            
            nbrs  = atom_list[i].neighbors
            intrs = atom_list[i].interactions
            new_atom = rf.atom( spinMag = S, spin = np.array([Sx,Sy,Sz]), neighbors = nbrs, interactions = intrs )
            fake_atom_list.append( new_atom )
        status = 0

        result = gen_Ham(fake_atom_list)
        print float(result)
        return [status, [result]]      

    #populate initial p list
    p0 = []
    for i in range(N_atoms):
        sx = atom_list[i].spin[0]
        sy = atom_list[i].spin[0]
        sz = atom_list[i].spin[0]
        s  = atom_list[i].spinMagnitude
        
        theta = arccos(sz/s)
        phi   = np.arctan2(sy,sx)
        
        p0.append(theta)
        p0.append(phi)

    parbase={'value':0., 'limited':[0,0], 'limits':[0.,0.]}
    parinfo=[]
    for i in range(len(p0)):
        parinfo.append(copy.deepcopy(parbase))
    for i in range(len(p0)):
        parinfo[i]['value']=p0[i]
        if(i%2 == 0):#even index means theta
            parinfo[i]['limits'][0]= -pi
            parinfo[i]['limits'][1] = pi
            parinfo[i]['limited'][0] = 1
            parinfo[i]['limited'][1] = 1
        else:#phi
            parinfo[i]['limits'][0]= 0
            parinfo[i]['limits'][1] = pi
            parinfo[i]['limited'][0] = 1
            parinfo[i]['limited'][1] = 1
        
    m = fmin_l_bfgs_b(hamiltonian, p0, bounds = parinfo['limits'])
    return m.params

if __name__ == '__main__':
    #print optimizeSpins('C:\\export.txt', 'C:\\spins.txt', 'C:\\spins2.txt')
    interfile = 'c:/test_montecarlo.txt'
    spinfile  = 'c:/test_Spins.txt'

    atom_list, jnums, jmats,N_atoms_uc=rf.readFiles(interfile,spinfile)
    N_atoms = len(atom_list)

    Ham = gen_Ham(atom_list)
    print Ham
    
    print local_optimizer(interfile, spinfile, 'C:\\spins_new.txt')
    
    
    