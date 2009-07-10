import sympy as sp
import sympy.matrices as spm
import numpy as N
from spinwaves.spinwavecalc.spinwave_calc_file import calculate_dispersion, calc_eigs
import spinwaves.spinwavecalc.readfiles as rf


def gen_chain(num, anti = False):
    ats = [rf.atom(pos = [i,0,0]) for i in range(num)]
    spin = rf.findmat(N.array([0,0,1]))
    
    if num == 2: 
        ints = [[0],[0]]
        nbrs = [[1],[0]]
        for i in range(num): 
            ats[i].interactions = ints[i]
            ats[i].neighbors = nbrs[i]
            ats[i].spinRmatrix = spin
            ats[i].spinMagnitude = 1
    elif num == 3:
        ints = [[0],[0,1],[1]]
        nbrs = [[1],[0,2],[1]]
        for i in range(num): 
            ats[i].interactions = ints[i]
            ats[i].neighbors = nbrs[i]
            ats[i].spinRmatrix = spin
            ats[i].spinMagnitude = 1
    if anti == True:
        spin = rf.findmat(N.array([0,0,-1]))
        for i in range(1,num,2): 
            ats[i].spinRmatrix = spin
            ats[i].spinMagnitude = 1
    for i in range(len(ats)):
        print ats[i].spinRmatrix
    return ats

def square(anti = False):
    ats = [rf.atom() for i in range(5)]
    spin = rf.findmat(N.array([0,0,1]))
    # Positions
    #        0        1        2        3        4
    pstn = [[0,0,0],[0,1,0],[1,0,0],[0,-1,0],[-1,0,0]]
    # Neighbors
    #        0        1    2   3   4
    nbrs = [[1,2,3,4],[0],[0],[0],[0]]
    # Interactions
    #        0        1    2   3   4
    ints = [[0,1,2,3],[0],[1],[2],[3]]

    for i in range(len(ats)):
        ats[i].pos = N.array(pstn[i])
        ats[i].neighbors = nbrs[i]
        ats[i].interactions = ints[i]
        ats[i].spinRmatrix = spin
        ats[i].spinMagnitude = 1
        
    if anti == True:
        spin = rf.findmat(N.array([0,0,-1]))
        ats[0].spinMagnitude = 1
        ats[0].spinRmatrix = spin
    return ats

def gen_cube(body = False, face = False, anti = False):
    ats = []
    spin = rf.findmat(N.array([0,0,1]))
    # Cubic
    if not body and not face:

        ats = [rf.atom() for i in range(7)]
        # Positions
        #          0        1         2        3        4         5        6
        pstn = [[0,0,0], [1,0,0], [-1,0,0], [0,1,0], [0,-1,0], [0,0,1], [0,0,-1]]

        # Neighbors
        #            0          1    2    3    4    5    6
        nbrs = [[1,2,3,4,5,6], [0], [0], [0], [0], [0], [0]]

        # Interactions
        #            0          1    2    3    4    5    6
        ints = [[0,1,2,3,4,5], [0], [1], [2], [3], [4], [5]]

        # Plug in values for positions and neighbors
        for i in range(len(ats)):
            ats[i].pos = N.array(pstn[i])
            ats[i].neighbors = nbrs[i]
            ats[i].interactions = ints[i]
            ats[i].spinRmatrix = spin
            ats[i].spinMagnitude = 1

    # BCC - Body Centered Cubic
    elif body and not face:   

        # Positions
        #          0            1               2               3               4              
        pstn = [[0,0,0], [-0.5,0.5,-0.5], [0.5,0.5,-0.5], [0.5,0.5,0.5], [-0.5,0.5,0.5],
        #               5                6                7               8               
                [-0.5,-0.5,-0.5], [0.5,-0.5,-0.5], [0.5,-0.5,0.5], [-0.5,-0.5,0.5]]

        # Neighbors
        #          0               1   2   3   4              
        nbrs = [[1,2,3,4,5,6,7,8], [], [], [], [],
        #        5   6   7   8               
                [], [], [], []]
        
        # Interactions
        #          0               1   2   3   4              
        ints = [[0,1,2,3,4,5,6,7], [], [], [], [],
        #        5   6   7   8               
                [], [], [], []]

        ats = [rf.atom() for i in range(len(pstn))]
        for i in range(len(ats)):
            ats[i].pos = N.array(pstn[i])
            ats[i].neighbors = nbrs[i]
            ats[i].interactions = ints[i]
            ats[i].spinRmatrix = spin
            ats[i].spinMagnitude = 1

    # FCC - Face Centered Cubic
    elif face and not body:
        # Positions
        #          0          1             2            3              4              
        pstn = [[0,0,0], [0.5,0,0.5], [-0.5,0,0.5], [0.5,0,-0.5], [-0.5,0,-0.5],
        #           5             6             7             8               9
                [0,0.5,0.5], [0,-0.5,0.5], [0,0.5,-0.5], [0,-0.5,-0.5], [0.5,0.5,0], 
        #             10            11            12
                [-0.5,0.5,0], [0.5,-0.5,0], [-0.5,-0.5,0]]
        # Neighbors
        #          0                               1             2            3              4              
        nbrs = [[1,2,3,4,5,6,7,8,9,10,11,12], [], [], [], [],
        #           5             6             7             8               9
                [], [], [], [], [], 
        #             10            11            12
                [], [], []]

        # Interactions
        #          0                               1             2            3              4              
        ints = [[0,1,2,3,4,5,6,7,8,9,10,11], [], [], [], [],
        #           5             6             7             8               9
                [], [], [], [], [], 
        #             10            11            12
                [], [], []]
        
        ats = [rf.atom() for i in range(len(pstn))]
        for i in range(len(ats)):
            ats[i].pos = N.array(pstn[i])
            ats[i].neighbors = nbrs[i]
            ats[i].interactions = ints[i]
            ats[i].spinRmatrix = spin
            ats[i].spinMagnitude = 1

    if anti == True:
        spin = rf.findmat(N.array([0,0,-1]))
        ats[0].spinMagnitude = 1
        ats[0].spinRmatrix = spin

    return ats


#-----------------------------------------------------------------------------------------------------

#calculate_dispersion(atom_list,N_atoms_uc,N_atoms,Jij,showEigs=False):
def run_dispersion(lattice, num_uc = 1, eigs = False):
    num = len(lattice)
    print num

    inters = [item for atom in lattice for item in atom.interactions]
    max_inter = max(inters)
    J = sp.Symbol('J', real = True)
    Jij = [N.matrix([[J,0,0],[0,J,0],[0,0,J]]) for i in range(max_inter+1)]
    Hsave = calculate_dispersion(lattice, num_uc, num, Jij, showEigs = eigs)

    if eigs == True:
        temp = 1
        x = sp.Symbol('x')
        for i in range(Hsave.cols): temp *= x-Hsave[i,i]
        poly = sp.simplify(temp.expand())
        print poly
    return Hsave[0]

def run_dispersion_from_file(interactionfile,spinfile, eigs=False):
    
    atom_list, jnums, jmats,N_atoms_uc=rf.readFiles(interactionfile,spinfile)
    N_atoms=len(atom_list)
    Hsave = calculate_dispersion(atom_list,N_atoms_uc,N_atoms,jmats,showEigs=eigs)

    if eigs == True:
        temp = 1
        x = sp.Symbol('x')
        for i in range(Hsave.cols): temp *= x-Hsave[i,i]
        poly = sp.simplify(temp.expand())
        print poly
    return sp.simplify(Hsave[0])

#------------------------------------------------------------------------------------------------------
# Run Dispersion from files
if 0:   
    interfile = 'c:\montecarlo.txt'
    spinfile = 'c:\Spins.txt'    
    
    run_dispersion_from_file(interfile,spinfile)

# Ferromagnetic Cases
if 0:
    chain2 = gen_chain(2)
    chain3 = gen_chain(3)
    sq = square()
    cube = gen_cube()
    bcc = gen_cube(body = True)
    fcc = gen_cube(face = True)
    
    Hchain2 = run_dispersion(chain2, num_uc = 1)#, eigs=True)
    Hchain3 = run_dispersion(chain3, num_uc = 1)#, eigs=True)
    Hsq = run_dispersion(sq, num_uc = 1)#, eigs=True)
    Hcube = run_dispersion(cube, num_uc = 1)#, eigs=True)
    Hbcc = run_dispersion(bcc, num_uc = 1)#, eigs=True)
    Hfcc = run_dispersion(fcc, num_uc = 1)#, eigs=True)
    
    print '\nchain2 =', Hchain2
    print '\nchain3 =', Hchain3
    print '\nsq     =', Hsq
    print '\ncube   =', Hcube
    print '\nbcc    =', Hbcc
    print '\nfcc    =', Hfcc

# Antiferromagnetic Cases
if 1:
    achain2 = gen_chain(2, anti=False)
    #achain3 = gen_chain(3, anti=True)
    #asq = square(anti=True)
    #acube = gen_cube(anti=True)
    #abcc = gen_cube(body = True, anti=True)
    #afcc = gen_cube(face = True, anti=True)
    
    aHchain2 = run_dispersion(achain2, num_uc = 1)#, eigs=True)
    #aHchain3 = run_dispersion(achain3, num_uc = 1)#, eigs=True)
    #aHsq = run_dispersion(asq, num_uc = 1)#, eigs=True)
    #aHcube = run_dispersion(acube, num_uc = 1)#, eigs=True)
    #aHbcc = run_dispersion(abcc, num_uc = 1)#, eigs=True)
    #aHfcc = run_dispersion(afcc, num_uc = 1)#, eigs=True)
    
    print '\nachain2 =', aHchain2
    #print '\nachain3 =', aHchain3
    #print '\nasq     =', aHsq
    #print '\nacube   =', aHcube
    #print '\nabcc    =', aHbcc
    #print '\nafcc    =', aHfcc

