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
#        # Positions
#        #            0        1         2        3
#        pstn = [[-1,-1,-1],[0,-1,-1],[1,-1,-1],[1,-1,0],
#        #           4       5         6         7
#                [1,-1,1],[0,1,-1],[-1,-1,1],[-1,-1,0],
#        #           8        9        10       11
#                [0,-1,0],[-1,0,-1],[0,0,-1],[1,0,-1],
#        #          12      13     14        15
#                [1,0,0],[1,0,1],[0,0,1],[-1,0,1],
#        #          16      17       18       19
#                [-1,0,0],[0,0,0],[-1,1,-1],[0,1,-1],
#        #          20      21     22      23
#                [1,1,-1],[1,1,0],[1,1,1],[0,1,1],
#        #           24      25      26
#                [-1,1,1],[-1,1,0],[0,1,0]]
#        # Neighbors
#        #          0        1         2         3
#        nbrs = [[1,7,9],[0,2,8,10],[1,3,11],[2,4,8,12],
#        #          4         5         6         7
#                [3,5,13],[4,6,8,14],[5,7,15],[0,6,8,16],
#        #            8            9            10            11
#                [1,3,5,7,17],[0,10,16,18],[1,9,11,17,19],[2,10,12,20],
#        #             12            13            14           15
#                [3,11,13,17,21],[4,12,14,22],[5,13,17,23],[6,14,16,24],
#        #            16               17             18          19
#                [7,9,15,17,25],[8,10,12,14,16,26],[9,19,25],[10,18,20,26],
#        #            20         21           22          23
#                [11,19,21],[12,20,22,26],[13,21,23],[14,22,24,26],
#        #           24           25             26
#                [15,23,25],[16,18,24,26],[17,19,21,23,25]]
#        # Interactions
#        #          0         1         2         3
#        ints = [[0,7,12],[0,1,9,13],[1,2,19],[2,3,8,10],
#        #          4         5          6         7
#                [3,4,16],[4,5,11,17],[5,6,18],[6,7,8,19],
#        #            8               9              10             11
#                [8,9,10,11,20],[12,21,28,33],[13,21,22,30,34],[14,22,23,35],
#        #             12             13              14             15
#                [8,23,24,31,36],[16,24,25,37],[17,25,26,32,38],[18,26,27,39],
#        #              16                 17             18          19
#                [19,27,28,29,40],[20,29,30,31,32,41],[33,42,49],[34,42,43,51],
#        #           20           21          22          23
#                [35,43,44],[36,44,45,52],[37,45,46],[38,46,47,53],
#        #           24         25            26
#                [39,47,48],[40,48,50],[41,50,51,52,53]]

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
#        # Positions
#        #          0        1         2        3        4              
#        pstn = [[0,0,0], [1,0,0], [1,0,1], [0,0,1], [0.5,0.5,0.5],
#        #          5        6        7        8
#                [0,1,0], [1,1,0], [1,1,1], [0,1,1]]
#
#        # Neighbors
#        #           0          1          2          3              4              
#        nbrs = [[1,3,4,5], [0,2,4,6], [1,3,4,7], [0,2,4,8], [0,1,2,3,5,6,7,8],
#        #           5          6          7          8
#                [0,4,6,8], [1,4,5,7], [2,4,6,8], [3,4,5,7]]
#
#        # Interactions
#        #           0           1           2            3                    4              
#        ints = [[0,3,8,12], [0,1,9,13], [1,2,10,17], [2,3,11,16], [12,13,14,15,16,17,18,19],
#        #           5           6           7            8
#                [4,7,8,15], [4,5,9,14], [5,6,10,18], [6,7,11,19]]

#        # Positions
#        #          0        1         2        3        4              
#        pstn = [[0,0,0], [1,0,0], [-1,0,0], [0,1,0], [0,-1,0],
#        #          5        6              7               8                9
#                [0,0,1], [0,0,-1], [-0.5,-0.5,0.5], [0.5,-0.5,0.5], [-0.5,-0.5,-0.5],
#        #              10               11             12             13              14
#                [0.5,-0.5,-0.5], [-0.5,0.5,0.5], [0.5,0.5,0.5], [-0.5,0.5,-0.5], [0.5,0.5,-0.5]]
#
#        # Neighbors
#        #                 0                  1            2              3            4              
#        nbrs = [[7,8,9,10,11,12,13,14], [8,10,12,14], [7,9,11,13], [11,12,13,14], [7,8,9,10],
#        #           5             6            7          8          9
#                [7,8,11,12], [9,10,13,14], [0,2,4,5], [0,1,4,5], [0,2,4,6],
#        #           10        11          12        13         14
#                [0,1,4,6], [0,2,3,5], [0,1,3,5], [0,2,3,6], [0,1,3,6]]
#
#        # Interactions
#        #                 0                  1            2              3            4              
#        ints = [[0,1,2,3,4,5,6,7], [17,30,15,14], [22,26,27,8], [11,10,13,29], [24,19,20,23],
#        #           5             6            7          8          9
#                [9,28,11,12], [20,21,31,18], [0,25,27,28], [1,11,15,24], [22,21,26,2],
#        #           10        11          12        13         14
#                [3,19,17,31], [4,8,10,9], [5,13,12,14], [26,29,6,20], [16,7,18,30]]

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

