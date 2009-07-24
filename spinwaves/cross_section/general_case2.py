from __future__ import division
import sys

import sympy as sp
import numpy as np
import sympy.matrices as spm
from sympy import I,pi,var,exp,oo,sqrt
from sympy.physics.paulialgebra import delta
from sympy.core.cache import *
from timeit import default_timer as clock
import matplotlib
matplotlib.use('WXAgg')
import pylab
from matplotlib._pylab_helpers import Gcf

from list_manipulation import *
from subin import sub_in
from printing import *
from spinwaves.spinwavecalc.readfiles import atom, readFiles
from spinwaves.spinwavecalc.spinwave_calc_file import calculate_dispersion, calc_eigs_direct
from periodictable import elements
sys.path.append('C:/tripleaxisproject-local/ tripleaxisproject/trunk/eclipse/src')
from rescalculator.lattice_calculator import Lattice, Orientation

# Computes the inner product with a metric tensor
def inner_prod(vect1,vect2,ten = spm.Matrix([[1,0,0],
                                             [0,1,0],
                                             [0,0,1]])):
    # For column vectors -  make sure vectors match eachother as well as # of rows in tensor
    if vect1.shape == vect2.shape == (3,1) == (ten.lines,1): 
        return (vect1.T * ten * vect2)[0]
    # For row vectors -  make sure vectors match eachother as well as # of cols in tensor
    elif vect1.shape == vect2.shape == (1,3) == (1,ten.cols): 
        return (vect1 * ten * vect2.T)[0]
    # Everything else
    else: 
        return None

#------------ CROSS SECTION CALC METHODS ---------------------------------------

# Lists of the b and b dagger operators
def generate_b_bd_operators(atom_list):
    """Generates b and b dagger operators"""
    b_list = []; bd_list = []
    N = len(atom_list)
    for i in range(N):
        b = sp.Symbol('b%i'%(i,), commutative = False)
        bd = sp.Symbol('bd%i'%(i,), commutative = False)

        b_list.append(b); bd_list.append(bd)
    print "Operators Generated: b, bd"
    return (b_list,bd_list)

# Generates the a and a dagger operators
#def generate_a_ad_operators(N, real_list, recip_list, b_list, bd_list):
def generate_a_ad_operators(atom_list, k, b_list, bd_list):
    """Generates a and a dagger operators"""
    a_list = []; ad_list = []
    N = len(atom_list)
    t = sp.Symbol('t', real = True)
    q = sp.Symbol('q', real = True)
    L = sp.Symbol('L', real = True)
    wq = sp.Symbol('wq', real = True)
    for i in range(N):
        temp = []; tempd = []
        
        for j in range(N):
#            q = spm.Matrix(atom_list[i].pos).T
            wj = sp.Symbol('w%i'%(j,), real = True)
            
            temp.append(exp(I*(q*L - wq*t)) * b_list[j])
            tempd.append(exp(-I*(q*L - wq*t)) * bd_list[j])
#            temp.append(exp(I*(inner_prod(q,k) - wj*t)) * b_list[j])
#            tempd.append(exp(-I*(inner_prod(q,k) - wj*t)) * bd_list[j])

        a = sp.Pow(sp.sqrt(N),-1) * sum(temp)
        ad = sp.Pow(sp.sqrt(N),-1) * sum(tempd)
        a_list.append(a); ad_list.append(ad)

    print "Operators Generated: a, ad"
    return (a_list,ad_list)

# Generates the Sp and Sm operators
def generate_Sp_Sm_operators(atom_list, a_list, ad_list):
    """Generates S+ and S- operators"""
    Sp_list = []; Sm_list = []
    N = len(atom_list)
    S = sp.Symbol('S', commutative = True)

    for i in range(N):
        Sp = sp.sqrt(2*S) * a_list[i]
        Sm = sp.sqrt(2*S) * ad_list[i]
        Sp_list.append(Sp) 
        Sm_list.append(Sm)
    print "Operators Generated: Sp, Sm"
    return (Sp_list,Sm_list)

def generate_Sa_Sb_Sn_operators(atom_list, Sp_list, Sm_list):
    """Generates Sa, Sb, Sn operators"""
    Sa_list = []; Sb_list = []; Sn_list = []
    N = len(atom_list)
    S = sp.Symbol('S', commutative = True)

    for i in range(N):
        Sa = ((1/2)*(Sp_list[i]+Sm_list[i])).expand()
        Sb = ((1/2)*(1/I)*(Sp_list[i]-Sm_list[i])).expand()
        Sn = (S - sp.Pow(2*S,-1) * Sm_list[i].expand() * Sp_list[i].expand()).expand()

        Sa_list.append(Sa)
        Sb_list.append(Sb)
        Sn_list.append(Sn)
    print "Operators Generated: Sa, Sb, Sn"
    return (Sa_list, Sb_list, Sn_list)

# Generates the Sx, Sy and Sz operators
def generate_Sx_Sy_Sz_operators(atom_list, Sa_list, Sb_list, Sn_list):
    """Generates Sx, Sy and Sz operators"""
    Sx_list = []; Sy_list = []; Sz_list = []
    N = len(atom_list)
    S = sp.Symbol('S', commutative = True)

    for i in range(N):
        rotmat = sp.Matrix(atom_list[i].spinRmatrix)
        loc_vect = spm.Matrix([Sa_list[i],Sb_list[i],Sn_list[i]])
        loc_vect = loc_vect.reshape(3,1)
        glo_vect = rotmat * loc_vect

        Sx = sp.powsimp(glo_vect[0].expand())
        Sy = sp.powsimp(glo_vect[1].expand())
        Sz = sp.powsimp(glo_vect[2].expand())

        Sx_list.append(Sx)
        Sy_list.append(Sy)
        Sz_list.append(Sz)
        
    #Unit vector markers
    kapxhat = sp.Symbol('kapxhat',commutative = False)#spm.Matrix([1,0,0])#
    kapyhat = sp.Symbol('kapyhat',commutative = False)#spm.Matrix([0,1,0])#
    kapzhat = sp.Symbol('kapzhat',commutative = False)#spm.Matrix([0,0,1])#
    
    Sx_list.append(kapxhat)
    Sy_list.append(kapyhat)
    Sz_list.append(kapzhat)
    print "Operators Generated: Sx, Sy, Sz"
    return (Sx_list,Sy_list,Sz_list)

# Generate Hamiltonian
def generate_Hamiltonian(atom_list, b_list, bd_list):
    """Generates the Hamiltonian operator"""
    # Ham = Ham0 + sum over q of hbar*omega_q * bdq * bq
    # Ham0 = - S^2 N sum over rho of J(rho)
    # hbar * omega_q = 2 S {cJ(0)-cJ(q)}
    # sum over rho of J(rho) = Sum J(l-lp) from lp 0 to N l fixed
    # cJ(0) = sum over rho of J(rho)
    # cJ(q) = cJ(0)*exp(I*q*(l-lp))
    N = len(atom_list)
    S = sp.Symbol('S', commutative = True)

    J = sp.Function('J')
    q = sp.Symbol('q', commutative = True)
    l = sp.Symbol('l', commutative = True)
    lp = sp.Symbol('lp', commutative = True)
    rho = sp.Symbol('p', commutative = True)
    rho = l - lp

    # Define Curly J function
    def cJ(N,q):
        temp = []
        for i in range(N):
            temp.append(J(0-i) * sp.exp(I * q * (0-i)))
        return sum(temp)

    # Define hbar*omega_q function
    def hbwq(N,q):
        return 2*S * (cJ(N,0) - cJ(N,q))

    Ham0 = -S**2 * N * cJ(N,0)

    # Sum over hbar*omega_q for all q
    temp2 = []
    for i in range(N):
        temp2.append(hbwq(N,i) * bd_list[i] * b_list[i])
    Ham_sum = sum(temp2)
    Ham = Ham0 + Ham_sum

    print "Generated: Hamiltonian"
    return Ham

# Define a method that generates the possible combinations of operators
def generate_possible_combinations(atom_list, alist):
    """This method returns the possible operator combinations from a list of operators"""
    # For a combination to be returned, the product must have an equal number of b
    # and b dagger operators. If not, they are rejected.
    op_list = []
    alista = []
    N = len(atom_list)
    t = sp.Symbol('t', commutative = True)
    L = sp.Symbol('L', real = True)
    q = sp.Symbol('q', real = True)
    qp = sp.Symbol('qp', real = True)
    wq = sp.Symbol('wq', real = True)
    wqp = sp.Symbol('wqp', real = True)

    alista = [[subelement.subs(t, 0) for subelement in element] for element in alist]
    for ele in alista:
        for sub in ele:
            sub = sub.subs(L,0)

    for i in range(len(alist)):
        for j in range(len(alist)):
            vect1 = alist[i][-1]
            vect2 = alist[j][-1]
            if cmp(vect1, vect2) == 0: delta = 1
            else: delta = 0

            allzerolist = [alista[i][0].subs(L,0) for k in range(len(alista[i])-1)]+[delta-vect1*vect2]
            otherlist = [alist[j][k].subs(q,qp).subs(wq,wqp) for k in range(len(alist[j])-1)]+[1]
            append_list = list_mult(allzerolist,otherlist)
            op_list.append(append_list)
    print "Generated: Possible Operator Combinations"
    return op_list

# 
def holstein(atom_list, arg):
    new = []
    N = len(atom_list)
    S = sp.Symbol('S', commutative = True)
    for k in range(len(arg)):
        temp = []
        orig = len(arg[k])
        for i in range(N):
            Snew = atom_list[i].spinMagnitude
            
            S2coeff = coeff(arg[k][i], S**2)
            Scoeff = coeff(arg[k][i], S)
            if S2coeff != None and Scoeff != None:
                temp.append((S2coeff*S**2 + Scoeff*S).subs(S,Snew))
            elif S2coeff != None and Scoeff == None:
                temp.append((S2coeff*S**2).subs(S,Snew))
            elif S2coeff == None and Scoeff != None:
                temp.append((Scoeff*S).subs(S,Snew))
        if temp != []:
            if temp[0] != 0:
                temp.append(arg[k][-1])
                new.append(temp)
    print "Applied: Holstein"
    return new

def reduce_options(atom_list, arg):
    """
    Further reduces possible operator combinations by removing combinations if
    they are the negative of another combination or they are not time dependent
    (i.e. elastic scattering)
    """
    new = []
    N = len(atom_list)
    t = sp.Symbol('t')
    for element in arg:
        if str(element[0]).find('t') > 0:
            new.append(element)

    for elementa in new:
        if elementa == 0:
            new.remove(elementa)
            break
        for elementb in new:
            if elementa[0].expand(deep = False) == (-1*elementb[0]).expand(deep = False):
                new.remove(elementa)
                new.remove(elementb)
                break
    print 'Applied: Possible Operator Reduction'
    return new

# Apply Commutation Relation
def apply_commutation(atom_list, arg):
    """Applies the commutation relation of [b_i, bd_j] = kronecker delta _ ij"""
    # [bi,bdj] = delta_ij
    # Thus commutator = 0 (THEY COMMUTE) for i != j
    # Thus commutator = 1 for i == j
        # Then just put '+1' after commutation
    # NOTE: This method will take bd*b*bd*b to bd*(bd*b+1)*b so
    # I have replace bd_b called first but implement it inside this method too.
    N = len(atom_list)
    if type(arg) == type([]):
        for k in range(len(arg)):
            for i in range(N):
                for j in range(N):
                    bj = sp.Symbol('b%i'%(j,), commutative = False)
                    bdj = sp.Symbol('bd%i'%(j,), commutative = False)
                    nj = sp.Symbol('n%i'%(j,), commutative = False)

                    for g in range(N):
                        bg = sp.Symbol('b%i'%(g,), commutative = False)
                        bdg = sp.Symbol('bd%i'%(g,), commutative = False)
                        
                        arg[k][i] = arg[k][i].subs(bg*bj,0)
                        arg[k][i] = arg[k][i].subs(bdg*bdj,0)
                        
                        if j == g:
                            arg[k][i] = arg[k][i].subs(bj*bdg, bdg*bj+1)
                        else:
                            arg[k][i] = arg[k][i].subs(bj*bdg, bdg*bj)
#
#                        if g == j:
##                            arg[k][i] = (arg[k][i].subs(bdg*bj, nj))
##                            arg[k][i] = (arg[k][i].subs(bj*bdg, (bdg*bj+1)))
#                            arg[k][i] = (arg[k][i].subs(bdg*bj, nj))
#                        elif g != j:
#                            arg[k][i] = (arg[k][i].subs(bg*bdj, 0))
#                            arg[k][i] = (arg[k][i].subs(bj*bdg, 0))
#                            arg[k][i] = (arg[k][i].subs(bdg*bj, 0))
#                            arg[k][i] = (arg[k][i].subs(bdj*bg, 0))

        print "Applied: Commutation"
        return arg

# Replaces expressions arranged by apply_commutation
def replace_bdb(atom_list, arg):
    """Replaces bdqbq with nq"""
    # Replaces bdq*bq' with nq when q = q'
    N = len(atom_list)
    for k in range(len(arg)):
        for i in range(N):
            for j in range(N):
                bj = sp.Symbol('b%i'%(j,), commutative = False)
                bdj = sp.Symbol('bd%i'%(j,), commutative = False)
                nj = sp.Symbol('n%i'%(j,), commutative = False)

                for g in range(N):
                    bg = sp.Symbol('b%i'%(g,), commutative = False)
                    bdg = sp.Symbol('bd%i'%(g,), commutative = False)

                    if j == g:
                        arg[k][i] = (arg[k][i].subs(bdg*bj, nj))

                    elif j != g:
                        arg[k][i] = (arg[k][i].subs((bdj*bg), 0))
                        arg[k][i] = (arg[k][i].subs((bdg*bj), 0))
                        

                    arg[k][i] = (arg[k][i].subs((bdj*bdg), 0))
                    arg[k][i] = (arg[k][i].subs((bj*bg), 0))
                    arg[k][i] = (arg[k][i].subs((bdg*nj), 0))
                    arg[k][i] = (arg[k][i].subs((bg*nj), 0))
            print '1', arg[k][i]
            q = sp.Symbol('q', real = True)
            qp = sp.Symbol('qp', real = True)
            wq = sp.Symbol('wq', real = True)
            wqp = sp.Symbol('wqp', real = True)
            arg[k][i] = arg[k][i].subs(qp,q).subs(wqp,wq)
            print '2', arg[k][i]
    print "Applied: bdq*bq Replacement"
    return arg

#   _        _
#  / \  |   |  \
# |   | |   |   |
#  \_/  |__ |_ /
# Inelastic Cross Section Equation
# (gamma r_0)^2 / 2 pi hbar * k'/k * N * {1/2 g F(kappa)}^2 sum over alpha,beta (delta_alpha,beta - kappa_alpha*kappa_beta)
#    * sum over l exp(i*kappa*l) X integral -oo to oo <exp{-i * kappa . u_0(0)} * exp{i * kappa . u_l(t)}>
#    <S^alpha_0(0) * S^beta_l(t)> * exp(-i omega t) dt
def generate_cross_section(atom_list, arg, q, real_list, recip_list):
    """Generates the Cross-Section Formula for the one magnon case"""
    N = len(atom_list)
    gam = 1.913 #sp.Symbol('gamma', commutative = True)
    r = sp.Symbol('r0', commutative = True)
    h = 1. # 1.05457148*10**(-34) #sp.Symbol('hbar', commutative = True)
    k = sp.Symbol('k', commutative = True)
    kp = sp.Symbol('kp', commutative = True)
    g = sp.Symbol('g', commutative = True)
    F = sp.Function('F')
    def FF(arg):
        F = sp.Function('F')
        if arg.shape == (3,1) or arg.shape == (1,3):
            return sp.Symbol("%r"%(F(arg.tolist()),),commutative = False)
    kap = spm.Matrix([sp.Symbol('kapx',commutative = False),sp.Symbol('kapy',commutative = False),sp.Symbol('kapz',commutative = False)])
    t = sp.Symbol('t', commutative = True)
    w = sp.Symbol('w', commutative = True)
    W = sp.Symbol('W', commutative = False)
    kappa = sp.Symbol('kappa', commutative = False)
    tau = sp.Symbol('tau', commutative = False)
    
    # Wilds for sub_in method
    A = sp.Wild('A',exclude = [0]); B = sp.Wild('B',exclude = [0]); C = sp.Wild('C'); D = sp.Wild('D'); 
    
    front_constant = (gam*r)**2/(2*pi*h)*(kp/k)*N
    front_func = (1./2.)*g#*F(k)
    vanderwaals = exp(-2*W)

    temp2 = []
    temp3 = []
    temp4 = []
    
    # Grabs the unit vectors from the back of the lists. 
    unit_vect = []
    kapx = sp.Symbol('kapxhat',)
    kapy = sp.Symbol('kapyhat',commutative = False)
    kapz = sp.Symbol('kapzhat',commutative = False)
    for i in range(len(arg)):
        unit_vect.append(arg[i].pop())
#    for ele in unit_vect:
#        ele = ele.subs(kapx,spm.Matrix([1,0,0]))
#        ele = ele.subs(kapy,spm.Matrix([0,1,0]))
#        ele = ele.subs(kapz,spm.Matrix([0,0,1]))
    # This is were the heart of the calculation comes in.
    # First the exponentials are turned into delta functions:
    for i in range(len(arg)):
        for j in range(N):
            arg[i][j] = sp.powsimp(arg[i][j], deep = True, combine = 'all')
            arg[i][j] = arg[i][j] * exp(-I*w*t) * exp(I*inner_prod(spm.Matrix(atom_list[j].pos).T,kap))
            arg[i][j] = sp.powsimp(arg[i][j], deep = True, combine = 'all')
            arg[i][j] = sub_in(arg[i][j],exp(I*t*A + I*t*B + C),sp.DiracDelta(A*t + B*t + C/I))#*sp.DiracDelta(C))
            arg[i][j] = sub_in(arg[i][j],sp.DiracDelta(A*t + B*t + C),sp.DiracDelta(A*h + B*h)*sp.DiracDelta(C + tau))
            arg[i][j] = sub_in(arg[i][j],sp.DiracDelta(-A - B),sp.DiracDelta(A + B))
            print arg[i][j]
    print "Applied: Delta Function Conversion"
    
#    for ele in arg:
#        for subele in ele:
#            temp2.append(subele)
#        temp3.append(sum(temp2))
#    
#    for i in range(len(temp3)):
#        temp4.append(unit_vect[i] * temp3[i])

    for k in range(len(arg)):
        temp4.append(arg[k][q])
    dif = (front_func**2 * front_constant * vanderwaals * sum(temp4))#.expand()#sp.simplify(sum(temp4))).expand()

    print "Complete: Cross-section Calculation"
    return dif

#def eval_cross_section(N, N_uc, atom_list, jmats, cross, qvals, temp, direction, lmin, lmax):
def eval_cross_section(interactionfile, spinfile, lattice, arg, 
                       tau_list, h_list, k_list, l_list, w_list, 
                       temperature, steps, eief, efixed = 14.7):
    """
    Calculates the cross_section given the following parameters:
    interactionfile, spinfile - files to get atom data
    lattice     - Lattice object from tripleaxisproject
    arg         - reduced list of operator combinations
    tau_list    - list of tau position
    w_list      - list of w's probed
    temp        - temperature
    kmin        - minimum value of k to scan
    kmax        - maximum value of k to scan
    steps       - number of steps between kmin, kmax
    eief        - True if fixed Ef
    efixed      - value of the fixed energy, either Ei or Ef
    """

    # Read files, get atom_list and such
    atom_list, jnums, jmats,N_atoms_uc=readFiles(interactionfile,spinfile)
    N_atoms = len(atom_list)
    N = N_atoms

    # TEMPORARY TO OVERRIDE ATOM_LIST ABOVE
#    atom1 = atom(pos = [0.00,0,0], neighbors = [1],   interactions = [0], int_cell = [0], 
#                 atomicNum = 26, valence = 3, spinRmatrix=spm.Matrix([[1, 0, 0],[0, 1, 0],[0, 0, 1]]))
#    atom2 = atom(pos = [0.25,0,0], neighbors = [0,2], interactions = [0], int_cell = [0], 
#                 atomicNum = 26, valence = 3, spinRmatrix=spm.Matrix([[1, 0, 0],[0, 1, 0],[0, 0, 1]]))
#    atom3 = atom(pos = [0.50,0,0], neighbors = [1,3], interactions = [0], int_cell = [0], 
#                 atomicNum = 26, valence = 3, spinRmatrix=spm.Matrix([[1, 0, 0],[0, 1, 0],[0, 0, 1]]))
#    atom4 = atom(pos = [0.75,0,0], neighbors = [2],   interactions = [0], int_cell = [0], 
#                 atomicNum = 26, valence = 3, spinRmatrix=spm.Matrix([[1, 0, 0],[0, 1, 0],[0, 0, 1]]))
#    atom_list,N_atoms_uc = ([atom1, atom2, atom3, atom4],1)
    N_atoms = len(atom_list)

    # Get Hsave to calculate its eigenvalues
    Hsave = calculate_dispersion(atom_list,N_atoms_uc,N_atoms,jmats,showEigs=False)
    print "Calculated: Dispersion Relation"

    # Generate kappa's from (h,k,l)
    kaprange = []
    kapvect = []
    #if len(h_list) == len(k_list) == len(l_list):
        #for i in range(len(h_list)):
            #kappa = lattice.modvec(h_list[i],k_list[i],l_list[i], 'latticestar')
            #kaprange.append(kappa[0])
            #kapvect.append(np.array([h_list[i],k_list[i],l_list[i]]))
##            kapvect.append(np.array([h_list[i]/kappa,k_list[i]/kappa,l_list[i]/kappa]))
    #else:
        #raise Exception('h,k,l not same lengths')
    # Generate q's from kappa and tau
    kaprange=lattice.modvec(h_list,k_list,l_list, 'latticestar')
    nqpts=len(kaprange)
    kapvect=np.empty((nqpts,3),'Float64')
    kapvect[:,0]=h_list
    kapvect[:,1]=k_list
    kapvect[:,2]=l_list
    #plusq=kappa-tau
    plusq=[]
    ones_list=np.ones((1,nqpts),'Float64')
    for tau in tau_list:
        taui=np.ones((nqpts,3),'Float64')
        taui[:,0]=ones_list*tau[0]
        taui[:,1]=ones_list*tau[1]
        taui[:,2]=ones_list*tau[2]
        kappa_minus_tau=kapvect-taui        
        plusq.append(kappa_minus_tau)
    #calculate kfki
    kfki=calc_kfki(w_list,eief,efixed)
    
    # Calculate w_q's using q
    qrange = []
#    krange = []
    wrange = []
    if 1: # Change this later
        eiglist=[]
        for q in plusq:
                eigs = calc_eigs_direct(Hsave,q[:,0],q[:,1],q[:,2])
                # Take only one set of eigs. Should probably have a flag here. 
                eiglist.append(eigs[:,0])
    print "Calculated: Eigenvalues"
    
    
    

    # Grab Form Factors
    ff_list = []
    s = sp.Symbol('s')
    for i in range(N_atoms):
        el = elements[atom_list[i].atomicNum]
        val = atom_list[i].valence
        if val != None:
            Mq = el.magnetic_ff[val].M_Q(kaprange)
        else:
            Mq = el.magnetic_ff[0].M_Q(kaprange)
        ff_list.append(Mq)
    print "Calculated: Form Factors"

    # Other Constants
    gamr0 = 2*0.2695*10**(-12) #sp.Symbol('gamma', commutative = True)
    hbar = 1.0 # 1.05457148*10**(-34) #sp.Symbol('hbar', commutative = True)
    g = 2.#sp.Symbol('g', commutative = True)
    # Kappa vector
    kap = sp.Symbol('kappa', real = True)#spm.Matrix([sp.Symbol('kapx',real = True),sp.Symbol('kapy',real = True),sp.Symbol('kapz',real = True)])
    t = sp.Symbol('t', real = True)
    w = sp.Symbol('w', real = True)
    W = sp.Symbol('W', real = True)
    tau = sp.Symbol('tau', real = True)
    q = sp.Symbol('q', real = True)
    L = sp.Symbol('L', real = True)
    boltz = 1.#1.3806503*10**(-23)     
    
    # Wilds for sub_in method
    A = sp.Wild('A',exclude = [0,t])
    B = sp.Wild('B',exclude = [0,t])
    C = sp.Wild('C')
    D = sp.Wild('D')
    K = sp.Wild('K')

    # First the exponentials are turned into delta functions:
    for i in range(len(arg)):
        for j in range(N):
            print '1', arg[i][j]
            arg[i][j] = sp.powsimp(arg[i][j])#sp.powsimp(arg[i][j], deep = True, combine = 'all')
            arg[i][j] = (arg[i][j] * exp(-I*w*t) * exp(I*kap*L)).expand()# * exp(I*inner_prod(spm.Matrix(atom_list[j].pos).T,kap))
            arg[i][j] = sp.powsimp(arg[i][j])#sp.powsimp(arg[i][j], deep = True, combine = 'all')
#            print '2', arg[i][j]
            arg[i][j] = sub_in(arg[i][j],exp(I*t*A + I*t*B + I*C + I*D + I*K),sp.DiracDelta(A*t + B*t + C + D + K))#*sp.DiracDelta(C))
#            print '3', arg[i][j]
            arg[i][j] = sub_in(arg[i][j],sp.DiracDelta(A*t + B*t + C*L + D*L + K*L),sp.DiracDelta(A*hbar + B*hbar)*sp.simplify(sp.DiracDelta(C + D + K + tau)))
#            print '4', arg[i][j]
            arg[i][j] = sub_in(arg[i][j],sp.DiracDelta(-A - B)*sp.DiracDelta(C),sp.DiracDelta(A + B)*sp.DiracDelta(C))
            print '5', arg[i][j]
    print "Applied: Delta Function Conversion"

    # Grabs the unit vectors from the back of the lists. 
    unit_vect = []
    kapxhat = sp.Symbol('kapxhat',commutative = False)
    kapyhat = sp.Symbol('kapyhat',commutative = False)
    kapzhat = sp.Symbol('kapzhat',commutative = False)
    for i in range(len(arg)):
        unit_vect.append(arg[i].pop())

    # Subs the actual values for kx,ky,kz, omega_q and n_q into the operator combos
    # The result is basically a nested list:
    #
    #    csdata     = [ op combos ]
    #    op combos  = [ one combo per atom ]
    #    1 per atom = [ evaluated exp ]
    #
    
    csdata = []
    for i in range(len(arg)):
        temp1 = []
        for j in range(len(arg[i])):
            temp2 = []
            for k in range(len(tau_list)):
                temp3 = []
                print i,j,k
                for g in range(len(qrange[k])):
                    pvalue = tau_list[k] + kapvect[g] + qrange[k][g]
                    mvalue = tau_list[k] + kapvect[g] - qrange[k][g]

                    if pvalue[0] == 0 and pvalue[1] == 0 and pvalue[2] == 0:
                        arg[i][j] = arg[i][j].subs(sp.DiracDelta(kap+tau+q),sp.DiracDelta(0))
                    else: arg[i][j] = arg[i][j].subs(sp.DiracDelta(kap+tau+q),0)
                    if mvalue[0] == 0 and mvalue[1] == 0 and mvalue[2] == 0:
                        arg[i][j] = arg[i][j].subs(sp.DiracDelta(kap+tau-q),sp.DiracDelta(0))
                    else: arg[i][j] = arg[i][j].subs(sp.DiracDelta(kap+tau-q),0)
#                    arg[i][j] = arg[i][j].subs(q,qrange[k][g])
#                    arg[i][j] = arg[i][j].subs(kap,kapvect[g])
#                    arg[i][j] = arg[i][j].subs(tau,tau_list[k])                   
                    
                    wq = sp.Symbol('wq', real = True)
                    nq = sp.Symbol('n%i'%(k,), commutative = False)

                    
                    arg[i][j] = arg[i][j].subs(wq,wrange[k][g])
                    arg[i][j] = arg[i][j].subs(w,w_calc[k][g])
                    n = sp.Pow( sp.exp(hbar*wrange[k][g]/boltz*temperature) - 1 ,-1)
                    arg[i][j] = arg[i][j].subs(nq,n)
                temp3.append(arg[i][j])
            temp2.append(temp3)
        temp1.append(temp2)
    csdata.append(temp1)
##            print arg[i][j]
#            for g in range(len(kaprange)):
#                arg[i][j] = arg[i][j].subs(kap, kaprange[g])
#            for g in range(len(krange)):
##                print 'calculating'
#                temp3 = []
#                wg = sp.Symbol('w%i'%(g,), real = True)
#                ng = sp.Symbol('n%i'%(g,), commutative = False)
##                kx = sp.Symbol('kx', real = True, commutative = True)
##                ky = sp.Symbol('ky', real = True, commutative = True)
##                kz = sp.Symbol('kz', real = True, commutative = True)
##                arg[i][j] = arg[i][j].subs(kx,krange[g][0])
##                arg[i][j] = arg[i][j].subs(ky,krange[g][1])
##                arg[i][j] = arg[i][j].subs(kz,krange[g][2])
#                arg[i][j] = arg[i][j].subs(wg,wrange[g])
#                arg[i][j] = arg[i][j].subs(w,w_calc[g])
#                nq = sp.Pow( sp.exp(hbar*wrange[g]/boltz*temp) - 1 ,-1)
#                arg[i][j] = arg[i][j].subs(ng,nq)
##                arg[i][j] = arg[i][j].subs(tau,tau_list[0])
#
#                temp3.append(arg[i][j])
##                print arg[i][j]
#            temp2.append(temp3)
##            print arg[i][j]
#        csdata.append(temp2)

    print csdata


    # Front constants and stuff for the cross-section
    front_constant = 1.0 # (gamr0)**2/(2*pi*hbar)
    front_func = (1./2.)*g#*F(k)
    vanderwaals = 1. #exp(-2*W)

    csrange = []
    if 1:
        temp1 = []
        temp2 = []
        for q in range(len(qrange)):
            print 'calculating'
            for ii in range(len(csdata)):
                for jj in range(len(csdata[ii])):
                    temp1.append(csdata[ii][jj])
            # Put on gamr0, 2pi hbar, form factor, kp/k, vanderwaals first
            dif = front_func**2 * front_constant*kfki# * kpk[q][0] * vanderwaals * ff_list[0][q]
#            for vect in unit_vect:
#                vect.subs(kapxhat,kapvect[q][0][0])
#                vect.subs(kapyhat,kapvect[q][1][0])
#                vect.subs(kapzhat,kapvect[q][2][0])
            print 'diff',dif
            print 'temp1',temp1
            print sum(temp1)
            dif = (dif * sum(temp1))# * sum(unit_vect))#.expand()
            csrange.append(dif)
    print "Calculated: Cross-section"
    csrange=np.real(csrange)
    csrange=np.array(csrange)

    xi=qrange
    yi=wrange
    zi=csrange
    Z=np.zeros((xi.shape[0],yi.shape[0])) 
    Z[range(len(xi)),range(len(yi))]=zi

    pylab.contourf(xi,yi,Z)
    pylab.colorbar()
    pylab.show()

#    pylab.contourf(qrange,wrange,csrange)
#    pylab.show()

    
def calc_kfki(w,eief,efixed):
    #eief==True if fixed Ef

    
    if eief==True:
        #fixed ef
        w_f=efixed*np.ones((len(w),1),'Float64')
        w_i=w-w_f
    else:
        #fixed ei
        w_i=efixed*np.ones((len(w),1),'Float64')
        w_f=w_i-w
        
    kfki=N.sqrt(wf/wi)
    return kfki
    
def run_cross_section(interactionfile, spinfile):
    start = clock()

    # Generate Inputs
    atom_list, jnums, jmats,N_atoms_uc=readFiles(interactionfile,spinfile)
    atom1 = atom(pos = [0.00,0,0], neighbors = [1],   interactions = [0], int_cell = [0], 
                 atomicNum = 26, valence = 3)
    atom2 = atom(pos = [0.25,0,0], neighbors = [0,2], interactions = [0], int_cell = [0], 
                 atomicNum = 26, valence = 3)
    atom3 = atom(pos = [0.50,0,0], neighbors = [1,3], interactions = [0], int_cell = [0], 
                 atomicNum = 26, valence = 3)
    atom4 = atom(pos = [0.75,0,0], neighbors = [2],   interactions = [0], int_cell = [0], 
                 atomicNum = 26, valence = 3)
#    atom_list,N_atoms_uc = ([atom1, atom2, atom3, atom4],1)
    N_atoms = len(atom_list)
    
    print N_atoms
    if 1:
        kx = sp.Symbol('kx', real = True)
        ky = sp.Symbol('ky', real = True)
        kz = sp.Symbol('kz', real = True)
        k = spm.Matrix([kx,ky,kz])
    
    (b,bd) = generate_b_bd_operators(atom_list)
#    list_print(b)
    (a,ad) = generate_a_ad_operators(atom_list, k, b, bd)
#    list_print(a)
    (Sp,Sm) = generate_Sp_Sm_operators(atom_list, a, ad)
#    list_print(Sp)
    (Sa,Sb,Sn) = generate_Sa_Sb_Sn_operators(atom_list, Sp, Sm)
#    list_print(Sa)
    (Sx,Sy,Sz) = generate_Sx_Sy_Sz_operators(atom_list, Sa, Sb, Sn)
#    list_print(Sx)
    print ''
    
    #Ham = generate_Hamiltonian(N_atoms, atom_list, b, bd)
    ops = generate_possible_combinations(atom_list, [Sx,Sy,Sz])
#    list_print(ops)
    ops = holstein(atom_list, ops)
#    list_print(ops)
    ops = apply_commutation(atom_list, ops)
#    list_print(ops)
    ops = replace_bdb(atom_list, ops)
#    list_print(ops)

    ops = reduce_options(atom_list, ops)
    list_print(ops)

    # Old Method
    #cross_sect = generate_cross_section(atom_list, ops, 1, real, recip)
    #print '\n', cross_sect
    
    if 1:
        aa = bb = cc = np.array([2.0*np.pi], 'Float64')
        alpha = beta = gamma = np.array([np.pi/2.0], 'Float64')
        vect1 = np.array([[1,0,0]])
        vect2 = np.array([[0,0,1]])
        lattice = Lattice(aa, bb, cc, alpha, beta, gamma, Orientation(vect1, vect2))
        data={}
        data['kx']=1.
        data['ky']=0.
        data['kz']=0.
        direction=data

        temperature = 100.0
        kmin = 0
        kmax = 2*sp.pi
        steps = 25

        tau_list = []
        for i in range(1):
            tau_list.append(np.array([0,0,0], 'Float64'))

        h_list = np.linspace(0,2,100)
        k_list = np.zeros(h_list.shape)
        l_list = np.zeros(h_list.shape)
        w_list = np.linspace(-4,4,100)

        efixed = 14.7 #meV
        eief = True
        eval_cross_section(interactionfile, spinfile, lattice, ops, 
                           tau_list, h_list, k_list, l_list, w_list,
                           temperature, steps, eief, efixed)

    end = clock()
    print "\nFinished %i atoms in %.2f seconds" %(N_atoms,end-start)
    
    
    #generate_output(cross_sect)


#---------------- MAIN --------------------------------------------------------- 


if __name__=='__main__':
    
    interfile = 'c:/test_montecarlo.txt'
    spinfile = 'c:/test_Spins.txt'    
    
    run_cross_section(interfile,spinfile)


    ### THINGS LEFT TO DO
    # - optimize for N_atoms > 2
    # - Fix definition of F function in cross_section