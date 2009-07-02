from __future__ import division
import sys

import sympy as sp
import numpy as np
import sympy.matrices as spm
from sympy import I,pi,var,exp,oo,sqrt
from sympy.physics.paulialgebra import delta
from timeit import default_timer as clock

import matplotlib as mpl
mpl.use('WxAgg')
import matplotlib.pyplot as plt

from list_manipulation import *
from subin import sub_in
import spinwaves.spinwavecalc.readfiles as rf


# Computes the inner product with a metric tensor
def inner_prod(vect1,vect2,ten = spm.Matrix([[1,0,0],
                                             [0,1,0],
                                             [0,0,1]])):
    # For column vectors -  make sure vectors match eachother as well as # of rows in tensor
    if vect1.shape == vect2.shape == (3,1) == (ten.lines,1): return (vect1.T * ten * vect2)[0]
    # For row vectors -  make sure vectors match eachother as well as # of cols in tensor
    elif vect1.shape == vect2.shape == (1,3) == (1,ten.cols): return (vect1 * ten * vect2.T)[0]
    # Everything else
    else: return None

#------------ CROSS SECTION CALC METHODS ---------------------------------------

# Generates the atom lists
def generate_atoms(N):
    """ Generates atoms """
    real_list = []; recip_list = []
    for i in range(N):
        qx = sp.Symbol('qx%i'%(i,),commutative = False)
        qy = sp.Symbol('qy%i'%(i,),commutative = False)
        qz = sp.Symbol('qz%i'%(i,),commutative = False)
        qvect = spm.Matrix([qx,qy,qz])

        lx = sp.Symbol('lx%i'%(i,),commutative = False)
        ly = sp.Symbol('ly%i'%(i,),commutative = False)
        lz = sp.Symbol('lz%i'%(i,),commutative = False)
        lvect = spm.Matrix([lx,ly,lz])

        real_list.append(lvect)
        recip_list.append(qvect)
    print "Atoms Generated!"
    return (real_list, recip_list)

# Lists of the b and b dagger operators
def generate_b_bd_operators(N):
    """Generates b and b dagger operators"""
    b_list = []; bd_list = []
    for i in range(N):
        b = sp.Symbol('b%i'%(i,), commutative = False)
        bd = sp.Symbol('bd%i'%(i,), commutative = False)

        b_list.append(b); bd_list.append(bd)
    print "Operators Generated: b, bd"
    return (b_list,bd_list)

# Generates the a and a dagger operators
def generate_a_ad_operators(N, real_list, recip_list, b_list, bd_list):
    """Generates a and a dagger operators"""
    a_list = []; ad_list = []
    t = sp.Symbol('t', commutative = True)
    for i in range(N):
        temp = []; tempd = []

        for j in range(N):
            wj = sp.Symbol('w%i'%(j,), commutative = True)
            qvect = real_list[i]
            lvect = recip_list[i]
            temp.append(exp(I*(inner_prod(qvect,lvect) - wj*t)) * b_list[j])
            tempd.append(exp(-I*(inner_prod(qvect,lvect) - wj*t)) * bd_list[j])

        a = sp.Pow(sp.sqrt(N),-1) * sum(temp)
        ad = sp.Pow(sp.sqrt(N),-1) * sum(tempd)
        a_list.append(a); ad_list.append(ad)

    print "Operators Generated: a, ad"
    return (a_list,ad_list)

# Generates the Sp and Sm operators
def generate_Sp_Sm_operators(N, atom_list, a_list, ad_list):
    """Generates S+ and S- operators"""
    Sp_list = []; Sm_list = []
    S = sp.Symbol('S', commutative = True)
    for i in range(N):
#        S = atom_list[i].spin
        Sp = sp.sqrt(2*S) * a_list[i]
        Sm = sp.sqrt(2*S) * ad_list[i]
        Sp_list.append(Sp); Sm_list.append(Sm)
    print "Operators Generated: Sp, Sm"
    return (Sp_list,Sm_list)

def generate_Sa_Sb_Sn_operators(N, atom_list, Sp_list, Sm_list):
    """Generates Sa, Sb, Sn operators"""
    Sa_list = []; Sb_list = []; Sn_list = []
    S = sp.Symbol('S', commutative = True)
    for i in range(N):
#        S = atom_list[i].spin
        Sa = ((1/2)*(Sp_list[i]+Sm_list[i])).expand()
        Sb = ((1/2)*(1/I)*(Sp_list[i]-Sm_list[i])).expand()
        Sn = (S - sp.Pow(2*S,-1) * Sm_list[i].expand() * Sp_list[i].expand()).expand()

        Sa_list.append(Sa); Sb_list.append(Sb); Sn_list.append(Sn)
    print "Operators Generated: Sa, Sb, Sn"
    return (Sa_list, Sb_list, Sn_list)

# Generates the Sx, Sy and Sz operators
def generate_Sx_Sy_Sz_operators(N, atom_list, Sa_list, Sb_list, Sn_list):
    """Generates Sx, Sy and Sz operators"""
    Sx_list = []; Sy_list = []; Sz_list = []
    S = sp.Symbol('S', commutative = True)
    rotmat = spm.eye(3)
    for i in range(N):
#        S = atom_list[i].spin
        rotmat = sp.Matrix(atom_list[i].spinRmatrix)
        loc_vect = spm.Matrix([Sa_list[i],Sb_list[i],Sn_list[i]])
        loc_vect = loc_vect.reshape(3,1)
        glo_vect = rotmat * loc_vect

        Sx = glo_vect[0].expand()
        Sy = glo_vect[1].expand()
        Sz = glo_vect[2].expand()

        Sx_list.append(Sx); Sy_list.append(Sy); Sz_list.append(Sz)
        
    #Unit vector markers
    kapx = sp.Symbol('kapx',commutative = False)
    kapy = sp.Symbol('kapy',commutative = False)
    kapz = sp.Symbol('kapz',commutative = False)

    Sx_list.append(kapx); Sy_list.append(kapy); Sz_list.append(kapz)
    print "Operators Generated: Sx, Sy, Sz"
    return (Sx_list,Sy_list,Sz_list)

# Generate Hamiltonian
def generate_Hamiltonian(N, b_list, bd_list):
    """Generates the Hamiltonian operator"""
    # Ham = Ham0 + sum over q of hbar*omega_q * bdq * bq
    # Ham0 = - S^2 N sum over rho of J(rho)
    # hbar * omega_q = 2 S {cJ(0)-cJ(q)}
    # sum over rho of J(rho) = Sum J(l-lp) from lp 0 to N l fixed
    # cJ(0) = sum over rho of J(rho)
    # cJ(q) = cJ(0)*exp(I*q*(l-lp))

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
def generate_possible_combinations(N, alist):
    """This method returns the possible operator combinations from a list of operators"""
    # For a combination to be returned, the product must have an equal number of b
    # and b dagger operators. If not, they are rejected.
    op_list = []
    alista = []
    t = sp.Symbol('t', commutative = True)
    alista = [[subelement.subs(t, 0) for subelement in element] for element in alist]

    for i in range(len(alist)):
        for j in range(len(alist)):
            vect1 = alist[i][-1]
            vect2 = alist[j][-1]

            allzerolist = [alista[i][0] for k in range(len(alista[i])-1)] + [delta(vect1,vect2)-vect1*vect2]
            otherlist = [alist[j][k] for k in range(len(alist[j])-1)] + [1]

            append_list = list_mult(allzerolist,otherlist)
            op_list.append(append_list)
    print "Generated: Possible Operator Combinations"
    return op_list

def reduce_options(N, arg):
    """
    Further reduces possible operator combinations by removing combinations if
    they are the negative of another combination or they are not time dependent
    (i.e. elastic scattering)
    """
    new = []
    t = sp.Symbol('t')
    for element in arg:
        if str(element[0]).find('t') > 0:
            new.append(element)

    for elementa in new:
        for elementb in new:
            if elementa[0].expand() == (-1*elementb[0]).expand():
                new.remove(elementa)
                new.remove(elementb)
                break
    print 'Applied: Possible Operator Reduction'
    return new

# Apply Commutation Relation
def apply_commutation(N, arg):
    """Applies the commutation relation of [b_i, bd_j] = kronecker delta _ ij"""
    # [bi,bdj] = delta_ij
    # Thus commutator = 0 (THEY COMMUTE) for i != j
    # Thus commutator = 1 for i == j
        # Then just put '+1' after commutation
    # NOTE: This method will take bd*b*bd*b to bd*(bd*b+1)*b so
    # I have replace bd_b called first but implement it inside this method too.
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

                        if g == j:
                            arg[k][i] = (arg[k][i].subs(bdg*bj, nj)).expand()
                            arg[k][i] = (arg[k][i].subs(bj*bdg, bdg*bj+1)).expand()
                            arg[k][i] = (arg[k][i].subs(bdg*bj, nj)).expand()

                        elif g != j:
                            arg[k][i] = (arg[k][i].subs(bg*bdj, 0)).expand()
                            arg[k][i] = (arg[k][i].subs(bj*bdg, 0)).expand()
                            arg[k][i] = (arg[k][i].subs(bdg*bj, 0)).expand()
                            arg[k][i] = (arg[k][i].subs(bdj*bg, 0)).expand()

                        arg[k][i] = (arg[k][i].subs((bdj*bdg), 0)).expand()
                        arg[k][i] = (arg[k][i].subs((bj*bg), 0)).expand()
                        arg[k][i] = (arg[k][i].subs((bdg*nj), 0)).expand()
                        arg[k][i] = (arg[k][i].subs((bg*nj), 0)).expand()
        print "Applied: Commutation"
        return arg

# Replaces expressions arranged by apply_commutation
def replace_bdb(N, arg):
    """Replaces bdqbq with nq"""
    # Replaces bdq*bq' with nq when q = q'
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
                        arg[k][i] = (arg[k][i].subs(bdg*bj, nj)).expand()

                    elif j != g:
                        arg[k][i] = (arg[k][i].subs((bdj*bg), 0)).expand()
                        arg[k][i] = (arg[k][i].subs((bdg*bj), 0)).expand()

                    arg[k][i] = (arg[k][i].subs((bdj*bdg), 0)).expand()
                    arg[k][i] = (arg[k][i].subs((bj*bg), 0)).expand()
                    arg[k][i] = (arg[k][i].subs((bdg*nj), 0)).expand()
                    arg[k][i] = (arg[k][i].subs((bg*nj), 0)).expand()
    print "Applied: bdq*bq Replacement"
    return arg

# 
def clean_up(N, arg, atom_list):
    S = sp.Symbol('S', commutative = True)
    for k in range(len(arg)):
        for i in range(N):
            S2coeff = coeff(arg[k][i], S**2)
            Scoeff = coeff(arg[k][i], S)
            if S2coeff != None and Scoeff != None:
                arg[k][i] = S2coeff*S**2 + Scoeff*S
            elif S2coeff != None and Scoeff == None:
                arg[k][i] = S2coeff*S**2
            elif S2coeff == None and Scoeff != None:
                arg[k][i] = Scoeff*S
            else: arg[k][i] == 0
    print "Applied: Final Clean up"
    return arg

# Inelastic Cross Section Equation
# (gamma r_0)^2 / 2 pi hbar * k'/k * N * {1/2 g F(kappa)}^2 sum over alpha,beta (delta_alpha,beta - kappa_alpha*kappa_beta)
#    * sum over l exp(i*kappa*l) X integral -oo to oo <exp{-i * kappa . u_0(0)} * exp{i * kappa . u_l(t)}>
#    <S^alpha_0(0) * S^beta_l(t)> * exp(-i omega t) dt
def generate_cross_section(N, arg, real_list, recip_list):
    """Generates the Cross-Section Formula for the one magnon case"""
    gam = sp.Symbol('gamma', commutative = True)
    r = sp.Symbol('r0', commutative = True)
    h = sp.Symbol('hbar', commutative = True)
    k = sp.Symbol('k', commutative = True)
    kp = sp.Symbol('kp', commutative = True)
    g = sp.Symbol('g', commutative = True)
    def FF(arg):
        F = sp.Function('F')
        if arg.shape == (3,1) or arg.shape == (1,3):
            return sp.Symbol("%r"%(F(arg.tolist()),),commutative = False)
    kap = spm.Matrix([sp.Symbol('kap1',commutative = False),sp.Symbol('kap2',commutative = False),sp.Symbol('kap3',commutative = False)])
    t = sp.Symbol('t', commutative = True)
    w = sp.Symbol('w', commutative = True)
    W = sp.Symbol('W', commutative = False)
    dif = sp.Symbol('diff', commutative = False)

    # Wilds for sub_in method
    A = sp.Wild('A',exclude = [0]); B = sp.Wild('B',exclude = [0]); C = sp.Wild('C',exclude = [0]); D = sp.Wild('D',exclude = [0])

    front_constant = (gam*r)**2/(2*pi*h)*(kp/k)*N
    front_func = (1./2.)*g*FF(kap)
    vanderwaals = exp(-2*W)

    temp2 = []
    temp3 = []
    temp4 = []
    
    # Grabs the unit vectors from the back of the lists. 
    unit_vect = []
    for i in range(len(arg)):
        unit_vect.append(arg[i].pop())

    # This is were the heart of the calculation comes in.
    # First the exponentials are turned into delta functions:
    #   exp(I(wq*t - w*t)) ---> delta(wq-w)
    #   exp(I(wq*t - w*t)+I*(q-qp)*l) ---> delta(wq*t-w*t+q*l-qp*l) ---> delta(wq-w)*delta(q*l-qp*l)        # NEEDS REVIEW
    for i in range(len(arg)):                                                                               # _
        for j in range(N):                                                                                  # ^
            arg[i][j] = (arg[i][j] * exp(-I*w*t)).expand()                                                  # |
            arg[i][j] = sub_in(arg[i][j],exp(A*I*t + B*I*t),sp.DiracDelta(A + B))                           # |
            arg[i][j] = sub_in(arg[i][j],exp(I*t*A + I*t*B + C),sp.DiracDelta(A*t + B*t + C))               # |
            arg[i][j] = sub_in(arg[i][j],sp.DiracDelta(A*t + B*t + C),sp.DiracDelta(A + B)*sp.DiracDelta(C))# |
            temp2.append(arg[i][j])                                                                         # |
        temp3.append(sum(temp2))                                                                            # |
    print "Applied: Delta Function Conversion"                                                              # |
    for i in range(len(temp3)):                                                                             # |
        temp4.append(unit_vect[i] * temp3[i])                                                               # V
    dif = front_func**2 * front_constant * vanderwaals * sum(temp4)#((sp.simplify(sum(temp4).expand())))#.expand()     # _

    print "Complete: Cross-section Calculation"
    return dif

def eval_cross_section(cross, qxval, qyval, qzval, etc):
    qx = sp.Symbol('qx', commutative = False)
    qy = sp.Symbol('qy', commutative = False)
    qz = sp.Symbol('qz', commutative = False)
    return cross.subs([(qx,qxval),(qy,qyval),(qy,qyval)])


def run_cross_section(interactionfile, spinfile):
    start = clock()

    # Generate Inputs
    atom_list, jnums, jmats,N_atoms_uc=rf.readFiles(interactionfile,spinfile)
    N_atoms = N_atoms_uc

    real, recip = generate_atoms(N_atoms)
    (b,bd) = generate_b_bd_operators(N_atoms)
    (a,ad) = generate_a_ad_operators(N_atoms, real, recip, b, bd)
    (Sp,Sm) = generate_Sp_Sm_operators(N_atoms, atom_list, a, ad)
    (Sa,Sb,Sn) = generate_Sa_Sb_Sn_operators(N_atoms, atom_list, Sp, Sm)
    list_print(Sa)
    list_print(Sb)
    list_print(Sn)    
    (Sx,Sy,Sz) = generate_Sx_Sy_Sz_operators(N_atoms, atom_list, Sa, Sb, Sn)
    list_print(Sx)
    list_print(Sy)
    list_print(Sz)
    print ''
    
    Ham = generate_Hamiltonian(N_atoms, b, bd)
    ops = generate_possible_combinations(N_atoms, [Sx,Sy,Sz])
#    list_print(ops)
    ops = replace_bdb(N_atoms, ops)
#    list_print(ops)
    ops = apply_commutation(N_atoms, ops)
#    list_print(ops)
    ops = clean_up(N_atoms, ops, atom_list)
#    list_print(ops)
    ops = reduce_options(N_atoms, ops)
#    list_print(ops)
    cross_sect = generate_cross_section(N_atoms, ops, real, recip)
#    list_print(ops)
    print '\n', cross_sect
    cross_sect = eval_cross_section(cross_sect, 0, 0, 0, 0)


    end = clock()
    print "\nFinished %i atoms in %.2f seconds" %(N_atoms,end-start)

##    fig = plt.figure(figsize = (10,7),facecolor = 'w')
##    str = repr(cross_sect)
##    str = str.replace('**','^').replace('DiracDelta','\delta').replace('I','\imath').replace('*',' ')#.split('+')
##    print str
##    #plt.title(str)
###    for s in str: s+'+\n'
###    fig.suptitle('$'+sum(str)+'$')
##    fig.suptitle('$'+str+'$')
##    #fig.savefig('test.pdf')
#    
#
#    end = clock()
#    print "\nFinished %i atoms in %.2f seconds" %(N_atoms,end-start)
#
##    plt.show()


#---------------- MAIN --------------------------------------------------------- 

# Will get rid of contents of main after integration is complete and I can run run_cross_section method
if __name__=='__main__':
    
    interfile = 'c:\montecarlo.txt'
    spinfile = 'c:\Spins.txt'
    
    run_cross_section(interfile,spinfile)


    ### THINGS LEFT TO DO
    # - optimize for N_atoms > 2
    # - Fix definition of F function in cross_section