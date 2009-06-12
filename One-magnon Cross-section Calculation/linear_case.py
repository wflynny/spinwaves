from __future__ import division
import find as fd
import sympy as sp
from sympy import I,pi,var,exp,oo
from sympy.physics.paulialgebra import delta
import numpy as np


# Method to grab coefficients (Ondrej Certik)
def coeff(expr, term):
    expr = sp.collect(expr, term)
    symbols = list(term.atoms(sp.Symbol))
    w = sp.Wild("coeff", exclude=symbols)
    m = expr.match(w*term+sp.Wild("rest"))
    if m:
        return m[w] 


#################### LIST ALGEBRA ##################################

def list_print(lista):
    print 'printing...'
    for element in lista:
        print element
    print ''

# Define a method for multiplying two lists together
def list_mult(lista,listb):
    "Defines a way to multiply two lists of the same length"
    if len(lista) != len(listb):
        print "lists not same length"
        return []
    else:
        temp = []
        for i in range(len(lista)):
            temp.append((lista[i].expand()*listb[i].expand()).expand())
        return temp

# Define a method for adding two lists together
def list_sum(lista,listb):
    "Defines a way to add two lists of the same length"
    if len(lista) != len(listb):
        print "lists not same length"
        return []
    else:
        temp = []
        for i in range(len(lista)):
            temp.append((lista[i].expand()+listb[i].expand()).expand())
        return temp

# Define a method for multiplying a list by a scalar
def scalar_mult(scalar,alist):
    "Defines a way to multiply a list by a scalar"
    temp = []
    for i in range(len(alist)):
        temp.append(scalar*alist[i])
    return temp


############# CROSS SECTION CALC METHODS #################


# Generates the atom array (1-D)  
def generate_atoms(N):
    atom_list=[]
    for i in range(N):
        Si = sp.Symbol('L%i'%(i,), commutative = True, real = True)
        atom_list.append(Si)
    print "Atoms Generated!"
    return atom_list

# Lists of the b and b dagger operators
def generate_b_bd_operators(N):
    "Generates b and b dagger operators"
    b_list = []
    bd_list = []
    for i in range(N):
        bi = sp.Symbol('b%i'%(i,), commutative = False, real = True)
        bdi = sp.Symbol('bd%i'%(i,), commutative = False, real = True)        
        b_list.append(bi)
        bd_list.append(bdi)   
    print "b,bd operators generated!"
    return (b_list,bd_list)

# Generates the a and a dagger operators
def generate_a_ad_operators(N,atom_list,b_list,bd_list):
    "Generates a and a dagger operators"  
    a_list=[]
    ad_list=[]
    t = sp.Symbol('t', commutative = True, real = True)   
    for i in range(N):
        ai = sp.Symbol("a%i"%(i,), commutative = False, real = True)
        adi = sp.Symbol('ad%i'%(i,), commutative = False, real = True)
        temp = []
        tempd = []
        for j in range(N):      
            wj = sp.Symbol('w%i'%(j,),commutative = False, real = True)
            temp.append(sp.exp(I*(atom_list[i]*atom_list[j] - wj*t)) * b_list[j])
            tempd.append(sp.exp(-I*(atom_list[i]*atom_list[j] - wj*t)) * bd_list[j])
        ai = sp.Pow(sp.sqrt(N),-1) * sum(temp)
        adi = sp.Pow(sp.sqrt(N),-1) * sum(tempd)       
        a_list.append(ai)
        ad_list.append(adi)
    print "a,ad operators generated!"                       
    return (a_list,ad_list)

# Generates the Sp and Sm operators
def generate_Sp_Sm_operators(N,a_list,ad_list):
    "Generates S+ and S- operators"
    Sp_list=[]
    Sm_list=[]
    S = sp.Symbol('S', commutative = True, real = True)
    for i in range(N):
        Spi = sp.Symbol('Sp%i'%(i,), commutative = False, real = True)
        Smi = sp.Symbol('Sm%i'%(i,), commutative = False, real = True)
        Spi = sp.sqrt(2*S)*a_list[i]    
        Smi = sp.sqrt(2*S)*ad_list[i]
        Sp_list.append(Spi)
        Sm_list.append(Smi)
    print "Sp and Sm operators generated!"
    return (Sp_list,Sm_list)

# Generates the Sx, Sy and Sz operators
def generate_Sx_Sy_Sz_operators(N,Sp_list,Sm_list):
    "Generates Sx, Sy and Sz operators"
    Sx_list=[]
    Sy_list=[]
    Sz_list=[]
    S = sp.Symbol('S', commutative = True, real = True)
    for i in range(N):
        Sxi = sp.Symbol('Sm%i'%(i,), commutative = False, real = True)
        Syi = sp.Symbol('Sm%i'%(i,), commutative = False, real = True)
        Szi = sp.Symbol('Sm%i'%(i,), commutative = False, real = True)
        Sxi = (1/2)*(Sp_list[i]+Sm_list[i]).expand()
        Syi = (1/2)*(1/I)*(Sp_list[i]-Sm_list[i]).expand()
        Szi = (S - sp.Pow(2*S,-1)*Sm_list[i].expand()*Sp_list[i].expand()).expand()
        Sx_list.append(Sxi)
        Sy_list.append(Syi)
        Sz_list.append(Szi)
    print "Sx, Sy, and Sz operators generated!"
    return (Sx_list,Sy_list,Sz_list)

# Generate Hamiltonian
def generate_Hamiltonian(N,b_list,bd_list):
    "Generates the Hamiltonian operator"
    # Ham = Ham0 + sum over q of hbar*omega_q * bdq * bq
    # Ham0 = - S^2 N sum over rho of J(rho)
    # hbar * omega_q = 2 S {cJ(0)-cJ(q)}
    # sum over rho of J(rho) = Sum J(l-lp) from lp 0 to N l fixed
    # cJ(0) = sum over rho of J(rho)
    # cJ(q) = cJ(0)*exp(I*q*(l-lp))

    S = sp.Symbol('S', commutative = True, real = True)
    J = sp.Function('J')
    q = sp.Symbol('q', commutative = True, real = True)
    l = sp.Symbol('l', commutative = True, real = True)
    lp = sp.Symbol('lp', commutative = True, real = True)
    rho = sp.Symbol('p', commutative = True, real = True)
    rho = l - lp

    # Define Curly J function
    def cJ(N,q):
        temp = []    
        for i in range(N):
            temp.append(J(0-i)*sp.exp(I*q*(0-i)))
        return sum(temp)
    
    # Define hbar*omega_q function
    def hbwq(N,q):
        return 2*S*(cJ(N,0)-cJ(N,q))

    Ham0 = -S**2 * N * cJ(N,0)

    # Sum over hbar*omega_q for all q
    temp2 = []
    for i in range(N):
        temp2.append(hbwq(N,i)*bd_list[i]*b_list[i])
    Ham_sum = sum(temp2)
    Ham = Ham0 + Ham_sum        
    
    print "Hamiltonian generated!"
    return Ham


def b_scanner(list1,i):
    b = 0
    bd = 0
    s = str(list1[i])
    indexbd = 0
    while indexbd < len(s):
        indexbd = s.find('bd',indexbd+1)
        if indexbd > 0:
            bd = bd + 1
        if indexbd < 0:
            break
    indexb = 0
    while indexb < len(s):
        indexb = s.find('b',indexb+1)
        if indexb > 0:
            b = b + 1
        if indexb < 0:
            break           
    return (b-bd, bd)


def generate_possible_combinations(N,alist):


    op_list = []
    alista = [] 
    t = sp.Symbol('t', commutative = True, real = True)
    alista = [[subelement.subs(t,0) for subelement in element] for element in alist]

    for i in range(len(alist)):                 
        for j in range(len(alist)):
            list1 = list_mult(alista[i],alist[j])
            (b1,b2) = b_scanner(list1,0)  
            if b1 == b2:
                op_list.append(list_mult(alista[i],alist[j]))    
    print "Possible Operator Combinations Generated!"
    return op_list

# Apply Commutation Relation
def apply_commutation(arg,title,N):
    "Applies the commutation relation of [b_i, bd_j] = kronecker delta _ ij"
    # [bi,bdj] = delta_ij
    # Thus commutator = 0 (THEY COMMUTE) for i != j
    # Thus commutator = 1 for i == j
        # Then just put '+1' after commutation
    if type(arg) == type([]):
        for k in range(len(arg)):
            for i in range(N):
                for j in range(N):
                    bj = sp.Symbol('b%i'%(j,), commutative = False, real = True)
                    bdj = sp.Symbol('bd%i'%(j,), commutative = False, real = True)
                    for g in range(N):
                        bg = sp.Symbol('b%i'%(g,), commutative = False, real = True)
                        bdg = sp.Symbol('bd%i'%(g,), commutative = False, real = True)
                        #if g == j:
                        #    arg[k][i] = arg[k][i].subs(bj*bdg,bdg*bj+1)
                        arg[k][i] = (arg[k][i].subs(bg*bdg,bdg*bg+1)).expand()
                        if g != j:
                            arg[k][i] = (arg[k][i].subs(bg*bdj,bdj*bg)).expand()
                            arg[k][i] = (arg[k][i].subs(bj*bdg,bdg*bj)).expand()
        print "Commutation applied on %r!"%title
        return arg

# Replaces expressions arranged by apply_commutation
def replace_bdb(arg,N):
    "Replaces bdqbq with nq"
    # Replaces bdq*bq' with nq when q = q'
    for k in range(len(arg)):
        for i in range(N):
            for j in range(N):
                bj = sp.Symbol('b%i'%(j,), commutative = False, real = True)
                bdj = sp.Symbol('bd%i'%(j,), commutative = False, real = True)
                nj = sp.Symbol('n%i'%(j,), commutative = False, real = True)
                arg[k][i] = arg[k][i].subs((bdj*bj) , nj)
                arg[k][i] = arg[k][i].subs((bj*bdj) , nj+1)
                for g in range(N):
                    bg = sp.Symbol('b%i'%(g,), commutative = False, real = True)
                    bdg = sp.Symbol('bd%i'%(g,), commutative = False, real = True)
                    if j != g:
                        arg[k][i] = (arg[k][i].subs((bdj*bg), 0)).expand()
                        arg[k][i] = (arg[k][i].subs((bdg*bj), 0)).expand()
                    arg[k][i] = (arg[k][i].subs((bdj*bdg), 0)).expand()
                    arg[k][i] = (arg[k][i].subs((bj*bg), 0)).expand()
                    arg[k][i] = (arg[k][i].subs((bdg*nj),0)).expand()
                    arg[k][i] = (arg[k][i].subs((bg*nj),0)).expand()                           

    print "bdq*bq replacement applied!"
    return arg

def reduce_options(arg,N):
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
    return new
  
def convert_to_delta(expr):
    return expr.subs(sp.exp,sp.Lambda(x,sp.DiracDelta(x/I)))        
     

# Cross Section Equation
# (gamma r_0)^2 / 2 pi hbar * k'/k * N * {1/2 g F(kappa)}^2 sum over alpha,beta (delta_alpha,beta - kappa_alpha*kappa_beta)
#    * sum over l exp(i*kappa*l) X integral -oo to oo <exp{-i * kappa . u_0(0)} * exp{i * kappa . u_l(t)}>
#    <S^alpha_0(0) * S^beta_l(t)> * exp(-i omega t) dt
def generate_cross_section(N,arg):
    "Generates the Cross-Section Formula for the one magnon case"
    S = sp.Symbol('S', commutative = True, real = True)
    gam = sp.Symbol('gamma', commutative = True, real = True)
    r = sp.Symbol('r0', commutative = True, real = True)
    h = sp.Symbol('hbar', commutative = True, real = True)
    k = sp.Symbol('k', commutative = True, real = True)
    kp = sp.Symbol('kp', commutative = True, real = True)
    g = sp.Symbol('g', commutative = True, real = True)
    F = sp.Function('F')
    kap = sp.Symbol('kappa', commutative = True, real = True)
    kapx = sp.Symbol('kappax', commutative = True, real = True)
    kapy = sp.Symbol('kappay', commutative = True, real = True)
    w = sp.Symbol('omega', commutative = True, real = True)
    t = sp.Symbol('t', commutative = True, real = True)
    dif = sp.Symbol('diff', commutative = False, real = True)
    
    front_constant = (gam*r)**2/(2*pi*h)*(kp/k)*N
    front_func = (1./2.)*g*F(kap)

    temp2=[]
    temp3=[]
    temp4=[]
    #temp5 = [sum([exp(I*kap*element.index(subelement)) * sp.Integral(subelement * exp(-I*w*t) ,(t,-oo,oo)) for subelement in element ]) for element in arg]
    
    for i in range(len(arg)):
        for j in range(N):
            arg[i][j] * exp(-I*w*t)
            temp2.append(exp(I*kap*j) * sp.Integral(arg[i][j] * exp(-I*w*t) ,(t,-oo,oo)))
        temp3.append(sum(temp2))
    for i in range(len(temp3)):
        temp4.append((1-kapx**2)*temp3[i])
    dif = front_constant * front_func**2*(sum(temp4))
    
    print "Cross-section calculated!"
    return dif


############### MAIN ###################

if __name__=='__main__':
    # Call Methods
    N_atoms = 2
    # CAUTION!! DO NOT SET N_atoms > 15 as of 6/3/09
    atom = generate_atoms(N=N_atoms)
    (b,bd) = generate_b_bd_operators(N=N_atoms)
    (a,ad) = generate_a_ad_operators(N_atoms,atom,b,bd)
    (Sp,Sm) = generate_Sp_Sm_operators(N_atoms,a,ad)
    (Sx,Sy,Sz) = generate_Sx_Sy_Sz_operators(N_atoms,Sp,Sm)
    print ''
    
    Ham = generate_Hamiltonian(N_atoms,b,bd)
    ops = generate_possible_combinations(N_atoms,[Sx,Sy,Sz]) 
    ops = apply_commutation(ops,'ops',N_atoms)
    ops = replace_bdb(ops,N_atoms)
    ops = reduce_options(ops,N_atoms)
    cross_sect = generate_cross_section(N_atoms,ops)   
    print ''
    
    list_print(ops)
    print ''
    
    #print "Cross-section =",cross_sect
    #print "Cross-section.evalf() =",cross_sect.evalf()
    
    #print ops[0][0] == ops[3][0]
    #print ops[1][0] == (-1*ops[2][0]).expand()
    
    #print ops[1][0] + ops[2][0]

    ### THINGS LEFT TO DO
    # - Tokenize the expression to substitute dirac deltas info for exps
    # - figure out a better way to replace bd*b -> n for expressions with more than 1 bd*b term
    # - optimize for N_atoms > 2