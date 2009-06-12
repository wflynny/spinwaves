import sympy as sp
from sympy import I,pi,var,exp,oo
from sympy.physics.paulialgebra import delta
import numpy as np

dependvar = sp.symbols('dependvar')
op_dict = dict([('op',(dependvar,[]))])

# Method to grab coefficients (Ondrej Certik)
def coeff(expr, term):
    expr = sp.collect(expr, term)
    symbols = list(term.atoms(sp.Symbol))
    w = sp.Wild("coeff", exclude=symbols)
    m = expr.match(w*term+sp.Wild("rest"))
    if m:
        return m[w] 

# Define a method for multiplying two lists together
def list_mult(lista,listb):
    "Defines a way to multiply two lists of the same length"
    if len(lista) != len(listb):
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

# Generates the atom array (1-D)  
def generate_atoms(N):
    atom_list=[]
    for i in range(N):
        Si = sp.Symbol('L%i'%(i,),commutative=True,real=True)
        atom_list.append(Si)
    print "Atoms Generated!"
    return atom_list

# Lists of the b operators and b dagger operators

# Generates the b operators
def generate_b_operators(N):
    "Generates operator b"
    b_list = []
    for i in range(N):
        bi = sp.Symbol('b%i'%(i,), commutative = False, real = True)
        b_list.append(bi)
    print "b operators generated!"
    return b_list

# Generates the b dagger operators
def generate_bd_operators(N):
    "Generates operator b dagger"
    bd_list = []
    for i in range(N):
        bdi = sp.Symbol('bd%i'%(i,), commutative = False, real = True)
        bd_list.append(bdi)   
    print "bd operators generated!"
    return bd_list


# Generates the a operators
def generate_a_operators(N,atom_list,b_list):
    "Generates operator a"  
    a_list=[]
    t = sp.Symbol('t',commutative=True,real=True)   
    for i in range(N):
        ai = sp.Symbol("a%i"%(i,),commutative=False,real=True)

        temp = []
        for j in range(N):      
            wj = sp.Symbol('w%i'%(j,),commutative=True,real=True)
            if j == 0:
                temp.append(sp.exp(I*(atom_list[i]*atom_list[j] - wj*t)) * b_list[j])
            if j > 0:
                temp[0] = temp[0] + sp.exp(I*(atom_list[i]*atom_list[j] - wj*t)) * b_list[j]
        ai = sp.Pow(N,-1./2.) * temp[0]       
        a_list.append(ai)
    print "a operators generated!"                       
    return a_list

# Generates the a operators
def generate_ad_operators(N,atom_list,bd_list):
    "Generates operator a dagger"
    ad_list=[]
    t = sp.Symbol('t',commutative=True,real=True)   
    for i in range(N):
        adi = sp.Symbol('ad%i'%(i,), commutative = False, real = True)

        tempd = []
        for j in range(N):
            wj = sp.Symbol('w%i'%(j,),commutative=True,real=True)
            if j == 0:
                tempd.append(sp.exp(-I*(atom_list[i]*atom_list[j] - wj*t)) * bd_list[j])
            if j > 0:
                tempd[0] = tempd[0] + sp.exp(-I*(atom_list[i]*atom_list[j] - wj*t)) * bd_list[j]
        adi = sp.Pow(N,-1./2.) * tempd[0]
        ad_list.append(adi)
    print "ad operators generated!"                       
    return ad_list


# Generates the Sp operators
def generate_Sp_operators(N,a_list):
    "Generates operator S+"
    Sp_list=[]
    S = sp.Symbol('S', commutative=True, real=True)
    for i in range(N):
        Spi = sp.Symbol('Sp%i'%(i,), commutative=False, real=True)
        Spi = sp.sqrt(2*S)*a_list[i]
        Sp_list.append(Spi)
    print "Sp operators generated!"
    return Sp_list

# Generates the Sm operators
def generate_Sm_operators(N,a_list):
    "Generates operator S-"
    Sm_list=[]
    S = sp.Symbol('S', commutative=True, real=True)
    for i in range(N):
        Smi = sp.Symbol('Sm%i'%(i,), commutative=False, real=True)
        Smi = sp.sqrt(2*S)*a_list[i]
        Sm_list.append(Smi)
    print "Sm operators generated!"
    return Sm_list

# Generates the Sz operators
def generate_Sz_operators(N,a_list,ad_list):
    "Generates operator Sz"
    Sz_list=[]
    S = sp.Symbol('S', commutative=True, real=True)
    for i in range(N):
        Szi = sp.Symbol('Sm%i'%(i,), commutative=False, real=True)
        Szi = (S - sp.Pow(2*S,-1)*ad_list[i].expand()*a_list[i].expand()).expand()
        #Szi = (S - sp.Pow(2*S,-1)*Sm_list[i].expand()*Sp_list[i].expand()).expand()
        Sz_list.append(Szi)
    print "Sz operators generated!"
    return Sz_list