import sympy
from sympy import sum, exp

#var('N_atoms N_zone S w t')

N_atoms = 7#Symbol("N_atoms", commutative=True, real=True)
N_zone = 4#Symbol("N_zone", commutative=True, real=True)
S = 1./4.#Symbol("S", commutative=True, real=True)
w = pi#Symbol("w", commutative=True, real=True)
t = 2#Symbol("t", commutative=True, real=True)


# Define Integer Property 
def integerp(x):
    return sympify(x).is_integer


# Generate Operators
b_list=[]
bd_list=[]
def generate_b_operators(N_atoms):
    for i in range(N_atoms):
        bi = Symbol('b%i'%(i,), commutative=False, real=True)
        bdi = Symbol('bd%i'%(i,), commutative=False, real=True)
        
        b_list.append(bi)
        bd_list.append(bdi)
    return b_list
    return bd_list


a_list=[]
ad_list=[]
def generate_a_operators(N_atoms):
    for i in range(N_atoms):
        ai, adi = symbols('a%i'%(i,), 'ad%i'%(i,), commutative=False, real=True)
        
        q = symbols('q', dummy=True, integer=True)
        ai = pow(float(N_atoms),-1/2) * sum( exp(q*i*1J - w*t) * b_list[q] , (q,0,N_zone) )
        adi = pow(float(N_atoms),-1/2) * sum( exp(q*i*1J - w*t) * bd_list[q] , (q,0,N_zone) )
        
        a_list.append(ai)
        ad_list.append(adi)
    return a_list
    return ad_list

SPlus_list=[]
SMinus_list=[]
Sz_list=[]
def generate_S_operators(N_atoms):
    for i in range(N_atoms):
        SPlusi = Symbol('SPlus%i'%(i,), commutative=False, real=True)
        SMinusi = Symbol('SMinus%i'%(i,), commutative=False, real=True)
        Szi = Symbol('Sz%i'%(i,), commutative=False, real=True)
        
        SPlusi = sqrt(2*S)*a_list[i]
        SMinusi = sqrt(2*S)*ad_list[i]
        Szi = S - ad_list[i]*a_list[i]
        
        SPlus_list.append(SPlusi)
        SMinus_list.append(SMinusi)
        Sz_list.append(Szi)
        
    return SPlus_list
    return SMinus_list
    return Sz_list



generate_b_operators(N_atoms)
print "b list: ", b_list, "\n", "bd list: ", bd_list
generate_a_operators(N_atoms)
print "a list: ", a_list, "\n", "ad list: ", ad_list
generate_S_operators(N_atoms)
print "SPlus list: ", SPlus_list, "\n", "SMinus list: ", SMinus_list, "Sz list: ", Sz_list



# 5 Define Hamiltonian in terms of Splus and Sminus
# 4 Define Splus and Sminus in terms of a and ad
# 3 Define a and ad in terms of b and bd 


#def generate_hamiltonian(N_atoms):
#    for l in range(N_atoms):
#        for lp in range(N_atoms):
#            H = - sympy.sum(sympy.sum(J*(l-lp)*(Splusl*Sminuslp+Szl*Szlp)))