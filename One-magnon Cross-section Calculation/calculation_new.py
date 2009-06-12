import sympy
import numpy as N
from sympy import I


# Define Variables
t = sympy.Symbol('t', commutative = True, real = True)


# Create Vector Class
class Vector:
    def __init__(self,listt):
         self.x_dim = listt[0]
         self.y_dim = listt[1]
         self.z_dim = listt[2]
         self.magnitude = sympy.sqrt(self.x_dim**2 + self.y_dim**2 + self.z_dim**2)
    def get_vector(self):
        return [self.x_dim,self.y_dim,self.z_dim]


# Define Dot Product Operation
def dot(vect1, vect2):
    return vect1[0]*vect2[0] + vect1[1]*vect2[1] + vect1[2]*vect2[2]

# Define Vector Subtraction
def vect_sub(vect1,vect2):
    return Vector([vect1[0]-vect2[0],vect1[1]-vect2[1],vect1[2]-vect2[2]]).get_vector()

# Define Scalar Mult
def scal_mult(scal, vect):
    return Vector(scal*vect[0],scal*vect[1],scal*vect[2]).get_vector()

# Test Vector Creation and Dot Product
"""
v = Vector(1,0,0)
w = Vector(0,4,0)
print v.magnitude
print w.magnitude
print dot(v.get_vector(),w.get_vector())
"""

# List Multiplication Method
def list_mult(list1,list2):
    if len(list1) == len(list2):
        temp = []
        for i in range(len(list1)):
            for j in range (len(list1)):
                temp.append(list1[i]*list2[j])
        return temp
    else: return list1,list2


# Create Atom Class
class Atom:
    def __init__(self,label,position,ang_freq,spin_vector=Vector([sympy.var('x_dim'),sympy.var('y_dim'),sympy.var('z_dim')]),neighbors=None,interactions=None):
        self.label = label
        self.position = position
        self.spin_vector = spin_vector
        self.ang_freq = ang_freq
        if neighbors == None:
            self.neighbors = []
        if interactions == None:
            self.interactions = []
    # Returns the atoms spin vector
    def get_spin_vector(self):
        return self.spin_vector.get_vector()


# Generates the atom array (1-D)
atom_list = []      
def generate_atoms(N_atoms):
    for i in range(N_atoms):
        atomi = Atom(i,i,sympy.Symbol("w%i"%(i,), commutative = True, real = True),Vector([i,0,0]))
        atom_list.append(atomi)
    print "Atoms Generated!"
    return atom_list


# Lists of the b operators and b dagger operators
b_list = []
bd_list = []
# Generates the b operators
def generate_b_operators(N_int):
    "Generates operators b and b dagger"
    S = sympy.Symbol('S', commutative = True, real = True)
    
    for i in range(N_int):
        bi = sympy.Symbol('b%i'%(i,), commutative = False, real = True)
        bdi = sympy.Symbol('bd%i'%(i,), commutative = False, real = True)
        
        b_list.append(bi)
        bd_list.append(bdi)   
    print "b operators generated!"
    return b_list,bd_list
   

# Lists of the a operators and a dagger operators
a_list=[]
ad_list=[]
# Generates the a operators
def generate_a_operators(N_atoms, N_int):
    "Generates operators a and a dagger"
    
    for i in range(N_atoms):
        ai = sympy.Symbol('a%i'%(i,), commutative = False, real = True)
        adi = sympy.Symbol('ad%i'%(i,), commutative = False, real = True)
        
        temp = []
        tempd = []
        for q in range(N_int):
            if q == 0:
                temp.append(sympy.exp(I*(dot(atom_list[i].get_spin_vector(),atom_list[q].get_spin_vector()) - atom_list[q].ang_freq * t))* b_list[q])
                tempd.append(sympy.exp(-I*(dot(atom_list[i].get_spin_vector(),atom_list[q].get_spin_vector()) - atom_list[q].ang_freq * t)) * bd_list[q])
            if q > 0:
                temp[0] = temp[0] + sympy.exp(I*(dot(atom_list[i].get_spin_vector(),atom_list[q].get_spin_vector()) - atom_list[q].ang_freq * t)) * b_list[q]
                tempd[0] = tempd[0] + sympy.exp(-I*(dot(atom_list[i].get_spin_vector(),atom_list[q].get_spin_vector()) - atom_list[q].ang_freq * t)) * bd_list[q]
        ai = sympy.Pow(N_atoms,-1./2.) * temp[0]
        adi = sympy.Pow(N_atoms,-1./2.) * tempd[0]
        
        a_list.append(ai)
        ad_list.append(adi)
    print "a operators generated!"                       
    return a_list, ad_list


# Lists of the S+ operators, S- and Sz operators
SPlus_list=[]
SMinus_list=[]
Sz_list=[]   
# Generates the S operators
def generate_S_operators(N_atoms):
    "Generates operators SPlus, SMinus, and Sz"
    S = sympy.Symbol('S', commutative = True, real = True)
    
    for i in range(N_atoms):
        SPlusi = sympy.Symbol('SPlus%i'%(i,), commutative=False, real=True)
        SMinusi = sympy.Symbol('SMinus%i'%(i,), commutative=False, real=True)
        Szi = sympy.Symbol('Sz%i'%(i,), commutative=False, real=True)
        
        SPlusi = sympy.sqrt(2*S)*a_list[i]
        SMinusi = sympy.sqrt(2*S)*ad_list[i]
        Szi = S - sympy.Mul(ad_list[i],a_list[i]).expand()
        
        SPlus_list.append(SPlusi)
        SMinus_list.append(SMinusi)
        Sz_list.append(Szi)
    
    print "S operators generated!"    
    return SPlus_list, SMinus_list, Sz_list


# The Energy Exchange Function
def J(stuff):
    J = sympy.Symbol('J',commutative=False)
    return J(Vector(stuff).magnitude)


# Generates the Hamiltonian 
def generate_hamiltonian(N_atoms,N_int):
    "Generates the Hamiltonian operator"
    Ham = sympy.Symbol("Ham", commutative = False, real = True)
    S = sympy.Symbol('S', commutative = True, real = True)
    
    temp = []
    for i in range(N_atoms):
        if i == 0:
            print J(vect_sub(atom_list[0].get_spin_vector(),atom_list[i].get_spin_vector()))
            temp.append(J(vect_sub(atom_list[0].get_spin_vector(),atom_list[i].get_spin_vector()))) 
        elif i > 0:
            print "I ran",i,"time(s)"
            print temp[0]

            temp[0] = temp[0] + J(vect_sub(atom_list[0].get_spin_vector(),atom_list[i].get_spin_vector()))
    Ham0_sum = temp[0]
    print Ham0_sum
    Ham0 = - S**2 * N_atoms * Ham0_sum
    
    q = sympy.Symbol('q', Integer=True)
    temp1=[]
    for i in range(N_atoms):
        if i ==0:
            temp1.append(sympy.sum(J(vect_sub(atom_list[0].get_spin_vector(),atom_list[i].get_spin_vector()))*sympy.exp(I*dot(atom_list[q].get_spin_vector(),vect_sub(atom_list[0].get_spin_vector(),atom_list[i].get_spin_vector())))))
        elif i > 0:
            temp1[0] = temp1[0] + sympy.sum(J(vect_sub(atom_list[0].get_spin_vector(),atom_list[i].get_spin_vector()))*sympy.exp(I*dot(atom_list[q].get_spin_vector(),vect_sub(atom_list[0].get_spin_vector(),atom_list[i].get_spin_vector()))))
    homega_intsum = temp1[0]
    homega = 2*S*(Ham0_sum - homega_intsum)
    
    temp2=[]
    for i in range(N_int):
        if i == 0:
            temp2.append(  ) 
        elif i > 0:
            temp2[0] = temp2[0] + 1
    
    Ham_sum = temp2[0]
    
    Ham = Ham0 + sympy.sum(k,(i,0,N_int))

    
# Call all generation methods     
N_atoms = 10
N_int = 5
generate_atoms(N_atoms)
generate_b_operators(N_int)
generate_a_operators(N_atoms,N_int)
generate_S_operators(N_atoms)
generate_hamiltonian(N_atoms,N_int)



#for i in range(N_int):
#    print sympy.Mul(SPlus_list[0].subs(t,0),SPlus_list[i]).expand()



# Makes certain combinations of the S operators 0 and multiplies the others out
def multiply_S_operators(list1, list2):
    print "Applying commutation"
    new_list = []
    op = ""
    if cmp(list1,list2) == 0:
        if cmp(list1,Sz_list) == 0:
            op = "Sz * Sz"
            for i in range(len(list1)):
                new_list.append(sympy.Mul(list1[0].subs(t,0),list2[i]).expand())
            pf = True
        else: new_list.append(0); pf = False
    elif cmp(list1,SPlus_list) == 0:
        if cmp(list2, SMinus_list) == 0:
            op = "S+ * S-"
            for i in range(len(list1)):
                new_list.append(sympy.Mul(list1[0].subs(t,0),list2[i]).expand())
            pf = True
        else: new_list.append(0); pf = False
    elif cmp(list1,SMinus_list) == 0:
        if cmp(list2, SPlus_list) == 0:
            op = "S- * S+"
            for i in range(len(list1)):
                new_list.append(sympy.Mul(list1[0].subs(t,0),list2[i]).expand())
            pf = True
        else: new_list.append(0); pf = False
    else: new_list.append(0); pf = False
    
    print "Operators are: ",op
    if pf == True: print "Commutation applied!"
    else: print "Commutation failed :("
    
#    print "Printing..."
#    for i in range(len(new_list)):
#        print new_list[i]
#    print "Done printing"
    return new_list

# Rearranges the b and b dagger operators




