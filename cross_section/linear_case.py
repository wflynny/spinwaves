from __future__ import division
import sympy as sp
from sympy import I,pi,var,exp,oo
from sympy.physics.paulialgebra import delta
from sympy.matrices import Matrix as spMatrix
#from sympy.matrices import matrix_multiply as spmat_mult
from sub_in import sub_in

#-------------------------------------------------------------------------------

# Method to grab coefficients (Ondrej Certik)
def coeff(expr, term):
    expr = sp.collect(expr, term)
    symbols = list(term.atoms(sp.Symbol))
    w = sp.Wild("coeff", exclude=symbols)
    m = expr.match(w*term+sp.Wild("rest"))
    if m:
        return m[w] 

#--------------------- LIST ALGEBRA --------------------------------------------


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


#------------ CROSS SECTION CALC METHODS ---------------------------------------


# Generates the atom array (1-D)
def generate_atoms(N):
    """ Generates atoms """
    atom_list=[]
    for i in range(N):
        Si = sp.Symbol('S%i'%(i,), commutative = False)
        atom_list.append(Si)
    print "Atoms Generated!"
    return atom_list

# Lists of the b and b dagger operators
def generate_b_bd_operators(N):
    """Generates b and b dagger operators"""
    b_list = []; bd_list = []
    for i in range(N):
        bi = sp.Symbol('b%i'%(i,), commutative = False)
        bdi = sp.Symbol('bd%i'%(i,), commutative = False)
        b_list.append(bi); bd_list.append(bdi)
    print "b,bd operators generated!"
    return (b_list,bd_list)

# Generates the a and a dagger operators
def generate_a_ad_operators(N,atom_list,b_list,bd_list):
    """Generates a and a dagger operators"""
    a_list = []; ad_list = []
    t = sp.Symbol('t', commutative = True)
    for i in range(N):
        temp = []
        tempd = []
        for j in range(N):
            wj = sp.Symbol('w%i'%(j,),commutative = True)
            temp.append(sp.exp(I*(atom_list[i]*atom_list[j] - wj*t)) * b_list[j])
            tempd.append(sp.exp(-I*(atom_list[i]*atom_list[j] - wj*t)) * bd_list[j])
        ai = sp.Pow(sp.sqrt(N),-1) * sum(temp)
        adi = sp.Pow(sp.sqrt(N),-1) * sum(tempd)
        a_list.append(ai); ad_list.append(adi)
    print "a,ad operators generated!"
    return (a_list,ad_list)

# Generates the Sp and Sm operators
def generate_Sp_Sm_operators(N,a_list,ad_list):
    """Generates S+ and S- operators"""
    Sp_list=[]; Sm_list=[]
    S = sp.Symbol('S', commutative = True)
    for i in range(N):
        Spi = sp.sqrt(2*S)*a_list[i]
        Smi = sp.sqrt(2*S)*ad_list[i]
        Sp_list.append(Spi); Sm_list.append(Smi)
    print "Sp and Sm operators generated!"
    return (Sp_list,Sm_list)

# Generates the Sx, Sy and Sz operators
def generate_Sx_Sy_Sz_operators(N,Sp_list,Sm_list):
    """Generates Sx, Sy and Sz operators"""
    Sx_list = []; Sy_list = []; Sz_list = []
    S = sp.Symbol('S', commutative = True)
    for i in range(N):
        Sxi = (1/2)*(Sp_list[i]+Sm_list[i]).expand()
        Syi = (1/2)*(1/I)*(Sp_list[i]-Sm_list[i]).expand()
        Szi = (S - sp.Pow(2*S,-1)*Sm_list[i].expand()*Sp_list[i].expand()).expand()
        Sx_list.append(Sxi); Sy_list.append(Syi); Sz_list.append(Szi)
    print "Sx, Sy, and Sz operators generated!"
    return (Sx_list,Sy_list,Sz_list)

# Generate Hamiltonian
def generate_Hamiltonian(N,b_list,bd_list):
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

# Define a method that returns the # of b/bd operators in an expression
def b_scanner(expr):
    """Finds the number of b and b dagger operators in an expression"""
    b = 0
    bd = 0
    s = str(expr)
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

# Define a method that generates the possible combinations of operators
def generate_possible_combinations(N,alist):
    """This method returns the possible operator combinations from a list of operators"""
    # For a combination to be returned, the product must have an equal number of b
    # and b dagger operators. If not, they are rejected.
    op_list = []
    alista = []
    t = sp.Symbol('t', commutative = True)
    alista = [[subelement.subs(t,0) for subelement in element] for element in alist]
    for i in range(len(alist)):
        for j in range(len(alist)):
            list1 = (alista[i][0].expand()*alist[j][0].expand()).expand()
            (b1,b2) = b_scanner(list1)
            if b1 == b2:
                app_list = list_mult([alista[i][0] for k in range(len(alista[i]))],alist[j])
                op_list.append(app_list)
    print "Possible Operator Combinations Generated!"
    return op_list

# Define method to further reduce possible combinations
# This is called after all others
def reduce_options(arg,N):
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
    print 'Possible Operator Combinations Reduced'
    return new

# Apply Commutation Relation
def apply_commutation(arg,title,N):
    """Applies the commutation relation of [b_i, bd_j] = kronecker delta _ ij"""
    # [bi,bdj] = delta_ij
    # Thus commutator = 0 (THEY COMMUTE) for i != j
    # Thus commutator = 1 for i == j
        # Then just put '+1' after commutation
    # NOTE: This method will take bd*b*bd*b and take it to bd*(bd*d+1)*d so 
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
                            arg[k][i] = (arg[k][i].subs(bdg*bj,nj)).expand()
                            arg[k][i] = (arg[k][i].subs(bj*bdg,bdg*bj+1)).expand()
                            arg[k][i] = (arg[k][i].subs(bdg*bj,nj)).expand()     
                        elif g != j:
                            arg[k][i] = (arg[k][i].subs(bg*bdj,0)).expand()
                            arg[k][i] = (arg[k][i].subs(bj*bdg,0)).expand()
                            arg[k][i] = (arg[k][i].subs(bdg*bj,0)).expand()
                            arg[k][i] = (arg[k][i].subs(bdj*bg,0)).expand()
                        arg[k][i] = (arg[k][i].subs((bdj*bdg), 0)).expand()
                        arg[k][i] = (arg[k][i].subs((bj*bg), 0)).expand()
                        arg[k][i] = (arg[k][i].subs((bdg*nj),0)).expand()
                        arg[k][i] = (arg[k][i].subs((bg*nj),0)).expand()      
        print "Commutation applied on %r!"%title
        return arg

# Replaces expressions arranged by apply_commutation
def replace_bdb(arg,N):
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
                        arg[k][i] = (arg[k][i].subs(bdg*bj,nj)).expand()
                    elif j != g:
                        arg[k][i] = (arg[k][i].subs((bdj*bg), 0)).expand()
                        arg[k][i] = (arg[k][i].subs((bdg*bj), 0)).expand()
                    arg[k][i] = (arg[k][i].subs((bdj*bdg), 0)).expand()
                    arg[k][i] = (arg[k][i].subs((bj*bg), 0)).expand()
                    arg[k][i] = (arg[k][i].subs((bdg*nj),0)).expand()
                    arg[k][i] = (arg[k][i].subs((bg*nj),0)).expand()   
    print "bdq*bq replacement applied!"
    return arg


# Cross Section Equation
# (gamma r_0)^2 / 2 pi hbar * k'/k * N * {1/2 g F(kappa)}^2 sum over alpha,beta (delta_alpha,beta - kappa_alpha*kappa_beta)
#    * sum over l exp(i*kappa*l) X integral -oo to oo <exp{-i * kappa . u_0(0)} * exp{i * kappa . u_l(t)}>
#    <S^alpha_0(0) * S^beta_l(t)> * exp(-i omega t) dt
# We want:
# (gamma r_0)^2 / 2 pi hbar * k'/k * N * {1/2 g F(kappa)}^2
#
#
def generate_cross_section(N,arg,atom_list):
    """Generates the Cross-Section Formula for the one magnon case"""
    S = sp.Symbol('S', commutative = True)
    gam = sp.Symbol('gamma', commutative = True)
    r = sp.Symbol('r0', commutative = True)
    h = sp.Symbol('hbar', commutative = True)
    k = sp.Symbol('k', commutative = True)
    kp = sp.Symbol('kp', commutative = True)
    g = sp.Symbol('g', commutative = True)
    F = sp.Function('F')
    kap = sp.Symbol('kappa', commutative = True)
    kapx = sp.Symbol('kappax', commutative = True)
    kapy = sp.Symbol('kappay', commutative = True)
    w = sp.Symbol('w', commutative = True)
    W = sp.Symbol('W', commutative = False)
    t = sp.Symbol('t', commutative = True)
    dif = sp.Symbol('diff', commutative = False)
    
    A = sp.Wild('A',exclude=[0]); B = sp.Wild('B',exclude=[0]); C = sp.Wild('C',exclude=[0]); D = sp.Wild('D',exclude=[0])
    
    front_constant = (gam*r)**2/(2*pi*h)*(kp/k)*N
    front_func = (1./2.)*g*F(kap)*exp(-2*W)

    temp2=[]
    temp3=[]
    temp4=[]

    # This is were the heart of the calculation comes in.
    # First the exponentials are turned into delta functions:
    #   exp(I(wq*t - w*t)) ---> delta(wq-w)
    #   exp(I(wq*t - w*t)+I*(q-qp)*l) ---> delta(wq*t-w*t+q*l-qp*l) ---> delta(wq-w)*delta(q*l-qp*l)        # NEEDS REVIEW
    for i in range(len(arg)):                                                                               # _
        for j in range(N):                                                                                  # ^
            arg[i][j] = (arg[i][j] * exp(-I*w*t)).expand()                                                  # |
            arg[i][j] = sub_in(arg[i][j],exp(A*I*t+B*I*t),sp.DiracDelta(A+B))                               # |
            arg[i][j] = sub_in(arg[i][j],exp(I*t*A+I*t*B+I*C+I*D),sp.DiracDelta(A*t+B*t+C+D))               # |
            arg[i][j] = sub_in(arg[i][j],sp.DiracDelta(A*t+B*t+C+D),sp.DiracDelta(A+B)*sp.DiracDelta(C+D))  # |
            temp2.append(exp(I*kap*atom_list[j]) * arg[i][j])                                               # |
        temp3.append(sum(temp2))                                                                            # |
    print "Converted to Delta Functions!"                                                                   # |
    for i in range(len(temp3)):                                                                             # |
        temp4.append((1-kapx**2)*temp3[i])                                                                  # V
    dif = front_constant * front_func**2*(sum(temp4))                                                       # _
    
    print "Cross-section calculated!"
    return dif

#---------------- MAIN --------------------------------------------------------- 

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
    ops = replace_bdb(ops,N_atoms)
    ops = apply_commutation(ops,'ops',N_atoms)
    ops = reduce_options(ops,N_atoms)
    cross_sect = generate_cross_section(N_atoms,ops,atom)   
    print ''
    
    list_print(ops)
    
    print "Cross-section =",cross_sect
    

    ### THINGS LEFT TO DO
    # - optimize for N_atoms > 2