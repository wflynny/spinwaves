# Import Math Stuff
from sympy import pi,exp,I,cos,sin
import sympy

#########BEGIN DOCUMENT########################################################

#Create the constants/variables
#    S - Spin
#    N - Number of Atoms in Crystal

S = sympy.Symbol('s', commutative=True, real=True)
N = 10


# Create the operators
#    splus    - Raising Spin Operator
#    sminus   - Lowering Spin Operator
#    bq       - Annihilation Operator for wave q
#    bdq      - Creation Operator for wave q
splus  = sympy.Symbol('splus', commuatative=False, real=True)
sminus = sympy.Symbol('sminus', commuatative=False, real=True)



#Define value for operators
for q in range(N):
    bq     = sympy.Symbol('b%q'%(q,),commutative=False, real=True)
    bdq    = sympy.Symbol('bd%q'%(q,),commutative=False, real=True)
    splus  = sympy.Pow(2*S/N, 1/2)*sympy.Sum(exp(I*(q*l - w*t))*bq,q)
    sminus = sympy.Pow(2*S/N, 1/2)*sympy.Sum(exp(-I*(q*l - w*t))*bdq,q)