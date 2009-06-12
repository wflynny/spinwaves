import sympy
from sympy import I

# Parameters
N = sympy.Symbol("N",commutative=True,real=True)
q = sympy.Symbol("q",commutative=True,real=True)
l = sympy.Symbol("l",commutative=True,real=True)
wq = sympy.Symbol("wq",commutative=True,real=True)
t = sympy.Symbol("t",commutative=True,real=True)
S = sympy.Symbol("S",commutative=True,real=True)

# Operators
bq = sympy.Symbol("bq",commutative=False,real=True)
bdq = sympy.Symbol("bdq",commutative=False,real=True)
al = sympy.Symbol("al",commutative=False,real=True)
adl = sympy.Symbol("adl",commutative=False,real=True)

al = sympy.Pow(N,-1./2.) * sympy.exp(I*(q*l - wq*t))*bq
adl = sympy.Pow(N,-1./2.) * sympy.exp(-I*(q*l - wq*t))*bdq

Sp = sympy.Symbol("Sp",commutative=False,real=True)
Sm = sympy.Symbol("Sm",commutative=False,real=True)
Sz = sympy.Symbol("Sz",commutative=False,real=True)

Sp = sympy.sqrt(2*S)*al
Sm = sympy.sqrt(2*S)*adl
Sz = S - sympy.Pow(2*S,-1)*Sm*Sp

