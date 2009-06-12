import sympy as sp
from sympy import oo,I,cos,sin,exp

x = sp.Symbol("X",commutative=True,real=True)
y = sp.Symbol("Y",commutative=False,real=True)

print x*y

zlist=[]

for i in range(10):
    zi = sp.Symbol("z%i"%(i,),commutative=False,real=True)
    zlist.append(zi)
    
print zlist

print x**(1/2)

print exp(x)
n=sp.Symbol("n",commutative=True,real=True)
print sp.sum(6*n**2 + 2**n, (n, 0, 6))


print sp.solve(x**2+x+1,x)

#Convolution
tau = sp.Symbol('tau')
t,T = sp.symbols('tT')
print sp.integrate(tau**2*sp.DiracDelta(tau-(t-T)),(tau,-oo,oo))
