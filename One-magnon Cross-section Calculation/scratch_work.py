from sympy import *
import sympy

var('x n a b')
print sum(6*n**2 + 2**n, (n, a, b))

print (1/cos(x)).series(x, 0, 6)

print diff(cos(x**2)**2 / (1+x), x)

print exp(I*x).subs(x,pi).evalf()

for i in range(10):
    ci=sympy.Symbol("c%d%s"%(i,''),commutative=False)
    cdi=sympy.Symbol("cd%d%s"%(i,''),commutative=False)
    print ci,cdi