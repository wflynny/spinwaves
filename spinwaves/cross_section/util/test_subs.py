import sympy as sp
from sympy import exp,log,sin,cos,oo,Sum,sympify,Symbol,symbols,Wild,simplify,Rational

from sub_in import sub_in

from sympy import legendre, Symbol, hermite, chebyshevu, chebyshevt, \
        chebyshevt_root, chebyshevu_root, assoc_legendre, Rational,  \
        roots, sympify, S

from sympy import Rational, sqrt, symbols, sin, exp, log, sinh, cosh, cos, pi, \
        I, S, erf, tan, asin, asinh, Function, Derivative, diff, trim
from sympy.integrals.risch import heurisch, components
from sympy.utilities.pytest import XFAIL, skip

x, y, z = symbols('xyz')
f = Function('f')

def test_chebyshev():
    assert chebyshevt(0, x) == 1
    assert chebyshevt(1, x) == x
    assert chebyshevt(2, x) == 2*x**2-1
    assert chebyshevt(3, x) == 4*x**3-3*x
    for n in range(1, 4):
        for k in range(n):
            z = chebyshevt_root(n, k)
            assert simplify(chebyshevt(n, z)) == 0
    for n in range(1, 4):
        for k in range(n):
            z = chebyshevu_root(n, k)
            assert simplify(chebyshevu(n, z)) == 0
            
test_chebyshev()
print 'pass'

def test_heurisch_radicals():
    assert heurisch(x**Rational(-1,2), x) == 2*x**Rational(1,2)
    assert heurisch(x**Rational(-3,2), x) == -2*x**Rational(-1,2)
    assert heurisch(x**Rational(3,2), x) == 2*x**Rational(5,2) / 5

    assert heurisch(sin(x)*sqrt(cos(x)), x) == -2*cos(x)**Rational(3,2) / 3
    assert heurisch(sin(y*sqrt(x)), x) == 2*y**(-2)*sin(y*x**S.Half) - \
                                          2*x**S.Half*cos(y*x**S.Half)/y

#test_heurisch_radicals()

def test_heurisch_hacking():
    assert heurisch(sqrt(1 + 7*x**2), x, hints=[]) == \
        x*sqrt(1+7*x**2)/2 + sqrt(7)*asinh(sqrt(7)*x)/14
    assert heurisch(sqrt(1 - 7*x**2), x, hints=[]) == \
        x*sqrt(1-7*x**2)/2 + sqrt(7)*asin(sqrt(7)*x)/14

    assert heurisch(1/sqrt(1 + 7*x**2), x, hints=[]) == \
        sqrt(7)*asinh(sqrt(7)*x)/7
    assert heurisch(1/sqrt(1 - 7*x**2), x, hints=[]) == \
        sqrt(7)*asin(sqrt(7)*x)/7

    assert heurisch(exp(-7*x**2),x,hints=[]) == \
        sqrt(7*pi)*erf(sqrt(7)*x)/14

    assert heurisch(1/sqrt(9 - 4*x**2), x, hints=[]) == \
        asin(2*x/3)/2

    assert heurisch(1/sqrt(9 + 4*x**2), x, hints=[]) == \
        asinh(2*x/3)/2

#test_heurisch_hacking()

def test_issue953():
    f = S(1)/2*asin(x) + x*(1 - x**2)**(S(1)/2)/2

    assert integrate(cos(asin(x)), x) == f
    assert integrate(sin(acos(x)), x) == f
    
#test_issue953()

def subs_test():
    # Define symbols
    a,b,c,d,K = sp.symbols('abcdK', commutative = True)
    w,x,y,z,Q = sp.symbols('wxyzQ', commutative = False)
   
    """ COMMUTATIVE TESTS """
    print '\n','COMMUTATIVE TESTS'
    assert (a*b    ).subs(a*b,K) == K,'Failed'; print '.'
    assert (a*b*a*b).subs(a*b,K) == K**2,'Failed'; print '.'
    assert (a*b*c*d).subs(a*b*c,K) == d*K,'Failed'; print '.'
    assert (a*b**c ).subs(a,K) == K*b**c, 'Failed'; print '.'
    assert (a*b**c ).subs(b,K) == a*K**c, 'Failed'; print '.'
    assert (a*b**c ).subs(c,K) == a*b**K, 'Failed'; print '.'
    assert (a*b*c*b*a  ).subs(a*b,K) == c*K**2,'Failed'; print '.'
    assert (a**3*b**2*a).subs(a*b,K) == a**2*K**2,'Failed'; print '.'
    

    """ NONCOMMUTATIVE TESTS """
    print '\n','NONCOMMUTATIVE TESTS'
    assert (x*y    ).subs(x*y,Q) == Q,'Failed'; print '.'
    assert (w*y*x  ).subs(x*y,Q) == w*y*x,'Failed'; print '.'
    assert (w*x*y*z).subs(x*y,Q) == w*Q*z,'Failed'; print '.'
    assert (x*y*x*y).subs(x*y,Q) == Q**2,'Failed'; print '.'
    assert (x*x*y  ).subs(x*y,Q) == x*Q,'Failed'; print '.'
    assert (x*x*y*y).subs(x*y,Q) == x*Q*y,'Failed'; print '.'
    assert (w*x*y  ).subs(x*y*z,Q) == w*x*y,'Failed'; print '.'
    assert (x*y**z     ).subs(x,Q) == Q*y**z, 'Failed'; print '.'
    assert (x*y**z     ).subs(y,Q) == x*Q**z, 'Failed'; print '.'
    assert (x*y**z     ).subs(z,Q) == x*y**Q, 'Failed'; print '.'
    assert (w*x*y*z*x*y).subs(x*y*z,Q) == w*Q*x*y,'Failed'; print '.'
    assert (w*x*y*y*w*x*x*y*x*y*y*x*y).subs(x*y,Q) == w*Q*y*w*x*Q**2*y*Q,'Failed'; print '.'
    
    """ OTHER OPERATION TESTS"""
    print '\n','OTHER OPERATION TESTS'
    assert (x+y  ).subs(x+y,Q) == Q,'Failed'; print '.'
    assert (x-y  ).subs(x-y,Q) == Q,'Failed'; print '.'
    assert (x/y  ).subs(x,Q) == Q/y,'Failed'; print '.'
    assert (x**y ).subs(x,Q) == Q**y,'Failed'; print '.'
    assert (x**y ).subs(y,Q) == x**Q,'Failed'; print '.'
    assert ( (a-c)/b  ).subs(b,K) == (a-c)/K,'Failed';print '.'
    assert (exp(x*y-z)).subs(x*y,Q) == exp(Q-z),'Failed'; print '.'
    assert (a*exp(x*y-w*z)+b*exp(x*y+w*z)).subs(z,0) == a*exp(x*y)+b*exp(x*y), 'Failed'; print '.'
    assert ((a-b)/(c*d-a*b)).subs(c*d-a*b,K) == (a-b)/K,'Failed'; print '.'
    assert (w*exp(a*b-c)*x*y/4).subs(x*y,Q) == w*exp(a*b-c)*Q/4,'Failed'; print '.'
    #assert (a/(b*c)).subs(b*c,K) == a/K,'Failed'; print '.' #FAILS DIVISION
    
    g = a*x + b/y - K**Q
    print '\nPrinting test\n', g,'\n', str(g)

    """ WILD TESTS """
    R = sp.Wild('R'); S = sp.Wild('S'); T = sp.Wild('T'); U = sp.Wild('U')
    print '\n','WILD TESTS'
    assert (R*S ).subs(R*S,T) == T,'Failed'; print '.'
    assert (S*R ).subs(R*S,T) == T,'Failed'; print '.'
    assert (R+S ).subs(R+S,T) == T,'Failed'; print '.'
    assert (R**S).subs(R,T) == T**S,'Failed'; print '.'
    assert (R**S).subs(S,T) == R**T,'Failed'; print '.'
    assert (R*S**T).subs(R,U) == U*S**T, 'Failed'; print '.'
    assert (R*S**T).subs(S,U) == R*U**T, 'Failed'; print '.'
    assert (R*S**T).subs(T,U) == R*S**U, 'Failed'; print '.'

    """ MIXED TESTS """ 
    print '\n','MIXED TESTS'
    assert (    a*x*y     ).subs(x*y,Q) == a*Q,'Failed'; print '.'
    assert (  a*b*x*y*x   ).subs(x*y,Q) == a*b*Q*x,'Failed'; print '.'
    assert (R*x*y*exp(x*y)).subs(x*y,Q) == R*Q*exp(Q),'Failed'; print '.'
    assert (     a*x*y*y*x-x*y*z*exp(a*b)   ).subs(x*y,Q) == a*Q*y*x-Q*z*exp(a*b),'Failed'; print '.'
    assert (c*y*x*y*x**(R*S-a*b)-T*(a*R*b*S)).subs(x*y,Q).subs(a*b,K).subs(R*S,U) == c*y*Q*x**(U-K)-T*(U*K),'Failed'; print '.'
    
    
if __name__ == '__main__':

    subs_test()
    

    