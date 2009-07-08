from sympy import *
from sympy.printing import *

#-------------------------------------------------------------------------------

def subin(expression, pattern, replacement, match=True, *vars):
    # Might take out match optionality because will always need it.
    if match:
        check = lambda expression, pattern: expression.match(pattern)
    else:
        check = lambda expression, pattern: isinstance(expression,pattern)
    new = _walk_it(expression, pattern, replacement, check, vars)
    if new != None: return new
    else: return None

def _walk_it(expression, pattern, replacement, check, vars):
    op = expression.__class__
    if isinstance(type(expression),FunctionClass):
        new = [expression]; v = []
        if check(expression,pattern) != None:
            ch = list(check(expression,pattern).iteritems())
            for i in ch: v.append(i)
            new.insert(new.index(expression),replacement.subs(v))
            new.remove(expression)
        return Mul(*new)
    elif expression.args:
        new = [subexpression for subexpression in expression.args]; v = []
        for sub in new:
            if check(sub,pattern) != None:
                ch = list(check(sub,pattern).iteritems())
                for i in ch: v.append(i)
                new.insert(new.index(sub),replacement.subs(v))
                new.remove(sub)
            else: 
                new.insert(new.index(sub),_walk_it(sub, pattern, replacement, check, vars))
                new.remove(sub)
        return op(*new)
    else: return expression

#-------------------------------------------------------------------------------

def test_subin():
    a,b,c,d = symbols('abcd', commmutative = True)
    t,x,y,z = symbols('txyz', commutative = False)
    j = Wild('j'); k = Wild('k'); l = Wild('l');
    F = WildFunction('f')
    
    assert subin(a*exp(x*y), exp(j*k), DiracDelta(j-k)) == a*DiracDelta(x-y), 'Failed'; print '.'
    assert subin(a*exp(x*y) + b*exp(t*z), exp(j*k), cos(j-k)) == a*cos(x-y) + b*cos(t-z), 'Failed'; print '.'
    assert subin(a*exp(x*y*z - y*x*t)*cos(x), exp(j*z-k*t), DiracDelta(j-k)) == a*DiracDelta(x*y-y*x)*cos(x), 'a*exp(x*y*z - y*x*t)*cos(x) != a*DiracDelta(x*y-y*x)*cos(x)'; print '.'

#-------------------------------------------------------------------------------

if __name__ == '__main__':

    test_subin()



