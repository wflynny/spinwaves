from sympy import *
from sympy.printing import *


#-------------------------------------------------------------------------------


"""
NOTICE:
PLEASE USE AS FEW WILDS AS NEEDED. THE MORE WILDS, THE (MUCH) MORE TIME NEEDED.
For example, 

e = exp(I*x*t + I*y*t - I*x*y + I*z + I*x) --> (x*y)*t + rest
you only need 3 wilds

A,B,C be Wilds

sub_in (  e  , exp(A*I*t + B*I*t + C)  , (A*B)*t + C )

"""

def sub_in(expression, pattern, replacement, match = True):
    """ Searches through the expression looking for the pattern and will
    replace the pattern with replacement if it finds it. Use wilds to write
    your pattern and replacement. Ex:
      x,y = symbols('xy'); j = Wild('j'); k = Wild('k')
      sub_in(exp(x*y), exp(j*k), cos(j-k))
    """

    # Might take out match optionality because will always need it.
    if match:
        check = lambda expression, pattern: expression.match(pattern)
    else:
        check = lambda expression, pattern: isinstance(expression, pattern)
    new = _walk_it(expression, pattern, replacement, check)
    if new != None: 
        return new
    else: return None

def _walk_it(expression, pattern, replacement, check):
    """ Helper method for sub_in """
    # Grab the operation of the expression
    expression = expression
    op = expression.__class__

    """ If the length of the expression and pattern are the same, then
    just scan through to make sure everything in the pattern is present in the
    expression. If so, then replace the wilds in the replacement with the
    variables they correspond to in the expression"""
    if len(expression.args) == len(pattern.args):
        if isinstance(type(expression), FunctionClass):
            new = [expression]; v = []
            if check(expression, pattern) != None:
                ch = list(check(expression, pattern).iteritems())
                for i in ch: i = list(i); i[1] = simplify(i[1]); v.append(i)
                new.insert(new.index(expression), replacement.subs(v))
                new.remove(expression)
                return Mul(*new)
        else:
            expr_terms = [sub for sub in expression.args]
            pat_terms = [sub for sub in pattern.args]
            new = [expression]; v = []
            allThere = True
            for i in range(len(expr_terms)):
                if check(expr_terms[i], pat_terms[i]) == None: allThere = False; break
            if allThere == True and check(expression, pattern) != None:
                ch = list(check(expression, pattern).iteritems())
                for i in ch: i = list(i); i[1] = simplify(i[1]); v.append(i)
                new.insert(new.index(expression), replacement.subs(v))
                new.remove(expression)
                return Mul(*new)

#    """ If the expression is just a function (generally just like exp(blah blah))
#    then check if the whole expression matches the pattern. If so, replace the
#    wilds in the replacement with the variables they correspond to in the
#    expression"""
    if isinstance(type(expression),FunctionClass):
        new = [expression]; v = []
        if check(expression,pattern) != None:
            ch = list(check(expression, pattern).iteritems())
            for i in ch: i = list(i); i[1] = simplify(i[1]); v.append(i)
            new.insert(new.index(expression), replacement.subs(v))
            new.remove(expression)
            return Mul(*new)
        elif expression.args:
            new = [subexpression for subexpression in expression.args]; v = []
            for sub in new:
                if check(sub,pattern) != None:
                    ch = list(check(sub, pattern).iteritems())
                    for i in ch: i = list(i); i[1] = simplify(i[1]); v.append(i)
                    new.insert(new.index(sub), replacement.subs(v))
                    new.remove(sub)
                else:
                    new.insert(new.index(sub), _walk_it(sub, pattern, replacement, check))
                    new.remove(sub)
            return op(*new)

#    """ Else if the expression has multiple arguments, scan through each. """
    elif expression.args:
        new = [subexpression for subexpression in expression.args]; v = []
        for sub in new:
            if check(sub, pattern) != None:
                ch = list(check(sub, pattern).iteritems())
                for i in ch: i = list(i); i[1] = simplify(i[1]); v.append(i)
                new.insert(new.index(sub), replacement.subs(v))
                new.remove(sub)
            else:
                new.insert(new.index(sub), _walk_it(sub, pattern, replacement, check))
                new.remove(sub)
        return op(*new)
    else: return expression


#-------------------------------------------------------------------------------


def test_sub_in():
    # Define symbols, wilds, etc. 
    a,b,c,d = symbols('abcd', commmutative = True)
    t,x,y,z = symbols('txyz', commutative = False)
    j = Wild('j'); k = Wild('k'); l = Wild('l');
    F = WildFunction('f')

#    # Demonstration
#    e = a*exp(x*y*z-y*x*t) + b*exp(x**2*z-y**2*t) + c*exp(y**2*z-x**2*t) + d*exp(y*x*z-x*y*t)
#    print e
#    print sub_in(e,exp(j*z-k*t),DiracDelta(j-k))
#    print sub_in(e,exp(j*z-k*t),DiracDelta(j*z-k*t))
#    print


    # Tests
    assert sub_in(a*exp(x), exp(j), DiracDelta(j)) == a*DiracDelta(x), 'Failed'; print 'Pass'
    assert sub_in(a*exp(x)+exp(y), exp(j), DiracDelta(j)) == a*DiracDelta(x)+DiracDelta(y), 'Failed'; print 'Pass'
    assert sub_in(a*exp(x*y), exp(j*k), DiracDelta(j-k)) == a*DiracDelta(x-y), 'Failed'; print 'Pass'
    assert sub_in(a*exp(x*y*z), exp(j*z-k*t), DiracDelta(j-k)) == a*exp(x*y*z), 'Failed'; print 'Pass'
    assert sub_in(a*exp(x*y) + b*exp(t*z), exp(j*k), cos(j-k)) == a*cos(x-y) + b*cos(t-z), 'Failed'; print 'Pass'
    assert sub_in(a*exp(x*y*z - y*x*t)*cos(x), exp(j*z-k*t), DiracDelta(j-k)) == a*DiracDelta(x*y-y*x)*cos(x), 'Failed'; print 'Pass'
    # Trig Tests
    assert sub_in(sin(x)/cos(x), sin(j)/cos(j), tan(j)) == tan(x), 'Failed'; print 'Pass'
    assert sub_in(sin(x)*cos(x), sin(j)/cos(j), tan(j)) != tan(x), 'Failed'; print 'Pass'
    assert sub_in(cos(x)/sin(x), cos(k)/sin(k), cot(k)) == cot(x), 'Failed'; print 'Pass'
    assert sub_in(cos(x)*sin(x), cos(k)/sin(k), cot(k)) != cot(x), 'Failed'; print 'Pass'
    assert sub_in(sin(x)/cos(x)+cos(x), sin(j)/cos(j), tan(j)) == tan(x)+cos(x), 'Failed'; print 'Pass'
    assert sub_in(cos(x-y), j-k, j) != cos(x), 'Passed'; print 'Failed'
   

#-------------------------------------------------------------------------------


# Main
if __name__ == '__main__':

    test_sub_in()