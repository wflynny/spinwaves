import sympy as sp
from sympy import exp

def subs_test():
    # Define symbols
    t,w,x,y,z = sp.symbols('twxyz',commutative=False)
    # Start off easy
    a = x*y
    b = w*y*x
    # Add some noise
    c = w*x*y*z
    # Double
    d = x*y*x*y
    # Hard
    e = x*x*y
    # Hardest
    f = x*x*y*y
    # New Hardest
    g = w*x*y
    
    # Check other functions still work
    h = x+y
    i = x-y
    j = exp(x*y-z)
    k = x/y
    
    assert a.subs(x*y,t) == t; print 'test A passed'
    assert b.subs(x*y,t) == w*y*x; print 'test B passed'
    assert c.subs(x*y,t) == w*t*z; print 'test C passed'
    assert d.subs(x*y,t) == t**2; print 'test D passed'
    assert e.subs(x*y,t) == x*t; print 'test E passed'
    assert f.subs(x*y,t) == x*t*y; print 'test F passed'
    assert g.subs(x*y*z,t) == w*x*y; print 'test G passed'
    assert h.subs(x+y,t) == t; print 'test H passed'
    assert i.subs(x-y,t) == t; print 'test I passed'
    assert j.subs(x*y,t) == exp(t-z); print 'test J passed'
    assert k.subs(x,t) == t/y; print 'test K passed'
    print a,h; print 'printing functional' 
    
subs_test()

    