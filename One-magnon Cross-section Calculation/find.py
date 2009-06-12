from sympy import *

def find(expr, pattern, match=True, _found=None):
    """Find all instances of pattern in expr.
    If match is True (default) perform full pattern matching else use only
    isinstance()."""
    if match:
        check = lambda expr, pattern: expr.match(pattern)
    else:
        check = lambda expr, pattern: isinstance(expr, pattern)
    found = []
    _walk_expr(expr, pattern, found, check)
    return found

def _walk_expr(expr, pattern, found, check):
    """Helper function for find()"""
    # iterate over all subexpressions
    if expr.args:
        for subexpr in expr.args:
            if check(subexpr, pattern):
                found.append(subexpr)
            _walk_expr(subexpr, pattern, found, check)
    return found

def test_find():
    x = Symbol('x')
    e = (log(tan(x**2)) - 1)/(tan(tan(x - 1)) + cos(tan(1)))
    assert set(find(e, tan, match=False)) == set([tan(tan(1 - x)), tan(1 - x),
                                                  tan(1), tan(x**2)])
    a = Wild('a')
    print "find test",find(e, log(a))

test_find()