import sympy as sp

# Method for printing lists neatly
def list_print(lista):
    print 'printing...'
    for element in lista:
        print element
    print ''

# Method for multiplying two lists together
def list_mult(lista, listb):
    "Defines a way to multiply two lists of the same length"
    if len(lista) != len(listb):
        print "lists not same length"
        return []
    else:
        temp = []
        for i in range(len(lista)):
            if isinstance(lista[i], int):
                temp.append((lista[i] * listb[i].expand()).expand())
            elif isinstance(listb[i], int):
                temp.append((lista[i].expand() * listb[i]).expand())
            else: temp.append((lista[i].expand() * listb[i].expand()).expand())
        return temp

# Method for adding two lists together
def list_sum(lista, listb):
    "Defines a way to add two lists of the same length"
    if len(lista) != len(listb):
        print "lists not same length"
        return []
    else:
        temp = []
        for i in range(len(lista)):
            temp.append((lista[i].expand() + listb[i].expand()).expand())
        return temp

# Method for multiplying a list by a scalar
def scalar_mult(scalar, alist):
    "Defines a way to multiply a list by a scalar"
    temp = []
    for i in range(len(alist)):
        temp.append(scalar * alist[i])
    return temp

# Method that returns the # of b/bd operators in an expression
def b_scanner(expr):
    """Finds the number of b and b dagger operators in an expression"""
    b = 0; bd = 0
    s = str(expr)
    indexbd = 0
    while indexbd < len(s):
        indexbd = s.find('bd', indexbd + 1)
        if indexbd > 0: bd = bd + 1
        if indexbd < 0: break
    indexb = 0
    while indexb < len(s):
        indexb = s.find('b', indexb + 1)
        if indexb > 0: b = b + 1
        if indexb < 0: break
    return (b - bd, bd)

# Method to grab coefficients (Ondrej Certik)
def coeff(expr, term):
    expr = sp.collect(expr, term)
    symbols = list(term.atoms(sp.Symbol))
    w = sp.Wild("coeff", exclude=symbols)
    m = expr.match(w*term+sp.Wild("rest"))
    if m:
        return m[w] 

