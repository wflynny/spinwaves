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
                temp.append(sp.powsimp((lista[i] * listb[i]).expand(deep=False)))
            elif isinstance(listb[i], int):
                temp.append(sp.powsimp((lista[i] * listb[i]).expand(deep=False)))
            else:
                temp.append(sp.powsimp((lista[i] * listb[i]).expand(deep=False)))
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

## Method that returns the # of b/bd operators in an expression
#def b_scanner(expr):
#    """Finds the number of b and b dagger operators in an expression"""
#    b = 0; bd = 0
#    s = str(expr)
#    indexbd = 0
#    while indexbd < len(s):
#        indexbd = s.find('bd', indexbd + 1)
#        if indexbd > 0:
#            marker = s[indexbd+3:indexbd+5]
#            if marker == '**':
#                expnum = s[indexbd+5]
#                b = b + int(expnum)
#            else: bd = bd + 1
#        if indexbd < 0: break
#    indexb = 0
#    while indexb < len(s):
#        indexb = s.find('b', indexb + 1)
#        if indexb > 0:
#            marker = s[indexb+2:indexb+4]
#            if marker == '**':
#                expnum = s[indexb+4]
#                b = b + int(expnum)
#            else: b = b + 1        
#        if indexb < 0: break
#    return (b - bd, bd)

# Method to grab coefficients (Ondrej Certik)
def coeff(expr, term):
    if isinstance(expr, int):
        return 0
    expr = sp.collect(expr, term)
    #print 'expr',expr
    symbols = list(term.atoms(sp.Symbol))
    #print 'symbols',symbols
    w = sp.Wild("coeff", exclude = symbols)
    #print 'w',w
    m = expr.match(w * term + sp.Wild("rest"))
    #print 'm',m
    m2 = expr.match(w * term)
    #print 'm2',m2
    res = False
    if m2 != None:
        #print 'm2[w]',m2[w]
        res = m2[w] * term == expr
    if m and res!= True:
        return m[w]
    #added the next two lines
    elif m2:
        return m2[w]

