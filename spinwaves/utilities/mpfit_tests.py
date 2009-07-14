from numpy.testing import *
import numpy as N
import copy

from mpfit import mpfit

def Flin(x,p):
    y =  p[0] -p[1]*x 
    return y

def Fquad(x,p):
    y = p[0] + p[1]*x + p[2]*x*x
    return y

def Fgauss(x,p):
    xc = x - p[2]
    sig2 = p[3]*p[3]
    y = p[0] + p[1] * N.exp(-0.5*xc*xc / sig2)
    return y

def myfunctlin(p, fjac=None, x=None, y=None, err=None):
    # Parameter values are passed in "p"
    # If fjac==None then partial derivatives should not be
    # computed.  It will always be None if MPFIT is called with default
    # flag.
    model = Flin(x, p)
    # Non-negative status value means MPFIT should continue, negative means
    # stop the calculation.
    status = 0
    return [status, (y-model)/err]

def myfunctquad(p, fjac=None, x=None, y=None, err=None):
    model = Fquad(x,p)
    status = 0
    return [status, (y-model)/err]

def myfunctgauss(p, fjac=None, x=None, y=None, err=None):
    model = Fgauss(x, p)
    status = 0
    return [status, (y-model)/err]

def test_linfit():
    x=N.array([-1.7237128E+00,1.8712276E+00,-9.6608055E-01,
        -2.8394297E-01,1.3416969E+00,1.3757038E+00,
        -1.3703436E+00,4.2581975E-02,-1.4970151E-01,
        8.2065094E-01])
    y=N.array([1.9000429E-01,6.5807428E+00,1.4582725E+00,
        2.7270851E+00,5.5969253E+00,5.6249280E+00,
        0.787615,3.2599759E+00,2.9771762E+00,
        4.5936475E+00])
    ey=0.07*N.ones(y.shape,dtype='float64')
    p0=N.array([1.0,1.0],dtype='float64')  #initial conditions
    pactual=N.array([3.2,1.78]) #actual values used to make data
    parbase={'value':0., 'fixed':0, 'limited':[0,0], 'limits':[0.,0.]}
    parinfo=[]
    for i in range(len(pactual)):
        parinfo.append(copy.deepcopy(parbase))
    for i in range(len(pactual)): 
        parinfo[i]['value']=p0[i]
    fa = {'x':x, 'y':y, 'err':ey}
    m = mpfit(myfunctlin, p0, parinfo=parinfo,functkw=fa)
    if (m.status <= 0): 
        print 'status = ', m.status
        print 'params = ', m.params
    assert N.allclose(m.params,N.array([ 3.20996572, -1.7709542 ],dtype='float64'))
    assert N.allclose(m.perror,N.array([ 0.02221018,  0.01893756],dtype='float64'))
    chisq=(myfunctlin(m.params, x=x, y=y, err=ey)[1]**2).sum()
    
    assert N.allclose(N.array([chisq],dtype='float64'),N.array([2.756284983],dtype='float64'))
    assert m.dof==8
    return
#    assert zzz()=='Hello from zzz'

def test_quadfit():

    x = N.array([-1.7237128E+00,1.8712276E+00,-9.6608055E-01,
          -2.8394297E-01,1.3416969E+00,1.3757038E+00,
          -1.3703436E+00,4.2581975E-02,-1.4970151E-01,
          8.2065094E-01])
    y = N.array([2.3095947E+01,2.6449392E+01,1.0204468E+01,
          5.40507,1.5787588E+01,1.6520903E+01,
          1.5971818E+01,4.7668524E+00,4.9337711E+00,
          8.7348375E+00])
    ey=0.2*N.ones(y.shape,dtype='float64')
    p0=N.array([1.0,1.0,1.0],dtype='float64')  #initial conditions
    pactual=N.array([4.7, 0.0, 6.2]) #actual values used to make data
    parbase={'value':0., 'fixed':0, 'limited':[0,0], 'limits':[0.,0.]}
    parinfo=[]
    for i in range(len(pactual)):
        parinfo.append(copy.deepcopy(parbase))
        parinfo[i]['value'] = p0[i]
    fa = {'x': x, 'y': y, 'err': ey}
    m = mpfit(myfunctquad, p0, parinfo = parinfo, functkw = fa)
    if (m.status <= 0):
        print 'status = ', m.status
    print m.perror
    assert N.allclose(m.params, N.array([ 4.703829, 0.062586,  6.163087], dtype='float64'))
    assert N.allclose(m.perror, N.array([ 0.097512, 0.054802,  0.054433], dtype='float64'))
    chisq = (myfunctquad(m.params, x=x, y=y, err=ey)[1]**2).sum()
    assert N.allclose(N.array([chisq], dtype='float64'), N.array([5.679323], dtype='float64'))
    assert m.dof == 7
    return

def test_quadfit_fix():
    x = N.array([-1.7237128E+00,1.8712276E+00,-9.6608055E-01,
        -2.8394297E-01,1.3416969E+00,1.3757038E+00,
        -1.3703436E+00,4.2581975E-02,-1.4970151E-01,
        8.2065094E-01])
    y = N.array([2.3095947E+01,2.6449392E+01,1.0204468E+01,
        5.40507,1.5787588E+01,1.6520903E+01,
        1.5971818E+01,4.7668524E+00,4.9337711E+00,
        8.7348375E+00])
    ey=0.2*N.ones(y.shape,dtype='float64')
    p0=N.array([1.0, 0.0, 1.0],dtype='float64')  #initial conditions
    pactual=N.array([4.7, 0.0, 6.2]) #actual values used to make data
    parbase={'value':0., 'fixed':0, 'limited':[0,0], 'limits':[0.,0.]}
    parinfo=[]
    for i in range(len(pactual)):
        parinfo.append(copy.deepcopy(parbase))
        parinfo[i]['value'] = p0[i]
    
    parinfo[1]['fixed'] = 1

    fa = {'x': x, 'y': y, 'err': ey}
    m = mpfit(myfunctquad, p0, parinfo = parinfo, functkw = fa)
    if (m.status <= 0):
        print 'status = ', m.status
    assert N.allclose(m.params, N.array([ 4.696254, 0.      ,  6.172954], dtype='float64'))
    assert N.allclose(m.perror, N.array([ 0.097286, 0.      ,  0.053743], dtype='float64'))
    chisq = (myfunctquad(m.params, x=x, y=y, err=ey)[1]**2).sum()
    assert N.allclose(N.array([chisq], dtype='float64'), N.array([6.983588], dtype='float64'))
    assert m.dof == 8
    return

def test_gaussfit():
    x = N.array([-1.7237128E+00,1.8712276E+00,-9.6608055E-01,
        -2.8394297E-01,1.3416969E+00,1.3757038E+00,
        -1.3703436E+00,4.2581975E-02,-1.4970151E-01,
        8.2065094E-01]);
    y = N.array([-4.4494256E-02,8.7324673E-01,7.4443483E-01,
        4.7631559E+00,1.7187297E-01,1.1639182E-01,
        1.5646480E+00,5.2322268E+00,4.2543168E+00,
        6.2792623E-01]);
    ey=0.5*N.ones(y.shape,dtype='float64')
    p0=N.array([0.0, 1.0, 1.0, 1.0],dtype='float64')  #initial conditions
    pactual=N.array([0.0, 4.70, 0.0, 0.5]) #actual values used to make data
    parbase={'value':0., 'fixed':0, 'limited':[0,0], 'limits':[0.,0.]}
    parinfo=[]
    for i in range(len(pactual)):
        parinfo.append(copy.deepcopy(parbase))
        parinfo[i]['value'] = p0[i]
    fa = {'x': x, 'y': y, 'err': ey}
    m = mpfit(myfunctgauss, p0, parinfo = parinfo, functkw = fa)
    if (m.status <= 0):
        print 'status = ', m.status
    assert N.allclose(m.params, N.array([ 0.480443, 4.550752, -0.062562, 0.397472], dtype='float64'))
    assert N.allclose(m.perror, N.array([ 0.232235, 0.395434,  0.074715, 0.089996], dtype='float64'))
    chisq = (myfunctgauss(m.params, x=x, y=y, err=ey)[1]**2).sum()
    assert N.allclose(N.array([chisq], dtype='float64'), N.array([10.350032], dtype='float64'))
    assert m.dof == 6
    return

def test_gaussfit_fixed():

    x = N.array([-1.7237128E+00,1.8712276E+00,-9.6608055E-01,
        -2.8394297E-01,1.3416969E+00,1.3757038E+00,
        -1.3703436E+00,4.2581975E-02,-1.4970151E-01,
        8.2065094E-01])
    y = N.array([-4.4494256E-02,8.7324673E-01,7.4443483E-01,
        4.7631559E+00,1.7187297E-01,1.1639182E-01,
        1.5646480E+00,5.2322268E+00,4.2543168E+00,
        6.2792623E-01])
    ey=0.5*N.ones(y.shape,dtype='float64')
    p0=N.array([0.0, 1.0, 0.0, 0.1],dtype='float64')  #initial conditions          
    pactual=N.array([0.0, 4.70, 0.0, 0.5]) #actual values used to make data
    parbase={'value':0., 'fixed':0, 'limited':[0,0], 'limits':[0.,0.]}
    parinfo=[]
    for i in range(len(pactual)):
        parinfo.append(copy.deepcopy(parbase))
        parinfo[i]['value'] = p0[i]
    
    parinfo[0]['fixed'] = 1
    parinfo[2]['fixed'] = 1

    fa = {'x': x, 'y': y, 'err': ey}
    m = mpfit(myfunctgauss, p0, parinfo = parinfo, functkw = fa)
    if (m.status <= 0):
        print 'status = ', m.status
    assert N.allclose(m.params, N.array([ 0., 5.059244,  0.,  0.479746 ], dtype='float64'))
    assert N.allclose(m.perror, N.array([ 0., 0.329307,  0.,  0.053804 ], dtype='float64'))
    chisq = (myfunctgauss(m.params, x=x, y=y, err=ey)[1]**2).sum()
    assert N.allclose(N.array([chisq], dtype='float64'), N.array([15.516134], dtype='float64'))
    assert m.dof == 8
    return


if __name__ == "__main__":
    run_module_suite()
