import copy
import numpy as np
from numpy import pi
from mpfit import mpfit
import spinwaves.vtkModel.Parameter_Manager as PM
from spinwaves.vtkModel.wxGUI.Session import Session

def fitting(session, spinwave_domain, spinwave_range, size=3, k = 100, tFactor = .95):
    """
    DOCUMENTATIONS
    Probably all wrong
    """

    # Create fitter
    fitter = PM.Fitter(session, spinwave_domain, size, k, tFactor)
    
    # Helper function for GetResult
    def myfunc(p, fjac=None, y=None, err=None):
        fitter.fit_list = p
        model = fitter.GetResult()
        status = 0
        print 'y:\n', y, '\n\nmodel:\n', model
        #print '\n\ny-model:\n', (y-model)
        #print '\n\nerr:\n', err
        result = (y-model)/err
        print '\n\nresult:\n', result
        return [status, result]

    # Function Keywords
    
    #spinwave_domain is a list of 3D lists:
    #-[(x1,y1,z1),(x2,y2,z2),...]
    #x must be a 1D numpy array:
    #[x1,y1,z1,x2,t2,z2,...]
    #xList = []
    #for point in fit.spinwave_domain:
    #    xList.append(point[0])
    #    yList.append(point[1])
    #    zList.append(point[2])
    #x = np.array(xList)
    #y = myfunc(fit)
    y = spinwave_range
    errVal = .001
    err = [errVal]*len(y)
    #fa = {'x':x, 'y':y, 'err':err}
    fa = {'y':y, 'err':err}
    
    # Set parinfo with the limits and such. Probably don't need
    parbase={'value':0., 'limited':[0,0], 'limits':[0.,0.]}
    parinfo=[]
    p0 = fitter.fit_list
    for i in range(len(p0)):
        parinfo.append(copy.deepcopy(parbase))
    for i in range(len(p0)):
        parinfo[i]['value']=p0[i]
        if(fitter.min_range_list[i] != np.NINF):
            parinfo[i]['limits'][0]=fitter.min_range_list[i]
            parinfo[i]['limited'][0] = 1
        else:
            parinfo[i]['limited'][0] = 0
        if fitter.max_range_list[i] != np.PINF:
            parinfo[i]['limits'][1] = fitter.max_range_list[i]
            parinfo[i]['limited'][1] = 1
        else:
            parinfo[i]['limited'][1] = 0

    # Run mpfit on fitlist with all the jazz. 
    m = mpfit(myfunc, p0, parinfo=parinfo, functkw = fa)
    # Bare Bones Run for when the above (especially keywords x) breaks instantly
    #m = mpfit(myfunc, fit.fitlist)
    
    #Return the parameters
    return (m.status, m.params, m.perror)
    
if __name__ == '__main__':
      print 'start'
      sess = Session()
      domain = [(0,0,0,),
                (0,0,.25*pi),
                (0,0,.5*pi),
                (0,0,.75 *pi),
                (0,0, pi),
                (0,0,1.25 *pi),
                (0,0,1.5 *pi),
                (0,0,1.75 * pi)]
      sw_range = []
      for point in domain:
          sw_range.append(4 - 4*np.cos(point[2]))
      sw_range = np.array(sw_range)
      print '\n\nspinwave range:\n', sw_range
      #sw_range = (4 - 4*np.cos(domain[))
      sess.openXMLSession('C:\\testsess.xml')
      print 'fitting...'
      print fitting(sess, domain, sw_range, size = 3, k = 1000)
#===============================================================================
#    def simpleFunc(p, fjac = None, y = None, err = None):
#        model = (np.cos((p**.34)/15 + 4))
#        result = (y-model)/err
#        return [0, result]
#    
#    
#    p0 = np.array([0,0,0])
#    x = np.array([1.1,2.034,3.247])
#    y = (np.cos((x**.34)/15 + 4))
#    print 'y:',y
#    err = 3*[1E-4]
#    fa = {'y':y, 'err':err}
#    
#    m = mpfit(simpleFunc, p0, functkw = fa)
# 
#    print (m.status, m.params, m.perror)
#===============================================================================

    
    
