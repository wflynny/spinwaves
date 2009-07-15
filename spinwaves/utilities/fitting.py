import vtkModel.Parameter_Manager as PM
from mpfit import mpfit
import numpy as np

def fitting(session, spinwave_domain = [], size=2, k = 100, tFactor = .95):
    """
    DOCUMENTATIONS
    Probably all wrong
    """

    # Create fitter
    fit = PM.Fitter(session, spinwave_domain, size, k, tFactor)
    
    # Helper function for GetResult
    def myfunc(fitter, fjac=None, x=None, y=None, err=None):
        model = fitter.GetResult().T
        status = 0
        return [status, (y-model)/err]

    # Function Keywords
    
    #spinwave_domain is a list of 3D lists:
    #-[(x1,y1,z1),(x2,y2,z2),...]
    #x must be a 1D numpy array:
    #[x1,y1,z1,x2,t2,z2,...]
    xList = []
    for point in fit.spinwave_domain:
        xList.append(point[0])
        yList.append(point[1])
        zList.append(point[2])
    x = np.array(xList)
    y = myfunc(fit)
    err = 0.01*y
    fa = {'x':x, 'y':y, 'err':err}
    
    # Set parinfo with the limits and such. Probably don't need
    parbase={'value':0., 'limited':[0,0], 'limits':[0.,0.]}
    parinfo=[]
    p0 = fit.fitlist
    for i in range(len(p0)):
        parinfo.append(copy.deepcopy(parbase))
    for i in range(len(p0)):
#        parinfo[i]['value']=x[i]
        if(fit.min_range_list[i] != np.NINF):
            parinfo[i]['limits'][0]=fit.min_range_list[i]
            parinfo[i]['limited'][0] = 1
        else:
            parinfo[i]['limited'][0] = 0
        if fit.max_range_list[i] != np.PINF:
            parinfo[i]['limits'][1] = fit.max_range_list[i]
            parinfo[i]['limited'][1] = 1
        else:
            parinfo[i]['limited'][1] = 0

    # Run mpfit on fitlist with all the jazz. 
    m = mpfit(myfunc, fit.fitlist, parinfo=parinfo, functkw = fa)
    # Bare Bones Run for when the above (especially keywords x) breaks instantly
    #m = mpfit(myfunc, fit.fitlist)
    
    #Return the parameters
    return (m.params, m.perror)
    
if '__name__' == '__main__':
    print 'poop'
