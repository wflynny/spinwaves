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
    x = fit.spinwave_domain
    y = myfunc(fit)
    err = 0.01*y
    fa = {'x':x, 'y':y, 'err':err}
    
    # Set parinfo with the limits and such. Probably don't need
    parbase={'value':0., 'limited':[0,0], 'limits':[0.,0.]}
    parinfo=[]
    for i in range(len(p0)):
        parinfo.append(copy.deepcopy(parbase))
    for i in range(len(p0)):
#        parinfo[i]['value']=x[i]
        parinfo[i]['limits'][0]=fit.min_range_list[i]
        parinfo[i]['limits'][1]=fit.max_range_list[i]
#        parinfo[i]['limited'][0]=fit.min_range_list[i]
#        parinfo[i]['limited'][1]=fit.min_range_list[i]

    # Run mpfit on fitlist with all the jazz. 
    m = mpfit(myfunc, fit.fitlist, parinfo=parinfo, functkw = fa)
    # Bare Bones Run for when the above (especially keywords x) breaks instantly
    #m = mpfit(myfunc, fit.fitlist)
    
    #Return the parameters
    return (m.params, m.perror)
    
if '__name__' == '__main__':
    print 'poop'
