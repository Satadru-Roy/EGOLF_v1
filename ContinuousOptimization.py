# from __future__ import print_function
from openmdao.api import IndepVarComp, Component, Problem, Group, ScipyOptimizer
from Initialization import initialize_cont_test
import numpy as np
from scipy.optimize import minimize
import sys
import time

class cont_opt(Component):
    """ This is the Step 2 of EGOLF. Performs optimization only w.r.t
    the continuous design variables.The integer/discrete type design variables
    are supplied as parameters in this step. """

    def __init__(self,num_des,prob):
        super(cont_opt, self).__init__()

        self.add_param('xI', val=np.zeros([1,num_des/2]))
        self.add_param('xC', val=np.zeros([1,num_des/2]))
        self.add_output('f', val=0.0)
        self.nfev = 0
        self.num_des = num_des
        self.prob = prob

    def solve_nonlinear(self, params, unknowns, resids):
        """ Define the function f(xI, xC)
        Here xI is integer and xC is continuous"""
        self.nfev += 1
        if self.prob == 1:  #Branin function
            x0 = params['xI']
            x1 = params['xC']
            x = np.array([x0,x1])
            a = 1.0; b = 5.1/(4.0*np.pi**2); c = 5.0/np.pi; d = 6.0; e = 10.0; f = 1.0/(8.0*np.pi)
            unknowns['f'] = a*(x[1] - b*x[0]**2 + c*x[0] - d)**2 + e*(1-f)*np.cos(x[0]) + e
        elif self.prob == 2:    #Griewank function
            xI = params['xI']
            xC = params['xC']
            f1C=0.0; f1I=0.0; f2C=1.0; f2I=1.0
            f1I = np.sum((xI**2/4000.0))
            f1C = np.sum((xC**2/4000.0))
            for ii in xrange(len(xC)):
                f2C *= np.cos(xC[ii]/np.sqrt(ii+1.))
            for ii in xrange(len(xI)):
                f2I *= np.cos(xI[ii]/np.sqrt(ii+1.))
            unknowns['f'] = ((f1C+f1I)-(f2C*f2I) + 1.)[0]
        elif self.prob == 3:  #Rosenbrock's function
            x = np.array([params['xI'],params['xC']])
            f = 0.0
            for ii in xrange(self.num_des-1):
                f += 100*(x[ii+1] - x[ii]**2)**2 + (x[ii] - 1)**2
            unknowns['f'] = f

    def linearize(self, params, unknowns, resids):
        """ Provide the Jacobian"""
        J = {}
        if self.prob == 1:
            x0 = params['xI']
            x1 = params['xC']
            x = np.array([x0,x1])
            a = 1.; b = 5.1/(4.*np.pi**2); c = 5.0/np.pi; d = 6.; e = 10.0; f = 1./(8.*np.pi)
            J['f', 'xC'] = 2.0*a*(x[1] - b*x[0]**2 + c*x[0] - d)
            # J['f', 'xI'] = 2.0*a*(x[1] - b*x[0]**2 + c*x[0] - d)*(-2.*b*x[0] + c) - e*(1.-f)*np.sin(x[0])
        return J

def continuous_optimization_test(x0I, M, num_des, prob):

    n = x0I.shape[0]
    num_xI = x0I.shape[1]
    [xC_lb, xC_ub] = initialize_cont_test(num_des, prob)
    num_xC = len(xC_lb)
    xC0 = (xC_lb + 0.0*(xC_ub-xC_lb)).reshape(1,num_xC)
    xC_opt = np.zeros([n,num_xC])
    obj = np.zeros([n,1])
    funCount = np.ones([n,1])
    eflag = np.ones([n,1])

    for ii in xrange(n):
        x0I_val = x0I[ii,:].reshape(1,num_xC)
        top = Problem()
        root = top.root = Group()
        root.add('Inp1', IndepVarComp('xI', x0I_val))
        root.add('Inp2', IndepVarComp('xC', xC0))
        root.add('copt', cont_opt(num_des, prob))

        root.connect('Inp1.xI', 'copt.xI')
        root.connect('Inp2.xC', 'copt.xC')

        top.driver = ScipyOptimizer()
        top.driver.options['optimizer'] = 'SLSQP'

        top.driver.add_desvar('Inp2.xC', lower=xC_lb, upper=xC_ub)
        top.driver.add_objective('copt.f')

        top.root.deriv_options['type'] = 'fd'

        top.setup()
        top.run()
        # data = top.check_partial_derivatives(out_stream=sys.stdout)
        print "Minimum found f = ",top['copt.f']
        print "at xC = ",top['copt.xC']
        print "for the given xI = ", top['copt.xI']
        xC_opt[ii] = top['copt.xC']
        obj[ii] = top['copt.f']
        eflag[ii] = top.driver.exit_flag
        funCount[ii] = root.copt.nfev

    g = []
    return xC_opt, obj, g, eflag, funCount

#Testing module to perform optimization outside OpenMDAO
def contopt_test_outOpenMDAO(x0I, M, num_des, prob):
    n = x0I.shape[0]
    [xC_lb, xC_ub] = initialize_cont_test(num_des, prob)
    xC0 = xC_lb + 0.5*(xC_ub-xC_lb)
    num_xC = len(xC_lb)
    xC_opt = np.zeros([n,num_xC])
    obj = np.zeros([n,1])
    eflag = np.ones([n,1])
    g = []
    funCount = np.zeros([n,1])
    bnds = [(xC_lb[ii], xC_ub[ii]) for ii in xrange(num_xC)]
    for ii in xrange(n):
        print "=============================================="
        print "Continuous optimization:"
        optResult = minimize(obj_func,xC0,\
        args=(x0I[ii,:],num_des, prob), method='SLSQP',\
        bounds=bnds, options={})
        xC_opt[ii] = optResult.x.reshape(1,num_xC)
        obj[ii] = optResult.fun
        funCount[ii] = optResult.nfev
        if not optResult.success:
            eflag[ii]=0
        print("Minimum found f = %f" %(obj[ii]))
        print "at xC = ", xC_opt[ii,:]
        print "for the given xI = ", x0I[ii,:]
        print "=============================================="
    return xC_opt, obj, g, eflag, funCount

def obj_func(xC, xI, num_des, probs):
    if probs == 1: #Branin
        x = np.array([xI,xC]).reshape(num_des,1)
        a = 1.0; b = 5.1/(4.0*np.pi**2); c = 5.0/np.pi; d = 6.0; e = 10.0; f = 1.0/(8.0*np.pi)
        fval = a*(x[1] - b*x[0]**2 + c*x[0] - d)**2 + e*(1.0-f)*np.cos(x[0]) + e
    elif probs == 2: #Griewank
        f1C=0.0; f1I=0.0; f2C=1.0; f2I=1.0
        f1I = np.sum((xI**2/4000.0))
        f1C = np.sum((xC**2/4000.0))

        for ii in xrange(len(xC)):
            f2C *= np.cos(xC[ii]/np.sqrt(ii+1.))
        for ii in xrange(len(xI)):
            f2I *= np.cos(xI[ii]/np.sqrt(ii+1.))

        fval = (f1C+f1I)-(f2C*f2I) + 1.
    return fval

if __name__ =="__main__":
    xC = np.array([[0.0],[0.0]])
    xI = np.array([[0.0],[0.0]])
    fval = obj_func(xC,xI,4,2)
    print fval
