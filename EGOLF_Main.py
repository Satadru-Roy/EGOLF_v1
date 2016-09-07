''' This is the EGO-Like Framework (EGOLF) for the
simultaneous design-mission-allocation optimization problem.
Handles mixed-integer/discrete type design variables in a
computationally efficient manner and finds a near-global solution
to the above MINLP/MDNLP problem.

Developed by Satadru Roy
Purdue University, West Lafayette, IN
Copyright, July 2016'''

import numpy as np
import copy as cp
from Initialization import initialize_test, ModelInfo
from ContinuousOptimization import contopt_test_outOpenMDAO, continuous_optimization_test
from openmdao.api import KrigingSurrogate
from solver_MINLP import MINLP_BB
from random import uniform
import time

num_des = 2 #Total number of design variables
prob = 2 # Problem type: 1. Branin 2. Griewank 3. Rosenbrock 4. 3 Bar Truss problem
################################################################################
#Step 0: Initialize
iter = 1
ter_crit = 0;Tot_pt_prev = 0;y_opt = 0.0;num_pt = 0;fea_pt = 0;Tot_newpt_added = 0
Tot_FunCount = 0;ec2 = 0;ei_max = 0.0
ei_tol_per = 0.001
ei_tol_abs = 0.001
new_sol = [];comp = [];FEA_obj = [];FEA_xopt = []; ModelInfo_obj=[]
y_opt = np.inf

[xI_lb,xI_ub,M] = initialize_test(num_des, prob)
num_xI = len(xI_lb)
max_pt_lim = 10*num_xI  # Max infill points to be added
n = num_xI + 2 # Number of starting points
num_xC = num_des - num_xI

ModelInfo_obj = ModelInfo(xI_lb, xI_ub, num_xI)
ModelInfo_g=[[]]*M
if M>0:
    for mm in xrange(M):
        ModelInfo_g[mm] = ModelInfo(xI_lb, xI_ub, num_xI)

################################################################################
# Step 1: Generate a set of initial integer points
# Use Latin Hypercube Sampling to generate the initial points
# User supplied (in future use LHS). Provide num_xI+2 starting points
if (prob == 1 or prob == 2 or prob ==3) and num_des == 2:
    x0I_01 = np.array([[0.0],[0.76],[1.0]]) #2-D Problem with 1 int var.
elif prob == 4:
    x0I_01 = np.array([[0.0, 0.51, 0.75],[0.75, 1.0, 1.0],[0.0, 0.75, 0.0],[0.25, 0.0, 0.25],[1.0,0.25,0.51]]) #3bar Truss (obtained using LHS of matlab)

# Randomly generated initial integer points (Use LHS in the future release)
# x0I_01 = np.zeros([n,num_xI])
# for ii in xrange(n):
#     for jj in xrange(num_xI):
#         x0I_01[ii,jj] = uniform(0,1)

#####################################################
#Complete enumeration for verification
# x0I_01 = np.array([[0.0],[0.07],[0.13],[0.2],[0.27],[0.33],[0.4],[0.47],[0.53],[0.6],[0.67],[0.73],[0.8],[0.87],[0.93],[1.]])
# n=16
#####################################################

x0I = np.zeros([n,num_xI])
x0I_hat = np.zeros([n,num_xI])

for ii in xrange(n):
    x0I[ii] = np.round(xI_lb + np.array([x0I_01[ii]]).T*(xI_ub - xI_lb)).reshape(1,num_xI)
    x0I_hat[ii] = ((np.array([x0I[ii]]).T - xI_lb)/(xI_ub - xI_lb)).reshape(1,num_xI)

print "x0I", x0I
# print "x0I_hat", x0I_hat
time.sleep(3)
while ter_crit == 0:
    ############################################################################
    # Step 2: Perform the optimization w.r.t continuous design variables
    print "======================ContinuousOptimization-Start====================================="
    [xC_opt, obj, g, eflag, funCount] = continuous_optimization_test(x0I, M, num_des, prob)
    # [xC_opt, obj, g, eflag, funCount] = contopt_test_outOpenMDAO(x0I, M, num_des, prob)
    print "======================ContinuousOptimization-End======================================="
    ############################################################################
    # Step 3: Build the surrogate models
    for nonNAN in xrange(x0I.shape[0]):
        # Put a check here to ensure obj is never NaN or imaginary
        Tot_FunCount += funCount[nonNAN,0]
        num_pt += 1
        Tot_newpt_added += 1

        # Surrogate data for the objective function
        # if num_pt==5:
        # ModelInfo_obj.X_hat = np.append(ModelInfo_obj.X_hat,x0I_hat[nonNAN]).reshape(num_pt,num_xI)
        ModelInfo_obj.X_org = np.append(ModelInfo_obj.X_org,x0I[nonNAN]).reshape(num_pt,num_xI)
        ModelInfo_obj.xC = np.append(ModelInfo_obj.xC,xC_opt[nonNAN]).reshape(num_pt,num_xC)
        ModelInfo_obj.y = np.append(ModelInfo_obj.y,obj[nonNAN]).reshape(num_pt,1)
        ModelInfo_obj.eflag = np.append(ModelInfo_obj.eflag,eflag[nonNAN]).reshape(num_pt,1)
        # Surrogate data for the constraint functions
        if M>0:
            for mm in xrange(M):
                # Put a check here to ensure cons are never NaN or imaginary
                # ModelInfo_g[mm].X_hat = np.append(ModelInfo_g[mm].X_hat,x0I_hat[nonNAN]).reshape(num_pt,num_xI)
                ModelInfo_g[mm].X_org = np.append(ModelInfo_g[mm].X_org,x0I[nonNAN]).reshape(num_pt,num_xI)
                ModelInfo_g[mm].y = np.append(ModelInfo_g[mm].y,g[nonNAN,mm]).reshape(num_pt,1)


        if eflag[nonNAN] >= 1:
            fea_pt += 1
            FEA_obj = np.append(FEA_obj,obj[nonNAN]).reshape(fea_pt,1)
            FEA_xopt = np.append(FEA_xopt,np.concatenate((x0I[nonNAN,:],xC_opt[nonNAN,:]))).reshape(fea_pt,num_xI+num_xC)
    # Call the surrogate building function
    surrogate = KrigingSurrogate() #Use ModelInfo_obj in the future release
    surrogate.train(ModelInfo_obj.X_org, ModelInfo_obj.y)
    ModelInfo_obj.X = surrogate.X
    # ModelInfo_obj.ynorm = surrogate.Y
    ModelInfo_obj.thetas = surrogate.thetas
    ModelInfo_obj.mu = np.mean(surrogate.Y) #This value should always be 0.0
    ModelInfo_obj.SigmaSqr = surrogate.sigma2/np.square(surrogate.Y_std) #This value should always be 1.0
    ModelInfo_obj.c_r = surrogate.alpha
    ModelInfo_obj.R_inv = surrogate.Vh.T.dot(np.einsum('i,ij->ij', surrogate.S_inv, surrogate.U.T))
    ModelInfo_obj.Y_mean = surrogate.Y_mean
    ModelInfo_obj.Y_std = surrogate.Y_std
    ModelInfo_obj.X_std = surrogate.X_std.reshape(num_xI,1)
    ModelInfo_obj.X_mean = surrogate.X_mean.reshape(num_xI,1)
    print "Surrogate building of the objective is complete..."

    # Call the surrogate for the constraints
    if M>0:
        for mm in xrange(M):
            surrogate = KrigingSurrogate()
            surrogate.train(ModelInfo_g[mm].X_org, ModelInfo_g[mm].y)
            ModelInfo_g[mm].X = surrogate.X
            ModelInfo_g[mm].thetas = surrogate.thetas
            ModelInfo_g[mm].mu = np.mean(surrogate.Y) #This value should always be 0.0
            ModelInfo_g[mm].SigmaSqr = surrogate.sigma2/np.square(surrogate.Y_std) #This value should always be 1.0
            ModelInfo_g[mm].c_r = surrogate.alpha
            ModelInfo_g[mm].R_inv = surrogate.Vh.T.dot(np.einsum('i,ij->ij', surrogate.S_inv, surrogate.U.T))
            ModelInfo_g[mm].Y_mean = surrogate.Y_mean
            ModelInfo_g[mm].Y_std = surrogate.Y_std
            ModelInfo_g[mm].X_std = surrogate.X_std.reshape(num_xI,1)
            ModelInfo_g[mm].X_mean = surrogate.X_mean.reshape(num_xI,1)
        print "Surrogate building of the constraints are complete..."

    if len(FEA_obj) >= 1:
        y_opt = np.min(FEA_obj)
        min_ind = FEA_obj.argmin()
        ModelInfo_obj.y_best = y_opt
        x_opt = FEA_xopt[min_ind]

    # Save all the date here
    # print "foobar-Print data for Matlab"
    # print xI_lb, xI_ub
    # print M
    # print surrogate.Y
    # print ModelInfo_g[0].X
    # print ModelInfo_g[0].thetas
    # print ModelInfo_g[0].mu
    # print ModelInfo_g[0].SigmaSqr
    # print ModelInfo_g[0].c_r
    # print ModelInfo_g[0].R_inv
    # print ModelInfo_g[0].Y_mean
    # print ModelInfo_g[0].Y_std
    # print ModelInfo_g[0].X_std
    # print ModelInfo_g[0].X_mean
    # exit()
    ############################################################################
    # Step 4: Maximize the expected improvement function to obtain an integer infill points
    # Choose the solver: 1. MINLP BB, 2. GA, 3. Both MINLP BB & GA
    app_step4 = 1
    print "EGOLF-Iter: %d" % iter
    int_con = range(len(xI_lb))
    if Tot_newpt_added != Tot_pt_prev:
        if app_step4 == 1: # Uses MINLP BB
            print "======================MINLPBB-Start====================================="
            x_new, ei_min, eflag_MINLPBB = MINLP_BB(xI_lb, xI_ub, ModelInfo_obj,ModelInfo_g)
            print "======================MINLPBB-End======================================="
        # elif app_step4 == 2: # Uses GA
        #
        # elif app_step4 == 3: # Using Both MINLP BB and GA (for comparison purpose)
        print "New xI = ", x_new
        print "EI_min = ", ei_min
        print "Eflag = ", eflag_MINLPBB
        # exit()
        if eflag_MINLPBB >= 1:
            new_sol = np.append(new_sol,ei_min)
            x0I = cp.deepcopy(x_new)
            x0I_hat = (x0I - xI_lb)/(xI_ub - xI_lb)
            ei_max = -ei_min
            Tot_pt_prev = Tot_newpt_added

            # Prevent the correlation matrix being close singular
            # No point allowed within the pescribed hypersphere of any existing point
            rad = 0.5
            cc = 0
            for ii in xrange(np.shape(ModelInfo_obj.X_org)[0]):
                dist = np.sum((ModelInfo_obj.X_org[ii] - x0I)**2)**0.5
                if dist <= rad:
                    print "Point already exists!"
                    ec2 = 1
                    break
        else:
            ec2 = 1
    else:
        ec2 = 1

    if np.abs(y_opt)<= 1e-6:
        term = ei_tol_abs
    else:
        term =  np.min(np.array([np.abs(ei_tol_per*y_opt),ei_tol_abs]))

    ############################################################################
    # Step 5: Check for termination
    if ei_max <= term or ec2 == 1 or Tot_newpt_added >= max_pt_lim:
        ter_crit = 1
        if ei_max <= term:
            print "No Further improvement expected! Terminating algorithm.."
        elif ec2 == 1:
            print "No new point found that improves the surrogate. Terminating algorithm.."
        elif Tot_newpt_added >= max_pt_lim:
            print "Maximum allowed sampling limit reached! Terminating algorithm.."

    iter += 1
################################################################################
print "\n===================Result Summary===================="
print "The best objective: %0.4f" % y_opt
print "Total number of continuous minimization: %d" % len(ModelInfo_obj.y)
print "Total number of objective function evaluation: %d" % Tot_FunCount
print "Best Integer designs: ", x_opt[:num_xI]
print "Corresponding continuous designs: ", x_opt[num_xI:]
print "====================================================="
