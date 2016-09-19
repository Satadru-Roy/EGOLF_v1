import numpy as np
import copy as cp
from scipy.special import erf
from Initialization import ModelInfo
from scipy.optimize import minimize
import time
from random import uniform

def MINLP_BB(xI_lb, xI_ub, ModelInfo_obj, ModelInfo_g):
    ''' This is the branch and bound algorithm that maximizes the
    constrained expected improvement function and returns an intger infill point.
    The algorithm uses the relaxation techniques proposed by
    Jones et.al. on their paper on EGO,1998. This enables the algorithm
    to use any gradient-based approach to obtain a global solution. Also, to satisfy
    the integer constraints, a new branching scheme has been implemented.

    Inputs:
    xI_lb - Lower bound integer/discrete type design variables (original design space)
    xI_ub - Lower bound integer/discrete type design variables (original design space)
    ModelInfo_obj - Objective surrogate model data
    ModelInfo_g - Constraint surrogate data

    Outputs:
    xopt - Optimal integer solution
    fopt - Optimal objective value
    eflag_MINLPBB - Exit flag MINLP BB

    Developed by Satadru Roy
    School of Aeronautics & Astronautics
    Purdue University, West Lafayette, IN 47906
    Copyright July, 2016 '''

    # Step 1: Initialize
    iter_=1
    tol_ec=1e-6
    term=0
    LBD=-np.inf
    node_num=0
    fopt=np.inf
    can_pt=0
    canX=[];canF=[];UBD_iter=[]
    Aset=[]; #Note: Active set fields: Aset = [[NodeNumber, lb, ub, LBD, UBD],[],..] #Each node is a list
    xL_iter=cp.deepcopy(xI_lb);xU_iter=cp.deepcopy(xI_ub)
    par_node=0
    LBD_prev=-np.inf
    fC = 0.0
    heur_search=0 # Set this to 1 for good starting point using heuristic search;
    # otherwise set this to 0

    print "Starting the MINLP BB...."
    num_des = len(xI_lb)

    if heur_search == 1:
        # Heuristic search goes here ...
        print "Heuristic search in on....."
    else:
        #Randomly generate an integer point
        xopt = np.round(xI_lb + uniform(0,1)*(xI_ub - xI_lb)).reshape((1,num_des))
        # xopt = np.array([[2.0,2.0,2.0]])
        fopt = combined_obj(xopt,ModelInfo_obj,ModelInfo_g,[],0)
        fC+=1.0
        UBD = 1.0*fopt
        eflag_MINLPBB = 1.0
    print "====================================================================================="
    print "%19s%12s%14s%21s" %("Global","Parent","Child1","Child2")
    print "%s%8s%10s%8s%9s%11s%10s%11s%11s" %("Iter","LBD","UBD","Node", \
    "Node1","LBD1","Node2","LBD2","Flocal")
    print "====================================================================================="

    # print "check before  while loop"
    # print xopt, fopt, UBD, fC
    # exit()
    while term == 0:
        con_fac=[] #concave_factor(xL_iter,xU_iter)
        con_flag = 0
        #Step 2: Obtain a local solution
        app_step2 = 1 #1-Uses gradient based approach
        loc_search = 0 #0-No local search, 1-local search
        if app_step2 == 1:
            xloc_iter = np.round(xL_iter + 0.49*(xU_iter - xL_iter)) #Put this to 0.49 for rounding direction towards left-bottom
            floc_iter = combined_obj(xloc_iter,ModelInfo_obj,ModelInfo_g,con_fac,con_flag)
            fC+=1.0
            efloc_iter = 1.0
            # print "Check with local sampling"
            # print xloc_iter
            # print floc_iter
            # exit()
            if loc_search == 1: #Perform the local search
                if np.abs(floc_iter) > 1.0e-6:
                    xC_iter = xloc_iter
                    bnds = [(xL_iter[ii], xU_iter[ii]) for ii in xrange(num_des)]
                    optResult = minimize(combined_obj,xC_iter,\
                    args=(ModelInfo_obj,ModelInfo_g,con_fac,con_flag), method='SLSQP',\
                    bounds=bnds, options={'ftol':1e-12})
                    xloc_iter = np.round(optResult.x.reshape(num_des,1))
                    floc_iter = combined_obj(xloc_iter,ModelInfo_obj,ModelInfo_g,con_fac,con_flag)
                    fC+= 1.0+optResult.nfev
                    if not optResult.success:
                        efloc_iter=0
                        floc_iter = np.inf
                    else:
                        efloc_iter=1.0
        #elif app_step2 == 2: Other methods goes here

        # Step 3: Partition the current rectangle as per the new branching scheme
        child_info = np.zeros([2,3])
        dis_flag = [' ',' ']
        l_iter = (xU_iter - xL_iter).argmax()
        if xloc_iter[l_iter]<xU_iter[l_iter]:
            delta = 0.5 #0<delta<1
        else:
            delta = -0.5 #-1<delta<0
        for ii in xrange(2):
            lb = cp.deepcopy(xL_iter)
            ub = cp.deepcopy(xU_iter)
            if ii == 0:
                ub[l_iter] = np.floor(xloc_iter[l_iter]+delta)
            elif ii == 1:
                lb[l_iter] = np.ceil(xloc_iter[l_iter]+delta)

            # print lb, ub
            if np.linalg.norm(ub-lb) > 1e-6: #Not a point
                # Step 4: Obtain an LBD of f in the newly created node
                S4_fail = 0
                x_comL, x_comU, Ain_hat, bin_hat = gen_coeff_bound(lb,ub,ModelInfo_obj)
                # print "foobar-step4 LBD"
                # print lb, ub
                # print x_comL, x_comU, Ain_hat, bin_hat
                # exit()
                sU, eflag_sU = maximize_S(x_comL, x_comU, Ain_hat, bin_hat, ModelInfo_obj)
                # print "foobar-Smax"
                # print sU, eflag_sU
                # exit()
                if eflag_sU >= 1:
                    yL, eflag_yL = minimize_y(x_comL, x_comU, Ain_hat, bin_hat, ModelInfo_obj)
                    # print "foobar-ymin"
                    # print yL, eflag_yL
                    # exit()
                    if eflag_yL >= 1:
                        NegEI = calc_conEI_norm([],ModelInfo_obj,sU,yL)
                        # print "foobar-NegEI after yL calc"
                        # print NegEI
                        # if ii==1:
                        #     exit()
                        M = len(ModelInfo_g)
                        EV = np.zeros([M,1])
                        if M>0:
                            # Expected violation goes here
                            for mm in xrange(M):
                                x_comL, x_comU, Ain_hat, bin_hat = gen_coeff_bound(lb,ub,ModelInfo_g[mm])
                                sU_g, eflag_sU_g = maximize_S(x_comL, x_comU, Ain_hat, bin_hat, ModelInfo_g[mm])
                                sL_g = -1.*sU_g
                                if eflag_sU_g >= 1:
                                    yL_g, eflag_yL_g = minimize_y(x_comL, x_comU, Ain_hat, bin_hat, ModelInfo_g[mm])
                                    if eflag_yL_g >= 1:
                                        EV[mm] = calc_conEV_norm([],ModelInfo_g[mm],sL_g,yL_g)
                                    else:
                                        S4_fail = 1.0
                                        break
                                else:
                                    S4_fail = 1.0
                                    break
                    else:
                        # print "Cannot solve Min y_hat problem!"
                        S4_fail = 1.0
                else:
                    # print "Cannot solve Max S problem!"
                    S4_fail = 1.0

                if S4_fail == 1: #Convex approximation failed!
                    if efloc_iter >= 1:
                        LBD_NegConEI = LBD_prev
                    else:
                        LBD_NegConEI = np.inf
                    dis_flag[ii] = 'F'
                else:
                    LBD_NegConEI = (NegEI/(1.0+np.sum(EV)))

                # Step5: Store any new node inside the active set that has LBD lower than the UBD
                if (LBD_NegConEI) < UBD-1.0e-6:
                    node_num += 1
                    newNode = [node_num, lb, ub, LBD_NegConEI, floc_iter]
                    Aset.extend([newNode])
                    child_info[ii] = np.array([node_num, LBD_NegConEI, floc_iter])
                else:
                    child_info[ii] = np.array([par_node, LBD_NegConEI, floc_iter])
                    dis_flag[ii] = 'X' #Flag for child created but not added to active set (fathomed)
            else:
                if ii == 1:
                    xloc_iter = ub
                    floc_iter = combined_obj(xloc_iter,ModelInfo_obj,ModelInfo_g,con_fac,con_flag)
                    fC+=1.0
                child_info[ii] = np.array([par_node, np.inf, floc_iter])
                dis_flag[ii] = 'x' #Flag for No child created

            #Update the active set
            can_pt+=1
            canX = np.reshape(np.append(canX,xloc_iter),(can_pt,num_des))
            canF = np.reshape(np.append(canF,floc_iter),(can_pt,1))
            if floc_iter < UBD: # Better integer solution found
                UBD = 1.0*floc_iter
                fopt = 1.0*UBD
                xopt = cp.deepcopy(xloc_iter).reshape(1,num_des)
                if len(Aset) >= 1:
                    del_flag = []
                    for aaa in xrange(len(Aset)):
                        if Aset[aaa][3] >= UBD:
                            del_flag.extend([aaa])
                    Aset = [ii for jj, ii in enumerate(Aset) if jj not in del_flag]
                    # print del_flag
                UBD_iter = np.append(UBD_iter,UBD)

        # print "foobar-loop end check"
        # print xloc_iter, floc_iter
        # print Aset, S4_fail
        # exit()
        if iter_ % 25 == 0:
            # Display output in a tabular format
            print "====================================================================================="
            print "%19s%12s%14s%21s" %("Global","Parent","Child1","Child2")
            print "%s%8s%10s%8s%9s%11s%10s%11s%11s" %("Iter","LBD","UBD","Node", \
            "Node1","LBD1","Node2","LBD2","Flocal")
            print "====================================================================================="
        print "%3d%10.2f%10.2f%6d%8d%1s%13.2f%8d%1s%13.2f%9.2f" %(iter_,LBD,UBD,\
        par_node,child_info[0,0],dis_flag[0],child_info[0,1],\
        child_info[1,0],dis_flag[1],child_info[1,1],child_info[1,2])

        if len(Aset) >= 1:
            # Update LBD and select the current rectangle
            # a. Set LBD as lowest in the active set
            LBD = min([Aset[ii][3] for ii in xrange(len(Aset))])
            ind_LBD = [Aset[ii][3] for ii in xrange(len(Aset))].index(LBD)
            LBD_prev = LBD
            # b. Select the lowest LBD node as the current node
            xL_iter = Aset[ind_LBD][1]
            xU_iter = Aset[ind_LBD][2]
            par_node = Aset[ind_LBD][0]
            iter_+=1
            # c. Delete the selected node from the Active set of nodes
            del Aset[ind_LBD]
            #Step 7: Check for convergence
            diff = np.abs(UBD - LBD)
            if diff < tol_ec:
                term = 1
                print "====================================================================================="
                print "Terminating! Absolute difference between the upper and lower bound is below the tolerence limit."
        else:
            term = 1
            print "====================================================================================="
            print "Terminating! No new node to explore."
    print "EI fun count: ", fC
    return xopt, fopt, eflag_MINLPBB
#End MINLP_BB
################################################################################
# Supporting modules
################################################################################
def combined_obj(xI,*param):
    ModelInfo_obj=param[0];ModelInfo_g=param[1];con_fac=param[2];flag=param[3]
    X = ModelInfo_obj.X
    k = np.shape(X)[1]
    lb_org = ModelInfo_obj.lb_org
    ub_org = ModelInfo_obj.ub_org
    lb = ModelInfo_obj.lb
    ub = ModelInfo_obj.ub

    if xI.shape[1] > 1:
        xI = xI.T
    xval = (xI - lb_org)/(ub_org - lb_org)     # Normalize to a unit hypercube

    # xval = (xI - ModelInfo_obj.X_mean.T)/ModelInfo_obj.X_std.T # Normalized as per the convention in kriging of openmdao
    NegEI = calc_conEI_norm(xval,ModelInfo_obj)

    M=len(ModelInfo_g)
    EV = np.zeros([M,1])
    if M>0:
        for mm in xrange(M):
            EV[mm] = calc_conEV_norm(xval,ModelInfo_g[mm])

    conNegEI = NegEI/(1.0+np.sum(EV))
    P=0.0
    if flag == 1: #Locally makes ei concave to get rid of flat objective space
        for ii in xrange(k):
            P += con_fac[ii]*(lb[ii] - xval[ii])*(ub[ii] - xval[ii])

    f = conNegEI + P

    return f
################################################################################
def calc_conEI_norm(xval,ModelInfo_obj,*param):
    '''This modules evaluates the expected improvement in the normalized
    design sapce'''
    # y_min = (ModelInfo_obj.y_best - ModelInfo_obj.Y_mean)/ModelInfo_obj.Y_std

    y_min = ModelInfo_obj.y_best
    X = ModelInfo_obj.X
    if len(param) == 0:
        c_r = ModelInfo_obj.c_r
        thetas = ModelInfo_obj.thetas
        SigmaSqr = ModelInfo_obj.SigmaSqr
        R_inv = ModelInfo_obj.R_inv
        mu = ModelInfo_obj.mu
        p = ModelInfo_obj.p
        n = np.shape(X)[0]
        one = np.ones([n,1])
        r = np.ones([n,1])
        for ii in xrange(n):
            r[ii] = np.exp(-np.sum(thetas.T*(xval.T - X[ii])**p))

        y_hat = mu + np.dot(r.T,c_r)
        SSqr = SigmaSqr*(1.0 - r.T.dot(np.dot(R_inv,r)) + \
        ((1.0 - one.T.dot(np.dot(R_inv,r)))**2)/(one.T.dot(np.dot(R_inv,one))))
    else:
        SSqr = param[0]
        y_hat = param[1]

    if SSqr<=0.0:
        NegEI = 0.0
    else:
        ei1 = (y_min-y_hat)*(0.5+0.5*erf((1/np.sqrt(2))*((y_min-y_hat)/np.sqrt(abs(SSqr)))))
        ei2 = np.sqrt(np.abs(SSqr))*(1.0/np.sqrt(2.0*np.pi))*np.exp(-0.5*((y_min-y_hat)**2/np.abs(SSqr)))
        NegEI = -(ei1 + ei2)

    return NegEI
################################################################################
def calc_conEV_norm(xval,ModelInfo_g,*param):
    '''This modules evaluates the expected improvement in the normalized
    design sapce'''
    g_min = 0.0
    X = ModelInfo_g.X
    if len(param) == 0:
        c_r = ModelInfo_g.c_r
        thetas = ModelInfo_g.thetas
        SigmaSqr = ModelInfo_g.SigmaSqr
        R_inv = ModelInfo_g.R_inv
        mu = ModelInfo_g.mu
        p = ModelInfo_g.p
        n = np.shape(X)[0]
        one = np.ones([n,1])
        r = np.ones([n,1])

        for ii in xrange(n):
            r[ii] = np.exp(-np.sum(thetas.T*(xval.T - X[ii])**p))

        g_hat = mu + np.dot(r.T,c_r)
        gSSqr = SigmaSqr*(1.0 - r.T.dot(np.dot(R_inv,r)) + \
        ((1.0 - one.T.dot(np.dot(R_inv,r)))**2)/(one.T.dot(np.dot(R_inv,one))))
    else:
        gSSqr = param[0]
        g_hat = param[1]

    if gSSqr<=0:
        EV = 0.0
    else:
        # Calculate expected violation
        ei1 = (g_hat-g_min)*(0.5+0.5*erf((1.0/np.sqrt(2.0))*((g_hat-g_min)/np.sqrt(np.abs(gSSqr)))))
        ei2 = np.sqrt(np.abs(gSSqr))*(1.0/np.sqrt(2.0*np.pi))*np.exp(-0.5*((g_hat-g_min)**2/np.abs(gSSqr)))
        EV = (ei1 + ei2)

    return EV
################################################################################
def gen_coeff_bound(xI_lb, xI_ub, ModelInfo):
    '''This modules generates the upper and lower bound of the artificial variable
    r and the coefficients for the linearized under estimator constraints.
    The version accepts design bound in the original design space, converts it
    to normalized design space'''

    #Normalized to 0-1 hypercube
    xL_hat = (xI_lb - ModelInfo.lb_org)/(ModelInfo.ub_org - ModelInfo.lb_org)
    xU_hat = (xI_ub - ModelInfo.lb_org)/(ModelInfo.ub_org - ModelInfo.lb_org)

    #Normalized as per openmdao krigging model
    # xL_hat = (xI_lb - ModelInfo.X_mean)/ModelInfo.X_std
    # xU_hat = (xI_ub - ModelInfo.X_mean)/ModelInfo.X_std

    rL, rU = interval_analysis(xL_hat, xU_hat, ModelInfo)

    # Combined design variables for supbproblem
    num = len(xL_hat) + len(rL)
    x_comL = np.append(xL_hat,rL).reshape(num,1)
    x_comU = np.append(xU_hat,rU).reshape(num,1)

    # Coefficients of the linearized constraints of the subproblem
    Ain_hat, bin_hat = lin_underestimator(x_comL, x_comU, ModelInfo)

    return x_comL, x_comU, Ain_hat, bin_hat #Checked
################################################################################
def interval_analysis(lb_x, ub_x, ModelInfo):
    ''' The module predicts the lower and upper bound of the artificial variable
    'r' from the bounds of the design variable x
    r is related to x by the following equation
    r_i = exp(-sum(theta_h*(x_h - x_h_i)^2))'''

    X = ModelInfo.X
    thetas = ModelInfo.thetas
    p = ModelInfo.p
    n = np.shape(X)[0]
    k = np.shape(X)[1]

    t1L = np.zeros([n,k]);t1U = np.zeros([n,k])
    t2L = np.zeros([n,k]);t2U = np.zeros([n,k])
    t3L = np.zeros([n,k]);t3U = np.zeros([n,k])
    t4L = np.zeros([n,1]);t4U = np.zeros([n,1])
    lb_r = np.zeros([n,1]);ub_r = np.zeros([n,1])
    eterm = 1
    if p % 2 == 0:
        for i in xrange(n):
            for h in xrange(k):
                t1L[i,h] = lb_x[h] - X[i,h]
                t1U[i,h] = ub_x[h] - X[i,h]

                t2L[i,h] = np.max(np.array([0,np.min(np.array([t1L[i,h]*t1L[i,h],t1L[i,h]*t1U[i,h],t1U[i,h]*t1U[i,h]]))]))
                t2U[i,h] = np.max(np.array([0,np.max(np.array([t1L[i,h]*t1L[i,h],t1L[i,h]*t1U[i,h],t1U[i,h]*t1U[i,h]]))]))

                t3L[i,h] = np.min(np.array([-thetas[h]*t2L[i,h],-thetas[h]*t2U[i,h]]))
                t3U[i,h] = np.max(np.array([-thetas[h]*t2L[i,h],-thetas[h]*t2U[i,h]]))

            t4L[i] = np.sum(t3L[i,:])
            t4U[i] = np.sum(t3U[i,:])

            lb_r[i] = np.exp(t4L[i])
            ub_r[i] = np.exp(t4U[i])
    else:
        print "\nWarning! Value of p should be 2. Cannot perform interval analysis"
        print "\nReturing global bound of the r variable"
        eterm =0
    return lb_r, ub_r #Checked
################################################################################
def lin_underestimator(lb, ub, ModelInfo):
    X = ModelInfo.X
    thetas = ModelInfo.thetas
    p = ModelInfo.p
    n = np.shape(X)[0]
    k = np.shape(X)[1]

    lb_x = lb[:k]; ub_x = ub[:k]
    lb_r = lb[k:]; ub_r = ub[k:]

    a1 = np.zeros([n,n]);a3 = np.zeros([n,n])
    a1_hat = np.zeros([n,n]);a3_hat = np.zeros([n,n])
    a2 = np.zeros([n,k]);a4 = np.zeros([n,k])
    b2 = np.zeros([n,k]);b4 = np.zeros([n,k])
    b1 = np.zeros([n,1]);b3 = np.zeros([n,1])
    b1_hat = np.zeros([n,1]);b3_hat = np.zeros([n,1])

    for i in xrange(n):
        #T1: Linearize under-estimator of ln[r_i] = a1[i,i]*r[i] + b1[i]
        if ub_r[i] < 1.0e-323 or (ub_r[i] - lb_r[i]) < 1.0e-308:
            # a1[i,i] = 0.
            # b1[i] = -np.inf
            a1_hat[i,i] = 0.0 #a1[i,i]*(ub_r[i]-lb_r[i])
            b1_hat[i] = -np.inf #a1[i,i]*lb_r[i] + b1[i]
        elif ub_r[i] <= lb_r[i]:
            # a1[i,i] = 0.0
            # b1[i] = np.log(ub_r[i])
            a1_hat[i,i] = 0.0 #a1[i,i]*(ub_r[i]-lb_r[i])
            b1_hat[i] = np.log(ub_r[i]) #a1[i,i]*lb_r[i] + b1[i]
        elif lb_r[i] < 1.0e-323:
            # a1[i,i] = np.inf
            # b1[i] = -np.inf
            a1_hat[i,i] = np.inf #a1[i,i]*(ub_r[i]-lb_r[i])
            b1_hat[i] = -np.inf #b1[i]
        else:
            a1[i,i] = ((np.log(ub_r[i]) - np.log(lb_r[i]))/(ub_r[i] - lb_r[i]))
            b1[i] = np.log(ub_r[i]) - a1[i,i]*ub_r[i]
            a1_hat[i,i] = a1[i,i]*(ub_r[i]-lb_r[i])
            b1_hat[i] = a1[i,i]*lb_r[i] + b1[i]

        #T3: Linearize under-estimator of -ln[r_i] = a3[i,i]*r[i] + b3[i]
        if ub_r[i] < 1.0e-323:
            a3_hat[i,i] = 0.0
            b3_hat[i] = np.inf
        else:
            r_m_i = (lb_r[i] + ub_r[i])/2.0
            if r_m_i < 1.0e-308:
                a3_hat[i,i] = -np.inf
                b3_hat[i] = np.inf
            else:
                a3[i,i] = -1.0/r_m_i
                b3[i] = -np.log(r_m_i) - a3[i,i]*r_m_i
                a3_hat[i,i] = a3[i,i]*(ub_r[i] - lb_r[i])
                b3_hat[i] = a3[i,i]*lb_r[i] + b3[i]

        for h in xrange(k):
            #T2: Linearize under-estimator of thetas_h*(x_h - X_h_i)^2 = a4[i,h]*x_h[h] + b4[i,h]
            x_m_h = (ub_x[h] + lb_x[h])/2.0
            a2[i,h] = p*thetas[h]*(x_m_h - X[i,h])**(p-1.0)
            yy = thetas[h]*(x_m_h - X[i,h])**p
            b2[i,h] = -a2[i,h]*x_m_h + yy

            #T4: Linearize under-estimator of -theta_h*(x_h - X_h_i)^2 = a4[i,h]*x_h[h] + b4[i,h]
            yy2 = -thetas[h]*(ub_x[h] - X[i,h])**p
            yy1 = -thetas[h]*(lb_x[h] - X[i,h])**p

            if ub_x[h] <= lb_x[h]:
                a4[i,h] = 0.0
            else:
                a4[i,h] = (yy2 - yy1)/(ub_x[h] - lb_x[h])

            b4[i,h] = -a4[i,h]*lb_x[h] + yy1

    Ain1 = np.concatenate((a2,a4), axis=0)
    Ain2 = np.concatenate((a1_hat,a3_hat), axis=0)
    Ain_hat = np.concatenate((Ain1,Ain2), axis=1)
    bin_hat = np.concatenate((-(b1_hat + np.sum(b2,axis=1).reshape(n,1)),-(b3_hat + np.sum(b4,axis=1).reshape(n,1))),axis=0)

    return Ain_hat, bin_hat #Checked
################################################################################
def maximize_S(x_comL,x_comU,Ain_hat,bin_hat,ModelInfo):
    '''This modules finds an upper bound to the SigmaSqr Error.
    The module scales up 'r' to provide a smooth design space for gradient-based
    approach.'''
    X = ModelInfo.X
    n = np.shape(X)[0]
    k = np.shape(X)[1]
    R_inv = ModelInfo.R_inv
    SigmaSqr = ModelInfo.SigmaSqr
    one = np.ones([n,1])

    xhat_comL = cp.deepcopy(x_comL)
    xhat_comU = cp.deepcopy(x_comU)
    xhat_comL[k:] = np.zeros([n,1])
    xhat_comU[k:] = np.ones([n,1])

    # Calculate the convexity factor alpha
    rL = x_comL[k:]
    rU = x_comU[k:]

    dr_drhat = np.zeros([n,n])
    for ii in xrange(n):
        dr_drhat[ii,ii] = rU[ii,0] - rL[ii,0]

    T2_num = np.dot(np.dot(R_inv,one),np.dot(R_inv,one).T)
    T2_den = np.dot(one.T,np.dot(R_inv,one))
    d2S_dr2 = 2.0*SigmaSqr*(R_inv - (T2_num/T2_den))
    H_hat = np.dot(np.dot(dr_drhat,d2S_dr2),dr_drhat.T)

    # Use Gershgorin's circle theorem to find a lower bound of the
    # min eigen value of the hessian
    eig_lb = np.zeros([n,1])
    for ii in xrange(n):
        dia_ele = H_hat[ii,ii]
        sum_rw = 0.0;sum_col=0.0
        for jj in xrange(n):
            if ii != jj:
                sum_rw += np.abs(H_hat[ii,jj])
                sum_col += np.abs(H_hat[jj,ii])

            eig_lb[ii] = dia_ele - np.min(np.array([sum_rw, sum_col]))
    eig_min = np.min(eig_lb)
    alpha = np.max(np.array([0.0,-0.5*eig_min]))
    ModelInfo.alpha = alpha

    # Maximize S
    x0 = 0.5*(xhat_comL + xhat_comU)
    bnds = [(xhat_comL[ii], xhat_comU[ii]) for ii in xrange(len(xhat_comL))]
    #Note: Python defines constraints like g(x) >= 0
    cons = [{'type' : 'ineq','fun' : lambda x : -np.dot(Ain_hat[ii,:],x) + bin_hat[ii,0],'jac': lambda x:-Ain_hat[ii,:]} for ii in xrange(2*n)]
    optResult = minimize(calc_SSqr_convex,x0,\
    args=(ModelInfo,x_comL,x_comU,xhat_comL,xhat_comU),method='SLSQP',\
    constraints=cons,bounds=bnds,options={'ftol':1e-4,'maxiter':100})
    Neg_sU = optResult.fun
    if not optResult.success:
        eflag_sU=0.0
    else:
        eflag_sU=1.0
        for ii in xrange(2*n):
            if np.dot(Ain_hat[ii,:],optResult.x) >  (bin_hat[ii,0] + 1.0e-6):
                eflag_sU=0.0
                break
    sU = - Neg_sU
    # print "foobar-Max S check"
    # print alpha, x0
    # print sU, eflag_sU
    # exit()
    return sU, eflag_sU
################################################################################
def calc_SSqr_convex(x_com,*param):
    ModelInfo=param[0];x_comL=param[1];x_comU=param[2]
    xhat_comL=param[3];xhat_comU=param[4]

    X = ModelInfo.X
    n = np.shape(X)[0]
    k = np.shape(X)[1]
    R_inv = ModelInfo.R_inv
    SigmaSqr = ModelInfo.SigmaSqr
    alpha = ModelInfo.alpha
    one = np.ones([n,1])
    rL = x_comL[k:]
    rU = x_comU[k:]
    rhat = np.array([x_com[k:]]).reshape(n,1)
    r = rL + rhat*(rU - rL)
    rhat_L = xhat_comL[k:]
    rhat_U = xhat_comU[k:]
    term1 = -SigmaSqr*(1.0-r.T.dot(np.dot(R_inv,r)) + \
    ((1.0-one.T.dot(np.dot(R_inv,r)))**2/(one.T.dot(np.dot(R_inv,one)))))

    term2 = alpha*np.sum((rhat-rhat_L)*(rhat-rhat_U))
    S2 = term1 + term2
    return S2[0,0]
################################################################################
def minimize_y(x_comL, x_comU, Ain_hat, bin_hat, ModelInfo):
    app = 1 # 1- Formulates y_hat as LP (weaker bound) 2- Uses non-convex relaxation technique (stronger bound) [Future release]
    X = ModelInfo.X
    n = np.shape(X)[0]
    k = np.shape(X)[1]

    xhat_comL = cp.deepcopy(x_comL)
    xhat_comU = cp.deepcopy(x_comU)
    xhat_comL[k:] = np.zeros([n,1])
    xhat_comU[k:] = np.ones([n,1])

    # print "check Min y"
    # print xhat_comL, xhat_comU, Ain_hat, bin_hat
    # exit()
    if app == 1:
        x0 = 0.5*(xhat_comL+xhat_comU)
        bnds = [(xhat_comL[ii], xhat_comU[ii]) for ii in xrange(len(xhat_comL))]
        cons = [{'type' : 'ineq','fun' : lambda x : -np.dot(Ain_hat[ii,:],x) + bin_hat[ii,0],'jac': lambda x:-Ain_hat[ii,:]} for ii in xrange(2*n)]
        optResult = minimize(calc_y_hat_convex,x0,\
        args=(x_comL,x_comU,ModelInfo),method='SLSQP',\
        constraints=cons,bounds=bnds,options={'ftol':1e-4,'maxiter':10000})
        yL = optResult.fun
        if not optResult.success:
            eflag_yL=0.0
        else:
            eflag_yL=1.0
            for ii in xrange(2*n):
                # print np.dot(Ain_hat[ii,:],optResult.x) - bin_hat[ii,0]
                # print np.dot(Ain_hat[ii,:],optResult.x.T) - bin_hat[ii,0]
                if np.dot(Ain_hat[ii,:],optResult.x.T) >  (bin_hat[ii,0] + 1.0e-6):
                    eflag_yL=0.0
                    break
    # print "foobar-Min y check"
    # print x0
    # print yL, eflag_yL
    # print
    # exit()
    return yL, eflag_yL
################################################################################
def calc_y_hat_convex(x_com,*param):
    x_comL=param[0];x_comU=param[1];ModelInfo=param[2]
    X = ModelInfo.X
    n = np.shape(X)[0]
    k = np.shape(X)[1]
    c_r = ModelInfo.c_r
    mu = ModelInfo.mu

    rL = x_comL[k:]
    rU = x_comU[k:]
    rhat = np.array([x_com[k:]]).reshape(n,1) #TODO remove scaling of 'r'
    r = rL + rhat*(rU - rL)

    y_hat = mu + np.dot(r.T,c_r)
    # print "foobar-calc_y check"
    # print x_com, rhat
    # print mu, r, c_r
    # print r
    # print x_com
    # print y_hat
    # exit()
    return y_hat[0,0]
################################################################################
def concave_factor(xL,xU):
    K = len(xL)
    per_htm = 5.0
    con_fac = np.zeros([K,1])
    for k in xrange(K):
        if np.abs(xL[k]-xU[k]) < 1e-10:
            con_fac[k] = 0.0
        else:
            h_req = (per_htm/100.0)*(xU[k] - xL[k])
            xm = (xL[k]+xU[k])/2.0
            h_act = (xm[k]-xL[k])*(xm[k]-xU[k])
            con_fac[k] = h_req/h_act
    return con_fac
################################################################################

if __name__ == "__main__":
    # Test MINLP_BB with dummy inputs
    xI_lb = np.array([[1.0],[1.0],[1.0]])
    xI_ub = np.array([[4.0],[4.0],[4.0]])
    M = 3
    n=5
    num_xI = len(xI_lb)
    #Objective data
    ModelInfo_obj = ModelInfo(xI_lb, xI_ub, num_xI)
    ModelInfo_obj.lb_org = xI_lb.copy()
    ModelInfo_obj.ub_org = xI_ub.copy()
    ModelInfo_obj.X = np.array([[0.0,0.66666666666667,0.66666666666667],\
                                [0.66666666666667,1.0,1.0],\
                                [0.0,0.66666666666667,0.0],\
                                [0.33333333333333,0.0,0.33333333333333],\
                                [1.0,0.33333333333333,0.66666666666667]])
    ModelInfo_obj.X_org = np.array([[1.0,3.0,3.0],[3.0,4.0,4.0],[1.0,3.0,1.0],[2.0,1.0,2.0],[4.0,2.0,3.0]])
    ModelInfo_obj.xC = np.array([[10.000000000000000,7.611934707614728,2.884713617365940],\
       [6.139509615250542,9.850382457843574,0.000000000100000],\
       [10.000000000000000,7.639024631926064,8.526550379139559],\
       [9.251987139688117,10.000000000000000,1.527127041337365],\
       [10.227954337559884,3.968105909184939,-0.227954337459885]])
    ModelInfo_obj.y = np.array([[15.626344347459352],[18.603914134812975],\
    [15.705089842796930],[11.471743511194173],[17.112772214006487]])
    ModelInfo_obj.eflag = np.array([[5.0],[1.0],[5.0],[5.0],[-2.0]])
    ModelInfo_obj.y_best = np.min(ModelInfo_obj.y)
    ModelInfo_obj.thetas = np.array([[10.0**1.022038094255710],[10.0**0.569609718087874],[10.0**-2.99999023616931]]) #10.**np.array([[0.465979016183845]])
    ModelInfo_obj.p = 2
    ModelInfo_obj.R_inv = np.array([[1125.44012840546,-0.00907704868570349,-1124.93817728487,-0.0298893649269364,0.000708142016463804],\
    [-0.00907704868570349,1.00365956919566,0.00330061256788934,-0.00690020760577544,-0.0598505705497161],\
    [ -1124.93817728487,0.00330061256788934,1125.44005717367,-0.0299745724787454,-1.21069549498704e-05],\
    [ -0.0298893649269364,-0.00690020760577544,-0.0299745724787454,1.00366035784801,-0.00577790767390208],\
    [0.000708142016463804,-0.0598505705497161,-1.21069549498704e-05,-0.00577790767390208,1.00360720056043]])
    ModelInfo_obj.SigmaSqr = 7.024300865720791
    ModelInfo_obj.mu = 15.716396564175804
    ModelInfo_obj.c_r = ModelInfo_obj.R_inv.dot(ModelInfo_obj.y-(np.ones([n,1])*ModelInfo_obj.mu))

    #constraints
    ModelInfo_g=[[]]*M
    #Constraint 1 data
    ModelInfo_g[0] = ModelInfo(xI_lb, xI_ub, num_xI)
    ModelInfo_g[0].lb_org = xI_lb.copy()
    ModelInfo_g[0].ub_org = xI_ub.copy()
    ModelInfo_g[0].X = np.array([[0.0,0.66666666666667,0.66666666666667],\
                                [0.66666666666667,1.0,1.0],\
                                [0.0,0.66666666666667,0.0],\
                                [0.33333333333333,0.0,0.33333333333333],\
                                [1.0,0.33333333333333,0.66666666666667]])
    ModelInfo_g[0].X_org = np.array([[1.0,3.0,3.0],[3.0,4.0,4.0],[1.0,3.0,1.0],[2.0,1.0,2.0],[4.0,2.0,3.0]])
    # ModelInfo_g[0].xC = np.array([[10.000000000000000   7.611934707614728   2.884713617365940],\
    #    [6.139509615250542   9.850382457843574   0.000000000100000],\
    #    [10.000000000000000   7.639024631926064   8.526550379139559],\
    #    [9.251987139688117  10.000000000000000   1.527127041337365],\
    #    [10.227954337559884   3.968105909184939  -0.227954337459885]])
    ModelInfo_g[0].y = np.array([[5.78500136683147e-10],[-0.590082052917028],\
    [1.66796954026438e-09],[-0.430705218486039],[0.105060456476318]])

    # ModelInfo_g[0].eflag = np.array([[5.0],[1.0],[5.0],[5.0],[-2.0]])
    # ModelInfo_g[0].y_best = np.min(ModelInfo_g.y)
    ModelInfo_g[0].thetas = np.array([[10.0**1.914806493588213],[10.0**1.834555468897605],[10.0**-2.999998679277249]]) #10.**np.array([[0.465979016183845]])
    ModelInfo_g[0].p = 2
    ModelInfo_g[0].R_inv = np.array([[1125.49665286775,-1.03579143745391e-19,-1124.99654175630,-3.51331249730258e-18,9.69285761486067e-37],\
    [-1.03579143745391e-19,1.00000000000000,3.45263789749402e-20,-2.29347906849811e-34,-7.02506386458047e-18],\
    [-1124.99654175630,3.45263789749402e-20,1125.49665286775,-3.51331249730304e-18,1.07824296558663e-40],\
    [-3.51331249730258e-18,-2.29347906849811e-34,-3.51331249730304e-18,1.00000000000000,-6.90681064649094e-20],\
    [9.69285761486067e-37,-7.02506386458047e-18,1.07824296558663e-40,-6.90681064649094e-20,1.00000000000000]])
    ModelInfo_g[0].SigmaSqr = 0.067022852107014
    ModelInfo_g[0].mu = -0.228918985690923
    ModelInfo_g[0].c_r = ModelInfo_g[0].R_inv.dot(ModelInfo_g[0].y-(np.ones([n,1])*ModelInfo_g[0].mu))

    #Constraint 2 data
    ModelInfo_g[1] = ModelInfo(xI_lb, xI_ub, num_xI)
    ModelInfo_g[1].lb_org = xI_lb.copy()
    ModelInfo_g[1].ub_org = xI_ub.copy()
    ModelInfo_g[1].X = np.array([[0.0,0.66666666666667,0.66666666666667],\
                                [0.66666666666667,1.0,1.0],\
                                [0.0,0.66666666666667,0.0],\
                                [0.33333333333333,0.0,0.33333333333333],\
                                [1.0,0.33333333333333,0.66666666666667]])
    ModelInfo_g[1].X_org = np.array([[1.0,3.0,3.0],[3.0,4.0,4.0],[1.0,3.0,1.0],[2.0,1.0,2.0],[4.0,2.0,3.0]])
    # ModelInfo_g[1].xC = np.array([[10.000000000000000   7.611934707614728   2.884713617365940],\
    #    [6.139509615250542   9.850382457843574   0.000000000100000],\
    #    [10.000000000000000   7.639024631926064   8.526550379139559],\
    #    [9.251987139688117  10.000000000000000   1.527127041337365],\
    #    [10.227954337559884   3.968105909184939  -0.227954337459885]])
    ModelInfo_g[1].y = np.array([[-0.624313944607018],[-1.23824395181771e-09],\
    [-0.625643064976861],[9.17266262945304e-12],[0.176556341141486]])

    # ModelInfo_g[1].eflag = np.array([[5.0],[1.0],[5.0],[5.0],[-2.0]])
    # ModelInfo_g[1].y_best = np.min(ModelInfo_g.y)
    ModelInfo_g[1].thetas = np.array([[10.0**1.75267291844906],[10.0**1.94803457741727],[10.0**-2.99999709235029]]) #10.**np.array([[0.465979016183845]])
    ModelInfo_g[1].p = 2
    ModelInfo_g[1].R_inv = np.array([[1125.49254209951,-9.40437651295121e-16,-1124.99243098766,-6.97368882673638e-21,-1.39876957796595e-29],\
    [-9.40437651295121e-16,1.00000000000000,3.13479196459711e-16,1.74869039997435e-35,-1.39442789041313e-20],\
    [-1124.99243098766,3.13479196459711e-16,1125.49254209951,-6.97368882673579e-21,2.03469344569694e-39],\
    [-6.97368882673638e-21,1.74869039997435e-35,-6.97368882673579e-21,1.00000000000000,-6.27097748898778e-16],\
    [-1.39876957796595e-29,-1.39442789041313e-20,2.03469344569694e-39,-6.27097748898778e-16,1.00000000000000]])
    ModelInfo_g[1].SigmaSqr = 0.074709184637465
    ModelInfo_g[1].mu = -0.112134032769722
    ModelInfo_g[1].c_r = ModelInfo_g[1].R_inv.dot(ModelInfo_g[1].y-(np.ones([n,1])*ModelInfo_g[1].mu))

    #Constraint 3 data
    ModelInfo_g[2] = ModelInfo(xI_lb, xI_ub, num_xI)
    ModelInfo_g[2].lb_org = xI_lb.copy()
    ModelInfo_g[2].ub_org = xI_ub.copy()
    ModelInfo_g[2].X = np.array([[0.0,0.66666666666667,0.66666666666667],\
                                [0.66666666666667,1.0,1.0],\
                                [0.0,0.66666666666667,0.0],\
                                [0.33333333333333,0.0,0.33333333333333],\
                                [1.0,0.33333333333333,0.66666666666667]])
    ModelInfo_g[2].X_org = np.array([[1.0,3.0,3.0],[3.0,4.0,4.0],[1.0,3.0,1.0],[2.0,1.0,2.0],[4.0,2.0,3.0]])
    # ModelInfo_g[2].xC = np.array([[10.000000000000000   7.611934707614728   2.884713617365940],\
    #    [6.139509615250542   9.850382457843574   0.000000000100000],\
    #    [10.000000000000000   7.639024631926064   8.526550379139559],\
    #    [9.251987139688117  10.000000000000000   1.527127041337365],\
    #    [10.227954337559884   3.968105909184939  -0.227954337459885]])
    ModelInfo_g[2].y = np.array([[-0.799251531258173],[9.09372777080364e-09],\
    [-0.649335674057793],[-0.905410643421811],[-0.204916110619887]])

    # ModelInfo_g[2].eflag = np.array([[5.0],[1.0],[5.0],[5.0],[-2.0]])
    # ModelInfo_g[2].y_best = np.min(ModelInfo_g.y)
    ModelInfo_g[2].thetas = np.array([[10.0**0.846781351101145],[10.0**0.111475071389736],[10.0**-0.648543251374813]]) #10.**np.array([[0.465979016183845]])
    ModelInfo_g[2].p = 2
    ModelInfo_g[2].R_inv = np.array([[5.54536361397759,-0.0426244147112794,-4.98268790802283,-0.137228297191174,0.0150273667062081],\
    [-0.0426244147112794,1.08047590606653,0.0337130928183813,-0.110770492761965,-0.267625325787606],\
    [-4.98268790802283,0.0337130928183813,5.54474147692432,-0.145084931027714,-0.00311100078459050],\
    [-0.137228297191174,-0.110770492761965,-0.145084931027714,1.08406724836207,-0.0122473123768546],\
    [0.0150273667062081,-0.267625325787606,-0.00311100078459050,-0.0122473123768546,1.06775736060134]])
    ModelInfo_g[2].SigmaSqr = 0.126924622798639
    ModelInfo_g[2].mu = -0.459719857128503
    ModelInfo_g[2].c_r = ModelInfo_g[2].R_inv.dot(ModelInfo_g[2].y-(np.ones([n,1])*ModelInfo_g[2].mu))

    # ModelInfo_obj.X_mean = 0.;ModelInfo_obj.X_std=1
    # ModelInfo_obj.Y_mean = 0.;ModelInfo_obj.Y_std=1
    # # xopt, fopt, eflag = MINLP_BB(xI_lb, xI_ub, ModelInfo_obj, [])
    # ModelInfo_obj.alpha = 2.75835168164
    # x_comL= np.array([[ 0.        ],[ 0.65641594],[ 0.65641593],[ 0.14471593]])
    # x_comU= np.array([[ 0.46666667],[ 1.        ],[ 1.        ],[ 0.57705012]])
    # xhat_comL= np.array([[ 0.],[ 0.],[ 0.],[ 0.]])
    # xhat_comU= np.array([[ 0.46666667],[ 1.        ],[ 1.        ],[ 1.        ]])
    # x_com = np.array([[0.233292711176682],[0.349656479929390],[0.349258214033471],[ 0.474123962414235]])
    # _ = calc_SSqr_convex(x_com,ModelInfo_obj,x_comL,x_comU,xhat_comL,xhat_comU)
    _ = MINLP_BB(xI_lb, xI_ub, ModelInfo_obj, ModelInfo_g)
    # print xopt
    # print fopt
    # print eflag
