# Initialization Module of EGOLF
import numpy as np
import copy

def initialize_test(num_des, prob):
    ''' Define the bounds of the integer type design variables'''
    #Please kindly put ##.0 while defining bounds(even for integer type design variables) or you dare!
    if prob == 1:
        # Branin Function
        xI_lb = np.array([[-5.0]]); xI_ub = np.array([[10.0]])
        M = 0
    elif prob == 2:
        # Griewank Function
        xI_lb = -5.0*np.ones([num_des/2,1]); xI_ub = 5.0*np.ones([num_des/2,1])
        M = 0
    elif prob == 3:
        # Rosenbrock function
        xI_lb = -5.0*np.ones([num_des/2,1]); xI_ub = 10.0*np.ones([num_des/2,1])
        M = 0
    elif prob == 4:
        # 3bar Truss Problem
        xI_lb = 1.0*np.ones([3,1]); xI_ub = 4.0*np.ones([3,1]) #Discrete material choice
        M = 3
    return xI_lb, xI_ub, M

def initialize_cont_test(num_des, prob):
    ''' Define the bounds of the continuous type design variables'''
    if prob == 1:
        # Branin Function
        xC_lb = np.array([[0.0]]); xC_ub = np.array([[15.0]])
    elif prob == 2:
        # Griewank Function
        xC_lb = -5.0*np.ones([num_des/2,1]); xC_ub = 5.0*np.ones([num_des/2,1])
    elif prob == 3:
        # Rosenbrock function
        xC_lb = -5.0*np.ones([num_des/2,1]); xC_ub = 10.0*np.ones([num_des/2,1])
    elif prob == 4:
        # 3bar Truss Problem
        xC_lb = 1.0e-10*np.ones([3,1]); xC_ub = 10.0*np.ones([3,1]) #in [cm]

    return xC_lb, xC_ub

class ModelInfo:
    def __init__(self,lb_org,ub_org,num):
        self.lb_org = copy.copy(lb_org)
        self.ub_org = copy.copy(ub_org)
        self.lb = np.zeros([num,1])
        self.ub = np.ones([num,1])
        self.X_hat =[]; self.X_org =[]; self.xC =[]; self.y =[]; self.eflag =[]; self.y_best = []
        self.thetas=[];self.R_inv=[];self.c_r=[];self.SigmaSqr=[];self.mu=[];self.alpha=[]
        self.p=2;
        self.X = []; self.ynorm = []

if __name__ == "__main__":
    [xC_lb,xC_ub] = initialize_cont_test(2,1)
    print xC_lb, xC_ub
