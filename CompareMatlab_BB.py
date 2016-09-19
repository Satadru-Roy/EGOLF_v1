""" Standalone vs matlab."""

import numpy as np

from openmdao.api import IndepVarComp, Group, Problem, ExecComp
from openmdao.drivers.branch_and_bound import Branch_and_Bound
from openmdao.test.three_bar_truss import ThreeBarTruss, ThreeBarTrussVector
from openmdao.test.util import assert_rel_error

xI_lb = np.array([[1.0],[1.0],[1.0]])
xI_ub = np.array([[4.0],[4.0],[4.0]])

M = 3
n=5
num_xI = len(xI_lb)

class ModelInfo(object):
    def __init__(self, a, b, c):
        pass

#Objective data
ModelInfo_obj = ModelInfo(xI_lb, xI_ub, num_xI)
ModelInfo_obj.lb_org = xI_lb.copy()
ModelInfo_obj.ub_org = xI_ub.copy()
ModelInfo_obj.X = np.array([[0.0,0.66666666666667,0.66666666666667],\
                            [0.66666666666667,1.0,1.0],\
                            [0.0,0.66666666666667,0.0],\
                            [0.33333333333333,0.0,0.33333333333333],\
                            [1.0,0.33333333333333,0.66666666666667]])
ModelInfo_obj.X_mean = np.mean(ModelInfo_obj.X)
ModelInfo_obj.X_std = np.std(ModelInfo_obj.X)
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

# Constraints
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
ModelInfo_g[0].X_mean = np.mean(ModelInfo_g[0].X)
ModelInfo_g[0].X_std = np.std(ModelInfo_g[0].X)
ModelInfo_g[0].X_org = np.array([[1.0,3.0,3.0],[3.0,4.0,4.0],[1.0,3.0,1.0],[2.0,1.0,2.0],[4.0,2.0,3.0]])
ModelInfo_g[0].y = np.array([[5.78500136683147e-10],[-0.590082052917028],\
[1.66796954026438e-09],[-0.430705218486039],[0.105060456476318]])
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
ModelInfo_g[1].X_mean = np.mean(ModelInfo_g[1].X)
ModelInfo_g[1].X_std = np.std(ModelInfo_g[1].X)
ModelInfo_g[1].X_org = np.array([[1.0,3.0,3.0],[3.0,4.0,4.0],[1.0,3.0,1.0],[2.0,1.0,2.0],[4.0,2.0,3.0]])
ModelInfo_g[1].y = np.array([[-0.624313944607018],[-1.23824395181771e-09],\
[-0.625643064976861],[9.17266262945304e-12],[0.176556341141486]])
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
ModelInfo_g[2].X_mean = np.mean(ModelInfo_g[2].X)
ModelInfo_g[2].X_std = np.std(ModelInfo_g[2].X)
ModelInfo_g[2].X_org = np.array([[1.0,3.0,3.0],[3.0,4.0,4.0],[1.0,3.0,1.0],[2.0,1.0,2.0],[4.0,2.0,3.0]])
ModelInfo_g[2].y = np.array([[-0.799251531258173],[9.09372777080364e-09],\
[-0.649335674057793],[-0.905410643421811],[-0.204916110619887]])
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

# Build the model

prob = Problem()
root = prob.root = Group()

root.add('xc_a1', IndepVarComp('area1', 5.0), promotes=['*'])
root.add('xc_a2', IndepVarComp('area2', 5.0), promotes=['*'])
root.add('xc_a3', IndepVarComp('area3', 5.0), promotes=['*'])
root.add('xi_m1', IndepVarComp('mat1', 1), promotes=['*'])
root.add('xi_m2', IndepVarComp('mat2', 1), promotes=['*'])
root.add('xi_m3', IndepVarComp('mat3', 1), promotes=['*'])
root.add('comp', ThreeBarTruss(), promotes=['*'])

bb = prob.driver = Branch_and_Bound()
#prob.driver.options['disp'] = False
root.deriv_options['type'] = 'fd'

#bb.add_desvar('area1', lower=0.0005, upper=10.0)
#bb.add_desvar('area2', lower=0.0005, upper=10.0)
#bb.add_desvar('area3', lower=0.0005, upper=10.0)
bb.add_desvar('mat1', lower=1, upper=4)
bb.add_desvar('mat2', lower=1, upper=4)
bb.add_desvar('mat3', lower=1, upper=4)
bb.add_objective('mass')
bb.add_constraint('stress', upper=1.0)

# Hack in our Matlab input

bb.obj_surrogate = ModelInfo_obj
bb.con_surrogate = ModelInfo_g
bb.xI_lb = xI_lb
bb.xI_ub = xI_ub
bb.standalone = False

# Run

prob.setup(check=False)

from time import time
t0 = time()
prob.run()
print("Elapsed time:", time()-t0)

print(bb.xopt, bb.fopt)
