%Dummy inputs from AMIEGO-Python version
% function [xI_lb,xI_ub,ModelInfo_obj,ModelInfo_g] = DummyInputs()
clc;clear;
xI_lb = [-5];
xI_ub = [10];
ModelInfo_obj = [];
ModelInfo_obj.p = 2;
ModelInfo_obj.lb_org = [-5];
ModelInfo_obj.ub_org = [10];
ModelInfo_obj.lb = [-5];
ModelInfo_obj.ub = [10];
ModelInfo_obj.y = [0.5941507;0.81433758;-1.40848828];
ModelInfo_obj.y_best = min(ModelInfo_obj.y);
ModelInfo_obj.X = [-1.36652966; 0.36791183; 0.99861783];
ModelInfo_obj.Theta = log10(34.50110524);
ModelInfo_obj.mu = -2.22044604925e-16;
ModelInfo_obj.SigmaSqr = 1.00000084;
ModelInfo_obj.c_r = [0.5941507;0.81433912;-1.40848918];
ModelInfo_obj.R_inv = [[  1.00000000e+00   8.41215593e-46   1.52292631e-84]
 [ -8.41215593e-46   1.00000000e+00  -1.09561297e-06]
 [  9.21646716e-52  -1.09561297e-06   1.00000000e+00]];
ModelInfo_obj.Psi = inv(ModelInfo_obj.R_inv);
ModelInfo_obj.Y_mean = [ 12.89036781];
ModelInfo_obj.Y_std = [ 7.7723239];
ModelInfo_obj.X_std = [ 6.3420992];
ModelInfo_obj.X_mean = [ 3.66666667];

ModelInfo_g=[];
[x_new, ei_min, eflag_MINLP_EI] = MINLP_BB(xI_lb,xI_ub,ModelInfo_obj,ModelInfo_g);