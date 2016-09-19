%% Dummy file to compare Matlab and Python
clc;clear;
load 3bar_Inputs2BB
[xopt,fopt,eflag_MINLPBB] = MINLP_BB(xI_lb,xI_ub,ModelInfo_obj,ModelInfo_g)