% This function generates the upper and lower bound of the artificial
% variable r and the coefficients for the linearized under estimator constraints
function [x_comL, x_comU, Ain_hat, bin_hat] = gen_coeff_bound2(xI_lb,xI_ub,ModelInfo)
%Version 2: Accepts the design bound in the original design space, convert
%it to the normalized design space
xL_hat = (xI_lb - ModelInfo.lb_org)./(ModelInfo.ub_org- ModelInfo.lb_org);
xU_hat = (xI_ub - ModelInfo.lb_org)./(ModelInfo.ub_org- ModelInfo.lb_org);

[rL,rU] = interval_analysis(xL_hat,xU_hat,ModelInfo);
%Combined design variables x_com = [x;r]
x_comL = [xL_hat;rL];
x_comU = [xU_hat;rU]; 
% Linearize the constraints and obtain the coefficient matrix
[Ain_hat, bin_hat] = lin_underestimator2(x_comL, x_comU, ModelInfo);
% keyboasrd
% %% Checkpost: 
% fprintf('\n%s','Step 4 (gen Coeff/bound) successfully completed')