% This is the objective function that evaluates the SSqr error associated
% with the surrogate model
function f = calc_SSqr_convex2(x_com, ModelInfo,x_comL,x_comU,xhat_comL,xhat_comU)
X=ModelInfo.X;
R_inv=ModelInfo.R_inv; %U=ModelInfo.U;
SigmaSqr = ModelInfo.SigmaSqr;
alpha = ModelInfo.alpha;
n = size(X,1);
k = size(X,2);
one = ones(n,1); 

rL = x_comL(k+1:end);
rU = x_comU(k+1:end);
rhat = x_com(k+1:end);

r = rL + rhat.*(rU-rL);
rhat_L = xhat_comL(k+1:end);
rhat_U = xhat_comU(k+1:end);

term1 = -SigmaSqr*(1 - (r'*(R_inv*r)) + ((1-one'*(R_inv*r))^2/(one'*(R_inv*one))));
term2 = alpha*sum((rhat-rhat_L).*(rhat-rhat_U));
f = term1 + term2;
% x_com
% f
% % pause(3)
% keyboard