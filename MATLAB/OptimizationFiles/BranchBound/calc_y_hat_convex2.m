% This function evaluates the linear y_hat w.r.t r
function y_hat = calc_y_hat_convex2(x_com,x_comL,x_comU,ModelInfo)
X=ModelInfo.X;
R_inv=ModelInfo.R_inv; %U=ModelInfo.U;
y=ModelInfo.y;
mu = ModelInfo.mu;
n = size(X,1);
k = size(X,2);
one = ones(n,1);
rL = x_comL(k+1:end);
rU = x_comU(k+1:end);
rhat = x_com(k+1:end);
r = rL + rhat.*(rU-rL);

c = (R_inv*(y-one*mu));
y_hat = mu + r'*c;
% r
% x_com
% y_hat
% keyboard