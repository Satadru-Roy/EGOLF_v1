% This function finds a lower bound of y_hat
function [yL, eflag_yL] = minimize_y2(x_comL, x_comU, Ain_hat, bin_hat, ModelInfo)
% Version 2: Ensures r is always locally bounded between [0,1] even though it has
% a global bound of [0,1]. This provides a smooth function of SQP. Without
% this scaling bound of r can be as low as [1.2e-13,1.1e-23]

app = 1;
X=ModelInfo.X;
n = size(X,1);
k = size(X,2);

xhat_comL = x_comL;
xhat_comU = x_comU;
xhat_comL(k+1:end) = zeros(n,1);
xhat_comU(k+1:end) = ones(n,1);
% keyboard
if app==1 % This solves the LP that minimizes the y_hat
    alg = 1;
    if alg == 1 % Use a gradient-based approach to solve this
        x0 = 0.5*(xhat_comL+xhat_comU);
        options = optimset('Display','none','algorithm','active-set');
        [x_com,yL,eflag_yL] = fmincon(@(x_com) calc_y_hat_convex2(x_com,x_comL,x_comU,ModelInfo),...
            x0,Ain_hat,bin_hat,[],[],xhat_comL,xhat_comU,[],options);
    elseif alg ==2 % Use LP solver to solve this
        R=ModelInfo.Psi; %U=ModelInfo.U;
        y=ModelInfo.y;
        mu = ModelInfo.mu;
        one = ones(n,1);
        
        rL = x_comL(k+1:end);
        rU = x_comU(k+1:end);

        c = (R\(y-one*mu)).*(rU-rL);
        fcoef = [zeros(k,1);c];
        [~,fopt,eflag_yL] = linprog(fcoef,Ain_hat,bin_hat,[],[],xhat_comL,xhat_comU);
        yL = mu + (c'*rL) + fopt;
    end
elseif app == 2 %This approach provides a tighter convex approximation of y_hat
    %Define this part
end
% keyboard

