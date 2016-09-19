% This function finds an upper bound to the SigmaSqr Error
function [sU, eflag_sU] = maximize_S2(x_comL,x_comU,Ain_hat,bin_hat,ModelInfo)
% Version 2: Ensures r is locally bounded between [0,1] even though it has
% a global bound of [0,1]. This provides a smooth function of SQP. Without
% this scaling bound of r can be as low as [1.2e-13,1.1e-23]

X=ModelInfo.X;
n = size(X,1);
k = size(X,2);
R_inv = ModelInfo.R_inv;
SigmaSqr = ModelInfo.SigmaSqr;
one=ones(n,1);

xhat_comL = x_comL;
xhat_comU = x_comU;
xhat_comL(k+1:end) = zeros(n,1);
xhat_comU(k+1:end) = ones(n,1);

%% Calculate alpha
rL = x_comL(k+1:end);
rU = x_comU(k+1:end);
dr_drhat = zeros(n,n);
for ii = 1:n
    dr_drhat(ii,ii) = (rU(ii)-rL(ii));
end
d2S_dr2 = 2*SigmaSqr*(R_inv - (((R_inv*one)*(R_inv*one)')/(one'*(R_inv*one))));
H_hat = (dr_drhat*d2S_dr2)*dr_drhat';
% eig_val = eig(H_hat);
% eig_min = min(real(eig_val));
%% Use Gershgorin's circle theorem to find a lower bound of the min eigen value
eig_lb = zeros(n,1);
for ii = 1:n
    dia_ele = H_hat(ii,ii);
    sum_rw=0;sum_col=0;
    for jj = 1:n
        if ii~=jj
            sum_rw = sum_rw + abs(H_hat(ii,jj));
            sum_col = sum_col + abs(H_hat(jj,ii));
        end
        eig_lb(ii) = dia_ele - min(sum_rw,sum_col);
    end
end
eig_min = min(real(eig_lb));

alpha = max(0,-0.5*eig_min);
if alpha<0 || isnan(alpha)==1
    fprintf('\n%s\n','Error: Alpha is negative or NAN!!')
end
% if imag(eig_val(min_ind))~= 0
%     fprintf('\n%s','Minimum eigen value has imaginary component! Considering only the real part..')
% end
ModelInfo.alpha = alpha;

% %% test with the Matlab calculated hessian
% r_hat = sym('rh',[n,1]);
% r = rL + r_hat.*(rU-rL);
% f = -SigmaSqr*(1 - (r'*(R\r)) + (1 - (one'*(R\r)))^2/(one'*(R\one)));
% Hmat = hessian(f,r_hat);
% Hmat_val = subs(Hmat,r_hat,ones(n,1));
% eig_mat = eig(Hmat_val)

x0 = 0.5*(xhat_comL + xhat_comU);
options = optimset('MaxFunEvals',100000,'Display','none','algorithm','active-set');
[x_com, Neg_sU,eflag_sU] = fmincon(@(x_com) calc_SSqr_convex2(x_com,...
    ModelInfo,x_comL,x_comU,xhat_comL,xhat_comU),x0, Ain_hat, bin_hat, [],[],xhat_comL,xhat_comU,[],options);
sU = -Neg_sU; %Upper bound on sU
% keyboard
