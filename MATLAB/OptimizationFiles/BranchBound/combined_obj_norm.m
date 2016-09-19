% This is the objective function that minimizes the negative of the 
% constrained expected improvement function
% This is in the normalized [0,1] design space
function f = combined_obj_norm(xval,ModelInfo_obj,ModelInfo_g,con_fac,flag)

X = ModelInfo_obj.X;
k = size(X,2);
lb = ModelInfo_obj.lb;
ub = ModelInfo_obj.ub;
NegEI = calc_conEI_norm2(xval,ModelInfo_obj);
M = length(ModelInfo_g);
EV = zeros(M,1);
if M>0
    for mm = 1:M
        EV(mm) = calc_conEV_norm2(xval,ModelInfo_g{mm});
    end
end
conNegEI = NegEI/(1+sum(EV));

P = 0;
if flag==1%% Locally makes the ei concave to get rid of the flat objective space
    for ii = 1:k
        P = P + con_fac(ii)*(lb(ii) - xval(ii))*(ub(ii) - xval(ii));
    end
end
f = conNegEI + P;