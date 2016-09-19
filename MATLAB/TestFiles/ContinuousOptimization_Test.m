% Performs the optimization with continuous variables
function [xC_opt,obj,g,eflag,funCount] = ContinuousOptimization_Test(x0I,M,num_des,prob)

[xC_lb,xC_ub] = Initialize_cont_Test(num_des,prob);

g=[];
[n,k] = size(x0I);
xC_opt = zeros(n,length(xC_lb));obj = zeros(n,1);eflag = zeros(n,1);

x0C = (xC_lb + xC_ub)./2;
options=optimset('Algorithm','active-set');
for ii = 1:n
    [xopt,fopt,ef,Output] = fmincon(@(xopt) objfun_Test(xopt,x0I(ii,:),1,prob),x0C,[],[],[],[],xC_lb,xC_ub,@(xopt) confun_Test(xopt,x0I(ii,:),1,prob),options);
    if size(xopt,1)>1
        xopt=xopt';
    end 
    funCount=Output.funcCount;
    xC_opt(ii,:) = xopt;obj(ii) = fopt; eflag(ii) = ef;
    if M>0
        g(ii,:) = confun_Test(xopt,x0I(ii,:),1,prob);
    end
end