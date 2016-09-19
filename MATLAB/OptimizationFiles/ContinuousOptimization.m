% Performs the optimization with continuous variables
function [xC_opt,obj,g,eflag,funCount] = ContinuousOptimization(x0I,M,num_route)

%% Aircarft uncertain parameters
FileName = 'AC_new';
acinputs.CD0_mult          = 1;
acinputs.e_mult            = 1;
acinputs.TSFC_mult         = 1;
acinputs.pax_wt             = 165;
acinputs.check=0;
%%
[n,xI_num] = size(x0I);
lookup = lookup_ref();

% network = networkdata(num_route);
% num_route = length(network.range);
existAC = existingAC();
num_ac = length(existAC.AC_name)+1;

% xI_num_alloc = num_route*num_ac;

%% This adds a top level profile which would be design variable in the
%% future release
x0I_toplev = 10*ones(n,1);
x0I = [x0I_toplev,x0I];
%%

g=[];
xC_opt = []; %zeros(n,length(xC_lb));
obj = zeros(n,1);eflag = zeros(n,1);
load Model
options=optimset('Algorithm','sqp','Display','iter','MaxFunEvals',5000);
for ii = 1:n
    trip = x0I(ii,2:end);
    code_profile = x0I(ii,1);
    [xC_lb, xC_ub, xC0, num_des] = Intialize_cont(code_profile, num_route, FileName, acinputs); 
    lookup.xC_lb=xC_lb;lookup.xC_ub=xC_ub;lookup.xC0 = xC0;
    lookup.num_des = num_des;
    lookup.num_route = num_route;
    [Ain_hat,bin_hat,Ain,bin] = lincon_DesAlloc(trip,lookup);
    xC0_hat = (xC0 - xC_lb)./(xC_ub-xC_lb);
    %%
%     keyboard
    xC_hat_lb = zeros(length(xC_lb),1);
    xC_hat_ub = ones(length(xC_ub),1);
    fprintf('\n%s\n','Optimizing the aircraft for the given discrete inputs')
    [xopt_hat,fopt,ef,Output] = fmincon(@(xopt_hat) objfun_DesAlloc(xopt_hat,x0I(ii,:),lookup,Model,acinputs),...
        xC0_hat,Ain_hat,bin_hat,[],[],xC_hat_lb,xC_hat_ub,@(xopt_hat) con_DesAlloc(xopt_hat,x0I(ii,:),lookup,Model,acinputs),options);
       
    %% Convert back to the original data
    if size(xopt_hat,1)>1
        xopt = xC_lb + xopt_hat.*(xC_ub - xC_lb);
        xopt_hat = xopt_hat';
    else
        xopt = xC_lb + xopt_hat'.*(xC_ub - xC_lb);
    end
    xopt = xopt';
    %%
    
    funCount=Output.funcCount;
    xC_opt = [xC_opt;xopt]; obj(ii) = fopt; eflag(ii) = ef;
    if M>0
        g1 = con_DesAlloc(xopt_hat,x0I(ii,:),lookup,Model,acinputs);
        g2 = Ain*xopt' - bin;
        g(ii,:) = [g1',g2'];
    end
%     keyboard
end
