%% Expected Improvement
% This function evaluates the expected improvement in the normalized design
% space
function [f] = calc_conEV_norm2(xval,ModelInfo_g,gSSqr_inp,g_hat_inp)
%% Inputs
% xval - Design variables in the normalized design space
% y_best - Best objective value so far from the sampled points
% ModelInfo_obj - Details about the objective surrogate model
% con_fac - Concavity factor as to what extent the flat surface needs to be
% concaved.
% flag - [0/1] 0 - means no concavity, returns the orginal function, 1 -
% returns the function value accounting for the concavity factor
% SSqr_inp - User supplied SSqr
% y_hat_inp - User supplied y_hat
% If the last two inputs are provided, function returns the EI value evaluated using
% the user input SSqr and y_hat
%% Output
% f - returns the expected improvement function value evaluated at xI
% if flag=0, returns the actual expected improvement value
% if flag=1, returns the expected improvement value with the concavity
% factor added to it.
%%
g_min = 0;
X=ModelInfo_g.X;
if nargin <= 2
    y=ModelInfo_g.y;
    theta=10.^ModelInfo_g.Theta;
    p=ModelInfo_g.p;  
    R=ModelInfo_g.Psi; %U=ModelInfo.U;
    SigmaSqr = ModelInfo_g.SigmaSqr;
    mu = ModelInfo_g.mu;
    n=size(X,1);
    
    one=ones(n,1);
    r = ones(n,1);
    if size(xval,2)<2
        xval = xval';
    end
    x = xval;
    % fill r vector
    for i=1:n
        r(i,1)=exp(-sum(theta.*(x-X(i,:)).^p));
    end
    %Calculate prediction and error
    g_hat = mu + r'*(R\(y-one*mu));
    gSSqr = SigmaSqr*(1 - (r'*(R\r)) + ((1-one'*(R\r))^2/(one'*(R\one))));
else
    g_hat = g_hat_inp;
    gSSqr = gSSqr_inp;
end

if gSSqr <= 0
    EV = 0;
else
    %Calculate the expected violation
    ei_1 = (g_hat-g_min)*(0.5 + 0.5*erf((1/sqrt(2))*((g_hat - g_min)/sqrt(abs(gSSqr)))));
    ei_2 = sqrt(abs(gSSqr))*(1/sqrt(2*pi))*exp(-(1/2)*((g_hat - g_min)^2/abs(gSSqr)));
    EV = (ei_1 + ei_2);
end

f = EV;
if isnan(f)==1
    fprintf('\n%s','Warning! Expected improvement is NAN')
    f=0;
end