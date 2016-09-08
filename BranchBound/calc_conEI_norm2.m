%% Expected Improvement
% This function evaluates the expected improvement in the normalized design
% space
function [f] = calc_conEI_norm2(xval,ModelInfo_obj,SSqr_inp,y_hat_inp)
%% Inputs
% xval - Design variables in the normalized design space
% y_best - Best objective value so far from the sampled points
% ModelInfo_obj - Details about the objective surrogate model
% SSqr_inp - User supplied SSqr
% y_hat_inp - User supplied y_hat
% If the last two inputs are provided, function returns the EI value evaluated using
% the user input SSqr and y_hat
%% Output
% f - returns the expected improvement function value evaluated at xval
%%
y_min = ModelInfo_obj.y_best;
X=ModelInfo_obj.X;
if nargin <= 2
    y=ModelInfo_obj.y;
    theta=10.^ModelInfo_obj.Theta;
    p=ModelInfo_obj.p;  
    R=ModelInfo_obj.Psi; %U=ModelInfo.U;
    SigmaSqr = ModelInfo_obj.SigmaSqr;
    mu = ModelInfo_obj.mu;
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
    y_hat = mu + r'*(R\(y-one*mu));
    SSqr = SigmaSqr*(1 - (r'*(R\r)) + ((1-one'*(R\r))^2/(one'*(R\one))));
else
    y_hat = y_hat_inp;
    SSqr = SSqr_inp;
end

if SSqr == 0
    NegEI = 0;
else   
    %Calculate the expected improvement
    ei_1 = (y_min - y_hat)*(0.5 + 0.5*erf((1/sqrt(2))*((y_min - y_hat)/sqrt(abs(SSqr)))));
    ei_2 = sqrt(abs(SSqr))*(1/sqrt(2*pi))*exp(-(1/2)*((y_min - y_hat)^2/abs(SSqr)));
    NegEI = -(ei_1 + ei_2);
end
f = NegEI;
% keyboard
if isnan(f)==1
    fprintf('\n%s','Warning! Expected improvement is NAN')
    f=0;
    keyboard
end