%% This function estimates the concavity factor for the expected improvment function
% This factor help make the ei concave and provides the needed
% gradient information for the SQP algorithm.
function [con_fac] = concave_factor(xL,xU)
K = length(xL);
per_htm = 5; % in percentage [%]
con_fac = zeros(K,1);        
for k = 1:K
    if abs(xL(k)-xU(k)) < 1e-10
        con_fac(k,1) = 0;
    else
        h_req = (per_htm/100)*(xU(k)-xL(k)); 
        xm = (xL(k)+xU(k))/2;
        h_act = (xm(k)-xL(k))*(xm(k)-xU(k));
        con_fac(k,1) = h_req/h_act;
    end
end