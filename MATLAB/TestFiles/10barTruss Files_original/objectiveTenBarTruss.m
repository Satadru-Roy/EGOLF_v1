function f = objectiveTenBarTruss(x)
% Calculates the weight of the ten bar truss with 4 materials.
% The input argument must have the form:
%   [MaterialElement1 AreaElement1 MaterialElement2 AreaElement2...
%   ... MaterialElement10 AreaElement 10]
% where 'MaterialElementX' is an integer that represents the material for
% element X and 'AreaElementX' represents the area of element X.
% The following materials are avalable:
% 1 : Aluminum
% 2: titanium
% 3: Steel
% 4: Nickel

if(length(x)~=20)
    error('Vector needs 20 inputs arguments');
end
rhoInit=[0.10 0.284 0.318 0.162]; % Note: some people use 0.101 for the 
                                  %  density of aluminum
i = 1; 
for count = 1:2:20 
   rho(i) =  rhoInit(x(count));
   A(i) = x(count+1);
   i=i+1;
end
f = calculateWeight(A,rho);