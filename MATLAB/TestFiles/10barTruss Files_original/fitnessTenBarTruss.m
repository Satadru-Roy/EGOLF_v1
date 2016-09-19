function f = fitnessTenBarTruss(x)



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

% Set the penalty multiplier

if(length(x)~=20)
    error('Vector needs 20 inputs arguments');
end
penaltyMultiplier = 10000;

[c] = constraintsTenBarTruss(x);
f =  objectiveTenBarTruss(x);
for i=1:length(c)
   if c(i)>1e-6      
      f = f + penaltyMultiplier*c(i); 
   end
end