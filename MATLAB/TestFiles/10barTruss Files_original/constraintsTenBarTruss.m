function [c ceq] = constraintsTenBarTruss(x)
% Calculates the constraints for the ten bar truss problem.
% These are the 10 stresses and the max. displacement of a node
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

%Sets the maximum displacement
displacementMax = 2;
%Sets the maximum stresses
stress_max=[25 20 13.3 73.3]*10^3;


Einit=[1e7 2.8e7 3e7 1.65e7];

i = 1; 
for count = 1:2:20
   E(i) =  Einit(x(count));  
   A(i) = x(count+1);
   i=i+1;
end
[displacement,stress] = tenBarTruss(A,E);
for i=1:10
    c(i) =  abs(stress(i))./stress_max(x(2*i-1))-1;
end
c(11) = max(abs(displacement))/displacementMax-1;
ceq = [];