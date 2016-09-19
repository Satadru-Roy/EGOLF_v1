function [mass,displacement,stress] = tenBarTrussObjective(x)
% This function computes the mass, displacement and the stress for a
% three bar truss problem. [Imperial Units]. the actual calculation of
% the stresses and displacements are done in matlab functions that
% accompy the book "A first course in Finite Elements" by
% Belytschko and Fish. This code is available at
% http://bcs.wiley.com/he-bcs/Books?action=index&itemId=0470035803&bcsId=36
% 25
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

Einit=[1e7 2.8e7 3e7 1.65e7];
rhoInit=[0.10 0.284 0.318 0.162]; % Note: some people use 0.101 for the 
                                  %  density of aluminum
i = 1; 
for count = 1:2:20
   E(i) =  Einit(x(count));
   rho(i) =  rhoInit(x(count));
   A(i) = x(count+1);
   i=i+1;
end

[displacement,stress] = tenBarTruss(A,E);
mass =  calculateWeight(A,rho);
