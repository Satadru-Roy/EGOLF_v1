function [displacement,stress] = tenBarTruss(A,E)
if(length(A)~=10)
   error('Need 10 areas for the 10 elements'); 
end
if(length(E)~=10 )
   error('Need 10 Youngs moduli for the 10 elements'); 
end

l=360.*[1 1 1 1 1 1 sqrt(2) sqrt(2) sqrt(2) sqrt(2)];
[displacement,f_E,stress]=truss(A,E);


