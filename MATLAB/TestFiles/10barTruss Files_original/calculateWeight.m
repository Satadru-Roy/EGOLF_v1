function weight = calculateWeight(A,rho)
% This function caculates the weight of a ten bar truss according to
% the input areas (A) and the input densities (rho)

l=360.*[1 1 1 1 1 1 sqrt(2) sqrt(2) sqrt(2) sqrt(2)];

weight =  sum(l.*A.*rho);