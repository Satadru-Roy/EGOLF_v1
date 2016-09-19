function f = GATrussFunc_sqp(x,rho)

% Area Variables [m^2]
A = x(1:3)/1e4;

L(1) = sqrt(1.2^2+1.2^2);	% length of element 1 [m]
L(2) = 1.2;               % length of element 2 [m]
L(3) = sqrt(1.2^2+1.2^2);   % length of element 3 [m]

% Objective Functions
f = sum(rho.*A.*L);