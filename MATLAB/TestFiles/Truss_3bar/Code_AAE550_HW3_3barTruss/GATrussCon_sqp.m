function [g,h] = GATrussCon_sqp(x,E,sigma_y)

% Area Variables [m^2]
A = x(1:3)/1e4;

% Stress Values
sigma = stressHW3(A,E);

% Constraints
g = abs(sigma) ./ sigma_y - 1;

h = [];