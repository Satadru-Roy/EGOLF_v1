function phi = NMfunc(x)
% objective function: Langermann function
% linear exterior penalty to enforce bounds

m = 10;

d = 2;

c = [ 3 10  2  6  1 10  1  3  5  7]';

A = [ 8  7  6  1  4  6  6  9  5  1;
      4  6 10  7  4  4  5  5  9 10]';

f = 0;
for i = 1:m
    x_A = 0;
    for j = 1:d
        x_A = x_A + (x(j)-A(i,j))^2;
    end
    f = f + c(i) * exp(-x_A/pi) * cos(pi*x_A);
end

g(1) = -x(1);  % enforces lower bound
g(2) = -x(2);
g(3) = x(1)/10 - 1; % enforces upper bound
g(4) = x(2)/10 - 1;

P = 0;    % initialize penalty function
for i = 1:4
    P = P + 10 * max(0,g(i));  % use c_j = 10 for all bounds
end
phi = f + P;