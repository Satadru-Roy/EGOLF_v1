function [phi f g] = GATrussFunc(x)

% Material Variables
material = x(4:6);

% Possible material properties for each element
for ii = 1:3
    if material(ii) == 1      % Aluminum
        E(ii) = 68.9e9;         % [Pa]
        rho(ii) = 2700;         % [kg/m^3]      
        sigma_y(ii) = 55.2e6;   % [Pa]
    elseif material(ii) == 2  % Titanium
        E(ii) = 116e9;          % [Pa]
        rho(ii) = 4500;         % [kg/m^3]
        sigma_y(ii) = 140e6;    % [Pa]
    elseif material(ii) == 3  % Steel
        E(ii) = 205e9;          % [Pa]
        rho(ii) = 7872;         % [kg/m^3]
        sigma_y(ii) = 285e6;    % [Pa]
    elseif material(ii) == 4  % Nickel
        E(ii) = 207e9;          % [Pa]
        rho(ii) = 8800;         % [kg/m^3]
        sigma_y(ii) = 59.0e6;   % [Pa]
    end
end

% Area Variables [m^2]
A = x(1:3)/1e4;

% Lengths of Elements
L(1) = sqrt(1.2^2+1.2^2);	% length of element 1 [m]
L(2) = 1.2;               % length of element 2 [m]
L(3) = sqrt(1.2^2+1.2^2);   % length of element 3 [m]

% Objective Functions
f = sum(rho.*A.*L);

% Stress Values
sigma = stressHW3(A,E);

% Constraints
g = abs(sigma) ./ sigma_y - 1;

%step-linear penalty function
P = 0;
for i = 1:length(g)
    if g(i) >= 0
        P = P + 1e6*(1+g(i));
    end
end

% Fitness Function
phi = f + P;