%script for HW 3 Part IV - solution with SQP (fmincon)
close all; clear all; clc;

vlb = [1e-10 1e-10 1e-10];   % Lower bound of area variable
vub = [10 10 10];   % Upper bound of area variables [cm^2]

Aeq = [];
beq = [];

A = [];
b = [];

x0 = [0.1 0.1 0.1];

options = optimset('LargeScale','off','Algorithm','sqp','Display','off');

fbest = 1000;
for i = 1:4
    for j = 1:4
        for k = 1:4
            material = [i j k];
            
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

            [x,fval,exitflag] = fmincon(@(x) GATrussFunc_sqp(x,rho),x0,A,b,Aeq,beq,vlb,vub,@(x) GATrussCon_sqp(x,E,sigma_y),options);
            
            if fval < fbest
                fbest = fval;
                xstar = x;
                best_material = material;
                Ebest = E;
                rhobest = rho;
                sigma_y_best = sigma_y;
                
                Abest = xstar./1e4;
                sigma = stressHW3(Abest,Ebest);
                gbest = abs(sigma)./sigma_y_best - 1;
            end
        end
    end
end

fbest
xstar
best_material
gbest

%% Answer from 2008
% disp('Answer from 2008')
% material = [1 1 3]
% for ii = 1:3
%     if material(ii) == 1      % Aluminum
%         E(ii) = 68.9e9;         % [Pa]
%         rho(ii) = 2700;         % [kg/m^3]
%         sigma_y(ii) = 55.2e6;   % [Pa]
%     elseif material(ii) == 2  % Titanium
%         E(ii) = 116e9;          % [Pa]
%         rho(ii) = 4500;         % [kg/m^3]
%         sigma_y(ii) = 140e6;    % [Pa]
%     elseif material(ii) == 3  % Steel
%         E(ii) = 205e9;          % [Pa]
%         rho(ii) = 7872;         % [kg/m^3]
%         sigma_y(ii) = 285e6;    % [Pa]
%     elseif material(ii) == 4  % Nickel
%         E(ii) = 207e9;          % [Pa]
%         rho(ii) = 8800;         % [kg/m^3]
%         sigma_y(ii) = 59.0e6;   % [Pa]
%     end
% end
% x=[2.0572 0.0001 0.0562]
% f=GATrussFunc_sqp(x,rho);
% f
% A=x/1e4;
% sigma = stressHW3(A,E);
% g = abs(sigma)./sigma_y - 1;
% g

