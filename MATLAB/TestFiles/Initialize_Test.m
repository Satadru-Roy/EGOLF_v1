% Initialize the test function
function [xI_lb,xI_ub,M] = Initialize_Test(num_des,prob)

if prob == 1 %Branin
    %% Set the bounds for the design variables
    xC_lb = 0;xC_ub=15;
    xI_lb=-5;xI_ub=10;
elseif prob == 2 % Griewank
    xI_num = num_des/2;
    xC_num = num_des/2;
    xC_lb=-5*ones(xC_num,1);
    xC_ub=5*ones(xC_num,1);
    xI_lb=-5*ones(xI_num,1);
    xI_ub = 5*ones(xI_num,1);
elseif prob ==3 % 3bar truss
    xC_lb = [1e-10;1e-10;1e-10];
    xC_ub = [10;10;10]; %[cm^2]
    xI_lb = [1;1;1];
    xI_ub = [4;4;4];
elseif prob == 4 % 10bar truss
    xC_lb = 0.1*ones(num_des/2,1);
    xC_ub = 40*ones(num_des/2,1);
    xI_lb = ones(num_des/2,1);
    xI_ub = 4*ones(num_des/2,1);
elseif prob ==5 % Goldstein Function
    xC_lb = [-2];
    xC_ub = [2]; 
    xI_lb = [-2];
    xI_ub = [2];
elseif prob ==6 % 3Hump Camel Function
    xC_lb = [-5];
    xC_ub = [5]; 
    xI_lb = [-5];
    xI_ub = [5];
elseif prob ==7 % Michalewicz Function
    xC_lb = [0];
    xC_ub = [pi]; 
    xI_lb = [0];
    xI_ub = [pi];
elseif prob ==8 % Styblinski-Tang Function
    xC_lb = [-5;-5];
    xC_ub = [5;5]; 
    xI_lb = [-5;-5];
    xI_ub = [5;5];
end

%% Number of constraints that are function of integer/discrete type design variables only
g = confun_Test(xC_lb,xI_lb,1,prob);
M = length(g);
%% Plot 

