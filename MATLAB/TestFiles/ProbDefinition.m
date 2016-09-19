%% Choose your test problem
function [prob,num_des,FileNam] = ProbDefinition()
%%
% 1. Branin Function (1-int, 1-cont) - 2 D
% 2. Griewank Function (N/2-int, N/2-cont) - N D
% 3. 3-Bar Truss problem (3-int, 3-cont) - 6 D
% 4. 10-Bar Truss problem (10-discrete, 10-cont) - 20 D
% 5. Goldstein-Price Function (1-int, 1-cont) - 2 D
% 6. 3Hump camel Function (1-int, 1-cont) - 2 D
% 7. Michalewicz Function (1-int, 1-cont) - 2 D
%%
prob = 3;
switch prob
    case 1
        num_des = 2;
        FileNam = 'TestCase_Branin';
    case 2
        num_des = 2;
        FileNam = 'TestCase_Griewank';
    case 3
        num_des = 6;
        FileNam = 'TestCase_3barTruss';
    case 4
        num_des = 20;
        FileNam = 'TestCase_10barTruss';
    case 5
        num_des = 2;
        FileNam = 'TestCase_Goldstein';
    case 6
        num_des = 2;
        FileNam = 'TestCase_3HumpCamel';
    case 7 
        num_des = 2;
        FileNam = 'TestCase_Michalewicz';
    case 8
        num_des = 4;
        FileNam = 'TestCase_Styblinski-Tang';
end
