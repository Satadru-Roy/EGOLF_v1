clear all
close all
%This script optimizes a ten bar truss. The problem statement is similar to
%Koth,Berke Venkayya, "Comparison  of Optimality Criteria Algorithms
%for Minimum Weight Design of Structures", AIAA 78-469, AIAA Journal, Vol
%17(2), 1978. In detail, the result is equal (or sufficiently close) to 
%case 2 in this paper.

fun = @aluminumObjective;
con = @aluminumConstraints;
x0 = 15* ones(1,10);
lb = 0.1*ones(1,10);
ub = 40 *ones(1,10);
options = optimset;
options.Algorithm = 'active-set';
options.Display = 'iter';
[xOut,f] = fmincon(fun,x0,[],[],[],[],lb,ub,con,options)




