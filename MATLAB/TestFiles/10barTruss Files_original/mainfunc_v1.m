clc;
clear all
close all

%Solved as NLP, material supplied
%     fun = @aluminumObjective;
%     con = @aluminumConstraints;
%     x0 = 15* ones(1,10);
%     lb = 0.1*ones(1,10);
%     ub = 40 *ones(1,10);
%     options = optimset;
%     options.Algorithm = 'active-set';
%     options.Display = 'iter';
%     options.MaxFunEvals = 5000;
%     [xOut,f] = fmincon(fun,x0,[],[],[],[],lb,ub,con,options)

%Solved as MDNLP using GA
    for ii = 1:2:19
        lb(ii) = 1; 
        lb(ii+1) = 0.1;
        ub(ii) = 4; 
        ub(ii+1) = 40;
    end
    [xOut_GA,f_GA,eflag_GA,Output] = ga(@objectiveTenBarTruss,length(lb),[],[],[],[],...
    lb,ub,@constraintsTenBarTruss,[1:2:19])
    g = constraintsTenBarTruss(xOut_GA)
    Output.funccount