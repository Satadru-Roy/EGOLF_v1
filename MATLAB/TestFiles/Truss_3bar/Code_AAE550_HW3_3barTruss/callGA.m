% this file provides input variables to the genetic algorithm
% upper and lower bounds, and number of bits chosen for "egg-crate" problem
% Modified on 11/10/09 by Bill Crossley.
clc
% close all;
clear all;
format long

options = goptions([]);
options(1) = 0;
options(14) = 200;

vlb = [0 0];	%Lower bound of each gene - all variables
vub = [10 10];	%Upper bound of each gene - all variables

nbits = [5 15 30];
k=0;

for i = 3:3
    fbest=500;
    rng('shuffle')
    k=0;
%     fprintf('\n%s%d','Number of bits: ',nbits(i))
    bits = nbits(i)*[1 1];	%number of bits describing each gene - all variables
    L = sum(bits);
    popsize = 4*L;
    pmut = (L+1)/(2*L*popsize);
    options(11) = popsize;
    options(13) = pmut;
    for j = 1:3     % run three times
        fprintf('\n%s%d\n','run ', j)
%         fprintf('\n')
        [xstar,fstar,stats,nfit,fgen,lgen,lfit] = GA550('GAfunc',[ ],options,vlb,vub,bits);
        if fstar <= -13.4
            k=k+1;
        end
        ngen = size(stats,1)-1;
%         fprintf('\n%s%d','      Number of bits: ',nbits(i))
%         fprintf('\n%s%d','     Population size: ',popsize)
%         fprintf('\n%s%d','       Mutation rate: ',pmut)
%         fprintf('\n%s%d','  No. of generations: ',ngen)
%         fprintf('\n%s%d','Function evaluations: ',nfit)
        xstar;
        fstar;
        if fstar <= fbest
            xbest = xstar;
            fbest = fstar;
        end
    end
    fprintf('\n%s%d\n','Number of bits: ',nbits(i))
    k
end
xbest
fbest