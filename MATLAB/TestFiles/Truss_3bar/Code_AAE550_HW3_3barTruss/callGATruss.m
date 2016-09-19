%script for HW 3 Part II
close all; clear all; clc;

options = goptions([]);

vlb = [1e-10 1e-10 1e-10 1 1 1];    % Lower bound of each gene - all variables
vub = [5 5 5 4 4 4 ];   % Upper bound of each gene - all variables
bits =[20 20 20 2 2 2]; % Number of bits describing each gene - all variables

l = sum(bits);
pop = 4*l
mutation = (l+1)/(2*pop*l)

options(2) = 0.95;          % Bit affinity
options(11) = pop;          % Population Size
options(12) = 0.5;          % Probability of Crossover
options(13) = mutation;     % Probability of Mutation
options(14) = 1000;         % Max number of generations

[x fbest stats fcount]= GA550('GATrussFunc',[ ],options,vlb,vub,bits);
x
fcount
ngen = size(stats,1);
[phi, f, g] = GATrussFunc(x)


