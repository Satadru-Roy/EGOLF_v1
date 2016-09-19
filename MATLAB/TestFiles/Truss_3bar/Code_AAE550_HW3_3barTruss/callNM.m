% this file provides input for and calls fminsearch to use the Nelder-Mead
% Simplex method
% Modified on 11/05/07 by Bill Crossley.
clc
close all;
clear all;
format long

options = optimset('Display','none');


%% Part I-(a)-(1)
disp('Part 1-(a)-(1)')
x0 = [1;1];
disp('runs 1')
[xstar,fstar,exitflag,output] = fminsearch('NMfunc',x0,options)

disp('runs 2')
[xstar,fstar,exitflag,output] = fminsearch('NMfunc',xstar,options)


%% Part 1-(a)-(2)
disp('Part 1-(a)-(2)')
x0 = [10;5];
disp('runs 1')
[xstar,fstar,exitflag,output] = fminsearch('NMfunc',x0,options)

disp('runs 2')
[xstar,fstar,exitflag,output] = fminsearch('NMfunc',xstar,options)


%% Part 1-(a)-(3)
disp('Part 1-(a)-(3)')
x0 = [10;10];
disp('runs 1')
[xstar,fstar,exitflag,output] = fminsearch('NMfunc',x0,options)

disp('runs 2')
[xstar,fstar,exitflag,output] = fminsearch('NMfunc',xstar,options)