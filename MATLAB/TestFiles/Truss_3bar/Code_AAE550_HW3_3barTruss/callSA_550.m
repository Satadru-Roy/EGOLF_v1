% this file provides input variables and parameters for simulated annealing
% upper and lower bounds, and initial design chosen for "egg-crate" problem
% Modified on 11/12/02 by Bill Crossley.
clc
close all;
clear all;
format long

bounds = [ 0 10;
           0 10];	% upper and lower bounds for each of the two variables

%% Part I-(b)-(1)
disp('Part 1-(b)-(1)')
X0 = [ 10; 5];	% initial design NOTE: this is a column vector

options = zeros(1,9);		% set up options array for non-default inputs

options(1) = 50;			% initial temperature (default = 50)
options(6) = 0.5;		% cooling rate (default = 0.5)

for i=1:2
    [xstar,fstar,count,accept,oob]=SA_550('SAfunc',bounds,X0,options);
    xstar
    fstar
    count
end


%% Part I-(b)-(2)
disp('Part 1-(b)-(2)')
X0 = [ 10; 5];	% initial design NOTE: this is a column vector

options = zeros(1,9);		% set up options array for non-default inputs

options(1) = 100;			% initial temperature (default = 50)
options(6) = 0.5;		% cooling rate (default = 0.5)

[xstar,fstar,count,accept,oob]=SA_550('SAfunc',bounds,X0,options);
xstar
fstar
count

options(1) = 5;			% initial temperature (default = 50)
[xstar,fstar,count,accept,oob]=SA_550('SAfunc',bounds,X0,options);
xstar
fstar
count


%% Part I-(b)-(3)
disp('Part 1-(b)-(3)')
X0 = [ 10; 5];	% initial design NOTE: this is a column vector

options = zeros(1,9);		% set up options array for non-default inputs

options(1) = 100;			% initial temperature (default = 50)
options(6) = 0.1;		% cooling rate (default = 0.5)

[xstar,fstar,count,accept,oob]=SA_550('SAfunc',bounds,X0,options);
xstar
fstar
count

options(6) = 0.9;		% cooling rate (default = 0.5)
[xstar,fstar,count,accept,oob]=SA_550('SAfunc',bounds,X0,options);
xstar
fstar
count