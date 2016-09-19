%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2D Truss (Chapter 2)        %
% Haim Waisman, Rensselaer    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clear all;
%close all; 
function [d,f_E,stress]= truss(A,Ein) 
% include global variables
include_flags;  
Ein=Ein;
% Preprocessor Phase 
 [K,f,d]	= preprocessor(A,Ein);

% Calculation and assembly of element matrices
for e = 1:nel
    ke	= trusselem(e);
    K	= assembly(K,e,ke);
end

% Solution Phase
 [d,f_E]	= solvedr(K,f,d);
 
% Postprocessor Phase 
stress=postprocessor(d);
