function  [K,f,d] = preprocessor(A,Ein)
include_flags;
A=A;
Ein=Ein;
% input file to include all variables 
  ten_bar_truss;
%  input_file_example2_8;

% Generate LM array 

for e = 1:nel
    for j = 1:nen
        for m = 1:ndof
            ind = (j-1)*ndof + m;
            LM(ind,e) = ndof*IEN(j,e) - ndof + m;
        end
    end
end
