% This function defines the constraints of the problem
function [g,h] = confun_Test(x, xI, flag,prob)
g=[];h=[];
if prob == 1 % General Constraints
    if flag == 1
        x1 = xI;
        x2 = x;
    else
        x1 = x(1);
        x2 = x(2);
    end
%     g(1,1) = x1*x2;
%     g(1,1) = 1 - ((x1+1)^2 + (x2-10)^2)/3^2;
%     g(1,2) = 1 - (((x1-2)/2.5)^2 + ((x2-2.5)/0.5)^2);
%     g(1,3) = 1 - (((x1-10)/2.5)^2 + ((x2-5)/4)^2);
elseif prob == 2
    if flag==1
        
%     g(1,1) = 1 - ((x1+1)^2 + (x2-10)^2)/3^2;
%     g(1,2) = 1 - (((x1-2)/2.5)^2 + ((x2-2.5)/0.5)^2);
%     g(1,3) = 1 - (((x1-10)/2.5)^2 + ((x2-5)/4)^2);
    end
elseif prob ==3 %3 bar truss
   [g,h] = TrussCon_sqp(x,xI); 
elseif prob == 4 %10 bar truss
        cc=1;x_combined=zeros(20,1);
    for ii = 1:10
        x_combined(cc) = xI(ii);
        x_combined(cc+1) = x(ii);
        cc=cc+2;
    end
    [g,h] = constraintsTenBarTruss(x_combined);
end
end
    

