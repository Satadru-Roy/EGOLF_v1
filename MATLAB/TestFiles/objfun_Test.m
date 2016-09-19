%Define the objective function here
function f = objfun_Test(x,xI,flag,prob)

if prob == 1 
    %Branin function (2D)
    if flag==1 %EGO Framework
        X1 = xI;
        X2 = x;
    elseif flag==2 %Other tools
        X1 = x(1);
        X2 = x(2);
    end
    a = 1;b = 5.1/(4*pi^2);c = 5/pi;d = 6;e = 10;ff = 1/(8*pi); 
    f = (a*( X2 - b*X1^2 + c*X1 - d )^2 + e*(1-ff)*cos(X1) + e); 
elseif prob==2
    % Griewank's function
    if flag == 1
        xC_num = length(x);
        xI_num = length(xI);
        f1C=0;f2C=1;
        for ii = 1:xC_num
            f1C = f1C + x(ii)^2/4000;
            f2C = f2C*cos(x(ii)/sqrt(ii));
        end
        f1I=0;f2I=1;
        for ii = 1:xI_num
            f1I = f1I + xI(ii)^2/4000;
            f2I = f2I*cos(xI(ii)/sqrt(ii));
        end
        f = 1*((f1C+f1I)-(f2C*f2I) + 1);
    elseif flag==2
        f1=0;f2=1;
        for ii=1:length(x)
            f1 = f1 + x(ii)/4000;
            f2 = f2*cos(x(ii)/sqrt(ii));
        end
        f = 1*(f1 - f2 + 1);
    end
elseif prob == 3 % 3 bar truss
   f = TrussFunc_sqp(x,xI); 
elseif prob == 4 % 10 bar truss
    cc=1;x_combined=zeros(20,1);
    for ii = 1:10
        x_combined(cc) = xI(ii);
        x_combined(cc+1) = x(ii);
        cc=cc+2;
    end
   f= objectiveTenBarTruss(x_combined);
elseif prob == 5 %Goldstein Function
    if flag==1 %EGO Framework
        X1 = xI;
        X2 = x;
    elseif flag==2 %Other tools
        X1 = x(1);
        X2 = x(2);
    end 
    f = (1 + (X1 + X2 + 1)^2*(19 - 14*X1 + 3*X1^2 - 14*X2 + 6*X1*X2 + 3*X2^2))*...
        (30 + (2*X1 - 3*X2)^2*(18 - 32*X1 + 12*X1^2 + 48*X2 - 36*X1*X2 + 27*X2^2));
elseif prob == 6 %3Hump Camel Function
    if flag==1 %EGO Framework
        X1 = xI;
        X2 = x;
    elseif flag==2 %Other tools
        X1 = x(1);
        X2 = x(2);
    end 
    f = 2*X1^2 - 1.05*X1^4 + (X1^6/4) + X1*X2 + X2^2;
elseif prob == 7 %Michalewicz Function
    if flag==1 %EGO Framework
        X1 = xI;
        X2 = x;
    elseif flag==2 %Other tools
        X1 = x(1);
        X2 = x(2);
    end 
    m=10;
    f = -(sin(X1)*(sin(1*X1^2/pi))^(2*m) + sin(X2)*(sin(2*X2^2/pi))^(2*m));
elseif prob == 8 %Stylblinski-Tang Function
    if flag == 1
        xC_num = length(x);
        xI_num = length(xI);
        fI=0;fC=0;
        for ii = 1:xI_num
            fI = fI + (xI(ii)^4 - 16*xI(ii)^2 + 5*xI(ii));
        end
        for ii = 1:xC_num
            fC = fC + (x(ii)^4 - 16*x(ii)^2 + 5*x(ii));
        end
        f = 0.5*(fI + fC);
    end
end
    