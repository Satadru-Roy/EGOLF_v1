function f = aluminumObjective(x)
mat = 4;
f = objectiveTenBarTruss([ mat x(1) mat x(2) mat x(3) mat x(4) mat x(5) ...
    mat x(6) mat x(7) mat x(8) mat x(9) mat x(10)]); 


end