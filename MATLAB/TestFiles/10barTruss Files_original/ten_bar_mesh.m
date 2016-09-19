function  ten_bar_mesh
include_flags;


% Node:  1      2      3       4       5       6  (origin placed at node 2) 
%---------------------------------------------------
x   =  [  0.0    0.0  360.0   360.0   720.0 720.0  ];     % X coordinate  
y   =  [360.0    0.0  360.0     0.0   360.0   0.0  ];     % Y coordinate

% connectivity array
IEN =  [1   3   2   4   4   6   2   1   3   4   
        3   5   4   6   3   5   3   4   6   5];     

%plot truss
%plottruss;

