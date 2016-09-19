% Input Data for Example 2.2 
nsd 	= 2;	      % Number of space dimensions 
ndof 	= 2;     	  % Number of degrees-of-freedom per node
nnp 	= 6;    	  % Number of nodal points
nel 	= 10;     	  % Number of elements
nen 	= 2;     	  % Number of element nodes
 
neq 	= ndof*nnp;	  % Number of equations
 
f 	= zeros(neq,1);   % Initialize force vector
d 	= zeros(neq,1);   % Initialize displacement matrix
K 	= zeros(neq);     % Initialize stiffness matrix
 
% Element properties
CArea 	= A;   	% Elements area  
leng  	= 360.*[1 1 1 1 1 1 sqrt(2) sqrt(2) sqrt(2) sqrt(2)];   	% Elements length
phi   	= [0 0 0 0 90 90 45 -45 -45 45  ];   	% Angle
E     	= Ein;   	% Young’s Modulus 
 
% prescribed displacements
% displacement     d1x    d1y    d2x    d2y
d           = [0      0      0      0]';
nd 	= 4; 	% Number of prescribed displacement degrees-of-freedom
 
% prescribed forces
f(8)	= -10^5;	   % Force at node 2 in the y-direction
f(12)	= -10^5;       % Force at node 3 in the y-direction
 
% output plots
plot_truss 	= 'yes';
plot_nod	= 'yes';

% mesh generation
ten_bar_mesh;




