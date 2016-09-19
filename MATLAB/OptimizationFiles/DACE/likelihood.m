function [NegLnLike,Psi,U]=likelihood(x,ModelInfo)

% Inputs:
%	x - vector of log(theta) parameters
%
% Global variables used:
%	ModelInfo.X - n x k matrix of sample locations
%	ModelInfo.y - n x 1 vector of observed data
%
% Outputs:
%	NegLnLike - concentrated log-likelihood *-1 for minimising
%	Psi - correlation matrix
%	U - Choleski factorisation of correlation matrix
%
if size(x,1)>1
    x=x';
end

% global ModelInfo
X=ModelInfo.X;
y=ModelInfo.y;
theta=10.^x;
p=2;  % added p definition (February 10)
n=size(X,1);
one=ones(n,1);

% Pre-allocate memory
Psi=zeros(n,n);
% Build upper half of correlation matrix
for i=1:n
	for j=i+1:n
		Psi(i,j)=exp(-sum(theta.*abs(X(i,:)-X(j,:)).^p)); 
	end
end

% Add upper and lower halves and diagonal of ones plus 
% small number to reduce ill-cokonditioning
Psi=Psi+Psi'+eye(n)+eye(n).*eps; 

% Cholesky factorisation
[U,p]=chol(Psi);

% % Use penalty if ill-conditioned
% if p>0
%     NegLnLike=1e4;
% else
    
%     Sum lns of diagonal to find ln(abs(det(Psi)))
LnDetPsi=2*sum(log(abs(diag(U))));

%     Use back-substitution of Cholesky instead of inverse
mu=(one'*(U\(U'\y)))/(one'*(U\(U'\one)));
SigmaSqr=((y-one*mu)'*(U\(U'\(y-one*mu))))/n;
NegLnLike=-1*(-(n/2)*log(SigmaSqr) - 0.5*LnDetPsi);
% end
% theta
% mu
% SigmaSqr
% U\(U'\(y-one*mu))
% keyboard