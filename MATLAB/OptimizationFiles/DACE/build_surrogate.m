% Build Surrogate - Using Krigging
function [ModelInfo] = build_surrogate(k,ModelInfo)
 
% global ModelInfo
p = 2;

%Set upper and lower bounds for the search of log theta
UpperTheta = ones(1,k).*2;
LowerTheta = ones(1,k).*-3;

%Run GA search of likelihood
%GA options
% options=[];
% options = gaoptimset('PopulationSize',250,'EliteCount',3,...
%     'Generations',500); %,'UseParallel','always');
% [Theta, MinNegLnLikelihood] = ga(@(Theta) likelihood(Theta,ModelInfo),k,[],[],[],[],LowerTheta, UpperTheta, [], options);

% Using fmincon
x0 = LowerTheta + 0.5*(UpperTheta - LowerTheta);
[Theta] = fmincon(@(Theta) likelihood(Theta,ModelInfo),x0, [],[],[],[],LowerTheta,UpperTheta);

[~, R, ModelInfo.U] = likelihood(Theta, ModelInfo);
ModelInfo.Theta = 10.^(Theta');
ModelInfo.p = p;
%Calculate mu and Sigma^2
y = ModelInfo.y;
X = ModelInfo.X;
n=size(X,1);
one=ones(n,1);
ModelInfo.R_inv = inv(R);
mu=(one'*(ModelInfo.R_inv*y))/(one'*(ModelInfo.R_inv*one));
SigmaSqr = ((y-one*mu)'*(ModelInfo.R_inv*(y-one*mu)))/n;

ModelInfo.SigmaSqr = SigmaSqr;
ModelInfo.mu = mu;
% keyboard