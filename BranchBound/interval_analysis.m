% Interval Analysis
function [lb_r, ub_r] = interval_analysis(lb_x,ub_x,ModelInfo)
% The function predicts the lower and the upper bound of the variable 'r'
% from the bounds of the design variable x
% r is related to x by the following equation
% r_i = exp(-sum(theta_h*(x_h - x_h_i)^2))
% For the derivation of each terms see notes (PhD Thesis: Roy)

X=ModelInfo.X;
theta=10.^ModelInfo.Theta;
p = ModelInfo.p;

n=size(X,1); %No. of sample points
k = size(X,2);  % No. of design variables
t1L = zeros(n,k);t1U = zeros(n,k);
t2L = zeros(n,k);t2U = zeros(n,k);
t3L = zeros(n,k);t3U = zeros(n,k);
t4L = zeros(n,1);t4U = zeros(n,1);
t5L = zeros(n,1);t5U = zeros(n,1);
lb_r = zeros(n,1);ub_r = zeros(n,1);
eterm = 1;
if mod(p,2) == 0
    for i = 1:n
        for h = 1:k
            t1L(i,h) = lb_x(h) - X(i,h);
            t1U(i,h) = ub_x(h) - X(i,h);

            t2L(i,h) = max([0,min([t1L(i,h)*t1L(i,h),t1L(i,h)*t1U(i,h),t1U(i,h)*t1U(i,h)])]);
            t2U(i,h) = max([0,max([t1L(i,h)*t1L(i,h),t1L(i,h)*t1U(i,h),t1U(i,h)*t1U(i,h)])]);

            t3L(i,h) = min([(-theta(h)*t2L(i,h)),(-theta(h)*t2U(i,h))]);
            t3U(i,h) = max([(-theta(h)*t2L(i,h)),(-theta(h)*t2U(i,h))]);
        end

        t4L(i,1) = sum(t3L(i,:));
        t4U(i,1) = sum(t3U(i,:));

%         t5L(i,1) = min(-t4L(i,1),-t4U(i,1));
%         t5U(i,1) = max(-t4L(i,1),-t4U(i,1));

        lb_r(i,1) = exp(t4L(i,1));
        ub_r(i,1) = exp(t4U(i,1));
    end
else
    fprintf('\n%s','Warning! Value of p should be 2. Cannot perform interval analysis')
    fprintf('\n%s','Returning global bound of the r variable')
    eterm = 0;
%             t2L(i,h) = min([t1L(i,h)*t1L(i,h),t1L(i,h)*t1U(i,h),t1U(i,h)*t1U(i,h)]);
end
if eterm == 0
    lb_r = zeros(n,1);
    ub_r = ones(n,1);
end