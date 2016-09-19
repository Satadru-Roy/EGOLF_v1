%% Linearized under-estimator
% This version scales up artificial variable r to [0,1] space
function [Ain, bin] = lin_underestimator2(lb, ub, ModelInfo_obj)
X=ModelInfo_obj.X;
% theta=10.^ModelInfo_obj.Theta;
theta = ModelInfo_obj.Theta;
p=ModelInfo_obj.p;

n=size(X,1); %No. of sample points
k = size(X,2);  % No. of design variables

lb_x = lb(1:k);
ub_x = ub(1:k);
lb_r = lb(k+1:end);
ub_r = ub(k+1:end);

% %To avoid any numerical issues
% aa = find(abs(lb_r-0)<1e-300); lb_r(aa) = 1e-300;
% aa = find(abs(ub_r-0)<1.1e-300); ub_r(aa) = 1.1e-300;
% bb = find(abs(lb_x-0)<1e-300); lb_x(bb) = 1e-300;
% bb = find(abs(ub_x-0)<1.1e-300); ub_x(bb) = 1.1e-300;
% Constraints of the MaxSSqr problem
%   T1                    T2
% ln(r_i) + sum_h(theta_h*(x_h - x_h_i)^p) <= 0 ...........(1)
%                       [assume p=2 don't change p]
%   T3                    T4
% -ln(r_i) - sum_h(theta_h*(x_h - x_h_i)^p) <= 0 ..........(2)
a1 = zeros(n,n);a3 = zeros(n,n);
a1_hat = zeros(n,n);a3_hat = zeros(n,n);
a2 = zeros(n,k);a4 = zeros(n,k);
b1 = zeros(n,1);b3 = zeros(n,1);
b1_hat = zeros(n,1);b3_hat = zeros(n,1);
b2 = zeros(n,k);b4 = zeros(n,k);

for i = 1:n

    %T1: Linearize under-estimator of ln(r_i) = a1(i,i)*r(i) + b1(i);
    if ub_r(i) <= lb_r(i)
        a1(i,i) = 0;
    else
        a1(i,i) = ((log(ub_r(i)) - log(lb_r(i)))/(ub_r(i) - lb_r(i)));
    end
    a1_hat(i,i) = a1(i,i)*(ub_r(i)-lb_r(i));
    b1(i,1) = log(ub_r(i)) - a1(i,i)*ub_r(i); 
    b1_hat(i,1) = a1(i,i)*lb_r(i) + b1(i,1);
    
    %T3: Linearize under-estimator of -ln(r_i) = a3(i,i)*r(i) + b3(i,1);
    r_m_i = (lb_r(i) + ub_r(i))/2;
    a3(i,i) = -1/r_m_i; 
    a3_hat(i,i) = a3(i,i)*(ub_r(i)-lb_r(i));
    b3(i,1) = -log(r_m_i) - a3(i,i)*r_m_i;
    b3_hat(i,1) = a3(i,i)*lb_r(i) + b3(i,1);
    
    for h = 1:k      
        %T2: Linearize under-estimator of theta_h*(x_h - x_h_i)^2 = a4(i,h)*x_h(h) + b4(i,h);
        x_m_h = (ub_x(h) + lb_x(h))/2;
        a2(i,h) = p*theta(h)*(x_m_h - X(i,h))^(p-1);
        yy = theta(h)*(x_m_h - X(i,h))^p;
        b2(i,h) = -a2(i,h)*x_m_h + yy;
    
        %T4: Linearize under-estimator of -theta_h*(x_h - x_h_i)^2 = a4(i,h)*x_h(h) + b4(i,h);
        yy2 = -theta(h)*(ub_x(h) - X(i,h))^p;
        yy1 = -theta(h)*(lb_x(h) - X(i,h))^p;
        if ub_x(h) <= lb_x(h)
            a4(i,h) = 0;
        else
            a4(i,h) = (yy2 - yy1)/(ub_x(h) - lb_x(h));
        end
        b4(i,h) = -a4(i,h)*lb_x(h) + yy1;
    end
end
Ain = [[a2;a4],[a1_hat;a3_hat]];
bin = [-(b1_hat + sum(b2,2));-(b3_hat + sum(b4,2))];   

% %% plot the original constraints and their underestimators
% SSE = 0;
% num_int = 5;
% fig_num=10;
% for ii = 1:n
%     %Plot T1 and T3
%     figure(fig_num+ii)
%     r_res = (ub_r(ii) - lb_r(ii))/num_int;
%     rran = [lb_r(ii):r_res:ub_r(ii)]';
%     f_act_T1 = zeros(length(rran),1);
%     f_act_T3 = zeros(length(rran),1);
%     f_app_T1 = zeros(length(rran),1);
%     f_app_T3 = zeros(length(rran),1);
%     err_T1 = zeros(length(rran),1);
%     err_T3 = zeros(length(rran),1);
%     for yy = 1:length(rran)
%         f_act_T1(yy,1) = log(rran(yy));
%         f_app_T1(yy,1) = a1(ii,ii)*rran(yy) + b1(ii);
%         err_T1(yy,1) = (f_act_T1(yy,1) - f_app_T1(yy,1))^2;
%         f_act_T3(yy,1) = -log(rran(yy));
%         f_app_T3(yy,1) = a3(ii,ii)*rran(yy) + b3(ii);
%         err_T3(yy,1) = (f_act_T3(yy,1) - f_app_T3(yy,1))^2;
%         if isnan(f_act_T1(yy,1)) == 1
%             err_T1(yy,1) = 0;
%         end
%         if isnan(f_act_T3(yy,1)) == 1
%             err_T3(yy,1) = 0;
%         end
%     end
%     subplot(k+1,2,1),plot(rran,f_act_T1,'k-',rran,f_app_T1,'r-')  
%     xlabel('r_1'); ylabel('log(r_1)')
%     subplot(k+1,2,2),plot(rran,f_act_T3,'k-',rran,f_app_T3,'r-') 
%     xlabel('r_1'); ylabel('- log(r_1)')
%     SSE = SSE + sum(err_T1) + sum(err_T3);
%     
%     %Plot T2 & T4  
%     for jj = 1:k
%         x_res = (ub_x(jj) - lb_x(jj))/num_int;
%         xran = [lb_x(jj):x_res:ub_x(jj)]';
%         f_act_T2 = zeros(length(xran),k);
%         f_act_T4 = zeros(length(xran),k);
%         f_app_T2 = zeros(length(xran),k);
%         f_app_T4 = zeros(length(xran),k);
%         err_T2 = zeros(length(xran),k);
%         err_T4 = zeros(length(xran),k);
%         for yy = 1:length(xran)
%             f_act_T2(yy,jj) = theta(jj)*(xran(yy) - X(ii,jj))^p;
%             f_app_T2(yy,jj) = a2(ii,jj)*xran(yy) + b2(ii,jj);
%             err_T2(yy,jj) = (f_act_T2(yy,jj) - f_app_T2(yy,jj))^2;
%             f_act_T4(yy,jj) = -theta(jj)*(xran(yy) - X(ii,jj))^p;
%             f_app_T4(yy,jj) = a4(ii,jj)*xran(yy) + b4(ii,jj);
%             err_T4(yy,jj) = (f_act_T4(yy,jj) - f_app_T4(yy,jj))^2;
%         end
%         subplot(k+1,2,3),plot(xran,f_act_T2(:,jj),'k-',xran,f_app_T2(:,jj),'r-') 
%         xlabel('x'); ylabel('Theta_1(x - x^(^1^))^2')
%         subplot(k+1,2,4),plot(xran,f_act_T4(:,jj),'k-',xran,f_app_T4(:,jj),'r-') 
%         xlabel('x'); ylabel('- Theta_1(x - x^(^1^))^2')
%         SSE = SSE + sum(err_T2(:,jj)) + sum(err_T4(:,jj));
%     end
% end 
% RMSE = sqrt(SSE);
% fig_num = fig_num+n;

% %% Plot the constraint (design varaibles: x1, r1, r2, r3)
% fig_num =10;
% num_int = 11;
% for ii = 1:n
%     %Plot T1 and T3
%     r_res = (ub_r(ii) - lb_r(ii))/num_int;
% %     r_res = 0.05;
%     rran = [lb_r(ii):r_res:ub_r(ii)];
%     
%     x_res = (ub_x - lb_x)/num_int;
% %     x_res = 0.05;
%     xran = [lb_x:x_res:ub_x];
%     
%     
%     for ry = 1:length(rran)
%         for xy = 1:length(xran)
%             g1_act(xy,ry) = log(rran(ry)) + (theta*(xran(xy) - X(ii,1))^p);
%             g1_app(xy,ry) = (a1(ii,ii)*rran(ry) + b1(ii,1)) + (a2(ii,1)*xran(xy) + sum(b2(ii),2));
%             g2_act(xy,ry) = -log(rran(ry)) - (theta*(xran(xy) - X(ii,1))^p);
%             g2_app(xy,ry) = (a3(ii,ii)*rran(ry) + b3(ii,1)) + (a4(ii,1)*xran(xy) + sum(b4(ii),2));       
%         end
%     end
%     fprintf('\n%s%d\n','Sample point ',ii)
%     fprintf('%s%f%s%f\n','Lowest actual g1 ',min(min(g1_act)),' and lowest approx g1 ',min(min(g1_app)))
%     fprintf('%s%f%s%f\n','Lowest actual g2 ',min(min(g2_act)),' and lowest approx g1 ',min(min(g2_app)))
%     fig_num = fig_num + 1;
%     figure(fig_num);
%     subplot(1,2,1);surf(rran,xran,g1_act);hold on
%     subplot(1,2,1);surf(rran,xran,g1_app);hold on 
%     subplot(1,2,2);surf(rran,xran,g2_act);hold on
%     subplot(1,2,2);surf(rran,xran,g2_app);hold on 
%     xlabel('r_i')
%     ylabel('x')
% end 

% keyboard