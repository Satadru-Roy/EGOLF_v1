% This function defines the problem and plot the objective function
function [ModelInfo] = ProblemFile(k,prob,run)
figure(10);
if k==1
    %% 1D Problem 
    X = [0;0.5;1];
    lb=0;ub=1;
    y = ((6.*X - 2).^2).*sin(12.*X -4);
    
    ModelInfo.X = X;ModelInfo.y=y;
    ModelInfo.lb = lb;ModelInfo.ub = ub;
    ModelInfo = build_surrogate(k,ModelInfo);
%     alpha = calc_alpha(ModelInfo);
%     ModelInfo.alpha = alpha;
    ModelInfo.y_best = min(y);
    %% Plot the EI
    xx = 0:.05:1;
    for ii = 1:length(xx)
        fval(ii) = calc_conEI_norm(xx(ii),ModelInfo,[],0);
    end
    plot(xx,fval)
    ylabel('Expected Improvement')
    xlabel('x')
elseif k==2
    %% 2D Problem
    if prob ==1
        %% Branin function
        lb = [-5;0];
        ub = [10;15];
        if run == 1
            n = 15*k;
            X = bestlh(n,k,50,20);
            for ii = 1:n
                xval = lb + (X(ii,:)'.*(ub-lb));
                y(ii,1) = objfun_Test(xval,[],2,prob);
            end
            ModelInfo.X = X;ModelInfo.y=y;
            ModelInfo.lb = zeros(k,1);
            ModelInfo.ub = ones(k,1);
            ModelInfo.lb_org = lb;
            ModelInfo.ub_org = ub;
            ModelInfo = build_surrogate(k,ModelInfo);
%             alpha = calc_alpha(ModelInfo);
%             ModelInfo.alpha = alpha;
            ModelInfo.y_best = min(y);
            cd DataTestFiles
            save('Data_2DBranin_30pts.mat','ModelInfo')
            cd ..
        else
            load Data_2DBranin_30pts
        end
    elseif prob == 2
        %% Griewank function
        lb = [0;0];
        ub = [5;5];
        if run == 1
            n = 5*k;
            X = bestlh(n,k,50,20);
            for ii = 1:n
                xval = lb + (X(ii,:)'.*(ub-lb));
                y(ii,1) = objfun_Test(xval,[],2,prob);
            end
            ModelInfo.X = X;ModelInfo.y=y;
            ModelInfo.lb = zeros(k,1);
            ModelInfo.ub = ones(k,1);
            ModelInfo.lb_org = lb;
            ModelInfo.ub_org = ub;
            ModelInfo = build_surrogate(k,ModelInfo);
%             alpha = calc_alpha(ModelInfo);
%             ModelInfo.alpha = alpha;
            ModelInfo.y_best = min(y);
            cd DataTestFiles
            save('Data_2DGriewank_10pts.mat','ModelInfo')
            cd ..
        else
            load Data_2DGriewank_10pts
        end
    end 

    % Plot the surface
    xx1 = 0:0.01:1;
    xx2 = 0:0.01:1;
    con_fac = concave_factor(ModelInfo.lb,ModelInfo.ub);
    for ii = 1:length(xx1)
        for jj = 1:length(xx2)
%                 xval = lb + ([xx1(ii);xx2(jj)].*(ub-lb));
            xval = [xx1(ii);xx2(jj)];
            fval1(jj,ii) = calc_conEI_norm(xval,ModelInfo,con_fac,1);
            fval2(jj,ii) = calc_conEI_norm(xval,ModelInfo,[],0);
        end
    end

    x1 = lb(1) + xx1.*(ub(1)-lb(1));
    x2 = lb(2) + xx2.*(ub(2)-lb(2));
    subplot(1,2,1);surf(x1,x2,fval2)
    title('Original')
    zlabel('Expected Improvement')
    xlabel('x1')
    ylabel('x2')
    subplot(1,2,2);surf(xx1,xx2,fval1)
    title('Concaved (Normalized)')
    zlabel('Expected Improvement')
    xlabel('x1')
    ylabel('x2')
elseif k>2
    if prob == 2
        %% Griewank function
        lb = -5.*ones(k,1);
        ub = 5.*ones(k,1);
        if run == 1
            n = 3*k;
            X = bestlh(n,k,50,20);
            for ii = 1:n
                xval = lb + (X(ii,:)'.*(ub-lb));
                y(ii,1) = objfun_Test(xval,[],2,prob);
            end
            ModelInfo.X = X;ModelInfo.y=y;
            ModelInfo.lb = zeros(k,1);
            ModelInfo.ub = ones(k,1);
            ModelInfo.lb_org = lb;
            ModelInfo.ub_org = ub;
            ModelInfo = build_surrogate(k,ModelInfo);
%             alpha = calc_alpha(ModelInfo);
%             ModelInfo.alpha = alpha;
            ModelInfo.y_best = min(y);
            cd DataTestFiles
            save('Data_4DGriewank_12pts.mat','ModelInfo')
            cd ..
        else
            load Data_4DGriewank_12pts
        end
    end 
end
