%% This is the EGO Fraemwork2 for the simultaneous optimization
% The version 2 normalizes the design variables within [0,1]
clear;close all;clc
tst_EGOLF= tic;
% addpath(genpath(pwd));
warning('off','all'); 
Run_no = 1;
%% Test problems
% Activate line 10,30,94 to run test problems and deactivate line 13,29,93
[prob,num_des,FileName] = ProbDefinition(); %Define this file to select the test problems
%% Simultaneous optimization
% Activate line 13,15,29,93 to run test problems and deactivate line 10,15,30,94
% num_route=3; 

% FileName = 'DesAlloc';
iter = 1;
gen_lh = 1;
FileName_mat = [FileName,'.mat'];
%% 
ter_crit  = 0;
Tot_pt_prev=0;
y_opt=0;
if iter == 1 % If iter = 0 start from a existing result file
    %% Step 0: Initialize
    num_pt = 0;Tot_newpt_added=0;Tot_FunCount=0;
    ec2 = 0;ei_max=0;ei_tol_per=0.001;ei_tol_abs=0.001;
    new_sol=[];y_opt=Inf;comp=[];

%     [xI_lb,xI_ub,M] = Initialize(num_route);
    [xI_lb,xI_ub,M] = Initialize_Test(num_des,prob);
    
    num_xI = length(xI_lb);
    ModelInfo_obj.lb_org = xI_lb;
    ModelInfo_obj.ub_org = xI_ub;
    ModelInfo_obj.lb = zeros(num_xI,1);
    ModelInfo_obj.ub = ones(num_xI,1);
    ModelInfo_obj.X=[];ModelInfo_obj.X_org=[];ModelInfo_obj.xC=[];
    ModelInfo_obj.y=[];ModelInfo_obj.eflag=[];
    ModelInfo_obj.y_best = 0;  
    PerforData.BBtime=[];PerforData.GAtime=[];
    if M > 0
        for mm = 1:M
            ModelInfo_g{mm}.X = [];
            ModelInfo_g{mm}.y = [];
            ModelInfo_g{mm}.lb_org = xI_lb;
            ModelInfo_g{mm}.ub_org = xI_ub;
        end
    else
        ModelInfo_g=[];
    end
    FEA_obj = [];
    FEA_xopt =[];
%     plot2DIprob(xC_lb,xC_ub,ModelInfo_obj,ModelInfo_g,1);
    %% Step 1: Generate a set of initial integer points
    % x0I=[-5;3;10];
    max_pt_lim = 10*length(xI_lb); %Max computational budget
    % n = length(xI_lb); %length(x0I);
    n = floor(max_pt_lim/3); %Number of starting points
    %Override n
    n = length(xI_lb)+2;
    if num_xI<2
        x0I_01 = [0.0;0.33;.66;1.0];
    else
        if gen_lh == 1
%             x0I_01 = bestlh(n,num_xI,50,20);
%             cd Output_Test_nBB
%             save('InitialPoints.mat','x0I_01');
%             cd ..
%             x0I_01 = [0.0,0.25,0;1.0,0.75,0.75;0.51,1.0,0.25;0.25,0.51,1.0; 0.75,0.0,0.51];
              x0I_01 = [[1.0, 0.25, 0.75],
                            [0.0, 0.75, 0.0],
                            [0.75, 0.0, 0.25],
                            [0.75, 1.0, 0.49],
                            [0.25, 0.49, 1.0]];
        else
            cd Output_Test_nBB
            load InitialPoints
            cd ..
%             x0I_01=[];
        end
    end

    for ii= 1:n 
        x0I(ii,:) = [round(xI_lb + x0I_01(ii,:)'.*(xI_ub-xI_lb))]';
        x0I_hat(ii,:) = (x0I(ii,:) - xI_lb')./(xI_ub - xI_lb)';  %Normalized to [0,1] space
    end
end
x0I
x0I_01
% x0I_hat1_mean = mean(x0I_hat1);
% x0I_hat1_std = std(x0I_hat1);
% x0I_hat = (x0I_hat1 - x0I_hat1_mean)./x0I_hat1_std

% keyboard 
while ter_crit == 0
    if iter == 0
        load(FileName);
        xI_lb = ModelInfo_obj.lb_org;
        xI_ub = ModelInfo_obj.ub_org;
    else
        %% Step 2: Perform the simultaneous optimization given  x0I_alloc
        tstart_con = tic;
%         [xC_opt,obj,g,eflag,funCount] = ContinuousOptimization(x0I,M,num_route);
        [xC_opt,obj,g,eflag,funCount] = ContinuousOptimization_Test(x0I,M,num_des,prob)
        tend_con = toc(tstart_con);
        Tot_FunCount = Tot_FunCount + funCount;
        %% Step 3: Build the response surface
        for nonNAN = 1:size(x0I,1)
            if (isnan(obj(nonNAN))~=1) && (imag(obj(nonNAN)) == 0)  %Do not include NaN, imaganary solution
                % For Objective function (obj)
                num_pt = num_pt + 1;
                ModelInfo_obj.X = [ModelInfo_obj.X;x0I_hat(nonNAN,:)];
                ModelInfo_obj.X_org = [ModelInfo_obj.X_org;x0I(nonNAN,:)];
                ModelInfo_obj.xC = [ModelInfo_obj.xC;xC_opt(nonNAN,:)];
                ModelInfo_obj.y = [ModelInfo_obj.y;obj(nonNAN)];
                ModelInfo_obj.eflag = [ModelInfo_obj.eflag;eflag(nonNAN)];
                Tot_newpt_added = Tot_newpt_added + 1;
                %% For the Constraint functions that function of integer/discrete type design variables)
%                 flag_fea_xI = 1;
                if M>0
                    % g_i(x) <= 0, i = 1,2,...,M
                    % g -> n x M matrix
                    for mm = 1:M
                        if isnan(g(nonNAN,mm))~=1 && imag(g(nonNAN,mm))==0 %Do not include NaN and imaganary numbers
                            ModelInfo_g{mm}.X = [ModelInfo_g{mm}.X;x0I_hat(nonNAN,:)];
                            ModelInfo_g{mm}.y = [ModelInfo_g{mm}.y;g(nonNAN,mm)];
                        end
%                         % Check for feasibility
%                         if g(nonNAN,mm) > 1e-6
%                             flag_fea_xI = 0;
%                         end
                    end
                end
                if eflag(nonNAN)>=1
                    FEA_obj = [FEA_obj;obj(nonNAN)];
                    FEA_xopt = [FEA_xopt;[x0I(nonNAN,:),xC_opt(nonNAN,:)]];
                end
%                 keyboard
            end
        end
        ModelInfo_obj = build_surrogate(num_xI,ModelInfo_obj);
        fprintf('\nSurroagte building of the objective fun is complete...\n')
        if M>0
            for mm = 1:M
                ModelInfo_g{mm} = build_surrogate(num_xI, ModelInfo_g{mm});
                fprintf('\n%s%d%s\n','Surroagte building of g',mm,' constraint is complete...') 
            end
        end
        if length(FEA_obj)>=1
            [y_opt,min_ind] = min(FEA_obj);
            ModelInfo_obj.y_best = y_opt;
            x_opt = FEA_xopt(min_ind,:);
        end

        cd Output_Test_nBB %Provide path for the output files
        save(FileName_mat,'ModelInfo_obj','ModelInfo_g',...
        'FEA_obj','FEA_xopt','y_opt','Tot_newpt_added','new_sol','comp','PerforData');
        fprintf('\n%s%d%s\n','Saving data complete...')
        cd ..
    end
    plot2DIprob(0,15,ModelInfo_obj,ModelInfo_g,2);
    drawnow
%     if iter == 1
%         keyboard
%         [ModelInfo_obj, ModelInfo_g] = manual_override(ModelInfo_obj,ModelInfo_g,iter);
%     end
%% Step 4: Maximize the expected improvement function
    app_step4 = 3; % 1-MINLP BB., 2-GA,  3-Both
    fprintf('\n%s%d','Run number: ',Run_no)
    fprintf('\n%s%d\n','EGO-Iter: ',iter)
    %% Solve to maximize the expected improvement
    int_con =[1:length(xI_lb)]';
    if Tot_newpt_added~=Tot_pt_prev
        if app_step4 == 1 % Using nested branch and bound
            %Call nested BB function here
            tst_MINLPBB = tic;
%             keyboard
            [x_new, ei_min, eflag_MINLP_EI] = MINLP_BB(xI_lb,xI_ub,ModelInfo_obj,ModelInfo_g)
%             keyboasrd
            tend_MINLPBB = toc(tst_MINLPBB)
            PerforData.BBtime = [PerforData.BBtime;tend_MINLPBB];
        elseif app_step4 == 2 % Using GA
            tst_GA = tic;
            options=[];
%             options = gaoptimset('PopulationSize',500,'Generations',500); %,'UseParallel','always');
            [x_new,ei_min,eflag_MINLP_EI] = ga(@(x_new) objGA_EGOLFStep4(...
                x_new,ModelInfo_obj,ModelInfo_g,[],0),num_xI,[],[],[],[],...
            xI_lb,xI_ub,[],int_con,options)
            tend_GA = toc(tst_GA)
            PerforData.GAtime = [PerforData.GAtime;tend_GA];
        elseif app_step4 == 3 % Compare both the results
            tst_MINLPBB = tic;
            [x_new, ei_min, eflag_MINLP_EI] = MINLP_BB(xI_lb,xI_ub,ModelInfo_obj,ModelInfo_g)
            tend_MINLPBB = toc(tst_MINLPBB)
            PerforData.BBtime = [PerforData.BBtime;tend_MINLPBB];
%             keyboard
            tst_GA = tic;
            options=[];
%             options = gaoptimset('PopulationSize',10); %,'UseParallel','always');
            [x_new_GA,ei_min_GA,eflag_GA] = ga(@(x_new) objGA_EGOLFStep4(...
                x_new,ModelInfo_obj,ModelInfo_g,[],0),num_xI,[],[],[],[],...
            xI_lb,xI_ub,[],int_con,options)
            tend_GA = toc(tst_GA)
            PerforData.GAtime = [PerforData.GAtime;tend_GA];
            comp = [comp;[ei_min,ei_min_GA]];
            % Plot the expected improvement plot over the iterations
            figure(20)
            subplot(1,2,1);plot(iter,y_opt,'k*-','MarkerSize',10)
            hold on
            xlabel('Iterations')
            ylabel('Best Objective')
            subplot(1,2,2);plot(iter,ei_min,'r*-',iter,ei_min_GA,'g*-','MarkerSize',6)
            hold on
            xlabel('Iterations')
            ylabel('nBB vs GA')
            drawnow
        end
%         keyboard
%         %Overrride for iter 2
%         x_new = [2,4,4]
        if eflag_MINLP_EI>=1
            new_sol = [new_sol;ei_min];
            x0I = x_new;
            if size(x0I,1)>1
                x0I = x0I';
            end
            x0I_hat = (x0I - xI_lb')./(xI_ub' - xI_lb');
            ei_max = -ei_min; 
            Tot_pt_prev = Tot_newpt_added;
            %% Prevent Correlation matrix being close to singular.
            % No point to be added within a defined hypersphere around an
            % existing point
            rad = 0.5; %Pick any positive value between 0<rad<1
            cc = 0;
            for ii = 1:size(ModelInfo_obj.X_org,1)
                dist = sqrt(sum((ModelInfo_obj.X_org(ii,:)-x0I).^2));
                if dist<=rad
                    fprintf('\n%s','Point already exists!')
                    ec2=1;
                    break;
                end
            end
        else
            ec2=1;
        end
    else
        ec2 = 1;
    end
    if abs(y_opt)<=1e-6
        term = ei_tol_abs;
    else
        term = min(abs(ei_tol_per*y_opt),ei_tol_abs);
    end
    if ei_max <= term || ec2==1 || Tot_newpt_added>=max_pt_lim
        if ei_max<=term
            fprintf('\n%s\n','No further improvement expected! Terminating Algorithm....')
        elseif ec2==1
            fprintf('\n%s\n','No new point found/added that improves the surrogate! Terminating Algorithm....')
        elseif Tot_newpt_added>=max_pt_lim
            fprintf('\n%s\n','Maximum allowed sampling limit reached! Terminating Algorithm....')
        end
        ter_crit = 1;
    end
    iter = iter + 1;
end
Tot_ConMin = length(ModelInfo_obj.y);
fprintf('\n%s%f','The best objective: ',y_opt)
fprintf('\n%s%d', 'Total number of continuous minimization: ',Tot_ConMin)
fprintf('\n%s%d', 'Total number of objective function evalaution: ',Tot_FunCount)
fprintf('\n%s','Best Integer design: [')
for ii = 1:length(ModelInfo_obj.lb)
    if ii < length(ModelInfo_obj.lb)
        fprintf('%f%s',x_opt(ii),' , ')
    else
        fprintf('%f%s\n',x_opt(ii),']')
    end
end  

%% Time plot
% if app_step4 == 3
%     xx = 1:length(PerforData.GAtime);
%     figure(100)
%     plot(xx,PerforData.GAtime','rx-',xx,PerforData.BBtime,'bo-');
%     xlabel('Iterations (EGOLF)')
%     ylabel('time[s]')
%     legend('GA','MINLP BB')
%     title('Time Comparison - Solve MINLP (Max EI)')
% end
% rmpath(genpath(pwd));
tend_EGOLF = toc(tst_EGOLF);
