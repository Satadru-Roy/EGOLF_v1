%% This is the Branch and Bound algorithm that maximizes the conEI function
% This combines the outer and inner BB of the previous version into a
% hybrid BB
% Problem formulation
% Minimize f
% Subject to: xI_L <= xI <= xI_U
% xI belongs to original integer space
function [xopt,fopt,eflag_MINLPBB] = MINLP_BB(xI_lb, xI_ub, ModelInfo_obj,ModelInfo_g)
%% Inputs
% xI_lb - Lower bound integer type design variables (orginal design space)
% xI_ub - Upper bound integer type design variables (orginal design space)
% ModelInfo_obj - Surrogate model data (objective) 
% ModelInfo_g - Surrogate model data (constraints)
%% Outputs
% x_best - Optimal integer solution
% f_best -  Optimal objective value
% eflag_MINLPBB -  Exit flag MINLP BB

%% Step 1: Initialize
iter = 1; 
tol_ec = 1e-6;
term = 0;
LBD = -Inf;
node_num = 0;
canX=[];canF=[];
xL_iter = xI_lb;xU_iter = xI_ub;
Aset = []; % Active Set fields: Aset = [NodeNumber, lb, ub, LBD, UBD]
UBD_iter = [];
par_node = 0;
heur_search = 0;
LBD_prev=-Inf;
fC = 0;
fprintf('%s\n','Starting the MINLP BB...')
num_des = length(xI_lb);
if heur_search == 1
    [xopt, fopt, eflag_heur] = ga(@(xopt)...
        combined_obj(xopt,ModelInfo_obj,ModelInfo_g,[],0),num_des,...
        [],[],[],[],xL_iter,xU_iter);  
    UBD=fopt; 
    xopt = xopt';
    fC = fC + outputheur.funcCount;
    if eflag_heur>=1
        eflag_MINLPBB=1;
        fprintf('\n%s\n','Good starting point found!...')
    end
else
%     UBD = Inf;
    %Randomly generates an integer point
    xopt = round(xI_lb + rand(length(xI_lb),1).*(xI_ub-xI_lb));
%     xopt =[2;2;2]; 
    UBD = combined_obj(xopt,ModelInfo_obj,ModelInfo_g,[],0);
    fC = fC+1;
    fopt = UBD;
    eflag_MINLPBB = 1;
end
% keyboard
fprintf('%s','====================================================================================')
fprintf('\n%19s%12s%14s%21s','Global','Parent','Child 1','Child 2')
fprintf('\n%s%8s%10s%8s%9s%11s%10s%11s%11s\n','Iter','LBD','UBD','Node',...
    'Node1','LBD1','Node2','LBD2','FLocal')
fprintf('%s\n','====================================================================================')
% keyboard
while term == 0
%     fprintf('\n%s%d','Iteration: ',iter)
%     con_fac=concave_factor(xL_iter, xU_iter);
    con_fac=[];
    con_flag = 0;
    %% Step 2: Obtain a local solution of f (-conEI) and update global
    app_step2 = 1; 
    loc_search = 0;
%         con_fac = concave_factor(xL_iter, xU_iter);
    if app_step2 == 1 %Approach 1: Use a gradient-based approach 
%         options = optimset();
        xloc_iter = round(0.49.*(xU_iter+xL_iter));
        floc_iter = combined_obj(xloc_iter,ModelInfo_obj,ModelInfo_g,con_fac,con_flag);
        efloc_iter = 1;
        fC=fC+1;
%         keyboard
        if loc_search == 1 %Perform the local search
            if abs(floc_iter) > 1e-6 %Do local search only if the starting region is NOT flat
                options = optimset('Display','none');
                xC_iter = xloc_iter;
                [xloc_iter, ~, efloc_iter,outputloc] = fmincon(@(xloc_iter)...
                    combined_obj(xloc_iter,ModelInfo_obj,ModelInfo_g,con_fac,con_flag),xC_iter,...
                    [],[],[],[] ,xL_iter,xU_iter,[],options);
                %Rounding to nearest integer solution
                xloc_iter = round(xloc_iter);
                floc_iter = combined_obj(xloc_iter,ModelInfo_obj,ModelInfo_g,con_fac,con_flag); 
                fC = fC + outputloc.funcCount+1;
            end
        end
%             keyboard
    elseif app_step2 == 2 %Approach 2: Use GA
%         options = gaoptimset('PopulationSize',50,'Generations',200); %,'UseParallel','always');
        [xloc_iter, ~, efloc_iter] = ga(@(xloc_iter)...
            combined_obj(xloc_iter,ModelInfo_obj,ModelInfo_g,con_fac,1),num_des,...
            [],[],[],[],xL_iter,xU_iter);  
            %Rounding to nearest integer solution
            xloc_iter = round(xloc_iter);
            floc_iter = combined_obj(xloc_iter,ModelInfo_obj,ModelInfo_g,con_fac,con_flag); 
    end
  
    %% Step 3: Partition current rectangle
    child_info = zeros(2,3);
    dis_flag = [' ';' '];
    for ii = 1:2
        lb = xL_iter;
        ub = xU_iter;
        [~,l_iter] = max(xU_iter - xL_iter); %Branch at variable with the largest edge
        if xloc_iter(l_iter)<ub(l_iter)
            delta = 0.5; % 0<delta<1 
        else
            delta = -0.5; % -1<delta<0
        end
        if ii == 1
            ub(l_iter) = floor(xloc_iter(l_iter)+delta); 
        elseif ii == 2
            lb(l_iter) = ceil(xloc_iter(l_iter)+delta);
        end
        if norm(ub - lb)>1e-6
            %% Step 4: Obtain an LBD of f ('-conEI') in both the newly created rectangles
            S4_fail = 0;
            % Expected Improvement function
            [x_comL, x_comU, Ain_hat, bin_hat] = gen_coeff_bound2(lb,ub,ModelInfo_obj);
%             keyboard
            [sU,eflag_sU] = maximize_S2(x_comL, x_comU, Ain_hat, bin_hat, ModelInfo_obj);
%             keyboard
            if eflag_sU>=1
                [yL,eflag_yL] = minimize_y2(x_comL, x_comU, Ain_hat, bin_hat, ModelInfo_obj);
%                 keyboard
                if eflag_yL>=1
                    NegEI = calc_conEI_norm2([],ModelInfo_obj,sU,yL);
                    % Expected Violation of each constraints
                    M = length(ModelInfo_g);
                    EV = zeros(M,1);
                    if M>0
                        for mm = 1:M
                            [x_comL, x_comU, Ain_hat, bin_hat] = gen_coeff_bound2(lb,ub,ModelInfo_g{mm});
                            [sU_g,eflag_sU_g] = maximize_S2(x_comL, x_comU, Ain_hat, bin_hat, ModelInfo_g{mm});
                            sL_g = -sU_g;
    %                             keyboard
                            if eflag_sU_g >=1
                                [yL_g,eflag_yL_g] = minimize_y2(x_comL, x_comU, Ain_hat, bin_hat, ModelInfo_g{mm});
    %                                 keyboard
                                if eflag_yL_g>=1
                                    EV(mm) = calc_conEV_norm2([],ModelInfo_g{mm},sL_g,yL_g);
                                else
                                    S4_fail = 1;
                                    break;
                                end
                            else
                                S4_fail=1;
                                break;
                            end
                        end
                    end
                else
    %                 fprintf('\n%s\n','Warning! Cannot solve Min y_hat problem!')
                    S4_fail = 1;
                end
            else
    %             fprintf('\n%s\n','Warning! Cannot solve Max S problem!')
                S4_fail = 1;
            end
            if S4_fail == 1 %If convex approximation is not possible, 
                % assign the parent's LBD value, if local search was successful
    %                 LBD_NegConEI = Inf;
                  if efloc_iter>=1
                      LBD_NegConEI = LBD_prev;
                  else
                      LBD_NegConEI = Inf;
                  end
                  dis_flag(ii)='F';
            else
                LBD_NegConEI = NegEI/(1+sum(EV)); 
            end  
    %             keyboard
            %% Step 5: Store any rectangle inside the active set that has LBD lower than the global UBD
            if (LBD_NegConEI) < UBD
                node_num = node_num + 1;
                newNode = [node_num, lb', ub', LBD_NegConEI, floc_iter];
                Aset = [Aset;newNode];
                child_info(ii,:) = [node_num,LBD_NegConEI,floc_iter];
                dis_flag(ii)=' ';
            else
                child_info(ii,:) = [par_node,LBD_NegConEI,floc_iter];
                dis_flag(ii)='X'; %Not added to active set
            end
        else
            xloc_iter = lb;
            floc_iter = combined_obj(xloc_iter,ModelInfo_obj,ModelInfo_g,con_fac,con_flag);  
            fC=fC+1;
            child_info(ii,:) = [par_node,Inf,floc_iter]; %No child created
            dis_flag(ii)='x'; %Flag for no child created
        end
        %% Update the active set
        canX = [canX;xloc_iter];
        canF = [canF;floc_iter];
        if floc_iter < UBD % Better integer solution found
            UBD = floc_iter;
            fopt = UBD;
            xopt = xloc_iter;  
            eflag_MINLPBB = 1;
            % Update the Active Set
            if size(Aset,1)>=1
                del_flag=[];
                for aaa = 1:size(Aset,1)
                    if Aset(aaa,end-1)>=UBD
                        del_flag=[del_flag;aaa];
                    end
                end
                Aset(del_flag,:)=[];
            end
        end
    end
%       keyboard
    %% Display the output in a tabular format
    if mod(iter,25)==0
        fprintf('%s','====================================================================================')
        fprintf('\n%19s%12s%14s%21s','Global','Parent','Child 1','Child 2')
        fprintf('\n%s%8s%10s%8s%9s%11s%10s%11s%11s\n','Iter','LBD','UBD','Node',...
            'Node1','LBD1','Node2','LBD2','FLocal')
        fprintf('%s\n','====================================================================================')
    end
    fprintf('%3d%10.2f%10.2f%6d%8d%1s%13.2f%7d%1s%13.2f%9.2f\n',iter,LBD,UBD,par_node,child_info(1,1),dis_flag(1),child_info(1,2),...
        child_info(2,1),dis_flag(2),child_info(2,2),child_info(2,3));
%     keyboard
    if isempty(Aset)==0
        %% Step 6: Update LBD and select the current rectangle
        % a. Set LBD as the lowest sol in the active set.
        [LBD,ind_LBD] = min(Aset(:,end-1));
%         if LBD<LBD_prev
%             fprintf('\n%s\n','Error: Previous LBD is greater than the present LBD. LBD should be non-decreasing!!')
%             keyboard
%         end
        LBD_prev = LBD;
        % b. Select this as the current rectangle
        xL_iter = Aset(ind_LBD,2:num_des+1)';
        xU_iter = Aset(ind_LBD,num_des+2:end-2)';
        par_node = Aset(ind_LBD,1);
        % c. Delete this rectangle from the set of active node
        Aset(ind_LBD,:)=[];
        iter = iter + 1;
        %% Step 7: Check for convergence
        diff = (UBD - LBD);
        if diff <= tol_ec % 1. Meets the tolerance gap criteria
            term = 1;
            fprintf('%s\n','====================================================================================')
            fprintf('%s\n','Terminating! Absolute difference between the upper and lower bound is below the tolerance limit.')
        end
    else
        term = 1;
        fprintf('%s\n','====================================================================================')
        fprintf('%s\n','Terminating! No new node to explore.')
    end
end
fprintf('%s%d\n','Function Count: ',fC)

