% This function plots the original function/Expected improvement function
function [] = plot2DIprob(xC_lb,xC_ub,ModelInfo_obj,ModelInfo_g,plot_flag)

xI_lb = ModelInfo_obj.lb_org;
xI_ub = ModelInfo_obj.ub_org;
M = length(ModelInfo_g);
kI = length(xI_lb);
kC = length(xC_lb);
num_des = kI+kC;
if num_des >= 1 % Any Dimension Problem
    if kI == 1 && kC == 0 % 1 integer, no continuous problem
    elseif kI==1 && kC == 1 %1 integer, 1 continuous problem
        if plot_flag == 1 % Plot the original function
            figure(10)
            xI = xI_lb:.1:xI_ub;
            xC = xC_lb:.1:xC_ub;
            for ii = 1:length(xI)
                for jj = 1:length(xC)
                    f(jj,ii) = objfun_Test(xC(jj),xI(ii),1,1);
                    infea_flag=0;
                    if M>0
                        g = confun_Test(xC(jj), xI(ii), 1,1);
                        for mm = 1:M  
                            if g(mm)>1e-6
                                infea_flag = 1;
                                break;
                            end
                        end
                    end
                    if infea_flag == 0
                        fcon(jj,ii) = objfun_Test(xC(jj),xI(ii),1,1);
                    else
                        fcon(jj,ii) = nan;
                    end          
                end
            end
            subplot(1,2,1);surf(xI,xC,f);
            xlabel('xI');ylabel('xC');zlabel('f(x)-Unconstrained');
            hold on
            subplot(1,2,2);surf(xI,xC,fcon);
            xlabel('xI');ylabel('xC');zlabel('f(x)-Constrained');
            drawnow;
        elseif plot_flag == 2 % Plot the expected improvement function/s
            figure(11)
            xI = xI_lb:.1:xI_ub;
            for ii = 1:length(xI)
                xval = (xI(ii) - xI_lb)/(xI_ub - xI_lb);
                ei(ii) = calc_conEI_norm2(xval,ModelInfo_obj);
                conei(ii) = combined_obj_norm(xval,ModelInfo_obj,ModelInfo_g,[],0);
                ev = 0;
                for mm = 1:M
                    ev = ev+calc_conEV_norm2(xval,ModelInfo_g{mm});
                end
                Fx(ii) = 1 + ev;
            end
%             subplot(2,2,1);plot(xI,ei,'k-');
%             plot(xI,ei,'k-');
%             xlabel('xI');ylabel('Exp Imp')
%             subplot(2,2,2);plot(xI,conei,'k-');
%             xlabel('xI');ylabel('Con Exp Imp')
%             subplot(2,2,3);plot(xI,Fx)
%             xlabel('xI');ylabel('Fx(Total Exp Violation)')
%             keyboard
        end
    elseif kI==1 && kC>=2 
    elseif kI == 2 && kC == 0 % 2 integer, 0 continuous problem
    elseif kI == 2 && kC>=1 % 2integer and Any continuous design problem
        if plot_flag == 1
            fprintf('\n%s\n','Cannot plot the original function! Problem dimensions greater than 2')
        elseif plot_flag == 2 % Plot the expected improvement function/s
            figure(11)
            xI1 = xI_lb(1):.1:xI_ub(2);
            xI2 = xI_lb(2):.1:xI_ub(2);
            for ii = 1:length(xI1)
                for jj =1:length(xI2)
                    x = [xI1(ii);xI2(jj)];
                    xval = (x - xI_lb)./(xI_ub - xI_lb);
                    ei(jj,ii) = calc_conEI_norm2(xval,ModelInfo_obj);
                    conei(jj,ii) = combined_obj_norm(xval,ModelInfo_obj,ModelInfo_g,[],0);
                    ev= 0;
                    for mm = 1:M
                        ev = ev + calc_conEV_norm2(xval,ModelInfo_g{mm});
                    end
                    Fx(jj,ii) = 1+ev;
                end
            end
            subplot(2,2,1);surf(xI1,xI2,ei);
            xlabel('x_I1');ylabel('x_I2');zlabel('Exp Imp')
            subplot(2,2,2);surf(xI1,xI2,conei);
            xlabel('x_I1');ylabel('x_I2');zlabel('Con Exp Imp')
            subplot(2,2,3);surf(xI1,xI2,Fx)
            xlabel('x_I1');ylabel('x_I2');zlabel('Fx(Total Exp Violation)')
        end
    else
        fprintf('\n%s\n','Cannot plot. Number of integer variables > 2')
    end
else
    fprintf('\n%s\n','Cannot be a zero dim problem');
end
