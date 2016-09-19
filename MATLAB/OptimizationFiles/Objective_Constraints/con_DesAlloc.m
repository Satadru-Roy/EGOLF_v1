function [g,ceq] = con_DesAlloc(xC_val,xI,lookup,Model,acinputs)

if size(xC_val,1)==1
        xC_val=xC_val';
end
xC = lookup.xC_lb + xC_val.*(lookup.xC_ub - lookup.xC_lb);

num_con = length(xC);
% for ii = 1:num_con
%     if xC(ii)<lookup.xC_lb(ii) || xC(ii)>lookup.xC_ub(ii)
%         fprintf('\n%s\n','Design variable outside its bound. Cannot continue!')
%         keyboard
%     end
% end
%%
% if size(xC,1)==1
%     xC=xC';
% end
% options2=optimset('Algorithm','active-set','Display','iter','MaxFunEvals',40);

num_route = lookup.num_route;
network = networkdata(num_route);
% dem = network.demand;
% price = network.price;
range = network.range;

existAC = existingAC();
num_ac = length(existAC.AC_name)+1;
AC_num = [lookup.AC_num_new;existAC.AC_num];
%% Integer/ Discrete data extract
% xI_num_alloc = num_route*num_ac;
trip = xI(2:end);
prof_lookup = xI(1);
top_lev = lookup.profile(prof_lookup,:);
design_range = top_lev(1);
seat_cap = top_lev(2);
Mc= top_lev(3);
%% Discrete design variable (See the discrete technology reference file to interpret the values)
% Engine settings (default) - For advance engine designs
BPR=0; % [0-10]
TET=0; % [0-500] Rankine
OPR=0; % [0-20]
FPR=0; % [0-0.1]   

% Discrete technologies (default)
x_dis=[ 2; %Laminar flow tech on wing [2-4]
        0; %Laminar flow tech on nascelle [0-1]
        0; %Laminar flow tech on tail [0-1]
        3; %Engine position [3-8]
        0; %Engine type [0-4]
        0; %Composite wing [0/1]
        0; %Composite fuselage [0/1]
        0; %Composite nascelle [0/1]
        0];%Composite tail [0/1]
    
x_dis = [x_dis;design_range;seat_cap;Mc];    
%% Continuous data extract
% xC_num_alloc = num_route*num_ac;
xC_num_des = lookup.num_des;
pax = xC(xC_num_des+1:end);
MTOW = xC(1);
TOW = xC(2:num_route+1);
wing_var = xC(num_route+2:xC_num_des-1);
eng_var = xC(xC_num_des);
x_con = [wing_var;BPR;TET;OPR;FPR;eng_var];

%% Size the aircraft for the given design range
ACdata = ReadInput2(MTOW,x_con,x_dis,acinputs);
[SizedAC,OutMission] = AircraftSizing(MTOW,ACdata,Model,acinputs);

scale_fac = 1e6;
con_tol = 100; %[lbs]

g=zeros(6+(2*num_route)+num_ac,1);
%Aircraft Sizing constraints
g(1,1) = (((SizedAC.W0 - MTOW)^2/con_tol) - 1)/scale_fac;
%Aircraft Performance constraints
g(2,1)=(OutMission.TOD/8500)-1; %Take-off distance
g(3,1)=(OutMission.LDD/7000)-1; %Landing distance
% g(4,1)=(SizedAC.LG_len*12/150)-1; %Landing gear length
% g(5,1)=(SizedAC.Wf_fuse/28800)-1; %Fuselage fuel capacity
% g(6,1)=(OutMission.y_min/OutMission.y_climb_SS)-1; %Second segment climb gradient requirement
g(4,1) = 1 - (OutMission.y_climb_SS/OutMission.y_min);
% g(7,1)=((100/60)/OutMission.Ps)-1; %Top of climb requirment 100ft/min
g(5,1) = 1 - (OutMission.Ps/(100/60));
g(6,1) = (OutMission.sin_yclimb_act/sind(15))-1;
%% Fly missions to obtain the performance data of all the aircraft across
%% the routes
cc=6; % If you change this, then also change line 74 above and Intialize.m (line 27) 
for kk = 1:num_ac
    con_val = 0;
    for jj = 1:num_route
        range_kj = range(jj);
        pax_kj = pax((kk-1)*num_route + jj);
        route_info = [range_kj;pax_kj];
        x_kj = trip((kk-1)*num_route + jj);
        if kk == 1 %New aircraft
            [routedata] = FlyMissions(TOW(jj),ACdata,Model,SizedAC,route_info,acinputs);
            g(cc+1,1) = (((routedata.W0 - TOW(jj))^2/con_tol) - 1)/scale_fac;
            g(cc+2,1) = (routedata.W0/SizedAC.W0)-1;
            cc=cc+2;
%             keyboard
        else
            
        end
        con_val = con_val + x_kj*(routedata.BH*(1+ACdata.RDTE.MH_FH(kk)) + 1);  
    end
    g(cc+1,1) = (con_val/(12*AC_num(kk)))-1;
    cc = cc+1;
end
% keyboard
% for ii = 1:length(g)
%     if isnan(g(ii))==1 || imag(g(ii))~=0
%         fprintf('\n%s\n','Constraint is NaN or Complex!')
% %         keyboard
%     end
% end
ceq = 0;