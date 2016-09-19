%% Design Allocation Objective function
function obj = objfun_DesAlloc(xC_val, xI, lookup, Model, acinputs)

if size(xC_val,1)==1
        xC_val=xC_val';
end

xC = lookup.xC_lb + xC_val.*(lookup.xC_ub - lookup.xC_lb);

num_con = length(xC);
% for ii = 1:num_con
%     if xC(ii)<lookup.xC_lb(ii) || xC(ii)>lookup.xC_ub(ii)
%         fprintf('\n%s\n','Design variable outside its bound. Cannot continue!')
% %         keyboard
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
price = network.price;
range = network.range;

existAC = existingAC();
num_ac = length(existAC.AC_name)+1;

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

profit= 0;
%% Size the aircraft for the given design range
ACdata = ReadInput2(MTOW,x_con,x_dis,acinputs);
[SizedAC] = AircraftSizing(MTOW,ACdata,Model,acinputs);

%% Fly missions to obtain the performance data of all the aircraft across
%% the routes
for kk = 1:num_ac
    for jj = 1:num_route
        range_kj = range(jj);
        pax_kj = pax((kk-1)*num_route + jj);
        route_info = [range_kj;pax_kj];
        x_kj = trip((kk-1)*num_route + jj);
        if kk == 1 % New aircraft
            [routedata] = FlyMissions(TOW(jj),ACdata,Model,SizedAC,route_info,acinputs);
            cost_kj = routedata.cost.Tot_cost_dep;
        else

        end
        profit = profit + ((price(jj)*pax_kj*x_kj) - (cost_kj*x_kj));
%         keyboard
    end
end
    
obj = -profit/1e3;
if isnan(obj)==1
%     fprintf('\n%s\n','Objective is NaN or Complex')
%     keyboard
end
% keyboard