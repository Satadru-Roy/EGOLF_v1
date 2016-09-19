%% Constraint function
%Generates the constraint coefficient for the problem
function [Ain_hat,bin_hat,Ain,bin] = lincon_DesAlloc(trip,lookup)
num_route = lookup.num_route;
network = networkdata(num_route);
dem = network.demand;
existAC = existingAC();
num_ac = length(existAC.AC_name)+1;
num_des = lookup.num_des;
%% Linear Inequality constraints (Ain*x <= bin)
% Upper demand constraint
Ain = zeros(num_route,length(lookup.xC_lb));
bin = dem;
for jj = 1:num_route
  for kk = 1:num_ac
      col1 = num_des + ((kk-1)*num_route + jj);
      col2 = (kk-1)*num_route + jj;
      Ain(jj,col1) = trip(col2);
  end
end
arr = [(lookup.xC_ub-lookup.xC_lb)';(lookup.xC_ub-lookup.xC_lb)';(lookup.xC_ub-lookup.xC_lb)'];
Ain_hat = Ain.*arr;
bin_hat = bin - Ain*lookup.xC_lb;
