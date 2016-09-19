clc
clear all
close all


%% Langermann function
x=0:0.05:10;
y=0:0.05:10;

c =[
     3
    10
     2
     6
     1
    10
     1
     3
     5
     7
    ];

m=size(c,1);

A =[
     8     4
     7     6
     6    10
     1     7
     4     4
     6     4
     6     5
     9     5
     5     9
     1    10
     ];

%% Calculate function values
f=zeros(length(y),length(x));
for j = 1:length(x)
    for i = 1:length(y)
        xx=[x(j),y(i)];
%         [f(i,j),m,c,A] = langer(xx);
        [f(i,j),m,c,A] = langer(xx,m,c,A);
    end
end

figure
surf(x,y,f)
% surfc(x,y,f)

Markersize=48;
Fontsize=25;        % 25 for 1 column, 35 for 2 columns

set(gca,'XLim',[0 10],'XTick',0:2:10,'FontSize',Fontsize)
set(gca,'YLim',[0 10],'YTick',0:2:10,'FontSize',Fontsize)
xlabel('x_1','FontSize',Fontsize)
ylabel('x_2','FontSize',Fontsize)
zlabel('f(x_1,x_2)','FontSize',Fontsize)
set(gcf,'PaperPositionMode','auto','Position',[0 0 1600 1000])
title('Langermann Function')


%% Plot the global minimum point
% hold on
% x1 = 6.5:0.01:7.5;
% x2 = 4:0.01:5;
% y = -13*ones(size(x1,2),size(x2,2));
% surf(x1,x2,y)
% plot3(6.945443187712963,4.275153689318769,-13.550738262018468,'r.','Markersize',Markersize)