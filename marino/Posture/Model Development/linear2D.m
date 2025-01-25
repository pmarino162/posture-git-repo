clear;
clc;
%close all

%% Define System 
% a = 1;
% b = 0;
% tau = 0;
% d = tau - a;
% detA = -1;
% delta = tau^2 - 4*detA;
% c = (a*d - detA)/b;
% if b == 0
%     detA = a*d;
%     delta = tau^2 - 4*detA;
%     c = 0;
% end
% 
% A = [a b; c d];
% I = [0;0];
% 
A = [0 -1; -1 0];
I = [3;3];
detA = det(A);
tau = trace(A);
delta = tau^2 - 4*detA;
%% Calculate Flow Field
x1Lim = [-10,10];
x2Lim = [-10,10];
dx = 1;
x1Pts = x1Lim(1,1):dx:x1Lim(1,2);
x2Pts = x2Lim(1,1):dx:x2Lim(1,2);
numPts = size(x1Pts,2)*size(x2Pts,2);
xdot = zeros(numPts,4);
FPs = nan(numPts,2);
i = 0;
for x1 = x1Pts
    for x2 = x2Pts
        i = i + 1;
        x = [x1,x2]';
        %Linear
        xdot(i,1:2) = x';
        xdot(i,3:4) = (A*x)' + I';
        if sum(abs(xdot(i,3:4)))==0
            FPs(i,1:2) = x';
        end
    end
end
FPs = FPs(~isnan(FPs(:,1)),:);

%% Plot Flow Field
f = figure
quiver(xdot(:,1),xdot(:,2),xdot(:,3),xdot(:,4))
hold on
plot(FPs(:,1),FPs(:,2),'r.','MarkerSize',10)
axis square
axis equal
xlabel('x_1')
ylabel('x_2')
title(['A = ','[',num2str(A(1,:)),';',num2str(A(2,:)),']',newline,'\tau = ',num2str(tau),...
    ' detA = ',num2str(detA),' \Delta = ',num2str(delta), ...
    newline, ' I = [',num2str(I(1,1)),';',num2str(I(2,1)),']']);
titleStr = ['A=','[',num2str(A(1,:)),';',num2str(A(2,:)),']','I=[',num2str(I(1,1)),';',num2str(I(2,1)),']'];
dirStr = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Figures\Model Development\Linear 2D\';
 saveas(gcf,[dirStr,titleStr,'.jpg'])  
%set(t, 'horizontalAlignment', 'left')
% set(t, 'units', 'normalized')
% h1 = get(t, 'position')
% set(t, 'position', [1 h1(2) h1(3)])
% for i=1:4
%     s=subplot(1,4,i)
%     system = Results(i).System;
%     flowfield = Results(i).FlowField;
%     quiver(flowfield(:,1),flowfield(:,2),flowfield(:,3),flowfield(:,4))
% 
%     axis square
%     axis equal
%     s.XLim = x1Lim;
%     s.YLim = x2Lim;
%     title(system)
% end
% 
% f.Position = [0 0 1700 500];
