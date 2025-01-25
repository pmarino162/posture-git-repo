clear; clc; clf; close all;

%% Setup saveFig   
    saveFig = false;
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Paper\20220907\Figure 4';
    set(0, 'DefaultFigureRenderer', 'painters');
    
%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\orli.mat')
    pcmap = orli;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat')
    tcmap = customRainbow;

%%  Create base trajectory
    t = 0:100:100;
    numPts = length(t);
    baseTraj = zeros(numPts,3);
    baseTraj(:,1) = t./100;
    baseTraj(:,2) = 0;
    baseTraj(:,3) = 0;
        
%% Set up posture and target lists
    postureList = 1:5;
    targetList = 1:4;
    
%% Create trajStruct
    trajStruct = struct('posture',[],'target',[],'traj',[]);
    trajStructInd = 1;
    for posture = postureList
        for target = targetList
            %Get posture offset
            postureOffset = [0,0,(posture-1)./2];
            %Get rotation matrix
            rotMat = rotZ(90*(target-1));
            %Apply stretch
            stretchFac = 1;
            stretchMat = stretchFac.*eye(3);
            %Fill trajStruct
            trajStruct(trajStructInd).posture = posture;
            trajStruct(trajStructInd).target = target;
            trajStruct(trajStructInd).traj = baseTraj*stretchMat*rotMat + postureOffset;
            trajStructInd = trajStructInd + 1;
        end
    end
       
%% Plot Trajectories
fs = 14; %Font size
ds = 25; %Dot Size
as = 10; %Arrow Size
lw = 3.5; %LineWidth
az = 40; el = 15;

%All trajectories 
figure; hold on;
xDim = 1; yDim = 2; zDim = 3;
for posture = postureList
    for target = targetList
        traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).traj;
        plot3(traj(:,xDim),traj(:,yDim),traj(:,zDim),'Color',pcmap(posture,:),'LineWidth',lw)
        plot3(traj(1,xDim),traj(1,yDim),traj(1,zDim),'.','MarkerEdgeColor',pcmap(posture,:),'MarkerSize',ds);
        plot3(traj(end,xDim),traj(end,yDim),traj(end,zDim),'v','MarkerEdgeColor',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:),'MarkerSize',as);
    end
end
grid on
xlabel('Neural Dim 1'); ylabel('Neural Dim 2'); zlabel('Neural Dim 3')
xticklabels({}); yticklabels({}); zticklabels({});
view(az,el)
axis equal
set(gca,'fontname','arial')
set(gca,'fontsize',fs)
if saveFig
    saveas(gcf,fullfile(saveDir,'schemAllTraj.svg'));
end
%Get axis lims from previous plot
ax = gca; xlims = ax.XLim; ylims = ax.YLim; zlims = ax.ZLim;

%All trajectories except p5t1
figure; hold on;
xDim = 1; yDim = 2; zDim = 3;
for posture = postureList
    for target = targetList
        if posture == 5 && target == 1
        else
            traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).traj;
            plot3(traj(:,xDim),traj(:,yDim),traj(:,zDim),'Color',pcmap(posture,:),'LineWidth',lw)
            plot3(traj(1,xDim),traj(1,yDim),traj(1,zDim),'.','MarkerEdgeColor',pcmap(posture,:),'MarkerSize',ds);
            plot3(traj(end,xDim),traj(end,yDim),traj(end,zDim),'v','MarkerEdgeColor',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:),'MarkerSize',as);
        end
    end
end
grid on
xlabel('Neural Dim 1'); ylabel('Neural Dim 2'); zlabel('Neural Dim 3')
xticklabels({}); yticklabels({}); zticklabels({});
view(az,el)
xlim(xlims); ylim(ylims); zlim(zlims);
axis equal
set(gca,'fontname','arial')
set(gca,'fontsize',fs)
if saveFig
    saveas(gcf,fullfile(saveDir,'schemAllTrajExceptP5T1.svg'));
end

%Only P1 t1; purple
figure; hold on;
xDim = 1; yDim = 2; zDim = 3;
for posture = 1
    for target = 1
        traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).traj;
        plot3(traj(:,xDim),traj(:,yDim),traj(:,zDim),'Color',tcmap(8,:),'LineWidth',lw)
        plot3(traj(1,xDim),traj(1,yDim),traj(1,zDim),'.','MarkerEdgeColor',tcmap(8,:),'MarkerSize',ds);
        plot3(traj(end,xDim),traj(end,yDim),traj(end,zDim),'v','MarkerEdgeColor',tcmap(8,:),'MarkerFaceColor',tcmap(8,:),'MarkerSize',as);
    end
end
grid on
xlabel('Neural Dim 1'); ylabel('Neural Dim 2'); zlabel('Neural Dim 3')
xticklabels({}); yticklabels({}); zticklabels({});
view(az,el)
axis equal
xlim(xlims); ylim(ylims); zlim(zlims);
set(gca,'fontname','arial')
set(gca,'fontsize',fs)
if saveFig
    saveas(gcf,fullfile(saveDir,'schemOnlyP1T1.svg'));
end

%Only IC
figure; hold on;
xDim = 1; yDim = 2; zDim = 3;
for posture = 5
    for target = 1
       traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).traj;
       plot3(traj(1,xDim),traj(1,yDim),traj(1,zDim),'.','MarkerEdgeColor',tcmap(8,:),'MarkerSize',ds);
    end
end
grid on
xlabel('Neural Dim 1'); ylabel('Neural Dim 2'); zlabel('Neural Dim 3')
xticklabels({}); yticklabels({}); zticklabels({});
view(az,el)
axis equal
xlim(xlims); ylim(ylims); zlim(zlims);
set(gca,'fontname','arial')
set(gca,'fontsize',fs)
if saveFig
    saveas(gcf,fullfile(saveDir,'schemOnlyP5IC.svg'));
end

%Only P1 t1; and IC 
figure; hold on;
xDim = 1; yDim = 2; zDim = 3;
for posture = 1
    for target = 1
        traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).traj;
        plot3(traj(:,xDim),traj(:,yDim),traj(:,zDim),'Color',tcmap(8,:),'LineWidth',lw)
        plot3(traj(1,xDim),traj(1,yDim),traj(1,zDim),'.','MarkerEdgeColor',tcmap(8,:),'MarkerSize',ds);
        plot3(traj(end,xDim),traj(end,yDim),traj(end,zDim),'v','MarkerEdgeColor',tcmap(8,:),'MarkerFaceColor',tcmap(8,:),'MarkerSize',as);
    end
end
for posture = 5
    for target = 1
       traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).traj;
       plot3(traj(1,xDim),traj(1,yDim),traj(1,zDim),'.','MarkerEdgeColor',tcmap(8,:),'MarkerSize',ds);
    end
end
grid on
xlabel('Neural Dim 1'); ylabel('Neural Dim 2'); zlabel('Neural Dim 3')
xticklabels({}); yticklabels({}); zticklabels({});
view(az,el)
axis equal
xlim(xlims); ylim(ylims); zlim(zlims);
set(gca,'fontname','arial')
set(gca,'fontsize',fs)
if saveFig
    saveas(gcf,fullfile(saveDir,'schemP1T1_andIC.svg'));
end

%Final predicted trajectory
figure; hold on;
xDim = 1; yDim = 2; zDim = 3;
for posture = 5
    for target = 1
        traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).traj;
        plot3(traj(:,xDim),traj(:,yDim),traj(:,zDim),'Color',tcmap(8,:),'LineWidth',lw)
        plot3(traj(1,xDim),traj(1,yDim),traj(1,zDim),'.','MarkerEdgeColor',tcmap(8,:),'MarkerSize',ds);
        plot3(traj(end,xDim),traj(end,yDim),traj(end,zDim),'v','MarkerEdgeColor',tcmap(8,:),'MarkerFaceColor',tcmap(8,:),'MarkerSize',as);
    end
end
grid on
xlabel('Neural Dim 1'); ylabel('Neural Dim 2'); zlabel('Neural Dim 3')
xticklabels({}); yticklabels({}); zticklabels({});
view(az,el)
axis equal
xlim(xlims); ylim(ylims); zlim(zlims);
set(gca,'fontname','arial')
set(gca,'fontsize',fs)
if saveFig
    saveas(gcf,fullfile(saveDir,'schemPredTraj.svg'));
end

%% Local Functions for Rotation Matrices
function [Rx] = rotX(theta)
    Rx = [1 0 0; ...
        0 cosd(theta) -1*sind(theta);...
        0 sind(theta) cosd(theta)];
end

function [Ry] = rotY(theta)
    Ry = [cosd(theta) 0 sind(theta); ...
        0 1 0;...
        -sind(theta) 0 cosd(theta)];
end

function [Rz] = rotZ(theta)
    Rz = [cosd(theta) -1*sind(theta) 0; ...
        sind(theta) cosd(theta) 0;...
        0 0 1];
end