clear; clc; clf; close all;

%% Setup saveFig   
    saveFig = true;
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Paper\20220907\Figure 4\Curved Traj Schematic';
    set(0, 'DefaultFigureRenderer', 'painters');
    
%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\orli.mat')
    pcmap = orli;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat')
    tcmap = customRainbow;

%%  Create base trajectory
    t = 0:10:100;
    numPts = length(t);
    baseTraj = zeros(numPts,3);
    baseTraj(:,1) = 0;
    baseTraj(:,2) = -cos(pi*t/100)+1;
    baseTraj(:,3) = sin(pi*t/100);
        
%% Set up posture and target lists
    postureList = 1:3;
    targetList = 1:2;
    
%% Create trajStruct
    trajStruct = struct('posture',[],'target',[],'traj',[]);
    trajStructInd = 1;
    for posture = postureList
        for target = targetList
            %Get posture offset
            postureOffset = [posture,0,0];
            %Get rotation matrix
            rotMat = rotX(180*(target-1));
            if target == 2
                rotMat = rotMat*rotY(180);
            end
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
az = 40; el = 20;

pcmap2 = pcmap; 
pcmap2(2,:) = pcmap(3,:);
pcmap2(3,:) = pcmap(5,:);


%All trajectories 
figure; hold on;
xDim = 3; yDim = 2; zDim = 1;
for posture = postureList
    for target = targetList
        traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).traj;
        if posture == 3 && target == 1
            plot3(traj(:,xDim),traj(:,yDim),traj(:,zDim),'--','Color',pcmap2(posture,:),'LineWidth',lw);
        else
            plot3(traj(:,xDim),traj(:,yDim),traj(:,zDim),'Color',pcmap2(posture,:),'LineWidth',lw);
        end
        plot3(traj(1,xDim),traj(1,yDim),traj(1,zDim),'.','MarkerEdgeColor',pcmap2(posture,:),'MarkerSize',ds);
        plot3(traj(end,xDim),traj(end,yDim),traj(end,zDim),'v','MarkerEdgeColor',pcmap2(posture,:),'MarkerFaceColor',pcmap2(posture,:),'MarkerSize',as);
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

%All trajectories except p3t1
figure; hold on;
xDim = 3; yDim = 2; zDim = 1;
for posture = postureList
    for target = targetList
        if posture == 3 && target == 1
        else
            traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).traj;
            plot3(traj(:,xDim),traj(:,yDim),traj(:,zDim),'Color',pcmap2(posture,:),'LineWidth',lw)
            plot3(traj(1,xDim),traj(1,yDim),traj(1,zDim),'.','MarkerEdgeColor',pcmap2(posture,:),'MarkerSize',ds);
            plot3(traj(end,xDim),traj(end,yDim),traj(end,zDim),'v','MarkerEdgeColor',pcmap2(posture,:),'MarkerFaceColor',pcmap2(posture,:),'MarkerSize',as);
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
    saveas(gcf,fullfile(saveDir,'schemAllTrajExceptP3T1.svg'));
end

%Only P1 t1; and IC 
figure; hold on;
xDim = 3; yDim = 2; zDim = 1;
for posture = 1
    for target = 1
        traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).traj;
        traj = traj + [0.5 0 0];
        plot3(traj(:,xDim),traj(:,yDim),traj(:,zDim),'Color',tcmap(8,:),'LineWidth',lw)
        plot3(traj(1,xDim),traj(1,yDim),traj(1,zDim),'.','MarkerEdgeColor',tcmap(8,:),'MarkerSize',ds);
        plot3(traj(end,xDim),traj(end,yDim),traj(end,zDim),'v','MarkerEdgeColor',tcmap(8,:),'MarkerFaceColor',tcmap(8,:),'MarkerSize',as);
    end
end
for posture = 3
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
xDim = 3; yDim = 2; zDim = 1;
for posture = 3
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