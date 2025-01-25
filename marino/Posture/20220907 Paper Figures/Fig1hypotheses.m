clear; clc; clf; close all;

%% Setup saveFig   
    saveFig = false;
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Paper\20220907\Figure 1';
    set(0, 'DefaultFigureRenderer', 'painters');
    
%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\orli.mat')
    pcmap = orli;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat')
    tcmap = customRainbow;

%%  Create base trajectory
    t = 0:1:100;
    numPts = length(t);
    baseTraj = zeros(numPts,3);
    baseTraj(:,1) = 0;
    baseTraj(:,2) = -cos(pi*t/100)+1;
    baseTraj(:,3) = sin(pi*t/100);
        
%% Set up posture and target lists
    postureList = 1:3;
    targetList = 1:2;
    
%% Create orgTrajStruct
    orgTrajStruct = struct('posture',[],'target',[],'traj',[]);
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
            %Fill orgTrajStruct
            orgTrajStruct(trajStructInd).posture = posture;
            orgTrajStruct(trajStructInd).target = target;
            orgTrajStruct(trajStructInd).traj = baseTraj*stretchMat*rotMat + postureOffset;
            trajStructInd = trajStructInd + 1;
        end
    end
    
%% Create disTrajStruct
    disTrajStruct = struct('posture',[],'target',[],'traj',[]);
    trajStructInd = 1;
    for posture = postureList
        randNum = rand;
        for target = targetList
            %Get posture offset
            %postureOffset = [posture,0,0];
            
            postureOffset = [posture+randNum,randNum,randNum];
            %Get rotation matrix
            thetaX = 360*rand;
            thetaY = 360*rand;
            thetaZ = 360*rand;
            rotMat = rotX(thetaX)*rotY(thetaY)*rotZ(thetaZ);
            %Apply stretch
            stretchFac = rand+1;
            stretchMat = stretchFac.*eye(3);
            %Fill orgTrajStruct
            disTrajStruct(trajStructInd).posture = posture;
            disTrajStruct(trajStructInd).target = target;
            disTrajStruct(trajStructInd).traj = baseTraj*stretchMat*rotMat + postureOffset;
            trajStructInd = trajStructInd + 1;
        end
    end
   
    
%% Plot Trajectories
fs = 14; %Font size
ds = 25; %Dot Size
as = 10; %Arrow Size
lw = 3.5; %LineWidth
az = 40; el = 20;

%Organized
figure; hold on;
xDim = 3; yDim = 2; zDim = 1;
for posture = postureList
    for target = targetList
        traj = orgTrajStruct([orgTrajStruct.posture]==posture & [orgTrajStruct.target]==target).traj;
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
    saveas(gcf,fullfile(saveDir,'orgHyp.svg'));
end

%Disorganized
figure; hold on;
xDim = 3; yDim = 2; zDim = 1;
for posture = postureList
    for target = targetList
        traj = disTrajStruct([disTrajStruct.posture]==posture & [disTrajStruct.target]==target).traj;
        plot3(traj(:,xDim),traj(:,yDim),traj(:,zDim),'Color',pcmap(posture,:),'LineWidth',lw)
        plot3(traj(1,xDim),traj(1,yDim),traj(1,zDim),'.','MarkerEdgeColor',pcmap(posture,:),'MarkerSize',ds);
        plot3(traj(end,xDim),traj(end,yDim),traj(end,zDim),'v','MarkerEdgeColor',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:),'MarkerSize',as);
    end
end
grid on
xlabel('X'); ylabel('Y'); zlabel('Z')
xlabel('Neural Dim 1'); ylabel('Neural Dim 2'); zlabel('Neural Dim 3')
xticklabels({}); yticklabels({}); zticklabels({});
view(az,el)
axis equal
set(gca,'fontname','arial')
set(gca,'fontsize',fs)
if saveFig
    saveas(gcf,fullfile(saveDir,'disOrgHyp.svg'));
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