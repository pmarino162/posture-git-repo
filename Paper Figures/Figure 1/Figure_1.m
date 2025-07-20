clear; clc; clf; close all;

%% Setup saveFig   
    saveFig = false;
    saveDir = "C:\Users\pmari\OneDrive\Documents\Posture\Paper\Reviewer responses\Analyses\Fig 1";
    set(0, 'DefaultFigureRenderer', 'painters');

%% Randomly generate disorganized trajectories (T/F)
    randomlyGenerateDisTraj = false;
    downsample = true;
    downsampleFactor = 2;
    
%% Setup colormaps 
    [pcmap,tcmap,rainbow] = getColorMaps(2);
    
%%  Create base trajectory
    t = 0:1:100;
    numPts = length(t);
    baseTraj = zeros(numPts,3);
    baseTraj(:,1) = 0;
    baseTraj(:,2) = -cos(pi*t/100)+1;
    baseTraj(:,3) = sin(pi*t/100);
        
%% Set up posture and target lists
    postureList = 1:2;
    targetList = 1:2;
    
%% Create orgTrajStruct
    orgTrajStruct = struct('posture',[],'target',[],'traj',[]);
    trajStructInd = 1;
    
    % Create a single rotation that applies to everything
    thetaXFrac = 0.01;
    thetaYFrac = 0.07;
    thetaZFrac = 0.08;
    thetaX = 360*thetaXFrac;
    thetaY = 360*thetaYFrac;
    thetaZ = 360*thetaZFrac;
    preRotMat = rotX(thetaX)*rotY(thetaY)*rotZ(thetaZ);
    
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
            orgTrajStruct(trajStructInd).traj = baseTraj*stretchMat*rotMat*preRotMat + postureOffset;
            trajStructInd = trajStructInd + 1;
        end
    end
    
%% Create disTrajStruct

    if randomlyGenerateDisTraj == true
        disTrajStruct = struct('posture',[],'target',[],'traj',[],'offsetX',[],'offsetY',[],'offsetZ',[],'thetaXFrac',[],'thetaYFrac',[],'thetaZFrac',[]);
        trajStructInd = 1;
        for posture = postureList
                % Randomly generate offset per posture
                offsetX = rand;
                offsetY = rand;
                offsetZ = rand;
            for target = targetList
                % Randomly generate rotation generate per target
                thetaXFrac = rand;
                thetaYFrac = rand;
                thetaZFrac = rand;
                % Get posture offset
                postureOffset = [posture+offsetX,offsetY,offsetZ];
                %Get rotation matrix
                thetaX = 360*thetaXFrac;
                thetaY = 360*thetaYFrac;
                thetaZ = 360*thetaZFrac;
                rotMat = rotX(thetaX)*rotY(thetaY)*rotZ(thetaZ);
                %Apply stretch
                %stretchFac = rand+1;
                stretchFac = 1;
                stretchMat = stretchFac.*eye(3);
                %Fill disTrajStruct
                disTrajStruct(trajStructInd).posture = posture;
                disTrajStruct(trajStructInd).target = target;
                disTrajStruct(trajStructInd).traj = baseTraj*stretchMat*rotMat + postureOffset;
                trajStructInd = trajStructInd + 1;
            end
        end
    else
        disTrajStruct = load("savedDisTrajStruct20250122.mat");
        disTrajStruct = disTrajStruct.disTrajStruct;

    end
    
%% Plot trajectories
fs = 14; %Font size
ds = 25; %Dot Size
as = 10; %Arrow Size
lw = 3.5; %LineWidth
az = 121.9112; el = 23.2414;

xLen = 2;
yLen = 5;
zLen = 4;

% Organized
figure; hold on;
xDim = 3; yDim = 2; zDim = 1;
for posture = [1,2]
    for target = targetList
        traj = orgTrajStruct([orgTrajStruct.posture]==posture & [orgTrajStruct.target]==target).traj;
        if downsample == true
            traj = traj(1:downsampleFactor:end,:);
        end
        if target == 1
            plot3(traj(:,xDim),traj(:,yDim),traj(:,zDim),'Color',pcmap(posture,:),'LineWidth',lw)
        elseif target == 2
            plot3(traj(1:10:end,xDim),traj(1:10:end,yDim),traj(1:10:end,zDim),'--','Color',pcmap(posture,:),'LineWidth',lw)
        end
        plot3(traj(1,xDim),traj(1,yDim),traj(1,zDim),'.','MarkerEdgeColor',pcmap(posture,:),'MarkerSize',ds);
        plot3(traj(end,xDim),traj(end,yDim),traj(end,zDim),'v','MarkerEdgeColor',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:),'MarkerSize',as);
    end
end
%xlabel('x');
%ylabel('y');
%zlabel('z');
xticklabels({}); yticklabels({}); zticklabels({});
ax = gca;

view(az,el)
axis equal
set(gca,'fontname','arial'); set(gca,'fontsize',fs)
ax.ZLim = [0.5 2.5];
orgXLims = ax.XLim; orgYLims = ax.YLim; orgZLims = ax.ZLim; 

overallMean = getOverallMean(orgTrajStruct);
orgXLims = [overallMean(xDim) - xLen/2, overallMean(xDim) + xLen/2];
orgYLims = [overallMean(yDim) - yLen/2, overallMean(yDim) + yLen/2];
orgZLims = [overallMean(zDim) - zLen/2, overallMean(zDim) + zLen/2];

ax.XLim = orgXLims;
ax.YLim = orgYLims;
ax.ZLim = orgZLims;

%Setup tickLists
xtickList = [orgXLims(1):.5:orgXLims(2)];
ytickList = [orgYLims(1):.5:orgYLims(2)];
ztickList = [orgZLims(1):.5:orgZLims(2)];
xticks(xtickList); yticks(ytickList); zticks(ztickList);
grid on
if saveFig
    saveas(gcf,fullfile(saveDir,'org.svg'));
end

% Disorganized
xLen = 2;
yLen = 5;
zLen = 5;
figure; hold on;
xDim = 3; yDim = 2; zDim = 1;
for posture = [1,2]
    for target = targetList
        traj = disTrajStruct([disTrajStruct.posture]==posture & [disTrajStruct.target]==target).traj;
        if downsample == true
            traj = traj(1:downsampleFactor:end,:);
        end
        if target == 1
            plot3(traj(:,xDim),traj(:,yDim),traj(:,zDim),'Color',pcmap(posture,:),'LineWidth',lw)
        elseif target == 2
            plot3(traj(1:10:end,xDim),traj(1:10:end,yDim),traj(1:10:end,zDim),'--','Color',pcmap(posture,:),'LineWidth',lw)
        end
        plot3(traj(1,xDim),traj(1,yDim),traj(1,zDim),'.','MarkerEdgeColor',pcmap(posture,:),'MarkerSize',ds);
        plot3(traj(end,xDim),traj(end,yDim),traj(end,zDim),'v','MarkerEdgeColor',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:),'MarkerSize',as);
    end
end
 xticklabels({}); yticklabels({}); zticklabels({});
ax = gca;
%xlabel('x');
%ylabel('y');
%zlabel('z');
view(az,el)
axis equal
set(gca,'fontname','arial'); set(gca,'fontsize',fs)
%ax.ZLim = [0.5 2.5];
% disXLims = ax.XLim; disYLims = ax.YLim; disZLims = ax.ZLim; 
% %Setup tickLists
% xtickList = [disXLims(1):.5:disXLims(2)];
% ytickList = [disYLims(1):.5:disYLims(2)];
% ztickList = [disZLims(1):.5:disZLims(2)];
% xticks(xtickList); yticks(ytickList); zticks(ztickList);


overallMean = getOverallMean(disTrajStruct);
disXLims = [overallMean(xDim) - xLen/2, overallMean(xDim) + xLen/2];
disYLims = [overallMean(yDim) - yLen/2, overallMean(yDim) + yLen/2];
disZLims = [overallMean(zDim) - zLen/2, overallMean(zDim) + zLen/2];

ax.XLim = disXLims;
ax.YLim = disYLims;
ax.ZLim = disZLims;

%Setup tickLists
xtickList = [orgXLims(1):.5:orgXLims(2)];
ytickList = [orgYLims(1):.5:orgYLims(2)];
ztickList = [orgZLims(1):.5:orgZLims(2)];
xticks(xtickList); yticks(ytickList); zticks(ztickList);
grid on
if saveFig
    saveas(gcf,fullfile(saveDir,'disorg.svg'));
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

function [overallMean] = getOverallMean(trajStruct)
    allData = [];
    for i = 1:size(trajStruct,2)
        allData = vertcat(allData,trajStruct(i).traj);
    end
    overallMean = mean(allData,1);
end