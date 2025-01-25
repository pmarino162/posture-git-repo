clear; clc; clf; close all;

%% Setup saveFig   
    saveFig = false;
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Paper\20231002\Figure 1';
    set(0, 'DefaultFigureRenderer', 'painters');
    
%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\orli.mat')
    pcmap = orli;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat')
    tcmap = customRainbow;
    posture2color = [5,1];
    
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
    

    
%% Plot Trajectories
fs = 14; %Font size
ds = 25; %Dot Size
as = 10; %Arrow Size
lw = 3.5; %LineWidth
az = 121.9112; el = 23.2414;


%Organized
figure; hold on;
xDim = 3; yDim = 2; zDim = 1;
for posture = [1,2]
    for target = targetList
        traj = orgTrajStruct([orgTrajStruct.posture]==posture & [orgTrajStruct.target]==target).traj;
        if target == 1
            plot3(traj(:,xDim),traj(:,yDim),traj(:,zDim),'Color',pcmap(posture2color(posture),:),'LineWidth',lw)
        elseif target == 2
            plot3(traj(1:10:end,xDim),traj(1:10:end,yDim),traj(1:10:end,zDim),'--','Color',pcmap(posture2color(posture),:),'LineWidth',lw)
        end
        plot3(traj(1,xDim),traj(1,yDim),traj(1,zDim),'.','MarkerEdgeColor',pcmap(posture2color(posture),:),'MarkerSize',ds);
        plot3(traj(end,xDim),traj(end,yDim),traj(end,zDim),'v','MarkerEdgeColor',pcmap(posture2color(posture),:),'MarkerFaceColor',pcmap(posture2color(posture),:),'MarkerSize',as);
    end
end
%xlabel('Neural Dim 1'); ylabel('Neural Dim 2'); zlabel('Neural Dim 3')
xticklabels({}); yticklabels({}); zticklabels({});
ax = gca;

view(az,el)
axis equal
set(gca,'fontname','arial'); set(gca,'fontsize',fs)
ax.ZLim = [0.5 2.5];
orgXLims = ax.XLim; orgYLims = ax.YLim; orgZLims = ax.ZLim; 
%Setup tickLists
xtickList = [orgXLims(1):.5:orgXLims(2)];
ytickList = [orgYLims(1):.5:orgYLims(2)];
ztickList = [orgZLims(1):.5:orgZLims(2)];
xticks(xtickList); yticks(ytickList); zticks(ztickList);
grid on
if saveFig
    saveas(gcf,fullfile(saveDir,'org.svg'));
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