clear; clc; clf; close all;
    
%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customSummer.mat');
    cmap = customSummer;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat');
    tcmap = customRainbow;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\rainbow2d.mat')
    pcmap = customRainbow;
    
%% Set up save 
saveFig = true;
dirStr = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Presentations\20211229 - adam meeting - controlling for kin and are traj parallel';

%% Load Data 
    [Data] = loadEarlData20210706_20211210;
    
%% Get trajStruct 
    trajFields = {'smoothFR'};
    trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
    condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
    trialInclStates(1).trialName = {'GridReaching'};
    binWidth = 25;

    %Reach
    trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-150},{'state','Target Hold','first',250}};
    reachTrajStruct = getTrajStruct20211210(Data,condFields,trajFields,trialInclStates,binWidth,'matchConditions',false);


%% Get posture and target lists
    trajStruct = reachTrajStruct;
    postureList = unique([trajStruct.posture]); numPostures = size(postureList,2);
    targetList = unique([trajStruct.target]); numTargets = size(targetList,2);

%% Do PCA on Condition-averaged reaching data
    avgSmoothFR = [reachTrajStruct.avgSmoothFR];
    allAvgs = vertcat(avgSmoothFR.traj);
    [coeff,score,latent,tsquared,explained,mu] = pca(allAvgs);    
    
%% Add Projections to trajStruct and modelStruct
    for i = 1:size(reachTrajStruct,2)
       reachTrajStruct(i).PCA = reachTrajStruct(i).avgSmoothFR.traj*coeff;
    end    
    
    
%% Plot PCA Trajectories
    for posture = [1,4]
       for target = [1,3,5,7]
           traj = reachTrajStruct([reachTrajStruct.posture]==posture & [reachTrajStruct.target]==target).PCA; 
           traj = traj - traj(1,:); 
           plot3(traj(:,1),traj(:,2),traj(:,3),'Color',tcmap(target,:));
           hold on
           plot3(traj(1,1),traj(1,2),traj(1,3),'.','MarkerSize',15,'Color',tcmap(target,:),'MarkerFaceColor',tcmap(target,:));
           plot3(traj(end,1),traj(end,2),traj(end,3),'>','MarkerSize',5,'Color',tcmap(target,:),'MarkerFaceColor',tcmap(target,:));
           
       end
    end
    xlabel('PC 1')
    ylabel('PC 2')
    zlabel('PC 3')
 
    if saveFig
        titleStr = '\ReachingTop3PC';
        saveas(gcf,[dirStr,titleStr,'.fig'])       
    end
    
%% Plot PCA Timecourses
    for target = 1:8
        figure
        hold on
        for posture = 1:7
            if ~isempty(reachTrajStruct([reachTrajStruct.posture]==posture & [reachTrajStruct.target]==target))
               time = reachTrajStruct([reachTrajStruct.posture]==posture & [reachTrajStruct.target]==target).avgSmoothFR.timestamps; 
               traj = reachTrajStruct([reachTrajStruct.posture]==posture & [reachTrajStruct.target]==target).PCA; 
               traj = traj - traj(1,:);
               for i = 1:10
                   subplot(2,5,i)
                   plot(time,traj(:,i),'Color',tcmap(posture,:),'LineWidth',1);
                   hold on
               end
               sgtitle(['Target ',num2str(target)])
            end
        end
        
            %Get axis ranges
    for plotInd = 1:10
       ax = subplot(2,5,plotInd);
       range(plotInd) = ax.YLim(1,2)-ax.YLim(1,1);
       center(plotInd) = (ax.YLim(1,2)+ax.YLim(1,1))/2;
    end
    maxRange = max(range);
    for plotInd = 1:10
       ax = subplot(2,5,plotInd);
       ax.YLim = [center(plotInd)-maxRange/2 center(plotInd)+maxRange/2];
       xlabel('time (ms)')
        ylabel(['PC ',num2str(plotInd)])
    end
       
    
    if saveFig
        titleStr = ['\ReachingAllPCTarget',num2str(target)];
        saveas(gcf,[dirStr,titleStr,'.fig'])       
    end
    
    end
        
        
        
 
    
    

