clear; clc; clf; close all

%% Setup saveFig   
    saveFig = false;
    saveDir = 'D:\Posture Paper\20230203\orthographic dpca pop traj';
    set(0, 'DefaultFigureRenderer', 'painters');
    
%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\orli.mat')
    pcmap = orli;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat')
    tcmap = customRainbow;
    
%% Get trajStruct
    %Load data
    dataset = 'E20200316';
    [Data,zScoreParams] = loadData(dataset);
    
    %Get trajStruct
    binWidth = 25;
    kernelStdDev = 25;
    trajFields = {'zSmoothFR'};
    trialInclStates = struct('trialName','','inclStates',[]);
    
    trialInclStates(1).trialName = {'GridTask_CO_Across_BC_ForceBar'};
    condFields = {{'target','targetData','target1ID'},{'posture','conditionData','postureID'}};
    
    %Execution
    trialInclStates(1).inclStates = {{'state','Step 1','first',-50},{'state','Step 2','first',250}};
    exTrajStruct = getTrajStruct20220419(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'matchConditions',false,'zScoreParams',zScoreParams);
    

%% Get numPts
    %Manually enter number of points to consider (execution)
   	numPts = 13;
           
    %Get posture and target lists
    postureList = unique([exTrajStruct.posture]);
    numPostures = size(postureList,2);
    targetList = unique([exTrajStruct.target]); 
    numTargets = size(targetList,2);
    numChannels = size(exTrajStruct(1).avgSmoothFR.traj,2);
    numConditions = size(exTrajStruct,2);

%% Marginalize neural responses
    %Form X
    trajStruct = exTrajStruct;
    X = NaN(numPts,numTargets,numPostures,numChannels);
    postureInd = 1;
    for posture = postureList
        targetInd = 1;
        for target = targetList
            traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgSmoothFR.traj;
            X(:,targetInd,postureInd,:) = traj(1:numPts,:); 
            targetInd = targetInd + 1;
        end
        postureInd = postureInd+1;
    end

    % Do marginalizations of X (Xdims: 1=time, 2=target, 3=posture, 4=channel)
    %Condition-invariant
    CIMargOffset = mean(X,[1 2 3],'omitnan');
    CIMargTraj = mean(X,[2 3],'omitnan') - CIMargOffset;
    %Posture and Target Traj
    tSig = mean(X,[3],'omitnan') - CIMargTraj  - CIMargOffset;
    pSig = mean(X,[2],'omitnan') - CIMargTraj - CIMargOffset;
tSigCI = std(X,0,3,'omitnan')./1.96;
pSigCI = std(X,0,2,'omitnan')./1.96;
%% Plot

     %Posture sig for select neurons
      time = exTrajStruct([exTrajStruct.posture]==1 & [exTrajStruct.target]==1).avgSmoothFR.timestamps(1:numPts);
    f = figure; fs = 14;
    f.Position = [100,-200,518,800];
    postureInd = 1;
    dimList = [3,49,78,76];
    postureInd = 1;
    for posture = postureList
        dimInd = 1;
        for dim = dimList
            subplot(4,1,dimInd); hold on
            shadedErrorBar(time,pSig(:,1,postureInd,dim),pSigCI(:,1,postureInd,dim),'lineprops',{'LineWidth',3,'Color',pcmap(posture,:)});
            dimInd = dimInd + 1;
        end
        postureInd = postureInd + 1;
    end

    dimInd = 1;
    for dim = dimList(1:4)
        subplot(4,1,dimInd); hold on
        ax = gca; xlimits = ax.XLim; ylimits = ax.YLim;
        ax.TickDir = 'out';
        set(gca,'fontname','arial'); set(gca,'fontsize',fs);
        if ismember(dimInd,[1,2,3])
            xticklabels({});
        end
        if dimInd == 4
            xticks([0,100,200])
            xlabel('time (ms)')
        end 

        dimInd = dimInd + 1;
    end
     
     
     
     
     %Target sig for select neurons
    f = figure; fs = 14;
    f.Position = [100,-200,518,800];
    targetInd = 1;
    dimList = [3,49,78,76];
    for target = targetList
        dimInd = 1;
        for dim = dimList
            subplot(4,1,dimInd); hold on
            shadedErrorBar(time,tSig(:,targetInd,1,dim),tSigCI(:,targetInd,1,dim),'lineprops',{'LineWidth',3,'Color',tcmap(targetInd,:)});
            dimInd = dimInd + 1;
        end
        targetInd = targetInd + 1;
    end
     
         dimInd = 1;
    for dim = dimList(1:4)
        subplot(4,1,dimInd); hold on
        ax = gca; xlimits = ax.XLim; ylimits = ax.YLim;
        ax.TickDir = 'out';
        set(gca,'fontname','arial'); set(gca,'fontsize',fs);
        if ismember(dimInd,[1,2,3])
            xticklabels({});
        end
        if dimInd == 4
            xticks([0,100,200])
            xlabel('time (ms)')
        end 

        dimInd = dimInd + 1;
    end
     
    

  