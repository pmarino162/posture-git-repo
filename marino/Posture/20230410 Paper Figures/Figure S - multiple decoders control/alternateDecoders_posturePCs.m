clear; clc; clf; close all

%% Setup saveFig   
    saveFig = false;
    saveDir = 'D:\Posture Paper\20230203\orthographic dpca pop traj';
    set(0, 'DefaultFigureRenderer', 'painters');
    
%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\orli.mat')
    pcmap = orli;
    
%% Load data; get trajStruct
    dataset = 'E20190729';
    [Data,zScoreParams] = loadData(dataset);

    %Get trajStruct
    binWidth = 25;
    kernelStdDev = 25;
    trajFields = {'zSmoothFR'};
    trialInclStates = struct('trialName','','inclStates',[]);
    condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'},{'decoder','conditionData','decoderPostureID'}};
    
    trialInclStates(1).trialName = {'CentetOut_BC_TouchBar'};
        trialInclStates(1).inclStates = {{'state','Step 1','first',0},{'state','Success with Reward','first',0}};
    trialInclStates(2).trialName = {'CenterOutCenter_BC_TouchBar'};
        trialInclStates(2).inclStates = {{'state','Step 1','first',0},{'state','Step 2','first',0}};
           
    trajStruct = getTrajStruct20220419(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'matchConditions',false,'zScoreParams',zScoreParams);
        
%% Fit posture PC on data with matched posture and decoder
    %Get posturePCs on blocks 1 and 2
    matchedTrajStruct = trajStruct([trajStruct.posture]==[trajStruct.decoder]); 
    [postureMargTraj,posturePCs] = getMarginalsAndPCs(matchedTrajStruct);

%% Add inidvidual trial projections to trajStruct
     for i = 1:size(trajStruct,2)
        numObs = size(trajStruct(i).allSmoothFR,2);
        pcProj = zeros(numObs,2);
        for j = 1:numObs
           obs = mean(trajStruct(i).allSmoothFR(j).traj); 
           pcProj(j,:) = obs*posturePCs(:,1:2);
        end
        trajStruct(i).pcProj = pcProj;
     end

%% Plot and save
    fs = 14;
    %pcmap = pcmap([1,5],:);
    figure; hold on;
    plotInd = 1;
    for posture = 1:2
        for decoder = 1:2
            tempTrajStruct = trajStruct([trajStruct.posture]==posture & [trajStruct.decoder]==decoder);
            plotData = vertcat(tempTrajStruct.pcProj);
            jitter = 0.5*(rand(1,length(plotData))-.5);
            plot(plotInd*ones(1,length(plotData))+jitter,plotData(:,1),'.','Color',pcmap(posture,:));
            plotInd = plotInd + 1;
        end
    end
    %ylabel('Projection of neural activity \n onto Posture Dim 1 (a.u.)')
    xlim([0.5 4.5])
    set(gca,'fontname','arial'); set(gca,'fontsize',fs)
    %             if saveFig
%                 saveas(gcf,fullfile(saveDir,[dataset,'_postureAxisVTrial.svg']));
%             end


%% Local functions
function [postureMargTraj,posturePCs] = getMarginalsAndPCs(trajStruct)
    %Get minimum number of timestamps in condition averages
    numTimestamps = [];
    for i = 1:size(trajStruct,2)
        numTimestamps(i) = length(trajStruct(i).avgSmoothFR.timestamps);
    end
    [minNumTimestamps,i] = min(numTimestamps);

    %Get posture and target lists
    postureList = unique([trajStruct.posture]);
    numPostures = size(postureList,2);
    targetList = unique([trajStruct.target]); 
    numTargets = size(targetList,2);
    numChannels = size(trajStruct(1).avgSmoothFR.traj,2);
    numConditions = size(trajStruct,2);

    %Form X, containing trial-averaged data for each condition
    X = NaN(minNumTimestamps,numTargets,numPostures,numChannels);
    postureInd = 1;
    for posture = postureList
        targetInd = 1;
        for target = targetList
            traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgSmoothFR.traj;
            X(:,targetInd,postureInd,:) = traj(1:minNumTimestamps,:); 
            targetInd = targetInd + 1;
        end
        postureInd = postureInd+1;
    end

    % Do marginalizations of X
    %Xdims: 1=time, 2=target, 3=posture, 4=channel
    %Condition-invariant
    CIMargOffset = mean(X,[1 2 3],'omitnan');
    CIMargTraj = mean(X,[2 3],'omitnan') - CIMargOffset;
    %Posture and Target Traj
    targetMargTraj = mean(X,[3],'omitnan') - CIMargTraj  - CIMargOffset;
    postureMargTraj = mean(X,[2],'omitnan') - CIMargTraj - CIMargOffset;

    %Compute PCs of marginalizations
    targetMargTraj = squeeze(targetMargTraj);
        targetMargTraj = reshape(targetMargTraj,[numTargets*minNumTimestamps,numChannels]);
        [targetPCs,~,~,~,~,targetMargMu] = pca(targetMargTraj); 
    postureMargTraj = squeeze(postureMargTraj);
        postureMargTraj = reshape(postureMargTraj,[numPostures*minNumTimestamps,numChannels]);
        [posturePCs,~,~,~,~,postureMargMu] = pca(postureMargTraj); 
        
    %Posture and Target Traj
    targetMargTraj = mean(X,[3],'omitnan') - CIMargTraj  - CIMargOffset;
    postureMargTraj = mean(X,[2],'omitnan') - CIMargTraj - CIMargOffset;
    postureMargTraj = squeeze(postureMargTraj);
end

