clear; clc; clf; close all

%% Parameters 
    numManPCs = 10;
    
%% Setup saveFig   
    saveFig = false;
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Paper\20220907\Error Reduction Figures';
    set(0, 'DefaultFigureRenderer', 'painters');

%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\orli.mat')
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat')
    pcmap = orli;
    scmap = customRainbow;
    
%% Setup resultStruct
    resultStruct = struct('animal',[],'dataset',[],'vPP',[],'vPT',[],'vTT',[],'vTP',[],'vPR',[],'vTR',[],'pAngle',[],'pAngleNull',[]);
    structInd = 1;
    
 %% Get traj struct   
    dataset = 'E20200316';
    [Data,zScoreParams] = loadData(dataset);

    %Get trajStruct
    binWidth = 25;
    kernelStdDev = 25;
    trajFields = {'zSmoothFR'};
    trialInclStates = struct('trialName','','inclStates',[]);
    switch dataset
        %BCI
        case {'E20200316','E20200317','E20200318'}
            trialInclStates(1).trialName = {'GridTask_CO_Across_BC_ForceBar'};
            condFields = {{'target','targetData','target1ID'},{'posture','conditionData','postureID'}};
            trialInclStates(1).inclStates = {{'state','Step 1','first',0},{'state','Step 2','first',0}};
    end
    trajStruct = getTrajStruct20220419(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'matchConditions',false,'zScoreParams',zScoreParams);


    % Get minimum number of timestamps in condition averages
    numTimestamps = [];
    for i = 1:size(trajStruct,2)
        numTimestamps(i) = length(trajStruct(i).avgSmoothFR.timestamps);
    end
    [minNumTimestamps,i] = min(numTimestamps);
    numPts = 10;    

    %Get posture and target lists
    postureList = unique([trajStruct.posture]);
    numPostures = size(postureList,2);
    targetList = unique([trajStruct.target]); 
    numTargets = size(targetList,2);
    numChannels = size(trajStruct(1).avgSmoothFR.traj,2);
    numConditions = size(trajStruct,2);


    %Project all data down to top PCs
    allTraj = NaN(numConditions*numPts,numChannels);
    j = 1;
    for i = 1:numConditions
       allTraj(j:j+numPts-1,:) = trajStruct(i).avgSmoothFR.traj(1:numPts,:);
       j = j + numPts;
    end
    [allPCs,~,~,~,explained,allMu] = pca(allTraj);
    for i = 1:size(trajStruct,2)
       trajStruct(i).avgSmoothFR.traj = (trajStruct(i).avgSmoothFR.traj-allMu)*allPCs(:,1:numManPCs); 
       for j = 1:size(trajStruct(i).allSmoothFR,2)
            trajStruct(i).allSmoothFR(j).traj = (trajStruct(i).allSmoothFR(j).traj-allMu)*allPCs(:,1:numManPCs); 
       end
    end
    
    %Get sigs
    [ciSig,pSig,tSig] = getMarg(trajStruct,numPts);
    pSigReshape = reshape(squeeze(pSig),[numPostures*numPts,numManPCs]);
    tSigReshape = reshape(squeeze(tSig),[numTargets*numPts,numManPCs]);
    
    %Plot
    %Target
    f = figure; f.Position = [50 50 1500 700];
    allTraj = targetMargTraj+targetMargOffset;
    targetInd = 1;
    for target = targetList
       traj = squeeze(allTraj(:,targetInd,:,:))*PCs;
       for i = 1:15
            subplot(3,5,i)
            plot(time,traj(:,i),'Color',tcmap(target,:))
            %shadedErrorBar(time,traj(:,i),targetErrBar(:,targetInd,i),'lineprops',{'LineWidth',2,'Color',tcmap(target,:)})
            ax = gca;
            %ax.YLim = [-250 250];
            hold on;
            xlabel('time (ms)')
            ylabel(['PC ',num2str(i)])
       end
       targetInd = targetInd + 1;
    end
    sgtitle('Target Marginalization')
    if saveFig
        saveas(gcf,fullfile(saveDir,task,'TargetMarginalization.jpg')); 
    end
    close all
    
    %% Local functions
    function [ciSig,pSig,tSig] = getMarg(trajStruct,numPts)
        %Get posture & target lists
        postureList = unique([trajStruct.posture]); numPostures = size(postureList,2);
        targetList = unique([trajStruct.target]); numTargets = size(targetList,2);
        numManPCs = size(trajStruct(1).avgSmoothFR.traj,2); 

        %Form X
        X = NaN(numPts,numTargets,numPostures,numManPCs);
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
        ciSig = mean(X,[2 3],'omitnan');
        %Posture and Target Traj
        tSig = mean(X,[3],'omitnan') - CIMargTraj  - CIMargOffset;
        pSig = mean(X,[2],'omitnan') - CIMargTraj - CIMargOffset;
    end