clear; clc; clf; close all

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
    resultStruct = struct('animal',[],'dataset',[],'vPP',[],'vPT',[],'vTT',[],'vTP',[],'vPR',[],'vTR',[]);
    structInd = 1;
    
%% Parameters
   %manVarThreshold = 90;
   numManPCs = 12;
   numIterations = 20;
   numPdims = 2;
   numTdims = 2;
   numChanceReps = 20;
   
%% Run loop for each dataset
    reachDatasetList = {'E20210706','E20210707','E20210708',...
        'N20190222','N20190226','N20190227','N20190228','N20190307'...
        'R20200221','R20200222'};
    
    bciDatasetList = {'E20200316','E20200317','E20200318','E20200319','N20171215','N20180221','R20201020','R20201021'};
    
    for datasetList = bciDatasetList%{'E20200316','E20200317','R20201020','R20201021'}% reachDatasetList%{'E20210707','N20190226','R20200221'}%{'E20200316','N20171215','R20201020'}%{ 'R20200221'}%bciDatasetList%{'E20210706','E20210707','E20210708','E20210709'}%reachDatasetList
        %Load data
        dataset = datasetList{1,1};
        [Data,zScoreParams] = loadData(dataset);
        %Get trajStruct
        binWidth = 25;
        kernelStdDev = 25;
        trajFields = {'zSmoothFR'};
        trialInclStates = struct('trialName','','inclStates',[]);
        switch dataset
            %BCI
            case {'E20200316','E20200317','E20200318','E20200319'}
                trialInclStates(1).trialName = {'GridTask_CO_Across_BC_ForceBar'};
                condFields = {{'target','targetData','target1ID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'state','Step 1','first',0},{'state','Step 2','first',0}};
                task = 'bci';
            case {'N20171215','N20180221'}
                trialInclStates(1).trialName = {'Nigel Posture BC Center Out'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'state','Cursor Freeze','first',0},{'state','Target Hold','first',0}};
                task = 'bci';
            case {'R20201020','R20201021'}
                trialInclStates(1).trialName = {'Rocky Posture BC Center Out'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'state','React','first',0},{'state','Hold','first',0}};
                task = 'bci';
            %Reach
            case {'E20210706','E20210707','E20210708','E20210709','E20210710'}
                trialInclStates(1).trialName = {'GridReaching'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',0}};
                task = 'reach';
            case {'R20200221','R20200222'}
                trialInclStates(1).trialName = {'Rocky Dissociation'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',0}};
                task = 'reach';
            case {'N20190222','N20190226','N20190227','N20190228','N20190301','N20190305','N20190306','N20190307'}
                trialInclStates(1).trialName = {'Nigel Dissociation'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',0}};
                task = 'reach';
        end
        
        %Remove any long or short trials - Reach
        if strcmpi(task,'reach')
            kinData = [Data.kinData];
            rxnTime = [kinData.rxnTime];
            reachTime = [kinData.reachTime];
            figure
                histogram(rxnTime)
                xlabel('Reaction Time')
                ylabel('Number of conditions')
            figure
                histogram(reachTime)
                xlabel('Reach Time')
                ylabel('Number of conditions')
           rxnCutoff = prctile(rxnTime,95);
           reachCutoff = prctile(reachTime,95);
           rmTrials =  unique([find(rxnTime > rxnCutoff), find(reachTime > reachCutoff)]);
           Data(rmTrials) = [];
       end

       trajStruct = getTrajStruct20220419(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'matchConditions',false,'zScoreParams',zScoreParams);
      
        %Remove long/short trials, BCI
        if strcmpi(task,'bci')
            reachTime = [];
            for i = 1:size(trajStruct,2)
               numTraj = size(trajStruct(i).allSmoothFR,2);
               for j = 1:numTraj
                  reachTime = [reachTime,size(trajStruct(i).allSmoothFR(j).traj,1)]; 
               end
            end     
            figure
                histogram(reachTime)
                xlabel('Reach Time')
                ylabel('Number of trials')
           reachCutoff = prctile(reachTime,95);
           for i = 1:size(trajStruct,2)
               numTraj = size(trajStruct(i).allSmoothFR,2);
               reachTime = [];
               for j = 1:numTraj
                  reachTime = [reachTime,size(trajStruct(i).allSmoothFR(j).traj,1)]; 
               end
               rmTrials =  find(reachTime > reachCutoff);
               trajStruct(i).allSmoothFR(rmTrials) = [];
            end  
        end
 
        % Get number of trials for each condition, numTimestamps in each
        % trail
        numCondTraj = [];
        for i = 1:size(trajStruct,2)
           numTraj = size(trajStruct(i).allSmoothFR,2);
           numCondTraj = [numCondTraj,numTraj];
        end        
        figure
        histogram(numCondTraj)
        xlabel('Number of trials')
        ylabel('Number of conditions')

        %Remove any conditions for which there weren't enough trials
        cutoffNumTraj = 10;
        trajStruct = trajStruct(numCondTraj >= cutoffNumTraj);
        
        %Keep only postures with all targets
        postureList = unique([trajStruct.posture]);
        targetList = unique([trajStruct.target]); 
        keepPosture = [];
        for posture = postureList
            tempTrajStruct = trajStruct([trajStruct.posture]==posture);
            postureTargetList = [tempTrajStruct.target];
            if isequal(postureTargetList,targetList)
                keepPosture = [posture,keepPosture];
            end
        end
        trajStruct = trajStruct(ismember([trajStruct.posture],keepPosture));
        
        %Get minimum number of trials and timestamps
        numTimestamps = [];
        numCondTraj = [];
        for i = 1:size(trajStruct,2)
           numTraj = size(trajStruct(i).allSmoothFR,2);
           numCondTraj = [numCondTraj,numTraj];
           numTimestamps = [numTimestamps,size(trajStruct(i).avgSmoothFR.timestamps,2)]; 
        end
        [minNumTimestamps,i] = min(numTimestamps);
        [minNumCondTraj,i] = min(numCondTraj);
        numPts = minNumTimestamps;
        
        %Manually enter number of points to consider
         switch dataset
            %BCI
            case {'E20200316','E20200317','E20200318','E20200319'}
                numPts = 8;
            case {'N20171215','N20180221'}
                numPts = 12;
            case {'R20201020','R20201021'}
                numPts = 10;
            %Reach
            case {'E20210706','E20210707','E20210708','E20210709','E20210710'}
                numPts = 8;
            case {'R20200221','R20200222'}
                numPts = 8;
            case {'N20190222','N20190226','N20190227','N20190228','N20190301','N20190305','N20190306','N20190307'}
                numPts = 8;
        end
        
        %Get posture and target lists
        postureList = unique([trajStruct.posture]);
        numPostures = size(postureList,2);
        targetList = unique([trajStruct.target]); 
        numTargets = size(targetList,2);
        numChannels = size(trajStruct(1).avgSmoothFR.traj,2);
        numConditions = size(trajStruct,2);

        %Project all data down to top PCs
        allTraj = NaN(numConditions*minNumTimestamps,numChannels);
        j = 1;
        for i = 1:numConditions
           allTraj(j:j+minNumTimestamps-1,:) = trajStruct(i).avgSmoothFR.traj(1:minNumTimestamps,:);
           j = j + minNumTimestamps;
        end
        [allPCs,~,~,~,explained,allMu] = pca(allTraj);
%         numManPCs = 1;
%         while sum(explained(1:numManPCs)) < manVarThreshold
%             numManPCs = numManPCs + 1;
%         end
        for i = 1:size(trajStruct,2)
           trajStruct(i).avgSmoothFR.traj = (trajStruct(i).avgSmoothFR.traj-allMu)*allPCs(:,1:numManPCs); 
           for j = 1:size(trajStruct(i).allSmoothFR,2)
                trajStruct(i).allSmoothFR(j).traj = (trajStruct(i).allSmoothFR(j).traj-allMu)*allPCs(:,1:numManPCs); 
           end
        end
        numDims = numManPCs;
        numChannels = numManPCs;

        %Do computation
        numSample = floor(minNumCondTraj/2); %Number of trajectories in each draw
        vPP = NaN(1,numIterations); 
        vPT = NaN(1,numIterations); 
        vTT = NaN(1,numIterations); 
        vTP = NaN(1,numIterations); 
        vPR = NaN(1,numIterations); 
        vTR = NaN(1,numIterations); 
        for i = 1:numIterations
            i
            %Split trajStruct into two groups
            trajStruct1 = trajStruct;
            trajStruct2 = trajStruct;
            for j = 1:size(trajStruct,2)
                numTraj = size(trajStruct(j).allSmoothFR,2);
                sampInd1 = randsample(numTraj,numSample);
                trajStruct1(j).allSmoothFR = trajStruct(j).allSmoothFR(sampInd1);
                [trajStruct1(j).avgSmoothFR.traj,trajStruct1(j).avgSmoothFR.timestamps] = getAvgTraj20211210(trajStruct1(j).allSmoothFR,binWidth);               
                numTrajRemaining = numTraj - numSample;
                sampInd2 = randsample(numTrajRemaining,numSample);
                remainingInd = setdiff(1:numTraj,sampInd1);
                sampInd2 = remainingInd(sampInd2);
                trajStruct2(j).allSmoothFR = trajStruct(j).allSmoothFR(sampInd2);
                [trajStruct2(j).avgSmoothFR.traj,trajStruct2(j).avgSmoothFR.timestamps] = getAvgTraj20211210(trajStruct2(j).allSmoothFR,binWidth);   
            end
            
            %Get P and T signals for each group
            [pSig1,tSig1] = getPandTsig(trajStruct1,numPts);
                pSig1Reshape = reshape(squeeze(pSig1),[numPostures*numPts,numChannels]);
                tSig1Reshape = reshape(squeeze(tSig1),[numTargets*numPts,numChannels]);
            [pSig2,tSig2] = getPandTsig(trajStruct2,numPts);
                pSig2Reshape = reshape(squeeze(pSig2),[numPostures*numPts,numChannels]);
                tSig2Reshape = reshape(squeeze(tSig2),[numTargets*numPts,numChannels]);
                
            %Get P and T dims for group 1
            [pDims,~,~,~,explainedP,pSigMu] = pca(pSig1Reshape); 
            [tDims,~,~,~,explainedT,tSigMu] = pca(tSig1Reshape); 
            
            %Get group 2 signal variances in each space
            CP = cov(pSig2Reshape);
            DP = pDims(:,1:numPdims);
            CT = cov(tSig2Reshape);
            DT = tDims(:,1:numTdims);
            
            vPP(i) = trace(DP'*CP*DP)./trace(CP);
            vPT(i) = trace(DT'*CP*DT)./trace(CP);
            vTT(i) = trace(DT'*CT*DT)./trace(CT);
            vTP(i) = trace(DP'*CT*DP)./trace(CT);
             
            %Compute chance distribution
            tempvPR = NaN(1,numChanceReps);
            tempvTR = NaN(1,numChanceReps);
            for j = 1:numChanceReps
                v1 = normrnd(0,1,numChannels,1);
                v2 = normrnd(0,1,numChannels,1);
                DR = orth([v1,v2]);
                tempvPR(j) = trace(DR'*CP*DR)/trace(CP);
                tempvTR(j) = trace(DR'*CT*DR)/trace(CT);
            end       
            vPR(i) = mean(tempvPR);
            vTR(i) = mean(tempvTR);
        end

        %Add results to resultStruct
        resultStruct(structInd).dataset = dataset;
        resultStruct(structInd).animal = dataset(1);
        resultStruct(structInd).vPP = mean(vPP);
        resultStruct(structInd).vPT = mean(vPT);
        resultStruct(structInd).vTT = mean(vTT);
        resultStruct(structInd).vTP = mean(vTP);
        resultStruct(structInd).vPR = mean(vPR);
        resultStruct(structInd).vTR = mean(vTR);
        structInd = structInd + 1;
    end
    
%% Plot results
    mcmap = [0.4*ones(1,3); 0.6*ones(1,3); 0.9*ones(1,3)];
    
    %Across datasets
    for i = 1:numel(resultStruct)
        resultMonkeyList{i} = resultStruct(i).animal;
    end
    monkeyList = unique(resultMonkeyList);
    numMonkeys = numel(monkeyList);
    
    %Plotting parameters
    fs = 12;
    offset = 1/7;   
    sessionOffset = 1/40;
    
    %Bar plot for posture signal
    f=figure; hold on;
    monkeyInd = 1;
    for monkey = monkeyList
       tempResultStruct = resultStruct(strcmp(resultMonkeyList,monkey{1,1}));
       vPP = [tempResultStruct.vPP]*100;
            vPPMean = mean(vPP);
            vPPStd = 1.96*std(vPP)./sqrt(length(vPP));
       vPT = [tempResultStruct.vPT]*100;
            vPTMean = mean(vPT);
            vPTStd = 1.96*std(vPT)./sqrt(length(vPT));
       vPR = [tempResultStruct.vPR]*100;
            vPRMean = mean(vPR);
            vPRStd = 1.96*std(vPR)./sqrt(length(vPR));
      
      bar(monkeyInd,vPPMean,'FaceColor',mcmap(monkeyInd,:),'EdgeColor',[0 0 0]);
      hold on;
            plot(monkeyInd,vPP,'.k','MarkerSize',15);
      bar(monkeyInd+3,vPTMean,'FaceColor',mcmap(monkeyInd,:),'EdgeColor',[0 0 0]);
            plot(monkeyInd+3,vPT,'.k','MarkerSize',15);

       monkeyInd = monkeyInd + 1; 
    end
          ax = gca;
      curXLims = ax.XLim;
      shadedErrorBar(curXLims,[vPRMean,vPRMean],vPRStd,'lineprops',{'--','LineWidth',1.5,'Color',[0.5 0.5 0.5]});

      %Bar plot for goal signal
      f=figure; hold on;
        monkeyInd = 1;
    for monkey = monkeyList
       tempResultStruct = resultStruct(strcmp(resultMonkeyList,monkey{1,1}));
       vTT = [tempResultStruct.vTT]*100;
            vTTMean = mean(vTT);
            vTTStd = 1.96*std(vTT)./sqrt(length(vTT));
       vTP = [tempResultStruct.vTP]*100;
            vTPMean = mean(vTP);
            vTPStd = 1.96*std(vTP)./sqrt(length(vTP));
       vTR = [tempResultStruct.vTR]*100;
            vTRMean = mean(vTR);
            vTRStd = 1.96*std(vTR)./sqrt(length(vTR));
      
      bar(monkeyInd,vTTMean,'FaceColor',mcmap(monkeyInd,:),'EdgeColor',[0 0 0]);
      hold on;
            plot(monkeyInd,vTT,'.k','MarkerSize',15);
      bar(monkeyInd+3,vTPMean,'FaceColor',mcmap(monkeyInd,:),'EdgeColor',[0 0 0]);
            plot(monkeyInd+3,vTP,'.k','MarkerSize',15);

       monkeyInd = monkeyInd + 1; 
    end
          ax = gca;
      curXLims = ax.XLim;
      shadedErrorBar(curXLims,[vTRMean,vTRMean],vTRStd,'lineprops',{'--','LineWidth',1.5,'Color',[0.5 0.5 0.5]});

      
      
      
      %     xlim([0.5,5.5])
%     fs = 14;
%     set(gca,'fontname','arial')
%     set(gca,'fontsize',fs)
%     ylabel('Variance Captured (%)')
%     yticks([0:25:100])
%     set(gca,'TickDir','out');
%     
%     %Bar plot for goal signal
%     f=figure; hold on;
%     monkeyInd = 1;
%     for monkey = monkeyList
%        tempResultStruct = resultStruct(strcmp(resultMonkeyList,monkey{1,1}));
%        vPP = [tempResultStruct.vPP]*100;
%             vPPMean = mean(vPP);
%             vPPStd = 1.96*std(vPP)./sqrt(length(vPP));
%        vPT = [tempResultStruct.vPT]*100;
%             vPTMean = mean(vPT);
%             vPTStd = 1.96*std(vPT)./sqrt(length(vPT));
%        vPR = [tempResultStruct.vPR]*100;
%             vPRMean = mean(vPR);
%             vPRStd = 1.96*std(vPR)./sqrt(length(vPR));
%        vTT = [tempResultStruct.vTT]*100;
%             vTTMean = mean(vTT);
%             vTTStd = 1.96*std(vTT)./sqrt(length(vTT));
%        vTP = [tempResultStruct.vTP]*100;
%             vTPMean = mean(vTP);
%             vTPStd = 1.96*std(vTP)./sqrt(length(vTP));
%        vTR = [tempResultStruct.vTR]*100;
%             vTRMean = mean(vTR);
%             vTRStd = 1.96*std(vTR)./sqrt(length(vTR));
%            
%        bar(
%        
% 
%        monkeyInd = monkeyInd + 1; 
%     end
%% Local functions

function [pSig,tSig] = getPandTsig(trajStruct,numPts)
    %Get posture & target lists
    postureList = unique([trajStruct.posture]); numPostures = size(postureList,2);
    targetList = unique([trajStruct.target]); numTargets = size(targetList,2);
    numChannels = size(trajStruct(1).avgSmoothFR.traj,2); 
            
    %Form X
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
end
