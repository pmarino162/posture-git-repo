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
    resultStruct = struct('animal',[],'dataset',[],'vPP',[],'vPT',[],'vTT',[],'vTP',[],'vPR',[],'vTR',[],'pAngle',[],'pAngleNull',[]);
    structInd = 1;
    
%% Parameters
   numManPCs = 10;     %Num PCs to project data into before analysis
   numPdims = 2;
   numTdims = 2;
   numBootReps = 10;   %Num bootstrap resamples
   numRandReps = 10; %Num random subspace draws per bootstrap
   numPts = 8;        %Num timePts per trajectory to include
   cutoffNumTraj = 10; %Num trials that must be present in a condition to keep it for analysis 
   
%% Run loop for each dataset    
    reachDatasetList = {'E20210706','E20210707','E20210708',...
        'N20190222','N20190226','N20190227','N20190228','N20190307'...
        'R20200221','R20200222'};
    
    bciDatasetList = {'E20200316','E20200317','E20200318','E20200319','N20171215','N20180221','R20201020','R20201021'};
           %NEED TO FIX 3/18 %{'E20200316','E20200317','E20200318','E20200319','N20171215','N20180221','R20201020','R20201021'};
    isoDatasetList = {'E20200116','E20200117','E20200120'};
    for datasetList = bciDatasetList
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
            %Iso
            case {'E20200116','E20200117','E20200120'}
                trialInclStates(1).trialName = {'IsometricForce_1D'};   
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'state','Target','first',0},{'state','Target','first',200}};
                task = 'iso';
            %Reach
            case {'E20210706','E20210707','E20210708','E20210709','E20210710'}
                trialInclStates(1).trialName = {'GridReaching'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',100}};
                task = 'reach';
            case {'R20200221','R20200222'}
                trialInclStates(1).trialName = {'Rocky Dissociation'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',100}};
                task = 'reach';
            case {'N20190222','N20190226','N20190227','N20190228','N20190301','N20190305','N20190306','N20190307'}
                trialInclStates(1).trialName = {'Nigel Dissociation'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',100}};
                task = 'reach';
        end 
        
       

       trajStruct = getTrajStruct20220419(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'matchConditions',false,'zScoreParams',zScoreParams);
      
        
         
        % Get number of trials for each condition, numTimestamps in each trial
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
        numCondTraj = [];
        for i = 1:size(trajStruct,2)
           numTraj = size(trajStruct(i).allSmoothFR,2);
           numCondTraj = [numCondTraj,numTraj];
        end
        [minNumCondTraj,i] = min(numCondTraj);

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

        %Preallocate 
        vPP = NaN(1,2*numBootReps); vPT = NaN(1,2*numBootReps); vPR = NaN(1,2*numBootReps);         
        vTT = NaN(1,2*numBootReps); vTP = NaN(1,2*numBootReps); vTR = NaN(1,2*numBootReps);
        pAngle = NaN(1,2*numBootReps);
 
        %Split data into two groups for cross-validation
        numSample = floor(minNumCondTraj/2); %Number of trajectories in each draw
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
        
        %Fit P and T spaces for each group
        [pSig1,tSig1] = getPandTsig(trajStruct1,numPts);
            pSig1Reshape = reshape(squeeze(pSig1),[numPostures*numPts,numManPCs]);
            tSig1Reshape = reshape(squeeze(tSig1),[numTargets*numPts,numManPCs]);
            [pDimsGroup1,~,~,~,explainedP,pSigMu] = pca(pSig1Reshape); 
            [tDimsGroup1,~,~,~,explainedT,tSigMu] = pca(tSig1Reshape); 
        [pSig2,tSig2] = getPandTsig(trajStruct2,numPts);
            pSig2Reshape = reshape(squeeze(pSig2),[numPostures*numPts,numManPCs]);
            tSig2Reshape = reshape(squeeze(tSig2),[numTargets*numPts,numManPCs]);
            [pDimsGroup2,~,~,~,~,~] = pca(pSig2Reshape); 
            [tDimsGroup2,~,~,~,~,~] = pca(tSig2Reshape); 
            
        %For each group, bootstrap resample, get P & T sig, then project into other group's space
        curInd = 1; %Index of projection vectors (e.g., vPP)
        for projGroup = 1:2
            if projGroup == 1
                projStruct = trajStruct1;
                pDims = pDimsGroup2;
                tDims = tDimsGroup2;
            elseif projGroup == 2
                projStruct = trajStruct2;
                pDims = pDimsGroup1;
                tDims = tDimsGroup1;
            end
            for bootRep = 1:numBootReps
                %Subsample with replacement
                possibleSampInd = 1:numSample;
                sampStruct = projStruct;
                for j = 1:size(sampStruct,2)
                    numTraj = size(projStruct(j).allSmoothFR,2);
                    sampInd = randsample(numTraj,numTraj,true); %True => with replacement
                    sampStruct(j).allSmoothFR = projStruct(j).allSmoothFR(sampInd);
                    [sampStruct(j).avgSmoothFR.traj,sampStruct(j).avgSmoothFR.timestamps] = getAvgTraj20211210(sampStruct(j).allSmoothFR,binWidth);               
                end
                %Get P and T signals
                [pSig,tSig] = getPandTsig(sampStruct,numPts);
                pSigReshape = reshape(squeeze(pSig),[numPostures*numPts,numManPCs]);
                tSigReshape = reshape(squeeze(tSig),[numTargets*numPts,numManPCs]);
                %Compute P and T projections
                CP = cov(pSigReshape);
                DP = pDims(:,1:numPdims);
                CT = cov(tSigReshape);
                DT = tDims(:,1:numTdims);
                vPP(curInd) = trace(DP'*CP*DP)./trace(CP);
                vPT(curInd) = trace(DT'*CP*DT)./trace(CP);
                vTT(curInd) = trace(DT'*CT*DT)./trace(CT);
                vTP(curInd) = trace(DP'*CT*DP)./trace(CT);
                
                %Get pAngle
                allPAngle = rad2deg(subspacea(DP,DT));
                %pAngle(curInd) = allPAngle(1);
                
                %Get pAngleNull
                v11 = normrnd(0,1,numManPCs,1);
                v21 = normrnd(0,1,numManPCs,1);
                DR1 = orth([v11,v21]);
                v12 = normrnd(0,1,numManPCs,1);
                v22 = normrnd(0,1,numManPCs,1);
                DR2 = orth([v12,v22]);
                %pAngleNull(curInd) = prinangle(DR1,DR2);
                
                %Compute random projections
                tempvPR = NaN(1,numRandReps);
                tempvTR = NaN(1,numRandReps);
                for j = 1:numRandReps
                    v1 = normrnd(0,1,numManPCs,1);
                    v2 = normrnd(0,1,numManPCs,1);
                    DR = orth([v1,v2]);
                    tempvPR(j) = trace(DR'*CP*DR)/trace(CP);
                    tempvTR(j) = trace(DR'*CT*DR)/trace(CT);

                end       
                vPR(curInd) = mean(tempvPR);
                vTR(curInd) = mean(tempvTR);
                
                %Update curInd
                curInd = curInd + 1;
            end
        end
        
        %Add results to resultStruct
        resultStruct(structInd).dataset = dataset; resultStruct(structInd).animal = dataset(1);
        resultStruct(structInd).vPP = vPP; resultStruct(structInd).vPT = vPT; resultStruct(structInd).vPR = vPR;
        resultStruct(structInd).vTT = vTT; resultStruct(structInd).vTP = vTP; resultStruct(structInd).vTR = vTR;
      %  resultStruct(structInd).pAngle = pAngle;
%        resultStruct(structInd).pAngleNull = pAngleNull;
        structInd = structInd + 1;
        
    end

    
%% Get monkey list 
    %Across datasets
    for i = 1:numel(resultStruct)
        resultMonkeyList{i} = resultStruct(i).animal;
    end
    monkeyList = unique(resultMonkeyList);
    numMonkeys = numel(monkeyList);
    
%% For each monkey, get overall distributions, CI's, t-tests
    tempStruct = struct();
    tStruct = struct('h',[],'p',[],'ci',[],'stats',tempStruct);
    plotStruct = struct('monkey',[],'vPP',[],'vPPCI',[],'vPT',[],'vPTCI',[],'vTT',[],'vTTCI',[],'vTP',[],...
        'vTPCI',[],'vPR',[],'vPRCI',[],'vTR',[],'vTRCI',[],'vPP_PTt',tStruct,'vPT_PRt',tStruct,'vTT_TPt',tStruct,'vTP_TRt',tStruct);
    CIpct = 95;
    edgeSize = (100-CIpct)/2;
    plotStructInd = 1;
    for monkey = monkeyList
        tempResultStruct = resultStruct(strcmp(resultMonkeyList,monkey{1,1}));
        
        vPP = [tempResultStruct.vPP]*100;
        vPPCI = [prctile(vPP,edgeSize) prctile(vPP,100-edgeSize)];
        vPT = [tempResultStruct.vPT]*100;
        vPTCI = [prctile(vPT,edgeSize) prctile(vPT,100-edgeSize)];
        vPR = [tempResultStruct.vPR]*100;
        vPRCI = [prctile(vPR,edgeSize) prctile(vPR,100-edgeSize)];
        vTT = [tempResultStruct.vTT]*100;
        vTTCI = [prctile(vTT,edgeSize) prctile(vTT,100-edgeSize)];
        vTP = [tempResultStruct.vTP]*100;
        vTPCI = [prctile(vTP,edgeSize) prctile(vTP,100-edgeSize)];
        vTR = [tempResultStruct.vTR]*100;
        vTRCI = [prctile(vTR,edgeSize) prctile(vTR,100-edgeSize)];
        
        plotStruct(plotStructInd).monkey = monkey{1,1};
        plotStruct(plotStructInd).vPP = vPP;
        plotStruct(plotStructInd).vPPCI = vPPCI;
        plotStruct(plotStructInd).vPT = vPT;
        plotStruct(plotStructInd).vPTCI = vPTCI;
        plotStruct(plotStructInd).vPR = vPR;
        plotStruct(plotStructInd).vPRCI = vPRCI;
        plotStruct(plotStructInd).vTT = vTT;
        plotStruct(plotStructInd).vTTCI = vTTCI;
        plotStruct(plotStructInd).vTP = vTP;
        plotStruct(plotStructInd).vTPCI = vTPCI;
        plotStruct(plotStructInd).vTR = vTR;
        plotStruct(plotStructInd).vTRCI = vTRCI;
        
        [h,p,ci,stats] = ttest2(vPP,vPT,'Vartype','unequal');
            plotStruct(plotStructInd).vPP_PTt.h = h;
            plotStruct(plotStructInd).vPP_PTt.p = p;
            plotStruct(plotStructInd).vPP_PTt.ci = ci;
            plotStruct(plotStructInd).vPP_PTt.stats = stats;

        [h,p,ci,stats] = ttest2(vPT,vPR,'Vartype','unequal');
            plotStruct(plotStructInd).vPT_PRt.h = h;
            plotStruct(plotStructInd).vPT_PRt.p = p;
            plotStruct(plotStructInd).vPT_PRt.ci = ci;
            plotStruct(plotStructInd).vPT_PRt.stats = stats;
            
        [h,p,ci,stats] = ttest2(vTT,vTP,'Vartype','unequal');
            plotStruct(plotStructInd).vTT_TPt.h = h;
            plotStruct(plotStructInd).vTT_TPt.p = p;
            plotStruct(plotStructInd).vTT_TPt.ci = ci;
            plotStruct(plotStructInd).vTT_TPt.stats = stats;
            
        [h,p,ci,stats] = ttest2(vTP,vTR,'Vartype','unequal');
            plotStruct(plotStructInd).vTP_TRt.h = h;
            plotStruct(plotStructInd).vTP_TRt.p = p;
            plotStruct(plotStructInd).vTP_TRt.ci = ci;
            plotStruct(plotStructInd).vTP_TRt.stats = stats;

        plotStructInd = plotStructInd + 1;
    end
    
%% Plot results
    %mcmap = [0.4*ones(1,3); 0.6*ones(1,3); 0.9*ones(1,3)];
    mcmap = summer(4);
    mcmap = mcmap(1:3,:);
    
    %Plotting parameters
    fs = 14;
    offset = 5;   
    
    for i = 1:numel(plotStruct)
        plotMonkeyList{i} = plotStruct(i).monkey;
    end
    
    %Bar plot for posture signal
    f=figure; hold on;
    monkeyInd = 1;
    for monkey = monkeyList
       tempPlotStruct = plotStruct(strcmp(plotMonkeyList,monkey{1,1}));
       vPPMean = mean(tempPlotStruct.vPP);
       vPPCI = tempPlotStruct.vPPCI;
       vPTMean = mean(tempPlotStruct.vPT);
       vPTCI = tempPlotStruct.vPTCI;
       
      bar(monkeyInd,vPPMean,'FaceColor',mcmap(monkeyInd,:),'EdgeColor',[0 0 0]);
      hold on;
      errorbar(monkeyInd,vPPMean,vPPMean-vPPCI(1),vPPCI(2)-vPPMean,'k','LineWidth',2);
      bar(monkeyInd+offset,vPTMean,'FaceColor',mcmap(monkeyInd,:),'EdgeColor',[0 0 0]);
      errorbar(monkeyInd+offset,vPTMean,vPTMean-vPTCI(1),vPTCI(2)-vPTMean,'k','LineWidth',2);
       monkeyInd = monkeyInd + 1; 
    end
    
    vPR = [resultStruct.vPR].*100;
    vPRMean = mean(vPR);
    vPRCI = [prctile(vPR,edgeSize) prctile(vPR,100-edgeSize)];
    ax = gca;
    curXLims = ax.XLim;
    shadedErrorBar(curXLims,[vPRMean,vPRMean],[vPRMean-vPRCI(1) vPRMean-vPRCI(1); vPRCI(2)-vPRMean vPRCI(2)-vPRMean],'lineprops',{'--','LineWidth',1.5,'Color',[0.3 0.3 0.3]});
    set(gca,'fontname','arial'); set(gca,'fontsize',fs);
    ax.TickDir = 'out';
    yticks([0:25:100]);
    xticks([]);
    ylabel('Variance Explained %');
    xlabel('Subspace');
    
    %Bar plot for goal signal
    f=figure; hold on;
    monkeyInd = 1;
    for monkey = monkeyList
       tempPlotStruct = plotStruct(strcmp(plotMonkeyList,monkey{1,1}));
       vTTMean = mean(tempPlotStruct.vTT);
       vTTCI = tempPlotStruct.vTTCI;
       vTPMean = mean(tempPlotStruct.vTP);
       vTPCI = tempPlotStruct.vTPCI;
       
      bar(monkeyInd,vTTMean,'FaceColor',mcmap(monkeyInd,:),'EdgeColor',[0 0 0]);
      hold on;
      errorbar(monkeyInd,vTTMean,vTTMean-vTTCI(1),vTTCI(2)-vTTMean,'k','LineWidth',2);
      bar(monkeyInd+offset,vTPMean,'FaceColor',mcmap(monkeyInd,:),'EdgeColor',[0 0 0]);
      errorbar(monkeyInd+offset,vTPMean,vTPMean-vTPCI(1),vTPCI(2)-vTPMean,'k','LineWidth',2);
       monkeyInd = monkeyInd + 1; 
    end
    
    vTR = [resultStruct.vTR].*100;
    vTRMean = mean(vTR);
    vTRCI = [prctile(vTR,edgeSize) prctile(vTR,100-edgeSize)];
    ax = gca;
    curXLims = ax.XLim;
    shadedErrorBar(curXLims,[vTRMean,vTRMean],[vTRMean-vTRCI(1) vTRMean-vTRCI(1); vTRCI(2)-vTRMean vTRCI(2)-vTRMean],'lineprops',{'--','LineWidth',1.5,'Color',[0.3 0.3 0.3]});
    set(gca,'fontname','arial'); set(gca,'fontsize',fs);
    ax.TickDir = 'out';
    yticks([0:25:100]);
    xticks([]);
    ylabel('Variance Explained %');
    xlabel('Subspace');

%% Local functions
    function [pSig,tSig] = getPandTsig(trajStruct,numPts)
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
        %Posture and Target Traj
        tSig = mean(X,[3],'omitnan') - CIMargTraj  - CIMargOffset;
        pSig = mean(X,[2],'omitnan') - CIMargTraj - CIMargOffset;
    end

