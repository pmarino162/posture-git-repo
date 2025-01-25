clear; clc; clf; close all

%% Setup saveFig   
    saveFig = false;
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Paper\20220907\Error Reduction Figures';
    set(0, 'DefaultFigureRenderer', 'painters');

%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\orli.mat')
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat')
    pcmap = orli;
    tcmap = customRainbow;
    
%% Setup resultStruct
    resultStruct = struct('animal',[],'dataset',[],'shiftR2',[],'shiftError',[],'noShiftR2',[],'noShiftError',[],'selfR2',[],'selfError',[]);
    structInd = 1;
    
%% Parameters
   %manVarThreshold = 90;
   numManPCs = 10;
   numIterations = 20;

%% Run loop for each dataset
    task = 'bci';
     reachDatasetList = {'E20210706','E20210707','E20210708',...
        'N20190222','N20190226','N20190227','N20190228','N20190307'...
        'R20200221','R20200222'};
    
    bciDatasetList = {'E20200316','E20200317','E20200318','E20200319','N20171215','N20180221','R20201020','R20201021'};
    
    for datasetList = reachDatasetList%{ 'R20200221'}%bciDatasetList%{'E20210706','E20210707','E20210708','E20210709'}%reachDatasetList
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
            case {'N20171215','N20180221'}
                trialInclStates(1).trialName = {'Nigel Posture BC Center Out'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'state','Cursor Freeze','first',0},{'state','Target Hold','first',0}};
            case {'R20201020','R20201021'}
                trialInclStates(1).trialName = {'Rocky Posture BC Center Out'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'state','React','first',0},{'state','Hold','first',0}};
            %Reach
            case {'E20210706','E20210707','E20210708','E20210709','E20210710'}
                trialInclStates(1).trialName = {'GridReaching'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',0}};
            case {'R20200221','R20200222'}
                trialInclStates(1).trialName = {'Rocky Dissociation'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',0}};
            case {'N20190222','N20190226','N20190227','N20190228','N20190301','N20190305','N20190306','N20190307'}
                trialInclStates(1).trialName = {'Nigel Dissociation'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',0}};
        end
        
        %Remove any long or short trials - Reach
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


       trajStruct = getTrajStruct20220419(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'matchConditions',false,'zScoreParams',zScoreParams);
      
%         %Remove long/short trials, BCI
%         reachTime = [];
%         for i = 1:size(trajStruct,2)
%            numTraj = size(trajStruct(i).allSmoothFR,2);
%            for j = 1:numTraj
%               reachTime = [reachTime,size(trajStruct(i).allSmoothFR(j).traj,1)]; 
%            end
%         end     
%         figure
%             histogram(reachTime)
%             xlabel('Reach Time')
%             ylabel('Number of trials')
%        reachCutoff = prctile(reachTime,95);
%        for i = 1:size(trajStruct,2)
%            numTraj = size(trajStruct(i).allSmoothFR,2);
%            reachTime = [];
%            for j = 1:numTraj
%               reachTime = [reachTime,size(trajStruct(i).allSmoothFR(j).traj,1)]; 
%            end
%            rmTrials =  find(reachTime > reachCutoff);
%            trajStruct(i).allSmoothFR(rmTrials) = [];
%         end  
%        
 
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
        
        %Keep only postures with all targets for earl
        postureList = unique([trajStruct.posture]);
        targetList = unique([trajStruct.target]); 
        if strcmpi(dataset(1),'E')
            keepPosture = [];
            for posture = postureList
                tempTrajStruct = trajStruct([trajStruct.posture]==posture);
                postureTargetList = [tempTrajStruct.target];
                if isequal(postureTargetList,targetList)
                    keepPosture = [posture,keepPosture];
                end
            end
            trajStruct = trajStruct(ismember([trajStruct.posture],keepPosture));
        end
        
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
        
        %Do computation
        numSample = floor(minNumCondTraj/2); %Number of trajectories in each draw
        
%         selfError = NaN(1,2);% NaN(1,numIterations*numPostures*numTargets);
%         selfR2 = NaN(1,2);%NaN(1,numIterations*numPostures*numTargets);       
%         noShiftError = NaN(1,2);%NaN(1,numIterations*numTargets*nchoosek(numPostures,2));
%         noShiftR2 = NaN(1,2);%NaN(1,numIterations*numTargets*nchoosek(numPostures,2));
%         shiftError = NaN(1,2);%NaN(1,numIterations*numTargets*nchoosek(numPostures,2));
%         shiftR2 = NaN(1,2);%NaN(1,numIterations*numTargets*nchoosek(numPostures,2));
                
        selfError = NaN(1,numIterations*numPostures*numTargets);
        selfR2 = NaN(1,numIterations*numPostures*numTargets);       
        noShiftError = NaN(1,numIterations*numTargets*nchoosek(numPostures,2));
        noShiftR2 = NaN(1,numIterations*numTargets*nchoosek(numPostures,2));
        shiftError = NaN(1,numIterations*numTargets*nchoosek(numPostures,2));
        shiftR2 = NaN(1,numIterations*numTargets*nchoosek(numPostures,2));
        
        selfCompInd = 1;
        acrossPostInd = 1;
        for i = 1:numIterations
            %Split trajStruct into 2 groups
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
            
            %Visualize all trajectories - is raw data ok? 
%             visualizeAllTraj(trajStruct,pcmap)
%             visualizeAllTraj(trajStruct1,pcmap)
%             visualizeAllTraj(trajStruct2,pcmap)

            %Do self-comparisons
            for condInd = 1:numConditions
               posture = trajStruct(condInd).posture;
               target = trajStruct(condInd).target;
               traj1 = trajStruct1(condInd).avgSmoothFR.traj;
               traj2 = trajStruct2(condInd).avgSmoothFR.traj;
               %numPts = min([size(traj1,1),size(traj2,1)]);
               unshiftTraj2 = traj2;
               %alpha = getOptimalAlpha(traj1,traj2,numPts);
               if min([size(traj1,1),size(traj2,1)]) < numPts
                    tempNumPts = min([size(traj1,1),size(traj2,1)]);
               else
                    tempNumPts = numPts;
               end
               cvShift = getCVShift(posture,posture,target,tempNumPts,numDims,numTargets,trajStruct1,trajStruct2,targetList);
               traj2 = traj2 + cvShift;                 
               selfError(selfCompInd) = getMeanDist(traj1,traj2,tempNumPts);
               selfR2(selfCompInd) = getR2(traj1,traj2,tempNumPts);
               selfCompInd = selfCompInd + 1;
               %Visualize Individual Comparisons - is procedure
               %working?
%                    cond1 = condInd; cond2 = condInd;
%                    coeff = NaN;
%                    dist = getMeanDist(traj1,traj2,numPts);
%                    alpha = vecnorm(cvShift);
%                    visualizeTraj(traj1,traj2,unshiftTraj2,trajStruct,cond1,cond2,coeff,tcmap,numPts,dist,alpha);
%                    close
            end
            
            %Do across-posture comparisons (all from trajStruct2)
            predStruct = trajStruct;
            for target = targetList %for every target, make every posture comparison & shuffle comaparison
               %targetPostureList = [trajStruct([trajStruct.target]==target).posture];
               for posture1Ind = 1:numPostures-1
                   posture1 = postureList(posture1Ind);
                   cond1 = find([trajStruct.target]==target & [trajStruct.posture]==posture1);
                   traj1 = trajStruct2(cond1).avgSmoothFR.traj;
                   predStruct(cond1).avgSmoothFR.traj = NaN;
                   for posture2Ind = posture1Ind+1:numPostures  
                        posture2 = postureList(posture2Ind);                       
                        cond2 = find([trajStruct.target]==target & [trajStruct.posture]==posture2);                      
                        traj2 = trajStruct2(cond2).avgSmoothFR.traj;
                        %Get average error between trajectories
                        %numPts = min([size(traj1,1),size(traj2,1)]);
                        unshiftTraj2 = traj2;
                        if min([size(traj1,1),size(traj2,1)]) < numPts
                            tempNumPts = min([size(traj1,1),size(traj2,1)]);
                        else
                            tempNumPts = numPts;
                        end
                        
                        
                        noShiftError(acrossPostInd) = getMeanDist(traj1,traj2,tempNumPts);
                        noShiftR2(acrossPostInd) = getR2(traj1,traj2,tempNumPts);
                        %Subtract mean from group 2 & shuffle
%                         alpha = getOptimalAlpha(traj1,traj2,numPts);
%                         traj2 = traj2 + alpha;
                        cvShift = getCVShift(posture1,posture2,target,tempNumPts,numDims,numTargets,trajStruct1,trajStruct2,targetList);
                        traj2 = traj2 + cvShift;
                        %Get remaining average error
                        shiftError(acrossPostInd) = getMeanDist(traj1,traj2,tempNumPts);
                        shiftR2(acrossPostInd) = getR2(traj1,traj2,tempNumPts);
                        acrossPostInd = acrossPostInd + 1;
                        %Visualize Individual Comparisons - is procedure working?
%                        coeff = NaN;
%                        dist = getMeanDist(traj1,traj2,numPts);
%                        alpha = vecnorm(cvShift);
%                        visualizeTraj(traj1,traj2,unshiftTraj2,trajStruct,cond1,cond2,coeff,tcmap,numPts,dist,alpha);
%                        close;
                        predStruct(cond1).avgSmoothFR.traj = traj2;
                        predStruct(cond1).avgSmoothFR.timestamps =  trajStruct2(cond2).avgSmoothFR.timestamps;
                   end
               end
            end
        end

        %Add results to resultStruct
        resultStruct(structInd).dataset = dataset;
        resultStruct(structInd).animal = dataset(1);
        resultStruct(structInd).selfError = selfError;
        resultStruct(structInd).selfR2 = selfR2;
        resultStruct(structInd).noShiftError = noShiftError;
        resultStruct(structInd).noShiftR2 = noShiftR2;
        resultStruct(structInd).shiftError = shiftError;
        resultStruct(structInd).shiftR2 = shiftR2;
        structInd = structInd + 1;
        
    end
    
%% Plot results
    %Across datasets
    for i = 1:numel(resultStruct)
        resultMonkeyList{i} = resultStruct(i).animal;
    end
    monkeyList = unique(resultMonkeyList);
    numMonkeys = numel(monkeyList);
    
    %Plotting parameters
    fs = 12;
    alpha = 0.3;
    offset = 1/7;   
    sessionOffset = 1/40;
    
    %Plot R2 for each dataset
    f=figure; hold on;
    monkeyInd = 1;
    for monkey = monkeyList
       tempResultStruct = resultStruct(strcmp(resultMonkeyList,monkey{1,1}));
       for session = 1:numel(tempResultStruct)
           %Plot NoShiftR2
           noShiftR2 = horzcat(tempResultStruct(session).noShiftR2);
           noShiftR2mean = mean(noShiftR2);
           noShiftR2std = std(noShiftR2);
           e = errorbar(monkeyInd - offset + (session-1)*sessionOffset,noShiftR2mean,noShiftR2std,'ok','LineWidth',2,'MarkerSize',2,'MarkerFaceColor','r');      
           %Plot Shift R2
           shiftR2 = horzcat(tempResultStruct(session).shiftR2);
           shiftR2mean = mean(shiftR2);
           shiftR2std = std(shiftR2);
           e = errorbar(monkeyInd + (session-1)*sessionOffset,shiftR2mean,shiftR2std,'ok','LineWidth',2,'MarkerSize',2,'MarkerFaceColor','g');      
           %Plot Self R2
           selfR2 = horzcat(tempResultStruct(session).selfR2);
           selfR2mean = mean(selfR2);
           selfR2std = std(selfR2);
           e = errorbar(monkeyInd + offset + (session-1)*sessionOffset,selfR2mean,selfR2std,'ok','LineWidth',2,'MarkerSize',2,'MarkerFaceColor','b');      
       end
       monkeyInd = monkeyInd + 1; 
    end
 
    %Plot R2
    f=figure; hold on;
    monkeyInd = 1;
    xtickList = []; xtickListInd = 1;
    for monkey = monkeyList
       tempResultStruct = resultStruct(strcmp(resultMonkeyList,monkey{1,1}));
       %Plot NoShiftR2
       noShiftR2 = horzcat(tempResultStruct.noShiftR2);
       noShiftR2mean = mean(noShiftR2);
       noShiftR2std = std(noShiftR2);
       e = errorbar(monkeyInd - offset,noShiftR2mean,noShiftR2std,'ok','LineWidth',2,'MarkerSize',2,'MarkerFaceColor','r');      
       %Plot Shift R2
       shiftR2 = horzcat(tempResultStruct.shiftR2);
       shiftR2mean = mean(shiftR2);
       shiftR2std = std(shiftR2);
       e = errorbar(monkeyInd,shiftR2mean,shiftR2std,'ok','LineWidth',2,'MarkerSize',2,'MarkerFaceColor','g');      
       %Plot Self R2
       selfR2 = horzcat(tempResultStruct.selfR2);
       selfR2mean = mean(selfR2);
       selfR2std = std(selfR2);
       e = errorbar(monkeyInd + offset,selfR2mean,selfR2std,'ok','LineWidth',2,'MarkerSize',2,'MarkerFaceColor','b');      
       %Update xtickList
       xtickList(xtickListInd) = monkeyInd - offset;
       xtickList(xtickListInd+1) = monkeyInd;
       xtickList(xtickListInd+2) = monkeyInd + offset;
       xtickListInd = xtickListInd + 3;
       monkeyInd = monkeyInd + 1;
    end
    set(gca,'fontname','arial'); set(gca,'fontsize',fs)
    ylabel('R2')
    if saveFig
         saveas(gcf,fullfile(saveDir,[dataset,'_ErrorRed.svg']));
    end
    
    
    %Plot Error
    f=figure; hold on;
    monkeyInd = 1;
    xtickList = []; xtickListInd = 1;
    for monkey = monkeyList
       tempResultStruct = resultStruct(strcmp(resultMonkeyList,monkey{1,1}));
       %Plot NoShift Error
       noShiftError = horzcat(tempResultStruct.noShiftError);
       noShiftErrormean = mean(noShiftError);
       noShiftErrorstd = std(noShiftError);
       e = errorbar(monkeyInd - offset,noShiftErrormean,noShiftErrorstd,'ok','LineWidth',2,'MarkerSize',2,'MarkerFaceColor','r');      
       %Plot Shift Error
       shiftError = horzcat(tempResultStruct.shiftError);
       shiftErrormean = mean(shiftError);
       shiftErrorstd = std(shiftError);
       e = errorbar(monkeyInd,shiftErrormean,shiftErrorstd,'ok','LineWidth',2,'MarkerSize',2,'MarkerFaceColor','g');             
       %Plot Self Error
       selfError = horzcat(tempResultStruct.selfError);
       selfErrormean = mean(selfError);
       selfErrorstd = std(selfError);
       e = errorbar(monkeyInd + offset,selfErrormean,selfErrorstd,'ok','LineWidth',2,'MarkerSize',2,'MarkerFaceColor','b');      
       %Update xtickList
       xtickList(xtickListInd) = monkeyInd - offset;
       xtickList(xtickListInd+1) = monkeyInd;
       xtickList(xtickListInd+2) = monkeyInd + offset;
       xtickListInd = xtickListInd + 3;
       monkeyInd = monkeyInd + 1;
    end
    set(gca,'fontname','arial'); set(gca,'fontsize',fs)
    ylabel('Error (std)')
    
%% Local function   
    %Get optimal alpha
    function alpha = getOptimalAlpha(traj1,traj2,numPts)
        numDims = size(traj1,2);
        alpha = NaN(1,numDims);
        for dim = 1:numDims
            alpha(dim) = (1/numPts)*sum(traj1(1:numPts,dim)-traj2(1:numPts,dim));
        end
    end
    
    %Get mean dist
    function dist = getMeanDist(traj1,traj2,numPts)
        dist = 0;
        for i = 1:numPts
            dist = vecnorm(traj1(i,:)-traj2(i,:)) + dist;
        end
        dist = dist/numPts;
    end
    
    %Get R2
    function r2 = getR2(traj1,traj2,numPts)
        SST = 0;
        SSR = 0;
        traj1mean = mean(traj1(1:numPts,:),1);
        numDims = size(traj1,2);
        for i = 1:numPts
            for j = 1:numDims
                SST = (traj1(i,j)-traj1mean(j)).^2 + SST;
                SSR = (traj1(i,j)-traj2(i,j)).^2 + SSR;
            end
        end
        r2 = 1-SSR/SST;
    end
    
    %Get cvShift
    function cvShift = getCVShift(posture1,posture2,curTarget,numPts,numDims,numTargets,trajStruct1,trajStruct2,targetList)
        %Find all targets (other than curTarget) which have data from both
        %postures
        
        %Get optimal alpha for every target besides current one
        alphaList = NaN(numTargets-1,numDims);
        listInd = 1;
        for target = setdiff(targetList,curTarget)
            traj1 = trajStruct1([trajStruct1.posture]==posture1 & [trajStruct1.target]==target).avgSmoothFR.traj;
            traj2 = trajStruct2([trajStruct2.posture]==posture2 & [trajStruct2.target]==target).avgSmoothFR.traj;
            if min([size(traj1,1),size(traj2,1)]) < numPts
                tempNumPts = min([size(traj1,1),size(traj2,1)]);
            else
                tempNumPts = numPts;
            end        
            alphaList(listInd,:) = getOptimalAlpha(traj1,traj2,tempNumPts);
            listInd = listInd + 1;
        end
        %Average together
        cvShift = mean(alphaList,1);
    end
    
    %Visualize all trajectories in struct
    function [] = visualizeAllTraj(trajStruct,pcmap,numManPCs,allTraj)
        % Get timestamps dist and min number timestamps
        numTimestamps = [];
        for i = 1:size(trajStruct,2)
            numTimestamps = [numTimestamps,size(trajStruct(i).avgSmoothFR.timestamps,2)]; 
        end
        
        % Get minimum number of trials and timestamps
        [minNumTimestamps,~] = min(numTimestamps);

        %Get posture and target lists
        postureList = unique([trajStruct.posture]);
        numPostures = size(postureList,2);
        targetList = unique([trajStruct.target]); 
        numTargets = size(targetList,2);
        numChannels = size(trajStruct(1).avgSmoothFR.traj,2);
        
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
    
        %dPCA Plot
        Xdpca = permute(X,[4,3,2,1]);
        combinedParams = {{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}};
        margNames = {'P', 'T', 'CI', 'PTI'};        
        
        [W,V,whichMarg] = dpca(Xdpca, numManPCs, ...
            'combinedParams', combinedParams);
        explVar = dpca_explainedVariance(Xdpca, W, V, ...
            'combinedParams', combinedParams);        
        
        % Plot Marginal dims
        postureMargID = find(strcmp(margNames,'P'));
        targetMargID = find(strcmp(margNames,'T'));
        CIMargID = find(strcmp(margNames,'CI'));

        postureDims = find(whichMarg==postureMargID,2);
        targetDims = find(whichMarg==targetMargID,2);
        CIDims = find(whichMarg==CIMargID,2);
        dimList = [postureDims(1),targetDims(1),CIDims(1),NaN(1,3)];

        if numel(postureDims) > 1
            dimList(4) = postureDims(2);
        end
        if numel(targetDims) > 1
            dimList(5) = targetDims(2);
        end
        if numel(CIDims) > 1
            dimList(6) = CIDims(2);
        end

        figure;
        if numel(postureList)>5
            postureList = postureList(1:5);
        end
        
        for i = 1:6
            subplot(2,3,i); hold on;
            if ~isnan(dimList(i))
                postureInd = 1;
                for posture = postureList
                    for target = targetList
                        if any([trajStruct.posture]==posture & [trajStruct.target]==target)
                            traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgSmoothFR.traj; 
                            traj = traj*W(:,dimList(i));
                            time = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgSmoothFR.timestamps; 
                            plot(time,traj,'Color',pcmap(postureInd,:),'LineWidth',2);
                        end
                    end
                    postureInd = postureInd + 1;
                end
            end
        end     
        
        %Orthonormalize pdpca1 and tdpca1-2
        pDPCA = W(:,postureDims);
        tDPCA = W(:,targetDims);
        [PTOrth,~] = qr([pDPCA(:,1),tDPCA(:,1:2)]); PTOrth = PTOrth(:,1:3);
        
        %Add Projections to trajStruct
        %totalVar = trace(cov(allTraj));
        for i = 1:size(trajStruct,2)
            %Add average traces to trajStruct
            trajStruct(i).PTOrth.traj = trajStruct(i).avgSmoothFR.traj*PTOrth;
            %Get VAF
            %trajStruct(i).PTOrth.VAF =  [explVar.componentVar(postureDims(1)),explVar.componentVar(targetDims(1)),explVar.componentVar(targetDims(2))];%100.*(diag(cov(allTraj*PTOrth))')./totalVar;
        end
        %100.*(diag(cov(allTraj*allPCA))')./totalVar;

        %Plot - Orthographic
        fs = 14;
        figure
        hold on
        timePts = 1:minNumTimestamps;
        xDim = 2; yDim = 3; zDim = 1;
        postureInd = 1;
        for posture = postureList
            for target = targetList
                if any([trajStruct.posture]==posture & [trajStruct.target]==target)
                    traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).PTOrth.traj(timePts,:); 
                    plot3(traj(:,xDim),traj(:,yDim),traj(:,zDim),'Color',pcmap(postureInd,:),'LineWidth',2);
                    plot3(traj(1,xDim),traj(1,yDim),traj(1,zDim),'.','MarkerSize',20,'Color',pcmap(postureInd,:),'MarkerFaceColor',pcmap(postureInd,:));
                    plot3(traj(end,xDim),traj(end,yDim),traj(end,zDim),'<','MarkerSize',5,'Color',pcmap(postureInd,:),'MarkerFaceColor',pcmap(postureInd,:));
                end
            end
            postureInd = postureInd +1;
        end
             
        
        grid on
        xticklabels({}); yticklabels({}); zticklabels({}); 
        %VAF = round(trajStruct(1).PTOrth.VAF);

%         xlabel(['Goal Dim 1 (',num2str(VAF(xDim)),'%)'])
%         ylabel(['Goal Dim 2 (',num2str(VAF(yDim)),'%)']) 
%         zlabel(['Posture Dim 1 (',num2str(VAF(zDim)),'%)'])
        xlabel('Goal Dim 1')
        ylabel('Goal Dim 2') 
        zlabel('Posture Dim 1')
        view([20 10])
        set(gca,'fontname','arial')
        set(gca,'fontsize',fs)

    end
    
    %Visualize individual comparison
    function [] = visualizeTraj(traj1,traj2,unshiftTraj2,trajStruct,cond1,cond2,coeff,tcmap,numPts,dist,alpha)
        traj1Proj = traj1;%traj1*coeff;
        traj2Proj = traj2;%traj2*coeff;
        unshiftTraj2Proj = unshiftTraj2;%unshiftTraj2*coeff;
        traj1AvgProj = trajStruct(cond1).avgSmoothFR.traj;%trajStruct(cond1).PCA;
        traj2AvgProj = trajStruct(cond2).avgSmoothFR.traj;%trajStruct(cond2).PCA;
        
        figure
        hold on
        %Traj 1
        plot3(traj1Proj(:,1),traj1Proj(:,2),traj1Proj(:,3),'Color','r');
        for pt = 1:numPts
            ptInd = mod(pt,8);
            if ptInd == 0
                ptInd = 8;
            end
            if pt==1
                plot3(traj1Proj(pt,1),traj1Proj(pt,2),traj1Proj(pt,3),'.','MarkerSize',20,'Color',tcmap(ptInd,:));
            else
                plot3(traj1Proj(pt,1),traj1Proj(pt,2),traj1Proj(pt,3),'.','MarkerSize',10,'Color',tcmap(ptInd,:));
            end
        end
        %Unshifted Traj 2
        plot3(unshiftTraj2Proj(:,1),unshiftTraj2Proj(:,2),unshiftTraj2Proj(:,3),'Color','b');
        %Shifted Traj 2
        plot3(traj2Proj(:,1),traj2Proj(:,2),traj2Proj(:,3),'Color','g');
        for pt = 1:numPts
            ptInd = mod(pt,8);
            if ptInd == 0
                ptInd = 8;
            end
            plot3(traj2Proj(pt,1),traj2Proj(pt,2),traj2Proj(pt,3),'.','MarkerSize',10,'Color',tcmap(ptInd,:));
        end
        %Avg Traj 1
        plot3(traj1AvgProj(:,1),traj1AvgProj(:,2),traj1AvgProj(:,3),'Color','r','LineWidth',2);
        plot3(traj1AvgProj(1,1),traj1AvgProj(1,2),traj1AvgProj(1,3),'.','MarkerSize',10,'Color','r');
        %Avg Traj 2
        plot3(traj2AvgProj(:,1),traj2AvgProj(:,2),traj2AvgProj(:,3),'Color','b','LineWidth',2);
        plot3(traj2AvgProj(1,1),traj2AvgProj(1,2),traj2AvgProj(1,3),'.','MarkerSize',10,'Color','b');

        %Get condition info 
        cond1P = trajStruct(cond1).posture;
        cond1T = trajStruct(cond1).target;
        cond2P = trajStruct(cond2).posture;
        cond2T = trajStruct(cond2).target;
        
        xlabel('PC 1')
        ylabel('PC 2')
        zlabel('PC 3')
        title(['P',num2str(cond1P),'T',num2str(cond1T),' vs ','P',num2str(cond2P),'T',num2str(cond2T)...
            newline,'Dist = ',num2str(dist),'Hz',newline,'Alpha = ',num2str(alpha)])
    end