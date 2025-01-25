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
    resultStruct = struct('animal',[],'dataset',[],'shiftR2',[],'shiftError',[],'noShiftR2',[],'noShiftError',[],'shufR2',[],'shufError',[],'shufShiftR2',[],'shufShiftError',[],'selfR2',[],'selfError',[]);
    structInd = 1;
        
%% Run loop for each dataset
    task = 'bci';
    for datasetList = {'E20200316','E20200317','N20171215','N20180221','R20201020'}
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
            case {'E20200316','E20200317','E20200318'}
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
        trajStruct = getTrajStruct20220419(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'matchConditions',false,'zScoreParams',zScoreParams);
        %Eliminate postures for which there are not trials to all targets
        switch dataset
            case {'E20210706'}
                trajStruct([trajStruct.posture]==5 | [trajStruct.posture]==6) = [];
        end
    
        % Get timestamps dist and min number timestamps
        numTimestamps = [];
        numCondTraj = [];
        for i = 1:size(trajStruct,2)
           numTraj = size(trajStruct(i).allSmoothFR,2);
           numCondTraj = [numCondTraj,numTraj];
           for j = 1:numTraj
              numTimestamps = [numTimestamps,size(trajStruct(i).allSmoothFR(j).timestamps,2)]; 
           end
        end
        figure
        histogram(numTimestamps)
        xlabel('Number of 25ms bins')
        ylabel('Number of trials')

        figure
        histogram(numCondTraj)
        xlabel('Number of trials')
        ylabel('Number of conditions')

        % Get minimum number of trials and timestamps
        [minNumTimestamps,i] = min(numTimestamps);
        [minNumCondTraj,i] = min(numCondTraj);

        %Get posture and target lists
        postureList = unique([trajStruct.posture]);
        numPostures = size(postureList,2);
        targetList = unique([trajStruct.target]); 
        numTargets = size(targetList,2);
        numChannels = size(trajStruct(1).avgSmoothFR.traj,2);
        numConditions = size(trajStruct,2);

        %Project all data down to top 10 PCs
        manVarThreshold = 90;
        %Form X
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

        %Do PCA; mean center; project down to manifold
        allTraj = reshape(X,[numTargets*numPostures*minNumTimestamps,numChannels]);
        [allPCs,~,~,~,explained,allMu] = pca(allTraj);
        numManPCs = 1;
        while sum(explained(1:numManPCs)) < manVarThreshold
            numManPCs = numManPCs + 1;
        end

        for i = 1:size(trajStruct,2)
           for j = 1:size(trajStruct(i).allSmoothFR,2)
                trajStruct(i).allSmoothFR(j).traj = (trajStruct(i).allSmoothFR(j).traj-allMu)*allPCs(:,1:numManPCs); 

           end
        end

        %Do computation
        %Parameters
        numIterations = 20;
        numSample = floor(minNumCondTraj/2); %Number of trajectories in each draw

        %Do self-comparisons
        selfError = NaN(1,numIterations*numPostures*numTargets);
        selfR2 = NaN(1,numIterations*numPostures*numTargets);
        index = 1;
        for target = targetList
           for posture = postureList
               cond1 = find([trajStruct.target]==target & [trajStruct.posture]==posture);
               numTraj1 = size(trajStruct(cond1).allSmoothFR,2);
               for i = 1:numIterations
                   %Split into even groups
                   sampInd1 = randsample(numTraj1,numSample);
                   traj1Struct = trajStruct(cond1).allSmoothFR(sampInd1);
                   numTrajRemaining = numTraj1 - numSample;
                   sampInd2 = randsample(numTrajRemaining,numSample);
                   remainingInd = setdiff(1:numTraj1,sampInd1);
                   sampInd2 = remainingInd(sampInd2);
                   traj2Struct = trajStruct(cond1).allSmoothFR(sampInd2);
                   %Get condition averages for each
                   traj1 = getAvgTraj20211210(traj1Struct,binWidth);
                   traj2 = getAvgTraj20211210(traj2Struct,binWidth);
                   %Get average error between trajectories
                   numPts = min([size(traj1,1),size(traj2,1)]);
                   alpha = getOptimalAlpha(traj1,traj2,numPts);
                   traj2 = traj2 + alpha;
                   selfError(index) = getMeanDist(traj1,traj2,numPts);
                   selfR2(index) = getR2(traj1,traj2,numPts);
                   index = index + 1;
               end
           end
        end

        %Get no shift error and shift error
        noShiftError = NaN(1,numIterations*numTargets*nchoosek(numPostures,2));
        noShiftR2 = NaN(1,numIterations*numTargets*nchoosek(numPostures,2));
        shiftError = NaN(1,numIterations*numTargets*nchoosek(numPostures,2));
        shiftR2 = NaN(1,numIterations*numTargets*nchoosek(numPostures,2));
        shufError = NaN(1,numIterations*numTargets*nchoosek(numPostures,2));
        shufR2 = NaN(1,numIterations*numTargets*nchoosek(numPostures,2));
        shufShiftError = NaN(1,numIterations*numTargets*nchoosek(numPostures,2));
        shufShiftR2 = NaN(1,numIterations*numTargets*nchoosek(numPostures,2));
        index = 1;
        for target = targetList %for every target, make every posture comparison & shuffle comaparison
           for posture1Ind = 1:numPostures-1
               posture1 = postureList(posture1Ind);
               for posture2Ind = posture1Ind+1:numPostures  
                   posture2 = postureList(posture2Ind);
                   cond1 = find([trajStruct.target]==target & [trajStruct.posture]==posture1);
                   numTraj1 = size(trajStruct(cond1).allSmoothFR,2);
                   cond2 = find([trajStruct.target]==target & [trajStruct.posture]==posture2);
                   numTraj2 = size(trajStruct(cond2).allSmoothFR,2);
                   condShuf = find([trajStruct.target]==targetList(randi(size(targetList,2))) & [trajStruct.posture]==posture2);
                   numTrajShuf = size(trajStruct(condShuf).allSmoothFR,2);
                   for i = 1:numIterations
                        %Select half the trials from each condition 
                        sampInd1 = randsample(numTraj1,numSample);
                        traj1Struct = trajStruct(cond1).allSmoothFR(sampInd1);
                        sampInd2 = randsample(numTraj2,numSample);
                        traj2Struct = trajStruct(cond2).allSmoothFR(sampInd2);
                        sampIndShuf = randsample(numTrajShuf,numSample);
                        trajShufStruct = trajStruct(condShuf).allSmoothFR(sampIndShuf);
                        %Form condition averages
                        traj1 = getAvgTraj20211210(traj1Struct,binWidth);
                        traj2 = getAvgTraj20211210(traj2Struct,binWidth);
                        trajShuf = getAvgTraj20211210(trajShufStruct,binWidth);
                        %Get average error between trajectories
                        numPts = min([size(traj1,1),size(traj2,1)]);
                        noShiftError(index) = getMeanDist(traj1,traj2,numPts);
                        noShiftR2(index) = getR2(traj1,traj2,numPts);
                        %Subtract mean from group 2 & shuffle
                        alpha = getOptimalAlpha(traj1,traj2,numPts);
                        traj2 = traj2 + alpha;
                        %Get remaining average error
                        shiftError(index) = getMeanDist(traj1,traj2,numPts);
                        shiftR2(index) = getR2(traj1,traj2,numPts);
                        numPtsShuf = min([size(traj1,1),size(trajShuf,1)]);
                        shufError(index) = getMeanDist(traj1,trajShuf,numPtsShuf);
                        shufR2(index) = getR2(traj1,trajShuf,numPtsShuf);
                        alphaShuf = getOptimalAlpha(traj1,trajShuf,numPtsShuf);
                        trajShuf = trajShuf + alphaShuf;
                        shufShiftError(index) = getMeanDist(traj1,trajShuf,numPtsShuf);
                        shufShiftR2(index) = getR2(traj1,trajShuf,numPtsShuf);
                        index = index + 1;
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
        resultStruct(structInd).shufError = shufError;
        resultStruct(structInd).shufR2 = shufR2;
        resultStruct(structInd).shufShiftError = shufShiftError;
        resultStruct(structInd).shufShiftR2 = shufShiftR2;
        structInd = structInd + 1;
    end
    
%% Plot results
    %Across datasets
    for i = 1:numel(resultStruct)
        resultMonkeyList{i} = resultStruct(i).animal;
    end
    monkeyList = unique(resultMonkeyList);
    numMonkeys = numel(monkeyList);

    fs = 12;
         alpha = 0.3;
         offset = (numMonkeys*2)/(5*numMonkeys-1);
    %Plot R2
    f=figure; hold on;
    monkeyInd = 1;
    xtickList = []; xtickListInd = 1;
    for monkey = monkeyList
       tempResultStruct = resultStruct(strcmp(resultMonkeyList,monkey{1,1}));

       %Plot ShuffleR2
       shufR2 = horzcat(tempResultStruct.shufR2);
       shufR2mean = mean(shufR2);
       shufR2std = std(shufR2);
       e = errorbar(monkeyInd*2-2*offset,shufR2mean,shufR2std,'ok','LineWidth',2,'MarkerSize',2,'MarkerFaceColor','r');      

       %Plot ShuffleShiftR2
       shufShiftR2 = horzcat(tempResultStruct.shufShiftR2);
       shufShiftR2mean = mean(shufShiftR2);
       shufShiftR2std = std(shufShiftR2);
       e = errorbar(monkeyInd*2-1*offset,shufShiftR2mean,shufShiftR2std,'ok','LineWidth',2,'MarkerSize',2,'MarkerFaceColor','r');      

       
       %Plot NoShiftR2
       noShiftR2 = horzcat(tempResultStruct.noShiftR2);
       noShiftR2mean = mean(noShiftR2);
       noShiftR2std = std(noShiftR2);
       e = errorbar(monkeyInd*2+0*offset,noShiftR2mean,noShiftR2std,'ok','LineWidth',2,'MarkerSize',2,'MarkerFaceColor','r');      

       %Plot Shift R2
       shiftR2 = horzcat(tempResultStruct.shiftR2);
       shiftR2mean = mean(shiftR2);
       shiftR2std = std(shiftR2);
       e = errorbar(monkeyInd*2+1*offset,shiftR2mean,shiftR2std,'ok','LineWidth',2,'MarkerSize',2,'MarkerFaceColor','g');      
       
       %Plot Self R2
       selfR2 = horzcat(tempResultStruct.selfR2);
       selfR2mean = mean(selfR2);
       selfR2std = std(selfR2);
       e = errorbar(monkeyInd*2+2*offset,selfR2mean,selfR2std,'ok','LineWidth',2,'MarkerSize',2,'MarkerFaceColor','b');      

       

       %Update xtickList
       xtickList(xtickListInd) = monkeyInd*2 - 2*offset;
       xtickList(xtickListInd+1) = monkeyInd*2 - 1*offset;
       xtickList(xtickListInd+2) = monkeyInd*2 + 0*offset;
       xtickList(xtickListInd+3) = monkeyInd*2 + 1*offset;
       xtickList(xtickListInd+4) = monkeyInd*2 + 2*offset;
       xtickListInd = xtickListInd + 5;
       monkeyInd = monkeyInd + 1;
    end
    
    set(gca,'fontname','arial'); set(gca,'fontsize',fs)
    ylabel('R2')

    if saveFig
         saveas(gcf,fullfile(saveDir,[dataset,'_ErrorRed.svg']));
    end
    
    
        %Plot R2
    f=figure; hold on;
    monkeyInd = 1;
    xtickList = []; xtickListInd = 1;
    for monkey = monkeyList
       tempResultStruct = resultStruct(strcmp(resultMonkeyList,monkey{1,1}));

       %Plot ShuffleR2
       shufError = horzcat(tempResultStruct.shufError);
       shufErrormean = mean(shufError);
       shufErrorstd = std(shufError);
       e = errorbar(monkeyInd*2-2*offset,shufErrormean,shufErrorstd,'ok','LineWidth',2,'MarkerSize',2,'MarkerFaceColor','r');      

              %Plot ShuffleShiftR2
       shufShiftError = horzcat(tempResultStruct.shufShiftError);
       shufShiftErrormean = mean(shufShiftError);
       shufShiftErrorstd = std(shufShiftError);
       e = errorbar(monkeyInd*2-1*offset,shufShiftErrormean,shufShiftErrorstd,'ok','LineWidth',2,'MarkerSize',2,'MarkerFaceColor','r');      

       
       %Plot NoShiftR2
       noShiftError = horzcat(tempResultStruct.noShiftError);
       noShiftErrormean = mean(noShiftError);
       noShiftErrorstd = std(noShiftError);
       e = errorbar(monkeyInd*2-0*offset,noShiftErrormean,noShiftErrorstd,'ok','LineWidth',2,'MarkerSize',2,'MarkerFaceColor','r');      

       %Plot Shift R2
       shiftError = horzcat(tempResultStruct.shiftError);
       shiftErrormean = mean(shiftError);
       shiftErrorstd = std(shiftError);
       e = errorbar(monkeyInd*2+1*offset,shiftErrormean,shiftErrorstd,'ok','LineWidth',2,'MarkerSize',2,'MarkerFaceColor','g');      
       
       %Plot Self R2
       selfError = horzcat(tempResultStruct.selfError);
       selfErrormean = mean(selfError);
       selfErrorstd = std(selfError);
       e = errorbar(monkeyInd*2+2*offset,selfErrormean,selfErrorstd,'ok','LineWidth',2,'MarkerSize',2,'MarkerFaceColor','b');      


       %Update xtickList
       xtickList(xtickListInd) = monkeyInd*2 - 2*offset;
       xtickList(xtickListInd+1) = monkeyInd*2 - 1*offset;
       xtickList(xtickListInd+2) = monkeyInd*2 + 0*offset;
       xtickList(xtickListInd+3) = monkeyInd*2 + 1*offset;
       xtickList(xtickListInd+4) = monkeyInd*2 + 2*offset;
       xtickListInd = xtickListInd + 5;
       monkeyInd = monkeyInd + 1;
    end
    
    set(gca,'fontname','arial'); set(gca,'fontsize',fs)
    ylabel('Error')
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