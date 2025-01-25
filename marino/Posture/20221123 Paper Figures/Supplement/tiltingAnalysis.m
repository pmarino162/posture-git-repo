clear; clc; clf; close all

%% Setup saveFig   
    saveFig = true;
    saveDir = 'D:\Posture Paper\20221123\Figures\Supplement\Tilting Analysis';
    set(0, 'DefaultFigureRenderer', 'painters');
    
%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\orli.mat')
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat')
    pcmap = vertcat(orli,customRainbow(4:5,:));
    tcmap = customRainbow;
    
%% Run Loop    
    for datasetList = {'E20210706'}%{'R20201020','E20200116','E20210706','N20171215','R20201020','N20190226'} % {'N20190226','R20200221'}
       clf; close all
        %{'E20200317','E20200116','E20210706','N20171215','R20201020','N20190226','R20200221'}
        % Load Data
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
            case 'E20210901'
                taskID = [Data.conditionData]; taskID = [taskID.taskID];
                Data = Data(taskID==1);
                trialInclStates(1).trialName = {'BCI Center Out'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'state','Step 1','first',0},{'state','Success with Reward','first',0}};
            %Iso
            case 'E20200116'
                trialInclStates(1).trialName = {'IsometricForce_1D'};   
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'state','Target','first',0},{'state','Target','first',250}};
            case 'E20211007'
                trialInclStates(1).trialName = {'IsometricForce_1D'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'state','Target','first',0},{'state','Target Hold','first',0}};
            %Reaching
            case 'E20210706'
                trialInclStates(1).trialName = {'GridReaching'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',50}};
            case{'N20190222','N20190226','N20190227','N20190228','N20190307'}
                trialInclStates(1).trialName = {'Nigel Dissociation'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',50}};
            case{'R20200221','R20200222'}
                trialInclStates(1).trialName = {'Rocky Dissociation'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',50}};
        end
        trajStruct = getTrajStruct20220419(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'matchConditions',false,'zScoreParams',zScoreParams);
         
        %Get posture and target lists
        postureList = unique([trajStruct.posture]);
        numPostures = size(postureList,2);
        targetList = unique([trajStruct.target]); 
        numTargets = size(targetList,2);
        numChannels = size(trajStruct(1).avgSmoothFR.traj,2);
        numConditions = size(trajStruct,2);
        
        %Keep only postures with every target for Earl
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
        
                
        %Update posture and list
        postureList = unique([trajStruct.posture]);
        numPostures = size(postureList,2);
        
        % Get minimum number of timestamps in condition averages
        numTimestamps = [];
        for i = 1:size(trajStruct,2)
            numTimestamps(i) = length(trajStruct(i).avgSmoothFR.timestamps);
        end
        [minNumTimestamps,i] = min(numTimestamps);


           
        %Get number of PCs required to capture 70% of variance
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
        
        varThreshold = 70;
        allTraj = reshape(X,[numTargets*numPostures*minNumTimestamps,numChannels]);
        [allPCs,~,~,~,explained,allMu] = pca(allTraj);
        numPCs = 1;
        while sum(explained(1:numPCs)) < varThreshold
            numPCs = numPCs + 1;
        end
        
        %Run cross-validated analysis
       numIterations = 20;
       result = zeros(numPostures); 
       posture1Ind = 1;
       for posture1 = postureList
            %Get posture 1 traj struct
            posture1TrajStruct = trajStruct([trajStruct.posture]==posture1);
            posture2Ind = 1;
            for posture2 = postureList
                %Get posture 2 traj struct
                posture2TrajStruct = trajStruct([trajStruct.posture]==posture2);
                tempResult = zeros(1,numIterations);
                for i = 1:numIterations
                    %Build up X1comp, X2comp, X1ref, X2ref
                    X1comp = []; X2comp = []; X1ref = []; X2ref = [];
                    for target = targetList
                        %Posture 1 (ref)
                        posture1TargetStruct = posture1TrajStruct([posture1TrajStruct.target]==target); 
                        numP1Trials = size(posture1TargetStruct.allSmoothFR,2);
                        sampInd1 = randsample(numP1Trials,floor(numP1Trials/2))';
                        sampInd2 = setdiff(1:numP1Trials,sampInd1);
                        traj1Struct = posture1TargetStruct.allSmoothFR(sampInd1);
                        traj1 = getAvgTraj20211210(traj1Struct,binWidth);
                        traj2Struct = posture1TargetStruct.allSmoothFR(sampInd2);
                        traj2 = getAvgTraj20211210(traj2Struct,binWidth);
                        X1ref = vertcat(X1ref,traj1);
                        X2ref = vertcat(X2ref,traj2);
                        
                        %Posture 2 (comp)
                        posture2TargetStruct = posture2TrajStruct([posture2TrajStruct.target]==target); 
                        numP2Trials = size(posture2TargetStruct.allSmoothFR,2);
                        sampInd1 = randsample(numP2Trials,floor(numP2Trials/2))';
                        sampInd2 = setdiff(1:numP2Trials,sampInd1);
                        traj1Struct = posture2TargetStruct.allSmoothFR(sampInd1);
                        traj1 = getAvgTraj20211210(traj1Struct,binWidth);
                        traj2Struct = posture2TargetStruct.allSmoothFR(sampInd2);
                        traj2 = getAvgTraj20211210(traj2Struct,binWidth);
                        X1comp = vertcat(X1comp,traj1);
                        X2comp = vertcat(X2comp,traj2);
                    end
                    %Get W
                    [W1ref,~,~,~,explained,~] = pca(X1ref);
                    [W2ref,~,~,~,explained,~] = pca(X2ref);
                    [W1comp,~,~,~,explained,~] = pca(X1comp);
                    [W2comp,~,~,~,explained,~] = pca(X2comp);
                    %Get result for current iteration
                    tempResult(i) = 0.5*(tiltVarExpl(X1comp,W2ref,numPCs)/tiltVarExpl(X1comp,W2comp,numPCs)...
                        +tiltVarExpl(X2comp,W1ref,numPCs)/tiltVarExpl(X2comp,W1comp,numPCs));        
                end
                result(posture1Ind,posture2Ind) = mean(tempResult);
                posture2Ind = posture2Ind + 1;
            end    
            posture1Ind = posture1Ind + 1;
       end
       
       for posture1 = 1:numPostures
               result(posture1,posture1) = 1;
       end
%        posture1Ind = 1;
%        for posture1 = postureList
%           %Fit Posture 1 PCs
%           posture1TrajStruct = trajStruct([trajStruct.posture]==posture1);
%           allPosture1Traj = vertcat(posture1TrajStruct.avgSmoothFR);
%           allPosture1Traj = vertcat(allPosture1Traj.traj);
%           [posture1PCs,~,~,~,explained,allMu] = pca(allPosture1Traj);
%           posture1PCs = posture1PCs(:,1:5);
%           posture2Ind = 1;
%           for posture2 = postureList
%               %Get Var Expl
%               posture2TrajStruct = trajStruct([trajStruct.posture]==posture2);
%               allPosture2Traj = vertcat(posture2TrajStruct.avgSmoothFR);
%               allPosture2Traj = vertcat(allPosture2Traj.traj);
%               cP2 = cov(allPosture2Traj);
%               varExpl = trace((posture1PCs'*cP2*posture1PCs))/trace(cP2);
%               result(posture1Ind,posture2Ind) = varExpl;
%               posture2Ind = posture2Ind + 1;
%           end    
%           posture1Ind = posture1Ind + 1;
%        end
        
       figure
       h = heatmap(result.*100);
       h.FontName = 'Arial'; h.FontSize = 14;
       ylabel('Posture Used for PCA')
       xlabel('Projected Posture')
        annotation('textarrow',[1,1],[.6,0.5],'string','Variance Captured (%)', ...
              'HeadStyle','none','LineStyle','none','HorizontalAlignment','center','TextRotation',90,'FontName','Arial');
        if saveFig
            saveas(gcf,fullfile(saveDir,[dataset,'tilting.svg']));
        end
        
        
    end
    
%% Local function for computing variance captured
    function [varExpl] = tiltVarExpl(X,W,numPCs)
        cX = cov(X);
        W = W(:,1:numPCs);
        varExpl = trace((W'*cX*W))/trace(cX);
    end