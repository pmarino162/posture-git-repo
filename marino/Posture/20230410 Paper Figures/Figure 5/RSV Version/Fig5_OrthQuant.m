clear; clc; clf; close all

%% Setup saveFig   
    saveFig = false;
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Paper\20230410\Figure 3\Projections';
    set(0, 'DefaultFigureRenderer', 'painters');
    
%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\orli.mat')
    pcmap = orli;
    
%% Parameters
   numManPCs = 10;     %Num PCs to project data into before analysis
   numPdims = 2;
   numBootReps = 10;   %Num bootstrap resamples
   numRandReps = 10; %Num random subspace draws per bootstrap
   numPts = 8;        %Num timePts per trajectory to include
   cutoffNumTraj = 10; %Num trials that must be present in a condition to keep it for analysis 
   
%% Setup resultStruct
    nanMat = NaN(1,2*numBootReps);
    resultStruct = struct('task',[],'g1pSpace',[],'g2pSpace',[],'g1numPts',[],'g2numPts',[],'vT1',nanMat,'vT2',nanMat,'vT3',nanMat,'vT4',nanMat,'vR',nanMat);
    resultStruct = repmat(resultStruct,1,3);
    structInd = 1;
    for task = 1:3
       resultStruct(task).task = task;      
    end
    
%% Get trajStruct
    %Load data
    dataset = 'E20200314';
    [Data,zScoreParams] = loadData(dataset);

    %Get trajStruct
    binWidth = 25; kernelStdDev = 25;
    trialInclStates = struct('trialName','','inclStates',[]);
    trajFields = {'zSmoothFR','markerVel','marker'};
    condFields = {{'task','conditionData','taskID'},{'target','targetData','targetID'},{'posture','conditionData','postureID'}};

    trialInclStates(1).trialName = {'GridTask_BC_ForceBar'};
        trialInclStates(1).inclStates = {{'state','Step 1 Freeze','first',0},{'state','Step 2','first',0}};
    trialInclStates(2).trialName = {'HC_CenterOut_ForceBar_20200314'};
        trialInclStates(2).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',100}};
    trialInclStates(3).trialName = {'IsometricForce_1D'};
        trialInclStates(3).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',100}};
    trajStruct = getTrajStruct20220419(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'matchConditions',false,'zScoreParams',zScoreParams);

    taskIDs = struct('ID',[],'task','');
    taskIDs(1).ID = 1; taskIDs(1).task = 'BC';
    taskIDs(2).ID = 2; taskIDs(2).task = 'HC';
    taskIDs(3).ID = 3; taskIDs(3).task = 'Iso';
    taskIDs(4).ID = 4; taskIDs(4).task = 'All';
    
%% Compute results    
    %Project all data down to top PCs
           % Get minimum number of timestamps in condition averages
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
        numPts = minNumTimestamps;

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
    
    %NEED TO UPDATE - ONLY 9 TRIALS FOR COND 31
    
    %Split data into two groups for cross-validaiton
    numCondTraj = [];
    for i = 1:size(trajStruct,2)
        numTraj = size(trajStruct(i).allSmoothFR,2);
        numCondTraj = [numCondTraj,numTraj];
    end
    [minNumCondTraj,i] = min(numCondTraj);
        
    numSample = 8; %floor(minNumCondTraj/2); %Number of trajectories in each draw
    trajStruct1 = trajStruct;
    trajStruct2 = trajStruct;
    for j = 1:size(trajStruct,2)
        numTraj = size(trajStruct(j).allSmoothFR,2);
        sampInd1 = randsample(numTraj,numSample);
        trajStruct1(j).allSmoothFR = trajStruct(j).allSmoothFR(sampInd1);
        [trajStruct1(j).avgSmoothFR.traj,trajStruct1(j).avgSmoothFR.timestamps] = getAvgTraj20211210(trajStruct1(j).allSmoothFR,binWidth);               
       
        if j == 31
            trajStruct2(j).allSmoothFR = trajStruct(j).allSmoothFR(sampInd1);
        else
            numTrajRemaining = numTraj - numSample;
            sampInd2 = randsample(numTrajRemaining,numSample);
            remainingInd = setdiff(1:numTraj,sampInd1);
            sampInd2 = remainingInd(sampInd2);
            trajStruct2(j).allSmoothFR = trajStruct(j).allSmoothFR(sampInd2);
        end
        
        [trajStruct2(j).avgSmoothFR.traj,trajStruct2(j).avgSmoothFR.timestamps] = getAvgTraj20211210(trajStruct2(j).allSmoothFR,binWidth);   
    end
    
    %For each group, get pSpace for all tasks
    for group = 1:2
        for task = 1:3
            if group == 1
                taskTrajStruct = trajStruct1([trajStruct1.task]==task);
            elseif group == 2
                taskTrajStruct = trajStruct2([trajStruct2.task]==task);
            end
            
           % Get minimum number of timestamps in condition averages
            numTimestamps = [];
            for i = 1:size(taskTrajStruct,2)
                numTimestamps(i) = length(taskTrajStruct(i).avgSmoothFR.timestamps);
            end
            [minNumTimestamps,i] = min(numTimestamps);

            if minNumTimestamps > 10
                numPts = 10;
            else
                numPts = minNumTimestamps;
            end

            %Get numPostures and numTargets
            postureList = unique([taskTrajStruct.posture]);
            numPostures = size(postureList,2);
            numChannels = size(taskTrajStruct(1).avgSmoothFR.traj,2);
            numConditions = size(taskTrajStruct,2);
                
            [pSig,~] = getPandTsig(taskTrajStruct,numPts);
            pSigReshape = reshape(squeeze(pSig),[numPostures*numPts,numChannels]);
            [pSpace,~,~,~,explainedP,pSigMu] = pca(pSigReshape); 
        
            if group == 1
                resultStruct(task).g1numPts = numPts;
                resultStruct(task).g1pSpace = pSpace;
            elseif group == 2
                resultStruct(task).g2numPts = numPts;
                resultStruct(task).g2pSpace = pSpace;
            end
        end
    end
    
    %For each group, bootstrap resample, get pSig for each task, then project into other task's space
    curInd = 1;
    for projGroup = 1:2
        if projGroup == 1
            projStruct = trajStruct1;
        elseif projGroup == 2
            projStruct = trajStruct2;
        end
        for bootRep = 1:numBootReps
            
            %Subsample w replacement 
            possibleSampInd = 1:numSample;
            sampStruct = projStruct;
            for j = 1:size(sampStruct,2)
                numTraj = size(projStruct(j).allSmoothFR,2);
                sampInd = randsample(numTraj,numTraj,true); %True => with replacement
                sampStruct(j).allSmoothFR = projStruct(j).allSmoothFR(sampInd);
                [sampStruct(j).avgSmoothFR.traj,sampStruct(j).avgSmoothFR.timestamps] = getAvgTraj20211210(sampStruct(j).allSmoothFR,binWidth);               
            end
            
            %For each task, get pSig, and compute projections for other spaces
            for task = 1:3
                taskTrajStruct = projStruct([projStruct.task]==task);
                
                if projGroup == 1
                    numPts = resultStruct(task).g1numPts;
                    DcurT = resultStruct(task).g2pSpace(:,1:numPdims);
                    DT1 = resultStruct(1).g2pSpace(:,1:numPdims);
                    DT2 = resultStruct(2).g2pSpace(:,1:numPdims);
                    DT3 = resultStruct(3).g2pSpace(:,1:numPdims);   
                else
                    numPts = resultStruct(task).g2numPts;
                    DcurT = resultStruct(task).g1pSpace(:,1:numPdims);
                    DT1 = resultStruct(1).g1pSpace(:,1:numPdims);
                    DT2 = resultStruct(2).g1pSpace(:,1:numPdims);
                    DT3 = resultStruct(3).g1pSpace(:,1:numPdims);  
                end
                
                [pSig,~] = getPandTsig(taskTrajStruct,numPts);
                pSigReshape = reshape(squeeze(pSig),[numPostures*numPts,numChannels]);
            
                CcurT = cov(pSigReshape); %Covariance of posture signal from current task
                
                vT1 = trace(DT1'*CcurT*DT1)./trace(CcurT);
                vT2 = trace(DT2'*CcurT*DT2)./trace(CcurT);
                vT3 = trace(DT3'*CcurT*DT3)./trace(CcurT);
                
                %Random spaces
                tempvR = NaN(1,numRandReps);
                for j = 1:numRandReps
                    v1 = normrnd(0,1,numManPCs,1);
                    v2 = normrnd(0,1,numManPCs,1);
                    DR = orth([v1,v2]);
                    tempvR(j) = trace(DR'*CcurT*DR)/trace(CcurT);
                end       
                
                %Fill resultStruct
                resultStruct(task).vT1(curInd) = vT1;
                resultStruct(task).vT2(curInd) = vT2;
                resultStruct(task).vT3(curInd) = vT3;
                resultStruct(task).vR(curInd) = mean(tempvR);
            end
            
            %Update curInd
            curInd = curInd + 1;
        end
    end
    
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
