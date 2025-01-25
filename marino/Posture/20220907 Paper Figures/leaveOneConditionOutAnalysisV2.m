function [resultStruct] = leaveOneConditionOutAnalysisV2(dataset,task)

%% Set up results struct
    resultStruct = struct('dataset','','task','','noModelDist',[],'CIDist',[],...
        'CIPDist',[],'fullDist',[],'actualTrajStruct',[],'predTrajStruct',[],'R2',[]);
    
%% Get trajStruct
    [Data,zScoreParams] = loadData(dataset);
    binWidth = 25; kernelStdDev = 25;
    trajFields = {'zSmoothFR'};
    trialInclStates = struct('trialName','','inclStates',[]);
    switch dataset
        %Earl BCI
        case {'E20200316','E20200317','E20200318','E20200319'}
            trialInclStates(1).trialName = {'GridTask_CO_Across_BC_ForceBar'};
            condFields = {{'target','targetData','target1ID'},{'posture','conditionData','postureID'}};
            trialInclStates(1).inclStates = {{'state','Step 1','first',25},{'state','Step 1','first',225}};
        %Earl Reach
        case {'E20210706','E20210707','E20210708','E20210709','E20210710'}
            trialInclStates(1).trialName = {'GridReaching'};
            condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
            trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',0}};
        %Earl Multiple tasks
        case {'E20200314'}
            allData = Data;
            taskID = [Data.conditionData]; taskID = [taskID.taskID];
            condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
            if strcmpi(task,'BCI')
                Data = allData(taskID==1);
                trialInclStates(1).trialName = {'GridTask_BC_ForceBar'};
                trialInclStates(1).inclStates = {{'state','Step 1 Freeze','first',25},{'state','Step 2','first',0}};
            elseif strcmpi(task,'Reach')
                Data = allData(taskID==2);
                trialInclStates(1).trialName = {'HC_CenterOut_ForceBar_20200314'};
                trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'state','moveOnsetTime','first',0}};
            elseif strcmpi(task,'Iso')
                Data = allData(taskID==3);
                trialInclStates(1).trialName = {'IsometricForce_1D'};
                trialInclStates(1).inclStates = {{'state','Target','first',25},{'state','Target Hold','first',0}};
            end
        %Nigel and Rocky Reaching
        case {'N20190222','N20190226','N20190227','N20190228','N20190301','N20190305',...
                'N20190306','N20190307','R20200221','R20200222'}
            if strcmpi(dataset(1),'N')
                trialInclStates(1).trialName = {'Nigel Dissociation'};
            elseif strcmpi(dataset(1),'R')
                trialInclStates(1).trialName = {'Rocky Dissociation'};
            end
            condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
            trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',0}};
        %Nigel BCI
        case {'N20171215','N20180221'}
            trialInclStates(1).trialName = {'Nigel Posture BC Center Out'};
            condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
            trialInclStates(1).inclStates = {{'state','Cursor Freeze','first',25},{'state','Cursor Freeze','first',225}};
        %Rocky BCI
        case {'R20201020','R20201021'}
            trialInclStates(1).trialName = {'Rocky Posture BC Center Out'};
            condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
            trialInclStates(1).inclStates = {{'state','React','first',25},{'state','React','first',225}};     
    end
    trajStruct = getTrajStruct20220419(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'matchConditions',false,'zScoreParams',zScoreParams);
 
    %Replace neurons with PCs
    allTraj = vertcat(trajStruct.avgSmoothFR);
    allTraj = vertcat(allTraj.traj);
    [PCs,score,latent,tsquared,explained,mu] = pca(allTraj); 
    manVarThreshold = 90;
    numPCs = 1;
    while sum(explained(1:numPCs)) < manVarThreshold
    numPCs = numPCs + 1;
    end
    PCs = PCs(:,1:numPCs);
    
    %Add Projections to trajStruct
    for i = 1:size(trajStruct,2)
        %Add average traces to trajStruct
        trajStruct(i).avgSmoothFR.traj = trajStruct(i).avgSmoothFR.traj*PCs;
    end  
    
    
    
    
%% Eliminate outlier trials and conditions, but ensure that things are balanced 
    switch dataset
        case 'E20210706'
            trajStruct([trajStruct.posture]==6) = [];
    end
    
%% Create model trajStruct
    %Get length of shortest mean
    numTimestamps = [];
    for i = 1:size(trajStruct,2)
        numTimestamps(i) = length(trajStruct(i).avgSmoothFR.timestamps);
    end
    figure
    histogram(numTimestamps)
    xlabel('Number of 25ms bins')
    ylabel('Number of trials')
    [minNumTimestamps,i] = min(numTimestamps);

%% Get parameters
    postureList = unique([trajStruct.posture]);
    numPostures = size(postureList,2);
    targetList = unique([trajStruct.target]); 
    numTargets = size(targetList,2);
    numChannels = size(trajStruct(1).avgSmoothFR.traj,2);
    numConditions = size(trajStruct,2);

%% For each condition, for each model, form and store prediction; also save actual trajectory
    %For each posture, form X using other postures.  Use this for CI and
    %target estimates. Then get posture estimate from current posture
    trajStructInd = 1;
    for posture = postureList
        %Form X (shorten to length of shortest traj) Xdims: 1=time, 2=target, 3=posture, 4=channel
        X = NaN(minNumTimestamps,numTargets,numPostures-1,numChannels);
        postureInd = 1;
        for tempPosture = setdiff(postureList,posture)
            targetInd = 1;
            for target = targetList
                traj = trajStruct([trajStruct.posture]==tempPosture & [trajStruct.target]==target).avgSmoothFR.traj;
                X(:,targetInd,postureInd,:) = traj(1:minNumTimestamps,:); 
                targetInd = targetInd + 1;
            end
            postureInd = postureInd+1;
        end
        %Get CI marg
        CITraj = mean(X,[2 3],'omitnan');
        %Get target marg
        targetMarg = mean(X,[3],'omitnan') - CITraj;
        
        targetTrajInd = 1;
        for target = targetList 
            %Get posture component by removing target components
            %from all traj from current posture, excluding held out
            %trajectory, then averaging to get CI1avg, then
            %taking C12avg on C1 component from other posture data
            CI1mat = zeros(minNumTimestamps,numTargets-1,numChannels);
            targetInd = 1;
            for tempTarget = setdiff(targetList,target)  
               CI1mat(:,targetInd,:) = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==tempTarget).avgSmoothFR.traj(1:minNumTimestamps,:); 
               CI1mat(:,targetInd,:) = squeeze(CI1mat(:,targetInd,:)) - squeeze(targetMarg(:,find(targetList==tempTarget),1,:));
               targetInd = targetInd + 1;
            end 
            CI1avg = squeeze(mean(CI1mat,[1 2]));
            CI2avg = squeeze(mean(CITraj,1));
            postureIC = (CI1avg-CI2avg)'; 

            %Get target traj for current target
            targetTraj = targetMarg(:,targetTrajInd,1,:);
            %Store predictions of LOCO model
            CITrajStruct(trajStructInd).posture = posture; CITrajStruct(trajStructInd).target = target;
                CITrajStruct(trajStructInd).traj = squeeze(CITraj);
            CIPTrajStruct(trajStructInd).posture = posture; CIPTrajStruct(trajStructInd).target = target;
                CIPTrajStruct(trajStructInd).traj = squeeze(CITraj) + postureIC;   
            fullTrajStruct(trajStructInd).posture = posture; fullTrajStruct(trajStructInd).target = target;
                fullTrajStruct(trajStructInd).traj = squeeze(CITraj) + postureIC + squeeze(targetTraj);
            actualTrajStruct(trajStructInd).posture = posture; actualTrajStruct(trajStructInd).target = target;
                actualTrajStruct(trajStructInd).traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgSmoothFR.traj(1:minNumTimestamps,:);
            trajStructInd = trajStructInd + 1;
            targetTrajInd = targetTrajInd + 1;
        end
    end

%% Fill in results struct
   %Get overal data mean
   grandMean = mean(vertcat(actualTrajStruct.traj),1);
   
   %For each model, get mean neural distance 
   noModelDist = 0; CIDist = 0; CIPDist = 0; fullDist = 0; totalNumPts = 0;
   
   for posture = postureList
       for target = targetList
           actualTraj = actualTrajStruct([actualTrajStruct.posture]==posture & [actualTrajStruct.target]==target).traj;
           CITraj = CITrajStruct([CITrajStruct.posture]==posture & [CITrajStruct.target]==target).traj;
           CIPTraj = CIPTrajStruct([CIPTrajStruct.posture]==posture & [CIPTrajStruct.target]==target).traj;
           fullTraj = fullTrajStruct([fullTrajStruct.posture]==posture & [fullTrajStruct.target]==target).traj;
           for i = 1:minNumTimestamps
               noModelDist = vecnorm(actualTraj(i,:)-grandMean) + noModelDist;
               %LOCO Model
               CIDist = vecnorm(actualTraj(i,:)-CITraj(i,:)) + CIDist;
               CIPDist = vecnorm(actualTraj(i,:)-CIPTraj(i,:)) + CIPDist;
               fullDist = vecnorm(actualTraj(i,:)-fullTraj(i,:)) + fullDist;
               totalNumPts = totalNumPts + 1;
           end
       end
   end
   
   noModelDist = noModelDist./totalNumPts;
   CIDist = CIDist./totalNumPts;
   CIPDist = CIPDist./totalNumPts; 
   fullDist = fullDist./totalNumPts; 
   
   %For each model, get R2
   SStot = 0; SSres = 0; CISSres = 0; CIPSSres = 0;
   for posture = postureList
       for target = targetList
           actualTraj = actualTrajStruct([actualTrajStruct.posture]==posture & [actualTrajStruct.target]==target).traj;
           CITraj = CITrajStruct([CITrajStruct.posture]==posture & [CITrajStruct.target]==target).traj;
           CIPTraj = CIPTrajStruct([CIPTrajStruct.posture]==posture & [CIPTrajStruct.target]==target).traj;
           fullTraj = fullTrajStruct([fullTrajStruct.posture]==posture & [fullTrajStruct.target]==target).traj;
           for i = 1:minNumTimestamps
               SStot = (actualTraj(i,:)-grandMean)*(actualTraj(i,:)-grandMean)' + SStot;
               SSres = (actualTraj(i,:)-fullTraj(i,:))*(actualTraj(i,:)-fullTraj(i,:))' + SSres;
               CISSres = (actualTraj(i,:)-CITraj(i,:))*(actualTraj(i,:)-CITraj(i,:))' + CISSres;
               CIPSSres = (actualTraj(i,:)-CIPTraj(i,:))*(actualTraj(i,:)-CIPTraj(i,:))' + CIPSSres;
           end
       end
   end
   R2 = 1-(SSres/SStot); 
   CIR2 = 1-(CISSres/SStot); 
   CIPR2 = 1-(CIPSSres/SStot); 
   

%% Store results
   resultStruct(1).dataset = dataset; 
   resultStruct(1).task = task;

   %Mean neural distance 
   resultStruct(1).noModelDist = noModelDist;
   resultStruct(1).CIDist = CIDist;
   resultStruct(1).CIPDist = CIPDist;
   resultStruct(1).fullDist = fullDist;
   
   %R2
   resultStruct(1).R2 = R2;
   resultStruct(1).CIR2 = CIR2;
   resultStruct(1).CIPR2 = CIPR2;
   
   %Actual and predicted trajectories
   resultStruct(1).actualTrajStruct = actualTrajStruct;
   resultStruct(1).predTrajStruct = fullTrajStruct;
   

end