function [resultStruct] = leaveOneConditionOutCompAnalysisReachVsBCI(dataset,task,epoch)

%% Set up results struct
    resultStruct = struct('dataset','','task','','epoch','','noModelDist',[],'CIDist',[],...
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
            trialInclStates(1).inclStates = {{'state','Step 1','first',0},{'state','Step 1','first',300}};
        %Earl Reach
        case {'E20210706','E20210707','E20210708','E20210709','E20210710'}
            trialInclStates(1).trialName = {'GridReaching'};
            condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
            trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',100}};
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
                trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'state','Hold','first',0}};
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
            trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',100}};

        %Nigel BCI
        case {'N20171215','N20180221'}
            trialInclStates(1).trialName = {'Nigel Posture BC Center Out'};
            condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
            trialInclStates(1).inclStates = {{'state','Cursor Release','first',0},{'state','Cursor Release','first',300}};
            
    end
    trajStruct = getTrajStruct20220419(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'matchConditions',false,'zScoreParams',zScoreParams);
 
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
    %Form XFull using all data
    XFull = NaN(minNumTimestamps,numTargets,numPostures,numChannels);
    postureInd = 1;
    for posture = postureList
        targetInd = 1;
        for target = targetList
            traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgSmoothFR.traj;
            XFull(:,targetInd,postureInd,:) = traj(1:minNumTimestamps,:); 
            targetInd = targetInd + 1;
        end
        postureInd = postureInd+1;
    end
    %Get CI marg
    CITrajFull = mean(XFull,[2 3],'omitnan');
    %Get target marg
    targetMargFull = mean(XFull,[3],'omitnan') - CITrajFull;

    %For each posture, form X using other postures.  Use this for CI and
    %target estimates. Then get posture estimate from current posture
    trajStructInd = 1;
    for posture = postureList
        
        %Get postureICFull (for model using all data)
        ICs = zeros(numTargets,numChannels);
        for target = targetList
           ICs(targetInd,:) = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgSmoothFR.traj(1,:); 
           targetInd = targetInd + 1;
        end
        postureICFull = mean(ICs,1)-squeeze(CITrajFull(1,:,:,:))';
        
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
            %Get posture IC excluding current target
            ICs = zeros(numTargets-1,numChannels);
            targetInd = 1;
            for tempTarget = setdiff(targetList,target)  
               ICs(targetInd,:) = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==tempTarget).avgSmoothFR.traj(1,:); 
               targetInd = targetInd + 1;
            end
            postureIC = mean(ICs,1)-squeeze(CITraj(1,:,:,:))';
            %Get target traj for current target
            targetTraj = targetMarg(:,targetTrajInd,1,:);
            targetTrajFull = targetMargFull(:,targetTrajInd,1,:);
            %Store predictions of LOCO model
            CITrajStruct(trajStructInd).posture = posture; CITrajStruct(trajStructInd).target = target;
                CITrajStruct(trajStructInd).traj = squeeze(CITraj);
            CIPTrajStruct(trajStructInd).posture = posture; CIPTrajStruct(trajStructInd).target = target;
                CIPTrajStruct(trajStructInd).traj = squeeze(CITraj) + postureIC;   
            fullTrajStruct(trajStructInd).posture = posture; fullTrajStruct(trajStructInd).target = target;
                fullTrajStruct(trajStructInd).traj = squeeze(CITraj) + postureIC + squeeze(targetTraj);
                
            %Store predictions of Full Model
            CITrajFullStruct(trajStructInd).posture = posture; CITrajFullStruct(trajStructInd).target = target;
                CITrajFullStruct(trajStructInd).traj = squeeze(CITrajFull);
            CIPTrajFullStruct(trajStructInd).posture = posture; CIPTrajFullStruct(trajStructInd).target = target;
                CIPTrajFullStruct(trajStructInd).traj = squeeze(CITrajFull) + postureICFull;   
            fullFullTrajStruct(trajStructInd).posture = posture; fullFullTrajStruct(trajStructInd).target = target;
                fullFullTrajStruct(trajStructInd).traj = squeeze(CITrajFull) + postureICFull + squeeze(targetTrajFull);
                
                
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
   noModelDist = 0;
   CIDist = 0; CIFullDist = 0;
   CIPDist = 0; CIPFullDist = 0;
   fullDist = 0; fullFullDist = 0;
   totalNumPts = 0;
   
   for posture = postureList
       for target = targetList
           actualTraj = actualTrajStruct([actualTrajStruct.posture]==posture & [actualTrajStruct.target]==target).traj;
           CITraj = CITrajStruct([CITrajStruct.posture]==posture & [CITrajStruct.target]==target).traj;
           CIPTraj = CIPTrajStruct([CIPTrajStruct.posture]==posture & [CIPTrajStruct.target]==target).traj;
           fullTraj = fullTrajStruct([fullTrajStruct.posture]==posture & [fullTrajStruct.target]==target).traj;
           CITrajFull = CITrajFullStruct([CITrajFullStruct.posture]==posture & [CITrajFullStruct.target]==target).traj;
           CIPTrajFull = CIPTrajFullStruct([CIPTrajFullStruct.posture]==posture & [CIPTrajFullStruct.target]==target).traj;
           fullFullTraj = fullFullTrajStruct([fullFullTrajStruct.posture]==posture & [fullFullTrajStruct.target]==target).traj;
           for i = 1:minNumTimestamps
               noModelDist = vecnorm(actualTraj(i,:)-grandMean) + noModelDist;
               %LOCO Model
               CIDist = vecnorm(actualTraj(i,:)-CITraj(i,:)) + CIDist;
               CIPDist = vecnorm(actualTraj(i,:)-CIPTraj(i,:)) + CIPDist;
               fullDist = vecnorm(actualTraj(i,:)-fullTraj(i,:)) + fullDist;
               %Full Model
               CIFullDist = vecnorm(actualTraj(i,:)-CITrajFull(i,:)) + CIFullDist;
               CIPFullDist = vecnorm(actualTraj(i,:)-CIPTrajFull(i,:)) + CIPFullDist;
               fullFullDist = vecnorm(actualTraj(i,:)-fullFullTraj(i,:)) + fullFullDist;
               
               totalNumPts = totalNumPts + 1;
           end
       end
   end
   
   noModelDist = noModelDist./totalNumPts;
   CIDist = CIDist./totalNumPts; CIFullDist = CIFullDist./totalNumPts;
   CIPDist = CIPDist./totalNumPts; CIPFullDist = CIPFullDist./totalNumPts;
   fullDist = fullDist./totalNumPts; fullFullDist = fullFullDist./totalNumPts;
   
   %For each model, get R2
   SStot = 0; SSres = 0; fullSSres = 0;
   CISSres = 0; CIFullSSres = 0;
   CIPSSres = 0; CIPFullSSres = 0;
   for posture = postureList
       for target = targetList
           actualTraj = actualTrajStruct([actualTrajStruct.posture]==posture & [actualTrajStruct.target]==target).traj;
           CITraj = CITrajStruct([CITrajStruct.posture]==posture & [CITrajStruct.target]==target).traj;
           CIPTraj = CIPTrajStruct([CIPTrajStruct.posture]==posture & [CIPTrajStruct.target]==target).traj;
           fullTraj = fullTrajStruct([fullTrajStruct.posture]==posture & [fullTrajStruct.target]==target).traj;
           CITrajFull = CITrajFullStruct([CITrajFullStruct.posture]==posture & [CITrajFullStruct.target]==target).traj;
           CIPTrajFull = CIPTrajFullStruct([CIPTrajFullStruct.posture]==posture & [CIPTrajFullStruct.target]==target).traj;
           fullFullTraj = fullFullTrajStruct([fullFullTrajStruct.posture]==posture & [fullFullTrajStruct.target]==target).traj;
           
           for i = 1:minNumTimestamps
               SStot = (actualTraj(i,:)-grandMean)*(actualTraj(i,:)-grandMean)' + SStot;
               SSres = (actualTraj(i,:)-fullTraj(i,:))*(actualTraj(i,:)-fullTraj(i,:))' + SSres;
               CISSres = (actualTraj(i,:)-CITraj(i,:))*(actualTraj(i,:)-CITraj(i,:))' + CISSres;
               CIPSSres = (actualTraj(i,:)-CIPTraj(i,:))*(actualTraj(i,:)-CIPTraj(i,:))' + CIPSSres;
               
               fullSSres = (actualTraj(i,:)-fullFullTraj(i,:))*(actualTraj(i,:)-fullFullTraj(i,:))' + fullSSres;
               CIFullSSres = (actualTraj(i,:)-CITrajFull(i,:))*(actualTraj(i,:)-CITrajFull(i,:))' + CIFullSSres;
               CIPFullSSres = (actualTraj(i,:)-CIPTrajFull(i,:))*(actualTraj(i,:)-CIPTrajFull(i,:))' + CIPFullSSres;
           end
       end
   end
   R2 = 1-(SSres/SStot); R2full = 1-(fullSSres/SStot); 
   CIR2 = 1-(CISSres/SStot); CIFullR2 = 1-(CIFullSSres/SStot);
   CIPR2 = 1-(CIPSSres/SStot); CIPFullR2 = 1-(CIPFullSSres/SStot);
   
%% Get Willett result
    %Form X using all data (shorten to length of shortest traj) Xdims: 1=time, 2=target, 3=posture, 4=channel
    %Xdims: 1=time, 2=target, 3=posture, 4=channel
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
    %Condition-invariant
    CIMargOffset = mean(X,[1 2 3],'omitnan');
    CIMargTraj = mean(X,[2 3],'omitnan') - CIMargOffset;
    %Posture and Target Offsets
    targetMargOffset = mean(X,[1 3],'omitnan') - CIMargOffset;
    postureMargOffset = mean(X,[1 2],'omitnan') - CIMargOffset;
    %Posture and Target Traj
    targetMargTraj = mean(X,[3],'omitnan') - CIMargTraj - targetMargOffset - CIMargOffset;
    postureMargTraj = mean(X,[2],'omitnan') - CIMargTraj - postureMargOffset - CIMargOffset;
    %Interaction
    intMargOffset = mean(X,[1],'omitnan') - targetMargOffset - postureMargOffset - CIMargOffset ;
    intMargTraj = X-intMargOffset-targetMargTraj-postureMargTraj-targetMargOffset-postureMargOffset-CIMargTraj-CIMargOffset;
    %Posture and Target Traj
    %targetMargTraj = mean(X,[3],'omitnan') - CIMargTraj  - CIMargOffset;
    %postureMargTraj = mean(X,[2],'omitnan') - CIMargTraj - CIMargOffset;
    %targetTrajNoCP = X - CIMargTraj - CIMargOffset - postureMargTraj - postureMargOffset;

    %Compute Variance accounted for by each
     %Get total var
    totalVar = 0;
    numObs = 0;
    postureInd = 1;
    for posture = postureList
        targetInd = 1;
        for target = targetList
            for i = 1:minNumTimestamps
                totalVar = totalVar + (X(i,targetInd,postureInd,:)-CIMargOffset).^2;
                numObs = numObs + 1;
            end
            targetInd = targetInd + 1;
        end
        postureInd = postureInd+1;
    end
    totalVar = sum(totalVar)./(numObs-1);
    
    %Get CI Var 
    CIVar = 0;
    numObs = 0;
    tempMean = mean(CIMargTraj,[1 2 3],'omitnan');
    postureInd = 1;
    for posture = postureList
        targetInd = 1;
        for target = targetList
            for i = 1:minNumTimestamps
%                 CIVar = CIVar + (CIMargTraj(i,1,1,:)+CIMargOffset-CIMargOffset).^2;
                CIVar = CIVar + (CIMargTraj(i,1,1,:)-tempMean).^2;
                numObs = numObs + 1;
            end
            targetInd = targetInd + 1;
        end
        postureInd = postureInd+1;
    end
    CIVar = sum(CIVar)./(numObs-1);
    CIVarPct = (CIVar/totalVar)*100;
    
    %Get Target Var
    targetVar = 0;
    numObs = 0;
    tempMean = mean(targetMargTraj+targetMargOffset,[1 2 3],'omitnan');
    postureInd = 1;
    for posture = postureList
        targetInd = 1;
        for target = targetList
            for i = 1:minNumTimestamps
%                 targetVar = targetVar + (targetMargTraj(i,target,1,:)+targetMargOffset(1,target,1,:)+CIMargTraj(i,1,1,:)+CIMargOffset-CIMargOffset).^2;
                targetVar = targetVar + (targetMargTraj(i,targetInd,1,:)+targetMargOffset(1,targetInd,1,:)-tempMean).^2;
                numObs = numObs + 1;
            end
            targetInd = targetInd + 1;
        end
        postureInd = postureInd+1;
    end
    targetVar = sum(targetVar)./(numObs-1);
    targetVarPct = (targetVar/totalVar)*100;
    
    
    %Get Target Offset Var
    targetOffsetVar = 0;
    numObs = 0;
    tempMean = mean(targetMargTraj+targetMargOffset,[1 2 3],'omitnan');
    postureInd = 1;
    for posture = postureList
        targetInd = 1;
        for target = targetList
            for i = 1:minNumTimestamps
                targetOffsetVar = targetOffsetVar + (targetMargOffset(1,targetInd,1,:)-tempMean).^2;
                numObs = numObs + 1;
            end
            targetInd = targetInd + 1;
        end
        postureInd = postureInd+1;
    end
    targetOffsetVar = sum(targetOffsetVar)./(numObs-1);
    targetOffsetVarPct = (targetOffsetVar/targetVar)*100;
    
    
    %Get Posture Var
    postureVar = 0;
    numObs = 0;
    tempMean = mean(postureMargTraj+postureMargOffset,[1 2 3],'omitnan');
    postureInd = 1;
    for posture = postureList
        targetInd = 1;
        for target = targetList
            for i = 1:minNumTimestamps
                %postureVar = postureVar + (postureMargTraj(i,1,posture,:)+postureMargOffset(1,1,posture,:)+CIMargTraj(i,1,1,:)+CIMargOffset-CIMargOffset).^2;
                postureVar = postureVar + (postureMargTraj(i,1,postureInd,:)+postureMargOffset(1,1,postureInd,:)-tempMean).^2;
                numObs = numObs + 1;
            end
            targetInd = targetInd + 1;
        end
        postureInd = postureInd+1;
    end
    postureVar = sum(postureVar)./(numObs-1);
    postureVarPct = (postureVar/totalVar)*100;
    
    %Get Posture Offset Var
    postureOffsetVar = 0;
    numObs = 0;
    tempMean = mean(postureMargTraj+postureMargOffset,[1 2 3],'omitnan');
    postureInd = 1;
    for posture = postureList
        targetInd = 1;
        for target = targetList
            for i = 1:minNumTimestamps
                postureOffsetVar = postureOffsetVar + (postureMargOffset(1,1,postureInd,:)-tempMean).^2;
                numObs = numObs + 1;
            end
            targetInd = targetInd + 1;
        end
        postureInd = postureInd+1;
    end
    postureOffsetVar = sum(postureOffsetVar)./(numObs-1);
    postureOffsetVarPct = (postureOffsetVar/postureVar)*100;
    
    %Get Interaction Var
    intVar = 0;
    numObs = 0;
    tempMean = mean(intMargTraj+intMargOffset,[1 2 3],'omitnan');
    postureInd = 1;
    for posture = postureList
        targetInd = 1;
        for target = targetList
            for i = 1:minNumTimestamps
%                 intVar = intVar + (intMargTraj(i,target,posture,:)+intMargOffset(1,target,posture,:)+targetMargTraj(i,target,1,:)+targetMargOffset(1,target,1,:)+postureMargTraj(i,1,posture,:)+postureMargOffset(1,1,posture,:)+CIMargTraj(i,1,1,:)+CIMargOffset-CIMargOffset).^2;
                intVar = intVar + (intMargTraj(i,targetInd,postureInd,:)+intMargOffset(1,targetInd,postureInd,:)-tempMean).^2;
                numObs = numObs + 1;
            end
            targetInd = targetInd + 1;
        end
        postureInd = postureInd+1;
    end
    intVar = sum(intVar)./(numObs-1);
    intVarPct = (intVar/totalVar)*100;
    
    %Get Int Offset Var Pct
    intOffsetVar = 0;
    numObs = 0;
    tempMean = mean(intMargTraj+intMargOffset,[1 2 3],'omitnan');
    postureInd = 1;
    for posture = postureList
        targetInd = 1;
        for target = targetList
            for i = 1:minNumTimestamps
%                 intVar = intVar + (intMargTraj(i,target,posture,:)+intMargOffset(1,target,posture,:)+targetMargTraj(i,target,1,:)+targetMargOffset(1,target,1,:)+postureMargTraj(i,1,posture,:)+postureMargOffset(1,1,posture,:)+CIMargTraj(i,1,1,:)+CIMargOffset-CIMargOffset).^2;
%                 intVar = intVar + (intMargTraj(i,targetInd,postureInd,:)+intMargOffset(1,targetInd,postureInd,:)-tempMean).^2;
                intOffsetVar = intOffsetVar + (intMargOffset(1,targetInd,postureInd,:)-tempMean).^2;
                numObs = numObs + 1;
            end
            targetInd = targetInd + 1;
        end
        postureInd = postureInd+1;
    end
    intOffsetVar = sum(intOffsetVar)./(numObs-1);
    intOffsetVarPct = (intOffsetVar/intVar)*100;

%% Store results
   resultStruct(1).dataset = dataset; 
   resultStruct(1).task = task;
   resultStruct(1).epoch = epoch;
   resultStruct(1).noModelDist = noModelDist;
   resultStruct(1).CIDist = CIDist;
   resultStruct(1).CIPDist = CIPDist;
   resultStruct(1).fullDist = fullDist;
   
   resultStruct(1).CIFullDist = CIFullDist;
   resultStruct(1).CIPFullDist = CIPFullDist;
   resultStruct(1).fullFullDist = fullFullDist;
   
   resultStruct(1).actualTrajStruct = actualTrajStruct;
   resultStruct(1).predTrajStruct = fullTrajStruct;
   
   resultStruct(1).R2 = R2;
   resultStruct(1).CIR2 = CIR2;
   resultStruct(1).CIPR2 = CIPR2;
   
   resultStruct(1).fullR2 = R2full;
   resultStruct(1).CIFullR2 = CIFullR2;
   resultStruct(1).CIPFullR2 = CIPFullR2;
   
   resultStruct(1).CIVarPct = CIVarPct;
   resultStruct(1).postureVarPct = postureVarPct;
   resultStruct(1).targetVarPct = targetVarPct;
   resultStruct(1).intVarPct = intVarPct;
   %figure
   %bar([noModelDist,CIDist,CIPDist,fullDist])
   
end