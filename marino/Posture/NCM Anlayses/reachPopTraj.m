clear; clc; clf; close all

%% Setup saveFig   
    saveFig = true;
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Conferences\NCM 2022\Figures\ReachPopTraj';
    set(0, 'DefaultFigureRenderer', 'painters');
    
%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\orli.mat')
    pcmap = orli;

%% Run loop for each dataset    
    holdEpoch = 'center';
    for datasetList = {'E20210706','N20190226'} 
        dataset = datasetList{1,1};
        
        switch dataset
            case 'E20210706'
                plotPostureList = [1:4];
            case 'N20190226'
                plotPostureList = [1:2];
        end
        
        [Data,zScoreParams] = loadData(dataset);
        
        centerHoldData = Data;
        targetHoldData = Data;
        
        %Get trajStruct
        binWidth = 25;
        kernelStdDev = 25;
        trajFields = {'zSmoothFR'};
        trialInclStates = struct('trialName','','inclStates',[]);
        switch dataset
            case 'E20210706'
                %All
                trialInclStates(1).trialName = {'GridReaching'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'state','Target Hold','first',500}};
                %Reach
                reachTrialInclStates(1).trialName = {'GridReaching'};
                reachCondFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                reachTrialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',-50}};
                %Center Hold
                centerHoldTrialInclStates(1).trialName = {'GridReaching'};
                centerHoldCondFields = {{'posture','conditionData','postureID'}};
                centerHoldTrialInclStates(1).inclStates = {{'state','Center Hold','first',0},{'state','Delay','first',50}};
                %Target Hold
                targetHoldTrialInclStates(1).trialName = {'GridReaching'};
                targetHoldCondFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                targetHoldTrialInclStates(1).inclStates = {{'state','Target Hold','first',0},{'state','Success with Reward','first',0}};
            case 'N20190226'
                %All
                trialInclStates(1).trialName = {'Nigel Dissociation'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'state','Target Hold','first',500}};
                %Reach
                reachTrialInclStates(1).trialName = {'Nigel Dissociation'};
                reachCondFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                reachTrialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',-50}};
                %Center Hold
                centerHoldTrialInclStates(1).trialName = {'Nigel Dissociation'};
                centerHoldCondFields = {{'posture','conditionData','postureID'}};
                centerHoldTrialInclStates(1).inclStates = {{'state','Center Hold','first',0},{'state','Reach','first',50}};
                %Target Hold
                targetHoldTrialInclStates(1).trialName = {'Nigel Dissociation'};
                targetHoldCondFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                targetHoldTrialInclStates(1).inclStates = {{'state','Target Hold','first',0},{'state','Success with Reward','first',0}};

        end
        trajStruct = getTrajStruct20220419(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'matchConditions',false,'zScoreParams',zScoreParams);
        reachTrajStruct = getTrajStruct20220419(Data,reachCondFields,trajFields,reachTrialInclStates,binWidth,kernelStdDev,'matchConditions',false,'zScoreParams',zScoreParams);
        centerHoldTrajStruct = getTrajStruct20220419(centerHoldData,centerHoldCondFields,trajFields,centerHoldTrialInclStates,binWidth,kernelStdDev,'matchConditions',false,'zScoreParams',zScoreParams);
        targetHoldTrajStruct = getTrajStruct20220419(targetHoldData,targetHoldCondFields,trajFields,targetHoldTrialInclStates,binWidth,kernelStdDev,'matchConditions',false,'zScoreParams',zScoreParams);
        
        
           %Clear unneccessary variables
        clearvars Data centerHoldData targetHoldData
        % Get minimum number of timestamps in condition averages
        %All 
            numTimestamps = [];
            for i = 1:size(trajStruct,2)
                numTimestamps(i) = length(trajStruct(i).avgSmoothFR.timestamps);
            end
            figure
            histogram(numTimestamps)
            xlabel('Number of 25ms bins')
            ylabel('Number of conditions')
            [minNumTimestamps,i] = min(numTimestamps);
        %Reach
            reachNumTimestamps = [];
            for i = 1:size(reachTrajStruct,2)
                reachNumTimestamps(i) = length(reachTrajStruct(i).avgSmoothFR.timestamps);
            end
            figure
            histogram(reachNumTimestamps)
            xlabel('Number of 25ms bins')
            ylabel('Number of conditions')
            [reachMinNumTimestamps,i] = min(reachNumTimestamps);
        %Center Hold 
            centerHoldNumTimestamps = [];
            for i = 1:size(centerHoldTrajStruct,2)
                centerHoldNumTimestamps(i) = length(centerHoldTrajStruct(i).avgSmoothFR.timestamps);
            end
            figure
            histogram(centerHoldNumTimestamps)
            xlabel('Number of 25ms bins')
            ylabel('Number of conditions')
            [centerHoldMinNumTimestamps,i] = min(centerHoldNumTimestamps);
        %Target Hold
            targetHoldNumTimestamps = [];
            for i = 1:size(targetHoldTrajStruct,2)
                targetHoldNumTimestamps(i) = length(targetHoldTrajStruct(i).avgSmoothFR.timestamps);
            end
            figure
            histogram(targetHoldNumTimestamps)
            xlabel('Number of 25ms bins')
            ylabel('Number of conditions')
            [targetHoldMinNumTimestamps,i] = min(targetHoldNumTimestamps);
            
         %Get posture and target lists 
           %All 
            postureList = unique([trajStruct.posture]);
            numPostures = size(postureList,2);
            targetList = unique([trajStruct.target]); 
            numTargets = size(targetList,2);
            %Reach
            reachPostureList = unique([reachTrajStruct.posture]);
            reachNumPostures = size(reachPostureList,2);
            reachTargetList = unique([reachTrajStruct.target]); 
            reachNumTargets = size(reachTargetList,2);
            %Center Hold
            centerHoldPostureList = unique([centerHoldTrajStruct.posture]);
            centerHoldNumPostures = size(centerHoldPostureList,2); 
            %Target Hold
            targetHoldPostureList = unique([targetHoldTrajStruct.posture]);
            targetHoldNumPostures = size(targetHoldPostureList,2);
            targetHoldTargetList = unique([targetHoldTrajStruct.target]); 
            targetHoldNumTargets = size(targetHoldTargetList,2);
            
        %num Channels and numConditions    
            numChannels = size(reachTrajStruct(1).avgSmoothFR.traj,2);
            numConditions = size(reachTrajStruct,2);
        
        %Remove postures for which not all targets were attempted
            rmPostures = [];
            for posture = reachPostureList
                tempTrajStruct = reachTrajStruct([reachTrajStruct.posture]==posture);
                postureTargetList = [tempTrajStruct.target];
                if ~isequal(postureTargetList,reachTargetList)
                    rmPostures = [rmPostures,posture];
                end
            end
            if ~isempty(rmPostures)
                %All
                rmInd = [trajStruct.posture]==rmPostures;
                trajStruct(rmInd) = [];
                postureList(postureList==rmPostures) = [];
                numPostures = length(postureList);
                %Reach
                rmInd = [reachTrajStruct.posture]==rmPostures;
                reachTrajStruct(rmInd) = [];
                reachPostureList(reachPostureList==rmPostures) = [];
                reachNumPostures = length(reachPostureList);
                %Target Hold
                rmInd = [targetHoldTrajStruct.posture]==rmPostures;
                targetHoldTrajStruct(rmInd) = [];
                targetHoldPostureList(targetHoldPostureList==rmPostures) = [];
                targetHoldNumPostures = length(targetHoldPostureList);              
            end
          
        
        %X
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
        
        %reachX
        reachX = NaN(reachMinNumTimestamps,reachNumTargets,reachNumPostures,numChannels);
        postureInd = 1;
        for posture = reachPostureList
            targetInd = 1;
            for target = reachTargetList
                traj = reachTrajStruct([reachTrajStruct.posture]==posture & [reachTrajStruct.target]==target).avgSmoothFR.traj;
                reachX(:,targetInd,postureInd,:) = traj(1:reachMinNumTimestamps,:); 
                targetInd = targetInd + 1;
            end
            postureInd = postureInd+1;
        end
        
        %centerHoldX
        centerHoldX = NaN(centerHoldMinNumTimestamps,centerHoldNumPostures,numChannels);
        postureInd = 1;
        for i=1:size(centerHoldTrajStruct,2)
            traj = centerHoldTrajStruct(i).avgSmoothFR.traj;
            centerHoldX(:,postureInd,:) = traj(1:centerHoldMinNumTimestamps,:); 
            postureInd = postureInd+1;
        end
        
        %targetHoldX
        targetHoldX = NaN(targetHoldMinNumTimestamps,targetHoldNumTargets,targetHoldNumPostures,numChannels);
        postureInd = 1;
        for posture = targetHoldPostureList
            targetInd = 1;
            for target = targetHoldTargetList
                traj = targetHoldTrajStruct([targetHoldTrajStruct.posture]==posture & [targetHoldTrajStruct.target]==target).avgSmoothFR.traj;
                targetHoldX(:,targetInd,postureInd,:) = traj(end-targetHoldMinNumTimestamps+1:end,:); 
                targetInd = targetInd + 1;
            end
            postureInd = postureInd+1;
        end
             
        % Do marginalizations of reachX
        %reachXdims: 1=time, 2=target, 3=posture, 4=channel
        %Condition-invariant
        reachCIMargOffset = mean(reachX,[1 2 3],'omitnan');
        reachCIMargTraj = mean(reachX,[2 3],'omitnan') - reachCIMargOffset;
        %Posture and Target Offsets
        reachTargetMargOffset = mean(reachX,[1 3],'omitnan') - reachCIMargOffset;
        reachPostureMargOffset = mean(reachX,[1 2],'omitnan') - reachCIMargOffset;
        %Posture and Target Traj
        reachTargetMargTraj = mean(reachX,[3],'omitnan') - reachCIMargTraj  - reachCIMargOffset;
        reachPostureMargTraj = mean(reachX,[2],'omitnan') - reachCIMargTraj - reachCIMargOffset;
        reachTargetTrajNoCP = reachX - reachCIMargTraj - reachCIMargOffset - reachPostureMargTraj - reachPostureMargOffset;
        reachNoCI = reachX - reachCIMargTraj - reachCIMargOffset;
        
        % Do marginalizations of targetHoldX
        %targetHoldXdims: 1=time, 2=target, 3=posture, 4=channel
        %Condition-invariant
        targetHoldCIMargOffset = mean(targetHoldX,[1 2 3],'omitnan');
        targetHoldCIMargTraj = mean(targetHoldX,[2 3],'omitnan') - targetHoldCIMargOffset;
        %%Posture and Target Offsets
        %targetHoldTargetMargOffset = mean(targetHoldX,[1 3],'omitnan') - targetHoldCIMargOffset;
        %targetHoldPostureMargOffset = mean(targetHoldX,[1 2],'omitnan') - targetHoldCIMargOffset;
        %%Posture and Target Traj
        %targetHoldTargetMargTraj = mean(targetHoldX,[3],'omitnan') - targetHoldCIMargTraj  - targetHoldCIMargOffset;
        %targetHoldPostureMargTraj = mean(targetHoldX,[2],'omitnan') - targetHoldCIMargTraj - targetHoldCIMargOffset;
        %targetHoldTargetTrajNoCP = targetHoldX - targetHoldCIMargTraj - targetHoldCIMargOffset - targetHoldPostureMargTraj - targetHoldPostureMargOffset;
        targetHoldNoCI = targetHoldX - targetHoldCIMargTraj - targetHoldCIMargOffset;
        
        % Do marginalizations of centerHoldX
        %centerHoldXdims: 1=time, 2=posture, 3=channel
        %Condition-invariant
        centerHoldCIMargOffset = mean(centerHoldX,[1 2],'omitnan');
        centerHoldCIMargTraj = mean(centerHoldX,[2],'omitnan') - centerHoldCIMargOffset;
        %Posture Offset
        centerHoldPostureMargOffset = mean(centerHoldX,[1 2],'omitnan') - centerHoldCIMargOffset;
        %Posture Traj
        centerHoldNoCI = centerHoldX - centerHoldCIMargTraj - centerHoldCIMargOffset; 
        
       
        
        %Do PCA on all condition averaged trajectories
        allReachTraj = reshape(reachX,[reachNumTargets*reachNumPostures*reachMinNumTimestamps,numChannels]);
        allTargetHoldTraj = reshape(targetHoldX,[targetHoldNumTargets*targetHoldNumPostures*targetHoldMinNumTimestamps,numChannels]);
        allCenterHoldTraj = reshape(centerHoldX,[centerHoldNumPostures*centerHoldMinNumTimestamps,numChannels]);
        
        if strcmp(holdEpoch,'center')
            allTraj = vertcat(allReachTraj,allCenterHoldTraj);
        elseif strcmp(holdEpoch,'target')
            allTraj = vertcat(allReachTraj,allTargetHoldTraj);
        end
        
        [allPCA,score,latent,tsquared,explained,mu] = pca(allTraj); 
        [reachPCA,score,latent,tsquared,explained,mu] = pca(allReachTraj); 
        [centerHoldPCA,score,latent,tsquared,explained,mu] = pca(allCenterHoldTraj); 
        [targetHoldPCA,score,latent,tsquared,explained,mu] = pca(allTargetHoldTraj); 
        Ca = cov(allTraj);
        [Ua,Sa] = eig(Ca);
 

        %Define CIMargTraj, postureMargTraj, targetMargTraj, targetTrajNoCP
        CIMargTraj = reachCIMargTraj;
        if strcmp(holdEpoch,'center')
            postureMargTraj = centerHoldNoCI;
        elseif strcmp(holdEpoch,'target')
            postureMargTraj = targetHoldNoCI;
        end
        targetMargTraj = reachTargetMargTraj;
        targetTrajNoCP = reachTargetTrajNoCP;
        
        CIMargTraj = squeeze(CIMargTraj);
        [CIPCA,score,latent,tsquared,explained,mu] = pca(CIMargTraj); 

        targetMargTraj = squeeze(targetMargTraj);
        targetMargTraj = reshape(targetMargTraj,[reachNumTargets*reachMinNumTimestamps,numChannels]);
        [Dt,score,latent,tsquared,explained,mu] = pca(targetMargTraj); 

        targetTrajNoCP = squeeze(targetTrajNoCP);
        targetTrajNoCP = reshape(targetTrajNoCP,[reachNumTargets*reachNumPostures*reachMinNumTimestamps,numChannels]);
        [targetNoCPPCA,score,latent,tsquared,explained,mu] = pca(targetTrajNoCP); 

        postureMargTraj = squeeze(postureMargTraj);
        if strcmp(holdEpoch,'center')
            postureMargTraj = reshape(postureMargTraj,[centerHoldNumPostures*centerHoldMinNumTimestamps,numChannels]);
        elseif strcmp(holdEpoch,'target')
            postureMargTraj = reshape(postureMargTraj,[targetHoldNumTargets*targetHoldNumPostures*targetHoldMinNumTimestamps,numChannels]);
        end
        [Dp,score,latent,tsquared,explained,mu] = pca(postureMargTraj); 

        
        
        %Peform LDA by posture on all data 
        obsStruct = struct('label',[],'numObs',[],'allObs',[]);
        structInd = 1;   
        for posture = plotPostureList
           obsStruct(structInd).label = posture;
           allObs = [];
           tempTrajStruct = reachTrajStruct([reachTrajStruct.posture]==posture);
           for i = 1:size(tempTrajStruct,2)
               for j = 1:size(tempTrajStruct(i).allSmoothFR,2)
                   timestamps = tempTrajStruct(i).allSmoothFR(j).timestamps;
                   traj = tempTrajStruct(i).allSmoothFR(j).traj;
                   if size(traj,1) > minNumTimestamps
                        traj = traj(1:minNumTimestamps,:);
                   end
                   allObs = vertcat(allObs,traj);
               end
           end
           obsStruct(structInd).allObs = allObs;
           obsStruct(structInd).numObs = size(obsStruct(structInd).allObs,1);
           structInd = structInd + 1;
        end
        [postureLDA] = doLDA(obsStruct);
        
        
         %Orthonormalize combinations of axes
        %[CPTOrth,~] = qr([postureLDA(:,1),targetNoCPPCA(:,1),CIPCA(:,1)]); CPTOrth = CPTOrth(:,1:3);
        %[CPOrth,~] = qr([postureLDA(:,1:2),CIPCA(:,1)]); CPOrth = CPOrth(:,1:3);
        
        targetPCDims = [1,3];
        switch dataset
            case 'E20210706'
                targetPCDims = [1,2];
            case 'N20190226'
                targetPCDims = [1,3];
        end
        [PTOrth,~] = qr([postureLDA(:,1),targetNoCPPCA(:,targetPCDims)]); PTOrth = PTOrth(:,1:3);

        
    
        %Add Projections to trajStruct
        totalVar = trace(cov(allTraj));
        for i = 1:size(trajStruct,2)
            %Add average traces to trajStruct
            trajStruct(i).PTOrth.traj = trajStruct(i).avgSmoothFR.traj*PTOrth;
            %Get VAF
            trajStruct(i).PTOrth.VAF =  100.*(diag(cov(allTraj*PTOrth))')./totalVar;
        end
        
        minNumTimestamps = 10;
        fs = 14;
        
        
        %Plot - Orthographic
        figure
        hold on
        timePts = 1:minNumTimestamps;
        xDim = 1; yDim = 2; zDim = 3;
        for posture = plotPostureList
            for target = targetList
                if any([trajStruct.posture]==posture & [trajStruct.target]==target)
                    traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).PTOrth.traj(timePts,:); 
                    plot3(traj(:,xDim),traj(:,yDim),traj(:,zDim),'Color',pcmap(posture,:),'LineWidth',2);
                    plot3(traj(1,xDim),traj(1,yDim),traj(1,zDim),'.','MarkerSize',20,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
                    plot3(traj(end,xDim),traj(end,yDim),traj(end,zDim),'<','MarkerSize',5,'Color',pcmap(posture,:),'MarkerFaceColor',pcmap(posture,:));
                end
            end
        end
        grid on
        xticklabels({}); yticklabels({}); zticklabels({}); 
        VAF = round(trajStruct(1).PTOrth.VAF);
        xlabel(['Posture LDA 1 (',num2str(VAF(xDim)),'%)'])
        ylabel(['Target PC 1 (',num2str(VAF(yDim)),'%)'])
        zlabel(['Target PC 2 (',num2str(VAF(zDim)),'%)']) 
        view([-70 15])
        set(gca,'fontname','arial')
        set(gca,'fontsize',fs)
        if saveFig
            saveas(gcf,fullfile(saveDir,[dataset,'_reachPopTraj.svg']));
            saveas(gcf,fullfile(saveDir,[dataset,'_reachPopTraj.fig']));
        end
    end
    
    
%% Local function for performing LDA
     function [LDAproj] = doLDA(obsStruct)
        %Preallocate
        minNumObs = min([obsStruct.numObs]);
        numClasses = size(obsStruct,2);
        numDims = size(obsStruct(1).allObs,2);
        obs = NaN(minNumObs*numClasses,numDims); labels =  NaN(minNumObs*numClasses,1); 

        %Fill
        k = 1;
        for i = 1:size(obsStruct,2)
            totalClassObs = size(obsStruct(i).allObs,1);
            obs(k:k+minNumObs-1,:) = obsStruct(i).allObs(randsample(totalClassObs,minNumObs),:);
            labels(k:k+minNumObs-1,1) = ones(minNumObs,1).*obsStruct(i).label;
            k = k+minNumObs;
        end

        LDAproj = fisherLDA(obs, labels);
        [LDAproj,~] = qr(LDAproj);
        LDAproj = LDAproj(:,1:numClasses-1);
     end  