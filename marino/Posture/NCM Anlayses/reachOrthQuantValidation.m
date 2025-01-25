clear; clc; clf; close all

%% Setup saveFig   
    saveFig = true;
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Conferences\NCM 2022\Figures\ReachOrthQuantValidation';
    set(0, 'DefaultFigureRenderer', 'painters');
    
%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat');
    tcmap = customRainbow;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\orli.mat')
    pcmap = orli;
    
%% Run loop for each dataset    
    for datasetList = {'E20210706','N20190226'} 
        %datasetList = {'N20190226'}
        %datasetList = {'E20210706'} 
        %Load data
        dataset = datasetList{1,1};
        [Data,zScoreParams] = loadData(dataset);

        %Look at distribution of center hold times; make holdData trials
        %bigger than 25th pctile
        figure
        kinData = [Data.kinData];
        centerHoldTimes = [kinData.centerHoldTime];
        histogram(centerHoldTimes);
        xlabel('Center Hold Time (ms)')
        ylabel('Count')
        centerHoldTimeThreshold = prctile(centerHoldTimes,10);
        centerHoldData = Data(centerHoldTimes > centerHoldTimeThreshold);
        if saveFig
            saveas(gcf,fullfile(saveDir,[dataset,'_CenterHoldDist.svg']));
        end
        
        %Look at distribution of target hold times; set threshold
        figure
        kinData = [Data.kinData];
        targetHoldTimes = [kinData.targetHoldTime];
        histogram(targetHoldTimes);
        xlabel('Target Hold Time (ms)')
        ylabel('Count')
        targetHoldTimeThreshold = prctile(targetHoldTimes,10);
        %targetHoldTimeThreshold = 800;
        targetHoldData = Data(targetHoldTimes > targetHoldTimeThreshold);
        if saveFig
            saveas(gcf,fullfile(saveDir,[dataset,'_TargetHoldDist.svg']));
        end
        
        %Get trajStruct
        binWidth = 25;
        kernelStdDev = 25;
        trajFields = {'zSmoothFR','marker','markerVel'};
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
                reachTrialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-300},{'kin','moveOnsetTime','first',300}};
                %Center Hold
                centerHoldTrialInclStates(1).trialName = {'GridReaching'};
                centerHoldCondFields = {{'posture','conditionData','postureID'}};
                centerHoldTrialInclStates(1).inclStates = {{'state','Center Hold','first',0},{'state','Delay','first',0}};
                %Target Hold
                targetHoldTrialInclStates(1).trialName = {'GridReaching'};
                targetHoldCondFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                targetHoldTrialInclStates(1).inclStates = {{'state','Target Hold','first',-200},{'state','Success with Reward','first',0}};
            case 'N20190226'
                %All
                trialInclStates(1).trialName = {'Nigel Dissociation'};
                condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'state','Target Hold','first',500}};
                %Reach
                reachTrialInclStates(1).trialName = {'Nigel Dissociation'};
                reachCondFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                reachTrialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-300},{'kin','moveOnsetTime','first',300}};
                %Center Hold
                centerHoldTrialInclStates(1).trialName = {'Nigel Dissociation'};
                centerHoldCondFields = {{'posture','conditionData','postureID'}};
                centerHoldTrialInclStates(1).inclStates = {{'state','Center Hold','first',0},{'state','Reach','first',0}};
                %Target Hold
                targetHoldTrialInclStates(1).trialName = {'Nigel Dissociation'};
                targetHoldCondFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
                targetHoldTrialInclStates(1).inclStates = {{'state','Target Hold','first',-200},{'state','Success with Reward','first',0}};
        end
        trajStruct = getTrajStruct20220419(Data,reachCondFields,trajFields,trialInclStates,binWidth,kernelStdDev,'matchConditions',false,'zScoreParams',zScoreParams);
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
        centerHoldPostureMargTraj = centerHoldX - centerHoldCIMargTraj - centerHoldCIMargOffset;        
        
        %Do PCA on all condition averaged trajectories
        allReachTraj = reshape(reachX,[reachNumTargets*reachNumPostures*reachMinNumTimestamps,numChannels]);
        allTargetHoldTraj = reshape(targetHoldX,[targetHoldNumTargets*targetHoldNumPostures*targetHoldMinNumTimestamps,numChannels]);
        allCenterHoldTraj = reshape(centerHoldX,[centerHoldNumPostures*centerHoldMinNumTimestamps,numChannels]);
        allTraj = vertcat(allReachTraj,allTargetHoldTraj);
        [allPCA,score,latent,tsquared,explained,mu] = pca(allTraj); 
        [reachPCA,score,latent,tsquared,explained,mu] = pca(allReachTraj); 
        [centerHoldPCA,score,latent,tsquared,explained,mu] = pca(allCenterHoldTraj); 
        [targetHoldPCA,score,latent,tsquared,explained,mu] = pca(allTargetHoldTraj); 
        Ca = cov(allTraj);
        [Ua,Sa] = eig(Ca);
        
        
        %Do PCA on Marginalizations
            %Define CIMargTraj, postureMargTraj, targetMargTraj, targetTrajNoCP
            CIMargTraj = reachCIMargTraj;
            postureMargTraj = targetHoldNoCI;
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
            postureMargTraj = reshape(postureMargTraj,[targetHoldNumTargets*targetHoldNumPostures*targetHoldMinNumTimestamps,numChannels]);
            [Dp,score,latent,tsquared,explained,mu] = pca(postureMargTraj); 

        %Add projections and Mean FR to trajStruct
            %All
            for i = 1:size(trajStruct,2)
                trajStruct(i).allPCA.traj = trajStruct(i).avgSmoothFR.traj*allPCA;
                trajStruct(i).reachPCA.traj = trajStruct(i).avgSmoothFR.traj*reachPCA;
                trajStruct(i).centerHoldPCA.traj = trajStruct(i).avgSmoothFR.traj*centerHoldPCA;
                trajStruct(i).targetHoldPCA.traj = trajStruct(i).avgSmoothFR.traj*targetHoldPCA;
                trajStruct(i).targetPCA.traj = trajStruct(i).avgSmoothFR.traj*Dt;
                trajStruct(i).posturePCA.traj = trajStruct(i).avgSmoothFR.traj*Dp;
                trajStruct(i).popMean.traj = mean(trajStruct(i).avgSmoothFR.traj,2);
            end
        
            %Reach
            for i = 1:size(reachTrajStruct,2)
                reachTrajStruct(i).allPCA.traj = reachTrajStruct(i).avgSmoothFR.traj*allPCA;
                reachTrajStruct(i).reachPCA.traj = reachTrajStruct(i).avgSmoothFR.traj*reachPCA;
                reachTrajStruct(i).centerHoldPCA.traj = reachTrajStruct(i).avgSmoothFR.traj*centerHoldPCA;
                reachTrajStruct(i).targetHoldPCA.traj = reachTrajStruct(i).avgSmoothFR.traj*targetHoldPCA;
                reachTrajStruct(i).targetPCA.traj = reachTrajStruct(i).avgSmoothFR.traj*Dt;
                reachTrajStruct(i).posturePCA.traj = reachTrajStruct(i).avgSmoothFR.traj*Dp;
                reachTrajStruct(i).popMean.traj = mean(reachTrajStruct(i).avgSmoothFR.traj,2);
            end
            
            %Center Hold
            for i = 1:size(centerHoldTrajStruct,2)
                centerHoldTrajStruct(i).allPCA.traj = centerHoldTrajStruct(i).avgSmoothFR.traj*allPCA;
                centerHoldTrajStruct(i).reachPCA.traj = centerHoldTrajStruct(i).avgSmoothFR.traj*reachPCA;
                centerHoldTrajStruct(i).centerHoldPCA.traj = centerHoldTrajStruct(i).avgSmoothFR.traj*centerHoldPCA;
                centerHoldTrajStruct(i).targetHoldPCA.traj = centerHoldTrajStruct(i).avgSmoothFR.traj*targetHoldPCA;
                centerHoldTrajStruct(i).targetPCA.traj = centerHoldTrajStruct(i).avgSmoothFR.traj*Dt;
                centerHoldTrajStruct(i).posturePCA.traj = centerHoldTrajStruct(i).avgSmoothFR.traj*Dp;
                centerHoldTrajStruct(i).popMean.traj = mean(centerHoldTrajStruct(i).avgSmoothFR.traj,2);
            end
            
            %Target Hold
            for i = 1:size(targetHoldTrajStruct,2)
                targetHoldTrajStruct(i).allPCA.traj = targetHoldTrajStruct(i).avgSmoothFR.traj*allPCA;
                targetHoldTrajStruct(i).reachPCA.traj = targetHoldTrajStruct(i).avgSmoothFR.traj*reachPCA;
                targetHoldTrajStruct(i).centerHoldPCA.traj = targetHoldTrajStruct(i).avgSmoothFR.traj*centerHoldPCA;
                targetHoldTrajStruct(i).targetHoldPCA.traj = targetHoldTrajStruct(i).avgSmoothFR.traj*targetHoldPCA;
                targetHoldTrajStruct(i).targetPCA.traj = targetHoldTrajStruct(i).avgSmoothFR.traj*Dt;
                targetHoldTrajStruct(i).posturePCA.traj = targetHoldTrajStruct(i).avgSmoothFR.traj*Dp;
                targetHoldTrajStruct(i).popMean.traj = mean(targetHoldTrajStruct(i).avgSmoothFR.traj,2);
            end
        
       
        %Plot Reach Marker Traj
        figure
        for posture = reachPostureList
            for target = reachTargetList
                traj = reachTrajStruct([reachTrajStruct.posture]==posture & [reachTrajStruct.target]==target).avgMarker.traj;
                plot(traj(:,1),traj(:,2),'Color',tcmap(target,:));
                hold on
            end
        end
        xlabel('x (mm)')
        ylabel('y (mm)')
        axis equal
        if saveFig
            saveas(gcf,fullfile(saveDir,[dataset,'_ReachMarkerTraj.svg']));
        end
        
        %Plot Reach Marker Time Course
        for epoch = {'reach','target hold'}
            if strcmp(epoch,'reach')
                tempTrajStruct = reachTrajStruct;
            elseif strcmp(epoch,'target hold')
                tempTrajStruct = targetHoldTrajStruct;
            end
            figure
            for posture = 2
                for target = reachTargetList
                    time = tempTrajStruct([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target).avgMarker.timestamps;
                    posTraj = tempTrajStruct([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target).avgMarker.traj;
                    velTraj = tempTrajStruct([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target).avgMarkerVel.traj;
                    subplot(4,1,1)
                        plot(time,posTraj(:,1),'Color',tcmap(target,:));
                        ylabel('x (mm)')
                        hold on
                    subplot(4,1,2)
                        plot(time,posTraj(:,2),'Color',tcmap(target,:));
                        ylabel('y (mm)')
                        hold on
                    subplot(4,1,3)
                        plot(time,velTraj(:,1),'Color',tcmap(target,:));
                        ylabel('vx (m/s)')
                        hold on
                    subplot(4,1,4)
                        plot(time,velTraj(:,2),'Color',tcmap(target,:));
                        ylabel('vy (m/s)')
                        hold on
                end
            end
            xlabel('t (ms)')
            sgtitle(epoch{1,1})
            if saveFig
                saveas(gcf,fullfile(saveDir,[dataset,'_',epoch{1,1},'_MarkerTimeCourse.svg']));
            end
        end
        
        epoch = {'center hold'};
        tempTrajStruct = centerHoldTrajStruct;
        figure
        for posture = centerHoldPostureList
            time = tempTrajStruct([tempTrajStruct.posture]==posture ).avgMarker.timestamps;
            posTraj = tempTrajStruct([tempTrajStruct.posture]==posture).avgMarker.traj;
            velTraj = tempTrajStruct([tempTrajStruct.posture]==posture).avgMarkerVel.traj;
            subplot(4,1,1)
                plot(time,posTraj(:,1),'Color',tcmap(target,:));
                ylabel('x (mm)')
                hold on
            subplot(4,1,2)
                plot(time,posTraj(:,2),'Color',tcmap(target,:));
                ylabel('y (mm)')
                hold on
            subplot(4,1,3)
                plot(time,velTraj(:,1),'Color',tcmap(target,:));
                ylabel('vx (m/s)')
                hold on
            subplot(4,1,4)
                plot(time,velTraj(:,2),'Color',tcmap(target,:));
                ylabel('vy (m/s)')
                hold on
        end
        xlabel('t (ms)')
        sgtitle(epoch{1,1})
        if saveFig
            saveas(gcf,fullfile(saveDir,[dataset,'_',epoch{1,1},'_MarkerTimeCourse.svg']));
        end
        
        %Plot Neuron Time Course
        for epoch = {'reach','target hold'}
            if strcmp(epoch,'reach')
                tempTrajStruct = reachTrajStruct;
            elseif strcmp(epoch,'target hold')
                tempTrajStruct = targetHoldTrajStruct;
            end
            f = figure;
            f.Position = [20 20 900 500];
            for posture = 2
                for target = reachTargetList
                    time = tempTrajStruct([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target).avgSmoothFR.timestamps;
                    traj = tempTrajStruct([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target).avgSmoothFR.traj;
                    numNeurons = size(traj,2);
                    for i = 1:25
                       subplot(5,5,i)
                           plot(time,traj(:,i),'Color',tcmap(target,:));
                           hold on
                           xlabel('t (ms)')
                           ylabel('FR (std)')
                    end
                end
            end
            sgtitle(epoch{1,1})
            if saveFig
                saveas(gcf,fullfile(saveDir,[dataset,'_',epoch{1,1},'_NeuronTimeCourse.svg']));
            end
        end
       
        %Plot Pop Mean Time Course
        for epoch = {'reach','target hold'}
            if strcmp(epoch,'reach')
                tempTrajStruct = reachTrajStruct;
            elseif strcmp(epoch,'target hold')
                tempTrajStruct = targetHoldTrajStruct;
            end
            figure
            for posture = 2
                for target = reachTargetList
                    time = tempTrajStruct([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target).avgSmoothFR.timestamps;
                    traj = tempTrajStruct([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target).popMean.traj;
                    numNeurons = size(traj,2);
                   plot(time,traj(:,1),'Color',tcmap(target,:));
                   hold on
                   xlabel('t (ms)')
                   ylabel('FR (std)')
                end
            end
            sgtitle(epoch{1,1})
            if saveFig
                saveas(gcf,fullfile(saveDir,[dataset,'_',epoch{1,1},'_PopMeanTimeCourse.svg']));
            end
        end
       
        %Plot Trajectories in Target Space
        for epoch = {'reach','target hold'}
            if strcmp(epoch,'reach')
                tempTrajStruct = reachTrajStruct;
            elseif strcmp(epoch,'target hold')
                tempTrajStruct = targetHoldTrajStruct;
            end
            f = figure;
            f.Position = [20 20 800 800];
            for posture = 2
                for target = reachTargetList
                    time = tempTrajStruct([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target).avgSmoothFR.timestamps;
                    traj = tempTrajStruct([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target).targetPCA.traj;
                    numPCs = 5;
                    for i = 1:numPCs
                        subplot(5,1,i)
                       plot(time,traj(:,i),'Color',tcmap(target,:));
                       hold on
                       xlabel('t (ms)')
                       ylabel(['Target PC ',num2str(i)])
                    end
                end
            end
            sgtitle(epoch{1,1})
             if saveFig
                saveas(gcf,fullfile(saveDir,[dataset,'_',epoch{1,1},'_TargetPCTimeCourse.svg']));
            end
        end
        
        %Plot Trajectories in Posture Space
        for epoch = {'reach','target hold'}
            if strcmp(epoch,'reach')
                tempTrajStruct = reachTrajStruct;
            elseif strcmp(epoch,'target hold')
                tempTrajStruct = targetHoldTrajStruct;
            end
            f = figure;
            f.Position = [20 20 800 800];
            for posture = 2
                for target = reachTargetList
                    time = tempTrajStruct([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target).avgSmoothFR.timestamps;
                    traj = tempTrajStruct([tempTrajStruct.posture]==posture & [tempTrajStruct.target]==target).posturePCA.traj;
                    numPCs = 5;
                    for i = 1:numPCs
                        subplot(5,1,i)
                       plot(time,traj(:,i),'Color',tcmap(target,:));
                       hold on
                       xlabel('t (ms)')
                       ylabel(['Posture PC ',num2str(i)])
                    end
                end
            end
            sgtitle(epoch{1,1})
            if saveFig
                saveas(gcf,fullfile(saveDir,[dataset,'_',epoch{1,1},'_PosturePCTimeCourse.svg']));
            end
        end
        
        
        
%         %Look at posture space
%         for i = 1:size(targetHoldTrajStruct,2)
%             allFR = targetHoldTrajStruct(i).allSmoothFR;
%             allFR = vertcat(allFR.traj);
%             targetHoldTrajStruct(i).condAvg = mean(allFR);
%         end
%         allCondAvg = vertcat(targetHoldTrajStruct.condAvg);
%         [targetHoldPCA,~,~,~,explained,~] = pca(allCondAvg);
%         figure
%         for posture = 1:5
%             for target = reachTargetList
%                 traj = targetHoldTrajStruct([targetHoldTrajStruct.posture]==posture & [targetHoldTrajStruct.target]==target).condAvg;
%                 traj = traj*targetHoldPCA;
%                 plot3(traj(:,1),traj(:,2),traj(:,3),'o','MarkerSize',10,'MarkerFaceColor',pcmap(posture,:),'MarkerEdgeColor',tcmap(target,:));
%                 hold on
%             end
%         end
%         
%         %Peform LDA by posture on all data 
%         obsStruct = struct('label',[],'numObs',[],'allObs',[]);
%         structInd = 1;   
%         for posture = reachPostureList
%             for target = [1,3,4,5,6,7]
%                obsStruct(structInd).label = str2num([num2str(posture),num2str(target)]);
%                allObs = [];
%                tempTrajStruct = targetHoldTrajStruct([targetHoldTrajStruct.target]==target & [targetHoldTrajStruct.posture]==posture);
%                for i = 1:size(tempTrajStruct,2)
%                    for j = 1:size(tempTrajStruct(i).allSmoothFR,2)
%                        timestamps = tempTrajStruct(i).allSmoothFR(j).timestamps;
%                        traj = tempTrajStruct(i).allSmoothFR(j).traj;
% %                        if size(traj,1) > minNumTimestamps
% %                             traj = traj(1:minNumTimestamps,:);
% %                        end
%                        allObs = vertcat(allObs,traj);
%                    end
%                end
%                obsStruct(structInd).allObs = allObs;
%                obsStruct(structInd).numObs = size(obsStruct(structInd).allObs,1);
%                structInd = structInd + 1;
%             end
%         end
%         [postureLDA] = doLDA(obsStruct);
        
%         for i = 1:size(obsStruct,2)
%            label = num2str(obsStruct(i).label); 
%            posture = str2num(label(1));
%            target = str2num(label(2));
%            obs = mean(obsStruct(i).allObs);
%            traj = obs*postureLDA;
%            plot3(traj(:,1),traj(:,2),traj(:,3),'o','MarkerSize',10,'MarkerFaceColor',pcmap(posture,:),'MarkerEdgeColor',tcmap(target,:));
%            hold on
%         end
        
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
         