clear; clc; clf; close all;

%% Setup saveFig   
    saveFig = false;
    %saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Presentations\20220906 - compositionality analysis';
    %saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Presentations\20220906 - compositionality analysis\compositionality';
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Paper\20220907\Figure 8';
    set(0, 'DefaultFigureRenderer', 'painters');
    
%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\orli.mat')
    pcmap = orli;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat')
    tcmap = customRainbow;

%% Get trajStruct
%'N20190226
%'R20200221'
    dataset = 'N20190226'; task = 'Reach';
    [Data,zScoreParams] = loadData(dataset);
    binWidth = 25; kernelStdDev = 25;
    trajFields = {'zSmoothFR'};
    trialInclStates = struct('trialName','','inclStates',[]);
    switch dataset
        %Earl Reach
        case {'E20210706','E20210707','E20210708','E20210709','E20210710'}
            trialInclStates(1).trialName = {'GridReaching'};
            condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
            trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'kin','moveOnsetTime','first',0}};
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
%         case 'N20190222'
%             trajStruct([trajStruct.target]==6) = [];
%             trajStruct([trajStruct.target]==7) = [];
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
            %trajectory, then averaging  averaging to get CI1avg, then
            %taking C12avg on C1 component from other posture data
            CI1mat = zeros(minNumTimestamps,numTargets-1,numChannels);
            targetInd = 1;
            for tempTarget = setdiff(targetList,target)  
               CI1mat(:,targetInd,:) = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==tempTarget).avgSmoothFR.traj; 
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
   
   %For each model, get R2
   ssStruct = struct('target',[],'posture',[],'SSres',[],'SStot',[],'R2',[]);
   structInd = 1;
   SStot = 0; SSres = 0;
   for posture = postureList
       for target = targetList
           SSresTemp = 0; SStotTemp = 0;
           actualTraj = actualTrajStruct([actualTrajStruct.posture]==posture & [actualTrajStruct.target]==target).traj;
           fullTraj = fullTrajStruct([fullTrajStruct.posture]==posture & [fullTrajStruct.target]==target).traj;
           for i = 1:minNumTimestamps
               SStot = (actualTraj(i,:)-grandMean)*(actualTraj(i,:)-grandMean)' + SStot;
               SSres = (actualTraj(i,:)-fullTraj(i,:))*(actualTraj(i,:)-fullTraj(i,:))' + SSres;
               SStotTemp = (actualTraj(i,:)-grandMean)*(actualTraj(i,:)-grandMean)' + SStotTemp;
               SSresTemp = (actualTraj(i,:)-fullTraj(i,:))*(actualTraj(i,:)-fullTraj(i,:))' + SSresTemp;
           end
           ssStruct(structInd).target = target;
           ssStruct(structInd).posture = posture;
           ssStruct(structInd).SSres = SSresTemp;
           ssStruct(structInd).SStot = SStotTemp;
           ssStruct(structInd).R2 = 1-(SSresTemp/SStotTemp);
           structInd = structInd + 1;
       end
   end
   R2 = 1-(SSres/SStot); 


%% Store results
   resultStruct(1).dataset = dataset; 
   resultStruct(1).task = task;
   
   %R2
   resultStruct(1).R2 = R2;
   
   %Actual and predicted trajectories
   resultStruct(1).actualTrajStruct = actualTrajStruct;
   resultStruct(1).predTrajStruct = fullTrajStruct;
   

%% Visualize results
    %Do PCA on model and actual results
    allActualTraj = vertcat(actualTrajStruct.traj);
    allFullTraj = vertcat(fullTrajStruct.traj);
    allTraj = vertcat(allActualTraj,allFullTraj);
    [PCs,score,latent,tsquared,explained,mu] = pca(allTraj); 
    
    %Add Projections to trajStruct
    for i = 1:size(actualTrajStruct,2)
        %Add average traces to trajStruct
        actualTrajStruct(i).PC = actualTrajStruct(i).traj*PCs;
        fullTrajStruct(i).PC = fullTrajStruct(i).traj*PCs;
    end  
    
    %Plot time courses
    for posture = postureList
        for target = targetList
            maxYLim = -100; minYLim = 100;
            f = figure;
            f.Position = [100 100 1000 500];
            if any([trajStruct.posture]==posture & [trajStruct.target]==target)
                actualTraj = actualTrajStruct([actualTrajStruct.posture]==posture & [actualTrajStruct.target]==target).PC; 
                predTraj = fullTrajStruct([fullTrajStruct.posture]==posture & [fullTrajStruct.target]==target).PC; 
                if numPCs < 10
                    numDims = numPCs;
                else
                    numDims = 10;
                end
                for dim = 1:numDims
                    subplot(2,5,dim);hold on
                    plot(actualTraj(:,dim),'Color',pcmap(posture,:),'LineWidth',2);
                    plot(predTraj(:,dim),'--','Color',pcmap(posture,:),'LineWidth',2);
                    yline(grandMean*PCs(:,dim));
                    ax = gca;
                    ylim = ax.YLim;
                    if ylim(1) < minYLim 
                        minYLim = ylim(1);
                    end
                    if ylim(2) > maxYLim 
                        maxYLim = ylim(2);
                    end

                end
                for dim = 1:numDims
                    subplot(2,5,dim)
                    if dim == 1
                        legend('Actual','Predicted')
                    end
                   ax = gca;
                   ax.YLim = [minYLim maxYLim];
                   xlabel('time (ms)')
                   ylabel(['PC ',num2str(dim)])
                end
                R2 = ssStruct([ssStruct.posture]==posture & [ssStruct.target]==target).R2;
                sgtitle(['P',num2str(posture),' T',num2str(target),' R2=',num2str(R2)])
            end
        end
    end
    
    
    R2 = resultStruct.R2
        