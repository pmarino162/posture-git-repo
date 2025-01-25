clear; clc; clf; close all;
    
%% Parameters
    numManPCs = 10;
%% Setup figure save
    saveFig = false;
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Presentations\20220121';
    
%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customSummer.mat');
    cmap = customSummer;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat');
    tcmap = customRainbow;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\haystacks.mat')
    pcmap = haystacks;
    
%% Load Data, subselect, get traj struct
    dataset = 'E20200316';
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
    end
    trajStruct = getTrajStruct20220419(Data,condFields,trajFields,trialInclStates,binWidth,kernelStdDev,'matchConditions',false,'zScoreParams',zScoreParams);


  %Get minimum number of timestamps in condition averages
    numTimestamps = [];
    for i = 1:size(trajStruct,2)
        numTimestamps(i) = length(trajStruct(i).avgSmoothFR.timestamps);
    end
    [minNumTimestamps,i] = min(numTimestamps);
    numPts = 10;    

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
    
%% Perform marginalizations and store them
    %Xdims: 1=time, 2=target, 3=posture, 4=channel
    
    %Condition-invariant
    CIMargOffset = mean(X,[1 2 3],'omitnan');
    CIMargTraj = mean(X,[2 3],'omitnan') - CIMargOffset;
    CIErrBar = NaN(numPts,numManPCs);
    for i = 1:numPts
        tempX = X(i,:,:,:)-CIMargOffset;
        tempX = reshape(tempX,[numTargets*numPostures,numManPCs]);
        CIErrBar(i,:) = 1.96.*std(tempX)./sqrt(size(tempX,1));
    end

    %Posture and Target Offsets
    targetMargOffset = mean(X,[1 3],'omitnan') - CIMargOffset;
    postureMargOffset = mean(X,[1 2],'omitnan') - CIMargOffset;
    
    %Posture and Target Traj
    targetMargTraj = mean(X,[3],'omitnan') - CIMargTraj - targetMargOffset - CIMargOffset;
    postureMargTraj = mean(X,[2],'omitnan') - CIMargTraj - postureMargOffset - CIMargOffset;
    
    targetErrBar = NaN(numPts,numTargets,numManPCs);
    for i = 1:numPts
        for targetInd = 1:numTargets
            tempX = X(i,targetInd,:,:)-CIMargOffset;
            tempX = reshape(tempX,[numPostures,numManPCs]);
            targetErrBar(i,targetInd,:) = 1.96.*std(tempX)./sqrt(size(tempX,1));
        end
    end
    
    postureErrBar = NaN(numPts,numPostures,numManPCs);
    for i = 1:numPts
        for postureInd = 1:numPostures
            tempX = X(i,:,postureInd,:)-CIMargOffset;
            tempX = reshape(tempX,[numTargets,numManPCs]);
            postureErrBar(i,postureInd,:) = 1.96.*std(tempX)./sqrt(size(tempX,1));
        end
    end
    
    %Interaction
    intMargOffset = mean(X,[1],'omitnan') - targetMargOffset - postureMargOffset - CIMargOffset ;
    intMargTraj = X-intMargOffset-targetMargTraj-postureMargTraj-targetMargOffset-postureMargOffset-CIMargTraj-CIMargOffset;
    
%% Get VAFs of marginalizations
    %Get total var
    totalVar = 0;
    numObs = 0;
    postureInd = 1;
    for posture = postureList
        targetInd = 1;
        for target = targetList
            for i = 1:numPts
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
            for i = 1:numPts
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
            for i = 1:numPts
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
            for i = 1:numPts
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
            for i = 1:numPts
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
            for i = 1:numPts
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
            for i = 1:numPts
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
            for i = 1:numPts
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
    
%% Pie chart
    figure
    labels = {['Condition Invariant ',num2str(round(CIVarPct)),'%'],['Target ',num2str(round(targetVarPct)),'%'],['Posture ',num2str(round(postureVarPct)),'%'],['PxT Interaction ',num2str(round(intVarPct)),'%']};
    pie([CIVarPct,targetVarPct,postureVarPct,intVarPct],[1 1 1 1],labels)
    if saveFig
       saveas(gcf,fullfile(saveDir,task,'Component Pie.jpg')); 
    end
    
    figure
    labels = {['Constant ',num2str(round(postureOffsetVarPct)),'%'],['Time Varying ',num2str(round(100-postureOffsetVarPct)),'%']};
    pie([postureOffsetVarPct,100-postureOffsetVarPct],[1 1],labels)
    if saveFig
       saveas(gcf,fullfile(saveDir,task,'Posture Pie.jpg')); 
    end
    close all
%% Visualize by projecting marginalizations into PC's
    
    time = trajStruct(1).avgSmoothFR.timestamps(1:numPts);
    
    %Condition-invariant
    f = figure; f.Position = [50 50 1500 700];
    traj = squeeze(CIMargTraj);
    for i = 1:10
        subplot(2,5,i)
        shadedErrorBar(time,traj(:,i),CIErrBar(:,i),'lineprops',{'LineWidth',2})
%         plot(time,traj(:,i))
        ax = gca;
        ax.YLim = [-4 4];
        xlabel('time (ms)')
        ylabel(['PC ',num2str(i)])
    end
    sgtitle('Time Marginalization')
    if saveFig
        saveas(gcf,fullfile(saveDir,task,'TimeMarginalization.jpg')); 
    end

    
    %Target
    f = figure; f.Position = [50 50 1500 700];
    allTraj = targetMargTraj+targetMargOffset;
    targetInd = 1;
    for target = targetList
       traj = squeeze(allTraj(:,targetInd,:,:));
       for i = 1:10
            subplot(2,5,i)
%             plot(time,traj(:,i),'Color',tcmap(target,:))
            shadedErrorBar(time,traj(:,i),targetErrBar(:,targetInd,i),'lineprops',{'LineWidth',2,'Color',tcmap(target,:)})
            ax = gca;
            ax.YLim = [-4 4];
            hold on;
            xlabel('time (ms)')
            ylabel(['PC ',num2str(i)])
       end
       targetInd = targetInd + 1;
    end
    sgtitle('Target Marginalization')
    if saveFig
        saveas(gcf,fullfile(saveDir,task,'TargetMarginalization.jpg')); 
    end
    
    %Posture
    f = figure; f.Position = [50 50 1500 700];
    allTraj = postureMargTraj+postureMargOffset;
    postureInd = 1;
    for posture = postureList
       traj = squeeze(allTraj(:,:,postureInd,:));
       for i = 1:10
            subplot(2,5,i)
%             plot(time,traj(:,i),'Color',pcmap(posture,:))
            shadedErrorBar(time,traj(:,i),postureErrBar(:,postureInd,i),'lineprops',{'LineWidth',2,'Color',pcmap(posture,:)})
            ax = gca;
            ax.YLim = [-4 4];
            hold on;
            xlabel('time (ms)')
            ylabel(['PC ',num2str(i)])
       end
       postureInd = postureInd + 1;
    end
    sgtitle('Posture Marginalization')
    if saveFig
        saveas(gcf,fullfile(saveDir,task,'PostureMarginalization.jpg')); 
    end
    
    %Interaction
    allTraj = intMargTraj+intMargOffset;
    postureInd = 1;
    for posture = postureList
        f = figure; f.Position = [50 50 1500 700];
        allTraj = intMargTraj+intMargOffset;
        targetInd = 1;
        for target = targetList
            traj = squeeze(allTraj(:,targetInd,postureInd,:));
            for i = 1:10
                subplot(2,5,i)
                plot(time,traj(:,i),'Color',tcmap(target,:))
                hold on;
                ax = gca;
                ax.YLim = [-4 4];
                hold on;
                xlabel('time (ms)')
                ylabel(['PC ',num2str(i)])
            end
            targetInd = targetInd + 1;
        end
        postureInd = postureInd + 1;
        sgtitle(['Int. Marg. Posture ',num2str(posture)])
        if saveFig
            saveas(gcf,fullfile(saveDir,task,['IntMargP',num2str(posture),'.jpg'])); 
        end

    end
    
% %% Estimate dimensionality of each marginalization
% 
