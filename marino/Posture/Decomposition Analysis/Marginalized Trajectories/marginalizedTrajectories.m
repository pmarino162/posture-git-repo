clear; clc; clf; close all;
    
%% Setup figure save
    saveFig = true;
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Presentations\20220121';
    
%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customSummer.mat');
    cmap = customSummer;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat');
    tcmap = customRainbow;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\haystacks.mat')
    pcmap = haystacks;
    
%% Load Data, subselect, get traj struct
    trajFields = {'smoothFR'};
    trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
    binWidth = 25;

    % BCI
%     task = 'BCI';
%     [Data,~,~,~,~,~,~,~,~,~,~,~] = loadEarlData20200317_20211210;
%     [Data] = subselectForMarginalization(Data,task);
%     trialInclStates(1).trialName = {'GridTask_CO_Across_BC_ForceBar'};
%     condFields = {{'target','targetData','target1ID'},{'posture','conditionData','postureID'}};
%     trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-150},{'state','Step 2','first',0}};
%     trajStruct = getTrajStruct20211210(Data,condFields,trajFields,trialInclStates,binWidth,'matchConditions',false);

%     %Reaching
%     task = 'Reaching';
%     [Data] = loadEarlData20210706_20211210; 
%     [Data] = subselectForMarginalization(Data,'reaching');
%     trialInclStates(1).trialName = {'GridReaching'};
%     condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
%     trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-200},{'state','Target Hold','first',0}};
%     trajStruct = getTrajStruct20211210(Data,condFields,trajFields,trialInclStates,binWidth,'matchConditions',false);
%     
%     
%     %Planning
%     [Data] = loadEarlData20210706_20211210; 
%     [Data] = subselectForMarginalization(Data,'planning');
%     task = 'Planning';
%     trialInclStates(1).trialName = {'GridReaching'};
%     condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
%     trialInclStates(1).inclStates = {{'state','Delay','first',0},{'state','Delay','first',250}};
%     trajStruct = getTrajStruct20211210(Data,condFields,trajFields,trialInclStates,binWidth,'matchConditions',false);
%     
    %Iso
    task = 'Iso';
        [Data] = loadEarlData20200116_20211210();
    [Data] = subselectForMarginalization(Data,'iso');
    trialInclStates(1).trialName = {'IsometricForce_1D'};
    condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
    trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-250},{'state','Target Hold','first',0}};
    trajStruct = getTrajStruct20211210(Data,condFields,trajFields,trialInclStates,binWidth,'matchConditions',false);
    
%% Get Postures, Targets, and Channels
    postureList = unique([trajStruct.posture]);
    targetList = unique([trajStruct.target]);
    numPostures = size(postureList,2);
    numTargets = size(targetList,2);
    numChannels = size(trajStruct(1).avgSmoothFR.traj,2);
    numConditions = size(trajStruct,2);
    
%% Get timestamps dist and min number timestamps
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
    [minNumCondTrials,i] = min(numCondTraj);
    [maxNumCondTrials,i] = max(numCondTraj);

%% Do PCA on condition averages
    allCondAvg = NaN(numConditions,numChannels);
    i = 1;
    for cond = 1:numConditions
            traj = trajStruct(cond).avgSmoothFR.traj;
            allCondAvg(i:i+minNumTimestamps-1,:) = traj(1:minNumTimestamps,:);
            i = i + minNumTimestamps;
    end
    [coeff,score,latent,tsquared,explained,mu] = pca(allCondAvg); 
    PCs = coeff(:,1:15);
    
%% Form X, containing trial-averaged data for each condition
    X = NaN(minNumTimestamps,numTargets,numPostures,maxNumCondTrials,numChannels);
    postureInd = 1;
    for posture = postureList
        targetInd = 1;
        for target = targetList
            numTrials = size(trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).allSmoothFR,2);
            for trial = 1:numTrials
                traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).allSmoothFR(trial).traj;
                X(:,targetInd,postureInd,trial,:) = traj(1:minNumTimestamps,:); 
            end
            targetInd = targetInd + 1;
        end
        postureInd = postureInd+1;
    end
    X = squeeze(mean(X,4,'omitnan'));
    %Clear Data from memory
    clearvars Data
    
%% Perform marginalizations and store them
    %Xdims: 1=time, 2=target, 3=posture, 4=channel
    
    %Condition-invariant
    CIMargOffset = mean(X,[1 2 3],'omitnan');
    CIMargTraj = mean(X,[2 3],'omitnan') - CIMargOffset;
    CIErrBar = NaN(minNumTimestamps,size(PCs,2));
    for i = 1:minNumTimestamps
        tempX = X(i,:,:,:)-CIMargOffset;
        tempX = reshape(tempX,[numTargets*numPostures,numChannels]);
        tempPC = tempX*PCs;
        CIErrBar(i,:) = 1.96.*std(tempPC)./sqrt(size(tempPC,1));
    end

    %Posture and Target Offsets
    targetMargOffset = mean(X,[1 3],'omitnan') - CIMargOffset;
    postureMargOffset = mean(X,[1 2],'omitnan') - CIMargOffset;
    
    %Posture and Target Traj
    targetMargTraj = mean(X,[3],'omitnan') - CIMargTraj - targetMargOffset - CIMargOffset;
    postureMargTraj = mean(X,[2],'omitnan') - CIMargTraj - postureMargOffset - CIMargOffset;
    
    targetErrBar = NaN(minNumTimestamps,numTargets,size(PCs,2));
    for i = 1:minNumTimestamps
        for targetInd = 1:numTargets
            tempX = X(i,targetInd,:,:)-CIMargOffset;
            tempX = reshape(tempX,[numPostures,numChannels]);
            tempPC = tempX*PCs;
            targetErrBar(i,targetInd,:) = 1.96.*std(tempPC)./sqrt(size(tempPC,1));
        end
    end
    
    postureErrBar = NaN(minNumTimestamps,numPostures,size(PCs,2));
    for i = 1:minNumTimestamps
        for postureInd = 1:numPostures
            tempX = X(i,:,postureInd,:)-CIMargOffset;
            tempX = reshape(tempX,[numTargets,numChannels]);
            tempPC = tempX*PCs;
            postureErrBar(i,postureInd,:) = 1.96.*std(tempPC)./sqrt(size(tempPC,1));
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
    
    time = trajStruct(1).avgSmoothFR.timestamps(1:minNumTimestamps);
    
    %Condition-invariant
    f = figure; f.Position = [50 50 1500 700];
    traj = squeeze(CIMargTraj)*PCs;
    for i = 1:15
        subplot(3,5,i)
        shadedErrorBar(time,traj(:,i),CIErrBar(:,i),'lineprops',{'LineWidth',2})
%         plot(time,traj(:,i))
        ax = gca;
        ax.YLim = [-200 200];
        xlabel('time (ms)')
        ylabel(['PC ',num2str(i)])
    end
    sgtitle('Time Marginalization')
    if saveFig
        saveas(gcf,fullfile(saveDir,task,'TimeMarginalization.jpg')); 
    end
    close all
    
    %Target
    f = figure; f.Position = [50 50 1500 700];
    allTraj = targetMargTraj+targetMargOffset;
    targetInd = 1;
    for target = targetList
       traj = squeeze(allTraj(:,targetInd,:,:))*PCs;
       for i = 1:15
            subplot(3,5,i)
%             plot(time,traj(:,i),'Color',tcmap(target,:))
            shadedErrorBar(time,traj(:,i),targetErrBar(:,targetInd,i),'lineprops',{'LineWidth',2,'Color',tcmap(target,:)})
            ax = gca;
            ax.YLim = [-250 250];
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
    close all
    
    %Posture
    f = figure; f.Position = [50 50 1500 700];
    allTraj = postureMargTraj+postureMargOffset;
    postureInd = 1;
    for posture = postureList
       traj = squeeze(allTraj(:,:,postureInd,:))*PCs;
       for i = 1:15
            subplot(3,5,i)
%             plot(time,traj(:,i),'Color',pcmap(posture,:))
            shadedErrorBar(time,traj(:,i),postureErrBar(:,postureInd,i),'lineprops',{'LineWidth',2,'Color',pcmap(posture,:)})
            ax = gca;
            ax.YLim = [-250 250];
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
    close all
    
    %Interaction
    allTraj = intMargTraj+intMargOffset;
    postureInd = 1;
    for posture = postureList
        f = figure; f.Position = [50 50 1500 700];
        allTraj = intMargTraj+intMargOffset;
        targetInd = 1;
        for target = targetList
            traj = squeeze(allTraj(:,targetInd,postureInd,:))*PCs;
            for i = 1:15
                subplot(3,5,i)
                plot(time,traj(:,i),'Color',tcmap(target,:))
                hold on;
                ax = gca;
                ax.YLim = [-250 250];
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
        close all
    end
    
% %% Estimate dimensionality of each marginalization
% 
%% Do dPCA
    %Xdims: 1=time, 2=target, 3=posture, 4=channel
    Xdpca = permute(X,[4,3,2,1]);
    
%     Xdpca(:,2:4,:,:) = [];
%     Xdpca(:,:,[2:4,6:8],:) = [];
%     Xdpca(:,:,[2,4,6,8],:) = [];
    %Xdpcadims: 1=neurons, 2=posture, 3=target, 4=time
    
    combinedParams = {{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}};
    margNames = {'P', 'T', 'CI', 'PTI'};
    margColours = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256;

    N = numChannels;
    P = numPostures;
    T = numTargets;
    timePts = minNumTimestamps;
    time = trajStruct(1).avgSmoothFR.timestamps(1:minNumTimestamps);
    timeEvents = [];
    
[W,V,whichMarg] = dpca(Xdpca, 10, ...
    'combinedParams', combinedParams);

explVar = dpca_explainedVariance(Xdpca, W, V, ...
    'combinedParams', combinedParams);

dpca_plot(Xdpca, W, V, @dpca_plot_default, ...
    'explainedVar', explVar, ...
    'marginalizationNames', margNames, ...
    'marginalizationColours', margColours, ...
    'whichMarg', whichMarg,                 ...
    'time', time,                        ...
    'timeEvents', timeEvents,               ...
    'timeMarginalization', 3, ...
    'legendSubplot', 1);

    if saveFig
        saveas(gcf,fullfile(saveDir,task,['dPCA.jpg'])); 
    end
    close all
    
%% Use dPCA decomposition to plot 
    figure
    postureInd = 1;
    for posture = postureList
        targetInd = 1;
        for target = targetList
            traj = squeeze(Xdpca(:,postureInd,targetInd,:));
            traj = traj'*V;
            for i = 1:10
                subplot(5,2,i)
                plot(time,traj(:,i),'Color',tcmap(target,:))
                hold on
            end
            targetInd = targetInd + 1;
        end
        postureInd = postureInd + 1;
    end
    

    figure
    postureInd = 1;
    for posture = postureList
        targetInd = 1;
        for target = targetList
            traj = squeeze(Xdpca(:,postureInd,targetInd,:));
            traj = traj'*V;
            for i = 1:10
                subplot(5,2,i)
                plot(time,traj(:,i),'Color',pcmap(posture,:))
                hold on
            end
            targetInd = targetInd + 1;
        end
        postureInd = postureInd + 1;
    end
    
    PTOrth = V(:,[1,2,5]);
    PTOrth = qr(PTOrth);
    
        figure
    postureInd = 1;
    for posture = postureList
        targetInd = 1;
        for target = targetList
            traj = squeeze(Xdpca(:,postureInd,targetInd,:));
            traj = traj'*PTOrth;
            plot3(traj(:,1),traj(:,2),traj(:,3),'Color',pcmap(posture,:),'LineWidth',2)
            hold on
            targetInd = targetInd + 1;
        end
        postureInd = postureInd + 1;
    end
    
