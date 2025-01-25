function [Data] = subselectForTrajDist(Data,dataType)

%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customSummer.mat');
    cmap = customSummer;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat');
    tcmap = customRainbow;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\rainbow2d.mat')
    pcmap = rainbow2d;
    
switch dataType
    %% BCI
    case 'BCI'
        %Reach time histogram (Exclude +/-1*std)
        kinData = [Data.kinData];
        reachTime = [kinData.movementTime];
        meanReachTime = mean(reachTime);
        stdReachTime = std(reachTime);
        figure; histogram(reachTime);
        hold on
        ax = gca;
        line([meanReachTime meanReachTime],[ax.YLim],'LineWidth',2,'Color','b')
        hold on
        acceptDist = .75*stdReachTime;
        line([meanReachTime+acceptDist meanReachTime+acceptDist],[ax.YLim],'LineWidth',2,'Color','r')
        line([meanReachTime-acceptDist meanReachTime-acceptDist],[ax.YLim],'LineWidth',2,'Color','r')
        numTrials = size(Data,2);
        rmTrials = reachTime < meanReachTime-acceptDist | reachTime > meanReachTime+acceptDist; 
        numTrialsRemaining = sum(numTrials - sum(rmTrials));
        text(0.95*ax.XLim(2),0.95*ax.YLim(2),[num2str(numTrialsRemaining),'/',num2str(numTrials),...
            ' (',num2str(round(100*numTrialsRemaining/numTrials)),'%) ',' remaining'],'HorizontalAlignment','right');
        xlabel('Reach Time (ms)')
        ylabel('Count')
        
        %Setup condFields
        if strcmpi(Data(1).trialName,'BCI Center Out')
            condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
        else strcmpi(Data(1).trialName,'GridTask_CO_Across_BC_ForceBar')
            condFields = {{'target','targetData','target1ID'},{'posture','conditionData','postureID'}};
        end
        
        %Compare effect of reach time cutoffs for each condition
        kinFields = {'movementTime'};
        [kinStruct] = getKinStruct(Data,condFields,kinFields);
        figure
        hold on
        for posture = 1:7
           for target = 1:8 
               if any([kinStruct.posture]==posture & [kinStruct.target]==target)
                    reachTimes = kinStruct([kinStruct.posture]==posture & [kinStruct.target]==target).allMovementTime;
                    pct25 = prctile(reachTimes,25);
                    pct75 = prctile(reachTimes,75);
                    med = median(reachTimes);
                    jitter = 0.7*rand(1,1)-0.35;
                    plot(med,posture+jitter,'.','MarkerSize',20,'Color',tcmap(target,:));
                    line([pct25 pct75],[posture+jitter posture+jitter],'Color',tcmap(target,:))
                    ax = gca;
                    line([meanReachTime+acceptDist meanReachTime+acceptDist],[ax.YLim],'LineWidth',2,'Color','r')
                    line([meanReachTime-acceptDist meanReachTime-acceptDist],[ax.YLim],'LineWidth',2,'Color','r')
                end
           end
        end
        
        %Exclude trials based on reach time
        Data(rmTrials) = [];
        
        %Examine remaining condition counts
        trialCounts = countConditionTrials(Data,condFields);
        numCondTrials = [trialCounts.numTrials];
        figure
        histogram(numCondTrials)
        xlabel('Number of Trials')
        ylabel('Number of Conditions')
        
        %Remove conditions based on condition counts
        rmCond = [trialCounts.numTrials] < 7;
        rmTrialNums = [trialCounts(rmCond).trialNum];
        trialNum = [Data.trialNum];
        Data(ismember(trialNum,rmTrialNums)) = [];
        close all
     
        
    %% Isometric force
    case 'iso'
        %Reach time histogram (Exclude +/-1.5*std)
        kinData = [Data.kinData];
        reachTime = [kinData.targetAcqTime];
        meanReachTime = mean(reachTime);
        stdReachTime = std(reachTime);
        figure; histogram(reachTime);
        hold on
        ax = gca;
        line([meanReachTime meanReachTime],[ax.YLim],'LineWidth',2,'Color','b')
        hold on
        acceptDist = 1.5*stdReachTime;
        line([meanReachTime+acceptDist meanReachTime+acceptDist],[ax.YLim],'LineWidth',2,'Color','r')
        line([meanReachTime-acceptDist meanReachTime-acceptDist],[ax.YLim],'LineWidth',2,'Color','r')
        numTrials = size(Data,2);
        rmTrials = reachTime < meanReachTime-acceptDist | reachTime > meanReachTime+acceptDist; 
        numTrialsRemaining = sum(numTrials - sum(rmTrials));
        text(0.95*ax.XLim(2),0.95*ax.YLim(2),[num2str(numTrialsRemaining),'/',num2str(numTrials),...
            ' (',num2str(round(100*numTrialsRemaining/numTrials)),'%) ',' remaining'],'HorizontalAlignment','right');
        xlabel('Reach Time (ms)')
        ylabel('Count')
        
        %Compare effect of reach time cutoffs for each condition
        condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
        kinFields = {'targetAcqTime'};
        [kinStruct] = getKinStruct(Data,condFields,kinFields);
        figure
        hold on
        for posture = 1:7
           for target = 1:8 
               if any([kinStruct.posture]==posture & [kinStruct.target]==target)
                    reachTimes = kinStruct([kinStruct.posture]==posture & [kinStruct.target]==target).allTargetAcqTime;
                    pct25 = prctile(reachTimes,25);
                    pct75 = prctile(reachTimes,75);
                    med = median(reachTimes);
                    jitter = 0.7*rand(1,1)-0.35;
                    plot(med,posture+jitter,'.','MarkerSize',20,'Color',tcmap(target,:));
                    line([pct25 pct75],[posture+jitter posture+jitter],'Color',tcmap(target,:))
                    ax = gca;
                    line([meanReachTime+acceptDist meanReachTime+acceptDist],[ax.YLim],'LineWidth',2,'Color','r')
                    line([meanReachTime-acceptDist meanReachTime-acceptDist],[ax.YLim],'LineWidth',2,'Color','r')
                end
           end
        end
        
        %Exclude trials based on reach time
        Data(rmTrials) = [];
        
        %Examine remaining condition counts
        condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
        trialCounts = countConditionTrials(Data,condFields);
        numCondTrials = [trialCounts.numTrials];
        figure
        histogram(numCondTrials)
        xlabel('Number of Trials')
        ylabel('Number of Conditions')
        
        %Remove conditions based on condition counts
        rmCond = [trialCounts.numTrials] < 14;
        rmTrialNums = [trialCounts(rmCond).trialNum];
        trialNum = [Data.trialNum];
        Data(ismember(trialNum,rmTrialNums)) = [];
        close all
        
    %Grid Reaching
    case 'reaching'
        %Reach time histogram (Exclude +/-1.5*std)
        kinData = [Data.kinData];
        reachTime = [kinData.reachTime];
        meanReachTime = mean(reachTime);
        stdReachTime = std(reachTime);
        figure; histogram(reachTime);
        hold on
        ax = gca;
        line([meanReachTime meanReachTime],[ax.YLim],'LineWidth',2,'Color','b')
        hold on
        acceptDist = 1.5*stdReachTime;
        line([meanReachTime+acceptDist meanReachTime+acceptDist],[ax.YLim],'LineWidth',2,'Color','r')
        line([meanReachTime-acceptDist meanReachTime-acceptDist],[ax.YLim],'LineWidth',2,'Color','r')
        numTrials = size(Data,2);
        rmTrials = reachTime < meanReachTime-acceptDist | reachTime > meanReachTime+acceptDist; 
        numTrialsRemaining = sum(numTrials - sum(rmTrials));
        text(0.95*ax.XLim(2),0.95*ax.YLim(2),[num2str(numTrialsRemaining),'/',num2str(numTrials),...
            ' (',num2str(round(100*numTrialsRemaining/numTrials)),'%) ',' remaining'],'HorizontalAlignment','right');
        xlabel('Reach Time (ms)')
        ylabel('Count')
        
        %Compare effect of reach time cutoffs for each condition
        condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
        kinFields = {'reachTime'};
        [kinStruct] = getKinStruct(Data,condFields,kinFields);
        figure
        hold on
        for posture = 1:7
           for target = 1:8 
               if any([kinStruct.posture]==posture & [kinStruct.target]==target)
                    reachTimes = kinStruct([kinStruct.posture]==posture & [kinStruct.target]==target).allReachTime;
                    pct25 = prctile(reachTimes,25);
                    pct75 = prctile(reachTimes,75);
                    med = median(reachTimes);
                    jitter = 0.7*rand(1,1)-0.35;
                    plot(med,posture+jitter,'.','MarkerSize',20,'Color',tcmap(target,:));
                    line([pct25 pct75],[posture+jitter posture+jitter],'Color',tcmap(target,:))
                    ax = gca;
                    line([meanReachTime+acceptDist meanReachTime+acceptDist],[ax.YLim],'LineWidth',2,'Color','r')
                    line([meanReachTime-acceptDist meanReachTime-acceptDist],[ax.YLim],'LineWidth',2,'Color','r')
                end
           end
        end
        
        %Exclude trials based on reach time
        Data(rmTrials) = [];
        
        %Examine remaining condition counts
        condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
        trialCounts = countConditionTrials(Data,condFields);
        numCondTrials = [trialCounts.numTrials];
        figure
        histogram(numCondTrials)
        xlabel('Number of Trials')
        ylabel('Number of Conditions')
        
        %Remove conditions based on condition counts
        rmCond = [trialCounts.numTrials] < 10;
        rmTrialNums = [trialCounts(rmCond).trialNum];
        trialNum = [Data.trialNum];
        Data(ismember(trialNum,rmTrialNums)) = [];
        close all
        
    %Grid Planning
    case 'planning'
        %Keep only trials with delay length >= 500ms
        kinData = [Data.kinData];
        delayLength = [kinData.delayLength];
        rmTrials = delayLength < 500;
        Data(rmTrials) = [];
        
        %Examine remaining condition counts
        condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
        trialCounts = countConditionTrials(Data,condFields);
        numCondTrials = [trialCounts.numTrials];
        figure
        histogram(numCondTrials)
        xlabel('Number of Trials')
        ylabel('Number of Conditions')
        
        %Remove conditions based on condition counts
        rmCond = [trialCounts.numTrials] < 8;
        rmTrialNums = [trialCounts(rmCond).trialNum];
        trialNum = [Data.trialNum];
        Data(ismember(trialNum,rmTrialNums)) = [];
        close all
end

    
% %% Look at trials
%     trialInd = min(find(step1AcqTime>600));
%     cursorTraj = Data(trialInd).Decoder.cursorTraj;
%     timestamps = Data(trialInd).Decoder.timestamps;
%     stateTransitions = Data(trialInd).stateData.stateTransitions;
%     step1Time = stateTransitions(2,find(stateTransitions(1,:)==4));
%     step2Time = stateTransitions(2,find(stateTransitions(1,:)==6));
%     
%     figure
%     subplot(2,1,1)
%     plot(timestamps,cursorTraj(:,1))
%     subplot(2,1,2)
%     plot(timestamps,cursorTraj(:,2))
end