clear; clc; clf; close all;

%% Load Data 
    %Nigel - 2/21/2018
    date = '20180221';
    Data = loadNigelData(date);
 
%% Get rid of channels that were always off (Nigel)
    zeroList = 1:384;
    for trial = 1:size(Data,2)
        allSortSpikeBins = Data(trial).spikes.allSortSpikeBins;
        trialZeroList = [];
        for sort = 1:384
            if sum(allSortSpikeBins(:,sort))==0
                trialZeroList = [trialZeroList,sort];
            end
        end
        zeroList = intersect(zeroList,trialZeroList);
    end

    for trial = 1:size(Data,2)
        allSortSpikeBins = Data(trial).spikes.allSortSpikeBins;
        allSortSpikeBins(:,zeroList) = [];
        Data(trial).spikes.allSortSpikeBins = allSortSpikeBins;
    end

%% Collect States of Interest into Traj Field
    trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);   
    %Nigel 
    trialInclStates(1).trialName = {'Nigel Posture BC Center Out'};
    trialInclStates(1).inclStates = {'Cursor Release','Center Exit','Target Hold'};
    trialInclStates(1).inclOccurrence = {'last','last','last'};
    trialInclStates(1).addTimeToBeginning = {-50,0,0};
    numTrials = size(Data,2);
    for trial = 1:numTrials
        [FR,FRTimestamps] = getStatesTraj(Data(trial),trialInclStates,'allSortSmoothedFR');
        Data(trial).Traj.FR = FR; Data(trial).Traj.FRTimestamps = FRTimestamps;
    end
    clearvars trialInclStates
        
%% Keep only trials with length w/in 2 stdDevs of mean length
    [Data] = excludeLengths(Data,'FR');
    
%% Create trajStruct
    %Nigel
    condFields = {{'posture','postureData','postureID'},{'target','targetData','targetID'}};
    trajFields = {'FR'};
    trajStruct = getTrajStruct(Data,condFields,trajFields,'addTimeToBeginning',-50);
    
%% Channel Grouping Struct
    channelGroups = struct('channels',[]);
    channelGroups(1).channels = 1:50; channelGroups(2).channels = 51:100; channelGroups(3).channels = 101:115;

%% Target Pairs Struct
    structInd = 1;
    targetPairs = struct('pair',[]);
    for target = 1:8
        %Get 90 degree target 
        t90 = target + 2;
        if t90 > 8
            t90 = t90-8;       
        end
        targetPairs(structInd).pair = [target,t90];
        structInd = structInd + 1;
        if target < 5
            %Get 180 degree target
            t180 = target + 4;
            if t180 > 8
                t180 = t180-8;        
            end
            targetPairs(structInd).pair = [target,t180];
            structInd = structInd + 1;
        end
    end
    numTargetPairs = size(targetPairs,2);
    
%% Colormap
    cmap = parula(5);
    cmap = cmap(2:4,:);

%% Plotting - try tight subplot     
close all
targetGroup = targetPairs(2).pair;
%Set up figure
f = figure; f.Position = [10 10 300 400];
% fs = 12;
[ha, pos] = tight_subplot(3,1,0.05,0.1,.15);
channelInd = 1;
maxTime = 0;
for channel = [102,13,101]
   axes(ha(channelInd));
   %Populate Figure
   for target = targetGroup
       for posture = 1:3
            traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgFR.traj;
            time = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgFR.timestamps;
            if time(end) > maxTime
                maxTime = time(end);
            end
            hold on
            if target == targetGroup(1)
                plot(time,traj(:,channel),'Color',cmap(posture,:),'LineWidth',2);
            elseif target == targetGroup(2)
                plot(time,traj(:,channel),'--','Color',cmap(posture,:),'LineWidth',2);
            end
            if posture == 3 && target == 5
                title(['Monkey N \newline Ch ',num2str(channel)])
            end
       end
   end
   if channelInd == 2
        ylabel('firing rate (Hz)')
   end
   channelInd = channelInd + 1;
   ax = gca;
   ax.TickDir = 'out';
   xticks([0 500 1000])
   ylimits = ylim;
   if ylimits(1) > 0
       ax.YLim(1) = 0;
   end
   yticks([0 floor(ylimits(2))])
   yticklabels({'0', num2str(floor(ylimits(2)))})
   set(gca,'fontname','arial')
%     set(gca,'fontsize',fs)
end
xlabel('time since go cue (ms)')
xticklabels({'0','500','1000'})
for channelInd = 1:3
    axes(ha(channelInd));
    ax = gca;
    ax.XLim = [-50 maxTime];
end
saveas(gcf,'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Figures\Paper Figures\PSTH\NigelPSTH.svg')