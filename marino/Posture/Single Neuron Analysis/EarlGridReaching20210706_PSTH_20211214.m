clear; clc; clf; close all;
    
%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customSummer.mat');
    cmap = customSummer;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat');
    tcmap = customRainbow;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\rainbow2d.mat')
    pcmap = rainbow2d;
    
%% Load Data 
% [44 87 88 77 78 71 67 69 118]
    load('D:\Animals\Earl\2021\07\20210706\Earl20210706_gridReaching_SI_translated.mat');
    Data = getDataStruct20211210(Data,'getSpikes',true,'getSorts',false,'exclCh',[],'exclZero',false,'getMarker',true,'centerMarker',false,'getKin',true);
    [Data] = cleanData20210706(Data);
    Data = Data([Data.trialStatus]==1);
    [Data,postureIDs] = labelPostures20210706(Data);
    
%% Get trajStruct 
    trajFields = {'smoothFR'};
    trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
    condFields = {{'target','targetData','targetID'},{'posture','conditionData','postureID'}};
    trialInclStates(1).trialName = {'GridReaching'};
    binWidth = 1;
    %Delay
    kinData = [Data.kinData]; delayLength = [kinData.delayLength];
    tempData = Data(delayLength >= 500);
    trialInclStates(1).inclStates = {{'state','Delay','first',-100},{'state','Target Acquire','first',0}};
    delayTrajStruct = getTrajStruct20211210(tempData,condFields,trajFields,trialInclStates,binWidth,'matchConditions',false);
    clearvars tempData
    %Reach
    trialInclStates(1).inclStates = {{'kin','moveOnsetTime','first',-100},{'state','Target Hold','first',0}};
    reachTrajStruct = getTrajStruct20211210(Data,condFields,trajFields,trialInclStates,binWidth,'matchConditions',false);
    %Target hold 
    trialInclStates(1).inclStates = {{'state','Target Hold','first',0},{'state','Success with Reward','first',0}};
    holdTrajStruct = getTrajStruct20211210(Data,condFields,trajFields,trialInclStates,binWidth,'matchConditions',false);
    
%% Get posture and target lists
    trajStruct = reachTrajStruct;
    postureList = unique([trajStruct.posture]); numPostures = size(postureList,2);
    targetList = unique([trajStruct.target]); numTargets = size(targetList,2);
   
%% Plot PSTH
    figure
      hold on
    channel = 5;
    cmap = customRainbow;
    reachTimeOffset = 1100;
    holdTimeOffset = 800;
    for target = 1:8
        for posture = 1
        %Delay
        traj = delayTrajStruct([delayTrajStruct.target]==target & [delayTrajStruct.posture]==posture).avgSmoothFR.traj(:,channel);
        timestamps = delayTrajStruct([delayTrajStruct.target]==target & [delayTrajStruct.posture]==posture).avgSmoothFR.timestamps;
        plot(timestamps,traj,'Color',cmap(target,:),'LineWidth',1)
        
        %Reach
        traj = reachTrajStruct([reachTrajStruct.target]==target & [reachTrajStruct.posture]==posture).avgSmoothFR.traj(:,channel);
        timestamps = reachTrajStruct([reachTrajStruct.target]==target & [reachTrajStruct.posture]==posture).avgSmoothFR.timestamps+reachTimeOffset;
        plot(timestamps,traj,'Color',cmap(target,:),'LineWidth',1)
        
        %Hold 
        traj = holdTrajStruct([holdTrajStruct.target]==target & [holdTrajStruct.posture]==posture).avgSmoothFR.traj(:,channel);
        timestamps = holdTrajStruct([holdTrajStruct.target]==target & [holdTrajStruct.posture]==posture).avgSmoothFR.timestamps+reachTimeOffset+holdTimeOffset;
        plot(timestamps,traj,'Color',cmap(target,:),'LineWidth',1)
        
        end
    end
    xlabel('time (ms)')
    ylabel('FR (Hz)')
    xlim([-250 reachTimeOffset+holdTimeOffset+600])
    xticks([0 500 reachTimeOffset+100 reachTimeOffset+600 reachTimeOffset+holdTimeOffset])
    xticklabels({'TO','500','MO','500','TH'})

%% All ch - target tuning
for channelGroup = 1:2
    if channelGroup == 1
        channelList = 1:64;
    else
        channelList = 65:128;
    end
    f = figure; f.Position = [5 5 1400 700];
    [ha, pos] = tight_subplot(8,8,0.05,0.05,.05);
    channelInd = 1;
    for channel = channelList
       axes(ha(channelInd));
       %Populate Figure
       for target = 1:6
           for posture = 5
                %Delay
                traj = delayTrajStruct([delayTrajStruct.target]==target & [delayTrajStruct.posture]==posture).avgSmoothFR.traj(:,channel);
                timestamps = delayTrajStruct([delayTrajStruct.target]==target & [delayTrajStruct.posture]==posture).avgSmoothFR.timestamps;
                plot(timestamps,traj,'Color',cmap(target,:),'LineWidth',1)
                hold on
                %Reach
                traj = reachTrajStruct([reachTrajStruct.target]==target & [reachTrajStruct.posture]==posture).avgSmoothFR.traj(:,channel);
                timestamps = reachTrajStruct([reachTrajStruct.target]==target & [reachTrajStruct.posture]==posture).avgSmoothFR.timestamps+reachTimeOffset;
                plot(timestamps,traj,'Color',cmap(target,:),'LineWidth',1)

                %Hold 
                traj = holdTrajStruct([holdTrajStruct.target]==target & [holdTrajStruct.posture]==posture).avgSmoothFR.traj(:,channel);
                timestamps = holdTrajStruct([holdTrajStruct.target]==target & [holdTrajStruct.posture]==posture).avgSmoothFR.timestamps+reachTimeOffset+holdTimeOffset;
                plot(timestamps,traj,'Color',cmap(target,:),'LineWidth',1)
           end
       end
       channelInd = channelInd + 1;
       ax = gca;
       ax.TickDir = 'out';
       
        xlim([-250 reachTimeOffset+holdTimeOffset+600])
        xticks([0 500 reachTimeOffset+100 reachTimeOffset+holdTimeOffset])
        xticklabels({'TO','500','MO','TH'})
       xtickangle(45)
       
       ylimits = ylim;
       if ylimits(1) > 0
           ax.YLim(1) = 0;
       end
       yticks([0 ylimits(2)])
       yticklabels({'0', [num2str(round(ylimits(2))),' Hz']})
       set(gca,'fontname','arial')
       title(['Ch ',num2str(channel)])
    end
end


%% All ch - posture tuning
for channelGroup = 1:2
    if channelGroup == 1
        channelList = 1:64;
    else
        channelList = 65:128;
    end
    f = figure; f.Position = [5 5 1400 700];
    [ha, pos] = tight_subplot(8,8,0.05,0.05,.05);
    channelInd = 1;
    for channel = channelList
       axes(ha(channelInd));
       %Populate Figure
       for target = 3
           for posture = [1:4]
                %Delay
                traj = delayTrajStruct([delayTrajStruct.target]==target & [delayTrajStruct.posture]==posture).avgSmoothFR.traj(:,channel);
                timestamps = delayTrajStruct([delayTrajStruct.target]==target & [delayTrajStruct.posture]==posture).avgSmoothFR.timestamps;
                plot(timestamps,traj,'Color',tcmap(posture,:),'LineWidth',1)
                hold on
                %Reach
                traj = reachTrajStruct([reachTrajStruct.target]==target & [reachTrajStruct.posture]==posture).avgSmoothFR.traj(:,channel);
                timestamps = reachTrajStruct([reachTrajStruct.target]==target & [reachTrajStruct.posture]==posture).avgSmoothFR.timestamps+reachTimeOffset;
                plot(timestamps,traj,'Color',tcmap(posture,:),'LineWidth',1)

                %Hold 
                traj = holdTrajStruct([holdTrajStruct.target]==target & [holdTrajStruct.posture]==posture).avgSmoothFR.traj(:,channel);
                timestamps = holdTrajStruct([holdTrajStruct.target]==target & [holdTrajStruct.posture]==posture).avgSmoothFR.timestamps+reachTimeOffset+holdTimeOffset;
                plot(timestamps,traj,'Color',tcmap(posture,:),'LineWidth',1)
           end
       end
       channelInd = channelInd + 1;
       ax = gca;
       ax.TickDir = 'out';
       
        xlim([-250 reachTimeOffset+holdTimeOffset+600])
        xticks([0 500 reachTimeOffset+100 reachTimeOffset+holdTimeOffset])
        xticklabels({'TO','500','MO','TH'})
       xtickangle(45)
       
       ylimits = ylim;
       if ylimits(1) > 0
           ax.YLim(1) = 0;
       end
       yticks([0 ylimits(2)])
       yticklabels({'0', [num2str(round(ylimits(2))),' Hz']})
       set(gca,'fontname','arial')
       title(['Ch ',num2str(channel)])
    end
end