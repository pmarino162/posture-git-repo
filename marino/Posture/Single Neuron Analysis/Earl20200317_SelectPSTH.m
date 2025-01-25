clear; clc; clf; close all;

%% Load Data 
    %Earl - 3/17/2020
    date = '20200317';
    exclCh =  [3 12 16 20 31 73 85 92 94];
    binWidth = 1;
    [Data,I30GPFAParams,I30DecoderParams,I15GPFAParams,I15DecoderParams,...
     N00GPFAParams,N00DecoderParams,E15GPFAParams,E15DecoderParams,...
     E30GPFAParams,E30DecoderParams] = loadEarlData20200317(binWidth);

 
%% Collect States of Interest into Traj Field
    trialInclStates = struct('trialName','','inclStates',[],'inclOccurrence',[]);
    trialInclStates(1).trialName = {'GridTask_CO_Across_BC_ForceBar'};
    trialInclStates(1).inclStates = {'Step 1'}; trialInclStates(1).inclOccurrence = {'last'};
    trialInclStates(1).addTimeToBeginning = {-500};
    
    numTrials = size(Data,2);
    for trial = 1:numTrials
        [FR,FRTimestamps] = getStatesTraj(Data(trial),trialInclStates,'allChannelSmoothedFR');
        Data(trial).Traj.FR = FR; Data(trial).Traj.FRTimestamps = FRTimestamps;
    end
    clearvars trialInclStates
    
    
%% Keep only trials with length w/in 2 stdDevs of mean length
    [Data] = excludeLengths(Data,'FR');
    
%% Create trajStruct
    condFields = {{'posture','postureData','postureID'},{'target','targetData','target1ID'}};
    trajFields = {'FR'};
    trajStruct = getTrajStruct(Data,condFields,trajFields,'addTimeToBeginning',-100);
    
%% Channel Grouping Struct
    channelGroups = struct('channels',[]);
    channelGroups(1).channels = 1:50; channelGroups(2).channels = 51:87;

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
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customSummer.mat');
    cmap = customSummer;
%% Plotting - try tight subplot     
close all
targetGroup = targetPairs(2).pair;
%Set up figure
f = figure; f.Position = [10 10 300 400];
% fs = 12;
[ha, pos] = tight_subplot(3,1,0.05,0.1,.15);
channelInd = 1;
maxTime = 0;
for channel = [3,34,67]
   axes(ha(channelInd));
   %Populate Figure
   for target = targetGroup
       for posture = 1:5
            traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgFR.traj;
            time = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgFR.timestamps;
            if time(end) > maxTime
                maxTime = time(end)
            end
            hold on
            if target == targetGroup(1)
                plot(time,traj(:,channel),'Color',cmap(posture,:),'LineWidth',2);
            elseif target == targetGroup(2)
                plot(time,traj(:,channel),'--','Color',cmap(posture,:),'LineWidth',2);
            end
            if posture == 5 && target == 5
                title(['Monkey E \newline Ch ',num2str(channel)])
            end
       end
   end
   if channelInd == 2
        ylabel('threshold crossings (Hz)')
   end
   channelInd = channelInd + 1;
   ax = gca;
   ax.TickDir = 'out';
   xticks([-100 0 250])
   ylimits = ylim;
   if ylimits(1) > 0
       ax.YLim(1) = 0;
   end
   yticks([0 ylimits(2)])
   yticklabels({'0', num2str(ylimits(2))})
   set(gca,'fontname','arial')
%     set(gca,'fontsize',fs)
end
xlabel('time since go cue (ms)')
xticklabels({'-100','0','250'})
for channelInd = 1:3
    axes(ha(channelInd));
    ax = gca;
    ax.XLim = [-100 maxTime];
end
saveas(gcf,'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Figures\Paper Figures\PSTH\EarlPSTH.svg')