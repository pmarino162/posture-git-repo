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
    trialInclStates(1).addTimeToBeginning = {-250};
    
    numTrials = size(Data,2);
    for trial = 1:numTrials
        [FR,FRTimestamps] = getStatesTraj(Data(trial),trialInclStates,'allChannelSmoothedFR');
    end
    clearvars trialInclStates
    
    
%% Keep only trials with length w/in 2 stdDevs of mean length
    [Data] = excludeLengths(Data,'FR');
    
%% Create trajStruct
    condFields = {{'posture','postureData','postureID'},{'target','targetData','target1ID'}};
    trajFields = {'FR'};
    trajStruct = getTrajStruct(Data,condFields,trajFields);
    
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
%% Plotting 
close all
for channelGroup = 1:2
   for target = 1
       %Get channel group
       channelList = channelGroups(channelGroup).channels;
       %Set up figure
       f = figure; f.Position = [0 0 1500 750];
       M = 5; N = 10;
       p = panel(); p.pack({1/14 [] 1/13}, {1/50 []}); p(2,2).pack(M,N); p.margin = [0 0 0 0];
       %Populate Figure
       for posture = 1:3
            traj = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgFR.traj;
            time = trajStruct([trajStruct.posture]==posture & [trajStruct.target]==target).avgFR.timestamps;
            for channel = channelList
                m = floor((channel-1)/N)+1;
                if m > 5
                    m = m-5;
                end
                n = channel-N*(m-1);
                if n > 50
                    n = n-50;
                end
                p(2,2,m,n).select();
                plot(time,traj(:,channel),'Color',cmap(posture,:),'LineWidth',2)
                box off            
                if posture == 5
                    ax = gca;
                    XLim = ax.XLim;
                    minXval = floor(XLim(1,1)); maxXval = ceil(XLim(1,end));
                    xRange = maxXval-minXval;
                    xMidPt = minXval + xRange./2;
                    YLim = ax.YLim;
                    minYval = floor(YLim(1,1)); maxYval = ceil(YLim(1,end));
                    yRange = maxYval-minYval;
                    if m==5
                        xTickLabelOffset = xRange*.2;
                        xlabel('time (ms)')
                        ax.XTick = [minXval,round(maxXval-xTickLabelOffset)];
                        ax.XTickLabel = {minXval,round(maxXval-xTickLabelOffset)};
                        xtickangle(45)
                    else
                        set(gca,'Xticklabel',[])
                    end
                    if mod(channel-1,10)==0
                        ylabel('Firing Rate (Hz)') 
                    end           
                    labelOffset = .1*yRange;
                    limitOffset =  .2*yRange;
                    ax.YLim = [minYval-limitOffset,maxYval+limitOffset];
                    ax.YTick = [round(minYval-labelOffset),round(maxYval+labelOffset)];
                    set(gca,'Yticklabel',[]) 
                    text(20,round(minYval-labelOffset),num2str(round(minYval-labelOffset)),'FontSize',9)
                    text(20,round(maxYval+labelOffset),num2str(round(maxYval+labelOffset)),'FontSize',9)                     
                    text(xMidPt-60,maxYval+labelOffset,['Ch ',num2str(channel)])
                end
                hold on
            end
       end
                titleStr = ['Target ',num2str(target),' Ch ',num2str(channelList(1,1)),'-',num2str(channelList(1,end))];
                sgtitle(titleStr);

%             filePath = 'Figures\PSTH\conditionPSTH_allPosturesForEachTarget\';
%             saveas(gcf,[filePath,titleStr,'.jpg']);
%             close all
   end
end