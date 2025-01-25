clear; clc; clf; close all; 

%% Load Data 
dateStr = '20210622';
switch dateStr
    %6/22/21
    case '20210622'
        load('D:\Animals\Earl\2021\06\20210619\Earl20210619_Delayed Center Out_SI_SW_translated.mat')
        tdtExists = arrayfun(@(x) ~isempty(x.TrialData.TDT.snippetInfo), Data);
        Data = Data(tdtExists);
        numTrials = size(Data,2);
        trialRange  = 1:10:numTrials-1;
        Data = Data(trialRange);
%         config = 'PMd';    
end

%Get Data Struct
    Data = getDataStruct(Data,'getSpikes',true,'getWaveforms',true); 
    numCh = size(Data(1).waveforms,2);

%% Preallocate
    for channel = 1:numCh
       allWaveforms(channel).waveforms = [];
    end
    
%% Collect All Waveforms from Trial Range  
    numTrials = size(Data,2);
    for trial = 1:numTrials
        for channel = 1:numCh
            if ~isempty(Data(trial).waveforms(channel).allSorts)
            channelWaveforms = [Data(trial).waveforms(channel).allSorts.waveforms]';
            end
            allWaveforms(channel).waveforms = vertcat(allWaveforms(channel).waveforms,channelWaveforms);
        end
    end
        
%% Get Vpp
    Vpp = zeros(1,numCh);
    meanWaveform = zeros(numCh,30);
    for channel = 1:numCh
        numWaveforms = size(allWaveforms(channel).waveforms,1);
        VppMat = zeros(1,numWaveforms);
        if numWaveforms > 250
            for i = 1:numWaveforms
                waveform = allWaveforms(channel).waveforms(i,:);
                VppMat(i) = (max(waveform)-min(waveform))./10^-6;
            end
            meanWaveform(channel,:) = mean(allWaveforms(channel).waveforms,1);
        else

        end
        Vpp(channel) = prctile(VppMat,95);
    end

%% Plot All Waveforms
    savePath = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Earl Array Characterization';
    f = figure
    for channel = 1:numCh
        channel
        subplot(10,10,channel)
        numWaveforms = size(allWaveforms(channel).waveforms,1);
        if numWaveforms >= 250
            minVal = 1000;
            maxVal = -1000;
            for i = 1:25:numWaveforms
                waveform = allWaveforms(channel).waveforms(i,:)./(10^-6);
                if min(waveform) < minVal
                    minVal = min(waveform);
                end
                if max(waveform) > maxVal
                    maxVal = max(waveform);
                end
                plot(waveform,'Color',[0.8 0.8 0.8]);
                hold on;
            end
            plot(meanWaveform(channel,:)./(10^-6),'k','LineWidth',1.5)       
        

        ylim([round(minVal) round(maxVal)])
        yticks([round(minVal) round(maxVal)])
        title(['Ch ', num2str(channel),', Vpp = ',num2str(round(Vpp(channel)))])
        end
                ylabel('\muV')
        xticklabels({})
        title(['Ch ', num2str(channel)])

    end
    f.Position = [0 0 1750 1000]
    sgtitle(dateStr)
%     saveas(gcf,[savePath,'Waveforms',dateStr,'.jpg'])
%    linkaxes
%% Plot Vpp by Channel
    figure
    bar(Vpp)
    xlabel('Channel')
    ylabel('Vpp (\muV)')
%     saveas(gcf,['Bar',saveString,'.jpg'])
%     figure
%     hist(Vpp)
%     saveas
    
%% Return Channels Below Threshold
    threshold = 38;
    lowAch = find(Vpp < threshold);
    highAch = find(Vpp > threshold);
    
%     close all