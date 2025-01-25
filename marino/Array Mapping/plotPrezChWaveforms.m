clear; clc; clf; close all; 

%% Load Data 
dateStr = '20210524 M1';
switch dateStr
    %5/24/21
    case '20210524 M1'
        load('D:\Animals\Prez\20210524\Prez20210524_20210524_SI_SW_translated.mat');
        tdtExists = arrayfun(@(x) ~isempty(x.TrialData.TDT.snippetInfo), Data);
        Data = Data(tdtExists);
        numTrials = size(Data,2);
        trialRange  = 1:10:numTrials;
        config = 'M1';

    %5/26/21
    case '20210526'
        load('D:\Animals\Prez\20210526\Prez20210526_20210526_SI_SW_translated_ms.mat')
        numTrials = size(Data,2);
        trialRange  = 1:40:numTrials;
        config = 'PMd';
%     %5/27/21
%     dateStr = '20210527';
%     load('D:\Animals\Prez\20210527\Prez20210527_20210527_SI_SW_translated.mat');
%     Data(1:35) = [];
    
%     %5/28/21
%     dateStr = '20210528';
%     load('D:\Animals\Prez\20210528\Prez20210528_20210528_SI_SW_translated.mat')
%     config = 'PMd';
%     numTrials = size(Data,2);
%     trialRange  = 1:40:numTrials;
    
% %5/28/21 - M1 Config
%     dateStr = '20210528 M1';
%     load('D:\Animals\Prez\20210528 M1 Config Test\Prez20210528_M1 config test 20210528_SI_SW_translated.mat')
%     config = 'M1';
%     trialRange = 1:8;
    
    
    
    case '20210530'
        load('D:\Animals\Prez\20210530\Prez20210530_20210530_SI_SW_translated.mat');
        tdtExists = arrayfun(@(x) ~isempty(x.TrialData.TDT.snippetInfo), Data);
        Data = Data(tdtExists);
        numTrials = size(Data,2);
        trialRange  = 1:10:numTrials-1;
        config = 'PMd';
        
     case '20210530 M1'
        load('D:\Animals\Prez\20210530 M1 Config Test\Prez20210530_M1 config test 20210530_SI_SW_translated.mat')
        tdtExists = arrayfun(@(x) ~isempty(x.TrialData.TDT.snippetInfo), Data);
        Data = Data(tdtExists);
        numTrials = size(Data,2);
        trialRange  = 1:numTrials-1;
        config = 'M1';
    
end
%Get Data Struct
    Data(end) = [];
    Data = getDataStruct(Data,'getSpikes',true,'getWaveforms',true); 
    

%% Preallocate
    for channel = 1:96
       allWaveforms(channel).waveforms = [];
    end
    
%% Collect All Waveforms from Trial Range  
    for trial = trialRange
        for channel = 1:96
            if ~isempty(Data(trial).waveforms(channel).allSorts)
            channelWaveforms = [Data(trial).waveforms(channel).allSorts.waveforms]';
            end
            allWaveforms(channel).waveforms = vertcat(allWaveforms(channel).waveforms,channelWaveforms);
        end
    end
        
%% Get Vpp
    Vpp = zeros(1,96);
    meanWaveform = zeros(96,30);
    for channel = 1:96
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
    savePath = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Choking\Prez Array Characterization';
    f = figure
    for channel = 1:96
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