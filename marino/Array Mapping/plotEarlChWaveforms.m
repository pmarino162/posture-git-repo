clear; clc; clf; close all; 

%% Setup Save
    savePath = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Presentations\20210624\';
    saveFig = true;
    
%% Load Data 
    dateStr = '20210622';
    switch dateStr
        %6/19/21
        case '20210619'
            load('D:\Animals\Earl\2021\06\20210619\Earl20210619_Delayed Center Out_SI_SW_translated.mat');
            config = 'earlRHM1';
        %6/22/21
        case '20210622'
            load('D:\Animals\Earl\2021\06\20210622\Earl20210622_Delayed Center Out_SI_SW_translated.mat');
            config = 'earlRHM1';
            Data = Data(1,500:end);
        %6/23/21
        case '20210623'
            load('D:\Animals\Earl\2021\06\20210623\Earl20210623_Delayed Center Out_SI_SW_translated.mat');
            config = 'earlRHPMd';
            Data = Data(1,500:end);
        %6/24/21
        case '20210624'
            load('D:\Animals\Earl\2021\06\20210624\Earl20210624_Delayed Center Out_SI_SW_translated.mat');
            config = 'earlRHM1';
            Data = Data(1,500:end);
    end
    %Get Data Struct
    tdtExists = arrayfun(@(x) ~isempty(x.TrialData.TDT.snippetInfo), Data);
    Data = Data(tdtExists);
    numTrials = size(Data,2);
    trialRange  = 1:10:numTrials-1;
    Data = Data(trialRange);
    Data = getDataStruct(Data,'getSpikes',true,'getWaveforms',true); 
    numChannels = size(Data(1).waveforms,2);

%% Load maps
    [tdt2adapterMap, adapter2cereportMap, cereport2arrayMap,array2anotomyMap,tdt2anatomyMap] = getMaps(config);

%% Concatenate waveforms across trials
    numTrials = size(Data,2);
    waveforms = struct('channel',[],'sortID',[],'waveforms',[]);
    structInd = 1;
    for channel = 1:numChannels
       allChWaveforms = [];
       numSorts = size(Data(end).waveforms(channel).sorts,2);
       for sort = 1:numSorts
           allWaveforms = [];
           waveforms(structInd).channel = channel;
           waveforms(structInd).sortID = Data(end).waveforms(channel).sorts(sort).sortID;
           for trial = 1:numTrials   
                allWaveforms = [allWaveforms,Data(trial).waveforms(channel).sorts(sort).waveforms];
                allChWaveforms = [allChWaveforms,Data(trial).waveforms(channel).sorts(sort).waveforms];
           end
           waveforms(structInd).waveforms = allWaveforms';
           structInd = structInd + 1;
       end
       waveforms(structInd).channel = channel;
       waveforms(structInd).sortID = 'all';
       waveforms(structInd).waveforms = allChWaveforms';
       structInd = structInd + 1;
    end
    clearvars Data
 
%% Get Vpp
    Vpp = zeros(1,numChannels);
    for channel = 1:numChannels
        chWaveforms = waveforms([waveforms.channel]==channel);
        chWaveforms = chWaveforms(end).waveforms;
        numWaveforms = size(chWaveforms,1);
        VppMat = zeros(1,numWaveforms);
        if numWaveforms > 250
            for i = 1:numWaveforms
                waveform = chWaveforms(i,:);
                VppMat(i) = (max(waveform)-min(waveform))./10^-6;
            end
        else

        end
        Vpp(channel) = prctile(VppMat,95);
    end

%% Define sort colormap
    sortcmap = struct('sortID',[],'color',[]);
    sortcmap(1).sortID = 1; sortcmap(1).color = [1 1 0];
    sortcmap(2).sortID = 2; sortcmap(2).color = [1 0 0];
    sortcmap(3).sortID = 3; sortcmap(3).color = [0 1 1];
    sortcmap(4).sortID = 4; sortcmap(4).color = [0 1 0];

%% Plot waveforms in anatomical locations
    PMdRows = [cell2mat(tdt2anatomyMap(:,2))<=8];
    PMdChList = [cell2mat(tdt2anatomyMap(PMdRows,1))];
    M1ChList = [cell2mat(tdt2anatomyMap(~PMdRows,1))];
    f = figure; f.Position = [10 10 700 800];
    
    [ha, pos] = tight_subplot(16,8,0,0,0);

    for array = {'PMd','M1'}
        if strcmpi(array,'PMd')
            channelList = PMdChList;
        else
            channelList = M1ChList;
        end
        for channel = channelList'
            chData = waveforms([waveforms.channel]==channel);
            sortIDs = [chData(1:end-1).sortID];
            row = tdt2anatomyMap{channel,2};
            col = tdt2anatomyMap{channel,3};
            axes(ha((row-1)*8 + col));
            for sortID = sortIDs
                if sortID ~= 0 && sortID ~= 31
                    %Get 10 evenly spaced waveform inds
                    sortWaveforms = chData([chData.sortID]==sortID).waveforms;
                    numWaveforms = size(sortWaveforms,1);
                    sortWaveforms = sortWaveforms(1:round(numWaveforms/10):end,:);
                    numWaveforms = size(sortWaveforms,1);
                    %Plot each one
                    color = sortcmap([sortcmap.sortID]==sortID).color;
                    for i = 1:numWaveforms
                        plot(sortWaveforms(i,:),'Color',color)
                        hold on
                    end

                else
                end
            end
            set(gca,'Color','k')
            yl = ylim;
            xl = xlim;
%             text(0,yl(1)+(yl(2)-yl(1))*.7,['Ch' ,num2str(channel),newline,'Vpp=',num2str(round(Vpp(channel)))],'Color','white','FontSize',7)
            text(0,yl(1)+(yl(2)-yl(1))*.9,['Ch' ,num2str(channel)],'Color','white','FontSize',7)
            text(xl(1)+(xl(2)-xl(1))*.5,yl(1)+(yl(2)-yl(1))*.9,['Vpp=',num2str(round(Vpp(channel)))],'Color','white','FontSize',7)
            xticks([])
            yticks([])
        end
    end
    if saveFig 
        set(gcf, 'InvertHardCopy', 'off');
        saveas(gcf,[savePath,'WaveformsAnatomical',dateStr,'.jpg'])
    end
%% Plot waveforms in TDT locations
    f = figure; f.Position = [10 10 1000 500];
    [ha, pos] = tight_subplot(10,10,0,0,0);

    for channel = 1:numChannels
        chData = waveforms([waveforms.channel]==channel);
        sortIDs = [chData(1:end-1).sortID];
        axes(ha(channel));
        for sortID = sortIDs
            if sortID ~= 0 && sortID ~= 31
                %Get 10 evenly spaced waveform inds
                sortWaveforms = chData([chData.sortID]==sortID).waveforms;
                numWaveforms = size(sortWaveforms,1);
                sortWaveforms = sortWaveforms(1:round(numWaveforms/10):end,:);
                numWaveforms = size(sortWaveforms,1);
                %Plot each one
                color = sortcmap([sortcmap.sortID]==sortID).color;
                for i = 1:numWaveforms
                    plot(sortWaveforms(i,:),'Color',color)
                    hold on
                end
            else
            end
        end
        set(gca,'Color','k')
        yl = ylim;
        text(0,yl(1)+(yl(2)-yl(1))*.9,['Vpp=',num2str(round(Vpp(channel)))],'Color','white','FontSize',7)
        xticks([])
        yticks([])
    end
    if saveFig 
        set(gcf, 'InvertHardCopy', 'off');
        saveas(gcf,[savePath,'WaveformsTDT',dateStr,'.jpg'])
    end
    
%% Make Vpp Array Map
    createArrayImage(config,Vpp,dateStr)
    if saveFig 
        saveas(gcf,[savePath,'ArrayMapVpp',dateStr,'.jpg'])
    end