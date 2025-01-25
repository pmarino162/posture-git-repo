function plotDrift(Data,inclStates,date)

% This function plots drifting over the course of an experiment.

%% Get parameters for this dataset
    %Get numTrials, numChannels, and binWidth
    numTrials = size(Data,2);
    i = 1;
    while isempty(Data(i).waveforms) | isempty(Data(i).spikes)
        i = i+1;
    end
    numChannels = size(Data(i).waveforms,2);
    binWidth = double(Data(i).spikes.binTimes(2) - Data(i).spikes.binTimes(1));
    
    %Get indices for inclStates
    numInclStates = size(inclStates,2);
    stateNames = Data(1).stateData.stateNames;
    for state = 1:numInclStates
        stateInd(state) = max(find([cellfun(@(x) strcmpi(x,inclStates{1,state}),stateNames)]==1));
    end 
    
%% Collect Relevant Data from Time Window(s) of Interest    
    %For each trial, collect firing rates during time window(s) of interest.
    %Add these to Data Struct    
    for trial = 1:numTrials
        stateTransitions = Data(trial).stateData.stateTransitions;  
        %Only collect data for trials that contain all inclStates and have
        %a non-empty spikes field
        if sum(ismember(stateInd,stateTransitions(1,:))) == numInclStates & ~isempty(Data(trial).spikes)
            %Get Onset and Offset Times for Window(s) of Interest. Store in
            %stateTimeWindows
            stateTimeWindows = zeros(numInclStates,2);
            for state = 1:numInclStates
                stateStartTime = double(stateTransitions(2,max(find(stateTransitions(1,:)==stateInd(state)))));
                stateEndTime = double(stateTransitions(2,max(find(stateTransitions(1,:)==stateInd(state)))+1));
                stateTimeWindows(state,:) = [stateStartTime,stateEndTime];
            end
            % Find binTimes and timestamps that fall in the window(s) of interest. 
            % Collect Firing Rates, waveforms, and P2P for trial
            for channel = 1:numChannels
                binTimes = Data(trial).spikes.binTimes;
                timestamps = Data(trial).waveforms(channel).allSorts.timestamps;
                binTimesMask = zeros(1,length(binTimes));
                timestampsMask = zeros(1,length(timestamps));
                for state = 1:numInclStates
                    stateStartTime = stateTimeWindows(state,1);
                    stateEndTime = stateTimeWindows(state,2);
                    binTimesMask = (binTimes >= stateStartTime & binTimes <= stateEndTime) | binTimesMask;
                    timestampsMask = (timestamps >= stateStartTime & timestamps <= stateEndTime) | timestampsMask;
                end
                %Firing Rates
                Data(trial).windowFR(channel).FR = Data(trial).spikes.channelSpikes(channel).allSorts(binTimesMask)./binWidth.*1000;
                %Waveforms
                Data(trial).windowWaveforms(channel).waveforms = Data(trial).waveforms(channel).allSorts.waveforms(:,timestampsMask);
                %P2P
                numWaveforms = sum(timestampsMask);
                P2P = zeros(1,numWaveforms);
                waveforms = Data(trial).windowWaveforms(channel).waveforms;
                for i = 1:numWaveforms
                    P2P(i) = max(waveforms(:,i))-min(waveforms(:,i));
                end
                Data(trial).windowP2P(channel).P2P = P2P;
            end
        end
    end

    
%% Collect Data for Plotting 
    %Create plotData Struct
    plotData = struct('meanFR',zeros(1,numTrials),'stdFR',zeros(1,numTrials),'meanP2P',zeros(1,numTrials),'stdP2P',zeros(1,numTrials));
        
    %For each trial, average data over trial window centered on current
    %trial (window size will be closest odd number to specified size) -
    %later, this should be modified to only get trials which have a spikes
    %field
    windowSize = 8;
    for trial = 1:numTrials
        trialList = trial-round(windowSize/2)+1:trial+round(windowSize/2)-1;
        trialListMask = trialList >= 1 & trialList <= numTrials;
        trialList = trialList(trialListMask);
        selectData = Data(trialList);
        for channel = 1:numChannels
            chFR = [];
            chP2P = [];
            for i = 1:size(selectData,2)
                if ~isempty(selectData(i).windowFR) & ~isempty(selectData(i).windowP2P)
                    chFR = [chFR,selectData(i).windowFR(channel).FR];
                    chP2P = [chP2P,selectData(i).windowP2P(channel).P2P];
                end
            end
            plotData(channel).meanFR(trial) = mean(chFR);
            plotData(channel).stdFR(trial) = std(chFR);
            plotData(channel).meanP2P(trial) = mean(chP2P);
            plotData(channel).stdP2P(trial) = std(chP2P);
        end     
    end
    
    

% Plot Results
early = 1;
middle = round(numTrials/2);
late = numTrials;

repTrials = [early,middle,late];

%Firing Rate
f=figure('Position',[0 0 1600 800])
f.Name = [date,'meanFR']
for channel = 1:96
   subplot(10,10,channel)
   meanFR = plotData(channel).meanFR;
   stdFR = plotData(channel).stdFR;
   shadedErrorBar([],meanFR,stdFR,'lineProps',{'LineWidth',3})
   title(['Ch ',num2str(channel)])
end
sgtitle([date,' Mean Firing Rate by Trial (Hz)'])
saveFigurePDF(f,'C:\Users\pmari\Documents\Energy Landscape\Presentations\7-19-19')

%P2P
f=figure('Position',[0 0 1600 800])
f.Name = [date,'meanP2P']
for channel = 1:96
   subplot(10,10,channel)
   meanP2P = plotData(channel).meanP2P;
   stdP2P = plotData(channel).stdP2P;
   shadedErrorBar([],meanP2P,stdP2P,'lineProps',{'LineWidth',3})
   title(['Ch ',num2str(channel)])
end
% linkaxes
sgtitle([date,' Mean P2P Voltage by Trial (\muV)'])
saveFigurePDF(f,'C:\Users\pmari\Documents\Energy Landscape\Presentations\7-19-19')

f=figure('Position',[0 0 1600 800])
f.Name = [date,'stdP2P']
for channel = 1:96
   subplot(10,10,channel)
   stdP2P = plotData(channel).stdP2P;
   plot(stdP2P);
   title(['Ch ',num2str(channel)])
end
sgtitle([date, ' Std Dev of P2P Voltage by Trial (\muV)'])
saveFigurePDF(f,'C:\Users\pmari\Documents\Energy Landscape\Presentations\7-19-19')
% %Waveforms
% for channel = 31
%    f=figure('Position',[0 0 1200 800])
%    f.Name = [date,'Ch',num2str(channel),'Waveforms']
%    for period = 1:3
%        trial = repTrials(period);
%        waveforms = plotData(channel).waveforms(trial).waveforms;
%        numWaveforms(period) = size(waveforms,2);
%    end
%    minNumWaveforms = min(numWaveforms);
%    for period = 1:3
%        subplot(1,3,period)
%        trial = repTrials(period);
%        waveforms = plotData(channel).waveforms(trial).waveforms;
%        plot(waveforms(:,1:minNumWaveforms))
%        ylabel('\muV')
%    end
%    linkaxes
%    sgtitle([date,' Ch ', num2str(channel)])
% end
% saveFigurePDF(f,'C:\Users\pmari\Documents\Energy Landscape\Presentations\7-5-19')


end