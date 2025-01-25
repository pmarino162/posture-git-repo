function [artifactTrials] = detectArtifact(Data)
% This function labels trials containing the artifact caused by the
% condition of Earl's array in 2019-2020.

%% Set threshold value and number of bins necessary for exlclusion
    threshold = 30;
    numBinsNec = 10;
    
%% Load data
    Data = Data.Data;

%% Preprocess data
    procData = preprocessDataMarino(Data,'CHECKPHASESPACESYNC',false,'CHECKDECODE',false,'REMOVESCREENFREEZETRIALS',false);
    numTrials = size(procData,2);
    
%% Get number of channels and sorts; create channelSpikes struct
    spikes = procData(1).TrialData.spikes;
    channelList = unique([spikes.channel]);
    numChannels = length(channelList);
    numSorts = zeros(1,numChannels);
    for channel = 1:numChannels
        chSpikes = spikes([spikes.channel]==channel);
        numSorts(channel) = size(chSpikes,2);   
    end
    
%% For each trial, check for artifact and forcebar failures.  Plot sum of spikes for each trial
    artifactTrials = [];
    for trial = 1:numTrials
        %Get Binned Spikes
        spikes = procData(trial).TrialData.spikes;
        stateTransitions = procData(trial).TrialData.stateTransitions;
        [binTimes,allSpikeBins] = getBinnedSpikesArtifactDetection(spikes,stateTransitions,channelList,numChannels,numSorts);
        %Sum Across The Array 
        spikeSum = sum(double(allSpikeBins),2);
% %         Plot Results
%         if mod(trial,50)==1
%             f=figure
%             f.Position = [0 0 1750 1000];
%         end
%         if mod(trial,50) == 0
%             subplot(ceil(sqrt(50)),ceil(sqrt(50)),50) 
%         else
%             subplot(ceil(sqrt(50)),ceil(sqrt(50)),mod(trial,50)) 
%         end
%         plot(spikeSum)
%         hold on
%         ax = gca;
%         plot(ax.XLim,[threshold,threshold],'r')
%         xlabel('time (ms)')
%         ylabel('# spikes')
%         ylim([0 50])
%         title(['Trial ',num2str(trial)])
        %Compare to threshold value and flag trials containing artifact
        if sum(spikeSum > threshold) >= numBinsNec
            artifactTrials = [artifactTrials,trial];
        end
    end
%     
    
    
end