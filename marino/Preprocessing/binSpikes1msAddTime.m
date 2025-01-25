function [Data] = binSpikes1msAddTime(Data)
%Replaces standard spikes struct (channels, sorts, timestamps) with a
%time x neurons matrix binned at 1ms.  Includes whole trial.

numTrials = size(Data,2);
for trial = 1:numTrials
   trial
   %Get trial spikes
   trialSpikes = Data(trial).spikes;
   %Set up bins
   stateTransitions = Data(trial).stateData.stateTransitions;
   startTime = stateTransitions(2,1);
   endTime = stateTransitions(2,end);
   timestamps = startTime+1:1:endTime-1;
   binEdges = [timestamps-0.5,timestamps(end)+0.5];
   %Bin spikes 
   [binnedSpikes] = getBinnedSpikes20211210(trialSpikes,binEdges);
   %Store
   Data(trial).spikes = binnedSpikes; 
   %Add time field
   Data(trial).time = timestamps;
end

 
end