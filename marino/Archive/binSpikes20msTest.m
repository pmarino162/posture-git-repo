function [Data] = binSpikes20msTest(Data)
%Replaces standard spikes struct (channels, sorts, timestamps) with a
%time x neurons matrix binned at 1ms.  Includes whole trial.

binWidth = 20;

numTrials = size(Data,2);
for trial = 1:numTrials
   trial
   %Get trial spikes
   trialSpikes = Data(trial).spikes;
   %Set up bins
   stateTransitions = Data(trial).stateData.stateTransitions;
   startTime = stateTransitions(2,1);
   endTime = stateTransitions(2,end);
   timestamps = startTime:binWidth:endTime;
   binEdges = [timestamps-binWidth/2,timestamps(end)+binWidth/2];
   %Bin spikes 
   [binnedSpikes] = getBinnedSpikes20211210(trialSpikes,binEdges);
   %Store
   Data(trial).spikes = binnedSpikes; 
   %Add time field
   Data(trial).time = timestamps;
end

 
end