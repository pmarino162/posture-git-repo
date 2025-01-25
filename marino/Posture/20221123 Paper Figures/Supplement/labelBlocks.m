function [Data] = labelBlocks(Data)

%Function for labeling separate experimental blocks for the posture/drift
%dissociation control

block = 1;
numTrials = size(Data,2);
curPosture = Data(1).conditionData.postureID;
for trial = 1:numTrials
   trialPosture = Data(trial).conditionData.postureID;
   if trialPosture ~= curPosture
       curPosture = trialPosture;
       block = block + 1;
   end
   Data(trial).block = block;
   
end


end