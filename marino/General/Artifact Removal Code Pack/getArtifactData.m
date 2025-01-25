function [artifactData] = getArtifactData(trialData)

%Create artifactData struct
artifactData = struct('timestamps',[],'neuralData',[]);

%Identify timestamps from spikeSum
spikeSum = trialData.spikeSum;
numPts = size(spikeSum,1);
i=1;
artifactInd = 1;
threshold = 30;
while i<numPts
    if spikeSum(i,1) > threshold
        %Save timestamps
        startTime = i;
        while i+51 <= numPts & sum(spikeSum(i+1:i+51) > threshold) > 0
            i = i+1;
        end
        if i+51 >= numPts
            i=numPts
        end
        endTime = i;
        i=i+1;
        artifactData(artifactInd).timestamps(1,1) = startTime;
        artifactData(artifactInd).timestamps(1,2) = endTime;
        %Save neural data
        decoderTime = trialData.Decoder.timestamps;
        timestampMask = decoderTime > startTime & decoderTime < endTime;
        artifactData(artifactInd).neuralData = trialData.Decoder.GPFATraj(timestampMask,:);
        artifactInd = artifactInd + 1;
    else
        i = i+1;
    end
end

end