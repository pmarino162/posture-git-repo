function Data = removeNoSpikeTrials(Data)
    noSpikeTrials = [];
    numTrials = size(Data,2);
    for trial = 1:numTrials
       if isempty(Data(trial).spikes)
           noSpikeTrials = [noSpikeTrials,trial];
       else
           spikeTimestamps = [Data(trial).spikes.timestamps];
           if size(spikeTimestamps,2) < 2
               noSpikeTrials = [noSpikeTrials,trial];
           end
       end
    end
    Data(noSpikeTrials) = [];


end