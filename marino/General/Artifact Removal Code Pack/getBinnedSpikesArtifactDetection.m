function [binTimes,allSpikeBins] = getBinnedSpikes(spikes,stateTransitions,channelList,numChannels,numSorts)
   
    %This function bins spikes for each sort on each channel.  It also
    %combines bins for all sorts on each channel.  It uses the histcounts
    %function for binning, so some extra steps have to be
    %taken to work around this function's limitations
    
    %% Set up bins
    startTime = double(stateTransitions(2,1));
    endTime = double(stateTransitions(2,end));
    binWidth = 1;
    binEdges = startTime:binWidth:endTime;
    binTimes = uint32(binEdges(2:end));
    numBins = size(binTimes,2);
    allSpikeBins = zeros(numBins,numChannels,'uint16');
    %Before running, add extra last bin so that histcounts won't include both
    %edges of final bin
    binEdges = [binEdges,binEdges(end)];
    %If you want to include the right bin edge, add 1 to each bin edge
    %(this works since timestamps are always integers).  Otherwise, bins
    %will be left-edged
    binEdges = binEdges + 1;
    
    %% For each channel, bin spikes for each sort and for all sorts 
    if ~isempty(spikes)
        for channel = 1:numChannels
            %Get channel data
            chSpikes = spikes([spikes.channel]==channel);
            totalChSpikes = zeros(1,numBins,'uint16');
            %Bin spikes for each sort.  Count up channel total
            for sortInd = 1:numSorts(channel)
               sortID = chSpikes(sortInd).sort;
               timestamps = double(chSpikes(sortInd).timestamps);
               binnedSpikeMat = uint16(histcounts(timestamps,binEdges));
               %Delete extra last bin
               binnedSpikeMat(end) = [];
               %Fill Struct
               totalChSpikes = totalChSpikes + binnedSpikeMat;
               %Spike Counts
               allSpikeBins(:,channel) = totalChSpikes';
            end
        end
    end
end